MODULE mod_chem_spack_rodas3_dyndt

  USE mem_chem1, ONLY: &
       chem1_vars        ! Type

  USE mem_spack, ONLY: &
       spack_type,     & ! Type
       spack_type_2d     ! Type

  USE Extras, ONLY: &
       ext3d             ! Type

!  USE solve_sparse, ONLY: &
!       Solve_linear      & ! Subroutine, OUT: spack_2d(ijk,inob)%DLk1

  USE mod_chem_spack_ros, ONLY: &
       get_always_zeros, &  ! Subroutine
       get_number_nonzeros  ! Subroutine

  USE mod_chem_spack_jacdchemdc, ONLY: &
       jacdchemdc        ! Subroutine

  USE mod_chem_spack_kinetic, ONLY: &
       kinetic           ! Subroutine

  USE mod_chem_spack_fexchem, ONLY: &
       fexchem           ! Subroutine


  use node_mod, only:  &
       mynum

  use mem_grid, only: &
       grid_g

  use chem1_list, only: &
      spc_name

  IMPLICIT NONE

  PRIVATE


  PUBLIC chem_rodas3_dyndt


CONTAINS

  !========================================================================================
  SUBROUTINE chem_rodas3_dyndt(m1,m2,m3,dtlt,pp,pi0,theta,rv,dn0,cosz,rcp, &
                               nspecies,nr,nr_photo,weight,PhotojMethod,maxnspecies, &
                               nspecies_chem_transported,nspecies_chem_no_transported, &
                               transp_chem_index,no_transp_chem_index,chem1_g,nob,maxblock_size, &
                               block_end,indexk,indexi,indexj,kij_index,last_accepted_dt, &
                               spack,spack_2d,atol,rtol,get_non_zeros,&
                               N_DYN_CHEM,split_method,cp,cpor,p00,na_extra3d,extra3d,&
                               jphoto,att)

use mem_grid, only: time
    IMPLICIT NONE

    INTEGER , INTENT(IN) :: m1
    INTEGER , INTENT(IN) :: m2
    INTEGER , INTENT(IN) :: m3

    ! mem_grid
    REAL             , INTENT(IN) :: dtlt

    ! mem_basic
    REAL             , INTENT(IN) :: pp(m1,m2,m3)
    REAL             , INTENT(IN) :: pi0(m1,m2,m3)
    REAL             , INTENT(IN) :: theta(m1,m2,m3)
    REAL             , INTENT(IN) :: rv(m1,m2,m3)
    REAL             , INTENT(IN) :: dn0(m1,m2,m3)

    ! mem_radiate
    REAL             , INTENT(IN) :: cosz(m2,m3)

    ! mem_micro
    REAL             , INTENT(IN) :: rcp(m1,m2,m3)

    ! chem1_list
    INTEGER          , INTENT(IN) :: nspecies
    INTEGER          , INTENT(IN) :: nr_photo
    INTEGER          , INTENT(IN) :: nr
    REAL             , INTENT(IN) :: weight(nspecies)
    CHARACTER(LEN=10), INTENT(IN) :: PhotojMethod

    ! mem_chem1
    INTEGER          , INTENT(IN) :: maxnspecies
    INTEGER          , INTENT(IN) :: nspecies_chem_transported
    INTEGER          , INTENT(IN) :: nspecies_chem_no_transported
    INTEGER          , INTENT(IN) :: transp_chem_index(maxnspecies)
    INTEGER          , INTENT(IN) :: no_transp_chem_index(maxnspecies)
    CHARACTER(LEN=20), INTENT(IN) :: split_method
    INTEGER          , INTENT(IN) :: N_DYN_CHEM

    ! mem_chem1
    TYPE (chem1_vars), INTENT(INOUT) :: chem1_g(nspecies)

    ! spack_utils
    INTEGER         ,   INTENT(IN)    :: nob
    INTEGER         ,   INTENT(IN)    :: maxblock_size
    INTEGER         ,   INTENT(IN)    :: block_end(:)
    INTEGER         ,   INTENT(IN)    :: indexk(:,:)
    INTEGER         ,   INTENT(IN)    :: indexi(:,:)
    INTEGER         ,   INTENT(IN)    :: indexj(:,:)
    INTEGER         ,   INTENT(IN)    :: kij_index(:,:)
    DOUBLE PRECISION,   INTENT(INOUT) :: last_accepted_dt(:)
!    INTEGER         ,   INTENT(IN)    :: block_end(nob)
!    INTEGER         ,   INTENT(IN)    :: indexk(maxblock_size,nob)
!    INTEGER         ,   INTENT(IN)    :: indexi(maxblock_size,nob)
!    INTEGER         ,   INTENT(IN)    :: indexj(maxblock_size,nob)
!    INTEGER         ,   INTENT(IN)    :: kij_index(maxblock_size,nob)
!    DOUBLE PRECISION,   INTENT(INOUT) :: last_accepted_dt(nob)

    ! mem_spack
    TYPE(spack_type)   , INTENT(INOUT) :: spack(1)
!    TYPE(spack_type_2d), INTENT(INOUT) :: spack_2d(maxblock_size,1)
    TYPE(spack_type_2d), INTENT(INOUT) :: spack_2d(:,:)
    DOUBLE PRECISION   , INTENT(IN)    :: atol(nspecies)
    DOUBLE PRECISION   , INTENT(IN)    :: rtol(nspecies)

    LOGICAL , INTENT(INOUT) :: get_non_zeros

    ! rconstants
    REAL, INTENT(IN) :: cp
    REAL, INTENT(IN) :: cpor
    REAL, INTENT(IN) :: p00

    ! extra
    INTEGER,     INTENT(IN)    :: na_extra3d
    TYPE(ext3d), INTENT(INOUT) :: extra3d(na_extra3d) ! Output JNO2

    ! FastJX
    DOUBLE PRECISION, INTENT(IN) :: jphoto(nr_photo,m1,m2,m3)

    ! uv_atten
    REAL, INTENT(IN) :: att(m2,m3)


!  INTEGER,PARAMETER :: block_size=2*32
!  INTEGER,PARAMETER :: block_size=4
!  INTEGER,PARAMETER :: block_size=1

   DOUBLE PRECISION,PARAMETER	 :: pmar=28.96D0
   DOUBLE PRECISION,PARAMETER	 :: Threshold1= 0.0D0
   DOUBLE PRECISION,PARAMETER	 :: Threshold2=-1.D-1

   DOUBLE PRECISION,PARAMETER	 :: Igamma = 0.5D0
   DOUBLE PRECISION,PARAMETER	 :: c43 = 4.D0/ 3.D0
   DOUBLE PRECISION,PARAMETER	 :: c83 = 8.D0/ 3.D0
   DOUBLE PRECISION,PARAMETER	 :: c56 = 5.D0/ 6.D0
   DOUBLE PRECISION,PARAMETER	 :: c16 = 1.D0/ 6.D0
   DOUBLE PRECISION,PARAMETER	 :: c112 = 1.D0/ 12.D0
   DOUBLE PRECISION dble_dtlt_i,fxc, Igamma_dtstep,dt_chem,dt_chem_i
   INTEGER :: i,ijk,n,j,k,ispc,Ji,Jj,k_,i_,j_,kij_,kij

  !integer,parameter ::  NOB_MEM=1   ! if 1 : alloc just one block (uses less memory)
  !                                  ! if 0 : alloc all blocks     (uses more memory)

  !- default number of allocatable blocks (re-use scratch arrays spack()% )
  integer,parameter :: inob=1

  !- parameters for dynamic timestep control
  integer all_accepted
  double precision :: dt_min,dt_max,dt_actual,dt_new
  double precision :: time_f,time_c
  double precision, parameter :: FacMin =0.2d0 ! lower bound on step decrease factor (default=0.2)
  double precision, parameter :: FacMax =6.0D0 ! upper bound on step increase factor (default=6)
  double precision, parameter :: FacRej =0.1D0 ! step decrease factor after multiple rejections
  double precision, parameter :: FacSafe=0.9D0 ! step by which the new step is slightly smaller
                                               ! than the predicted value  (default=0.9)
  double precision, parameter :: uround  = 1.d-15,   elo = 3.0D0, Roundoff = 1.d-8
  double precision Fac, Tol,err1,max_err1


  double precision, allocatable :: DLmat(:,:,:)
  integer, allocatable :: ipos(:)
  integer, allocatable :: jpos(:)

  integer :: NumberOfNonZeros, nz
  integer :: offset, offsetDLmat, offsetnz
  integer :: blocksize, sizeOfMatrix, maxnonzeros, blocknonzeros
  real :: start, finish
  real :: elapsed_time_solver,elapsed_time,elapsed_time_alloc,elapsed_time_dealloc,elapsed_time_copy


  INTEGER(KIND=8) :: matrix_Id
  INTEGER :: error
  INTEGER,PARAMETER :: complex_t=0
  INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:) :: element
  INTEGER(kind=8), EXTERNAL :: sfCreate_solve
  INTEGER(kind=8), EXTERNAL :: sfGetElement
  INTEGER, EXTERNAL :: sfFactor

  integer :: processor
 
  processor = mynum


  call write_input(time, processor, m1, m2, m3, nspecies, nr, nr_photo, maxnspecies, &
                      nspecies_chem_transported, nspecies_chem_no_transported, &
                      nob, maxblock_size, na_extra3d, dtlt, pp, pi0, theta , rv, dn0, rcp, &
                      cosz, weight, atol, rtol, jphoto, att, block_end, indexk, indexi, &
                      indexj, kij_index, transp_chem_index, no_transp_chem_index, &
                      extra3d, spack, spack_2d, last_accepted_dt, get_non_zeros, chem1_g, &
                      cp, cpor, n_dyn_chem, p00, PhotojMethod, split_method)



!!  if(.not. spack_alloc) then
!!  !- Determining number of blocks and masks, allocating spack type
!!   call allocIndex(block_size,N_DYN_CHEM)
!!   call alloc_spack(CHEMISTRY)
!!  end if

!  IF(.NOT. get_non_zeros) then
  !- Get always zero elements of Jacobian matrix to save
  !- computational time at solver:
  ! non-opt
  !get_non_zeros = .true.
  !call Prepare(nspecies)

   !  opt
   !!call get_always_zeros(nr_photo,nr,nspecies,get_non_zeros,spack(1)%rk(1,1:nr), &
   !                          spack(1)%jphoto(1,1:nr_photo),p00)

   call get_number_nonzeros(nr_photo,nr,nspecies,spack(1)%rk(1,1:nr), &
                             spack(1)%jphoto(1,1:nr_photo),p00,maxnonzeros)

  sizeOfMatrix = nspecies

  ALLOCATE(ipos(maxnonzeros));ipos=1
  ALLOCATE(jpos(maxnonzeros));jpos=1
!  ENDIF



  DO i=1,nob !- loop over all blocks

    !- copying structure from input to internal

    DO ijk=1,block_end(i) !index_g%block_end(i)

      kij_ =kij_index(ijk,i) !index_g%kij_index(ijk,i)!- kij brams tendency index
      k_=indexk(ijk,i) ! index_g%indexk(ijk,i) ! k brams index
      i_=indexi(ijk,i) ! index_g%indexi(ijk,i) ! i brams index
      j_=indexj(ijk,i) ! index_g%indexj(ijk,i) ! j brams index

      !- Exner function divided by Cp
      spack(inob)%press(ijk) = ( pp (k_,i_,j_) + pi0(k_,i_,j_) )/cp

      !- Air temperature (K)
       spack(inob)%temp(ijk)= theta(k_,i_,j_) * spack(inob)%press(ijk)

      !- transform from Exner function to pressure
       spack(inob)%press(ijk)= spack(inob)%press(ijk)**cpor*p00

      !- Water vapor
       spack(inob)%vapp(ijk) =  rv(k_,i_,j_)!&
      !                  * dn0(k_,i_,j_)*6.02D20/18.D0

      !- volmol e' usado para converter de ppbm para molecule/cm^3 (fator fxc incluido)
      ! spack(inob)%volmol(ijk)=(spack(inob)%press(ijk))/(8.314e0*spack(inob)%temp(ijk))
      ! fxc = (6.02e23*1e-15*spack(inob)%volmol(ijk))*pmar
      !
       spack(inob)%volmol(ijk)=(6.02D23*1D-15*pmar)*(spack(inob)%press(ijk))/(8.314D0*spack(inob)%temp(ijk))

      !- inverse of volmol to get back to ppbm at the end of routine
       spack(inob)%volmol_I(ijk)=1.0d0 / spack(inob)%volmol(ijk)

      !- get liquid water content
       spack(inob)%xlw(ijk) =  rcp(k_,i_,j_) * dn0(k_,i_,j_)*1.D-3


      !- convert from brams chem (ppbm) arrays to spack (molec/cm3)
      !- no transported species section
      DO ispc=1,nspecies_chem_no_transported

        !- map the species to NO transported ones
     	n=no_transp_chem_index(ispc)

        !- initialize no-transported species (don't need to convert, because these
        !- species are already saved using molecule/cm^3 units)
        spack(inob)%sc_p(ijk,n) = chem1_g(n)%sc_p(k_,i_,j_)
!       spack(inob)%sc_p_4(ijk,n) = chem1_g(n)%sc_p(k_,i_,j_)
      END DO

    END DO

    !- convert from brams chem (ppbm) arrays to spack (molec/cm3)
    !- transported species section
    IF(split_method == 'PARALLEL' .and. N_DYN_CHEM > 1) then
      DO ijk=1,block_end(i) !index_g%block_end(i)

	 kij_ =kij_index(ijk,i) !index_g%kij_index(ijk,i)!- kij brams tendency index
         k_=indexk(ijk,i) !index_g%indexk(ijk,i) ! k brams index
         i_=indexi(ijk,i) !index_g%indexi(ijk,i) ! i brams index
         j_=indexj(ijk,i) !index_g%indexj(ijk,i) ! j brams index

         DO ispc=1,nspecies_chem_transported

           !- map the species to transported ones
           n=transp_chem_index(ispc)

           !- conversion from ppbm to molecule/cm^3
	   !- get back the concentration at the begin of the chemistry timestep integration:

           spack(inob)%sc_p(ijk,n) =  (  chem1_g(n)%sc_p(k_,i_,j_) -                 &  ! updated mixing ratio
	                                 chem1_g(n)%sc_t_dyn(kij_)*N_DYN_CHEM*dtlt  )&  ! accumulated tendency
        			       * spack(inob)%volmol(ijk)/weight(n)
           !- for testing
           !spack(inob)%sc_p(ijk,n)   = max(0.,spack(inob)%sc_p(ijk,n))
           !spack(inob)%sc_p_4(ijk,n) =(chem1_g(n)%sc_p(k_,i_,j_) - beta*chem1_g(n)%sc_t_dyn(kij_)*N_DYN_CHEM*dtlt)& !dtlt \E9 o dinamico
           !				* spack(inob)%volmol(ijk)/weight(n)
         END DO
      END DO
    ELSE

      DO ijk=1,block_end(i) !index_g%block_end(i)

	 kij_ =kij_index(ijk,i) !index_g%kij_index(ijk,i)!- kij brams tendency index
         k_=indexk(ijk,i) !index_g%indexk(ijk,i) ! k brams index
         i_=indexi(ijk,i) !index_g%indexi(ijk,i) ! i brams index
         j_=indexj(ijk,i) !index_g%indexj(ijk,i) ! j brams index

         DO ispc=1,nspecies_chem_transported

           !- map the species to transported ones
           n=transp_chem_index(ispc)

           !- conversion from ppbm to molecule/cm^3
           spack(inob)%sc_p(ijk,n) = chem1_g(n)%sc_p(k_,i_,j_) * spack(inob)%volmol(ijk)/weight(n)
           spack(inob)%sc_p(ijk,n) = max(0.D0,spack(inob)%sc_p(ijk,n))

         END DO
      END DO
    ENDIF



    !- Photolysis section
    IF(trim(PhotojMethod) == 'FAST-JX' .or. trim(PhotojMethod) == 'FAST-TUV' ) THEN

     DO ijk=1,block_end(i) !index_g%block_end(i)

	k_=indexk(ijk,i) !index_g%indexk(ijk,i) ! k brams index
        i_=indexi(ijk,i) !index_g%indexi(ijk,i) ! i brams index
        j_=indexj(ijk,i) !index_g%indexj(ijk,i) ! j brams index

	DO n=1,nr_photo
          spack(inob)%jphoto(ijk,n)= jphoto(n,k_,i_,j_) !fast_JX_g%jphoto(n,k_,i_,j_)
        END DO

      ENDDO

     ELSEIF(trim(PhotojMethod) == 'LUT') then

      DO ijk=1,block_end(i) !index_g%block_end(i)

        i_=indexi(ijk,i) !index_g%indexi(ijk,i) ! i brams index
        j_=indexj(ijk,i) !index_g%indexj(ijk,i) ! j brams index

	!- UV  attenuation un function AOT
        spack(inob)%att(ijk)= att(i_,j_) !uv_atten_g%att(i_,j_)

	!- get zenital angle (for LUT PhotojMethod)
        spack(inob)%cosz(ijk)=  cosz(i_,j_)

      ENDDO
    ENDIF


!- compute kinetical and photochemical reactions (array: spack(inob)%rk)
    	CALL kinetic(nr_photo,spack(inob)%Jphoto  &
    			     ,spack(inob)%rk      &
  			     ,spack(inob)%temp    &
    			     ,spack(inob)%vapp    &
    			     ,spack(inob)%Press   &
    			     ,spack(inob)%cosz    &
    			     ,spack(inob)%att     &
			     ,1 &
                             ,block_end(i) & !index_g%block_end(i), &
                             ,maxblock_size,nr)


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     !tmp
!
!  !- srf: to use these arrays, be sure NA_EXTRAD3D RAMSIN is at least 8
!    DO ijk=1,block_end(i) !index_g%block_end(i)
!       k_=indexk(ijk,i) !index_g%indexk(ijk,i) ! k brams index
!       i_=indexi(ijk,i) !index_g%indexi(ijk,i) ! i brams index
 !      j_=indexj(ijk,i) !index_g%indexj(ijk,i) ! j brams index

 !  extra3d(7)%d3(k_,i_,j_) =spack(inob)%rk(ijk,1)
    !extra3d(8)%d3(k_,i_,j_) =spack(inob)%rk(ijk,2)
    !extra3d(9)%d3(k_,i_,j_) =spack(inob)%rk(ijk,3)
    !extra3d(10)%d3(k_,i_,j_)=spack(inob)%rk(ijk,4)
    !extra3d(11)%d3(k_,i_,j_)=spack(inob)%rk(ijk,5)
    !extra3d(12)%d3(k_,i_,j_)=spack(inob)%rk(ijk,6)
    !extra3d(13)%d3(k_,i_,j_)=spack(inob)%rk(ijk,7)
    !extra3d(14)%d3(k_,i_,j_)=spack(inob)%rk(ijk,8)
    !extra3d(15)%d3(k_,i_,j_)=spack(inob)%rk(ijk,9)
    !extra3d(16)%d3(k_,i_,j_)=spack(inob)%rk(ijk,10)
    !extra3d(17)%d3(k_,i_,j_)=spack(inob)%rk(ijk,11)
    !extra3d(18)%d3(k_,i_,j_)=spack(inob)%rk(ijk,12)
    !extra3d(19)%d3(k_,i_,j_)=spack(inob)%rk(ijk,13)
    !extra3d(20)%d3(k_,i_,j_)=spack(inob)%rk(ijk,14)
    !extra3d(21)%d3(k_,i_,j_)=spack(inob)%rk(ijk,16)
    !extra3d(22)%d3(k_,i_,j_)=spack(inob)%rk(ijk,17)
    !extra3d(23)%d3(k_,i_,j_)=spack(inob)%rk(ijk,18)
    !extra3d(24)%d3(k_,i_,j_)=spack(inob)%rk(ijk,19)
    !extra3d(25)%d3(k_,i_,j_)=spack(inob)%rk(ijk,20)
    !extra3d(26)%d3(k_,i_,j_)=spack(inob)%rk(ijk,21)
   ! extra3d(27)%d3(k_,i_,j_)=spack(inob)%rk(ijk,22)
  !  extra3d(28)%d3(k_,i_,j_)=spack(inob)%rk(ijk,23)

   ! ENDDO
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     !tmp






!- ROSENBROCK METHOD ----------------------------------------------------------------------
!- kml: Inicio do equivalente a roschem
!- srf: extending to RODAS 3 (4 stages, order 3)

  dt_chem   = last_accepted_dt(i) !index_g%last_accepted_dt(i) ! dble(dtlt)
  dt_min   = max(1.0D0,1.D-2*dble(dtlt*N_DYN_CHEM))
  dt_max   = dble(dtlt*N_DYN_CHEM)
  dt_new   = 0.0d0

  time_c   = 0.0d0
  time_f   = dble(dtlt*N_DYN_CHEM)

  run_until_integr_ends:   DO WHILE (time_c + Roundoff < time_f)

!-   Compute the Jacobian (DLRDC).
     CALL jacdchemdc (spack(inob)%sc_p  &
!     CALL jacdchemdc (spack(inob)%sc_p_4  &
	             ,spack(inob)%rk	&
      	             ,spack(inob)%DLdrdc& ! Jacobian matrix
		     ,nspecies,1 &
                     ,block_end(i) & !index_g%block_end(i), &
                     ,maxblock_size,nr)

!-   Compute chemical net production terms (actually, P-L) at initial time (array spack(inob)%DLr)
     CALL fexchem (spack(inob)%sc_p &
!     CALL fexchem (spack(inob)%sc_p_4 &
        	  ,spack(inob)%rk   &
        	  ,spack(inob)%DLr  & !production term
        	  ,nspecies,1 &
                  ,block_end(i) & !index_g%block_end(i), &
                  ,maxblock_size,nr)

     DO Ji=1,nspecies
       DO ijk=1,block_end(i) !index_g%block_end(i)
	   spack_2d(ijk,inob)%DLb1(Ji) = spack(inob)%DLr(ijk,Ji)
       ENDDO
     ENDDO

!-   compute matrix (1/(Igamma*dt)  - Jacobian)  where Jacobian = DLdrdc
!-   fill DLMAT with non-diagonal (fixed in the timestep) Jacobian
     DO Jj=1,nspecies
        DO Ji=1,nspecies
          DO ijk=1,block_end(i) !index_g%block_end(i)

           spack_2d(ijk,inob)%DLmat(Ji,Jj) = - spack(inob)%DLdrdc(ijk,Ji,Jj)

          ENDDO
        ENDDO
     ENDDO

     UntilAccepted: DO

!- fill DLMAT with diagonal Jacobian (which changes in time, because of dt_chem)
     Igamma_dtstep = 1.0D0/(Igamma * dt_chem)

     DO Jj=1,nspecies
        DO ijk=1,block_end(i) !index_g%block_end(i)
             spack_2d(ijk,inob)%DLmat(Jj,Jj) = Igamma_dtstep - spack(inob)%DLdrdc(ijk,Jj,Jj)
        ENDDO
     ENDDO

     if ( (i .eq. 1) ) then
        blocknonzeros = 0
        DO Ji=1,nspecies
          DO Jj=1,nspecies
              !IF(spack_2d(1,inob)%DLmat(Ji,Jj)<=0 .AND. spack_2d(1,inob)%DLmat(Ji,Jj)>=(-1)*0) CYCLE
              if ( (spack_2d(1,inob)%DLmat(Ji,Jj) .ne. 0.0d+00)  ) then
                 blocknonzeros = blocknonzeros + 1
                 ipos(blocknonzeros)=Ji
                 jpos(blocknonzeros)=Jj
              endif
          ENDDO
        ENDDO
        NumberOfNonZeros = blocknonzeros

        !create matrix
        matrix_Id = sfCreate_solve(sizeOfMatrix,complex_t,error)

        IF(ALLOCATED(element)) DEALLOCATE (element)
        ALLOCATE(element(NumberOfNonZeros))

        DO nz=1,NumberOfNonZeros
           Ji = ipos(nz)
           Jj = jpos(nz)
           element(nz)=sfGetElement(matrix_Id,Ji,Jj)
           !element(nz)=sfGetElement(matrix_Id,INT(DLmat(1,nz,1)),INT(DLmat(1,nz,2)))
        END DO

        CALL sfZero(matrix_Id)
        DO nz=1,NumberOfNonZeros
           Ji = ipos(nz)
           Jj = jpos(nz)
           CALL sfAdd1Real(element(nz),spack_2d(1,inob)%DLmat(Ji,Jj))
           !CALL sfAdd1Real(element(nz),DLmat(1,nz,3))
        END DO

        error = sfFactor(matrix_Id)
     ENDIF

   !@LNCC: begin points loop
   DO ijk=1,block_end(i) !index_g%block_end(i)
      CALL sfZero(matrix_Id)
      DO nz=1,NumberOfNonZeros
         Ji = ipos(nz)
         Jj = jpos(nz)
         CALL sfAdd1Real(element(nz),spack_2d(ijk,inob)%DLmat(Ji,Jj))
         !CALL sfAdd1Real(element(nz),DLmat(ijk,nz,3))
      END DO
      error = sfFactor(matrix_Id)

!---------------------------------------------------------------------------------------------------------------------
!    1- First step
!-   Compute DLk1 by Solving (1/(Igamma*dt) - DLRDC) DLk1=DLR

     !- solver sparse
     call sfSolve(matrix_Id, spack_2d(ijk,inob)%DLb1, spack_2d(ijk,inob)%DLk1)
!
!---------------------------------------------------------------------------------------------------------------------
!    2- Second step
!    compute   K2 by solving (1/0.5 h JAC)K2 = 4/h * K1 +  F(Yn)
!    compute DLK2 by solving (1/Igama*dt - DLRDC)DLK2 = 4/h * DLK1 +  DLb1
     dt_chem_i = 1.0d0/dt_chem
     DO Ji=1,nspecies
        spack_2d(ijk,inob)%DLb2(Ji) = (4.0D0*dt_chem_i)*spack_2d(ijk,inob)%DLk1(Ji)  + &
                                                        spack_2d(ijk,inob)%DLb1(Ji)
     ENDDO
     call sfSolve(matrix_Id, spack_2d(ijk,inob)%DLb2, spack_2d(ijk,inob)%DLk2)

!---------------------------------------------------------------------------------------------------------------------
!    3- Third step
!    a) update concentrations

     !dt_chem_i = 1.0d0/dt_chem
     DO Ji=1,nspecies
    	    spack(inob)%sc_p_new(ijk,Ji) = spack(inob)%sc_p(ijk,Ji) + 2.0D0* spack_2d(ijk,inob)%DLk1(Ji)

	    IF (spack(inob)%sc_p_new(ijk,Ji) .LT. Threshold1) THEN
	        spack(inob)%sc_p_new(ijk,Ji) = Threshold1
	        spack_2d(ijk,inob)%DLk1(Ji)  = 0.5D0*(spack(inob)%sc_p_new(ijk,Ji) - spack(inob)%sc_p(ijk,Ji))
     	    ENDIF
     ENDDO
!
!    b) update the net production term (DLr= P-L = F(Y3)) at this stage with the first-order
!       approximation with the new concentration
!
     CALL fexchem (spack(inob)%sc_p_new  &
     		  ,spack(inob)%rk        &
	  	  ,spack(inob)%DLr       &
     		  ,nspecies,ijk &
                  ,ijk & !index_g%block_end(i), &
                  ,maxblock_size,nr)

!    c) compute   K3 by solving (1 /(0.5 h) - JAC)K3 =   F(Y3) + 0.5 (K1-K2)
     DO Ji=1,nspecies
	   spack_2d(ijk,inob)%DLb3(Ji) = spack(inob)%DLr(ijk,Ji) + dt_chem_i* &
	                               ( spack_2d(ijk,inob)%DLk1(Ji) - spack_2d(ijk,inob)%DLk2(Ji))
     ENDDO
     call sfSolve(matrix_Id, spack_2d(ijk,inob)%DLb3, spack_2d(ijk,inob)%DLk3)

!---------------------------------------------------------------------------------------------------------------------
!    4- Fourth step
!    a) update concentrations
!       Y4 = Yn + 2 * k1 +  K3
     dt_chem_i = 1.0d0/dt_chem
     DO Ji=1,nspecies
    	    spack(inob)%sc_p_new(ijk,Ji) = spack(inob)%sc_p(ijk,Ji)          +  &
	                                 2.0D0 * spack_2d(ijk,inob)%DLk1(Ji) + &
					 	 spack_2d(ijk,inob)%DLk3(Ji)

	    IF (spack(inob)%sc_p_new(ijk,Ji) .LT. Threshold1) THEN
	        spack(inob)%sc_p_new(ijk,Ji) = Threshold1
	        spack_2d(ijk,inob)%DLk3(Ji) = (spack(inob)%sc_p_new(ijk,Ji) - spack(inob)%sc_p(ijk,Ji)) &
	   				      -2.0D0 * spack_2d(ijk,inob)%DLk1(Ji)
     	    ENDIF
     ENDDO
!    b) update the net production term (DLr= P-L = F(Y4) ) at this stage with the 3rd-order
!       approximation with the new concentration
!
     CALL fexchem (spack(inob)%sc_p_new  & ! Y4
     		  ,spack(inob)%rk        &
	  	  ,spack(inob)%DLr       & ! F(Y4)
     		  ,nspecies,ijk &
                  ,ijk & !index_g%block_end(i)
                  ,maxblock_size,nr)

!    c) compute   K4 by solving (1/(0.5 h)- JAC)K4 =    F(Y3) + K1/h -K2/h -8/3 K3/h
     DO Ji=1,nspecies
	   spack_2d(ijk,inob)%DLb4(Ji) =  spack(inob)%DLr (ijk,Ji) + dt_chem_i* & ! F(Y4)
	                                (        spack_2d(ijk,inob)%DLk1(Ji)    &
				         -       spack_2d(ijk,inob)%DLk2(Ji)    &
				         - c83 * spack_2d(ijk,inob)%DLk3(Ji)    )
     ENDDO
     call sfSolve(matrix_Id, spack_2d(ijk,inob)%DLb4, spack_2d(ijk,inob)%DLk4)

   ENDDO
   !@LNCC: end points loop


!---------------------------------------------------------------------------------------------------------------------
!   - the solution
     dt_chem_i = 1.0d0/dt_chem
     DO Ji=1,nspecies
         DO ijk=1,block_end(i) !index_g%block_end(i)

    	    spack(inob)%sc_p_new(ijk,Ji) = spack(inob)%sc_p(ijk,Ji) +  &
	                                    2.0D0 * spack_2d(ijk,inob)%DLk1(Ji)  &
					  +	    spack_2d(ijk,inob)%DLk3(Ji)  &
					  +         spack_2d(ijk,inob)%DLk4(Ji)
            spack(inob)%sc_p_new(ijk,Ji) = max (spack(inob)%sc_p_new(ijk,Ji),Threshold1)

	! IF (spack(inob)%sc_p_new(ijk,Ji) .LT. Threshold1) THEN
	!    spack(inob)%sc_p_new(ijk,Ji) = Threshold1
	!    spack_2d(ijk,inob)%DLk3(Ji) = (spack(inob)%sc_p_new(ijk,Ji) - spack(inob)%sc_p(ijk,Ji)) * dt_chem_i
     	! ENDIF
        !- keep non-negative values for concentrations
	!   spack(inob)%sc_p_new(ijk,Ji) = max(spack(inob)%sc_p_new(ijk,Ji),Threshold1)
        !- update DLK2, in case sc_p_new is changed by the statement above
	!   spack_2d(ijk,inob)%DLk2(Ji) = 2.0D0*( (spack(inob)%sc_p_new(ijk,Ji)-spack(inob)%sc_p(ijk,Ji)) &
	!                                         * dt_chem_i - 1.5D0 * spack_2d(ijk,inob)%DLk1(Ji) )
       ENDDO
     ENDDO


!
!-  Compute the error estimation : spack(inob)%err
     DO ijk=1,block_end(i) !index_g%block_end(i)

       spack(inob)%err(ijk) = 0.0d0
       DO Ji=1,nspecies

	 Tol = ATol(Ji) + RTol(Ji)*DMAX1( DABS(spack(inob)%sc_p(ijk,Ji)),DABS(spack(inob)%sc_p_new(ijk,Ji)) )

	 err1= spack_2d(ijk,inob)%DLk4(Ji)

         spack(inob)%err(ijk) = spack(inob)%err(ijk) + (err1/Tol)**2.0D0

       ENDDO

       spack(inob)%err(ijk) = DMAX1( uround, DSQRT( spack(inob)%err(ijk)/nspecies ) )
     ENDDO

     all_accepted = 1
     DO ijk=1,block_end(i) !index_g%block_end(i)
        if( spack(inob)%err(ijk) - Roundoff > 1.0D0) then
	   all_accepted = 0 ; exit
	endif
     ENDDO

     !- find the maximum error occurred
     max_err1 = maxval( spack(inob)%err(1:block_end(i))) !index_g%block_end(i) ) )

     !- use it to determine the new time step for all block elements
     !- new step size is bounded by FacMin <= Hnew/H <= FacMax
     Fac  = MIN(FacMax,MAX(FacMin,FacSafe/max_err1**(1.0D0/elo)))

     !- possible new timestep
     dt_new = dt_chem * Fac
     dt_new = max(dt_min,min(dt_max,dt_new))

     !- to reset the timestep resizing in function of the error estimation, use the statements below:
     ! all_accepted = 1; dt_new=dt_max

     if( all_accepted == 0 .and. dt_new > dt_min) then  ! current solution is not accepted

       !- resize the timestep and try again
       dt_chem = dt_new

     else    ! current solution is     accepted

      !- go ahead, updating spack(inob)%sc_p with the solution (spack(inob)%sc_p_new)
      !- next time
      time_c = time_c + dt_chem

      !- next timestep (dt_new but limited by the time_f-time_c, the rest of time integration interval)
      dt_chem = min(dt_new, time_f-time_c)

      !- save the accepted timestep for the next integration interval
      if(time_c < time_f)     last_accepted_dt(i) = dt_new !index_g%last_accepted_dt(i) =    dt_new

      !- pointer (does not work yet)
      ! spack(inob)%sc_p=>spack(inob)%sc_p_new     ! POINTER
      !- copy
        DO Ji=1,nspecies
           DO ijk=1,block_end(i) !index_g%block_end(i)
               spack(inob)%sc_p(ijk,Ji) = spack(inob)%sc_p_new(ijk,Ji)
              ! spack(inob)%sc_p_4(ijk,Ji) = spack(inob)%sc_p_new(ijk,Ji)
           ENDDO
        ENDDO

	EXIT UntilAccepted

     endif

   END DO UntilAccepted


  ENDDO run_until_integr_ends ! time-spliting


  !--------------------------------------------------------------------------------------------
  !- Restoring species tendencies OR updated mixing ratios from internal to brams structure

  !- transported species section

    if(split_method == 'PARALLEL' .and.  N_DYN_CHEM > 1) then

       DO ijk=1,block_end(i) !index_g%block_end(i)

         kij_ =kij_index(ijk,i) !index_g%kij_index(ijk,i)
         k_   =indexk(ijk,i) !index_g%indexk(ijk,i)
         i_   =indexi(ijk,i) !index_g%indexi(ijk,i)
         j_   =indexj(ijk,i) !index_g%indexj(ijk,i)

         DO ispc=1,nspecies_chem_transported

            !- map the species to transported ones
            n=transp_chem_index(ispc)

            !- include the chemical tendency at total tendency (convert to unit: ppbm/s)
            chem1_g(n)%sc_p(k_,i_,j_ ) =  chem1_g(n)%sc_t_dyn(kij_)*N_DYN_CHEM*dtlt +&
        				       spack(inob)%sc_p(ijk,n)*weight(n)*spack(inob)%volmol_I(ijk)

            chem1_g(n)%sc_p(k_,i_,j_ ) = max(0.,chem1_g(n)%sc_p(k_,i_,j_ ))

        END DO
      ENDDO

    ELSEIF(split_method == 'PARALLEL' .and.  N_DYN_CHEM == 1) then

       dble_dtlt_i=1.0d0 / dble( dtlt)
       DO ijk=1,block_end(i) !index_g%block_end(i)

         kij_ =kij_index(ijk,i) !index_g%kij_index(ijk,i)
         k_   =indexk(ijk,i) !index_g%indexk(ijk,i)
         i_   =indexi(ijk,i) !index_g%indexi(ijk,i)
         j_   =indexj(ijk,i) !index_g%indexj(ijk,i)

         DO ispc=1,nspecies_chem_transported

            !- map the species to transported ones
            n=transp_chem_index(ispc)

            !- include the chemical tendency at total tendency (convert to unit: ppbm/s)
!           chem1_g(n)%sc_t(kij_) =		    +  &! use this for update only chemistry (No dyn/emissions)
            chem1_g(n)%sc_t(kij_) = chem1_g(n)%sc_t(kij_)			    +  &! previous tendency
      			     ( spack(inob)%sc_p(ijk,n)*weight(n)*spack(inob)%volmol_I(ijk)  -  &! new mixing ratio
       			       chem1_g(n)%sc_p(k_,i_,j_ ) 		  )	       &! old mixing ratio
       			     * dble_dtlt_i							! inverse of timestep
         END DO
       END DO

     ELSE

       DO ijk=1,block_end(i) !index_g%block_end(i)

         kij_ =kij_index(ijk,i) !index_g%kij_index(ijk,i)
         k_   =indexk(ijk,i) !index_g%indexk(ijk,i)
         i_   =indexi(ijk,i) !index_g%indexi(ijk,i)
         j_   =indexj(ijk,i) !index_g%indexj(ijk,i)

	 DO ispc=1,nspecies_chem_transported

            !- map the species to transported ones
            n=transp_chem_index(ispc)

  	    chem1_g(n)%sc_p(k_,i_,j_ ) = spack(inob)%sc_p(ijk,n)*weight(n)*spack(inob)%volmol_I(ijk)
	    chem1_g(n)%sc_p(k_,i_,j_ ) = max(0.,chem1_g(n)%sc_p(k_,i_,j_ ))
        END DO
       END DO

     ENDIF


  !- no transported species section
     DO ijk=1,block_end(i) !index_g%block_end(i)

         kij_ =kij_index(ijk,i) !index_g%kij_index(ijk,i)
         k_   =indexk(ijk,i) !index_g%indexk(ijk,i)
         i_   =indexi(ijk,i) !index_g%indexi(ijk,i)
         j_   =indexj(ijk,i) !index_g%indexj(ijk,i)
         DO ispc=1,nspecies_chem_no_transported

	  !- map the species to no transported ones
          n=no_transp_chem_index(ispc)

	  !- save no-transported species (keep current unit : molec/cm3)
          chem1_g(n)%sc_p(k_,i_,j_ )  = max(0., real ( spack(inob)%sc_p(ijk,n) ))
         END DO

     END DO

 END DO ! enddo loop over all blocks

call write_output(time, m1, m2, m3, processor, chem1_g, nspecies)

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!DO n=1,nspecies
!print*,'max1=',spc_name(n),maxval(chem1_g(n)%sc_p(:,:,:)),maxloc(chem1_g(n)%sc_p(:,:,:))
!print*,'max2=',n,spc_name(n),maxval(spack(:)%sc_p(:,n)),maxloc(spack(:)%sc_p(:,n))
!call flush(6)
!END DO
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

END SUBROUTINE chem_rodas3_dyndt


subroutine write_input(time, processor, m1, m2, m3, nspecies, nr, nr_photo, maxnspecies, &
                      nspecies_chem_transported, nspecies_chem_no_transported, &
                      nob, maxblock_size, na_extra3d, dtlt, pp, pi0, theta , rv, dn0, rcp, &
                      cosz, weight, atol, rtol, jphoto, att, block_end, indexk, indexi, &
                      indexj, kij_index, transp_chem_index, no_transp_chem_index, &
                      extra3d, spack, spack_2d, last_accepted_dt, get_non_zeros, chem1_g, &
                      cp, cpor, n_dyn_chem, p00, PhotojMethod, split_method)
    
        integer, intent(in) :: processor 
        integer, intent(in) :: m1, m2, m3
        integer, intent(in) :: nspecies, nr, nr_photo, maxnspecies
        integer, intent(in) :: nspecies_chem_transported, nspecies_chem_no_transported
        integer, intent(in) :: nob, maxblock_size, na_extra3d, n_dyn_chem

        integer, intent(in) :: block_end(nob)
        integer, intent(in) :: indexk(maxblock_size,nob)
        integer, intent(in) :: indexi(maxblock_size,nob)
        integer, intent(in) :: indexj(maxblock_size,nob)
        integer, intent(in) :: kij_index(maxblock_size,nob)
        integer, intent(in) :: transp_chem_index(maxnspecies)
        integer, intent(in) :: no_transp_chem_index(maxnspecies)

        real, intent(in)    :: time
        real, intent(in) :: dtlt
        real, intent(in) :: cp, cpor, p00
        real, intent(in) :: pp    (m1,m2,m3)
        real, intent(in) :: pi0   (m1,m2,m3)
        real, intent(in) :: theta (m1,m2,m3)
        real, intent(in) :: rv    (m1,m2,m3)
        real, intent(in) :: dn0   (m1,m2,m3)
        real, intent(in) :: rcp   (m1,m2,m3)
        real, intent(in) :: cosz  (m2,m3)
        real, intent(in) :: weight(nspecies)
        real, intent(in) :: att(m2, m3)

        double precision, intent(in) :: last_accepted_dt(:)
        double precision, intent(in) :: atol(nspecies)
        double precision, intent(in) :: rtol(nspecies)
        double precision, intent(in) :: jphoto(nr_photo,m1,m2,m3)

        logical , intent(inout) :: get_non_zeros
        character(len=10), intent(in) :: PhotojMethod
        character(len=20), intent(in) :: split_method

        type (chem1_vars), intent(inout) :: chem1_g(nspecies)
        type(ext3d), intent(inout) :: extra3d(na_extra3d)
        type(spack_type)   , intent(inout) :: spack(1)
        type(spack_type_2d), intent(inout) :: spack_2d(:,:)

        character(len=64) :: fileName
        integer :: iunit, ispc, i, ii, ijk

        write(filename, fmt='("chem_ia/chem_inp_",I6.6,"-",I4.4,".bin")') int(time),processor
        print *,"File to write: ",filename
        open (newunit=iunit, file=filename, status="replace", form="unformatted", action="write")    
        write(iunit)  m1 !int scalar
        write(iunit)  m2 !int scalar
        write(iunit)  m3 !int scalar

        write(iunit)  nspecies   !int scalar
        write(iunit)  nr_photo   !int scalar
        write(iunit)  nr         !int scalar
        write(iunit)  maxnspecies !int scalar
        write(iunit)  nspecies_chem_transported !int scalar
        write(iunit)  nspecies_chem_no_transported !int scalar
        write(iunit)  nob        !int scalar
        write(iunit)  maxblock_size !int scalar
        write(iunit)  na_extra3d !int scalar
        write(iunit)  n_dyn_chem !int scalar

        write(iunit)  block_end ! int 1D (nob)
        write(iunit)  indexk ! int 2D (maxblock_size,nob)
        write(iunit)  indexi ! int 2D (maxblock_size,nob)
        write(iunit)  indexj ! int 2D (maxblock_size,nob)
        write(iunit)  kij_index ! int 2D (maxblock_size,nob)
        write(iunit)  transp_chem_index ! int 1D (maxnspecies)
        write(iunit)  no_transp_chem_index ! int 1D (maxnspecies)

        write(iunit)  dtlt   !real scalar
        write(iunit)  cp     !real scalar
        write(iunit)  cpor   !real scalar
        write(iunit)  p00    !real scalar
        write(iunit)  pp     !real 3D (m1,m2,m3) 
        write(iunit)  pi0    !real 3D (m1,m2,m3)
        write(iunit)  theta  !real 3D (m1,m2,m3)
        write(iunit)  rv     !real 3D (m1,m2,m3)
        write(iunit)  dn0    !real 3D (m1,m2,m3)
        write(iunit)  rcp    !real 3D (m1,m2,m3)
        write(iunit)  cosz   !real 3D (m1,m2,m3)
        write(iunit)  weight !real 1D (nspecies)
        write(iunit)  att    !real 2D (m2,m3)
        
        write(iunit)  last_accepted_dt !double 1D (:)
        write(iunit)  atol !double 1D (nspecies)
        write(iunit)  rtol !double 1D (nspecies)
        write(iunit)  jphoto !double 4D (nr_photo,m1,m2,m3)

        write(iunit)  get_non_zeros !scalar logical
        
        write(iunit)  PhotojMethod ! char 10 len
        write(iunit)  split_method ! char 20 len

        do ispc=1,nspecies
            write(iunit) size(chem1_g(ispc)%sc_t) !real
            write(iunit) size(chem1_g(ispc)%sc_t_dyn) !real
        end do

        write(iunit) grid_g(1)%glat !real 2d
        write(iunit) grid_g(1)%glon !real 2d

        write(iunit) spc_name ! char 1D(2D?) len 8 (nspecies)

        do ispc=1,nspecies
             write(iunit)  chem1_g(ispc)%sc_p     ! real 3D (n1,n2,n3)
             write(iunit)  chem1_g(ispc)%sc_t     ! real 1D (nmxp(ng)*nmyp(ng)*nmzp(ng))
             write(iunit)  chem1_g(ispc)%sc_t_dyn ! real 1D (n1*n2*n3)
        end do
        !Double precision reals
        ! spack size is (nob), nob=1 is (supposed to) the default
        write(iunit)  spack(1)%DLdrdc   !3D (1:maxblock_size,nspecies,nspecies))
        write(iunit)  spack(1)%jphoto   !2D (1:maxblock_size,nr_photo))
        write(iunit)  spack(1)%rk       !2D (1:maxblock_size,nr))
        write(iunit)  spack(1)%w        !2D (1:maxblock_size,nr))
        write(iunit)  spack(1)%sc_p     !2D (1:maxblock_size,nspecies))
        write(iunit)  spack(1)%sc_p_new !2D (1:maxblock_size,nspecies))
        write(iunit)  spack(1)%DLr	 	  !2D (1:maxblock_size,nspecies))

        write(iunit)  spack(1)%temp	    !1D (1:maxblock_size))
        write(iunit)  spack(1)%press    !1D (1:maxblock_size))
        write(iunit)  spack(1)%cosz	    !1D (1:maxblock_size))
        write(iunit)  spack(1)%att	    !1D (1:maxblock_size)) 
        write(iunit)  spack(1)%vapp     !1D (1:maxblock_size))
        write(iunit)  spack(1)%volmol   !1D (1:maxblock_size))
        write(iunit)  spack(1)%volmol_i !1D (1:maxblock_size))
        write(iunit)  spack(1)%xlw      !1D (1:maxblock_size))
        write(iunit)  spack(1)%err      !1D (1:maxblock_size))

        ! spack_2d size is (maxblock_size,nob)
        do i = 1, 1
            do ijk = 1, maxblock_size
              ! double precision reals
                write(iunit)  spack_2d(ijk,i)%DLmat !2D (nspecies, nspecies))
                write(iunit)  spack_2d(ijk,i)%DLb1	!1D (nspecies)
                write(iunit)  spack_2d(ijk,i)%DLb2	!1D (nspecies)
                write(iunit)  spack_2d(ijk,i)%DLk1  !1D (nspecies)
                write(iunit)  spack_2d(ijk,i)%DLk2  !1D (nspecies)
                write(iunit)  spack_2d(ijk,i)%DLb3	!1D (nspecies)
                write(iunit)  spack_2d(ijk,i)%DLb4	!1D (nspecies)
                write(iunit)  spack_2d(ijk,i)%DLk3  !1D (nspecies)
                write(iunit)  spack_2d(ijk,i)%DLk4  !1D (nspecies)
            end do
        end do 

        do i=1,na_extra3d
            write(iunit)  extra3d(1)%d3 !real 3D (m1,m2,m3)

        end do

      close(iunit)

end subroutine write_input


subroutine write_output(time, m1, m2, m3, processor, chem1_g, nspecies)
    
        integer, intent(in) :: processor, m1,m2,m3
        integer, intent(in) :: nspecies !, nob, maxblock_size, na_extra3d
        !integer, intent(in) :: block_end(nob)

        real, intent(in)    :: time

        !double precision, intent(in) :: last_accepted_dt(:)

        !logical , intent(inout) :: get_non_zeros

        type (chem1_vars), intent(inout) :: chem1_g(nspecies)
        !type(ext3d), intent(inout) :: extra3d(na_extra3d)
        !type(spack_type)   , intent(inout) :: spack(1)
        !type(spack_type_2d), intent(inout) :: spack_2d(:,:)

        integer :: iunit, ispc, i, ijk
        character(len=64) :: fileName

        if (.false.) then
        write (filename, fmt='("chem_ia/chem_out_",I6.6,"-",I4.4,".bin")') int(time),processor
        open (newunit=iunit, file=filename, status="replace", form="unformatted", action="write")
        
        write(iunit) m1,m2,m3
        write(iunit) grid_g(1)%glat
        write(iunit) grid_g(1)%glon
        
        do ispc=1,nspecies
            write(iunit)  chem1_g(ispc)%sc_p
            write(iunit)  chem1_g(ispc)%sc_t
            write(iunit)  chem1_g(ispc)%sc_t_dyn
        end do

        write(iunit) spc_name
!        write(iunit)  spack(1)%DLdrdc
!        write(iunit)  spack(1)%sc_p_new
!        !write(iunit,*)  spack(1)%sc_p_4
!        write(iunit)  spack(1)%DLr	 
!        !write(iunit,*)  spack(1)%DLr3	  
!        write(iunit)  spack(1)%jphoto 
!        write(iunit)  spack(1)%rk    
!        write(iunit)  spack(1)%w    
!        write(iunit)  spack(1)%sc_p
!        write(iunit)  spack(1)%temp	
!        write(iunit)  spack(1)%press 
!        write(iunit)  spack(1)%cosz	  
!        write(iunit)  spack(1)%att	  
!        write(iunit)  spack(1)%vapp 
!        write(iunit)  spack(1)%volmol 
!        write(iunit)  spack(1)%volmol_i 
!        write(iunit)  spack(1)%xlw 
!        write(iunit)  spack(1)%err    
!        do i = 1, nob
!          do ijk = 1, block_end(i)
!            write(iunit)  spack_2d(ijk,inob)%DLmat 
!            write(iunit)  spack_2d(ijk,inob)%DLb1	
!            write(iunit)  spack_2d(ijk,inob)%DLb2	
!            write(iunit)  spack_2d(ijk,inob)%DLb3	
!            write(iunit)  spack_2d(ijk,inob)%DLb4	
!            write(iunit)  spack_2d(ijk,inob)%DLk1  
!            write(iunit)  spack_2d(ijk,inob)%DLk2  
!            write(iunit)  spack_2d(ijk,inob)%DLk3  
!            write(iunit)  spack_2d(ijk,inob)%DLk4 
!          end do
!        end do 
!        write(iunit)  last_accepted_dt !(:) !IO !(nob)
!        write(iunit)  get_non_zeros !IO
!        do i=1,na_extra3d
!          write(iunit)  extra3d(i)%d3 !(na_extra3d) ! Output JNO2 !IO
!        end do

        close(iunit)
      end if
        call create_netcdf_file(time, processor, m1, m2, m3, nspecies, chem1_g)

end subroutine write_output


subroutine create_netcdf_file(time, processor, m1, m2, m3, nspecies, chem1_g)
  USE netcdf
  IMPLICIT NONE

  integer, intent(in) :: m1, m2, m3, nspecies, processor
  real, intent(in)    :: time
  type(chem1_vars), intent(inout) :: chem1_g(nspecies)

  integer :: ncid, status
  integer :: dim_m1, dim_m2, dim_m3, dim_nspecies, dim_strlen
  integer :: var_m1, var_m2, var_m3, var_nspecies
  integer :: var_glat, var_glon, var_spc_name
  integer :: var_scp, var_sct, var_sctdyn
  integer :: i, n, str_len
  character(len=64) :: filename

  ! Create filename
  write(filename, fmt='("chem_ia/chem_out_",I6.6,"-",I4.4,".nc4")') int(time), processor

  status = NF90_CREATE(TRIM(filename), OR(NF90_NETCDF4, NF90_CLOBBER), ncid)
  if (status /= NF90_NOERR) CALL handle_error(status, "Create File")

  ! Define dimensions
  status = NF90_DEF_DIM(ncid, 'level', m1, dim_m1)
  if (status /= NF90_NOERR) CALL handle_error(status, "def level")
  status = NF90_DEF_DIM(ncid, 'latitude', m2, dim_m2)
  if (status /= NF90_NOERR) CALL handle_error(status, "def lat")
  status = NF90_DEF_DIM(ncid, 'longitude', m3, dim_m3)
  if (status /= NF90_NOERR) CALL handle_error(status, "def lon")
  status = NF90_DEF_DIM(ncid, 'nspecies', nspecies, dim_nspecies)
  if (status /= NF90_NOERR) CALL handle_error(status, "def nspecies")
  status = NF90_DEF_DIM(ncid, 'strlen', 64, dim_strlen)
  if (status /= NF90_NOERR) CALL handle_error(status, "def strlen")

  ! Define Coordinates
  status = NF90_DEF_VAR(ncid, 'latitude', NF90_FLOAT, dim_m2, var_glat)
  if (status /= NF90_NOERR) CALL handle_error(status, "")
  status = NF90_PUT_ATT(ncid, var_glat, 'units', 'degrees_north')
  if (status /= NF90_NOERR) CALL handle_error(status, "")

  status = NF90_DEF_VAR(ncid, 'longitude', NF90_FLOAT,  dim_m3, var_glon)
  if (status /= NF90_NOERR) CALL handle_error(status, "")
  status = NF90_PUT_ATT(ncid, var_glon, 'units', 'degrees_east')
  if (status /= NF90_NOERR) CALL handle_error(status, "")

  ! Define Species Names Variable
  status = NF90_DEF_VAR(ncid, 'species_names', NF90_CHAR, (/dim_strlen, dim_nspecies/), var_spc_name)
  if (status /= NF90_NOERR) CALL handle_error(status, "")

  ! Define Main Variables WITH DEFLATE AND SHUFFLE
  status = NF90_DEF_VAR(ncid, 'concentration', NF90_FLOAT, (/dim_m1, dim_m2, dim_m3, dim_nspecies/), var_scp)
  if (status /= NF90_NOERR) CALL handle_error(status, "")
  status = NF90_DEF_VAR_DEFLATE(ncid, var_scp, shuffle=1, deflate=1, deflate_level=1)
  if (status /= NF90_NOERR) CALL handle_error(status, "")

  status = NF90_DEF_VAR(ncid, 'tendency', NF90_FLOAT, (/dim_m1, dim_m2, dim_m3, dim_nspecies/), var_sct)
  if (status /= NF90_NOERR) CALL handle_error(status, "")
  status = NF90_DEF_VAR_DEFLATE(ncid, var_sct, shuffle=1, deflate=1, deflate_level=1)
  if (status /= NF90_NOERR) CALL handle_error(status, "")

  status = NF90_DEF_VAR(ncid, 'dynamics', NF90_FLOAT, (/dim_m1, dim_m2, dim_m3, dim_nspecies/), var_sctdyn)
  if (status /= NF90_NOERR) CALL handle_error(status, "")
  status = NF90_DEF_VAR_DEFLATE(ncid, var_sctdyn, shuffle=1, deflate=1, deflate_level=1)
  if (status /= NF90_NOERR) CALL handle_error(status, "")

  ! End Definitions
  status = NF90_ENDDEF(ncid)
  if (status /= NF90_NOERR) CALL handle_error(status, "End Def")

  ! Write Grid Data
  status = NF90_PUT_VAR(ncid, var_glat, grid_g(1)%glat(1,:),start=[1],count=[m2])
  if (status /= NF90_NOERR) CALL handle_error(status, "put lat")
  status = NF90_PUT_VAR(ncid, var_glon, grid_g(1)%glon(:,1),start=[1],count=[m3])
  if (status /= NF90_NOERR) CALL handle_error(status, "put lon")

  ! Write Species Names
  DO n = 1, nspecies
    str_len = MAX(1, LEN_TRIM(spc_name(n)))
    status = NF90_PUT_VAR(ncid, var_spc_name, TRIM(spc_name(n)), start=(/1, n/), count=(/str_len, 1/))
  if (status /= NF90_NOERR) CALL handle_error(status, "put name")
  END DO

  ! Write 4D Species Data
  DO n = 1, nspecies
    status = NF90_PUT_VAR(ncid, var_scp, chem1_g(n)%sc_p, start=(/1, 1, 1, n/), count=(/m1, m2, m3, 1/))
  if (status /= NF90_NOERR) CALL handle_error(status, "put concentration")
    status = NF90_PUT_VAR(ncid, var_sct, chem1_g(n)%sc_t, start=(/1, 1, 1, n/), count=(/m1, m2, m3, 1/))
  if (status /= NF90_NOERR) CALL handle_error(status, "put tendency")
    status = NF90_PUT_VAR(ncid, var_sctdyn, chem1_g(n)%sc_t_dyn, start=(/1, 1, 1, n/), count=(/m1, m2, m3, 1/))
  if (status /= NF90_NOERR) CALL handle_error(status, "put dynamics")
  END DO

  status = NF90_CLOSE(ncid)
  if (status /= NF90_NOERR) CALL handle_error(status, "Close File")

CONTAINS
  subroutine handle_error(status, loc)
    integer, intent(in) :: status
    character(len=*), intent(in) :: loc
    IF (status /= NF90_NOERR) THEN
        PRINT *, 'Output NetCDF Error [', TRIM(loc), ']: ', TRIM(NF90_STRERROR(status))
        STOP
    END IF
  END subroutine handle_error
END subroutine create_netcdf_file

!Creates a 20-33% larger netcdf file, but each species type is its own variable
subroutine create_netcdf_file_per_species(time, processor, m1, m2, m3, nspecies, chem1_g)
  USE netcdf
  IMPLICIT NONE

  integer, intent(in) :: m1, m2, m3, nspecies, processor
  real, intent(in)    :: time
  type(chem1_vars), intent(inout) :: chem1_g(nspecies)

  integer :: ncid, status
  integer :: dim_m1, dim_m2, dim_m3
  integer :: var_glat, var_glon
  
  ! Arrays to hold the variable IDs for each species
  integer, allocatable :: var_scp(:), var_sct(:), var_sctdyn(:)
  integer :: i, n
  character(len=64) :: filename
  character(len=128) :: var_name

  ! Create filename
  write(filename, fmt='("chem_ia/chem_out_",I6.6,"-",I4.4,".nc4")') int(time), processor

  status = NF90_CREATE(TRIM(filename), OR(NF90_NETCDF4, NF90_CLOBBER), ncid)
    if (status /= NF90_NOERR) CALL handle_error(status, " ")

  ! Define dimensions
  status = NF90_DEF_DIM(ncid, 'level', m1, dim_m1)
    if (status /= NF90_NOERR) CALL handle_error(status, " ")
  status = NF90_DEF_DIM(ncid, 'latitude', m2, dim_m2)
    if (status /= NF90_NOERR) CALL handle_error(status, " ")
  status = NF90_DEF_DIM(ncid, 'longitude', m3, dim_m3)
    if (status /= NF90_NOERR) CALL handle_error(status, " ")

  ! Define Coordinates
  status = NF90_DEF_VAR(ncid, 'latitude', NF90_FLOAT, dim_m2, var_glat)
    if (status /= NF90_NOERR) CALL handle_error(status, " ")
  status = NF90_PUT_ATT(ncid, var_glat, 'units', 'degrees_north')
    if (status /= NF90_NOERR) CALL handle_error(status, " ")

  status = NF90_DEF_VAR(ncid, 'longitude', NF90_FLOAT,  dim_m3, var_glon)
    if (status /= NF90_NOERR) CALL handle_error(status, " ")
  status = NF90_PUT_ATT(ncid, var_glon, 'units', 'degrees_east')
    if (status /= NF90_NOERR) CALL handle_error(status, " ")

  ! Allocate arrays for species variable IDs
  allocate(var_scp(nspecies), var_sct(nspecies), var_sctdyn(nspecies))

  ! Define Variables per species WITH DEFLATE AND SHUFFLE
  DO n = 1, nspecies
    ! Define concentration
    var_name = TRIM(spc_name(n)) // '_concentration'
    status = NF90_DEF_VAR(ncid, TRIM(var_name), NF90_FLOAT, (/dim_m1, dim_m2, dim_m3/), var_scp(n))
    if (status /= NF90_NOERR) CALL handle_error(status, " ")
    status = NF90_DEF_VAR_DEFLATE(ncid, var_scp(n), shuffle=1, deflate=1, deflate_level=1)
    if (status /= NF90_NOERR) CALL handle_error(status, " ")

    ! Define tendency
    var_name = TRIM(spc_name(n)) // '_tendency'
    status = NF90_DEF_VAR(ncid, TRIM(var_name), NF90_FLOAT, (/dim_m1, dim_m2, dim_m3/), var_sct(n))
    if (status /= NF90_NOERR) CALL handle_error(status, " ")
    status = NF90_DEF_VAR_DEFLATE(ncid, var_sct(n), shuffle=1, deflate=1, deflate_level=1)
    if (status /= NF90_NOERR) CALL handle_error(status, " ")

    ! Define dynamics
    var_name = TRIM(spc_name(n)) // '_dynamics'
    status = NF90_DEF_VAR(ncid, TRIM(var_name), NF90_FLOAT, (/dim_m1, dim_m2, dim_m3/), var_sctdyn(n))
    if (status /= NF90_NOERR) CALL handle_error(status, " ")
    status = NF90_DEF_VAR_DEFLATE(ncid, var_sctdyn(n), shuffle=1, deflate=1, deflate_level=1)
    if (status /= NF90_NOERR) CALL handle_error(status, " ")
  END DO

  ! End Definitions
  status = NF90_ENDDEF(ncid)
  if (status /= NF90_NOERR) CALL handle_error(status, "End Def")

  ! Write Grid Data
  status = NF90_PUT_VAR(ncid, var_glat, grid_g(1)%glat(1,:), start=[1], count=[m2])
  if (status /= NF90_NOERR) CALL handle_error(status, "put lat")
  status = NF90_PUT_VAR(ncid, var_glon, grid_g(1)%glon(:,1), start=[1], count=[m3])
  if (status /= NF90_NOERR) CALL handle_error(status, "put lon")

  ! Write 3D Species Data
  DO n = 1, nspecies
    status = NF90_PUT_VAR(ncid, var_scp(n), chem1_g(n)%sc_p, start=(/1, 1, 1/), count=(/m1, m2, m3/))
    if (status /= NF90_NOERR) CALL handle_error(status, "put concentration")
    status = NF90_PUT_VAR(ncid, var_sct(n), chem1_g(n)%sc_t, start=(/1, 1, 1/), count=(/m1, m2, m3/))
    if (status /= NF90_NOERR) CALL handle_error(status, "put tendency")
    status = NF90_PUT_VAR(ncid, var_sctdyn(n), chem1_g(n)%sc_t_dyn, start=(/1, 1, 1/), count=(/m1, m2, m3/))
    if (status /= NF90_NOERR) CALL handle_error(status, "put dynamics")

  END DO

  status = NF90_CLOSE(ncid)
    if (status /= NF90_NOERR) CALL handle_error(status, " ")
  
  ! Clean up allocations
  deallocate(var_scp, var_sct, var_sctdyn)

CONTAINS
  subroutine handle_error(status, loc)
    integer, intent(in) :: status
    character(len=*), intent(in) :: loc
    IF (status /= NF90_NOERR) THEN
        PRINT *, 'Output NetCDF Error [', TRIM(loc), ']: ', TRIM(NF90_STRERROR(status))
        STOP
    END IF
  END subroutine handle_error
END subroutine create_netcdf_file_per_species

!---------------------------------------------------------------------------------------------------
END MODULE mod_chem_spack_rodas3_dyndt
