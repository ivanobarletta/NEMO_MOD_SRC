MODULE solpetsc
   !!======================================================================
   !!                     ***  MODULE  solfet
   !! Ocean solver :  preconditionned conjugate gradient solver
   !!=====================================================================

   !!----------------------------------------------------------------------
   !!   sol_petsc    : petsc krylov subspace solver
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables 
   USE sol_oce         ! ocean solver variables
   USE lib_mpp         ! distributed memory computing
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager  ! I/O manager
   USE lib_fortran     ! Fortran routines library
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sol_petsc_init, sol_petsc_finalize    ! 
   PUBLIC   petsc_set_mat , petsc_set_rhs, petsc_solve    ! 

   !! * Substitutions
#  include "vectopt_loop_substitute.h90"
   !! Petsc include files     
#  include "petsc/finclude/petscsys.h"
#  include "petsc/finclude/petscvec.h"
#  include "petsc/finclude/petscmat.h"
#  include "petsc/finclude/petscksp.h"
#  include "petsc/finclude/petscpc.h"
#  include "petsc/finclude/petscviewer.h"
#  include "petsc/finclude/petscvec.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: solpcg.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   PetscInt       ::   vecsize,nnz, matnnz, vecnnz, N
   PetscInt       ::   loc_matrows, glob_matrows, matcols, vecrows, veccols
   PetscInt       ::   i, xsize, m, one, col_shift
   PetscInt       ::   three, four, five, two
   PetscInt       ::   ji,jj, ishift, jshift, nx, ny
   PetscInt       ::   LI, col5(5), col4(4), col3(3), row
   PetscErrorCode ::   ierr
   PetscBool      ::   flg, print_out
   PetscScalar    ::   dot,ione, values5(5), values4(4), values3(3)
   PetscReal      ::   norm,rdot, zero, tol, minus_one
   Vec            ::   x,rhs,residual, mv
   Mat            ::   Amat
   KSP            ::   ksp
   PC             ::   pc
   PetscScalar, pointer ::  xx_v(:), rr_v(:), mv_v(:)
   PetscMPIInt    ::   rank, comm_size
   PetscDraw      ::   draw           
   PetscViewer    ::   vview, resview
   PetscViewer    ::   mview  ! DO NOT USE matview. IT LOOKS 
                              ! LIKE IS A PRIVATE PETSC WORD

CONTAINS
   
   SUBROUTINE sol_petsc_init
      CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)     

      ! global dimension of matrix
      N = jpiglo * jpjglo 
      ! # of local rows of matrix
      nx = nx_self(narea)
      ny = ny_self(narea)

      m = nx * ny

      !ALLOCATE( petsc_prealloc_d( m ) )  
      !ALLOCATE( petsc_prealloc_o( m ) )  

      write(6,*) 'count bmask ==0 ',count( bmask == 0 )   
      write(6,*) 'count bmask /=0 ',count( bmask /= 0 )   
      write(6,*) 'count gcp_1 ==0 ',count( gcp2(:,:,1) == 0 )   
      write(6,*) 'count gcp_1 /=0 ',count( gcp2(:,:,1) /= 0 )   
      write(6,*) 'count gcp_2 ==0 ',count( gcp2(:,:,2) == 0 )   
      write(6,*) 'count gcp_2 /=0 ',count( gcp2(:,:,2) /= 0 )   
      write(6,*) 'count gcp_3 ==0 ',count( gcp2(:,:,3) == 0 )   
      write(6,*) 'count gcp_3 /=0 ',count( gcp2(:,:,3) /= 0 )   
      write(6,*) 'count gcp_4 ==0 ',count( gcp2(:,:,4) == 0 )   
      write(6,*) 'count gcp_4 /=0 ',count( gcp2(:,:,4) /= 0 )   


      !petsc_prealloc_d = 5 
      !petsc_prealloc_o = 2  

      !DO jj = nldj, nlej
      !   DO ji = nldi, nlei
      !      LI = (jj-1-jshift)*nx + ji-ishift
      !      IF ( bmask(ji,jj) /= 0 ) THEN
      !         petsc_prealloc_d( LI ) = 5     
      !         petsc_prealloc_o( LI ) = 2     
      !      END IF    
      !   END DO
      !END DO

      !write(6,*) 'sum of prealloc diag  ', sum(petsc_prealloc_d)  
      !write(6,*) 'sum of prealloc off d ', sum(petsc_prealloc_o)  

      one   = 1  
      two   = 2    
      three = 3
      four  = 4    
      five  = 5

      ishift = nldi - 1
      jshift = nldj - 1

      !CALL MatCreate( mpi_comm_opa,A,ierr)
      CALL MatCreate( PETSC_COMM_WORLD,Amat,ierr)
      CALL MatSetSizes(Amat,m,m,N,N,ierr)
      CALL MatSetType(Amat, MATMPIAIJ,ierr)
      CALL MatSetFromOptions(Amat,ierr)

      CALL MatMPIAIJSetPreallocation(Amat,five,PETSC_NULL_INTEGER,two,PETSC_NULL_INTEGER,ierr)
      !CALL MatMPIAIJSetPreallocation(Amat,five,petsc_prealloc_d,two,petsc_prealloc_o,ierr)

      CALL VecCreate( PETSC_COMM_WORLD, rhs, ierr)
      CALL VecSetSizes(rhs, m, N, ierr)
      CALL VecSetFromOptions(rhs,ierr)
      CALL VecDuplicate(rhs,x, ierr)  

      CALL KSPCreate( PETSC_COMM_WORLD, ksp, ierr)
      CALL KSPSetOperators(ksp, Amat, Amat, ierr)
      !flg = PETSC_TRUE
      !CALL KSPSetInitialGuessNonzero(ksp, flg, ierr)
      !CALL KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED, ierr)
      CALL KSPSetFromOptions(ksp, ierr)

      WRITE(6,*) 'Petsc initialized' 
 
      911  FORMAT(a6,2i5,a5,i6,a6,5i6)
      9111 FORMAT(a6,2i5,a5,i6,a6,5i6,5e15.2)
      912  FORMAT(a6,2i5,a5,i6,a6,4i6)
      913  FORMAT(a6,2i5,a5,i6,a6,3i6)
      9134 FORMAT(f27.18)
      1985 FORMAT(a5,2i4,1x,f35.20)
      1984 FORMAT(2i4,1x,f38.20)
      9898 FORMAT(a5,2i3,a4,i4,6(a2,e14.7))
      9899 FORMAT(a5,2i3,a8,a14,5(a2,e14.7))
      !! --  
   END SUBROUTINE sol_petsc_init               

   SUBROUTINE sol_petsc_finalize
      !! Finalize Petsc Package  
      WRITE(6,*) 'Petsc finalized'  

      !DEALLOCATE( petsc_prealloc_d )  
      !DEALLOCATE( petsc_prealloc_o )  

      CALL MatDestroy(Amat, ierr)
      CALL VecDestroy(rhs, ierr)
      CALL VecDestroy(x, ierr)
      CALL KSPDestroy(ksp, ierr)
      CALL PetscFinalize(ierr)     

   END SUBROUTINE sol_petsc_finalize               

   !!=====================================================================

   SUBROUTINE petsc_set_mat  
      !! Set Matrix Values   
      911  FORMAT(a6,2i5,a5,i6,a6,5i6)
      912  FORMAT(a6,2i5,a5,i6,a6,4i6)

      nx = nx_self(narea)
      ny = ny_self(narea)

      one   = 1  
      two   = 2    
      three = 3
      four  = 4    
      five  = 5

      ishift = nldi - 1
      jshift = nldj - 1
  
      !! --------------------------------------- !!
      !!    TREATING INSIDE OF INNER DOMAINS     !!
      !! (BORDERS AND CORNERS ARE TREATED LATER) !!
      !! --------------------------------------- !!
  
      DO jj = nldj+1, nlej-1
         DO ji = nldi+1, nlei-1
            LI = (jj-1-jshift)*nx + ji-ishift
            row = RB(narea) + LI - 1
            IF ( bmask(ji,jj) /= 0 ) THEN    
               !                SOUTH         WEST          DIAG           EAST          NORTH            
               col5(1:5)    = (/row-nx       ,row-1        ,row           ,row+1        ,row+nx       /)
               values5(1:5) = (/gcp2(ji,jj,1),gcp2(ji,jj,2),gcdmat2(ji,jj),gcp2(ji,jj,3),gcp2(ji,jj,4)/)
               !WRITE(6,911) 'ji ,jj' ,ji, jj, 'row', row,'col5 ', col5
               CALL MatSetValues(Amat,one,row,five,col5, values5, INSERT_VALUES, ierr)
            ELSE
               CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
            END IF         
         END DO
      END DO

      !! ------------------------------------ !!
      !! TREATING BORDERS ( WITHOUT CORNERS ) !!
      !! ------------------------------------ !!

      ! < EAST BORDER >
      IF ( someone_east(narea) == .TRUE. ) THEN
         ji = nlei ! eastern col
         DO jj = nldj+1,nlej-1
            LI = (jj-1-jshift)*nx + ji-ishift
            row  = RB(narea) + LI - 1
            IF ( bmask(ji,jj) /= 0 ) THEN    
               !                SOUTH         WEST          DIAG           EAST                   NORTH            
               col5(1:5)    = (/row-nx       ,row-1        ,row           ,row+col_shift_east(jj),row+nx       /)
               values5(1:5) = (/gcp2(ji,jj,1),gcp2(ji,jj,2),gcdmat2(ji,jj),gcp2(ji,jj,3)         ,gcp2(ji,jj,4)/)
               !WRITE(6,911) 'ji ,jj' ,ji, jj, 'row', row,'col5 ', col5
               CALL MatSetValues(Amat,one,row,five,col5, values5, INSERT_VALUES, ierr) 
            ELSE
               CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
            END IF        
         END DO
      ELSE
         ! EASTERN EDGE OF GLOBAL DOMAIN
         ji = nlei ! eastern col
         DO jj = nldj+1,nlej-1
            LI = (jj-1-jshift)*nx + ji-ishift
            row  = RB(narea) + LI - 1
            IF ( bmask(ji,jj) /= 0 ) THEN    
               !                SOUTH         WEST          DIAG           NORTH
               col4(1:4)    = (/row-nx       ,row-1        ,row           ,row+nx       /)
               values4(1:4) = (/gcp2(ji,jj,1),gcp2(ji,jj,2),gcdmat2(ji,jj),gcp2(ji,jj,4)/)
               !WRITE(6,912) 'ji ,jj' ,ji, jj, 'row', row,'col4 ', col4
               CALL MatSetValues(Amat,one,row,four,col4, values4, INSERT_VALUES, ierr) 
            ELSE
               CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
            END IF        
         END DO  
      END IF

      ! < WEST BORDER >
      IF ( someone_west(narea) == .TRUE. ) THEN
         ji = nldi ! western col
         DO jj = nldj+1,nlej-1
            LI = (jj-1-jshift)*nx + ji-ishift
            row  = RB(narea) + LI - 1
            IF ( bmask(ji,jj) /= 0 ) THEN    
               !                SOUTH         WEST                   DIAG           EAST          NORTH
               col5(1:5)    = (/row-nx       ,row-col_shift_west(jj),row           ,row+1        ,row+nx       /)
               values5(1:5) = (/gcp2(ji,jj,1),gcp2(ji,jj,2)         ,gcdmat2(ji,jj),gcp2(ji,jj,3),gcp2(ji,jj,4)/)
               !WRITE(6,911) 'ji ,jj' ,ji, jj, 'row', row,'col5 ', col5
               CALL MatSetValues(Amat,one,row,five,col5, values5, INSERT_VALUES, ierr) 
            ELSE
               CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
            END IF        
         END DO  
      ELSE
         ! WESTERN EDGE OF GLOBAL DOMAIN 
         ji = nldi ! western col
         DO jj = nldj+1,nlej-1
            LI  = (jj-1-jshift)*nx + ji-ishift
            row = RB(narea) + LI - 1
            IF ( bmask(ji,jj) /= 0 ) THEN    
               !                SOUTH         DIAG           EAST          NORTH
               col4(1:4)    = (/row-nx       ,row           ,row+1        ,row+nx/)
               values4(1:4) = (/gcp2(ji,jj,1),gcdmat2(ji,jj),gcp2(ji,jj,3),gcp2(ji,jj,4)/)
               !WRITE(6,912) 'ji ,jj' ,ji, jj, 'row', row,'col4 ', col4
               !write(6,*) 'gcdmat2 western ', gcdmat2(ji,jj)    
               CALL MatSetValues(Amat,one,row,four,col4, values4, INSERT_VALUES, ierr) 
            ELSE
               CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
            END IF        
         END DO  
      END IF
 
      ! < NORTH BORDER >
      IF ( someone_north(narea) == .TRUE. ) THEN
         jj = nlej ! northern row
         DO ji = nldi+1, nlei-1
            LI = (jj-1-jshift)*nx + ji-ishift
            row  = RB(narea) + LI - 1
            IF ( bmask(ji,jj) /= 0 ) THEN    
               !                SOUTH         WEST          DIAG           EAST          NORTH
               col5(1:5)    = (/row-nx       ,row-1        ,row           ,row+1        ,row+col_shift_north(ji)/)
               values5(1:5) = (/gcp2(ji,jj,1),gcp2(ji,jj,2),gcdmat2(ji,jj),gcp2(ji,jj,3),gcp2(ji,jj,4)          /)
               !WRITE(6,911) 'ji ,jj' ,ji, jj, 'row', row,'col5 ', col5
               CALL MatSetValues(Amat,one,row,five,col5, values5, INSERT_VALUES, ierr) 
            ELSE
               CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
            END IF        
         END DO
      ELSE
         ! NORTHERN EDGE OF DOMAIN
         jj = nlej ! northern row
         DO ji = nldi+1, nlei-1
            LI = (jj-1-jshift)*nx + ji-ishift
            row  = RB(narea) + LI - 1
            IF ( bmask(ji,jj) /= 0 ) THEN    
               !                SOUTH         WEST          DIAG           EAST
               col4(1:4)    = (/row-nx       ,row-1        ,row           ,row+1        /)
               values4(1:4) = (/gcp2(ji,jj,1),gcp2(ji,jj,2),gcdmat2(ji,jj),gcp2(ji,jj,3)/)
               !WRITE(6,912) 'ji ,jj' ,ji, jj, 'row', row,'col4 ', col4
               !write(6,*) 'gcdmat2 northern ', gcdmat2(ji,jj)    
               CALL MatSetValues(Amat,one,row,four,col4, values4, INSERT_VALUES, ierr) 
            ELSE
               CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
            END IF        
         END DO
      END IF

      ! < SOUTH BORDER >
      IF ( someone_south(narea) == .TRUE. ) THEN
         jj = nldj ! southern row
         DO ji = nldi+1, nlei-1
            LI = (jj-1-jshift)*nx + ji-ishift
            row  = RB(narea) + LI - 1
            IF ( bmask(ji,jj) /= 0 ) THEN    
               !                SOUTH                   WEST          DIAG           EAST          NORTH
               col5(1:5)    = (/row-col_shift_south(ji),row-1        ,row           ,row+1        ,row+nx        /)
               values5(1:5) = (/gcp2(ji,jj,1)          ,gcp2(ji,jj,2),gcdmat2(ji,jj),gcp2(ji,jj,3),gcp2(ji,jj,4) /)
               !WRITE(6,911) 'ji ,jj' ,ji, jj, 'row', row,'col5 ', col5
               CALL MatSetValues(Amat,one,row,five,col5, values5, INSERT_VALUES, ierr) 
            ELSE
               CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
            END IF        
         END DO
      ELSE
         ! SOUTHERN EDGE OF GLOBAL DOMAIN
         jj = nldj ! southern row
         DO ji = nldi+1, nlei-1
            LI = (jj-1-jshift)*nx + ji-ishift
            row  = RB(narea) + LI - 1
            IF ( bmask(ji,jj) /= 0 ) THEN    
               !                WEST          DIAG           EAST          NORTH
               col4(1:4)    = (/row-1        ,row           ,row+1        ,row+nx        /)
               values4(1:4) = (/gcp2(ji,jj,2),gcdmat2(ji,jj),gcp2(ji,jj,3),gcp2(ji,jj,4) /)
               !WRITE(6,912) 'ji ,jj' ,ji, jj, 'row', row,'col4 ', col4
               !write(6,*) 'gcdmat2 southern ', gcdmat2(ji,jj)    
               CALL MatSetValues(Amat,one,row,four,col4, values4, INSERT_VALUES, ierr) 
            ELSE
               CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
            END IF        
         END DO
      END IF
     
      !! -----------------!!
      !! TREATING CORNERS !!
      !! -----------------!!

      ! < SOUTH-WEST CORNER >
      ji = nldi
      jj = nldj
      LI = (jj-1-jshift)*nx + ji-ishift
      row  = RB(narea) + LI - 1
      IF ( someone_south(narea) == .TRUE. .AND. someone_west(narea) == .TRUE. ) THEN
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                SOUTH                   WEST                   DIAG           EAST          NORTH
            col5(1:5)    = (/row-col_shift_south(ji),row-col_shift_west(jj),row           ,row+1        ,row+nx        /)
            values5(1:5) = (/gcp2(ji,jj,1)          ,gcp2(ji,jj,2)         ,gcdmat2(ji,jj),gcp2(ji,jj,3),gcp2(ji,jj,4) /)
            !WRITE(6,911) 'ji ,jj' ,ji, jj, 'row', row,'col5 ', col5
            CALL MatSetValues(Amat,one,row,five,col5, values5, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      ELSE IF ( someone_south(narea) == .TRUE. .AND. someone_west(narea) == .FALSE. ) THEN
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                SOUTH                   DIAG           EAST          NORTH
            col4(1:4)    = (/row-col_shift_south(ji),row           ,row+1        ,row+nx       /)
            values4(1:4) = (/gcp2(ji,jj,1)          ,gcdmat2(ji,jj),gcp2(ji,jj,3),gcp2(ji,jj,4)/)
            !WRITE(6,912) 'ji ,jj' ,ji, jj, 'row', row,'col4 ', col4
            CALL MatSetValues(Amat,one,row,four,col4, values4, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      ELSE IF ( someone_south(narea) == .FALSE. .AND. someone_west(narea) == .TRUE. ) THEN
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                WEST                   DIAG           EAST          NORTH
            col4(1:4)    = (/row-col_shift_west(jj),row           ,row+1        ,row+nx        /)
            values4(1:4) = (/gcp2(ji,jj,2)         ,gcdmat2(ji,jj),gcp2(ji,jj,3),gcp2(ji,jj,4) /)
            !WRITE(6,912) 'ji ,jj' ,ji, jj, 'row', row,'col4 ', col4
            CALL MatSetValues(Amat,one,row,four,col4, values4, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      ELSE ! there isn't any process either at south or west 
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                DIAG           EAST          NORTH
            col3(1:3)    = (/row           ,row+1        ,row+nx       /)
            values3(1:3) = (/gcdmat2(ji,jj),gcp2(ji,jj,3),gcp2(ji,jj,4)/)
            !WRITE(6,913) 'ji ,jj' ,ji, jj, 'row', row,'col3 ', col3
            CALL MatSetValues(Amat,one,row,three,col3, values3, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      END IF

      ! < SOUTH-EAST CORNER >
      ji = nlei
      jj = nldj
      LI = (jj-1-jshift)*nx + ji-ishift
      row  = RB(narea) + LI - 1
      IF ( someone_south(narea) == .TRUE. .AND. someone_east(narea) == .TRUE. ) THEN
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                SOUTH                   WEST          DIAG           EAST                   NORTH
            col5(1:5)    = (/row-col_shift_south(ji),row-1        ,row           ,row+col_shift_east(jj),row+nx       /)
            values5(1:5) = (/gcp2(ji,jj,1)          ,gcp2(ji,jj,2),gcdmat2(ji,jj),gcp2(ji,jj,3)         ,gcp2(ji,jj,4)/)
            !WRITE(6,911) 'ji ,jj' ,ji, jj, 'row', row,'col5 ', col5
            CALL MatSetValues(Amat,one,row,five,col5, values5, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      ELSE IF ( someone_south(narea) == .TRUE. .AND. someone_east(narea) == .FALSE. ) THEN
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                SOUTH                   WEST          DIAG           NORTH
            col4(1:4)    = (/row-col_shift_south(ji),row-1        ,row           ,row+nx       /)
            values4(1:4) = (/gcp2(ji,jj,1)          ,gcp2(ji,jj,2),gcdmat2(ji,jj),gcp2(ji,jj,4)/)
            !WRITE(6,912) 'ji ,jj' ,ji, jj, 'row', row,'col4 ', col4
            CALL MatSetValues(Amat,one,row,four,col4, values4, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      ELSE IF ( someone_south(narea) == .FALSE. .AND. someone_east(narea) == .TRUE. ) THEN
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                WEST          DIAG           EAST                   NORTH
            col4(1:4)    = (/row-1        ,row           ,row+col_shift_east(jj),row+nx       /) 
            values4(1:4) = (/gcp2(ji,jj,2),gcdmat2(ji,jj),gcp2(ji,jj,3)         ,gcp2(ji,jj,4)/) 
            !WRITE(6,912) 'ji ,jj' ,ji, jj, 'row', row,'col4 ', col4
            CALL MatSetValues(Amat,one,row,four,col4, values4, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      ELSE
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                WEST          DIAG           NORTH
            col3(1:3)    = (/row-1        ,row           ,row+nx       /)
            values3(1:3) = (/gcp2(ji,jj,2),gcdmat2(ji,jj),gcp2(ji,jj,4)/)
            !WRITE(6,913) 'ji ,jj' ,ji, jj, 'row', row,'col3 ', col3
            CALL MatSetValues(Amat,one,row,three,col3, values3, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      END IF

      ! < NORTH-WEST CORNER > 
      ji = nldi
      jj = nlej
      LI = (jj-1-jshift)*nx + ji-ishift
      row  = RB(narea) + LI - 1
      IF ( someone_north(narea) == .TRUE. .AND. someone_west(narea) == .TRUE. ) THEN
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                SOUTH         WEST                   DIAG           EAST          NORTH
            col5(1:5)    = (/row-nx       ,row-col_shift_west(jj),row           ,row+1        ,row+col_shift_north(ji)/)
            values5(1:5) = (/gcp2(ji,jj,1),gcp2(ji,jj,2)         ,gcdmat2(ji,jj),gcp2(ji,jj,3),gcp2(ji,jj,4)          /) 
            !WRITE(6,911) 'ji ,jj' ,ji, jj, 'row', row,'col5 ', col5
            CALL MatSetValues(Amat,one,row,five,col5, values5, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      ELSE IF ( someone_north(narea) == .FALSE. .AND. someone_west(narea) == .TRUE. ) THEN
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                SOUTH         WEST                   DIAG           EAST
            col4(1:4)    = (/row-nx       ,row-col_shift_west(jj),row           ,row+1        /)
            values4(1:4) = (/gcp2(ji,jj,1),gcp2(ji,jj,2)         ,gcdmat2(ji,jj),gcp2(ji,jj,3)/)
            !WRITE(6,912) 'ji ,jj' ,ji, jj, 'row', row,'col4 ', col4
            CALL MatSetValues(Amat,one,row,four,col4, values4, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      ELSE IF ( someone_north(narea) == .TRUE. .AND. someone_west(narea) == .FALSE. ) THEN
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                SOUTH         DIAG           EAST          NORTH
            col4(1:4)    = (/row-nx       ,row           ,row+1        ,row+col_shift_north(ji)/)
            values4(1:4) = (/gcp2(ji,jj,1),gcdmat2(ji,jj),gcp2(ji,jj,3),gcp2(ji,jj,4)          /)
            !WRITE(6,912) 'ji ,jj' ,ji, jj, 'row', row,'col4 ', col4
            CALL MatSetValues(Amat,one,row,four,col4, values4, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      ELSE
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                SOUTH         DIAG           EAST
            col3(1:3)    = (/row-nx       ,row           ,row+1        /)
            values3(1:3) = (/gcp2(ji,jj,1),gcdmat2(ji,jj),gcp2(ji,jj,3)/)
            !WRITE(6,913) 'ji ,jj' ,ji, jj, 'row', row,'col3 ', col3
            CALL MatSetValues(Amat,one,row,three,col3, values3, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      END IF

      ! < NORTH-EAST CORNER >
      ji = nlei
      jj = nlej
      LI = (jj-1-jshift)*nx + ji-ishift
      row  = RB(narea) + LI - 1
      IF ( someone_north(narea) == .TRUE. .AND. someone_east(narea) == .TRUE. ) THEN
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                SOUTH         WEST          DIAG           EAST                   NORTH
            col5(1:5)    = (/row-nx       ,row-1        ,row           ,row+col_shift_east(jj),row+col_shift_north(ji)/)
            values5(1:5) = (/gcp2(ji,jj,1),gcp2(ji,jj,2),gcdmat2(ji,jj),gcp2(ji,jj,3)         ,gcp2(ji,jj,4)          /)
            !WRITE(6,911) 'ji ,jj' ,ji, jj, 'row', row,'col5 ', col5
            CALL MatSetValues(Amat,one,row,five,col5, values5, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      ELSE IF ( someone_north(narea) == .TRUE. .AND. someone_east(narea) == .FALSE. ) THEN
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                 SOUTH         WEST          DIAG           NORTH
            col4(1:4)     = (/row-nx       ,row-1        ,row           ,row+col_shift_north(ji)/)
            values4(1:4)  = (/gcp2(ji,jj,1),gcp2(ji,jj,2),gcdmat2(ji,jj),gcp2(ji,jj,4)          /)
            !WRITE(6,912) 'ji ,jj' ,ji, jj, 'row', row,'col4 ', col4
            CALL MatSetValues(Amat,one,row,four,col4, values4, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      ELSE IF ( someone_north(narea) == .FALSE. .AND. someone_east(narea) == .TRUE. ) THEN
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                SOUTH         WEST          DIAG           EAST
            col4(1:4)    = (/row-nx       ,row-1        ,row           ,row+col_shift_east(jj)/)
            values4(1:4) = (/gcp2(ji,jj,1),gcp2(ji,jj,2),gcdmat2(ji,jj),gcp2(ji,jj,3)         /)
            !WRITE(6,912) 'ji ,jj' ,ji, jj, 'row', row,'col4 ', col4
            CALL MatSetValues(Amat,one,row,four,col4, values4, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      ELSE
         IF ( bmask(ji,jj) /= 0 ) THEN    
            !                SOUTH         WEST          DIAG 
            col3(1:3)    = (/row-nx       ,row-1        ,row           /)
            values3(1:3) = (/gcp2(ji,jj,1),gcp2(ji,jj,2),gcdmat2(ji,jj)/)
            !WRITE(6,913) 'ji ,jj' ,ji, jj, 'row', row,'col3 ', col3
            CALL MatSetValues(Amat,one,row,three,col3, values3, INSERT_VALUES, ierr) 
         ELSE
            CALL MatSetValue (Amat,row,row, gcdmat2(ji,jj), INSERT_VALUES, ierr)
         END IF        
      END IF

      CALL MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY,ierr)
      CALL MatAssemblyEnd  (Amat,MAT_FINAL_ASSEMBLY,ierr)

      !CALL PetscViewerASCIIOpen( PETSC_COMM_WORLD ,'ORCA2_PETSC_MAT',mview,ierr)
      !CALL MatView(Amat,mview,ierr)
      !CALL PetscViewerDestroy(mview,ierr)

   END SUBROUTINE petsc_set_mat

   SUBROUTINE petsc_set_rhs
      !! Set RHS   
      nx = nx_self(narea)

      !CALL VecSet(rhs,0.d0,ierr)  

      ishift = nldi - 1
      jshift = nldj - 1

      DO jj = nldj, nlej
         DO ji = nldi, nlei
            IF ( bmask(ji,jj) /= 0 ) THEN    
               LI = (jj-1-jshift)*nx + ji-ishift
               row = RB(narea) + LI - 1
               CALL VecSetValue(rhs,row,gcb2(ji,jj),INSERT_VALUES,ierr)
            END IF    
         END DO
      END DO

      CALL VecAssemblyBegin(rhs,ierr)
      CALL VecAssemblyEnd(rhs,ierr)

   END SUBROUTINE petsc_set_rhs          

   SUBROUTINE petsc_solve( kindic )
      INTEGER, INTENT(inout) ::   kindic   ! solver indicator, < 0 if the conver-
      !                                    ! gence is not reached: the model is stopped in step
      !                                    ! set to zero before the call of solver
      IF( nn_timing == 1 )  CALL timing_start('petsc_solve')

      ! < SET FIRST GUESS >
      CALL VecGetArrayF90(x,xx_v,ierr)
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            LI = (jj-1-jshift)*nx + ji-ishift
            row = LI - 1
            xx_v(row+1) = gcx(ji,jj)
         END DO
      END DO
      CALL VecRestoreArrayF90(x,xx_v,ierr)

      ! < SOLVE SYSTEM >
      CALL KSPSetInitialGuessNonzero(ksp, PETSC_TRUE, ierr)
      CALL KSPSolve(ksp, rhs, x, ierr)
      CALL KSPGetIterationNumber(ksp, niter, ierr)

      IF( niter == 0 ) kindic = -2

      ! < RETURN SOLUTION TO NEMO >
      CALL VecGetArrayReadF90(x,xx_v,ierr)
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            IF ( bmask(ji,jj) /= 0 ) THEN    
               LI = (jj-1-jshift)*nx + ji-ishift
               row = LI - 1
               gcx(ji,jj) = xx_v(row+1) 
            END IF      
         END DO
      END DO
      CALL VecRestoreArrayReadF90(x,xx_v,ierr)

      CALL lbc_lnk( gcx, c_solver_pt, 1. )      ! Output in gcx with lateral b.c. applied

      IF( nn_timing == 1 )  CALL timing_stop('petsc_solve')

   END SUBROUTINE petsc_solve          

        
END MODULE solpetsc
