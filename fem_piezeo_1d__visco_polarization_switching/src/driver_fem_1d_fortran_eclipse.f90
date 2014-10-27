! ============================================================================
! name        : driver_fem_1d_fortran_eclipse.f90
! author      : amir
! version     :
! copyright   : your copyright notice
! description : hello world in fortran
! ============================================================================
program driver_fem_1d_fortran_eclipse

use linearsolvers
use fem_geometry
use material_behavior
use fem_libs

implicit none
!     ================================================================
!                           the variables
!     ================================================================

!     ================================================================
!                      solution variables
!     ================================================================
      real(iwp),allocatable::glk(:,:)

      real(iwp),allocatable::glu(:) !the global solution
      real(iwp),allocatable::glp(:) !the previus degrees of freedom vector
      real(iwp),allocatable::glb(:) !the before previus degrees of freedom
!     ================================================================
!                      solver variables
!     ================================================================

      real(iwp),allocatable::glq(:) !internal force vector
      real(iwp),allocatable::glt(:) !external fource vector
      real(iwp),allocatable::glr(:) !total residual vector
!     ============================iteration variables
      real(iwp)::  eps,error,normal !the error tolerance and error
      integer:: itmax !maxumim number of iteration
!      integer:: iter ! iteration counter
      logical:: converge
!     ================================================================
!                      file unit
!     ================================================================
!      integer::curved_node !curved node for
!     ================================================================
!                      time variables
!     ================================================================
!     integer::ntime
!     real(iwp)::deltime
      real(iwp)::loadfactor    ! time
      real(iwp)::freq
!      real*4:: telapse ! timearray,
!      integer ( kind = 4 )start_stamp_v(8),end_stamp_v(8),values(8)
!      character*8::ctime(2)
      integer :: clck_counts_beg, clck_counts_end, clck_rate
      real ::  beg_cpu_time, end_cpu_time
!     ================================================================
!      integer::inode,jnode,ibond,jbond ,i,j,id,in
!     ================================================================
!      integer::igauss
!     =======================trivial meshing arrays the meshing seed parameters
      integer ::neldirectional(dimen)
      real(iwp)  ::length(dimen)
!     ================================================================
!     ===============================================================
!                       p r e p r o c e s s o r   u n i t
!     ===============================================================
      call system_clock ( clck_counts_beg, clck_rate)
      call cpu_time (beg_cpu_time)
      call timestamp()
!     ===============================================================
!     reading iteration and convergence critria's
!     ===============================================================
      itmax=30;eps=1e-2
!     ===============================================================
!     reading time increment varibales
!     ===============================================================
      dtime=0.01d0; freq=0.25d0
      ntime=int(24.0/dtime)

      length=10.0d0 ! mm
      neldirectional(1)=5
      call geometry_1D(neldirectional,length)
!     ===============================================================
!     define the solution parameters
!     ===============================================================
      nnm=size(coords,dim=1)
      nem=size(nod,dim=1)
      neq  = nnm*ndf;     nn=npe*ndf;
      ipdf = iel+1;    ngauss=(ipdf**dimen)*nem
!     ===============================================================
!     reading the boundary conditions
!     ===============================================================
      call form_history(ngauss)

      call bounday_1D()
      allocate(glk(neq,neq),glu(neq),glq(neq),glt(neq),glr(neq),glp(neq),glb(neq))
      glu=0.0d0;glt=0.0d0;glr=0.0d0;glp=0.0d0;glb=0.0d0
!     ===============================================================
!                        time increment starts here
!     ===============================================================
      vspvt=vspv
!      vssvt=vssv
     do itime=0,ntime
     time(1)=dtime;time(2)=itime*dtime

!     loadfactor=sawtooth(time(2)*freq)*100 ! * 100e3
     loadfactor=sin(2*3.14515*freq*time(2))*100
     vspv=loadfactor*vspvt
!     write(3,*)vspv
!    vssv=loadfactor*vssvt


      glu=0.0d0
!      glu(bnd_no_pr_vec)=0.0d0;
!      glu(bnd_no_pr_vec)=vspv;
      glt=0.0d0;
!     ===============================================================
!                 nonlinear solution iteration starts here
!     ===============================================================
      iter=0;CONVERGE=.false.;
130   iter=iter+1;
      if(iter.gt.itmax)then; write(*,330);goto 500;endif
!     ===============================================================
!                        forming the global matrices
!     ===============================================================
      glk=0.0d0;glq=0.0d0;!
      call glbmatrcs(glu,glk,glq,glp,glb)
!     ===============================================================
!                             newton raphson sprocedure
!     ===============================================================
      glt(bnd_no_se_vec) = vssv ;
      glr=glt-glq ;
!     ===============================================================
!                        solving the linear system of equations
!     ===============================================================
      call symmetric_primary_bounday(glk,glr)
      vspv=0.0d0
      call gaussian_elimination_solver(glk,glr)
!      call lapack_gesv( glk, glr )
!     ===============================================================
!                        error forming
!     ===============================================================
      normal=norm_vect(glu);if(normal.lt.default_smallest_pivot)then;normal=1;endif
      error=norm_vect(glr)/normal
!     ===============================================================
!                        updating the solution
!     ===============================================================
      glu=glu+glr

      write(out,*)'time',time
      write(out,*)'error',error
      write(out,*)'normal',normal
      write(out,*)'iter',iter
!
      if(error.gt.eps)then;converge=.false. ;goto 130;endif;
      converge=.true.

    call result_printer_1D(iter,glu)
    write(7,*)time(2),curn_electric_field(2,1),curn_polarization_function(2,1)

    call update_history()
    glb=glp;glp=glu

    enddo !itime=0,ntime


500   call system_clock ( clck_counts_end, clck_rate )
      write (*, *)'elapsed system clock=', &
      (clck_counts_end - clck_counts_beg) / real (clck_rate)
      call cpu_time (end_cpu_time)
      write (*, *)'elapsed cpu clock=', end_cpu_time - beg_cpu_time
      call timestamp ()

330   FORMAT(/,5X,'***** CONVERGENCE CRITERION IS NOT SATISFIED *****')
end program

