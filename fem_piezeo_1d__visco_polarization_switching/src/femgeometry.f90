      module fem_geometry
      use linearsolvers
      implicit none


      integer::npe !number of nodes in element
      integer::nem !number of elements
      integer::nnm !number of nodes
      integer::iel   !element type tag
      integer::ipdf  !element type tag
      integer::neq !number of equations
      integer::nn !number of degrees of freedom in an element
      integer::nhbw !half band width
      integer::nw
      integer::nbw !band width
      real(iwp),parameter  :: pi = 3.14159265 !the pie number
!     ================================================================
!                      geormetry variables
!     ================================================================
      real(iwp),allocatable  ::coords(:,:) !the nodes geometry array
      integer,allocatable ::nod(:,:) !the connectivity matrix
!     ================================================================
!                      boundary condition variables
!     ================================================================
      integer,allocatable::ispv(:,:)!the primary boundary condition array
      real(iwp),allocatable::vspv(:),vspvt(:)!the primary boundary condition array
      integer,allocatable::issv(:,:)!the secondry boundary condition
      real(iwp),allocatable::vssv(:),vssvt(:)!the secondry boundary conditi

      integer,allocatable::bnd_no_pr_vec(:) !the secondry boundary conditi
      real(iwp),allocatable::bnd_va_pr_vec(:) !the secondry boundary conditi

      integer,allocatable::bnd_no_se_vec(:) !the secondry boundary conditi
      real(iwp),allocatable::bnd_va_se_vec(:) !the secondry boundary conditi

      integer::nspv,np,nb  !number of specified primary boundary condition
      integer::nssv  !number of specified secondry boundary condition


      contains

    subroutine geometry_1D(neldirectional,length)
    integer::neldirectional(:);
    real(iwp)::length(:)
    real(iwp)::delta_l
    integer::i,in,id

    delta_l=length(1)/neldirectional(1)
!     ===============================================================
!     mesh and fem information number of elements and number of nodes
!     ===============================================================

!     ===============================================================
!     define the solution parameters
!     ===============================================================
      iel=1;npe=2;
      nem=neldirectional(1)
      nnm=(iel*neldirectional(1) + 1)

    allocate(nod(nem,npe),coords(nnm,dimen))
    id=1;
    do in=1,nnm
    coords(in,id)=(in-1)*delta_l;enddo
!     ===============================================================
!     the connectivity matrix
!     ===============================================================
     do i=1,nem;nod(i,1)=1+(i-1);nod(i,2)=2+(i-1);enddo


      end subroutine geometry_1D

      subroutine bounday_1D()
      nspv=3       !reading primary variables

    if(nspv.gt.0)then
    allocate(ispv(nspv,2),vspv(nspv))
    ispv(1,:)=[1,1];
    ispv(2,:)=[1,2];
    vspv=0.0d0;

    ispv(3,:)=[nnm,2];
    vspv(3)=1.0d0;
      endif


     nssv=0!reading secondry variables
    if (nssv.gt.0)then
    allocate(issv(nssv,2),vssv(nssv))
    issv(1,:)=[nnm,1];
    vssv=0.0d0;
     endif


      allocate(bnd_va_pr_vec(nspv),bnd_no_pr_vec(nspv))
      allocate(bnd_va_se_vec(nssv),bnd_no_se_vec(nssv))

         do np=1,nspv
            nb=(ispv(np,1)-1)*ndf+ispv(np,2)
            bnd_va_pr_vec(np)=vspv(np)
            bnd_no_pr_vec(np)=nb
            enddo

         do np=1,nssv
           nb=(issv(np,1)-1)*ndf+issv(np,2)
           bnd_va_se_vec(np)=vssv(np)
           bnd_no_se_vec(np)=nb
         enddo



      end subroutine bounday_1D


      end module fem_geometry
