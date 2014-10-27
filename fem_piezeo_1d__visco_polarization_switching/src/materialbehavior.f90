module material_behavior
    use linearsolvers
    !     ================================================================
    real(iwp),parameter::k00=0.8d0
    real(iwp),parameter::k01=0.2d0
    real(iwp),parameter::lambda_01=5.0
    !     ================================================================
    real(iwp),parameter::ec=1;
    real(iwp),parameter::pr_sat=26e-2;
    real(iwp),parameter::eps_s_r=1e-3;
    real(iwp),parameter::c=1.0;
    real(iwp),parameter::eta=10e-2;
    real(iwp),parameter::m=2;
    !     ================================================================
    real(iwp),parameter::c11 = 13.9e4;
    real(iwp),parameter::c12 = 7.78e4;
    real(iwp),parameter::c13 = 7.43e4;
    real(iwp),parameter::c33 = 11.3e4;
    real(iwp),parameter::c44 = 2.56e4;
    real(iwp),parameter::e31 = 6.98;
    real(iwp),parameter::e33 = 13.84;
    real(iwp),parameter::e15 = 13.44;
    real(iwp),parameter::eps_11 = 6e-3;
    real(iwp),parameter::eps_33 = 5.47e-3;
    logical::is_polarized
    real(iwp)::a


    real(iwp)::lambda
    real(iwp)::mu
    real(iwp)::alpha1
    real(iwp)::alpha2
    real(iwp)::alpha3
    real(iwp)::beta1
    real(iwp)::beta2
    real(iwp)::beta3
    real(iwp)::gamma1
    real(iwp)::gamma2

    !     =========================element center point
    real(iwp) :: el_center_coord(dimen)
    real(iwp) :: qtransform(dimen,dimen) ! this is the transformtion matrix
    !     =======================mechanical properties variables
    real(iwp)::ctens(3,3,3,3)
    !     ===========================materials tensors
    real(iwp),dimension(3,3,3):: epz
    real(iwp),dimension(3,3)::ktense
    !     ================================================================
    !      time variables
    !     ================================================================
    real(iwp)::time(2),dtime
    integer::ntime,itime  ! time
    !     ================================================================
    integer::noelem
    !     ================================================================
    !      history variables
    !     ================================================================
    real(iwp),allocatable::remanent_polarization_vector(:,:)

    real(iwp),allocatable::hist_polarization_function(:,:)
    real(iwp),allocatable::curn_polarization_function(:,:)

    real(iwp),allocatable::mechanical_hist_curentt(:,:,:,:,:)
    real(iwp),allocatable::mechanical_hist_previus(:,:,:,:,:)
    real(iwp),allocatable::mechanical_hist_beforep(:,:,:,:,:)

    real(iwp),allocatable::mechanical_strain_curentt(:,:,:)
    real(iwp),allocatable::mechanical_strain_previus(:,:,:)
    real(iwp),allocatable::mechanical_strain_beforep(:,:,:)

    real(iwp),allocatable::hist_electric_field(:,:)
    real(iwp),allocatable::curn_electric_field(:,:)

    real(iwp),allocatable::hist_electric_displ(:,:)
    real(iwp),allocatable::curn_electric_displ(:,:)

    real(iwp)::total_displacement(dimen)
    real(iwp)::total_polarization(dimen)
    real(iwp)::pr
contains


    subroutine stress_elect_displacement(ur,up,ub,der,sigma)
        implicit none
        !     ================================================================
        !  input variables
        !     ================================================================
        real(iwp),intent(in) :: ur(:,:),up(:,:),ub(:,:) ! values of field function n point
        !     ================================================================
        ! output variables
        !     ================================================================
        integer :: i , j, k , l !integer counters
        real(iwp)  :: der(:),sigma(:,:)
        !     ================================================================
        real(iwp)::strain(dimen,dimen),strain_p(dimen,dimen),strain_b(dimen,dimen);
        real(iwp)::d_strain(dimen,dimen),d_strain_p(dimen,dimen)
        real(iwp)::d_electric_field(3),d_electric_field_p(3)
        real(iwp)::strain_r(dimen,dimen),strain_elastic(dimen,dimen);
        real(iwp)::electric_field(3),electric_field_p(3),electric_field_b(3);
        !     ================================================================
        lambda  =c12;
        mu      =(c11-c12)/2.0d0;
        alpha1  =2*c44+c12-c11;
        alpha2  =(c11+c33)/2.0d0-2.0d0*c44-c
        alpha3  =c13-c12;
        beta1   =-e31;
        beta2   =-e33+2.0d0*e15+e31;
        beta3   =-2.0d0*e15;
        gamma1  =-eps_11/2.0d0;
        gamma2  =(eps_11-eps_33)/2.0d0;

        !     ================================================================
        strain=0.0d0;electric_field=0.0d0;
        strain_p=0.0d0;electric_field_p=0.0d0;
        strain_b=0.0d0;electric_field_b=0.0d0;

        strain(1:dimen,1:dimen)=0.5d0*(   ur(1:dimen,:)+transpose(ur(1:dimen,:))+      &
            matmul(  transpose(   ur(1:dimen,:)  ),ur(1:dimen,:))    );

        strain_p(1:dimen,1:dimen)=0.5d0*(   up(1:dimen,:)+transpose(up(1:dimen,:))+      &
            matmul(  transpose(   up(1:dimen,:)  ),up(1:dimen,:))    );

        strain_b(1:dimen,1:dimen)=0.5d0*(   ub(1:dimen,:)+transpose(ub(1:dimen,:))+      &
            matmul(  transpose(   ub(1:dimen,:)  ),ub(1:dimen,:))    );

        electric_field(1:dimen)=-ur(dimen+1,:)
        electric_field_p(1:dimen)=-up(dimen+1,:)
        electric_field_b(1:dimen)=-ub(dimen+1,:)

        d_electric_field=electric_field-electric_field_p
        d_electric_field_p=electric_field_p-electric_field_b
        !     ================================================================
        curn_electric_field(gauss_point_number,1)=electric_field(1)
        !     ================================================================
        call remanent_polarization(d_electric_field,electric_field)
        pr= curn_polarization_function(gauss_point_number,1)
!        pr=pr_sat
        call direction_polarization(pr,a)
        call material_properties()
        !     ===================================remanent strain
        strain_r=0.0d0
        strain_r(1,1)=eps_s_r*dabs(pr)/pr_sat
        strain_elastic=strain-strain_r
        !     ==============================================increment in strain
        mechanical_strain_curentt(gauss_point_number,:,:)=strain_elastic
        d_strain=mechanical_strain_curentt(gauss_point_number,:,:)-mechanical_strain_previus(gauss_point_number,:,:)
        d_strain_p=mechanical_strain_previus(gauss_point_number,:,:)-mechanical_strain_beforep(gauss_point_number,:,:)
        !     ==============================the history variable
        do i=1,dimen
        do j=1,dimen
        do k=1,dimen
        do l=1,dimen
        mechanical_hist_curentt(gauss_point_number,i,j,k,l)= &
        mechanical_hist_previus(gauss_point_number,i,j,k,l)*exp(-lambda_01*dtime)+ &
        ctens(i,j,k,l)*(k01*0.5d0)* &
        (  &
        exp(-lambda_01*dtime)*d_strain_p(k,l)+d_strain(k,l) &
        )
        end do;
        end do;
        end do;
        end do;




        do i=1,dimen
        do j=1,dimen
        do k=1,dimen
        do l=1,dimen
        mechanical_hist_curentt(gauss_point_number,i,j,k,l)= &
        mechanical_hist_previus(gauss_point_number,i,j,k,l)*exp(-lambda_01*dtime)+ &
        ctens(i,j,k,l)*(k01*0.5d0)* &
        (  &
        exp(-lambda_01*dtime)*d_strain_p(k,l)+d_strain(k,l) &
        )
        end do;
        end do;
        end do;
        end do;



        sigma =  0.0d0 ;
        der   =  0.0d0 ;

        do i=1,dimen
        do j=1,dimen
        der(i)=der(i)+ktense(i,j)*electric_field(j)
        do k=1,dimen
        der(i)=der(i)+epz(i,j,k)*( strain_elastic(i,j) )
        sigma(i,j)=sigma(i,j) + epz(i,j,k)*electric_field(k)
        do l=1,dimen
 sigma(i,j)=sigma(i,j)+k00*ctens(i,j,k,l)*( strain_elastic(k,l) )  &
 + mechanical_hist_curentt(gauss_point_number,i,j,k,l)
        end do;
        end do;
        end do;
        end do;

    end subroutine stress_elect_displacement

    !     ================================================================
    !                       polarization swithing function
    !    this will find the polarization with the given state of material
    !     ================================================================
    subroutine remanent_polarization(d_electric_field,electric_field)
        implicit none
        real(iwp)::d_electric_field(3)
        real(iwp)::electric_field(3)
        real(iwp)::el_1,el_0
        real(iwp)::del_el_1
        real(iwp)::Dc
!        real(iwp)::electric_yield_e,electric_yield_d
!        real(iwp)::D_0,D_1! ,del_D_1
!        real(iwp)::Dp_0,Dp_1,del_Dp_1
        real(iwp)::Dp_1,electric_drive
        real(iwp)::delta_time
        real(iwp)::kappa_e

        kappa_e=2*(gamma2+gamma1)
        dc=ec*kappa_e


!        Dp_0 = hist_polarization_function(gauss_point_number,1)
!        Dp_1 = curn_polarization_function(gauss_point_number,1)

        delta_time=time(1)

!        D_1=curn_electric_displ(gauss_point_number,1)
!        D_0=hist_electric_displ(gauss_point_number,1)

!        el_1=curn_electric_field(gauss_point_number,1)
!        el_0=hist_electric_field(gauss_point_number,1)

        del_el_1=d_electric_field(1)
        el_1=electric_field(1)
        el_0=electric_field(1)-d_electric_field(1)
        !     ================================================================
        !       checking for electric yeald
        !     ================================================================
!        electric_yield_e=el_1**2-ec**2
!        electric_yield_d=D_1**2-dc**2
        !     ====================================inside the yield surface
       Dp_1=pr_sat*(tanh(abs(el_1)-3*ec)+1)*signum(el_1)*0.5d0

        if(abs(el_1)  .gt. 8*ec) then
            is_polarized=.true.
        endif
        !
        if(is_polarized)then
            electric_drive=el_1-3*ec*signum(del_el_1)
            Dp_1=pr_sat*tanh(abs(electric_drive))*signum(electric_drive)
        endif

        curn_polarization_function(gauss_point_number,1)=Dp_1

    end subroutine remanent_polarization


    !     ===============================================================
    !   forming history variables
    !     ===============================================================
    subroutine  form_history(ngauss)
        integer::ngauss
        !     ===================allocate history variables
        allocate( remanent_polarization_vector(ngauss,dimen) ,&
            hist_polarization_function(ngauss,dimen), &
            curn_polarization_function(ngauss,dimen), &
            hist_electric_field (ngauss,dimen),&
            curn_electric_field (ngauss,dimen),&
            hist_electric_displ (ngauss,dimen),&
            curn_electric_displ (ngauss,dimen))

        allocate(    mechanical_hist_curentt(ngauss,dimen,dimen,dimen,dimen), &
                     mechanical_hist_previus(ngauss,dimen,dimen,dimen,dimen), &
                     mechanical_hist_beforep(ngauss,dimen,dimen,dimen,dimen)  )
allocate(  &
mechanical_strain_curentt(ngauss,dimen,dimen) ,&
mechanical_strain_previus(ngauss,dimen,dimen) ,&
mechanical_strain_beforep(ngauss,dimen,dimen)  &
)

        mechanical_hist_curentt=0.0d0
        mechanical_hist_previus=0.0d0
        mechanical_hist_beforep=0.0d0

mechanical_strain_curentt=0.0d0;
mechanical_strain_previus=0.0d0;
mechanical_strain_beforep=0.0d0;

        remanent_polarization_vector=0.0d0

        hist_polarization_function=0.0d0
        curn_polarization_function=0.0d0

        hist_electric_field=0.0d0
        curn_electric_field=0.0d0

        hist_electric_displ=0.0d0
        curn_electric_displ=0.0d0
        is_polarized=.false.


    end subroutine  form_history

    !     ===============================================================
    !   forming history variables
    !     ===============================================================
    subroutine  update_history()


        hist_polarization_function=curn_polarization_function
        hist_electric_field =curn_electric_field
        hist_electric_displ =curn_electric_displ


        mechanical_hist_beforep=mechanical_hist_previus
        mechanical_hist_previus=mechanical_hist_curentt

        mechanical_strain_beforep=mechanical_strain_previus
        mechanical_strain_previus=mechanical_strain_curentt

    end subroutine  update_history

    !     ===============================================================
    !   Material properties
    !     ===============================================================
subroutine material_properties()
        implicit none
        !     ================================================================
        !      material variables
        !     ================================================================

        !     ================================================================
        integer::eye(dimen,dimen)
        integer::i ! ,j,k,l,m,n
        !     ================================================================
        !     ================================================================
        !     ================================================================
        lambda  =c12;
        mu      =(c11-c12)/2.0d0;
        alpha1  =2*c44+c12-c11;
        alpha2  =(c11+c33)/2.0d0-2.0d0*c44-c
        alpha3  =c13-c12;
        beta1   =-e31;
        beta2   =-e33+2.0d0*e15+e31;
        beta3   =-2.0d0*e15;
        gamma1  =-eps_11/2.0d0;
        gamma2  =(eps_11-eps_33)/2.0d0;

        eye = 0 ; do i = 1,dimen; eye(i,i) = 1;enddo

        !    1D formulation
        ctens=0;
        ctens(1,1,1,1) =lambda+2*mu+2*a**2*alpha3+2*a**4*alpha2+2*a**2*alpha1
        epz=0;
        epz(1,1,1)=a*(beta3+a**2*beta2+beta1)*dabs(pr)/pr_sat

        ktense=0.0d0  ;
        ktense(1,1)=(a**2*gamma2+gamma1)

end subroutine material_properties




    !     ================================================================
    !                       polarization swithing function
    !    this will find the polarization with the given state of material
    !     ================================================================
    subroutine direction_polarization(pr,a)
        real(iwp)::pr,a
        a=0.0d0
        if (abs(pr).gt.0.0d0)then
            a=abs(pr)/pr
        endif

    end subroutine direction_polarization


!    function derivitaive_dissipation(e);
!        real(iwp)::derivitaive_dissipation
!        real(iwp)::e
!
!        derivitaive_dissipation=0.0d0
!
!        if (abs(e).gt.abs(ec))then
!            derivitaive_dissipation=abs(e)*(abs(e)/ec-1)**m /e/eta
!        endif
!
!
!    end function derivitaive_dissipation;


!    function direction_vector(e);
!        real(iwp)::direction_vector
!        real(iwp)::e
!
!        direction_vector=0.0d0
!
!        if (abs(e).gt.0)then
!            direction_vector=e/abs(e)
!        endif
!    end function direction_vector


!    function dissipation(e);
!        real(iwp)::dissipation
!        real(iwp)::e
!        dissipation=ec*(abs(abs(e)/abs(ec)-1)+abs(e)/abs(ec)-1)**(m+1)*2**(-m-1)/(eta&
!            &*(m+1))
!        if (e.eq.0)dissipation=0
!    end function dissipation


end module material_behavior
