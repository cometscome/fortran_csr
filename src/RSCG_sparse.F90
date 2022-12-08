module RSCG_sparse

    implicit none

    contains

    subroutine calculate_Green_function(H,ii,jj,vec_sigma,vec_theta,epsilon,eps,maximumsteps)
        use sparsevector
        use csrmodule
        implicit none
        type(CSRcomplex),intent(in)::H
        integer,intent(in)::ii,jj
        complex(8),intent(in)::vec_sigma(:)
        complex(8),intent(out)::vec_theta(:)
        real(8),intent(in)::epsilon,eps
        integer,intent(in)::maximumsteps
        integer::M,N
        type(Sparse_complex_vector)::b,x,p,Ap,r
        complex(8),parameter::ci=(0d0,1d0)
        complex(8)::alpha_m,beta_m,Lsigma
        complex(8)::rnorm,alpha,beta,alpha_kj,beta_kj
        complex(8),allocatable,dimension(:)::rho_k,rho_km,rho_kp,vec_Pi
        integer::k,j
        real(8)::hi2
        complex(8)::denom
        complex(8),parameter::cone = (1d0,0d0)

        M = ubound(vec_sigma,1)
        N = H%N
        b = Sparse_complex_vector(N,epsilon)
        x = Sparse_complex_vector(N,epsilon)
        p = Sparse_complex_vector(N,epsilon)
        Ap = Sparse_complex_vector(N,epsilon)
        r = Sparse_complex_vector(N,epsilon)
        !call H%print()

        call b%set(cone,jj)
        call r%set(cone,jj) !r = b
        call p%set(cone,jj) !p = b

        alpha_m = 1d0
        beta_m = 0d0

        Lsigma = b%get(ii)
        allocate(vec_Pi(1:M))
        vec_Pi = Lsigma
        vec_theta = 0d0
        allocate(rho_k(1:M))
        rho_k = 1d0
        allocate(rho_km(1:M))
        rho_km = 1d0
        allocate(rho_kp(1:M))
        rho_kp = 1d0

        do k=0,maximumsteps
            call Ap%matmul(H,p) !Ap = H*p
            rnorm = r%norm()**2

            alpha = -rnorm/p%dot(Ap)


            call x%add_a_axby(cone,alpha,p)
            call r%add_a_axby(cone,alpha,Ap)

            beta = r%norm()**2/rnorm
            !write(*,*) "p",beta
            call p%add_a_axby(beta,cone,r) !p = beta*p + r
            Lsigma =  r%get(ii)


            do j=1,M
                if (abs(rho_k(j)) > 1d-12) then
                    denom = rho_km(j)*alpha_m*(1d0+alpha*vec_sigma(j)) +alpha*beta_m*(rho_km(j)-rho_k(j))
                    rho_kp(j) = rho_k(j)*rho_km(j)*alpha_m/denom
                    alpha_kj = alpha*rho_kp(j)/rho_k(j)
                    vec_theta(j) = vec_theta(j) + alpha_kj*vec_Pi(j)
                    beta_kj = ((rho_kp(j)/rho_k(j))**2)*beta
                    vec_Pi(j) = rho_kp(j)*Lsigma+beta_kj*vec_Pi(j)
                    rho_km(j) = rho_k(j)
                    rho_k(j) = rho_kp(j)
                end if
            end do
            alpha_m = alpha
            beta_m = beta
            hi2 = dble(rnorm)*maxval(abs(rho_k)) 
            write(*,*) k,hi2,x%numofdata

            if (hi2 < eps) then
                write(*,*) k,hi2,x%numofdata
                return 
            end if

            if(mod(k,10)==0) then
                call x%reflesh()
                call r%reflesh()
                call p%reflesh()
            end if


        end do


        




    end subroutine

    subroutine calculate_Green_function_csr(H,ii,jj,vec_sigma,vec_theta,epsilon,eps,maximumsteps)
        use csrmodule
        implicit none
        type(CSRcomplex),intent(in)::H
        integer,intent(in)::ii,jj
        complex(8),intent(in)::vec_sigma(:)
        complex(8),intent(out)::vec_theta(:)
        real(8),intent(in)::epsilon,eps
        integer,intent(in)::maximumsteps
        integer::M,N
        complex(8),allocatable,dimension(:)::b,x,p,Ap,r
        complex(8),parameter::ci=(0d0,1d0)
        complex(8)::alpha_m,beta_m,Lsigma
        complex(8)::rnorm,alpha,beta,alpha_kj,beta_kj
        complex(8),allocatable,dimension(:)::rho_k,rho_km,rho_kp,vec_Pi
        integer::k,j,i
        real(8)::hi2
        complex(8)::denom
        complex(8),parameter::cone = (1d0,0d0)

        M = ubound(vec_sigma,1)
        N = H%N
        allocate(b(N),x(N),p(N),Ap(N),r(N))
        b = 0d0
        r = 0d0
        p = 0d0
        Ap = 0d0
        x = 0d0
        b(jj) = cone
        r(jj) = cone
        p(jj) = cone

        alpha_m = 1d0
        beta_m = 0d0

        Lsigma = b(ii)
        allocate(vec_Pi(1:M))
        vec_Pi = Lsigma
        vec_theta = 0d0
        allocate(rho_k(1:M))
        rho_k = 1d0
        allocate(rho_km(1:M))
        rho_km = 1d0
        allocate(rho_kp(1:M))
        rho_kp = 1d0

        do k=0,maximumsteps
            Ap = H*p!matmul(H,p)
            rnorm = dot_product(r,r)


            alpha = -rnorm/dot_product(p,Ap)
            x = x + alpha*p
            r = r + alpha*Ap


            beta = dot_product(r,r)/rnorm
            p = beta*p + r
            
            Lsigma =  r(ii)


            do j=1,M
                if (abs(rho_k(j)) > 1d-12) then
                    denom = rho_km(j)*alpha_m*(1d0+alpha*vec_sigma(j)) +alpha*beta_m*(rho_km(j)-rho_k(j))
                    rho_kp(j) = rho_k(j)*rho_km(j)*alpha_m/denom
                    alpha_kj = alpha*rho_kp(j)/rho_k(j)
                    vec_theta(j) = vec_theta(j) + alpha_kj*vec_Pi(j)
                    beta_kj = ((rho_kp(j)/rho_k(j))**2)*beta
                    vec_Pi(j) = rho_kp(j)*Lsigma+beta_kj*vec_Pi(j)
                    rho_km(j) = rho_k(j)
                    rho_k(j) = rho_kp(j)
                end if
            end do
            alpha_m = alpha
            beta_m = beta
            hi2 = dble(rnorm)*maxval(abs(rho_k)) 
            write(*,*) k,hi2

            if (hi2 < eps) then
                write(*,*) k,hi2
                return 
            end if


        end do


        




    end subroutine

    subroutine calculate_Green_function_dense(H,ii,jj,vec_sigma,vec_theta,epsilon,eps,maximumsteps)
        implicit none
        complex(8),intent(in)::H(:,:)
        integer,intent(in)::ii,jj
        complex(8),intent(in)::vec_sigma(:)
        complex(8),intent(out)::vec_theta(:)
        real(8),intent(in)::epsilon,eps
        integer,intent(in)::maximumsteps
        integer::M,N
        complex(8),allocatable,dimension(:)::b,x,p,Ap,r
        complex(8),parameter::ci=(0d0,1d0)
        complex(8)::alpha_m,beta_m,Lsigma
        complex(8)::rnorm,alpha,beta,alpha_kj,beta_kj
        complex(8),allocatable,dimension(:)::rho_k,rho_km,rho_kp,vec_Pi
        integer::k,j
        real(8)::hi2
        complex(8)::denom
        complex(8),parameter::cone = (1d0,0d0)

        M = ubound(vec_sigma,1)
        N = ubound(H,1)
        allocate(b(N),x(N),p(N),Ap(N),r(N))
        b = 0d0
        r = 0d0
        p = 0d0
        Ap = 0d0
        x = 0d0
        b(jj) = cone
        r(jj) = cone
        p(jj) = cone

        alpha_m = 1d0
        beta_m = 0d0

        Lsigma = b(ii)
        allocate(vec_Pi(1:M))
        vec_Pi = Lsigma
        vec_theta = 0d0
        allocate(rho_k(1:M))
        rho_k = 1d0
        allocate(rho_km(1:M))
        rho_km = 1d0
        allocate(rho_kp(1:M))
        rho_kp = 1d0

        do k=0,maximumsteps
            Ap = matmul(H,p)
            rnorm = dot_product(r,r)

            alpha = -rnorm/dot_product(p,Ap)
            x = x + alpha*p
            r = r + alpha*Ap

            beta = dot_product(r,r)/rnorm
            p = beta*p + r
            
            Lsigma =  r(ii)


            do j=1,M
                if (abs(rho_k(j)) > 1d-12) then
                    denom = rho_km(j)*alpha_m*(1d0+alpha*vec_sigma(j)) +alpha*beta_m*(rho_km(j)-rho_k(j))
                    rho_kp(j) = rho_k(j)*rho_km(j)*alpha_m/denom
                    alpha_kj = alpha*rho_kp(j)/rho_k(j)
                    vec_theta(j) = vec_theta(j) + alpha_kj*vec_Pi(j)
                    beta_kj = ((rho_kp(j)/rho_k(j))**2)*beta
                    vec_Pi(j) = rho_kp(j)*Lsigma+beta_kj*vec_Pi(j)
                    rho_km(j) = rho_k(j)
                    rho_k(j) = rho_kp(j)
                end if
            end do
            alpha_m = alpha
            beta_m = beta
            hi2 = dble(rnorm)*maxval(abs(rho_k)) 
            write(*,*) k,hi2

            if (hi2 < eps) then
                write(*,*) k,hi2
                return 
            end if


        end do


        




    end subroutine



end module