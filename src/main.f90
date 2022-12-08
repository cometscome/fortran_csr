program main
    use csrmodule
    call rscgtest()
    !call vectortest()
    !call matrixtest()
end program main

subroutine rscgtest()
    use sparsevector
    use csrmodule
    use RSCG_sparse

    implicit none
    type(Sparse_complex_vector)::x_s,y_s,ytemp2_s,ytemp3_s
    type(CSRcomplex)::H
    integer::N
    integer::Nx,Ny
    complex(8)::v,alpha,beta
    real(8)::mu,t
    complex(8),parameter::ci = (0d0,1d0)
    integer::ii,jj
    complex(8),allocatable::vec_sigma(:),vec_theta(:)
    real(8)::epsilon,eps
    integer::maximumsteps
    integer::M,i,ix,iy
    real(8)::emax,emin,de,eta
    real(8),parameter::pi = 3.141592653589793238d0
    real(8)::d
    integer::j
    complex(8),allocatable::Hdense(:,:)
    real(8)::time_begin,time_end
    


    v = -1d0
    Nx = 150
    Ny = 150
    N = Nx*Ny
    t = -1d0
    mu = -1.5d0
    d = 0.5d0



    call make_2dcsrmatrix(H,Nx,Ny,mu,t,d)
    !call H%print()

    !call make_csrmatrix(H,N,mu,t,d)
    !allocate(Hdense(2*N,2*N))
    !call H%convert_to_dense(Hdense)
    !call make_matrix(H,N,mu,t,d,Hdense)
    !do i=1,2*N
    !    do j=1,2*N
    !        if(Hdense(i,j) .ne. 0d0) then
    !            write(*,*) i,j,Hdense(i,j)
    !        end if
    !    end do
    !end do
    !stop

    ix = Nx/2
    iy = Ny/2
    ii = (iy-1)*Nx + ix
    jj = ii

    maximumsteps = 10000
    eps = 1d-7
    epsilon = 1d-12

    M = 1000
    allocate(vec_sigma(M),vec_theta(M))
    emax = 6d0
    emin = -emax
    eta = 0.05d0
    de = (emax-emin)/dble(M-1)
    do i=1,M
        vec_sigma(i) =dble(i-1)*de + emin + ci*eta
        !vec_sigma(i) = dble(i-1)*ci !dble(i-1)*de + emin + ci*eta
    end do

    write(*,*) "calculating Green's function..."
    call cpu_time(time_begin)
    call calculate_Green_function_csr(H,ii,jj,vec_sigma,vec_theta,epsilon,eps,maximumsteps)
    call cpu_time(time_end)
    write(*,*) time_end-time_begin,"[sec]"
    write(*,*) "calculating Green's function..."
    call cpu_time(time_begin)
    call calculate_Green_function(H,ii,jj,vec_sigma,vec_theta,epsilon,eps,maximumsteps)
    call cpu_time(time_end)
    write(*,*) time_end-time_begin,"[sec]"

    stop
    !call calculate_Green_function_dense(Hdense,ii,jj,vec_sigma,vec_theta,epsilon,eps,maximumsteps)
    !open(11,file="ldos.txt")
    open(11,file="mat.txt")
    do i=1,M
        write(11,*) dble(vec_sigma(i)),(-1d0/pi)*dimag(vec_theta(i))
        !write(*,*) dble(vec_sigma(i)),(-1d0/pi)*(vec_theta(i))
    end do
    close(11)

end subroutine

subroutine vectortest()
    use csrvector
    use sparsevector
    implicit none
    type(CSRcomplex_vector)::x,y,ytemp2
    type(Sparse_complex_vector)::x_s,y_s,ytemp2_s,ytemp3_s
    integer::N
    complex(8)::v,alpha,beta
    real(8)::mu,t
    integer::i
    complex(8),allocatable::xdense(:)
    complex(8),allocatable::ytemp(:)
    complex(8),allocatable::xdense2(:)
    complex(8),allocatable::xdense3(:)
    type(CSRcomplex)::H
    complex(8),allocatable::Hdense(:,:)
    complex(8),parameter::ci = (0d0,1d0)
    real(8)::d
    
    

    v = -1d0


    N = 10

    x = CSRcomplex_vector(2*N)
    call x%set(v,5)
    call x%set(ci,8)
    call x%set(0.2*ci,1)

    call x%print()

    x_s = Sparse_complex_vector(2*N)
    call x_s%set(v,5)
    call x_s%set(ci,8)
    call x_s%set(0.2*ci,1)

    call x_s%print()


    write(*,*) "dense"

    allocate(xdense(2*N))
    call x%convert_to_dense(xdense)
    do i=1,N
        if (xdense(i) .ne. 0d0) then
            write(*,*) i,xdense(i)
        end if
    end do

    write(*,*) "dense"

    call x_s%convert_to_dense(xdense)
    do i=1,N
        if (xdense(i) .ne. 0d0) then
            write(*,*) i,xdense(i)
        end if
    end do


    y = CSRcomplex_vector(xdense)
    write(*,*) "y"
    call y%print()

    alpha = y%get(8)
    write(*,*) alpha

    y_s = Sparse_complex_vector(xdense)
    write(*,*) "y_s"
    call y_s%print()

    alpha = y_s%get(8)
    write(*,*) alpha
    !alpha = y%get(9)
    !write(*,*) alpha
    

    allocate(ytemp(1:2*N))

    t = -1d0
    mu = -1.5d0
    d = 0.5d0

    allocate(Hdense(2*N,2*N))
    call make_matrix(H,N,mu,t,d,Hdense)
    ytemp = H*xdense
    

    xdense2 = matmul(Hdense,xdense)
    ytemp2 = H*y
    ytemp2_s = H*y_s

    write(*,*) y_s%dot(ytemp2_s)
    write(*,*) ytemp2_s%dot(ytemp2_s)
    write(*,*) ytemp2_s%norm()**2

    xdense3 = 0.2*matmul(Hdense,xdense2) + 0.1*xdense

    ytemp3_s = Sparse_complex_vector(2*N)
    alpha = 0.2d0
    beta = 0.1d0
    call ytemp3_s%matmul_add(alpha,H,ytemp2_s,beta,y_s)
    do i=1,2*N
        write(*,*) i,xdense3(i),ytemp3_s%get(i)
    end do
    


    stop
    do i=1,2*N
        write(*,*) i,ytemp(i),ytemp2%get(i),ytemp2_s%get(i)
    end do

    

    xdense = H*ytemp
    xdense3 = matmul(Hdense,xdense2)

    y = H*ytemp2
    y_s = H*ytemp2_s

    do i=1,2*N
        write(*,*) i,xdense(i),y%get(i),y_s%get(i)
    end do

    ytemp = H*xdense
    xdense2 = matmul(Hdense,xdense3)
    ytemp2 = H*y
    ytemp2_s = H*y_s

    do i=1,2*N
        write(*,*) i,ytemp(i),ytemp2%get(i),ytemp2_s%get(i)
    end do







end subroutine

subroutine make_2dcsrmatrix(H,Nx,Ny,mu,t,d)
    use csrmodule
    implicit none
    type(CSRcomplex)::H
    integer::i,j,ix,iy,jx,jy
    integer,intent(in)::Nx,Ny
    real(8)::mu,t
    complex(8)::v
    real(8)::d
    integer::N

    N = Nx*Ny
    H = CSRcomplex(2*N)
    do ix=1,Nx
        do iy=1,Ny
            i = (iy-1)*Nx + ix

            jx = ix
            jy = iy
            j = (jy-1)*Nx + jx
            v = -mu
            call H%set(v,i,j)
            call H%set(-v,i+N,j+N)


            jx = ix +1
            if(jx > Nx) jx = jx - Nx
            jy = iy
            j = (jy-1)*Nx + jx
            v = t
            call H%set(v,i,j)
            call H%set(-v,i+N,j+N)

            jx = ix -1
            if(jx < 1) jx = jx + Nx
            jy = iy
            j = (jy-1)*Nx + jx
            v = t
            call H%set(v,i,j)
            call H%set(-v,i+N,j+N)

            jx = ix 
            jy = iy +1
            if(jy >  Ny) jy = jy - Ny
            j = (jy-1)*Nx + jx
            v = t
            call H%set(v,i,j)
            call H%set(-v,i+N,j+N)

            jx = ix 
            jy = iy-1
            if(jy <  1) jy = jy + Ny
            j = (jy-1)*Nx + jx
            v = t
            call H%set(v,i,j)
            call H%set(-v,i+N,j+N)

            j = i + N
            v = d
            call H%set(v,i,j)
            call H%set(v,j,i)

        end do
    end do



end subroutine

subroutine make_csrmatrix(H,N,mu,t,d)
    use csrmodule
    implicit none
    type(CSRcomplex)::H
    integer::i,j
    integer,intent(in)::N
    real(8)::mu,t
    complex(8)::v
    real(8)::d



    

    H = CSRcomplex(2*N)
    do i = 1,N
        j = i
        v = -mu

        call H%set(v,i,j)
        call H%set(-v,i+N,j+N)

        v = t
        j = i+1
        if(j > N) j = j -N

        if (j > 0 .and. j < N+1) then
            call H%set(v,i,j)
            call H%set(-v,i+N,j+N)
        end if

        j = i-1
        if(j < 1) j = j +N
        v =t

        if (j > 0 .and. j < N+1) then
            call H%set(v,i,j)
            call H%set(-v,i+N,j+N)
        end if

        j = i+N
        v = d
        call H%set(v,i,j)
        call H%set(v,j,i)

    end do




end subroutine

subroutine make_matrix(H,N,mu,t,d,Hdense)
    use csrmodule
    implicit none
    type(CSRcomplex)::H
    integer::i,j
    integer,intent(in)::N
    real(8)::mu,t
    real(8)::d
    complex(8)::v
    complex(8)::Hdense(2*N,2*N)



    

    H = CSRcomplex(2*N)
    do i = 1,N
        j = i
        v = -mu

        call H%set(v,i,j)
        Hdense(i,j) = v
        call H%set(-v,i+N,j+N)
        Hdense(i+N,j+N) = -v

        v = t
        j = i+1
        if(j > N) j = j -N
        if (j > 0 .and. j < N+1) then
            call H%set(v,i,j)
            Hdense(i,j) = v

            call H%set(-v,i+N,j+N)
            Hdense(i+N,j+N) = -v

        end if

        j = i-1
        if(j < 1) j = j +N
        v =t

        if (j > 0 .and. j < N+1) then
            call H%set(v,i,j)
            Hdense(i,j) = v
            call H%set(-v,i+N,j+N)
            Hdense(i+N,j+N) = -v

        end if

        j = i+N
        v = d
        call H%set(v,i,j)
        Hdense(i,j) = v
        call H%set(v,j,i)
        Hdense(j,i) = v

    end do




end subroutine

subroutine matrixtest()
    use csrmodule
    implicit none
    type(CSRcomplex)::H
    integer::N
    complex(8)::v,alpha,beta
    integer::i,j
    real(8)::mu,t
    complex(8),allocatable::x(:)
    complex(8),allocatable::y(:)
    complex(8),allocatable::ytemp(:)
    complex(8),allocatable::z(:)
    complex(8),allocatable::Hdense(:,:)
    complex(8),allocatable::ytemp2(:)

    N = 10
    v = -1d0
    t = -1d0
    mu = -1.5d0
    allocate(Hdense(2*N,2*N))

    H = CSRcomplex(2*N)
    do i = 1,N
        j = i
        v = -mu

        call H%set(v,i,j)
        Hdense(i,j) = v
        call H%set(-v,i+N,j+N)
        Hdense(i+N,j+N) = -v

        v = t
        j = i+1
        if (j > 0 .and. j < N+1) then
            call H%set(v,i,j)
            Hdense(i,j) = v
            call H%set(v,j,i)
            Hdense(j,i) = v

            call H%set(-v,i+N,j+N)
            Hdense(i+N,j+N) = -v
            call H%set(-v,j+N,i+N)
            Hdense(j+N,i+N) = -v
        end if

        j = i-1
        v =t

        if (j > 0 .and. j < N+1) then
            call H%set(v,i,j)
            Hdense(i,j) = v
            call H%set(v,j,i)
            Hdense(j,i) = v

            call H%set(-v,i+N,j+N)
            Hdense(i+N,j+N) = -v
            call H%set(-v,j+N,i+N)
            Hdense(j+N,i+N) = -v
        end if

        j = i+N
        v = 0.5d0
        call H%set(v,i,j)
        Hdense(i,j) = v
        call H%set(v,j,i)
        Hdense(j,i) = v

    end do

    call H%print()

    do i=1,N
        j = i+N
        v = 0.6d0
        call H%update(v,i,j)
        Hdense(i,j) = v
        call H%update(v,j,i)
        Hdense(j,i) = v


    end do
    call H%print()

    allocate(x(1:2*N))
    x = 0d0
    allocate(y(1:2*N))
    allocate(ytemp(1:2*N))
    allocate(ytemp2(1:2*N))
    y = 0d0
    ytemp = 0d0
    allocate(z(1:2*N))
    z(5) = 10d0
    x(3) = 1d0
    call H%matmul(x,y)
    !write(*,*) "matmul(x,y) ", y


    ytemp = H*x
    write(*,*) "H*x", ytemp

    ytemp2 = matmul(Hdense,x)
    do i=1,2*N
        write(*,*) i,ytemp(i),ytemp2(i)
    end do
    write(*,*) "diff = ",dot_product(ytemp-ytemp2,ytemp-ytemp2)/N




    alpha = 2d0
    beta = 3d0
    call H%matmul2(x,y,z,alpha,beta)
    !write(*,*) "y = alpha*A*x+beta*z",y
end subroutine
