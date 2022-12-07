program main
    use csrmodule
    call vectortest()
    !call matrixtest()
end program main

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

    allocate(Hdense(2*N,2*N))
    call make_matrix(H,N,mu,t,v,Hdense)
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

subroutine make_matrix(H,N,mu,t,v,Hdense)
    use csrmodule
    implicit none
    type(CSRcomplex)::H
    integer::i,j
    integer,intent(in)::N
    real(8)::mu,t
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

    !call H%print()

    do i=1,N
        j = i+N
        v = 0.6d0
        call H%update(v,i,j)
        Hdense(i,j) = v
        call H%update(v,j,i)
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
