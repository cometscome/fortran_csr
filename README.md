# Fortran library for CSR format
This module can construct a sparce matrix with CSR format. The example is written in main.f90.
Please see CMakeList.txt if you want to know how to link it.


How to test it.
```
mkdir build
cd build
cmake ..
./test
```

# example 

```fortran
program main
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


end program main
```