module sparsevector
    use csrmodule
    implicit none
    public Sparse_complex_vector
    public :: operator(*)

    type Sparse_complex_vector
        integer::N 
        complex(8),allocatable::val(:)
        integer,allocatable::indices(:)
        integer::numofdata
        real(8)::epsilon

        contains
        procedure::set => set_c
        procedure::print => print_csr 
        procedure::convert_to_dense 
        procedure::matmul => matmul_axy  !y = A*x: matmul(x,y)
        procedure::get => get_c
        procedure::clear
        procedure::norm => norm_sp
        procedure::dot => dot_sp
        procedure::add_vec => add_c_axby
        procedure::add_c_axby ! c = a*x + b*y
        procedure::add_a_axby ! a = a*x + b*y
        procedure::matmul_add
        procedure::reflesh
    end type


    interface Sparse_complex_vector
        module procedure::init_Sparse_complex_vector
        module procedure::init_Sparse_complex_vector_0
        module procedure::init_Sparse_complex_vector_dense
        module procedure::init_Sparse_complex_vector_dense_0
    end interface Sparse_complex_vector


    interface operator(*)
        module procedure mult_sp
    end interface

    contains

    type(Sparse_complex_vector) function init_Sparse_complex_vector_0(N) result(x)
        implicit none
        integer,intent(in)::N    
        x = Sparse_complex_vector(N,0d0)
        return
    end function

    type(Sparse_complex_vector) function init_Sparse_complex_vector(N,epsilon) result(x)
        implicit none
        integer,intent(in)::N    
        real(8),intent(in)::epsilon
        x%N = N
        x%epsilon = epsilon
        allocate(x%indices(N))
        allocate(x%val(N)) 
        x%val = 0d0
        x%indices = 0
        x%numofdata = 0
        return
    end function

    type(Sparse_complex_vector) function init_Sparse_complex_vector_dense(xdense,epsilon) result(x)
        implicit none 
        complex(8),intent(in)::xdense(:)
        real(8),intent(in)::epsilon
        integer::i
        integer::N
        N = ubound(xdense,1)

        x = Sparse_complex_vector(N,epsilon)
        do i=1,N
            if(abs(xdense(i)) > epsilon) then
                call x%set(xdense(i),i)
            end if
        end do
    

        return
    end function    



    type(Sparse_complex_vector) function init_Sparse_complex_vector_dense_0(xdense) result(x)
        implicit none 
        complex(8),intent(in)::xdense(:)
        x = Sparse_complex_vector(xdense,0d0)
        return
    end function

    real(8) function norm_sp(self) result(v)
        implicit none
        class(Sparse_complex_vector)::self
        integer::ip,i

        v = 0d0

        do ip=1,self%numofdata
            i = self%indices(ip)
            v = v + abs(self%val(i))**2
        end do        
        v = sqrt(v)

        return
    end function

    subroutine matmul_add(self,alpha,matA,a,beta,b) !self = alpha*matA*a + beta*b
        implicit none
        class(Sparse_complex_vector)::self
        type(Sparse_complex_vector),intent(in)::a,b
        type(CSRcomplex),intent(in)::matA
        complex(8),intent(in)::alpha,beta
        integer::ip,i
        logical::isb
        logical,allocatable::check(:)
        complex(8)::v,vi,vtmp
        type(Sparse_complex_vector)::ytmp
        integer::nonzeros,k,k2,j

        ytmp = Sparse_complex_vector(self%N,self%epsilon)
        !y = Sparse_complex_vector(self%N,self%epsilon)
        call self%clear()
        call ytmp%clear()


        nonzeros = a%numofdata
        !write(*,*) "nonzeros",nonzeros
        do k=1,nonzeros
            i = a%indices(k)
            vi = a%val(i)
            do k2=matA%row(i), matA%row(i+1)-1
                j = matA%col(k2)
                !write(*,*) "i,j ",i,j
                
                v = conjg(matA%val(k2))*vi
                vtmp = ytmp%get(j)
                call ytmp%set(vtmp + v,j)
            end do
        end do


        isb = .false.
        if(ytmp%numofdata < b%numofdata) then
            isb = .true.
        end if
        allocate(check(self%N))
        check = .false.

        if (isb) then
            do ip=1,b%numofdata
                i = b%indices(ip)
                check(i) = .true.
                v = alpha*ytmp%val(i) + beta*b%val(i)
                call self%set(v,i)
            end do
        else
            do ip=1,ytmp%numofdata
                i = ytmp%indices(ip)
                check(i) = .true.
                v = alpha*ytmp%val(i) + beta*b%val(i)
                call self%set(v,i)
            end do
        end if


        if (isb) then
            do ip=1,ytmp%numofdata
                i = ytmp%indices(ip)
                if(check(i) .neqv. .true.) then
                    v = alpha*ytmp%val(i) + beta*b%val(i)
                    call self%set(v,i)
                end if
            end do
        else
            do ip=1,b%numofdata
                i = b%indices(ip)
                if(check(i) .neqv. .true.) then
                    v = alpha*ytmp%val(i) + beta*b%val(i)
                    call self%set(v,i)
                end if
            end do
        end if




    end subroutine

    subroutine add_c_axby(self,alpha,a,beta,b) !c = alpha*a + beta*b
        implicit none
        class(Sparse_complex_vector)::self
        type(Sparse_complex_vector),intent(in)::a,b
        complex(8),intent(in)::alpha,beta
        integer::ip,i
        logical::isb
        logical,allocatable::check(:)
        complex(8)::v

        isb = .false.
        if(a%numofdata < b%numofdata) then
            isb = .true.
        end if
        allocate(check(self%N))
        check = .false.

        call self%clear()
        if (isb) then
            do ip=1,b%numofdata
                i = b%indices(ip)
                check(i) = .true.
                v = alpha*a%val(i) + beta*b%val(i)
                call self%set(v,i)
            end do
        else
            do ip=1,a%numofdata
                i = a%indices(ip)
                check(i) = .true.
                v = alpha*a%val(i) + beta*b%val(i)
                call self%set(v,i)
            end do
        end if


        if (isb) then
            do ip=1,a%numofdata
                i = a%indices(ip)
                if(check(i) .neqv. .true.) then
                    v = alpha*a%val(i) + beta*b%val(i)
                    call self%set(v,i)
                end if
            end do
        else
            do ip=1,b%numofdata
                i = b%indices(ip)
                if(check(i) .neqv. .true.) then
                    v = alpha*a%val(i) + beta*b%val(i)
                    call self%set(v,i)
                end if
            end do
        end if





    end subroutine

    subroutine add_a_axby(self,alpha,beta,b) !a = alpha*a + beta*b
        implicit none
        class(Sparse_complex_vector)::self
        type(Sparse_complex_vector),intent(in)::b
        complex(8),intent(in)::alpha,beta
        integer::ip,i
        logical::isb
        logical,allocatable::check(:)
        complex(8)::v

        isb = .false.
        if(self%numofdata < b%numofdata) then
            isb = .true.
        end if
        allocate(check(self%N))
        check = .false.


        if (isb) then
            do ip=1,b%numofdata
                i = b%indices(ip)
                
                check(i) = .true.
                v = alpha*self%val(i) + beta*b%val(i)
                call self%set(v,i)
            end do
        else
            do ip=1,self%numofdata
                i = self%indices(ip)
                check(i) = .true.
                v = alpha*self%val(i) + beta*b%val(i)
                call self%set(v,i)
            end do
        end if


        if (isb) then
            do ip=1,self%numofdata
                i = self%indices(ip)
                
                if(check(i) .neqv. .true.) then
                    v = alpha*self%val(i) + beta*b%val(i)
                    call self%set(v,i)
                end if
            end do
        else
            do ip=1,b%numofdata
                i = b%indices(ip)
                if(check(i) .neqv. .true.) then
                    v = alpha*self%val(i) + beta*b%val(i)
                    call self%set(v,i)
                end if
            end do
        end if





    end subroutine


    complex(8) function dot_sp(self,b) result(v)
        implicit none
        class(Sparse_complex_vector)::self
        type(Sparse_complex_vector),intent(in)::b
        integer::ip,i
        logical::isb
        logical,allocatable::check(:)

        isb = .false.
        if(self%numofdata < b%numofdata) then
            isb = .true.
        end if
        allocate(check(self%N))
        check = .false.

        v = 0d0

        if (isb) then
            do ip=1,b%numofdata
                i = b%indices(ip)
                check(i) = .true.
                v = v + conjg(self%val(i))*b%val(i)

            end do
        else
            do ip=1,self%numofdata
                i = self%indices(ip)
                check(i) = .true.
                v = v + conjg(self%val(i))*b%val(i)
            end do
        end if

        if (isb) then
            do ip=1,self%numofdata
                i = self%indices(ip)
                if(check(i) .neqv. .true.) then
                    v = v + conjg(self%val(i))*b%val(i)
                end if
            end do
        else
            do ip=1,b%numofdata
                i = b%indices(ip)
                if(check(i) .neqv. .true.) then
                    v = v + conjg(self%val(i))*b%val(i)
                end if
            end do
        end if
        


    end

    function mult_sp(A,x) result(y)
        implicit none
        class(CSRcomplex),intent(in)::A
        type(Sparse_complex_vector),intent(in)::x
        type(Sparse_complex_vector)::y

        y = Sparse_complex_vector(x%N,x%epsilon)
        call y%clear()

        call y%matmul(A,x)


    end function

    subroutine matmul_axy(self,A,x) !self = A*x 
        implicit none
        type(CSRcomplex),intent(in)::A
        class(Sparse_complex_vector),intent(in)::self
        type(Sparse_complex_vector),intent(in)::x

        integer::i,j,k,k2
        integer::nonzeros
        type(Sparse_complex_vector)::ytmp
        complex(8)::v,vi,vtmp
        ytmp = Sparse_complex_vector(self%N,self%epsilon)
        !y = Sparse_complex_vector(self%N,self%epsilon)
        call self%clear()
        call ytmp%clear()


        nonzeros = x%numofdata
        !write(*,*) "nonzeros",nonzeros
        do k=1,nonzeros
            i = x%indices(k)
            vi = x%val(i)
            do k2=A%row(i), A%row(i+1)-1
                j = A%col(k2)
                !write(*,*) "i,j ",i,j
                
                v = conjg(A%val(k2))*vi
                vtmp = ytmp%get(j)
                call ytmp%set(vtmp + v,j)
            end do
        end do

        !write(*,*) ytmp%indices

        do k=1,ytmp%numofdata
            i = ytmp%indices(k)
            call self%set(ytmp%val(i),i)
        end do
            
    end subroutine

    subroutine clear(self)
        implicit none
        class(Sparse_complex_vector)::self
        integer::i,k
        !self%val(1:self%numofdata) = 0d0
        do k=1,self%numofdata
            i = self%indices(k)
            self%val(i) = 0d0
            self%indices(k) = 0 
        end do
        !self%indices(1:self%numofdata) = 0
        self%numofdata = 0

    end subroutine

    subroutine convert_to_dense(self,xdense)
        implicit none
        class(Sparse_complex_vector)::self
        complex(8),intent(out)::xdense(:)
        integer::i,j

        if(self%N .ne. ubound(xdense,1)) then
            write(*,*) "The length check is failed. The length of output dense vector is wrong"
            stop
        end if
        xdense = 0d0


        do i=1,self%numofdata
            j = self%indices(i)
            xdense(j) = self%val(j)
        end do

        return
    end subroutine convert_to_dense


    subroutine print_csr(self)
        implicit none
        class(Sparse_complex_vector)::self
        integer::i,j


        do i=1,self%numofdata
            j = self%indices(i)
            write(*,*) j,self%val(j)
        end do

        return
    end subroutine print_csr

    complex(8) function get_c(self,j) result(v)
        implicit none
        class(Sparse_complex_vector)::self
        integer,intent(in)::j

        v = self%val(j)
        return

    end function

    subroutine reflesh(self)
        implicit none
        class(Sparse_complex_vector)::self
        integer::ip,i
        type(Sparse_complex_vector)::tmp
        complex(8)::v

        tmp = Sparse_complex_vector(self%N,self%epsilon)

        do ip=1,self%numofdata
            i = self%indices(ip)
            v = self%val(i)
            if(abs(v) > self%epsilon) then
                call tmp%set(v,i)
            end if
        end do

        call self%clear()
        do ip=1,tmp%numofdata
            i = tmp%indices(ip)
            v = tmp%val(i)
            call self%set(v,i)
        end do

    end subroutine


    subroutine set_c(self,v,j)
        implicit none
        class(Sparse_complex_vector)::self
        complex(8),intent(in)::v
        integer,intent(in)::j
        integer::ip,i
        complex(8)::vi
        type(Sparse_complex_vector)::tmp
        complex(8)::vtmp
        vi = self%val(j)


        if(abs(vi) ==  0d0) then
            if(abs(v) > self%epsilon) then
                self%numofdata = self%numofdata + 1
                self%indices(self%numofdata) = j
                self%val(j) = v
            end if
        else
            self%val(j) = v
            !if(abs(v) > self%epsilon) then
            !    self%val(j) = v
            !end if
        end if
        return


        if(abs(v) .le. self%epsilon) then
            if(abs(vi) .ne.  0d0) then
                self%val(j) = 1d-18
            end if
            if(abs(vi) > self%epsilon) then
                !write(*,*) "nonzero",self%numofdata,v
                !tmp = Sparse_complex_vector(self%N,self%epsilon)
                !self%val(j) = 0d0
                !do ip=1,self%numofdata
                !    i = self%indices(ip)
                !    vtmp = self%val(i)
                !    if(abs(vtmp) > self%epsilon) then
                !        call tmp%set(vtmp,i)
                !    end if
                !end do
                !call self%clear()
                !do ip=1,tmp%numofdata
                !    i = tmp%indices(ip)
                !    vtmp = tmp%val(i)
                !    call self%set(vtmp,i)
                !end do
                !self%val(j) = self%epsilon+1d-16
                self%val(j) = 1d-18
                !write(*,*) "nonzero after",self%numofdata
            end if
            !return
        end if
        

        !write(*,*) v,j,self%numofdata
        if(abs(vi) .ne.  0d0) then
        !if(abs(vi) > self%epsilon) then
            self%val(j) = v
        else
            self%numofdata = self%numofdata + 1
            self%indices(self%numofdata) = j
            self%val(j) = v
        end if
        return
    end subroutine set_c






end module sparsevector