module csrvector
    use csrmodule
    implicit none
    public CSRcomplex_vector
    !public :: operator(*)

    type CSRcomplex_vector
        private
        integer::N 
        complex(8),allocatable::val(:)
        integer,allocatable::indices(:)
        integer::numofdata
        real(8)::epsilon

        contains
        procedure::set => set_c
        procedure::print => print_csr 
        procedure::convert_to_dense 
    end type


    interface CSRcomplex_vector
        module procedure::init_CSRcomplex_vector
        module procedure::init_CSRcomplex_vector_0
        module procedure::init_CSRcomplex_dense
        module procedure::init_CSRcomplex_dense_0
    end interface CSRcomplex_vector

    interface insert_element
        module procedure::insert_c
        module procedure::insert_int
    end interface

    contains

    type(CSRcomplex_vector) function init_CSRcomplex_vector_0(N) result(x)
        implicit none
        integer,intent(in)::N    
        x = CSRcomplex_vector(N,0d0)
        return
    end function

    type(CSRcomplex_vector) function init_CSRcomplex_vector(N,epsilon) result(x)
        implicit none
        integer,intent(in)::N    
        real(8),intent(in)::epsilon
        x%N = N
        x%epsilon = 0d0
        allocate(x%indices(N))
        allocate(x%val(N)) 
        x%val = 0d0
        x%indices = 0
        x%numofdata = 0
        return
    end function

    type(CSRcomplex_vector) function init_CSRcomplex_dense(xdense,epsilon) result(x)
        implicit none 
        complex(8),intent(in)::xdense(:)
        real(8),intent(in)::epsilon
        integer::i
        integer::N
        N = ubound(xdense,1)

        x = CSRcomplex_vector(N,epsilon)
        do i=1,N
            if(abs(xdense(i)) .ge. epsilon) then
                call x%set(xdense(i),i)
            end if
        end do
    

        return
    end function    

    type(CSRcomplex_vector) function init_CSRcomplex_dense_0(xdense) result(x)
        implicit none 
        complex(8),intent(in)::xdense(:)
        x = CSRcomplex_vector(xdense,0d0)
        return
    end function

    subroutine convert_to_dense(self,xdense)
        implicit none
        class(CSRcomplex_vector)::self
        complex(8),intent(out)::xdense(:)
        integer::i,j

        if(self%N .ne. ubound(xdense,1)) then
            write(*,*) "The length check is failed. The length of output dense vector is wrong"
            stop
        end if
        xdense = 0d0


        do i=1,self%numofdata
            j = self%indices(i)
            xdense(j) = self%val(i)
        end do

        return
    end subroutine convert_to_dense


    subroutine print_csr(self)
        implicit none
        class(CSRcomplex_vector)::self
        integer::i,j


        do i=1,self%numofdata
            j = self%indices(i)
            write(*,*) j,self%val(i)
        end do

        return
    end subroutine print_csr


    subroutine set_c(self,v,j)
        implicit none
        class(CSRcomplex_vector)::self
        complex(8),intent(in)::v
        integer,intent(in)::j
        integer::searchk
        integer::rowifirstk,rowilastk
        integer::nonzeros
        integer::m

        call boundcheck(self%N,j)
        nonzeros = self%numofdata

        rowifirstk = 1
        rowilastk = nonzeros
        if  (rowilastk == 0) then
            searchk = 1
        else
            call searchsortedfirst(self%indices(1:nonzeros),j,rowifirstk,rowilastk,searchk)
            

            if (searchk .le. rowilastk .and. self%indices(searchk) .eq. j) then
                self%val(searchk) = v
                return
            end if
        end if
        !write(*,*) "s",searchk,j,rowifirstk,rowilastk 

        if (abs(v) .ge. self%epsilon ) then

            call insert_element(self%indices,searchk,j,nonzeros)
            call insert_element(self%val,searchk,v,nonzeros)
            self%numofdata =  self%numofdata + 1

        end if

        return
    end subroutine set_c

    subroutine boundcheck(N,i) !配列外参照をチェック
        implicit none
        integer,intent(in)::N,i
        if (i < 1 .or. i > N) then
            write(*,*) "error! i should be [1:N]."
            stop
        end if
        return
    end subroutine

    subroutine searchsortedfirst(vec,i,istart,iend,searchk) 
        implicit none
        integer,intent(in)::vec(:)
        integer,intent(in)::i,istart,iend
        integer,intent(out)::searchk
        integer::k
        if (iend - istart < 0) then
            searchk = istart
            return
        end if

        searchk = iend+1!ubound(vec,1)+1
        do k=istart,iend
            if (vec(k) .eq. i) then
                searchk = k
                return
            end if
        end do

        return
    end subroutine searchsortedfirst

    subroutine insert_c(vec,pos,item,nz)
        implicit none
        complex(8),intent(inout),allocatable::vec(:)
        integer,intent(in)::pos,nz
        complex(8),intent(in)::item
        integer::vlength
        complex(8),allocatable::temp(:)

        vlength = ubound(vec,1)

        if (nz .ge. vlength) then !確保している配列valの長さが足りない時
            write(*,*) "error! the insert position is outside the vector"
            stop
        else
            
            vec(pos+1:nz+1) = vec(pos:nz)
            vec(pos) = item
        end if
        !write(*,*) vec

        return
    end subroutine insert_c

    subroutine insert_int(vec,pos,item,nz)
        implicit none
        integer,intent(inout),allocatable::vec(:)
        integer,intent(in)::pos,item,nz
        integer::vlength
        integer,allocatable::temp(:)

        vlength = ubound(vec,1)

        if (nz .ge. vlength) then !確保している配列valの長さが足りない時
            write(*,*) "error! the insert position is outside the vector"
            stop
        else
            vec(pos+1:nz+1) = vec(pos:nz)
            vec(pos) = item
        end if

        return
    end subroutine insert_int




end module csrvector