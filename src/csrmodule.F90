module csrmodule
    implicit none
    private

    public CSRcomplex
    public :: operator(*)

    type CSRcomplex
        integer::N !正方行列を仮定している。
        complex(8),allocatable::val(:)
        integer,allocatable::row(:)
        integer,allocatable::col(:)

        contains
        procedure::convert_to_dense => convert_to_dense_matrix
        procedure::set => set_c !(i,j)成分にvを代入 set(v,i,j)
        procedure::print => print_csr !行列の中身をprint
        procedure::nonzero => get_nonzeronum  !非ゼロ要素の数を数える
        procedure::update => update_c !(i,j)成分の値をアップデート。なければエラーで止まる。
        procedure::matmul => matmul_axy  !y = A*x: matmul(x,y)
        procedure::matmul2 => matmul_axy2 !y = alpha*A*x +beta*b: matmul2(x,y,b,alpha,beta)
        
    end type

    interface CSRcomplex
        module procedure::init_CSRcomplex
        module procedure::init_CSRcomplex_initialnum
    end interface CSRcomplex

    interface insert_element
        module procedure::insert_c
        module procedure::insert_int
    end interface

    interface operator(*)
        module procedure mult
    end interface


    contains 



    type(CSRcomplex) function init_CSRcomplex(N) result(A)!疎行列を初期化
        implicit none
        integer,intent(in)::N      
        A = init_CSRcomplex_initialnum(N,N) !配列の初期の長さをNとした。
        return
    end function

    type(CSRcomplex) function init_CSRcomplex_initialnum(N,initialnum) result(A)!疎行列を初期化。配列の初期の長さはinitialnumに設定。
        implicit none
        integer,intent(in)::N,initialnum

        A%N = N
        allocate(A%row(N+1))
        allocate(A%val(initialnum)) !長さinitialnumに初期化
        allocate(A%col(initialnum)) !長さinitialnumに初期化
        A%val = 0d0
        A%col = 0
        A%row = 1
    end function


    subroutine matrix_vector_c(N,val,row,col,x,y)
        implicit none
        integer,intent(in)::N
        complex(8),intent(in)::val(:)
        integer,intent(in)::row(:)
        integer,intent(in)::col(:)
        complex(8),intent(in)::x(N)
        complex(8),intent(out)::y(N)
        integer::i,j
        do i=1,N
            y(i) = 0.0d0
            do j=row(i), row(i+1)-1
                !write(*,*) j,col(j),val(j),x(col(j))
                y(i) = y(i)+val(j)*x(col(j))
            end do
            !write(*,*) "ij",i,y(i)
        end do
    end subroutine



    function mult(A,x) result(y)
        implicit none
        type(CSRcomplex),intent(in)::A
        complex(8),intent(in)::x(:)
        complex(8)::y(ubound(x,1))

        if (A%N .ne. ubound(x,1)) then
            write(*,*) "error in CSRmodules! size mismatch"
            stop
        end if
        y(1:ubound(x,1)) = 0d0

#ifdef MKL
        call mkl_zcsrgemv("N", ubound(x,1), A%val, A%row, A%col, x, y)
#else
        call matrix_vector_c(ubound(x,1), A%val, A%row, A%col, x, y)
#endif

    end function

    subroutine matmul_axy(self,x,y) 
        implicit none
        class(CSRcomplex),intent(in)::self
        complex(8),intent(in)::x(:)
        complex(8),intent(out)::y(*)
        !allocate(y(1:ubound(x,1)))
        if (self%N .ne. ubound(x,1)) then
            write(*,*) "error in CSRmodules! size mismatch"
            stop
        end if
        y(1:ubound(x,1)) = 0d0

#ifdef MKL
        call mkl_zcsrgemv("N", ubound(x,1), self%val, self%row, self%col, x, y)
#else
        call matrix_vector_c(ubound(x,1), self%val, self%row, self%col, x, y)
#endif


        return
    end subroutine

    subroutine matmul_axy2(self,x,y,b,alpha,beta) !y = alpha*A*x+ beta*b 
        implicit none
        class(CSRcomplex),intent(in)::self
        complex(8),intent(in)::x(:)
        complex(8),intent(in)::b(:)
        complex(8),intent(in)::alpha,beta
        complex(8),intent(out)::y(*)
        !allocate(y(1:ubound(x,1)))
        if (self%N .ne. ubound(x,1)) then
            write(*,*) "error in CSRmodules! size mismatch"
            stop
        end if
        y(1:ubound(x,1)) = 0d0
#ifdef MKL
        call mkl_zcsrgemv("N", ubound(x,1), self%val, self%row, self%col, x, y)
#else
        call matrix_vector_c(ubound(x,1), self%val, self%row, self%col, x, y)
#endif

        y(1:ubound(x,1)) = alpha*y(1:ubound(x,1)) + beta*b(1:ubound(x,1))

        return
    end subroutine    


    subroutine print_csr(self)
        implicit none
        class(CSRcomplex)::self
        integer::i,j,k


        do i=1,self%N
            do k=self%row(i),self%row(i+1) -1
                j = self%col(k)
                write(*,*) i,j,self%val(k)
            end do
        end do

        return
    end subroutine print_csr

    subroutine convert_to_dense_matrix(self,xdense)
        implicit none
        class(CSRcomplex)::self
        complex(8),intent(out)::xdense(:,:)
        integer::i,j,N

        N = self%N
        xdense = 0d0

        do i=1,N
            do j=self%row(i), self%row(i+1)-1
                xdense(i,self%col(j)) = self%val(j)
            end do
        end do
    end subroutine


    subroutine set_c(self,v,i,j)
        implicit none
        class(CSRcomplex)::self
        complex(8),intent(in)::v
        integer,intent(in)::i,j
        integer::rowifirstk,rowilastk
        integer::searchk
        integer::nonzeros
        integer::m

        call boundcheck(self%N,i,j)
        nonzeros = self%nonzero()
        rowifirstk = self%row(i)
        rowilastk = self%row(i+1)-1 
        call searchsortedfirst(self%col(1:nonzeros),j,rowifirstk,rowilastk,searchk)
        !write(*,*) "s",searchk,i,j,rowifirstk,rowilastk 

        if (searchk .le. rowilastk .and. self%col(searchk) .eq. j) then
            self%val(searchk) = v
            return
        end if

        if (abs(v) .ne. 0d0) then

            call insert_element(self%col,searchk,j,nonzeros)
            call insert_element(self%val,searchk,v,nonzeros)
            do m=i+1,self%N+1
                self%row(m) = self%row(m) + 1
            end do


        end if

        return
    end subroutine set_c

    subroutine insert_c(vec,pos,item,nz)
        implicit none
        complex(8),intent(inout),allocatable::vec(:)
        integer,intent(in)::pos,nz
        complex(8),intent(in)::item
        integer::vlength
        complex(8),allocatable::temp(:)

        vlength = ubound(vec,1)

        if (nz .ge. vlength) then !確保している配列valの長さが足りない時
            allocate(temp(vlength+1))
            if (pos > 1) then
                temp(1:pos-1)= vec(1:pos-1)
            end if
            temp(pos+1:nz+1) = vec(pos:nz)
            temp(pos) = item
            deallocate(vec)
            allocate(vec(vlength+vlength))
            vec(1:vlength+1) = temp(1:vlength+1)
            vec(vlength+2:vlength+vlength) = 0d0
        else
            vec(pos+1:nz+1) = vec(pos:nz)
            vec(pos) = item
        end if

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
            allocate(temp(vlength+1))
            if (pos > 1) then
                temp(1:pos-1)= vec(1:pos-1)
            end if
            temp(pos+1:nz+1) = vec(pos:nz)
            temp(pos) = item
            deallocate(vec)
            allocate(vec(vlength+vlength))
            vec(1:vlength+1) = temp(1:vlength+1)
            vec(vlength+2:vlength+vlength) = 0
        else
            vec(pos+1:nz+1) = vec(pos:nz)
            vec(pos) = item
        end if

        return
    end subroutine insert_int

    subroutine update_c(self,v,i,j)
        implicit none
        class(CSRcomplex)::self
        complex(8),intent(in)::v
        integer,intent(in)::i,j
        integer::rowifirstk,rowilastk
        integer::searchk

        call boundcheck(self%N,i,j)
        rowifirstk = self%row(i)
        rowilastk = self%row(i+1)-1 
        call searchsortedfirst(self%col,j,rowifirstk,rowilastk,searchk)


        if (searchk .le. rowilastk .and. self%col(searchk) .eq. j) then
            self%val(searchk) = v
            return
        else
            write(*,*) "error! in CSRmodules. There is no entry at",i,j
            stop
        end if
        return
    end subroutine update_c

    integer function get_nonzeronum(self) result(nonzeros)
        implicit none
        class(CSRcomplex)::self
        nonzeros = self%row(self%N+1)
        return
    end function

    subroutine boundcheck(N,i,j) !配列外参照をチェック
        implicit none
        integer,intent(in)::N,i,j
        if (i < 1 .or. i > N) then
            write(*,*) "error! i should be [1:N]."
            stop
        end if
        if (j < 1 .or. j > N) then
            write(*,*) "error! j should be [1:N]."
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

end module 