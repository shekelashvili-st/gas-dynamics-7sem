module field_mod
    implicit none 
    
    type:: Field
!��������� ���� ������������ ���������� ��������.
        character(:),allocatable:: long_name
        integer                 :: imax, jmax
        real,allocatable        :: fdata(:,:)
    contains
        procedure,private,pass(self):: Field_substract,  &
                                       Field_substractR, &
                                       Field_Rsubstract, &
                                       Field_substractRm,&
                                       Field_negate,     &
                                       Field_add,        &
                                       Field_addR,       &
                                       Field_Radd,       &
                                       Field_divideR,    &
                                       Field_multiplyR,  &
                                       Field_Rmultiply,  &
                                       Field_assignR,    &
                                       Field_Rassign,    &
                                       Field_assignRm
        generic:: operator(-)=>Field_substract, &
                               Field_substractR,Field_Rsubstract, &
                               Field_substractRm, Field_negate 
        generic:: operator(+)=>Field_add, Field_addR, Field_Radd
        generic:: operator(/)=>Field_divideR
        generic:: operator(*)=>Field_multiplyR, Field_Rmultiply
        generic:: assignment(=)=>Field_assignR, Field_Rassign, Field_assignRm
    end type Field

!������ �� ������������ ��� ������������ ����� ������
    interface Field
        module procedure:: init_Field
    end interface Field
    
!***************************************************************!    
    contains
    
    subroutine output_Fields(fields,filename)
!����� ������ � ���� � ������� Tecplot. ��������������, ��� ���� ����� �����������.
    type(Field), intent(in) :: fields(1:)
    character(*),intent(in) :: filename
    character(len=100)      :: formt
    integer                 :: iu, i, sz
    
    sz = size(fields)
    open(newunit=iu, file=filename)
    !����� ����� �����
    write(iu,'(a\)') 'VARIABLES = '
    do i=1,sz-1
        write(iu,'(a\)') '"' // fields(i).long_name // '", ' 
    end do
    write(iu,'(a)') '"' // fields(sz).long_name // '"' 
    write(iu,*) 'ZONE I=',fields(1).imax,', J=',fields(1).jmax, ', DATAPACKING=BLOCK'
    
    !����� �����, ������ DATAPACKING=BLOCK
    write(formt,*) fields(1).imax * fields(1).jmax
    formt = "(" // '100' // "es25.16)"
    do i=1,sz
        write(iu,fmt=formt) fields(i).fdata(1:fields(1).imax,1:fields(1).jmax)
    end do
    
    close(iu)
    end subroutine output_Fields
    
    
    type(Field) function init_Field(name,imax,jmax) result(res)
!����������� ��� ����
    integer,intent(in)       :: imax, jmax
    character(*),intent(in)  :: name
    res.imax = imax
    res.jmax = jmax
    res.long_name = name
    !���� ����������� �����
    allocate(res.fdata(0:imax+1,0:jmax+1))
    res.fdata = 0.0
    end function init_Field
    
    
!����� ������� ���� ��� ������� �������������� ��������
    pure type(Field) function Field_substract(self,x) result(res)
    class(Field),intent(in):: self, x
    res.fdata = self.fdata - x.fdata
    end function Field_substract
    
    
    pure type(Field) function Field_substractR(self,x) result(res)
    class(Field),intent(in):: self
    real,intent(in)        :: x
    res.fdata = self.fdata - x
    end function Field_substractR
    
    
    pure type(Field) function Field_substractRm(self,x) result(res)
    class(Field),intent(in):: self
    real,intent(in)        :: x(:,:)
    res.fdata = self.fdata - x
    end function Field_substractRm
    
    
    pure type(Field) function Field_Rsubstract(x,self) result(res)
    class(Field),intent(in):: self
    real,intent(in)        :: x(:,:)
    res.fdata = x - self.fdata 
    end function Field_Rsubstract

    
    pure type(Field) function Field_negate(self) result(res)
    class(Field),intent(in):: self
    res.fdata = -self.fdata
    end function Field_negate
    
    
    pure type(Field) function Field_add(self,x) result(res)
    class(Field),intent(in):: self, x
    res.fdata = self.fdata + x.fdata
    end function Field_add
    
    
    pure type(Field) function Field_addR(self,x) result(res)
    class(Field),intent(in):: self
    real,intent(in)        :: x
    res.fdata = self.fdata + x
    end function Field_addR
    

    pure type(Field) function Field_Radd(x,self) result(res)
    class(Field),intent(in):: self
    real,intent(in)        :: x
    res.fdata = self.fdata + x
    end function Field_Radd
    
    
    pure type(Field) function Field_divideR(self,x) result(res)
    class(Field),intent(in):: self
    real,intent(in)        :: x
    res.fdata = self.fdata / x
    end function Field_divideR

    
    pure type(Field) function Field_multiplyR(self,x) result(res)
    class(Field),intent(in):: self
    real,intent(in)        :: x
    res.fdata = self.fdata * x
    end function Field_multiplyR
    
    
    pure type(Field) function Field_Rmultiply(x,self) result(res)
    class(Field),intent(in):: self
    real,intent(in)        :: x
    res.fdata = self.fdata * x
    end function Field_Rmultiply
    
    
    subroutine Field_assignR(self,x)
    class(Field),intent(inout):: self
    real,intent(in)           :: x
    self.fdata = x
    end subroutine Field_assignR
    
    
    subroutine Field_assignRm(self,x)
    class(Field),intent(inout):: self
    real,intent(in)           :: x(:,:)
    self.fdata = x
    end subroutine Field_assignRm
    
    
    subroutine Field_Rassign(x,self)
    class(Field),intent(in):: self
    real,intent(inout)     :: x(:,:)
    x = self.fdata
    end subroutine Field_Rassign

    
    
    end module field_mod
    
    
    
    
    