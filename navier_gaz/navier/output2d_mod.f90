module output2d_mod
    implicit none
    
    contains
    
    subroutine output_tau_w(x,tau_w,filename)
    real,intent(in)         :: x(:), tau_w(:)
    character(*),intent(in) :: filename
    character(len=100)      :: formt
    integer                 :: iu, i
    
    open(newunit=iu, file=filename)
    
    !Вывод шапки файла
    write(iu,'(a\)') 'VARIABLES = '
    write(iu,'(a)') '"X", "tau_w"' 
    write(iu,*) 'ZONE I=',size(x),', J=',1, ', DATAPACKING=BLOCK'
    
    !Вывод полей, формат DATAPACKING=BLOCK
    write(formt,*) size(x)
    formt = "(" // '100' // "E25.16)"
    write(iu,fmt=formt) x
    write(iu,fmt=formt) tau_w
    end subroutine output_tau_w
    
end module output2d_mod