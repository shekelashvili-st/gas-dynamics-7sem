module output2d_mod
    implicit none
    
    contains
    
    subroutine output_Ct(x,Ct_th,Ct_ras,filename)
    real,intent(in)         :: x(:), Ct_th(:), Ct_ras(:)
    character(*),intent(in) :: filename
    character(len=100)      :: formt
    integer                 :: iu, i
    
    open(newunit=iu, file=filename)
    
    !Вывод шапки файла
    write(iu,'(a\)') 'VARIABLES = '
    write(iu,'(a)') '"X", "Ct_th", "Ct_ras"' 
    write(iu,*) 'ZONE I=',size(x),', J=',1, ', DATAPACKING=BLOCK'
    
    !Вывод полей, формат DATAPACKING=BLOCK
    write(formt,*) size(x)
    formt = "(" // trim(adjustl(formt)) // "E25.16)"
    write(iu,fmt=formt) x
    write(iu,fmt=formt) Ct_th
    write(iu,fmt=formt) Ct_ras
    end subroutine output_Ct
    
end module output2d_mod