  module processing_mod
    use field_mod
    implicit none
    
    contains
    
    subroutine tridiag(A,B,C,D)
!Решение системы с трехдиагональной матрицей, с диагональным
!преобладанием. Результат перезаписывается в D.
    real,intent(inout):: A(1:), B(1:), C(1:), D(1:)
    integer:: info,n
    external DDTTRFB
    external DDTTRSB
    
    n = size(B)
    !Computes the factorization of a diagonally dominant
    !tridiagonal matrix.
    call DDTTRFB(n,A,B,C,info)
    if (info /= 0) print*, 'DDTTRFB error, code:', info
    
!A is overwritten by the (n-1) multipliers that define the matrix L from
!the LU factorization
    
!B is overwritten by the n diagonal element reciprocals of the upper
!triangular matrix U from the factorization of A.
    
!C is overwritten by the (n-1) elements of the superdiagonal of U.
    
    !Solves a system of linear equations with a diagonally
    !dominant tridiagonal matrix, using LU factorization
    call DDTTRSB('N',n,1,A,B,C,D,n,info)
    if (info /= 0) print*, 'DDTTRSB error, code:', info
    
    end subroutine tridiag
    
    pure function calc_tau_w(U,mu,dy)
!Расчёт трения на поверхности пластины
    type(Field),intent(in) :: U
    real                   :: calc_tau_w(1:U.imax)
    real,intent(in)        :: dy, mu
    integer                :: i,jm
    jm = U.jmax
    calc_tau_w = [(-U.fdata(i,jm-2) + 4.0*U.fdata(i,jm-1) + &
            -3*U.fdata(i,jm), i=1,U.imax)] * mu/(2.0*dy)
    end function calc_tau_w
    
end module processing_mod   