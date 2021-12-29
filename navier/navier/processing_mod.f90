  module processing_mod
    use field_mod
    implicit none
    
    contains
    
    function diff(array,dm)
!Считаем среднюю разность по двум ячейкам, вдоль измерения dim
    real,intent(in)                 :: array(1:,1:)
    integer,intent(in)              :: dm
    real, allocatable               :: diff(:,:)
    integer                         :: imax, jmax
    
    imax = size(array,dim=1)
    jmax = size(array,dim=2)
    if (dm==1) then
        allocate(diff(imax-1,jmax))
        diff = array(2:imax,:) - array(1:imax-1,:)
    else if (dm==2) then
        allocate(diff(imax,jmax-1))
        diff = array(:,2:jmax) - array(:,1:jmax-1)
    end if
    end function diff
    
    
    function diff2(array,dm)
!Считаем среднюю разность по трем ячейкам, вдоль измерения dim  
    real,intent(in)                 :: array(0:,0:)
    integer,intent(in)              :: dm
    real, allocatable               :: diff2(:,:)
    integer                         :: imax, jmax
    
    imax = size(array,dim=1) - 2
    jmax = size(array,dim=2) - 2
    allocate(diff2(imax,jmax))
    if (dm==1) then
        diff2 = array(2:imax+1,1:jmax) - 2*array(1:imax,1:jmax) + array(0:imax-1,1:jmax)                            
    else if (dm==2) then
        diff2 = array(1:imax,2:jmax+1) - 2*array(1:imax,1:jmax) + array(1:imax,0:jmax-1)
    end if
    end function diff2
        

    subroutine b_cond(U_n,V_n,P_n,U_0,U_e,p_e)
!Постановка условий на границах
    real,intent(inout):: U_n(0:,0:), V_n(0:,0:), P_n(0:,0:)
    real,intent(in)   :: U_0, U_e, p_e
    integer           :: imax,jmax
    
    imax = size(U_n,1) - 2
    jmax = size(U_n,2) - 2
    
    U_n(0,:) = U_0
    V_n(0,:) = 0.0
    p_n(0,:) = p_n(1,:)
    
    U_n(imax+1,:) = U_n(imax,:)
    V_n(imax+1,:) = V_n(imax,:)
    p_n(imax+1,:) = 0.0
    
    U_n(:,0) = -U_n(:,1)
    V_n(:,0) = -V_n(:,1)
    p_n(:,0) = p_n(:,1)
    
    V_n(:,jmax+1) = V_n(:,jmax)
    where (V_n(:,jmax) >= 0.0 ) 
        p_n(:,jmax+1) = p_e
        U_n(:,jmax+1) = U_n(:,jmax)
    elsewhere
        p_n(:,jmax+1) = p_n(:,jmax)
        U_n(:,jmax+1) = U_e
    end where
    end subroutine b_cond
    
    pure function calc_tau_w(U,mu,dy)
!Расчёт трения на поверхности пластины
    type(Field),intent(in) :: U
    real                   :: calc_tau_w(1:U.imax)
    real,intent(in)        :: dy, mu
    integer                :: i
    calc_tau_w = [(U.fdata(i,1) - U.fdata(i,0), i=1,U.imax)] * mu/dy
    end function calc_tau_w
    
    
    
end module processing_mod   