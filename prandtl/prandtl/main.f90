program main
    use field_mod
    use processing_mod
    use output2d_mod
    
    implicit none
    
    real,parameter              :: pi = 4*atan(1.0_16)
    character(len=100),parameter:: input_file='params.nml'
    integer                     :: imax, jmax, S_max, iu
    real                        :: L, H, U_0, mu, rho, nu, dx, dy
    real                        :: eps, eps_u, eps_v
    type(Field)                 :: x, y, p, U, V
    real,allocatable            :: A(:), B(:), C(:), U_new(:), V_new(:)
    integer                     :: i, j, k, S
    real,allocatable            :: Ct_th(:), Ct_w(:)
    
    !Чтение параметров задачи из файла
    namelist /params/ imax,jmax,L,H,U_0,mu,rho, &
                      eps, S_max
    open(newunit=iu, file=input_file, action='read')
    read(iu, nml=params)
    close(iu)
    nu = mu/rho
    dx = L/(imax-1)
    dy = H/(jmax-1)
    print*, 'Re_L=', u_0*L/nu
    pause
    
    !Инициализация массивов
    p = Field('Pressure',imax,jmax)
    U = Field('U',imax,jmax)
    V = Field('V',imax,jmax)
    x = Field('X',imax,jmax)
    y = Field('Y',imax,jmax)
    do concurrent (k=1:jmax) local(i) shared(x,imax,L)
        x.fdata(:,k) = L/(imax-1) * [(i-1, i=1,imax)]
    end do
    do concurrent (k=1:imax) local(j) shared(y,jmax)
        y.fdata(k,:) = H/(jmax-1) * [(j-1, j=1,jmax)]
    end do
    
    !Начальные условия
    U.fdata(1,:) = U_0
    P = 0.0
    V = 1e-5
    
    !Итерационный процесс
    allocate(A(jmax),B(jmax),C(jmax),U_new(jmax), V_new(jmax))
    do i=2,imax
        U.fdata(i,:) = U.fdata(i-1,:)
        V.fdata(i,:) = V.fdata(i-1,:)
        eps_u = 0.0
        eps_v = 0.0
        do S=1,S_max
            !Коэффициенты для прогонки
            A = [0.0, &
                (-V.fdata(i,j-1)/(2*dy) - nu/dy**2, j=2,jmax-1), &
                0.0]
            B = [1.0, &
                (U.fdata(i,j)/dx + 2*nu/dy**2, j=2,jmax-1), &
                1.0]
            C = [0.0, &
                (V.fdata(i,j+1)/(2*dy) - nu/dy**2, j=2,jmax-1), &
                0.0]
            U_new = [0.0, &
                           (U.fdata(i-1,j)**2/dx,j=2,jmax-1), &
                           U_0]
            
            !Вызов прогонки
            call tridiag(A(2:),B,C(:jmax-1),U_new)
            V_new = 0.0
            do j=2,jmax
                V_new(j) = V_new(j-1) - dy/2 * 1/dx * &
                    (U_new(j)-U.fdata(i-1,j)+U_new(j-1)-U.fdata(i-1,j-1))
            end do
    
            !Расчёт погрешности
            eps_u = maxval(abs(U_new - U.fdata(i,:))) / maxval(abs(U_new)) 
            eps_v = maxval(abs(V_new - V.fdata(i,:))) / maxval(abs(V_new))
            write(*,*) 'Iteration:', s
            write(*,*) 'eps_u=', eps_u
            write(*,*) 'eps_v=', eps_v
            print*, ' '
            
            if (eps_u < eps .and. eps_v < eps) exit
            
            V.fdata(i,:) = V_new
            U.fdata(i,:) = U_new
        end do 
    end do
    
    !Работа с трением
    allocate(Ct_th(imax),Ct_w(imax))
    Ct_th = 0.664/sqrt(U_0*x.fdata(:,1)/nu)
    Ct_w = calc_tau_w(U,mu,dy) * 2 / (rho * U_0**2)


    !Вывод результатов
    call output_Fields([x, y, u, v ,p],'sol_new.dat')
    call output_Ct(x.fdata(2:,1),Ct_th(2:),Ct_w(2:),'tau_new.dat')
    
    
    end program main