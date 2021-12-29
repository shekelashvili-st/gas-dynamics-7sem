program main
    use field_mod
    use processing_mod
    use output2d_mod
    
    implicit none
    
    real,parameter              :: pi = 4*atan(1.0_16)
    character(len=100),parameter:: input_file='params.nml'
    integer                     :: imax, jmax, S_max, iu
    real                        :: L, H, U_0, mu, dx, dy, gamma, rho_0, p_0, C_0
    real                        :: eps, eps_u, eps_v, eps_p
    real                        :: G, G_new
    type(Field)                 :: x, y, p, U, V, rho
    real,allocatable            :: A(:), B(:), C(:), U_new(:), V_new(:), rho_new(:), p_new(:)
    integer                     :: i, j, k, S, S_p
    real,allocatable            :: tau_w(:)
    
    !Чтение параметров задачи из файла
    namelist /params/ imax,jmax,L,H,U_0,mu,rho_0, &
                      eps, S_max, gamma, p_0
    open(newunit=iu, file=input_file, action='read')
    read(iu, nml=params)
    close(iu)
    dx = L/(imax-1)
    dy = H/(jmax-1)
    C_0 = p_0/(rho_0**gamma)
    
    !Инициализация массивов
    p = Field('Pressure',imax,jmax)
    rho = Field('Density',imax,jmax)
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
    V = 0.0
    p = p_0
    rho = rho_0
    print*, 'Re = ', rho_0*U_0*H/mu
    pause
    allocate(A(jmax),B(jmax),C(jmax), &
             U_new(jmax),V_new(jmax),rho_new(jmax),p_new(jmax))    

    !Итерационный процесс по координатам
    do i=2,imax
        U.fdata(i,:) = U.fdata(i-1,:)
        V.fdata(i,:) = V.fdata(i-1,:)
        p.fdata(i,:) = p.fdata(i-1,:)
        rho.fdata(i,:) = rho.fdata(i-1,:)
        eps_u = 0.0
        eps_v = 0.0
        eps_p = 0.0
        p_new(:) = p.fdata(i,:)
        rho_new(:) = rho.fdata(i,:)
        G = 0.5*dy * sum(rho.fdata(i,2:jmax)*U.fdata(i,2:jmax) + &
                       rho.fdata(i,1:jmax-1)*U.fdata(i,1:jmax-1))
        print*, 'G=', G
!Итерационный процесс по давлению
do S_p=1,S_max
        do S=1,S_max
            !Коэффициенты для прогонки
            A = [0.0, &
                (-rho.fdata(i,j-1)*V.fdata(i,j-1)/(2*dy) - mu/dy**2, j=2,jmax-1), &
                0.0]
            B = [-1.0, &
                (rho.fdata(i,j)*U.fdata(i,j)/dx + 2*mu/dy**2, j=2,jmax-1), &
                1.0]
            C = [1.0, &
                (rho.fdata(i,j+1)*V.fdata(i,j+1)/(2*dy) - mu/dy**2, j=2,jmax-1), &
                0.0]
            U_new = [0.0, &
                    (rho.fdata(i-1,j)*U.fdata(i-1,j)**2/dx - &
                    (p_new(j) - p.fdata(i-1,j))/dx,j=2,jmax-1), &
                    0.0]
            !Вызов прогонки
            call tridiag(A(2:),B,C(:jmax-1),U_new)
            V_new = 0.0
            do j=2,jmax-1
                V_new(j) = rho.fdata(i,j-1)*V_new(j-1) - dy/2 * 1/dx * &
                    (rho.fdata(i,j)*U_new(j) - rho.fdata(i-1,j)*U.fdata(i-1,j) + &
                    rho.fdata(i,j-1)*U_new(j-1) - rho.fdata(i-1,j-1)*U.fdata(i-1,j-1))
                V_new(j) = V_new(j)/rho.fdata(i,j)
            end do
    
            !Расчёт погрешности
            eps_u = maxval(abs(U_new - U.fdata(i,:))) / maxval(abs(U_new)) 
            eps_v = maxval(abs(V_new - V.fdata(i,:))) / maxval(abs(V_new))
            write(*,*) 'Iteration(s_p,s_u):', S_p, S
            write(*,*) 'eps_u=', eps_u
            write(*,*) 'eps_v=', eps_v
            print*, ' '
            
            if (eps_u < eps .and. eps_v < eps) exit
            
            V.fdata(i,:) = V_new
            U.fdata(i,:) = U_new
        end do
    
    p.fdata(i,:) = p_new(:)
    G_new = 0.5*dy * sum(rho.fdata(i,2:jmax)*U.fdata(i,2:jmax) + &
                       rho.fdata(i,1:jmax-1)*U.fdata(i,1:jmax-1))
    p_new(:) = p.fdata(i,:) - G/(H**2 * rho.fdata(i,:)) * (G - G_new)
    rho.fdata(i,:) = rho_new(:)
    rho_new(:) = rho_0*(p_new(:)/p_0)**(1/gamma)
    
    eps_p = abs(p_new(2) - p.fdata(i,2))/abs(p_new(2))
    write(*,*) 'eps_p=', eps_p
    if (eps_p < eps) exit
end do
    end do
    
    !!Проверка расхода
    !do i=1, imax
    !    print*,'i=', i
    !    print*, 'G=', 0.5*dy * sum(rho.fdata(i,2:jmax)*U.fdata(i,2:jmax) + &
    !                  rho.fdata(i,1:jmax-1)*U.fdata(i,1:jmax-1))
    !enddo
    
    !Работа с трением
    allocate(tau_w(imax))
    tau_w = calc_tau_w(U,mu,dy)


    !Вывод результатов
    call output_Fields([x, y, u, v, p, rho],'sol.dat')
    call output_tauw(x.fdata(2:,1),tau_w(2:),'tau.dat')
    
    
    end program main