program main
    use field_mod
    use processing_mod
    use output2d_mod
    
    implicit none
    
    
    real,parameter              :: pi = 4*atan(1.0_16)
    character(len=100),parameter:: input_file='params.nml'
    character(len=10)           :: prompt
    integer                     :: imax, jmax, S_max, iu
    real                        :: L, H, U_0, mu, rho_0, p_0, gamma, nu, dx, dy
    real                        :: U_ref, CFL, eps, eps_u, eps_v, eps_p
    real                        :: dtau, C_0
    type(Field)                 :: x, y, p, U, V, rho
    real,allocatable            :: Ur_barx(:,:), Ur_bary(:,:),                     &
                                   Vr_bary(:,:), Vr_barx(:,:),                     &
                                   U_flxx(:,:), V_flxy(:,:),  p_flx(:,:),        &
                                   U_flxy(:,:), V_flxx(:,:),                     &
                                   U_n(:,:), V_n(:,:), p_n(:,:), rho_n(:,:)
    integer                     :: i, j, k, S
    real,allocatable            :: tau_w(:)
    
    !Чтение параметров задачи из файла
    namelist /params/ imax,jmax,L,H,U_0,mu,rho_0,p_0,gamma, &
                      eps, S_max, U_ref, CFL
    open(newunit=iu, file=input_file, action='read')
    read(iu, nml=params)
    close(iu)
    nu = mu/rho_0
    dx = L/imax
    dy = H/jmax
    
    !Инициализация массивов, вычисление параметров
    p = Field('Pressure',imax,jmax)
    rho = Field('Density',imax,jmax)
    U = Field('U',imax,jmax)
    V = Field('V',imax,jmax)
    x = Field('X',imax,jmax)
    y = Field('Y',imax,jmax)
    do concurrent (k=1:jmax) local(i) shared(x,imax,dx)
        x.fdata(1:imax,k) = 0.5*dx + dx * [(i-1, i=1,imax)]
    end do
    do concurrent (k=1:imax) local(j) shared(y,jmax,dy)
        y.fdata(k,1:jmax) = 0.5*dy + dy * [(j-1, j=1,jmax)]
    end do
    dtau = CFL * min(dx,dy) / U_0
    C_0 = p_0/rho_0**gamma
    
    print*, 'Time step=', dtau
    print*, 'Re=', U_0*H/nu
    
    !Начальные+граничные условия
    U = U_0
    P = p_0
    V = 1e-5
    rho = rho_0
    call b_cond(U.fdata,V.fdata,p.fdata,rho.fdata,U_0,p_0,rho_0,gamma)


    open(101,file='log.txt')
    write(101,*) 'iter eps_u eps_v eps_p'
    !Итерационный процесс
    allocate(Ur_barx(imax+1,jmax),Vr_bary(imax,jmax+1),                               &
             Ur_bary(imax,jmax+1),Vr_barx(imax+1,jmax),                               & 
             U_flxx(imax+1,jmax),V_flxy(imax,jmax+1),p_flx(imax+1,jmax+1),            &
             U_flxy(imax,jmax+1),V_flxx(imax+1,jmax),                                 &
             U_n(0:imax+1,0:jmax+1),V_n(0:imax+1,0:jmax+1),                           &
             p_n(0:imax+1,0:jmax+1), rho_n(0:imax+1,0:jmax+1))
    do s=1, s_max
        !Физические потоки через границы, индексация от единицы
        Ur_barx = 0.5 * (rho.fdata(1:imax+1,1:jmax) * U.fdata(1:imax+1,1:jmax) &
                        + rho.fdata(0:imax,1:jmax) * U.fdata(0:imax,1:jmax))     
        Vr_bary = 0.5 * (rho.fdata(1:imax,1:jmax+1) * V.fdata(1:imax,1:jmax+1) &
                        + rho.fdata(1:imax,0:jmax) * V.fdata(1:imax,0:jmax))
        
        Ur_bary = 0.5 * (rho.fdata(1:imax,1:jmax+1) * U.fdata(1:imax,1:jmax+1) &
                        + rho.fdata(1:imax,0:jmax) * U.fdata(1:imax,0:jmax))     
        Vr_barx = 0.5 * (rho.fdata(1:imax+1,1:jmax) * V.fdata(1:imax+1,1:jmax) &
                        + rho.fdata(0:imax,1:jmax) * V.fdata(0:imax,1:jmax))  
        
    !Уравнение неразрывности
        !Потоки согласно схеме
        where (Ur_barx >= 0.0) U_flxx = U.fdata(0:imax,1:jmax)*rho.fdata(0:imax,1:jmax)
        where (Ur_barx < 0.0)  U_flxx = U.fdata(1:imax+1,1:jmax)*rho.fdata(1:imax+1,1:jmax)
        where (Vr_bary >= 0.0) U_flxy = U.fdata(1:imax,0:jmax)
        where (Vr_bary < 0.0)  U_flxy = U.fdata(1:imax,1:jmax+1)
        
        where (Vr_bary >= 0.0) V_flxy = V.fdata(1:imax,0:jmax)*rho.fdata(1:imax,0:jmax)
        where (Vr_bary < 0.0)  V_flxy = V.fdata(1:imax,1:jmax+1)*rho.fdata(1:imax,0:jmax)
        where (Ur_barx >= 0.0) V_flxx = V.fdata(0:imax,1:jmax)
        where (Ur_barx < 0.0)  V_flxx = V.fdata(1:imax+1,1:jmax)
        V_flxy(:,jmax+1) = 0.0
        U_flxy(:,jmax+1) = 0.0

        
        p_n(1:imax,1:jmax) = p.fdata(1:imax,1:jmax) - dtau*U_ref**2 &
                             * (diff(U_flxx,1)/dx + diff(V_flxy,2)/dy)
        rho_n(1:imax,1:jmax) = rho_0*(p_n(1:imax,1:jmax)/p_0)**(1/gamma)
        

        where (Ur_barx >= 0.0) U_flxx = U.fdata(0:imax,1:jmax)
        where (Ur_barx < 0.0)  U_flxx = U.fdata(1:imax+1,1:jmax)
        where (Vr_bary >= 0.0) V_flxy = V.fdata(1:imax,0:jmax)
        where (Vr_bary < 0.0)  V_flxy = V.fdata(1:imax,1:jmax+1)
    !Уравнение движения, проекция на ось x
        !Потоки согласно схеме, индексация от единицы
        where (Ur_barx >= 0.0) p_flx(:,1:jmax) = p.fdata(1:imax+1,1:jmax)
        where (Ur_barx < 0.0)  p_flx(:,1:jmax) = p.fdata(0:imax,1:jmax)
        
        U_n(1:imax,1:jmax) = - diff(Ur_barx*U_flxx,1)/dx - diff(Vr_bary*U_flxy,2)/dy  & !Конвекция
            - diff(p_flx(:,1:jmax),1)/dx                                              & !Градиент давления
            + 4/3*mu/dx**2 * diff2(U.fdata,1)                                         & !Вязкость, x
            + mu/dy**2 * diff2(U.fdata,2)                                             & !Вязкость, y
            + 1/3*mu/(4*dy*dx) * diff11(V.fdata)                                        !Вязкость, xy
        U_n(1:imax,1:jmax) = U_n(1:imax,1:jmax)*dtau                                  &
                             + rho.fdata(1:imax,1:jmax)*U.fdata(1:imax,1:jmax)
        U_n(1:imax,1:jmax) = U_n(1:imax,1:jmax)/rho.fdata(1:imax,1:jmax)

        
    !Уравнение движения, проекция на ось y
        !Потоки согласно схеме, индексация от единицы
        where (Vr_bary >= 0.0) p_flx(1:imax,:) = p.fdata(1:imax,1:jmax+1)
        where (Vr_bary < 0.0)  p_flx(1:imax,:) = p.fdata(1:imax,0:jmax)
        
        
        V_n(1:imax,1:jmax) = - diff(Vr_bary*V_flxy,2)/dy - diff(Ur_barx*V_flxx,1)/dx  & !Конвекция
            - diff(p_flx(1:imax,:),2)/dy                                              & !Градиент давления
            + mu/dx**2 * diff2(V.fdata,1)                                             & !Вязкость, x
            + 4/3*mu/dy**2 * diff2(V.fdata,2)                                         & !Вязкость, y
            + 1/3*mu/(4*dy*dx) * diff11(U.fdata)   
        V_n(1:imax,1:jmax) = V_n(1:imax,1:jmax)*dtau                                  &
                             + rho.fdata(1:imax,1:jmax)*V.fdata(1:imax,1:jmax)
        V_n(1:imax,1:jmax) = V_n(1:imax,1:jmax)/rho.fdata(1:imax,1:jmax)
    
    !Граничные условия
        call b_cond(U_n,V_n,p_n,rho_n,U_0,p_0,rho_0,gamma) 
    !Расчёт невязок
        eps_u = maxval(abs(U_n - U.fdata)) / dtau
        eps_v = maxval(abs(V_n - V.fdata)) / dtau
        eps_p = maxval(abs(p_n - p.fdata)) / dtau
        write(*,*) 'Iteration:', s
        write(*,*) 'eps_u=', eps_u
        write(*,*) 'eps_v=', eps_v
        write(*,*) 'eps_p=', eps_p
        print*, ' '
        write(101,'(i0,1x,es23.16,1x,es23.16,1x,es23.16)') s, eps_u, eps_v, eps_p
        if (eps_u < eps .and. eps_v < eps .and. eps_p < eps) exit
        if (mod(s,100000)==0) then
            print*, 'Exit? y/n'
            read(*,*) prompt
            if (prompt=='y') exit
        end if
        
    !Переприсвоение массивов        
        U.fdata = U_n
        V.fdata = V_n
        p.fdata = p_n
        rho.fdata = rho_n
    end do
    
    close(101)

    !Работа с трением
    allocate(tau_w(imax))
    tau_w = calc_tau_w(U,mu,dy)


    !Вывод результатов
    call output_Fields([x, y, u, v ,p, rho],'sol.dat')
    call output_tau_w(x.fdata(1:imax,1),tau_w,'tau.dat')
    
    
    end program main