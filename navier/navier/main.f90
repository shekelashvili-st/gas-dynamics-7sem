program main
    use field_mod
    use processing_mod
    use output2d_mod
    
    implicit none
    
    
    real,parameter              :: pi = 4*atan(1.0_16)
    character(len=100),parameter:: input_file='params.nml'
    integer                     :: imax, jmax, S_max, iu
    real                        :: L, H, U_0, mu, rho, nu, dx, dy
    real                        :: U_ref, CFL, eps, eps_u, eps_v, eps_p
    real                        :: dtau, U_e, p_e
    type(Field)                 :: x, y, p, U, V
    real,allocatable            :: U_barx(:,:), U_bary(:,:),                     &
                                   V_bary(:,:), V_barx(:,:),                     &
                                   U_flxx(:,:), V_flxy(:,:),  p_flx(:,:),        &
                                   U_flxy(:,:), V_flxx(:,:),                     &
                                   U_n(:,:), V_n(:,:), p_n(:,:)
    integer                     :: i, j, k, S
    real,allocatable            :: Ct_th(:), Ct_w(:)
    
    !Чтение параметров задачи из файла
    namelist /params/ imax,jmax,L,H,U_0,mu,rho, &
                      eps, S_max, U_ref, CFL
    open(newunit=iu, file=input_file, action='read')
    read(iu, nml=params)
    close(iu)
    nu = mu/rho
    dx = L/imax
    dy = H/jmax
    
    !Инициализация массивов, вычисление параметров
    p = Field('Pressure',imax,jmax)
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
    p_e = 0.0
    U_e = U_0
    
    print*, 'Time step=', dtau
    print*, 'Re=', U_0*L/nu
    pause
    
    !Начальные+граничные условия
    U = U_0
    P = 0.0
    V = 1e-5
    call b_cond(U.fdata,V.fdata,p.fdata,U_0,U_e,p_e)


    open(101,file='log.txt')
    write(101,*) 'iter eps_u eps_v eps_p'    
    !Итерационный процесс
    allocate(U_barx(imax+1,jmax),V_bary(imax,jmax+1),                               &
             U_bary(imax,jmax+1),V_barx(imax+1,jmax),                               & 
             U_flxx(imax+1,jmax),V_flxy(imax,jmax+1),p_flx(imax+1,jmax+1),          &
             U_flxy(imax,jmax+1),V_flxx(imax+1,jmax),                               &
             U_n(0:imax+1,0:jmax+1),V_n(0:imax+1,0:jmax+1),p_n(0:imax+1,0:jmax+1))
    do s=1, s_max
        !Физические потоки через границы, индексация от единицы
        U_barx = 0.5 * (U.fdata(1:imax+1,1:jmax) + U.fdata(0:imax,1:jmax))     
        V_bary = 0.5 * (V.fdata(1:imax,1:jmax+1) + V.fdata(1:imax,0:jmax))
        
        U_bary = 0.5 * (U.fdata(1:imax,1:jmax+1) + U.fdata(1:imax,0:jmax))     
        V_barx = 0.5 * (V.fdata(1:imax+1,1:jmax) + V.fdata(0:imax,1:jmax))  
        
    !Уравнение неразрывности
        !Потоки согласно схеме
        where (U_barx >= 0.0) U_flxx = U.fdata(0:imax,1:jmax)
        where (U_barx < 0.0)  U_flxx = U.fdata(1:imax+1,1:jmax)
        where (V_bary >= 0.0) U_flxy = U.fdata(1:imax,0:jmax)
        where (V_bary < 0.0)  U_flxy = U.fdata(1:imax,1:jmax+1)
        
        where (V_bary >= 0.0) V_flxy = V.fdata(1:imax,0:jmax)
        where (V_bary < 0.0)  V_flxy = V.fdata(1:imax,1:jmax+1)
        where (U_barx >= 0.0) V_flxx = V.fdata(0:imax,1:jmax)
        where (U_barx < 0.0)  V_flxx = V.fdata(1:imax+1,1:jmax)
        V_flxy(:,1) = V_bary(:,1)
        U_flxy(:,1) = U_bary(:,1)

        
        p_n(1:imax,1:jmax) = p.fdata(1:imax,1:jmax) + &
                                 - dtau*U_ref**2 * (diff(U_flxx,1)/dx + diff(V_flxy,2)/dy)

        
    !Уравнение движения, проекция на ось x
        !Потоки согласно схеме, индексация от единицы
        where (U_barx >= 0.0) p_flx(:,1:jmax) = p_n(1:imax+1,1:jmax)
        where (U_barx < 0.0)  p_flx(:,1:jmax) = p_n(0:imax,1:jmax)
        
        U_n(1:imax,1:jmax) = - diff(U_barx*U_flxx,1)/dx - diff(V_bary*U_flxy,2)/dy  & !Конвекция
            - diff(p_flx(:,1:jmax),1)/dx                                            & !Градиент давления
            + nu/dx**2 * diff2(U.fdata,1)                                           & !Вязкость, x
            + nu/dy**2 * diff2(U.fdata,2)                                             !Вязкость, y
        U_n(1:imax,1:jmax) = U_n(1:imax,1:jmax)*dtau + U.fdata(1:imax,1:jmax)       
        
        
    !Уравнение движения, проекция на ось y
        !Потоки согласно схеме, индексация от единицы
        where (V_bary >= 0.0) p_flx(1:imax,:) = p_n(1:imax,1:jmax+1)
        where (V_bary < 0.0)  p_flx(1:imax,:) = p_n(1:imax,0:jmax)
        
        
        V_n(1:imax,1:jmax) = - diff(V_bary*V_flxy,2)/dy - diff(U_barx*V_flxx,1)/dx  & !Конвекция
            - diff(p_flx(1:imax,:),2)/dy                                            & !Градиент давления
            + nu/dx**2 * diff2(V.fdata,1)                                           & !Вязкость, x
            + nu/dy**2 * diff2(V.fdata,2)                                             !Вязкость, y
        V_n(1:imax,1:jmax) = V_n(1:imax,1:jmax)*dtau + V.fdata(1:imax,1:jmax)           
    
    !Граничные условия
        call b_cond(U_n,V_n,p_n,U_0,U_e,p_e)
    
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
        
    !Переприсвоение массивов        
        U.fdata = U_n
        V.fdata = V_n
        p.fdata = p_n
    end do
    

    !Работа с трением
    allocate(Ct_th(imax),Ct_w(imax))
    Ct_th = 0.664/sqrt(U_0*x.fdata(1:imax,1)/nu)
    Ct_w = calc_tau_w(U,mu,dy) * 2 / (rho * (maxval(U.fdata(1:imax,1:jmax),2))**2)


    !Вывод результатов
    call output_Fields([x, y, u, v ,p],'sol.dat')
    call output_Ct(x.fdata(1:imax,1),Ct_th,Ct_w,'tau.dat')
    
    
    end program main