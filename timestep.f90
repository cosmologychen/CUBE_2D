
subroutine timestep
    use variables
    implicit none
    save
    integer ntry,j
    real ra,da_1,da_2,a_next,z_next,ai

    real(16) ::  k1, k4

    dt_old=dt
    sim%timestep=sim%timestep+1
    call tic(1)
    print*, ''
    print*, '-------------------------------------------------------'
    print*, 'timestep    :',sim%timestep
    dt_e=dt_max
    ntry=0
    do
        ntry=ntry+1
        da = expansion(sim%a,dt_e)
        ra=da/(sim%a+da)
        ! print*,ntry,dt_e,ra,da,sim%a+da
        if (ra>ra_max) then
        dt_e=dt_e*(ra_max/ra)
        else
        exit
        endif
        if (ntry>10) exit
    enddo
    dt = min(dt_e,sim%dt_pm1,sim%dt_pm2,dt_refine*sim%dt_pp,sim%dt_vmax)
    da = expansion(sim%a,dt)
      
    checkpoint_step=.false.
    z_next=z_checkpoint(sim%cur_checkpoint)
    a_next=1.0/(1+z_next)
    if (da>=a_next-sim%a) then
        if (z_next==z_checkpoint(sim%cur_checkpoint)) then
        checkpoint_step=.true.
        if (sim%cur_checkpoint==n_checkpoint) final_step=.true.
        endif
        do while (abs((sim%a+da)/a_next-1)>=1e-6 .or. (sim%a+da) > 1)
        dt=dt*(a_next-sim%a)/da
        da = expansion(sim%a,dt)
        enddo
    endif

    ra=da/(sim%a+da)
    a_mid=sim%a+(da/2)

    tcat(41,istep)=sim%a
    tcat(42,istep)=a_mid
    tcat(43,istep)=sim%a+da
    dtau = 0

    print*, 'tau         :',sim%tau,sim%tau+dtau
    print*, 'z         :',1.0/sim%a-1.0,1.0/(sim%a+da)-1.0
    print*, 'a         :',sim%a,a_mid,sim%a+da
    print*, 'expansion :',ra
    print*, 'dt        :',dt
    print*, 'dt_e      :',dt_e
    print*, 'dt_pm1    :',sim%dt_pm1
    print*, 'dt_pm2    :',sim%dt_pm2
    ! print*, 'dt_pm3    :',sim%dt_pm3
    print*, 'dt_pp     :',sim%dt_pp
    print*, 'dt_vmax   :',sim%dt_vmax
    print*, 'cur_powerpoint :',sim%cur_powerpoint,z_powerpoint(sim%cur_powerpoint)
    print*, ''
    sim%tau=sim%tau+dtau
    sim%t=sim%t+dt
    sim%a=sim%a+da
    call toc(1)

    contains


    real function expansion(a0,dt0)
        use variables
        ! use parameters
        implicit none
        real(8) :: a_x,adot,t_x,tdoa,a8_0
        real(4) :: a0,dt0
        integer i1,i2

        a8_0=a0
        i1 = 1
        do while(s2a(i1+1)>a8_0 .and. i1 < istep_max)
            i1 = i1+1
        enddo

        tdoa = (stime(i1+1)-stime(i1))/(s2a(i1+1)-s2a(i1))
        t_x = stime(i1)+tdoa*(a8_0-s2a(i1))+dt0


        i2 = i1
        do while(stime(i2+1)<t_x .and. i2 > 1)
            i2 = i2-1
        enddo

        adot = (s2a(i2+1)-s2a(i2))/(stime(i2+1)-stime(i2))
        a_x = s2a(i2)+adot*(t_x-stime(i2))

        expansion=a_x-a8_0
    endfunction

endsubroutine timestep


