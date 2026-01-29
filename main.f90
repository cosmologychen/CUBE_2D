!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CUBE™ in Coarray Fortran  !
!   haoran@xmu.edu.cn         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program CUBE
    use omp_lib
    use variables
    implicit none
    save
    include 'fftw3.f'

    ! real :: x1(2)
    ! integer :: x_wrapped1(2)
    ! integer :: i1

    ! x1 = [0.000,1024.00000 ]
    ! do i1=1,2
    !   print*,i,'x(i)=',x1(i)
    !   x_wrapped1(i) = floor(wrap_position(x1(i)))+1
    !   print*,'  x_wrapped1(i)=',x_wrapped1(i)
    !   if (x_wrapped1(i) >= box) x_wrapped1(i) = x_wrapped1(i) - rng  
    !   if (x_wrapped1(i) <  0.0) x_wrapped1(i) = x_wrapped1(i) + rng
    !   print*,'  x_wrapped1d(i)=',x_wrapped1(i)
    ! enddo
    ! stop

    call initialize
    call particle_initialization
    sim%cur_checkpoint=sim%cur_checkpoint+1
    ! sim%cur_halofind=sim%cur_halofind+1
    open(77,file=output_name('vinfo'),access='stream',status='replace')

    print*, '---------- starting main loop ----------'
    call system_clock(ta1,tac)
    do istep=sim%timestep,istep_max
    ! do istep=sim%timestep,sim%timestep+50
        call system_clock(t_start,t_rate)
        call tic(100)
        ! print*, 'Write xpos into',output_name_step('xpos')
        ! open(11,file=output_name_step('xpos'),status='replace',access='stream')
        ! write(11) xp
        ! close(11)
        call timestep
        call drift
        call kick
        if (checkpoint_step) then
          dt_old=0
          call drift
          if (checkpoint_step) then
              call checkpoint
              sim%cur_checkpoint=sim%cur_checkpoint+1
          endif
          call print_header(sim)
          if (final_step) exit
          dt=0
        endif
        call system_clock(t_end,t_rate)
        call toc(100)
        print*, 'total elapsed time =',tcat(100,istep),real(t_end-t_start)/t_rate,'secs';
        ! stop
    enddo
    call system_clock(ta2,tac)
        print*, 'total time =',real(ta2-ta1)/tac/60,'mins';
    close(77)
    call sfftw_destroy_plan( plan)
    call sfftw_destroy_plan(iplan)
    do iteam=1,ncore
      call sfftw_destroy_plan( plan2(iteam))
      call sfftw_destroy_plan(iplan2(iteam))
    enddo
    print*,'end in ',istep
contains

  function Dgrow(scale_factor)
    implicit none
    real, parameter :: om=omega_m
    real, parameter :: ol=omega_l
    real scale_factor
    real Dgrow
    real g,ga,hsq,oma,ola
    hsq=om/scale_factor**3+(1-om-ol)/scale_factor**2+ol
    oma=om/(scale_factor**3*hsq)
    ola=ol/hsq
    g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
    ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))
    Dgrow=scale_factor*ga/g
  endfunction Dgrow

endprogram
