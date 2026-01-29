subroutine drift
    use variables
    implicit none

    integer l

    print*,'drift'
    call tic(2)
    call system_clock(t1,t_rate)
    dt_mid=(dt_old+dt)/2
    print*,'  dt_mid =',dt_mid

    !$omp paralleldo default(shared) schedule(dynamic)&
    !$omp& private(l)
    do l=1,sim%np
        xp(1,l) = wrap_position(xp(1,l) + dt_mid * vp(1,l))
        xp(2,l) = wrap_position(xp(2,l) + dt_mid * vp(2,l))
    enddo
    !$omp endparalleldo
    

    call system_clock(t2,t_rate)
    print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'
    print*, ''
    call toc(2)
endsubroutine
