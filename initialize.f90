    subroutine initialize
    use variables
    implicit none
    save
    include 'fftw3.f'

    logical,parameter :: read_Gks=.false.
    integer i,j,l
    istep=0; tictoc=0; tcat=0;

    print*, ''
    print*, 'CUBE run on',int(ncore,1),'cores'
    print*, '  call geometry'
    call omp_set_num_threads(ncore)
    call omp_set_max_active_levels(2)
    print*,'  omp_get_max_threads() =',omp_get_max_threads()
    print*,'  omp_get_num_procs()   =',omp_get_num_procs()
    print*,'  omp_set_threads()     =',ncore
    ! stop

    ! do i=1,nmax_redshift-1
    !   z_checkpoint(i) = i
    ! enddo
    ! print*,z_checkpoint(2:4)
    ! stop

    dt=0
    dt_old=0
    da=0 ! change for resuming checkpoints
    ! sim%cur_halofind=1
    z_checkpoint=-9999
    ! z_halofind=-9999
    checkpoint_step=.false.
    ! halofind_step=.false.
    final_step=.false.

    print*,'read z'

    open(16,file='z_checkpoint.txt',status='old')
    do i=1,nmax_redshift-1
        read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
    enddo
71 n_checkpoint=i-1
    close(16)
    if (n_checkpoint==0) stop 'z_checkpoint.txt empty'

    sim%tau=-3/sqrt(1./(1+z_checkpoint(sim%cur_checkpoint)))
    n_checkpoint=n_checkpoint
    print*,'mkdir -p '//opath
    call system('mkdir -p '//opath//'/runtime/')



    print*,'read s_a_tau_H'
    open(10,file=nupath//'s_a_tau_H.txt',form='formatted')
    read(10,*) stime
    read(10,*) s2a
    read(10,*) s2tau
    read(10,*) s2chi
    close(10)

    if (one_run_lightcone) then
        allocate(a_grid(ng,ng),D_grid(ng,ng))
    endif

    allocate(pid(ng**2));do l = 1, int(ng**2); pid(l) = l; enddo


    print*, ''
    print*, 'checkpoint information'
    print*, '  ',z_checkpoint(1),'< CDM initial conditions'
    do i=2,n_checkpoint
        print*, '  ',z_checkpoint(i)
    enddo

    print*,'  initialize Green''s functions'
    print*,'    nc,ngt,ng =',int(nc,2),int(ngt,2),int(ng,2)

    allocate(Gk1(nc/2+1 ,nc ))
    allocate(Gk2(ngt/2+1,ngt))

    if (read_Gks) then
        open(10,file=output_name_ng('Gk1'),access='stream')
        read(10) Gk1
        close(10)
        open(10,file=output_name_ng('Gk2'),access='stream')
        read(10) Gk2
        close(10)
    else
        call tic(31)
        call Green_2D(Gk1,nc,nc/2+1, nc, apm1c,   0., real(ratio_cs))
        ! call Green_2D_dyn(Gk1,nc,nc/2+1, nc, apm1c,   0., real(ratio_cs))
        call toc(31)
        open(10,file=output_name_ng('Gk1'),status='replace',access='stream')
        write(10) Gk1
        close(10)
        call tic(32)
        call Green_2D(Gk2,      ngt,      ngt/2+1,ngt, apm2,  apm1,1.)
        ! call Green_2D_dyn(Gk2,      ngt,      ngt/2+1,ngt, apm2,  apm1,1.)
        call toc(32)
        open(10,file=output_name_ng('Gk2'),status='replace',access='stream')
        write(10) Gk2
        close(10)
    endif
    ! Gk1 = Gk1*2
    ! Gk2 = Gk2*2
    print*,maxval(Gk1),minval(Gk1)
    print*,maxval(Gk2),minval(Gk2)
    ! stop


    print*, '  create fft plan'
    call tic(40)

    do iteam=1,ncore
        call sfftw_plan_dft_r2c_2d( plan2(iteam),ngt,ngt,rho2k(:,:,iteam),rho2k(:,:,iteam),FFTW_MEASURE)
        call sfftw_plan_dft_c2r_2d(iplan2(iteam),ngt,ngt,rho2k(:,:,iteam),rho2k(:,:,iteam),FFTW_MEASURE)
    enddo

    call sfftw_init_threads(l)
    print*, '    sfftw_init_threads status',l
    call sfftw_plan_with_nthreads(ncore)
    call sfftw_plan_dft_r2c_2d( plan,nw,nw,rho1k,rho1k,FFTW_MEASURE)
    call sfftw_plan_dft_c2r_2d(iplan,nw,nw,rho1k,rho1k,FFTW_MEASURE)
    call toc(40)
    print*, '    elapsed time =',tcat(5,0),'secs'

    print*,'  initialize PP neighbors'
    l=0
    do j=-nrange,-1
        do i=-nrange,nrange
            l=l+1
            ij(:,l)=[i,j]
        enddo
    enddo
    j=0
    do i=-nrange,-1
        l=l+1
        ij(:,l)=[i,j]
    enddo

    ! do i=1,n_neighbor
    !     print*,'  neighbor',i,':',ij(:,i)
    ! enddo
    ! stop

endsubroutine
