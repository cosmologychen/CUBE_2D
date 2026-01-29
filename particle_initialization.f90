subroutine particle_initialization
    use parameters
    use variables, only: xp,vp
    implicit none

    character(200) fn10,fn11,fn12,fn13,fn14
    call tic(41)

    print*, ''
    print*, 'particle_initialization'
    print*, '  at redshift',z_checkpoint(sim%cur_checkpoint)
    fn10=output_name('info')
    fn11=output_name('xp')
    fn12=output_name('vp')

    open(10,file=fn10,status='old',access='stream')
    read(10) sim
    close(10)
    sim%a=1./(1+z_checkpoint(sim%cur_checkpoint))
    ! sigma_vi=sim%sigma_vi

    if (sim%izipx/=izipx .or. sim%izipv/=izipv) then
        print*, '  zip format incompatable'
        stop
    endif

    allocate(xp(2,sim%np),vp(2,sim%np))
    
    !$omp parallelsections default(shared)
    !$omp section
        !omp workshare
        open(11,file=fn11,status='old',access='stream'); read(11) xp(:,:sim%np); close(11)
        !omp endworkshare
    !$omp section
        open(12,file=fn12,status='old',access='stream'); read(12) vp(:,:sim%np); close(12)
    print*,'  read',sim%np,' particles'
    !$omp endparallelsections
    call toc(41)
    print*,'  np           =', sim%np
    print*,'  mass_p_cdm   =', sim%mass_p_cdm
    print*,'  vsim2phys    =',sim%vsim2phys, ' (km/s)/(1.0)'
    ! print*,'  sigma_vi     =',sigma_vi,'(simulation unit)'
    print*,'  elapsed time =',tcat(6,0),'secs'
    print*,''
    ! stop
    call print_header(sim)



    ! !xp shift
    ! xp=mod(xp+ng/2+real(ng),real(ng))
endsubroutine
