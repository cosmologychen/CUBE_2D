subroutine checkpoint
    use omp_lib
    use variables
    implicit none
    save

    character(200) fn10,fn11,fn12,fn13

    print*, 'checkpoint'

    fn10=output_name('info')
    fn11=output_name('xp')
    fn12=output_name('vp')
    fn13=output_name('a_grid')

    sim%vsim2phys=(1.5/sim%a)*box*100.*sqrt(omega_m)/ng
    sim%sigma_vi=sigma_vi
    !$omp parallelsections default(shared)
    !$omp section
    open(10,file=fn10,status='replace',access='stream'); write(10) sim; close(10)
    !$omp section
    open(11,file=fn11,status='replace',access='stream'); write(11) xp(:,:sim%np); close(11)
    !$omp section
    open(12,file=fn12,status='replace',access='stream'); write(12) vp(:,:sim%np); close(12)
    !$omp endparallelsections
    if (one_run_lightcone) then
        open(13,file=fn13,status='replace',access='stream'); write(13) a_grid; close(13)
    endif
    print*,'  wrote',sim%np,'CDM particles'
endsubroutine
