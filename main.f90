!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   CUBE™ in Coarray Fortran  !
!   haoran@xmu.edu.cn         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program CUBE
   use omp_lib
   use variables
   use runtime_lightcone_module
   implicit none
   save
   include 'fftw3.f'

   sim%cur_checkpoint=1
   call initialize
   call particle_initialization
   sim%cur_checkpoint=sim%cur_checkpoint+1
   ! sim%cur_halofind=sim%cur_halofind+1
   open(77,file=output_name('vinfo'),access='stream',status='replace')

   ! runtime_lightcone: 初始化
   if (enable_runtime_lightcone) then
      call init_runtime_lightcone()
      print*, 'Runtime Lightcone enabled, mode:', lightcone_mode
   endif

   print*, '---------- starting main loop ----------'
   call system_clock(ta1,tac)
   do istep=sim%timestep,istep_max
      ! do istep=sim%timestep,sim%timestep+50
      call system_clock(t_start,t_rate)
      call tic(100)
      write(11) xp
      close(11)
      call timestep


      if (enable_runtime_lightcone) then
         xp_new(:, :) = xp(:, :)
         vp_new(:, :) = vp(:, :)
      endif

      call drift
      call kick

      if (checkpoint_step) then
         dt_old=0
         call drift
         if (checkpoint_step) then
            if ( .not. one_run_lightcone ) call checkpoint
            sim%cur_checkpoint=sim%cur_checkpoint+1
         endif
         call print_header(sim)
         if (final_step) then
            print*,'final_step',sim%cur_checkpoint,n_checkpoint
            exit
         endif
         dt=0
      endif
      ! runtime_lightcone: 检测光锥粒子
      if (enable_runtime_lightcone) then
         call check_lightcone_crossing()
      endif
      call system_clock(t_end,t_rate)
      call toc(100)
      print*, 'total elapsed time =',tcat(100,istep),real(t_end-t_start)/t_rate,'secs';
      print*, 'Write xpos into ',output_name_step('xpos')
      open(11,file=output_name_step('xpos'), status='replace',access='stream')
      ! stop
   enddo

   ! runtime_lightcone: 关闭文件并输出统计
   if (enable_runtime_lightcone) then
      call finalize_runtime_lightcone()
   endif
   if ( one_run_lightcone ) then 
      sim%cur_checkpoint = n_checkpoint
      call checkpoint
      open(10,file=trim(opath) // '/one_run_pid.bin',status='replace',access='stream'); write(10) pid; close(10)
      open(10,file=trim(opath) // '/one_run_xp.bin',status='replace',access='stream'); write(10) xp; close(10)
      open(10,file=trim(opath) // '/one_run_vp.bin',status='replace',access='stream'); write(10) vp; close(10)
   endif

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
endprogram
