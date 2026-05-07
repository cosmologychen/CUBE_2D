! ========================================================================
! lightcone.f90 - 标准化后处理光锥生成工具 (Strict Initialize Mode + Density)
! ========================================================================

program lightcone
   use parameters
   use omp_lib
   implicit none

   ! --- 全局变量声明 (同步 initialize.f90/ic.f90) ---
   character(200) :: arg_opath
   real stime(istep_max),s2a(istep_max),s2tau(istep_max),s2chi(istep_max)
   real(4), allocatable :: xp_i(:,:), vp_i(:,:)
   real(4), allocatable :: xp_n(:,:), vp_n(:,:)
   real(4), allocatable :: xp_lc(:,:), vp_lc(:,:)
   integer(izipi), allocatable :: pid_i(:), pid_n(:), pid_lc(:)
   integer(8) :: np_total, np_lc

   integer(8), allocatable :: pid_to_idx(:)

   ! --- 密度场相关 (参考 density_power.f90) ---
   real(4), allocatable :: rho_grid(:,:,:), rho_lc(:,:)
   real(8) :: rho_mean

   type(sim_header) :: sim_i, sim_n
   integer :: lc_mode, cur_checkpoint,start_checkpoint
   real :: chi_i, chi_n

   integer :: k, iter, l, iteam, idx1(2), idx2(2)
   integer(8) :: iostat, ip, idx_next
   character(200) :: arg, fn_xp, fn_vp, fn_pid, fn_rho
   character(1) :: mode_str
   real :: pos1(2), dx1(2), dx2(2)

   if ( one_run_lightcone ) stop
   ! ========================================================================
   ! 1. 初始化与数据读取
   ! ========================================================================
   call omp_set_num_threads(ncore)
   call system_clock(t1, t_rate)
   print*, '========================================'
   print*, '   Lightcone Generator (Strict Loading)'
   print*, '========================================'

   arg_opath = opath
   lc_mode = lightcone_mode
   ! observer_x = box / 2.0
   ! observer_y = box / 2.0

   if (command_argument_count() >= 2) then
      call get_command_argument(2, arg); read(arg, *) lc_mode
   endif
   print*, '  [Config] Mode:', lc_mode, ' Threads:', ncore

   ! A. 读取红移列表
   open(16, file='./z_checkpoint.txt', status='old')
   do l = 1, nmax_redshift
      read(16, end=71, fmt='(f8.4)') z_checkpoint(l)
   enddo
71 n_checkpoint = l - 1
   close(16)

   ! B. 严格按照 initialize.f90 的方式读取 TXT
   print*, '  [Setup] Reading s_a_tau_H.txt...'
   open(10, file=trim(nupath)//'s_a_tau_H.txt', form='formatted', status='old', iostat=iostat)
   if (iostat /= 0) stop '  [Error] Cannot open s_a_tau_H.txt'

   read(10, *) stime
   read(10, *) s2a
   read(10, *) s2tau
   read(10, *) s2chi
   close(10)

   print*, '  [Config] s_a_tau_H.txt'

   sim%cur_checkpoint = 1
   open(11, file=output_name('info'), access='stream'); read(11) sim_i; close(11)
   np_total = sim_i%np
   allocate(xp_i(ndim, np_total), vp_i(ndim, np_total), pid_i(np_total))
   allocate(xp_n(ndim, np_total), vp_n(ndim, np_total), pid_n(np_total))
   allocate(xp_lc(ndim, np_total*2), vp_lc(ndim, np_total*2), pid_lc(np_total*2))
   allocate(pid_to_idx(np_total))

   ! ========================================================================
   ! 2. 环带遍历逻辑
   ! ========================================================================

   do cur_checkpoint = 1, n_checkpoint - 1
      chi_n = find_chi(1/(1.0+z_checkpoint(cur_checkpoint+1)))
      if (chi_n < sqrt(2.0) * box / 2.0) exit
   enddo
   start_checkpoint = cur_checkpoint
   np_lc = 0

   chi_n = find_chi(1/(1.0+z_checkpoint(start_checkpoint)))
   sim%cur_checkpoint = start_checkpoint
   open(11, file=output_name('info'), access='stream'); read(11) sim_n; close(11)
   open(11, file=output_name('xp'),   access='stream'); read(11) xp_n;  close(11)
   open(11, file=output_name('vp'),   access='stream'); read(11) vp_n;  close(11)
   ! open(11, file=output_name('pid'),  access='stream', iostat=iostat)
   ! if (iostat == 0) then; read(11, iostat=iostat) pid_n; close(11); endif
   ! if (iostat /= 0) then; do l = 1, int(np_total); pid_n(l) = l; enddo; endif
   do l = 1, int(sim_n%np); pid_n(l) = l; enddo

   do cur_checkpoint = start_checkpoint, n_checkpoint - 1
      chi_i = chi_n
      chi_n = find_chi(1/(1.0+z_checkpoint(cur_checkpoint+1)))
      ! 预判逻辑：跳过超出盒���物理范围的环带
      if (chi_n > sqrt(2.0) * box / 2.0) cycle

      sim_i = sim_n; xp_i = xp_n; vp_i = vp_n; pid_i = pid_n
      print*, '---------------------------------------------------'
      print*, 'Band', cur_checkpoint, ': z =', z_checkpoint(cur_checkpoint), '->', z_checkpoint(cur_checkpoint+1)
      print*, '  Chi Radius:', chi_i, '->', chi_n


      sim%cur_checkpoint = cur_checkpoint + 1
      open(11, file=output_name('info'), access='stream'); read(11) sim_n; close(11)
      open(11, file=output_name('xp'),   access='stream'); read(11) xp_n;  close(11)
      open(11, file=output_name('vp'),   access='stream'); read(11) vp_n;  close(11)
      ! open(11, file=output_name('pid'),  access='stream', iostat=iostat)
      ! if (iostat == 0) then; read(11, iostat=iostat) pid_n; close(11); endif
      ! if (iostat /= 0) then; do l = 1, int(sim%np); pid_n(l) = l; enddo; endif
      do l = 1, int(sim_n%np); pid_n(l) = l; enddo
      print*, 'Band', cur_checkpoint, ': a =', sim_i%a, '->', sim_n%a

      pid_to_idx = 0
      !$omp parallel do private(ip)
      do ip = 1, np_total
         if (pid_n(ip) > 0 .and. pid_n(ip) <= np_total) pid_to_idx(pid_n(ip)) = ip
      enddo
      !$omp end parallel do

      call process_band_final(sim_i, sim_n, xp_i, vp_i, pid_i, xp_n, vp_n, chi_i, chi_n)
   enddo

   ! ========================================================================
   ! 3. 汇总密度场 (参考 density_power.f90)
   ! ========================================================================
   if (np_lc > 1e5) then
      print*, '---------------------------------------------------'
      write(*, '(a,i8.0,a,i8.0,a)') '  Generating Lightcone Density Field (nw =', nw, ', threads =', ncore, ')...'
      allocate(rho_grid(nw, nw, ncore), rho_lc(nw, nw))
      rho_grid = 0.0

      !$omp parallel do default(shared) private(ip, iteam, pos1, idx1, idx2, dx1, dx2)
      do ip = 1, np_lc
         iteam = omp_get_thread_num() + 1
         if (xp_lc(1, ip) /= xp_lc(1, ip) .or. xp_lc(2, ip) /= xp_lc(2, ip)) cycle

         pos1 = xp_lc(:, ip) * (real(nw)/real(ng)) - 0.5
         idx1 = floor(pos1); idx2 = idx1 + 1
         dx1 = idx1 + 1.0 - pos1; dx2 = 1.0 - dx1

         call assign_rho_robust(idx1(1)+1, idx1(2)+1, iteam, dx1(1)*dx1(2))
         call assign_rho_robust(idx1(1)+1, idx2(2)+1, iteam, dx1(1)*dx2(2))
         call assign_rho_robust(idx2(1)+1, idx1(2)+1, iteam, dx2(1)*dx1(2))
         call assign_rho_robust(idx2(1)+1, idx2(2)+1, iteam, dx2(1)*dx2(2))
      enddo
      !$omp end parallel do

      rho_lc = 0.0
      do iteam = 1, ncore
         rho_lc = rho_lc + rho_grid(:, :, iteam)
      enddo

      rho_mean = sum(real(rho_lc, 8)) / real(nw*nw, 8)
      if (rho_mean > 0) then
         rho_lc = rho_lc / real(rho_mean) - 1.0
      endif

      call save_all_results()
   endif

   deallocate(xp_i, vp_i, pid_i, xp_n, vp_n, pid_n, xp_lc, vp_lc, pid_lc, pid_to_idx)
   if (allocated(rho_grid)) deallocate(rho_grid, rho_lc)
   call system_clock(t2, t_rate)
   write(*, '(a,I8.0,a,f8.4,a)') '  [Done] Total collected:', int(sqrt(real(np_lc, 8))), '^2 Particles. Time:', real(t2-t1)/t_rate, 's'

contains

   subroutine assign_rho_robust(ix, iy, it, val)
      integer, intent(in) :: ix, iy, it
      real, intent(in) :: val
      integer :: iix, iiy
      iix = mod(ix - 1, nw); if (iix < 0) iix = iix + nw; iix = iix + 1
      iiy = mod(iy - 1, nw); if (iiy < 0) iiy = iiy + nw; iiy = iiy + 1
      rho_grid(iix, iiy, it) = rho_grid(iix, iiy, it) + val
   endsubroutine

   subroutine process_band_final(si, sn, xi, vi, pi, xn, vn, ci, cn)
      type(sim_header), intent(in) :: si, sn
      real(4), intent(in) :: xi(:,:), vi(:,:), xn(:,:), vn(:,:)
      integer(izipi), intent(in) :: pi(:)
      real, intent(in) :: ci, cn
      real :: r_i, r_n, d_i, d_n, frac, a_c, x_c(2), v_c(2), r_c, chi_c, da, dt_p, c_eff
      integer(8) :: ip, idx_n, count

      da = sn%a - si%a; dt_p = abs(sn%t - si%t)
      c_eff = - (si%a + sn%a)*0.5 * (cn - ci) / (dt_p + 1e-20)
      count = 0
      print*, '  ci:', ci, 'cn:', cn, 'da:', da, 'dt_p:', dt_p, 'c_eff:', c_eff
      !$omp parallel do default(shared) private(ip, idx_n, r_i, r_n, d_i, d_n, frac, a_c, x_c, v_c, r_c, chi_c, iter) reduction(+:count) !! num_threads(1)
      do ip = 1, np_total
         idx_n = pid_to_idx(pi(ip)); if (idx_n == 0) cycle
         r_i = sqrt((xi(1, ip)*grid2phys - observer_x)**2 + (xi(2, ip)*grid2phys - observer_y)**2)
         r_n = sqrt((xn(1, idx_n)*grid2phys - observer_x)**2 + (xn(2, idx_n)*grid2phys - observer_y)**2)
         d_i = ci - r_i; d_n = cn - r_n
         ! print*, '    ip:', ip, 'idx_n:', idx_n, 'r_i:', r_i, 'r_n:', r_n, 'd_i:', d_i, 'd_n:', d_n, 'frac:', frac, 'a_c:', a_c, 'r_c:', r_c, 'chi_c:', chi_c, 'c_eff:', c_eff
         if (d_i * d_n <= 0.0 .and. d_i /= d_n) then
            if (lc_mode == 0) then
               x_c = xi(:, ip); v_c = vi(:, ip)
            else
               frac = d_i / (d_i - d_n)
               if (lc_mode == 2 .and. c_eff > 0) then
                  do iter = 1, 2
                     a_c = si%a + frac * da; chi_c = find_chi(a_c)
                     x_c = xi(:, ip) + frac * (xn(:, idx_n) - xi(:, ip))
                     r_c = sqrt((x_c(1)*box/ng - observer_x)**2 + (x_c(2)*box/ng - observer_y)**2)
                     frac = frac + (chi_c - r_c) / (c_eff / a_c * dt_p + 1e-20)
                     frac = max(0.0, min(1.0, frac))
                  end do
               endif
               x_c = xi(:, ip) + frac * (xn(:, idx_n) - xi(:, ip))
               v_c = vi(:, ip) + frac * (vn(:, idx_n) - vi(:, ip))
            endif
            !$omp critical
            np_lc = np_lc + 1
            if (np_lc <= np_total*2) then
               xp_lc(:, np_lc) = x_c; vp_lc(:, np_lc) = v_c; pid_lc(np_lc) = pi(ip)
            endif
            !$omp end critical
            count = count + 1
         endif
      enddo
      !$omp end parallel do
      print*, '    Crossing count:', count
   endsubroutine

   subroutine save_all_results()
      write(mode_str, '(I1)') lc_mode
      fn_xp = trim(arg_opath) // '/lightcone_post_m' // mode_str // '_xp.bin'
      fn_vp = trim(arg_opath) // '/lightcone_post_m' // mode_str // '_vp.bin'
      fn_pid = trim(arg_opath) // '/lightcone_post_m' // mode_str // '_pid.bin'
      fn_rho = trim(arg_opath) // '/lightcone_post_m' // mode_str // '_delta.bin'
      open(10, file=fn_xp, status='replace', access='stream'); write(10) xp_lc(:, 1:np_lc); close(10)
      open(11, file=fn_vp, status='replace', access='stream'); write(11) vp_lc(:, 1:np_lc); close(11)
      open(12, file=fn_pid, status='replace', access='stream'); write(12) pid_lc(1:np_lc); close(12)
      open(13, file=fn_rho, status='replace', access='stream'); write(13) rho_lc; close(13)
      print*, '  [Save] Collection and Density saved to:', trim(arg_opath)
   endsubroutine

   real function find_chi(a_in)
      real, intent(in) :: a_in
      integer :: il, ir, imid
      il = 1; ir = istep_max
      do while (ir - il > 1)
         imid = (il + ir) / 2
         if (s2a(imid) > a_in) then; il = imid; else; ir = imid; endif
      enddo
      find_chi = s2chi(il) + (s2chi(il+1) - s2chi(il)) * (a_in - s2a(il)) / (s2a(il+1) - s2a(il) + 1d-20)
   endfunction

endprogram
