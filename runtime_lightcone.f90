module runtime_lightcone_module
   use omp_lib
   use parameters
   use variables, only: sim, stime, s2a, s2chi, istep_max, dt_mid, dt, da, box, ng, ncore, &
      xp, vp, pid, xp_new, vp_new, tic, toc, tcat
   implicit none

!   ! 观测者位置
!   real :: observer_x, observer_y

   ! 光锥粒子计数
   integer(8) :: np_lc = 0

   ! 文件单元号
   integer, parameter :: lc_file_xp = 111
   integer, parameter :: lc_file_vp = 112
   integer, parameter :: lc_file_pid = 113

   type type_lc_buffer
      integer :: count
      real(4), allocatable :: xp(:, :), vp(:, :)
      integer(8), allocatable :: pid(:)
   end type type_lc_buffer

   type(type_lc_buffer), allocatable :: thread_buffers(:)
   integer, parameter :: buffer_size_per_thread = 20000

contains

   ! %% 初始化光锥 (支持续跑追加)
   subroutine init_runtime_lightcone()
      implicit none
      character(250) :: fn_xp, fn_vp, fn_pid
      character(20) :: op_status, op_pos
      integer :: it

      ! observer_x = box / 2.0
      ! observer_y = box / 2.0

      if (.not. allocated(thread_buffers)) then
         allocate(thread_buffers(ncore))
         do it = 1, ncore
            allocate(thread_buffers(it)%xp(2, buffer_size_per_thread))
            allocate(thread_buffers(it)%vp(2, buffer_size_per_thread))
            allocate(thread_buffers(it)%pid(buffer_size_per_thread))
            thread_buffers(it)%count = 0
         end do
      end if

      fn_xp = trim(opath) // '/runtime_lightcone_xp.bin'
      fn_vp = trim(opath) // '/runtime_lightcone_vp.bin'
      fn_pid = trim(opath) // '/runtime_lightcone_pid.bin'

      if (sim%timestep > 0) then
         op_status = 'old'
         op_pos = 'append'
      else
         op_status = 'replace'
         op_pos = 'rewind'
      end if

      open(lc_file_xp, file=fn_xp, status=trim(op_status), position=trim(op_pos), access='stream')
      open(lc_file_vp, file=fn_vp, status=trim(op_status), position=trim(op_pos), access='stream')
      open(lc_file_pid, file=fn_pid, status=trim(op_status), position=trim(op_pos), access='stream')

      np_lc = 0
      print*, 'Runtime Lightcone Initialized (Optimized):'
      print*, '  Restart Status:', op_status, ' Position:', op_pos
   endsubroutine

   ! %% 穿越检测
   subroutine check_lightcone_crossing()
      implicit none
      real(8) :: a_i, a_f, chi_i, chi_f, t_i, t_f, da_8, dt_8, d_f_denom
      integer(8) :: ip
      integer :: iteam, iter
      real(8) :: x_i(2), x_f(2), v_i(2), v_f(2), r_i, r_f, d_i, d_f
      real(8) :: f, a_c, x_c(2), v_c(2), r_c, chi_c, epsilon
      logical :: converged

      call tic(15)

      a_f = real(sim%a, 8); a_i = a_f - real(da, 8)
      t_f = real(sim%t, 8); t_i = t_f - real(dt, 8)
      da_8 = real(da, 8); dt_8 = real(dt, 8)

      chi_i = real(find_chi_from_a_local(real(a_i, 4)), 8)
      chi_f = real(find_chi_from_a_local(real(a_f, 4)), 8)

      if (chi_f > (sqrt(2.0d0) * real(box, 8) / 2.0d0)) return

      ! grid2phys = real(box, 8) / real(ng, 8)

      !$omp parallel do default(shared) &
      !$omp& private(ip, iteam, x_i, x_f, v_i, v_f, r_i, r_f, d_i, d_f, f, a_c, x_c, v_c, r_c, chi_c,epsilon, iter, converged, d_f_denom)
      do ip = 1, sim%np
         iteam = omp_get_thread_num() + 1

         x_i = real(xp_new(:, ip), 8) * grid2phys; v_i = real(vp_new(:, ip), 8)
         x_f = real(xp(:, ip), 8) * grid2phys;     v_f = real(vp(:, ip), 8)

         r_i = sqrt((x_i(1) - observer_x)**2 + (x_i(2) - observer_y)**2)
         r_f = sqrt((x_f(1) - observer_x)**2 + (x_f(2) - observer_y)**2)

         d_i = chi_i - r_i; d_f = chi_f - r_f

         if (d_i * d_f <= 0.0 .and. d_i /= d_f) then
            f = d_i / (d_i - d_f)

            if (lightcone_mode == 2) then
               iter = 0; converged = .false.
               do while (.not. converged .and. iter <= 3)
                  iter = iter + 1
                  a_c = a_i + f * da_8
                  chi_c = real(find_chi_from_a_local(real(a_c, 4)), 8)
                  x_c = x_i + f * (x_f - x_i); v_c = v_i + f * (v_f - v_i)
                  r_c = sqrt((x_c(1) - observer_x)**2 + (x_c(2) - observer_y)**2)
                  epsilon = abs((r_c - chi_c) / (chi_c + 1e-20))
                  if (epsilon < 1.0d-7) then
                     converged = .true.
                  else
                     d_f_denom = (chi_f - chi_i) - (r_f - r_i)
                     f = f - (chi_c - r_c) / (d_f_denom + 1e-20)
                     f = max(0.0d0, min(1.0d0, f))
                  endif
               end do
            end if

            x_c = real(xp_new(:, ip), 8) + f * (real(xp(:, ip), 8) - real(xp_new(:, ip), 8))
            v_c = real(vp_new(:, ip), 8) + f * (real(vp(:, ip), 8) - real(vp_new(:, ip), 8))

            if (ip /= pid(ip)) then
               print*, 'Warning: pid(', ip, ') =', pid(ip), ' != ', ip
               stop
            end if

            call add_to_buffer(iteam, real(x_c, 4), real(v_c, 4), pid(ip))
         end if
      end do
      !$omp end parallel do

      call flush_lightcone_buffers()
      call toc(15)
   endsubroutine

   subroutine add_to_buffer(iteam, x, v, id)
      integer, intent(in) :: iteam
      real(4), intent(in) :: x(2), v(2)
      integer(8), intent(in) :: id
      integer :: old_size, new_size
      thread_buffers(iteam)%count = thread_buffers(iteam)%count + 1
      if (thread_buffers(iteam)%count > size(thread_buffers(iteam)%pid)) then
         old_size = size(thread_buffers(iteam)%pid)
         new_size = old_size * 2
         call resize_buffer(iteam, new_size)
      end if
      thread_buffers(iteam)%xp(:, thread_buffers(iteam)%count) = x
      thread_buffers(iteam)%vp(:, thread_buffers(iteam)%count) = v
      thread_buffers(iteam)%pid(thread_buffers(iteam)%count) = id
   endsubroutine

   subroutine resize_buffer(iteam, new_size)
      integer, intent(in) :: iteam, new_size
      real(4), allocatable :: tmp_xp(:, :), tmp_vp(:, :)
      integer(8), allocatable :: tmp_pid(:)
      integer :: n
      n = thread_buffers(iteam)%count - 1
      allocate(tmp_xp(2, new_size), tmp_vp(2, new_size), tmp_pid(new_size))
      tmp_xp(:, 1:n) = thread_buffers(iteam)%xp(:, 1:n)
      tmp_vp(:, 1:n) = thread_buffers(iteam)%vp(:, 1:n)
      tmp_pid(1:n) = thread_buffers(iteam)%pid(1:n)

      deallocate(thread_buffers(iteam)%xp, thread_buffers(iteam)%vp, thread_buffers(iteam)%pid)
      call move_alloc(tmp_xp, thread_buffers(iteam)%xp)
      call move_alloc(tmp_vp, thread_buffers(iteam)%vp)
      call move_alloc(tmp_pid, thread_buffers(iteam)%pid)
   endsubroutine

   subroutine flush_lightcone_buffers()
      integer :: it, l
      integer(8) :: total_this_step
      total_this_step = 0
      do it = 1, ncore
         l = thread_buffers(it)%count
         if (l > 0) then
            write(lc_file_xp) thread_buffers(it)%xp(:, 1:l)
            write(lc_file_vp) thread_buffers(it)%vp(:, 1:l)
            write(lc_file_pid) thread_buffers(it)%pid(1:l)
            total_this_step = total_this_step + l
            thread_buffers(it)%count = 0
         end if
      end do
      if (total_this_step > 0) then
         flush(lc_file_xp); flush(lc_file_vp); flush(lc_file_pid)
         np_lc = np_lc + total_this_step
         write(*, '(a,I0,a,I0,a)') '  Lightcone: +', total_this_step, ' pts (Collected: ', np_lc, ')'
      end if
   endsubroutine

   real function find_chi_from_a_local(a_in)
      use variables, only: s2a, s2chi, istep_max
      implicit none
      real, intent(in) :: a_in
      integer :: il, ir, imid
      real :: chi1, chi2, a1, a2
      il = 1; ir = istep_max
      do while (ir - il > 1)
         imid = (il + ir) / 2
         if (s2a(imid) > a_in) then; il = imid; else; ir = imid; endif
      enddo
      a1 = s2a(il); a2 = s2a(il+1); chi1 = s2chi(il); chi2 = s2chi(il+1)
      if (a2 > a1) then; find_chi_from_a_local = chi1 + (chi2 - chi1) * (a_in - a1) / (a2 - a1)
      else; find_chi_from_a_local = chi1; endif
   endfunction

   subroutine finalize_runtime_lightcone()
      implicit none
      close(lc_file_xp); close(lc_file_vp); close(lc_file_pid)
      print*, 'Runtime Lightcone Finalized. Total this run:', np_lc
   endsubroutine

endmodule runtime_lightcone_module
