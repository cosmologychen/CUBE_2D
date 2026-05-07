program test_lc_geom
   use parameters
   use omp_lib
   implicit none

   ! 1. 结构体拆解为连续数组 (SoA)，优化向量化与缓存
   real(4), allocatable :: pos_a(:,:), pos_b(:,:)
   integer(8), allocatable :: pid_a(:), pid_b(:)

   ! 输出数据缓冲区 (Batch I/O)
   integer(8), allocatable :: out_pid(:)
   real(4), allocatable :: out_dist(:), out_rproj(:), out_tproj(:), out_theta(:)

   integer(8) :: n_a, n_b, match_count
   character(250) :: file_axp, file_bxp, file_out, file_apid, file_bpid, file_info, dir_a

   type(sim_header) :: sim_info

   ! 2. 细网格链表变量 (HOC/LL)
   ! 使用 Ghost Zones (0 和 ng+1) 彻底消除搜索循环内的 mod 运算
   integer, allocatable :: hoc(:, :), ll(:)
   integer :: idx(2), ix, iy, i, j
   integer(8) :: ip, ip_b, idx_write

   ! 几何分析变量
   real(4) :: dx(2), dist, los(2), r_dist, r_proj, t_proj, cos_theta, theta
   real(4) :: obs(2)
   integer :: last_slash, io, idx_ext
   logical :: is_matched

   ! 查重与统计
   integer(1), allocatable :: seen_a(:), seen_b(:)
   integer(8) :: np_total, dup_a, dup_b, mis_a, mis_b

   ! --- A. 严格解析命令行参数 ---
   if (command_argument_count() < 3) then
      print *, "用法: ./test_lc_geom.x <A_xp.bin> <B_xp.bin> <output.bin>"
      print *, "注意: 程序将自动通过 xp.bin 推断同目录下的 pid.bin"
      stop
   end if

   call get_command_argument(1, file_axp)
   call get_command_argument(2, file_bxp)
   call get_command_argument(3, file_out)

   ! 自动定位 PID 文件 (xp.bin -> pid.bin)
   idx_ext = index(file_axp, 'xp.bin', back=.true.)
   if (idx_ext > 0) then
      file_apid = file_axp(1:idx_ext-1) // 'pid.bin'
   else
      print *, "错误: A 必须是 xp.bin 路径"
      stop
   end if

   idx_ext = index(file_bxp, 'xp.bin', back=.true.)
   if (idx_ext > 0) then
      file_bpid = file_bxp(1:idx_ext-1) // 'pid.bin'
   else
      print *, "错误: B 必须是 xp.bin 路径"
      stop
   end if

   ! 初始化基础参数
   obs = [real(ng)/2.0, real(ng)/2.0]
   np_total = ng**2

   ! 读取 A 所在目录的 info 文件
   last_slash = index(trim(file_axp), '/', back=.true.)
   if (last_slash > 0) then; dir_a = file_axp(1:last_slash)
   else; dir_a = './'; end if
   file_info = trim(dir_a) // '0.000_info.bin'

   open(10, file=file_info, access='stream', status='old', iostat=io)
   if (io == 0) then
      read(10) sim_info; close(10)
      print *, "Info Loaded: a=", sim_info%a, " np=", sim_info%np
      print *, "A path:", file_axp(1:idx_ext-1)//'*'
   else
      print *, "Warning: Cannot find info.bin, using default parameters."
   end if

   ! --- B. 读取数据并去重 (List A) ---
   n_a = get_file_size(file_apid) / 8
   print *, "A: Loaded", n_a, " particles."
   allocate(pid_a(n_a), pos_a(2, n_a), seen_a(np_total))
   seen_a = 0; dup_a = 0; match_count = 0


   open(10, file=file_axp,  access='stream', status='old')
   open(11, file=file_apid, access='stream', status='old')
   do i = 1, int(n_a)
      read(10) dx ! 借用 dx 变量暂存
      read(11) ip ! 借用 ip 变量暂存 pid
      if (ip > 0 .and. ip <= np_total) then
         if (seen_a(ip) == 1) then
            dup_a = dup_a + 1
         else
            seen_a(ip) = 1
            match_count = match_count + 1
            pid_a(match_count) = ip
            pos_a(:, match_count) = dx
         end if
      end if
   end do
   close(10); close(11)
   n_a = match_count
   print *, "List A: Loaded", n_a, " Unique particles. Dups:", dup_a

   ! --- C. 读取数据并去重 (List B) ---
   print *, "B path:", file_bxp(1:idx_ext-1)//'*'
   n_b = get_file_size(file_bpid) / 8
   print *, "B: Loaded", n_b, " particles."
   allocate(pid_b(n_b), pos_b(2, n_b), seen_b(np_total))
   seen_b = 0; dup_b = 0; match_count = 0

   open(10, file=file_bxp,  access='stream', status='old')
   open(11, file=file_bpid, access='stream', status='old')
   do i = 1, int(n_b)
      read(10) dx; read(11) ip
      if (ip > 0 .and. ip <= np_total) then
         if (seen_b(ip) == 1) then
            dup_b = dup_b + 1
         else
            seen_b(ip) = 1
            match_count = match_count + 1
            pid_b(match_count) = ip
            pos_b(:, match_count) = dx
         end if
      end if
   end do
   close(10); close(11)
   n_b = match_count
   print *, "List B: Loaded", n_b, " Unique particles. Dups:", dup_b

   ! 统计各版本独有粒子 (Missing Statistics)
   mis_a = 0; mis_b = 0
   do i = 1, int(np_total)
      if (seen_a(i) == 1 .and. seen_b(i) == 0) mis_b = mis_b + 1
      if (seen_a(i) == 0 .and. seen_b(i) == 1) mis_a = mis_a + 1
   end do
   print *, "Missing Statistics: A_only=", mis_b, " B_only=", mis_a

   ! --- D. 优化 HOC 构建 (Ghost Zones) ---
   allocate(hoc(nc,nc), ll(n_a))
   hoc = 0; ll = 0
   do i = 1, n_a
      idx = floor(pos_a(:, i)/ratio_cs) + 1 ! 转换为 1:ng 索引
      idx = modulo(idx-1,nc)+1
      ll(i) = hoc(idx(1), idx(2))
      hoc(idx(1), idx(2)) = i
   end do
   print *, "HOC with Ghost Zones built (ng=", ng, ")"

   ! --- E. OpenMP 并行匹配与几何分析 ---
   ! 预分配输出缓冲区 (Batch Buffer)
   allocate(out_pid(n_b), out_dist(n_b), out_rproj(n_b), out_tproj(n_b), out_theta(n_b))
   match_count = 0

   !$omp parallel do default(shared) &
   !$omp& private(ip_b, idx, ix, iy, j, is_matched, dx, dist, los, r_dist, r_proj, t_proj, cos_theta, theta, idx_write)
   do ip_b = 1, n_b
      idx = floor(pos_b(:, ip_b) / ratio_cs) + 1
      idx = modulo(idx-1,nc)+1

      is_matched = .false.
      
      do iy = idx(2)-1, idx(2)+1
         do ix = idx(1)-1, idx(1)+1
            j = hoc(modulo(ix-1,nc)+1, modulo(iy-1,nc)+1)
            do while (j > 0 .and. .not. is_matched)
               if (pid_a(j) == pid_b(ip_b)) then
                  is_matched = .true.

                  ! PBC 距离计算 (保持格点单位)
                  dx = pos_b(:, ip_b) - pos_a(:, j)
                  dx = modulo(dx + real(ng)/2.0, real(ng)) - real(ng)/2.0
                  dist = sqrt(sum(dx**2))

                  ! 视线方向 (LOS) 分解
                  los = obs - pos_a(:, j)
                  los = modulo(los + real(ng)/2.0, real(ng)) - real(ng)/2.0
                  r_dist = sqrt(sum(los**2))
                  if (r_dist > 0.0) los = los / r_dist

                  r_proj = sum(dx * los)
                  t_proj = dx(1)*los(2) - dx(2)*los(1)

                  if (dist > 0.0) then
                     cos_theta = max(-1.0, min(1.0, r_proj / dist))
                     theta = acos(cos_theta) * (180.0 / 3.1415926535)
                  else; theta = 0.0; end if

                  ! 安全获取写入索引 (Atomic Batch)
                  !$omp atomic capture
                  match_count = match_count + 1
                  idx_write = match_count
                  !$omp end atomic

                  out_pid(idx_write)   = pid_b(ip_b)
                  out_dist(idx_write)  = dist
                  out_rproj(idx_write) = r_proj
                  out_tproj(idx_write) = t_proj
                  out_theta(idx_write) = theta
               endif
               j = ll(j)
            end do
         end do
      end do
   end do
   !$omp end parallel do

   ! --- F. 批量整块写入 (Batch I/O) ---
   ! NOTE: 在三维模拟中，若 n_b 极其庞大，应考虑切块 (chunk) 保存以控制内存。
   open(20, file=file_out, access="stream", status="replace")
   write(20) match_count
   write(20) box
   write(20) ng
   write(20) n_a
   write(20) n_b
   write(20) mis_a
   write(20) mis_b

   if (match_count > 0) then
      write(20) out_pid(1:match_count)
      print *,'| dist: min=', minval(abs(out_dist(1:match_count))), ' max=', maxval(abs(out_dist(1:match_count))), ' mean=', sum(out_dist(1:match_count)) / real(match_count)
      write(20) out_dist(1:match_count)
      print *,'| rproj: min=', minval(abs(out_rproj(1:match_count))), ' max=', maxval(abs(out_rproj(1:match_count))), ' mean=', sum(out_rproj(1:match_count)) / real(match_count)
      write(20) out_rproj(1:match_count)
      print *,'| tproj: min=', minval(abs(out_tproj(1:match_count))), ' max=', maxval(abs(out_tproj(1:match_count))), ' mean=', sum(out_tproj(1:match_count)) / real(match_count)
      write(20) out_tproj(1:match_count)
      print *,'| theta: min=', minval(abs(out_theta(1:match_count))), ' max=', maxval(abs(out_theta(1:match_count))), ' mean=', sum(out_theta(1:match_count)) / real(match_count)
      write(20) out_theta(1:match_count)
   end if
   close(20)

   print *, "Geometric Analysis Complete!"
   print *, "Final Match Count:", match_count, " -> Saved to", trim(file_out)

contains

   function get_file_size(filename) result(res)
      character(*), intent(in) :: filename
      integer(8) :: res
      integer :: unit, io
      open(newunit=unit, file=filename, access='stream', status='old', iostat=io)
      if (io /= 0) then
         print *, "Error: Cannot open file ", trim(filename); stop
      end if
      inquire(unit, size=res); close(unit)
   end function get_file_size

end program test_lc_geom
