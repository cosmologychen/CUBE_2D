program DT_void
  use parameters
  implicit none

  real,parameter :: rate_in = 1!sqrt(2)
  real :: np_halo_min_r,  b_link_r
  integer :: halo_np,  nhalos, nhalos_r

  integer(8) i,j,l,cur_checkpoint,ip,ip1,ip2,np_void,idl(2),cm(2),rm,iq(2)
  real,allocatable :: xp_mean(:,:),xp(:,:),triangles(:,:,:)
  integer(8),allocatable :: hoc(:,:),ll(:)
  integer :: iunit, iostat
  character(len=400) :: cmd
  integer :: unit, filesize, num_reals, ntri ,t ,n_void

  real :: center(2),radius, a(2), b(2), c(2), r2, d2,xpos(2),qp_mean(2)
  integer :: tmax(2),tmin(2),current_pos


  open(16,file='../z_checkpoint.txt',status='old')
  do i=1,nmax_redshift-1
    read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
  enddo
  71 n_checkpoint=i-1
  close(16)
  if (n_checkpoint==0) stop 'z_checkpoint.txt empty'
  
  sim%cur_checkpoint = 5



  ! ! 打开文件
  ! open(11, file=output_name('halo_xp_mean'), access='stream', form='unformatted', status='old', action='read', iostat=iostat)
  ! if (iostat /= 0) then
  !   print *, '无法打开文件: ', output_name('halo_xp_mean')
  !   stop
  ! endif

  ! ! 读取文件头
  ! read(11) b_link_r, np_halo_min_r,nhalos_r
  ! print *, 'b_link = ', b_link_r
  ! print *, 'np_halo_min = ', np_halo_min_r
  ! print *, 'nhalos = ', nhalos_r
  ! allocate(xp_mean(2, nhalos_r))

  ! nhalos = 0
  ! do
  !   ! 读取 halo_np
  !   read(11, iostat=iostat) halo_np
  !   if (iostat /= 0) exit
    
  !   ! 检查结束标记
  !   if (halo_np == 0) then
  !     print *, '找到文件结束标记'
  !     exit
  !   endif
    
  !   ! 读取 xp_mean
  !   nhalos = nhalos + 1
  !   read(11, iostat=iostat) xp_mean(:, nhalos)
  !   if (iostat /= 0) exit
    
  !   print *, 'Halo ', nhalos, ': NP = ', halo_np, &
  !             ', X_mean = ', xp_mean(1, nhalos), ', Y_mean = ', xp_mean(2, nhalos)
  ! enddo

  ! close(11)

  ! if (nhalos /= nhalos_r) stop '读取的halo数量与文件头不符'

  ! print *, '总共读取了 ', nhalos, ' 个halo数据'

  ! xpos = [10,0]
  ! a = [0,0]
  ! b = [10,0]
  ! c = [0,10]
  ! print*,point_in_triangle(xpos, a, b, c)
  ! stop

  ! xpos = [10,11]
  ! center = [800,1022]
  ! ! print*,wrap_position2(xpos - center)
  ! ! xpos = xpos-center
  ! ! print*,xpos,rng,(rng/2)
  ! ! do i=1,2
  ! !   xpos(i) =  xpos(i) - floor(xpos(i) / (rng/2)) * (rng/2)
  ! !   print*,xpos
  ! !   if (xpos(i) >= ng ) xpos(i) = xpos(i) - (rng/2) 
  ! !   if (xpos(i) <  0.0) xpos(i) = xpos(i) + (rng/2)
  ! !   print*,xpos
  ! ! enddo
  ! ! print*,xpos

  ! print*,pbc_vec(xpos - center)
  ! xpos = [790,1024]
  ! center = [10,10]
  ! print*,pbc_vec(xpos - center)
  ! stop


  ! do ip = 1,ng*4,100
  !   iq(1)=(ip-1)/ng
  !   iq(2)=modulo(ip-1,int(ng,4))
  !   print*,ip,iq
  ! enddo
  ! stop

  !print*,output_name('info')
  open(11,file=output_name('info'),access='stream'); read(11) sim; close(11)
  print*, 'np =',sim%np
  allocate(xp(2,sim%np))
  open(11,file=output_name('xp'),access='stream'); read(11) xp; close(11)

  allocate(ll(sim%np),hoc(ng,ng))
  ll = 0 ;hoc = 0
  do ip=1,sim%np
    xpos=xp(:,ip)
    idl=xpos2mesh(xpos,ng)
    ll(ip)=hoc(idl(1),idl(2))
    hoc(idl(1),idl(2))=ip
  enddo ! ip


  print*, 'Reading DTFE...'

  ! 构建调用命令
  write(cmd, '(A,A)') "~/anaconda3/envs/camb/bin/python3 DTFE.py ", trim(output_name('halo_xp_mean_only'))

  ! 调用 Python 脚本
  print*, "Calling Python script..."
  print*, "Command: ", trim(cmd)
  call system(trim(cmd))

  inquire(file="triangles.bin", size=filesize)
  num_reals = filesize / 4
  ntri = num_reals / (2 * 3)  ! 每个三角形有 2*3 个坐标

  allocate(triangles(2,3,ntri))
  open(unit=20, file="triangles.bin", form="unformatted", access="stream")
  read(20) triangles
  close(20)

  print*, "ntri: ", ntri
  print*,maxval(triangles),minval(triangles)



  print*, "Output: ", output_name('DT_void')
  open(11,file=output_name('DT_void'),status='replace',access='stream')
  n_void = 0
  do t = 1, ntri
    a = triangles(:,1,t)
    b = triangles(:,2,t)
    c = triangles(:,3,t)  
    call calculate_circumcircle(a, b, c, center, radius)
    if (center(1) .lt. 0.0 .or. center(1) .ge. real(ng) .or. &
        center(2) .lt. 0.0 .or. center(2) .ge. real(ng)) then
        ! print*, "Outside: ", t, a, b, c, center, radius
        cycle
    endif
    if (radius .gt. 5.0) then
      n_void = n_void + 1
      write(11) center, radius, a, b, c
    endif
  enddo
  close(11)
  print*,n_void
  stop



  open(12,file=output_name('DT_void_r'),status='replace',access='stream')
  open(13,file=output_name('DT_void_t'),status='replace',access='stream')
  write(12) n_void
  write(13) n_void

  n_void = 0
  do t = 1, ntri
    a = triangles(:,1,t)
    b = triangles(:,2,t)
    c = triangles(:,3,t)  
    call calculate_circumcircle(a, b, c, center, radius)
    if (center(1) .lt. 0.0 .or. center(1) .ge. real(ng) .or. &
        center(2) .lt. 0.0 .or. center(2) .ge. real(ng)) then
        ! print*, "Outside: ", t, a, b, c, center, radius
        cycle
    endif
    if (radius .gt. 20.0) then
      write(12) center, radius, a, b, c
      ! print*,center, radius, a, b, c

      ! r
      np_void = 0
      write(12) np_void
      r2 = radius**2
      qp_mean = 0
      cm=xpos2mesh(center,ng)
      rm=floor(rate_in*radius)+1
      do i=-rm,rm
      do j=-rm,rm
        idl(1) = mod(cm(1)+i+ng-1,ng)+1
        idl(2) = mod(cm(2)+j+ng-1,ng)+1
        if (idl(1) < 1 .or. idl(1) > ng .or. idl(2) < 1 .or. idl(2) > ng) then
          print*,idl
          print*,cm
          print*,i,j
          print*,cm(1)+i+ng,cm(2)+j+ng
          print*,mod(cm(1)+i+ng,ng)+1,mod(cm(2)+j+ng,ng)+1
          stop 'out of range'
        endif
        ip = hoc(idl(1),idl(2))
        do while (ip /= 0)
          xpos=pbc_vec(xp(:,ip) - center)
          d2 = sum(xpos**2)

          if (d2 < r2) then
            np_void = np_void + 1
            iq(1)=(ip-1)/ng
            iq(2)=modulo(ip-1,int(ng,4))
            xpos=iq+0.5
            qp_mean = qp_mean + pbc_vec(xpos-center)
            write(12) ip
            write(12) xpos
          endif
          ip = ll(ip)
        enddo
      enddo
      enddo
      inquire(12, pos=current_pos)  ! 获取当前指针位置
      write(12, pos=(current_pos - np_void * 8 * 2 - 8 )) np_void
      qp_mean = qp_mean / np_void + center
      write(12, pos=(current_pos)) qp_mean
      ! print*,center,qp_mean, radius, np_void
      np_void = 0
      write(12) np_void




      !t
      
      write(13) center, radius, a, b, c
      write(13) np_void
      np_void=0
      qp_mean = 0
      tmax(1) = min(1024,floor(max(max(a(1),b(1)),c(1)))+1)
      tmax(2) = min(1024,floor(max(max(a(2),b(2)),c(2)))+1)
      tmin(1) = max(1,floor(min(min(a(1),b(1)),c(1)))-1)
      tmin(2) = max(1,floor(min(min(a(2),b(2)),c(2)))-1)
      ! print*, tmax
      ! print*, tmin
      do i=tmin(1),tmax(1)
      do j=tmin(2),tmax(2)
        ip = hoc(i,j)
        do while (ip /= 0)
          if  (point_in_triangle(xp(:,ip), a, b, c)) then
            np_void = np_void + 1
            iq(1)=(ip-1)/ng
            iq(2)=modulo(ip-1,int(ng,4))
            xpos=iq+0.5
            qp_mean = qp_mean + pbc_vec(xpos-center)
            write(13) ip!xp(:,ip)
            write(13) xpos
          endif
          ip = ll(ip)
        enddo
      enddo
      enddo
      inquire(13, pos=current_pos)  ! 获取当前指针位置
      write(13, pos=(current_pos - np_void * 8 * 2 - 8 )) np_void
      qp_mean = qp_mean / np_void+center
      write(13, pos=(current_pos)) qp_mean
      ! print*,center,qp_mean, radius, np_void
      np_void = 0
      write(13) np_void

      n_void = n_void + 1
    endif
    ! if (n_void == 2) then
    !   close(12)
    !   close(13)
    !   stop '停止'
    ! endif
  enddo
  close(12)
  close(13)


  print*, "Output: ", output_name('DT_void_in')
  open(11,file=output_name('DT_void_in'),status='replace',access='stream')
  n_void = 0
  do t = 1, ntri
    a = triangles(:,1,t)
    b = triangles(:,2,t)
    c = triangles(:,3,t)  
    call calculate_incircle(a, b, c, center, radius)
    if (center(1) .lt. 0.0 .or. center(1) .ge. real(ng) .or. &
        center(2) .lt. 0.0 .or. center(2) .ge. real(ng)) then
        ! print*, "Outside: ", t, a, b, c, center, radius
        cycle
    endif
    if (radius .gt. 1.0) then
      n_void = n_void + 1
      write(11) center, radius, a, b, c
    endif
  enddo
  close(11)
  print*, "n_void =", n_void

  
  open(12,file=output_name('DT_void_in_r'),status='replace',access='stream')
  write(12) n_void

  n_void = 0
  do t = 1, ntri
    a = triangles(:,1,t)
    b = triangles(:,2,t)
    c = triangles(:,3,t)  
    call calculate_incircle(a, b, c, center, radius)
    if (center(1) .lt. 0.0 .or. center(1) .ge. real(ng) .or. &
        center(2) .lt. 0.0 .or. center(2) .ge. real(ng)) then
        ! print*, "Outside: ", t, a, b, c, center, radius
        cycle
    endif
    if (radius .gt. 20.0) then
      write(12) center, radius, a, b, c

      ! r
      np_void = 0
      write(12) np_void
      r2 = radius**2
      qp_mean = 0
      cm=xpos2mesh(center,ng)
      rm=floor(rate_in*radius)+1
      do i=-rm,rm
      do j=-rm,rm
        idl(1) = mod(cm(1)+i+ng-1,ng)+1
        idl(2) = mod(cm(2)+j+ng-1,ng)+1
        if (idl(1) < 1 .or. idl(1) > ng .or. idl(2) < 1 .or. idl(2) > ng) then
          print*,idl
          print*,cm
          print*,i,j
          print*,cm(1)+i+ng,cm(2)+j+ng
          print*,mod(cm(1)+i+ng,ng)+1,mod(cm(2)+j+ng,ng)+1
          stop 'out of range'
        endif
        ip = hoc(idl(1),idl(2))
        do while (ip /= 0)
          xpos=pbc_vec(xp(:,ip) - center)
          d2 = sum(xpos**2)

          if (d2 < r2) then
            np_void = np_void + 1
            iq(1)=(ip-1)/ng
            iq(2)=modulo(ip-1,int(ng,4))
            xpos=iq+0.5
            qp_mean = qp_mean + pbc_vec(xpos-center)
            write(12) ip
            write(12) xpos
          endif
          ip = ll(ip)
        enddo
      enddo
      enddo
      inquire(12, pos=current_pos)  ! 获取当前指针位置
      write(12, pos=(current_pos - np_void * 8 * 2 - 8 )) np_void
      qp_mean = qp_mean / np_void + center
      write(12, pos=(current_pos)) qp_mean
      print*,center,qp_mean, radius, np_void
      np_void = 0
      write(12) np_void
    endif

  enddo
  close(12)



  contains

  subroutine calculate_circumcircle(a, b, c, center, radius)
    real, intent(in) :: a(2), b(2), c(2)
    real, intent(out) :: center(2), radius
    
    real :: d, ux, uy, uz, dx, dy
    real :: denom
    
    ! 计算外接圆圆心
    d = 2.0 * (a(1)*(b(2)-c(2)) + b(1)*(c(2)-a(2)) + c(1)*(a(2)-b(2)))
    
    if (abs(d) < 1e-10) then
      ! 三点共线，无法形成外接圆
      center = 0.0
      radius = 0.0
      return
    endif
    
    ux = a(1)*a(1) + a(2)*a(2)
    uy = b(1)*b(1) + b(2)*b(2)
    uz = c(1)*c(1) + c(2)*c(2)
    
    center(1) = (ux*(b(2)-c(2)) + uy*(c(2)-a(2)) + uz*(a(2)-b(2))) / d
    center(2) = (ux*(c(1)-b(1)) + uy*(a(1)-c(1)) + uz*(b(1)-a(1))) / d
    
    ! 计算半径
    dx = a(1) - center(1)
    dy = a(2) - center(2)
    radius = sqrt(dx*dx + dy*dy)
    if (radius**2 > triangle_area(a, b, c)) then
      center = 0.0
      radius = 0.0
      return
    endif
  end subroutine calculate_circumcircle

  subroutine calculate_incircle(a, b, c, center, radius)
    implicit none
    real, intent(in) :: a(2), b(2), c(2)  ! 三角形的三个顶点坐标
    real, intent(out) :: center(2)         ! 内切圆圆心坐标
    real, intent(out) :: radius            ! 内切圆半径
    
    real :: ab_len, bc_len, ca_len         ! 三角形边长
    real :: perimeter                      ! 三角形周长
    real :: area                           ! 三角形面积
    
    ! 计算三角形边长
    ab_len = sqrt((b(1)-a(1))**2 + (b(2)-a(2))**2)
    bc_len = sqrt((c(1)-b(1))**2 + (c(2)-b(2))**2)
    ca_len = sqrt((a(1)-c(1))**2 + (a(2)-c(2))**2)
    
    ! 计算三角形周长
    perimeter = ab_len + bc_len + ca_len
    
    ! 使用海伦公式计算三角形面积
    area = 0.25 * sqrt((ab_len+bc_len+ca_len) * &
                      (-ab_len+bc_len+ca_len) * &
                      (ab_len-bc_len+ca_len) * &
                      (ab_len+bc_len-ca_len))
    
    ! 计算内切圆半径
    radius = 2.0 * area / perimeter
    
    ! 计算内切圆圆心坐标
    center(1) = (a(1)*bc_len + b(1)*ca_len + c(1)*ab_len) / perimeter
    center(2) = (a(2)*bc_len + b(2)*ca_len + c(2)*ab_len) / perimeter
    if (radius**2*8 < triangle_area(a, b, c)) then
      center = 0.0
      radius = 0.0
      return
    endif
end subroutine calculate_incircle



  function point_in_triangle(xp, a, b, c) result(inside)
    implicit none
    real, intent(in) :: xp(2), a(2), b(2), c(2)
    logical :: inside
    
    real :: area, alpha, beta, gamma
    real :: eps = 1e-10
    
    ! 计算整个三角形的面积（有符号面积）
    area = triangle_area(a, b, c)
    
    ! 计算三个子三角形的面积
    alpha = triangle_area(xp, b, c) / area
    beta = triangle_area(a, xp, c) / area
    gamma = triangle_area(a, b, xp) / area
    
    ! 判断点是否在三角形内
    inside = (alpha >= -eps) .and. (beta >= -eps) .and. (gamma >= -eps) .and. &
              (alpha + beta + gamma <= 1.0 + eps)
    
  end function point_in_triangle

  ! 计算三角形有符号面积的辅助函数
  function triangle_area(p1, p2, p3) result(area)
    real, intent(in) :: p1(2), p2(2), p3(2)
    real :: area
    
    area = 0.5 * ((p2(1)-p1(1))*(p3(2)-p1(2)) - (p3(1)-p1(1))*(p2(2)-p1(2)))
  end function triangle_area
  
    
end program DT_void