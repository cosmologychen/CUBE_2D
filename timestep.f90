
subroutine timestep
   use variables
   implicit none
   save
   integer ntry,j,i
   real ra,da_1,da_2,a_next,z_next,ai
   real D_center  ! one_run_lightcone 变量

   dt_old=dt
   sim%timestep=sim%timestep+1
   call tic(1)
   print*, ''
   print*, '-------------------------------------------------------'
   print*, 'timestep    :',sim%timestep
   dt_e=dt_max
   ntry=0
   do
      ntry=ntry+1
      da = expansion(dt_e)
      ra=da/(sim%a+da)
      ! print*,ntry,dt_e,ra,da,sim%a+da
      if (ra>ra_max) then
         dt_e=dt_e*(ra_max/ra)
      else
         exit
      endif
      if (ntry>10) exit
   enddo
   dt = 4*min(dt_e,sim%dt_pm1,sim%dt_pm2,dt_refine*sim%dt_pp,sim%dt_vmax)
   da = expansion(dt)

   checkpoint_step=.false.
   z_next=z_checkpoint(sim%cur_checkpoint)
   a_next=1.0/(1+z_next)
   if (da>=a_next-sim%a) then
      if (z_next==z_checkpoint(sim%cur_checkpoint)) then
         checkpoint_step=.true.
         if (sim%cur_checkpoint==n_checkpoint) final_step=.true.
      endif
      ntry=0
      do while (abs((sim%a+da)/a_next-1)>=1e-6)
         dt=dt*(a_next-sim%a)/da
         da = expansion(dt)
         j = j+1
         if (ntry>10) exit
      enddo
   endif

   ra=da/(sim%a+da)
   a_mid=sim%a+(da/2)

   tcat(41,istep)=sim%a
   tcat(42,istep)=a_mid
   tcat(43,istep)=sim%a+da
   dtau = 0

   ! one_run_lightcone: 更新a_grid
   ! 每个格点的a_grid根据其膨胀历史独立演化
   ! 使用find_a_from_t函数，根据格点当前位置计算新的尺度因子
   D_center = Dgrow(sim%a+da)
   if (one_run_lightcone) then
      !$omp parallel do default(shared) private(j,i)
      do j = 1, ng
         do i = 1, ng
            a_grid(i, j) = find_a_from_t(find_t_from_a(a_grid(i,j)) + dt)
            D_grid(i, j) = Dgrow(a_grid(i, j)) / D_center
         enddo
      enddo
      !$omp end parallel do
   endif

   print*, 'tau         :',sim%tau,sim%tau+dtau
   print*, 'z         :',1.0/sim%a-1.0,1.0/(sim%a+da)-1.0
   print*, 'a         :',sim%a,a_mid,sim%a+da
   if ( one_run_lightcone ) then
      print*, 'a_range   :',sim%a+da,minval(a_grid) / (sim%a+da), maxval(a_grid) / (sim%a+da)
      print*, 'D_range   :',D_center,minval(D_grid),maxval(D_grid)
      print*, 'a_grid    :',a_grid(1,1),a_grid(ng/2,ng/2)
      print*, 'D_grid    :',D_grid(1,1),D_grid(ng/2,ng/2)
   endif

   print*, 'expansion :',ra
   print*, 'dt        :',dt
   print*, 'dt_a      :',merge(100.,1e-2/(sim%a**2),sim%a<1)
   print*, 'dt_e      :',dt_e
   print*, 'dt_pm1    :',sim%dt_pm1
   print*, 'dt_pm2    :',sim%dt_pm2
   print*, 'dt_pp     :',sim%dt_pp
   print*, 'dt_vmax   :',sim%dt_vmax
   print*, 'cur_powerpoint :',sim%cur_powerpoint,z_powerpoint(sim%cur_powerpoint)
   print*, ''
   sim%tau=sim%tau+dtau
   sim%t=sim%t+dt
   sim%a=sim%a+da

   call toc(1)

contains


   real function expansion(dt0)
      use variables, only: stime,s2a,ia,sim,istep_max
      implicit none
      real(8) :: a_x,adot,t_x,tdoa,a8_0
      real(4) :: dt0
      integer i1,i2,il,ir,imid

      a8_0=sim%a

      ! Use binary search to find i1 such that s2a(i1) <= a8_0 < s2a(i1+1)
      il = 1
      ir = istep_max
      do while (ir - il > 1)
         imid = (il + ir) / 2
         if (s2a(imid) < a8_0) then
            ir = imid
         else
            il = imid
         endif
      enddo
      i1 = il
      ia = i1

      tdoa = (stime(i1+1)-stime(i1))/(s2a(i1+1)-s2a(i1))
      t_x = stime(i1)+tdoa*(a8_0-s2a(i1))+dt0

      ! Use binary search to find i2 such that stime(i2) <= t_x < stime(i2+1)
      il = 1
      ir = istep_max
      do while (ir - il > 1)
         imid = (il + ir) / 2
         if (stime(imid) < t_x) then
            ir = imid
         else
            il = imid
         endif
      enddo
      i2 = il

      adot = (s2a(i2+1)-s2a(i2))/(stime(i2+1)-stime(i2))
      a_x = s2a(i2)+adot*(t_x-stime(i2))

      expansion=a_x-a8_0
   endfunction


   real function find_a_from_t(t_in)
      implicit none
      real :: t_in,a1,a2,t_step,t1,t2
      integer il,ir,imid

      t_step = stime(2)
      il = floor(t_in/t_step)+1


      t1 = stime(il)
      t2 = stime(il+1)
      a1 = s2a(il)
      a2 = s2a(il+1)
      ! print*,stime(1:4)
      ! print*,il,t_in,t_step,t_in-t1,t_in-t2,a1,a2
      ! stop

      find_a_from_t = a1 + (a2 - a1) / (t2 - t1) * (t_in - t1)
   endfunction

   real function find_t_from_a(a_in)
      implicit none
      real :: a_in,a1,a2,t1,t2
      integer il,ir,imid

      il = 1
      ir = istep_max
      do while (ir - il > 1)
         imid = (il + ir) / 2
         if (s2a(imid) < a_in) then
            ir = imid
         else
            il = imid
         endif
      enddo

      t1 = stime(il)
      t2 = stime(ir)
      a1 = s2a(il)
      a2 = s2a(ir)

      find_t_from_a = t1 + (t2 - t1) * (a_in - a1) / (a2 - a1)
   endfunction

   function Dgrow(a) ! growth function
      implicit none
      real, parameter :: om=omega_m
      real, parameter :: ol=omega_l
      real :: a
      real :: Dgrow
      real :: g,ga,hsq,oma,ola
      hsq=om/a**3+(1-om-ol)/a**2+ol
      oma=om/(a**3*hsq)
      ola=ol/hsq
      g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
      ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))
      Dgrow=a*ga/g
   end function Dgrow

endsubroutine timestep


