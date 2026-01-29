
#define READ_SEED
#define USE_PKIC
! #define READ_NOISE

program initial_conditions
  use omp_lib
  use parameters
  use iso_fortran_env, only : int64
  implicit none
  save
  include 'fftw3.f'

  integer,parameter :: nk=npbin ! transfer function length
  logical,parameter :: norm_as_sigma8=.false.
  logical,parameter :: write_phik=.false.
  logical,parameter :: write_potential=.true.

  integer(8) plan,iplan

  real        rho1(nw+2,nw)
  complex     rho1k(nw/2+1,nw)
  equivalence(rho1,rho1k)
  

  integer istat
  
  integer nthreads

  real tf(13,nk)
  real kr,kx,ky,kmax,temp_r,temp_theta,pow
  real(8) norm_As,v8,j0kl

  integer(4) seedsize
  integer(4),allocatable :: iseed(:)

  complex     delta_k(nw/2+1,nw)
  real,allocatable :: phi(:,:)

  real grad_max(2),vmax(2),vf
  real diff_x(2,2),diff_y(2,2)
  real(4) svz(500,2),svr(100,2) ! for velocity conversion between integers and reals
  real(8) sigma_vc,sigma_vf ! coarse- and fine- grid based velocity dispersion
  integer(int64) time64


  integer i,j,ii,jj,ilayer,nlayer,itx,ity,imove,idx,iq(2)
  real lg

   !! zip format arrays
  integer,parameter :: n_buffer=1 ! depth in coarse cells
  integer,parameter :: nb=32 ! replace ngb
  integer(8),parameter :: npt=ng ! np (number of particle) / dimension (dim)
  integer(8),parameter :: npb=n_buffer*np_nc ! np / buffer depth
  integer(8),parameter :: npmax=2*(npt+2*npb)**3 ! maximum

  real(4),allocatable :: vp(:,:),xp(:,:)
  real(8) gradphi(2)
  




  call omp_set_num_threads(ncore)
  call system('mkdir -p '//opath)

  sim%cur_checkpoint=1 !! set current checkpoint to 1
  open(16,file='./z_checkpoint.txt',status='old') !! open redshift list to do checkpoint
  read(16,fmt='(f8.4)') z_checkpoint(sim%cur_checkpoint)
  close(16)

  print*, ''
  print*, 'CUBE_2D IC run on',int(ncore,1),'cores'
  print*, 'IC resolution: ngic   =', int(ngic,2)
  print*, 'CUBE resolution: ng   =', int(ng,2),npbin
  print*, 'To genterate np       =', int(ng,2),'^2'
  print*, 'at redshift           =', z_checkpoint(sim%cur_checkpoint)
  print*, 'Box size              =', box
  print*, 'output                :', opath
  ! print*, 'body_centered_cubic   =',body_centered_cubic
  print*, '-----------------------------------------'
  call system('mkdir -p '//opath//'code')
  call system('cp -r ./*.f* '//opath//'code/')
  call system('cp ./z_*.txt '//opath//'code/')
  print*, 'export PYTHONUNBUFFERED=1'
  call system('export PYTHONUNBUFFERED=1')
  print*, 'python ./neutrinos/Pk.py'
  call system('~/anaconda3/envs/camb/bin/python ./neutrinos/Pk.py')
  call system('cp /Users/chenbinghang/Library/CloudStorage/OneDrive-个人/陈冰航/天文/sim_2D/test2048/*seed* '//opath)

  ! initialize the simulation information in 'sim' structure
  sim%np=0 ! local (on this image) CDM particle number
  sim%a=1./(1+z_checkpoint(sim%cur_checkpoint)) ! scale factor of the universe
  sim%t=0 ! time
  sim%tau=0 ! conformal time

  sim%timestep=0 ! number of timesteps
  sim%dt_pm1=9999;  sim%dt_pm2=9999; sim%dt_pp=9999;  sim%dt_vmax=9999

  sim%cur_checkpoint=1 ! current checkpoint
  sim%box=box ! box size in Mpc/h
  sim%nnt=nnt ! number of tiles / image / dim
  sim%nt=nt ! number of coarse grid / tile / dim
  sim%ncb=ncb ! buffer depth in coarse cell; always 6
  sim%izipx=izipx ! integer format for CDM position
  sim%izipv=izipv ! integer format for CDM velocity

  sim%h0=h0 ! Hubble
  sim%omega_m=omega_m
  sim%omega_l=omega_l
  sim%s8=s8 ! sigma_8
  sim%vsim2phys=(140./sim%a)*box*sqrt(omega_m)/ng ! velocity unit
  sim%z_i=z_checkpoint(sim%cur_checkpoint) ! initial redshift
  sim%Mass_nu=Mass_nu
  sim%m_nu=m_nu
  sim%cur_powerpoint=1
  sim%calculate_PK=calculate_PK


  tf=0
  lg = real(nw)/real(box)
  print*,lg
  ! stop
  

  call sfftw_init_threads(istat)
  print*, '    sfftw_init_threads status',istat
  nthreads=omp_get_max_threads()
  print*, '    omp_get_max_threads() =',nthreads
  call sfftw_plan_with_nthreads(nthreads)


  ! transferfnc --------------------------------------
  print*,''
  print*,'Transfer function'
  call system_clock(t1,t_rate)


# ifdef USE_PKIC
   print*,'    read Pk_ic from '//nupath//'IC/Pcb_ic.txt'
   open(11,file=nupath//'IC/Pcb_ic.txt',form='formatted')
   read(11,*) tf(:2,:npbin)
   close(11)
   
   tf(2,:)= sqrt(tf(2,:)) * lg
# else
   write(str_z,'(f8.4)') z_checkpoint(sim%cur_checkpoint)
   print*,'    read TF_ic from IC/IC_rescaled_transfer_z'//trim(adjustl(str_z))//'.txt'
   open(11,file='IC/IC_rescaled_transfer_z'//trim(adjustl(str_z))//'.txt',form='formatted')
   read(11,*) tf
   close(11)
   print*,tf(:2,:)
   tf(2,:)= tf(1,:)**(3+n_s) * (tf(2,:)**2)**(2/3) / (2*pi**2)
   if (norm_as_sigma8) then
      print*,'  scale amplitude as sigma_8'
      print*,'  sigma_8 =',s8
      v8=0; tf(4,1)=tf(1,2)/2 ! calculate v8
      do k=2,nk-1
         tf(4)=(tf(1+1)-tf(1-1))/2
      enddo
      tf(4,nk)=tf(1,nk)-tf(1,nk-1)
      kmax=2*pi*sqrt(3.)*nyquist/box
      do k=1,nk
         if (tf(1)>kmax) exit
         
         v8=v8+tf(2)*tf(4)/tf(1)*merge(1d0,3*(sin((tf(1)*8)*1d0)-cos((tf(1)*8)*1d0)*(tf(1)*8)*1d0)/((tf(1)*8)*1d0)**3,(tf(1)*8)==0)**2
      enddo
      print*, '  v8, (s8^2/v8) =', v8, s8**2/v8
      tf(2,:)=tf(2,:)*(s8**2/v8)*Dgrow(sim%a)**2
      then
         open(11,file=output_dir()//'tf_s8.txt',status='replace',access='stream')
         write(11) tf
         close(11)
      endif
   else
      print*,'  scale as A_s'
      !norm_As=2.*pi**2*h0**4*(h0/0.05)**(n_s-1)
      norm_As=3.878
      print*,'  norm_As =',norm_As
      print*,'  A_s =',A_s
      tf(2,:)=tf(2,:)*A_s*norm_As!*Dgrow(sim%a)**2
      then
         open(11,file=output_dir()//'tf_as.txt',status='replace',access='stream')
         write(11) tf(:2,:162)
         close(11)
      endif
   endif
# endif

  ! noisemap -------------------------------------
  print*,''
  print*,'Generating random noise'
  call random_seed(size=seedsize) ! generate random seed
  print*,'  min seedsize =', seedsize
  seedsize=max(seedsize,12)
  allocate(iseed(seedsize))
#ifdef READ_SEED
   print*, '  Copy and read seeds from ../confings/'
   open(11,file=output_name('seed'),status='old',access='stream')
   read(11) iseed
   close(11)
   ! Input iseed
   call random_seed(put=iseed)
   print*, 'iseed', iseed
#else
# ifdef READ_NOISE
   open(11,file=output_name('noise'),access='stream')
   read(11) rho1(1:ngic,1:ngic)
   close(11)
   print*, '  READ IN NOISE MAP:', rho1(1,1), rho1(ngic,ngic)
# else
   ! Generate at least 12 seeds according to system clock
   call system_clock(time64)
   do i = 1, seedsize
      iseed(i) = lcg(time64) + 137
   enddo
   call random_seed(put=iseed) ! generate random seed
   open(11,file=output_name('seed'),status='replace',access='stream')
   write(11) iseed
   close(11)
# endif
#endif
  deallocate(iseed)

# ifndef READ_NOISE
    call random_number(rho1(1:ngic,1:ngic)) ! generate random numbers ! rho1 is defined in pencil_fft ! tophat [0,1) distribution
    open(11,file=output_name('noise'),status='replace',access='stream')
    write(11) rho1(1:ngic,1:ngic) ! write random numbers into a file
    close(11)
    print*, '  noise',int(image,1),rho1(1:2,1)
# endif

  ! Box-Muller transform ----------------------------------------------
  ! convert to standard Normal (Gaussian) distribution
  print*,'  Box-Muller transform'
  
  !$omp paralleldo&
  !$omp& default(shared) &
  !$omp& private(j,i,temp_theta,temp_r)
  do j=1,ngic
      do i=1,ngic,2
        temp_theta=2*pi*rho1(i,j)
        temp_r=sqrt(-2*log(1-rho1(i+1,j)))
        rho1(i  ,j)=temp_r*cos(temp_theta)
        rho1(i+1,j)=temp_r*sin(temp_theta)
      enddo
  enddo
  !$omp endparalleldo
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';

  ! rho1 = 0
  ! rho1(ngic/2,ngic/2) = 100
  ! ! rho1(ngic/4,ngic/4) = 100
  ! ! rho1(ngic/4,ngic/4*3) = 100
  ! ! rho1(ngic/4*3,ngic/4) = 100
  ! ! rho1(ngic/4*3,ngic/4*3) = 300

  ! delta_field ----------------------------------------------------
  ! apply transfer function in Fourier space (Wiener filter) and get delta_L
  print*, ''
  print*, 'delta field'
  call system_clock(t1,t_rate)
  print*, '  ftran'
  call sfftw_plan_dft_r2c_2d( plan,ngic,ngic,rho1k,rho1k,FFTW_MEASURE)
  call sfftw_plan_dft_c2r_2d(iplan,ngic,ngic,rho1k,rho1k,FFTW_MEASURE)
  print*, rho1(1:6,1),maxval(rho1),minval(rho1)
  call sfftw_execute( plan) ! Fourier transform
  ! print*, rho1k(1:6,1)


  print*, '  Wiener filter'
  !$omp paralleldo&
  !$omp& default(shared) &
  !$omp& private(j,i,kx,ky,kr)
  do j=1,ngic
      do i=1,nyquist+1
        ky=mod(j+nyquist-1,nw)-nyquist
        kx=i-1
        kr=sqrt(kx**2+ky**2)*(2*real(pi,kind=4)/box) ! kr is |k_n|
        rho1k(i,j)=rho1k(i,j) * interp_tf(kr)
      enddo
  enddo
  !$omp endparalleldo
  rho1k(1,1)=0 ! zero frequency
  print*,('   Wiener filter done')
  print*, rho1k(1:6,1)
  delta_k=rho1k ! backup (Fourier) k-space delta_L


  print*,'  btran'
  call sfftw_execute(iplan) ! inverse Fourier transform
  rho1 = rho1/real(ngic*ngic)
  print*, rho1(1:6,1),maxval(rho1)/Dgrow(sim%a),minval(rho1)/Dgrow(sim%a)


  ! rho1 = 0
  ! rho1(ngic/2,ngic/2) = 100

  print*,'  delta_L',rho1(1:4,1)
  print*,'  rms of delta',sqrt(sum(rho1**2*1.d0)/ngic/ngic)

  print*,'  write delta_L into file',output_name('delta_L')
  print*,'  growth factor Dgrow(',sim%a,') =',Dgrow(sim%a) ! growth factor
  open(11,file=output_name('delta_L'),status='replace',access='stream')
  do i=1,ngic
    write(11) rho1(1:ngic,i)/Dgrow(sim%a) ! write layer by layer to avoid bug
  enddo
  close(11)


  print*,'  write delta_ic into file',output_name('delta_ic')
  open(11,file=output_name('delta_ic'),status='replace',access='stream')
  do i=1,ngic
    write(11) rho1(1:ngic,i) ! write layer by layer to avoid bug
  enddo
  close(11)

  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';



  ! print*, ''
  ! print*, 'Potential field',nyquist
  ! call system_clock(t1,t_rate)
  ! !$omp paralleldo&
  ! !$omp& default(shared) &
  ! !$omp& private(j,i,ky,kx,kr,j0kl)
  ! do j=1,ngic
  !     ky=mod(j+nyquist-1,ngic)-nyquist
  !     do i=1,nyquist+1
  !       kx=i-1
  !       kr=kx**2+ky**2
  !       kr=max(kr,1.0) * ((2 * pi / box )**2)! avoid kr being 0
  !       rho1k(i,j)=-2*pi/kr
  !       ! call Bessel(0,sqrt(kr)*ngic,j0kl)
  !       ! rho1k(i,j)=-2*pi/kr*(1-j0kl) ! raw potential kernel
  !     enddo
  ! enddo
  ! !$omp endparalleldo
  ! rho1k(1,1)=0 ! zero frequency
  ! ! stop

  ! print*, ''
  ! print*, 'Potential field'
  ! call system_clock(t1,t_rate)
  ! !$omp paralleldo&
  ! !$omp& default(shared) &
  ! !$omp& private(j,i,ky,kx,kr,j0kl)
  ! do j=1,ngic
  !     ky=mod(j+nyquist-1,ngic)-nyquist
  !     ky=2*sin(pi*ky/nw)
  !     do i=1,nyquist+1
  !       kx=i-1
  !       kx=2*sin(pi*kx/nw)
  !       kr=kx**2+ky**2
  !       kr=max(kr,1.0/nw**2)! avoid kr being 0
  !       call Bessel(0,sqrt(kr)*ngic,j0kl)
  !       rho1k(i,j)=-2*pi/kr*(1-j0kl) ! raw potential kernel
  !     enddo
  ! enddo
  ! !$omp endparalleldo
  ! rho1k(1,1)=0 ! zero frequency

  print*, ''
  print*, 'Potential field'
  call system_clock(t1,t_rate)
  !$omp paralleldo&
  !$omp& default(shared) &
  !$omp& private(j,i,ky,kx,kr)
  do j=1,ngic
      do i=1,nyquist+1
        ky=mod(j+nyquist-1,ngic)-nyquist
        kx=i-1
        ky=2*sin(pi*ky/nw)
        kx=2*sin(pi*kx/nw)
        kr=kx**2+ky**2
        kr=max(kr,1.0/nw**2) ! avoid kr being 0
        rho1k(i,j)=-4*pi/kr ! raw potential kernel
      enddo
  enddo
  !$omp endparalleldo
  rho1k(1,1)=0 ! zero frequency
  ! rho1k=1
   
  rho1k=real(rho1k)*delta_k ! phi(k) = kernel(k) * delta_L(k)
  if (write_phik) then
    print*, '  write phi1 into file'
    open(11,file=output_name('phik'),status='replace',access='stream')
    write(11) rho1k
    close(11)
  endif
  call sfftw_execute(iplan) ! inverse Fourier transform
  rho1 = rho1/(ngic*ngic)

  allocate(phi(-nb:ngic+nb+1,-nb:ngic+nb+1))
  phi = 0
  phi(1:ngic,1:ngic)=rho1(1:ngic,1:ngic) ! phi1
  print*,'  phi',phi(1:4,1)
  if (write_potential) then
    print*, '  write phi1 into file'
    open(11,file=output_name('phi1'),status='replace',access='stream')
    write(11) rho1(1:ngic,1:ngic)
    close(11)
  endif

  ! buffer phi ---------------------------------------------------
  print*, '  buffer phi'
  phi(:0,:)=phi(ngic-nb:ngic,:)
  phi(ngic+1:,:)=phi(1:nb+1,:)
  phi(:,:0)=phi(:,ngic-nb:ngic)
  phi(:,ngic+1:)=phi(:,1:nb+1)

  print*, '  destroying FFT plans'
  call sfftw_destroy_plan( plan)
  call sfftw_destroy_plan(iplan); 
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs'; 


  ! zip checkpoints ------------------------------------------------
  print*, ''
  print*, 'zip checkpoints'
  vf=vfactor(sim%a) ! velocity factor
  print*, '  vf =',vf
  !! maximum gradient of phi
  grad_max(1)=maxval(abs(phi(-nb:ngic+nb-1,:)-phi(-nb+2:ngic+nb+1,:)))
  grad_max(2)=maxval(abs(phi(:,-nb:ngic+nb-1)-phi(:,-nb+2:ngic+nb+1)))
  vmax=grad_max/(4*nic*pi)*vf ! maximum velocity
  sim%dt_vmax=vbuf*16./maxval(abs(vmax)) ! constrain dt by maximum velocity
  sim%vz_max=vmax(2)
  nlayer=2*ceiling(grad_max(2)*np_nc/(4*nic*pi*ratio_cs))+3 ! for OpenMP
  
  print*, '  grad_max',grad_max
  print*, '  max dsp',grad_max/(4*nic*pi)
  print*, '  vmax',vmax
  print*, '  vz_max',sim%vz_max
  if (maxval(grad_max)/(4*nic*pi)>=nb) then
      print*, '  particle dsp > buffer' ! particle might move beyond buffer depth
      print*, maxval(grad_max)/(4*nic*pi),nb
      stop
  endif
  print*,'  Thread save nlayer =',nlayer
  print*,''

  ! velocity conversion as a function of redshift and scale
  open(11,file='./velocity_conversion/sigmav_z.bin',access='stream')
  read(11) svz
  close(11)
  open(11,file='./velocity_conversion/sigmav_r.bin',access='stream')
  read(11) svr
  close(11)

  sigma_vf=interp_sigmav(sim%a,box/ng) ! sigma(v) on scale of fine grid, in km/s
  sigma_vc=interp_sigmav(sim%a,box/nc) ! sigma(v) on scale of coarse grid, in km/s
  sim%sigma_vres=sqrt(sigma_vf**2-sigma_vc**2) ! sigma(v) residual, in km/s
  sim%sigma_vi=sim%sigma_vres/sim%vsim2phys/sqrt(2.) ! sigma(v_i) residual, in sim unit

  print*, ''
  print*, 'Read velocity dispersion prediction'
  print*,'sigma_vf(a=',sim%a,', r=',box/ng,'Mpc/h)=',real(sigma_vf,4),'km/s'
  print*,'sigma_vc(a=',sim%a,', r=',box/nc,'Mpc/h)=',real(sigma_vc,4),'km/s'
  print*,'sigma_vres=',real(sim%sigma_vres,4),'km/s'
  print*,'sigma_vi =',real(sim%sigma_vi,4),'(simulation unit)'


  ! create particles (no communication) ----------------------------
  print*,''
  print*, 'Create particles'
  diff_x(1,:)=-0.5; diff_x(2,:)=0.5
  diff_y(:,1)=-0.5; diff_y(:,2)=0.5
  call system_clock(t1,t_rate)

  allocate(xp(2,ng**2),vp(2,ng**2))

  !$omp paralleldo default(shared) num_threads(ncore) schedule(dynamic,1) private(j,i,jj,ii,gradphi,idx)
  do j=1,ng
    jj = j * nic-1
    do i=1,ng
      ii = i * nic-1
      gradphi(1)=sum(phi(ii-1:ii,jj-1:jj)*diff_x)
      gradphi(2)=sum(phi(ii-1:ii,jj-1:jj)*diff_y)

      idx = (i-1)*ng+j
      xp(1,idx)=wrap_position(real(i-0.5-gradphi(1)/(4*nic*pi),kind=4))
      xp(2,idx)=wrap_position(real(j-0.5-gradphi(2)/(4*nic*pi),kind=4))

      vp(:,idx)=-gradphi/(4*pi*nic)*vf
    enddo
  enddo ! j
  !$omp endparalleldo

  sim%np=ng**2

  open(11,file=output_name('xp'),status='replace',access='stream') ! position list
  write(11) xp
  close(11)

  open(12,file=output_name('vp'),status='replace',access='stream') ! velocitie list
  write(12) vp
  close(12)
  call system_clock(t2,t_rate)
  

  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';

  print*,'np',sim%np
  sim%mass_p_cdm=real(ng**2,kind=8)/sim%np ! particle mass in unit of fine cell
  print*,'sim%mass_p_cdm =',sim%mass_p_cdm
  

  call print_header(sim)

  
  open(10,file=output_name('info'),status='replace',access='stream')
  write(10) sim
  close(10)
  

  call system_clock(ttt2,t_rate)
  print*, 'total elapsed time =',real(ttt2-ttt1)/t_rate,'secs';
  print*, 'initial condition done'

contains

  real function interp_sigmav(aa,rr)
    implicit none
    integer(8) ii,i1,i2
    real aa,rr,term_z,term_r
    i1=1
    i2=500
    do while (i2-i1>1)
        ii=(i1+i2)/2
        if (aa>svz(ii,1)) then
          i1=ii
        else
          i2=ii
        endif
    enddo
    term_z=svz(i1,2)+(svz(i2,2)-svz(i1,2))*(aa-svz(i1,1))/(svz(i2,1)-svz(i1,1))
    i1=1
    i2=100
    do while (i2-i1>1)
        ii=(i1+i2)/2
        if (rr>svz(ii,1)) then
          i1=ii
        else
          i2=ii
        endif
    enddo
    term_r=svr(i1,2)+(svr(i2,2)-svr(i1,2))*(rr-svr(i1,1))/(svr(i2,1)-svr(i1,1))
    interp_sigmav=term_z*term_r
  endfunction

  real function interp_tf(kr) ! interpolation in log space
    implicit none
    integer(8) ii,i1,i2
    real kr,xx,yy,x1,x2,y1,y2
    i1=1
    i2=nk
    do while (i2-i1>1)
        ii=(i1+i2)/2
        if (kr>tf(1,ii)) then
          i1=ii
        else
          i2=ii
        endif
    enddo
    x1=log(tf(1,i1))
    y1=log(tf(2,i1))
    x2=log(tf(1,i2))
    y2=log(tf(2,i2))
    xx=log(kr)
    yy=y1+(y2-y1)*(xx-x1)/(x2-x1)
    interp_tf=exp(yy)
  endfunction

  function Dgrow(a) ! growth function
    implicit none
    real, parameter :: om=omega_m
    real, parameter :: ol=omega_l
    real a
    real Dgrow
    real g,ga,hsq,oma,ola
    hsq=om/a**3+(1-om-ol)/a**2+ol
    oma=om/(a**3*hsq)
    ola=ol/hsq
    g=2.5*om/(om**(4./7)-ol+(1+om/2)*(1+ol/70))
    ga=2.5*oma/(oma**(4./7)-ola+(1+oma/2)*(1+ola/70))
    Dgrow=a*ga/g
  end function Dgrow

  real function vfactor(a)
    implicit none
    real :: a
    real :: H,km,lm
    lm=omega_l/omega_m
    km=(1-omega_m-omega_l)/omega_m
    H=2/(3*sqrt(a**3))*sqrt(1+a*km+a**3*lm)
    vfactor=a**2*H
  endfunction

  function lcg(s) ! Linear congruential generator
    implicit none
    integer(4) :: lcg
    integer(int64) :: s
    if (s == 0) then
        s = 104729
    else
        s = mod(s, 4294967296_int64)
    end if
    s = mod(s * 279470273_int64, 4294967291_int64)
    lcg = int(mod(s, int(huge(0), int64)), kind(0))
  endfunction

  subroutine Bessel(a,x,j)
    use parameters
    use, intrinsic :: ieee_arithmetic
    implicit none
    real, intent(in)  :: x
    integer, intent(in)  :: a
    real(8), intent(out) :: j
    real(16) f,ti
    integer(8) mi,mn

    mn = 1000

      
    f = 2*pi/mn
    j=0
    ti=0
    do mi = 1, mn
      ti = ti+f
      j = j + cos(a*ti-x*sin(ti))*f
    end do
    j = j/2/pi
    if (.not. ieee_is_finite(j)) then
      print *, 'j 是 Inf 或 NaN'
      print *, 'a = ', a
      print *, 'x = ', x
      stop
    endif
  endsubroutine
end
