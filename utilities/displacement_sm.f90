#define Cpower
program displacement
  use omp_lib
  use parameters
  use powerspectrum
  implicit none
  save
  include 'fftw3.f'

  integer,parameter,dimension(6) :: smooth_ks=[200,20,15,13,12,10]*3
  integer :: ip,iq(2),i_dim,idx1(2),idx2(2),k_s
  real :: dx1(2),dx2(2),pdim(2)
  integer(8) np,istat,plan,iplan

  integer :: i,j,iteam,cur_checkpoint,nthreads
  real :: pos0(2),pos1(2),dpos(2),kx(2),xi(10,0:nbin)
  real,allocatable :: rho_grid(:,:,:)
  real,allocatable :: dsp(:,:,:),xp(:,:)
  complex,allocatable :: cdiv(:,:),cphi(:,:),rhok_L(:,:),rhok_R(:,:),rhok_N(:,:),grid_sk(:,:)
  real(8) xpos

  print*,64-nearest(real(64, kind(xpos)), -1.0)
  print*, 'Displacement field analysis on resolution:'
  print*, 'ng=',ng
  print*, 'checkpoint at:'
  open(16,file='./z_checkpoint.txt',status='old')
  do i=1,nmax_redshift
    read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
    print*, z_checkpoint(i)
  enddo
  71 n_checkpoint=i-1
  close(16)
  print*,''
  nthreads=omp_get_max_threads()
  print*, '    omp_get_max_threads() =',nthreads
  call omp_set_num_threads(nthreads)

  call sfftw_init_threads(istat)
  print*, '    sfftw_init_threads status',istat
  call sfftw_plan_with_nthreads(nthreads)
  call sfftw_plan_dft_r2c_2d( plan,ng,ng,rho1k,rho1k,FFTW_MEASURE)
  call sfftw_plan_dft_c2r_2d(iplan,ng,ng,rho1k,rho1k,FFTW_MEASURE)


  ! do cur_checkpoint= 1,6
  do cur_checkpoint= 5,5!n_checkpoint,n_checkpoint !2,-1
    print*, ''
    print*,'==========================================================='

    print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))

    k_s = smooth_ks(cur_checkpoint)
    allocate(grid_sk(ng/2+1,ng))

    grid_sk(1,1) = 1
    !$omp paralleldo default(shared) &
    !$omp& private(i,j,kx,pdim)
    do j=1,ng
    do i=1,ng/2+1
      if (j == 1 .and. i == 1) cycle
      kx=modulo([i,j]+ng/2-1,ng)-ng/2 !k
      grid_sk(i,j)=exp(-0.5D0 * sum(kx**2) / (k_s**2))
    enddo
    enddo
    !$omp endparalleldo
    sim%cur_checkpoint=cur_checkpoint

    ! call dsp_halo()
    ! call Divergence('sD')
    call Divergence('H')
    ! call decompose_Mesh_D('sD')
    ! call decompose_Mesh_D('sH')
    ! call dep2delta_e('sD','E_sD',4,.false.,istat)
    call dep2delta_e('H','E_H',4,.false.,istat)
    stop
    call Divergence('sD')
    ! call decompose_Mesh_FFT('sD')
    call dep2delta_e('sD','k_sD',4,.false.,istat)
    stop
    call gassim_file('u_D_q')
    call Divergence('D')
    call Divergence('sD')
    call dep2delta_e('sD','k_sD',4,.false.,istat)
    ! call dep2delta_e('sD','E_D',4,.false.,istat)
    call dep2delta_e('sD','E_sD',4,.false.,istat)

    stop



    open(11,file=output_name('info'),access='stream'); read(11) sim; close(11)
    np = sim%np
    allocate(xp(2,sim%np))
    open(11,file=output_name('xp'),access='stream'); read(11) xp(:,:sim%np); close(11)

    allocate(dsp(2,ng,ng))
    dsp=0

    !$omp paralleldo default(shared) &
    !$omp& private(ip,iq,pos0,pos1,dpos)
    do ip=1,np
      iq(1)=(ip-1)/ng
      iq(2)=modulo(ip-1,int(ng,4))
      pos0=iq+0.5
      pos1=xp(:,ip)
      dpos=pos1-pos0
      dpos=modulo(dpos+ng/2,real(ng))-ng/2
      dsp(:,iq(1)+1,iq(2)+1)=dpos
    enddo
    !$omp endparalleldo
    deallocate(xp)

    do i_dim=1,2
      print*, '   dsp: dimension',int(i_dim,1),'min,max values ='
      print*, '   ',minval(dsp(i_dim,:,:)), maxval(dsp(i_dim,:,:))
    enddo

    print*,'    Write dsp into file:'
    print*,'      save:',output_name('dsp_D')
    open(15,file=output_name('dsp_D'),status='replace',access='stream')
    write(15) dsp
    close(15)


    print*,''
    print*,'    Smoothing dsp'




    do i_dim=1,2
      call gassim(dsp(i_dim,1:ng,1:ng))
      print*, '   dsp: dimension',int(i_dim,1),'min,max values ='
      print*, '   ',minval(dsp(i_dim,:,:)), maxval(dsp(i_dim,:,:))
    enddo ! i_dim
    
    print*,'    Write Smoothed dsp into file:'
    print*,'      save:',output_name('dsp_sD')
    open(15,file=output_name('dsp_sD'),status='replace',access='stream')
    write(15) dsp
    close(15)
    deallocate(dsp)
    ! stop

    call Divergence('D')
    call Divergence('sD')
    
    call decompose_Mesh_D('D')
    call decompose_Mesh_D('sD')
    call dep2delta_e('D','A_D',4,.false.,istat)
    call dep2delta_e('sD','l1_sD',4,.true.,istat)
    call dep2delta_e('sD','l2_sD',4,.true.,istat)
    call dep2delta_e('sD','A_sD',4,.false.,istat)
    ! call decompose_Mesh_FFT('D')

    deallocate(grid_sk)

  enddo
  print*,'displacement done'

  contains

  subroutine Divergence(dsp_name)
    implicit none
    character(len=*), intent(in) :: dsp_name
    integer :: i_dim
    real,allocatable :: dsp(:,:,:)
    complex,allocatable :: cdiv(:,:)

    allocate(dsp(2,ng,ng))
    dsp = 0
    print*,'  read:'
    print*,'    ',output_name('dsp_'//dsp_name)
    open(15,file=output_name('dsp_'//dsp_name),status='old',access='stream')
    read(15) dsp
    close(15)

    print*,''
    print*,'    Start computing Divergence '//dsp_name
    allocate(cdiv(ng/2+1,ng))
    cdiv=0
    do i_dim=1,2
      rho1(1:ng,1:ng)=dsp(i_dim,1:ng,1:ng)
      call sfftw_execute( plan) ! Fourier transform

      !$omp paralleldo default(shared) &
      !$omp& private(i,j,kx,pdim)
      do j=1,ng
      do i=1,ng/2+1
        if (j == 1 .and. i == 1) cycle
        kx=modulo([i,j]+ng/2-1,ng)-ng/2 !k
        pdim=sin(2*pi*kx/ng)
        cdiv(i,j)=cdiv(i,j)+(0,1)*rho1k(i,j)*pdim(i_dim) !c means complex 
      enddo
      enddo
      !$omp endparalleldo
    enddo ! i_dim

    cdiv(1,1)=0

    rho1k=cdiv

    call sfftw_execute(iplan)
    rho1(1:ng,1:ng) = -real(rho1(1:ng,1:ng))/real(ng**2)
    print*,''
    print*,'    write Divergence '
    print*,'    ',minval(rho1(1:ng,1:ng)),maxval(rho1(1:ng,1:ng)),sum(rho1(1:ng,1:ng)*1d0)
    open(15,file=output_name('E_'//dsp_name//'_q'),status='replace',access='stream')
    write(15) rho1(1:ng,1:ng)
    close(15)
  endsubroutine

  subroutine dep2delta_e(dsp_name,rho_name,n_min,sm,state)
    implicit none
    logical, intent(in) :: sm
    integer, intent(in)  :: n_min     ! 至少需要 n 个非零邻居
    character(len=*), intent(in) :: dsp_name,rho_name
    character(len=100) :: out_name
    real rhoe(ng,ng)
    integer(8) , intent(inout) :: state

    integer :: ilayer, i, j, i_dim, j_dim,i_n(4),j_n(4), di, dj, ii, jj, n, count_bad, count_bad_prev,c_count
    real :: sum_d
    real(8) :: A(2,2)
    real,allocatable :: A_mesh(:,:,:,:),delta(:,:),rho_grid(:,:,:),dsp(:,:,:),dsp_t(:,:,:)




    if (n_min > 8 .or. n_min < 2) then
        print *, "bad n_min = ",n_min, "exiting subroutine."
        return
    end if

    print*,''
    print*, 'CIC interpolation'//rho_name//' by dsp_'//dsp_name

    allocate(dsp(2,ng,ng))
    print*,size(rho1k(:,1)),size(rho1k(1,:))
    dsp = 0
    print*,'  read:'
    print*,'    ',output_name('dsp_'//dsp_name)
    print*,'    ',output_name(rho_name//'_q')
    open(15,file=output_name('dsp_'//dsp_name),status='old',access='stream')
    read(15) dsp
    close(15)

    rho1=0
    open(15,file=output_name(rho_name//'_q'),status='old',access='stream')
    read(15) rho1(1:ng,1:ng)
    close(15)




    allocate(delta(ng, ng),rho_grid(0:ng+1,0:ng+1,2*nthreads))
    rho_grid=0
    print*,'init '
    !$omp paralleldo default(shared) &
    !$omp& private(i,j,iteam,pos1,idx1,idx2,dx1,dx2,A,j_n,i_n)
    do j=1,ng
      iteam=omp_get_thread_num()+1
      ! print*,iteam,' ',j
      if (iteam > nthreads) error stop 'thread number out of range'
      do i=1,ng
        pos1=[i,j]-0.5+dsp(:,i,j)
        pos1=wrap_position2(pos1)
        idx1=floor(pos1)+1
        idx2=idx1+1
        dx1=idx1-pos1
        dx2=1-dx1

        if (maxval(idx1) > ng  .or. minval(idx1) < 0) then
          print*, 'xp out of range in kick ll',i,j
          print*, 'dsp=',dsp(:,i,j)
          print*, 'xpos =',pos1
          print*, 'fxpos =',floor(pos1)
          print*, 'idl =',idx1
          stop 
        endif

        rho_grid(idx1(1),idx1(2),iteam) = rho_grid(idx1(1),idx1(2),iteam)+dx1(1)*dx1(2)*rho1(i,j)
        rho_grid(idx2(1),idx1(2),iteam) = rho_grid(idx2(1),idx1(2),iteam)+dx2(1)*dx1(2)*rho1(i,j)
        rho_grid(idx1(1),idx2(2),iteam) = rho_grid(idx1(1),idx2(2),iteam)+dx1(1)*dx2(2)*rho1(i,j)
        rho_grid(idx2(1),idx2(2),iteam) = rho_grid(idx2(1),idx2(2),iteam)+dx2(1)*dx2(2)*rho1(i,j)


        rho_grid(idx1(1),idx1(2),iteam+nthreads) = rho_grid(idx1(1),idx1(2),iteam+nthreads)+dx1(1)*dx1(2)
        rho_grid(idx2(1),idx1(2),iteam+nthreads) = rho_grid(idx2(1),idx1(2),iteam+nthreads)+dx2(1)*dx1(2)
        rho_grid(idx1(1),idx2(2),iteam+nthreads) = rho_grid(idx1(1),idx2(2),iteam+nthreads)+dx1(1)*dx2(2)
        rho_grid(idx2(1),idx2(2),iteam+nthreads) = rho_grid(idx2(1),idx2(2),iteam+nthreads)+dx2(1)*dx2(2)
      enddo
    enddo
    !$omp endparalleldo
    deallocate(dsp)

    delta = 0
    rhoe  = 0
    do iteam=1,nthreads
      rho_grid(1 ,:,iteam) = rho_grid(1 ,:,iteam) + rho_grid(ng+1,:,iteam)
      rho_grid(ng,:,iteam) = rho_grid(ng,:,iteam) + rho_grid(0   ,:,iteam)
      rho_grid(: ,1,iteam) = rho_grid(: ,1,iteam) + rho_grid(:,ng+1,iteam)
      rho_grid(:,ng,iteam) = rho_grid(:,ng,iteam) + rho_grid(:   ,0,iteam)

      rho_grid(1 ,:,iteam+nthreads) = rho_grid(1 ,:,iteam+nthreads) + rho_grid(ng+1,:,iteam+nthreads)
      rho_grid(ng,:,iteam+nthreads) = rho_grid(ng,:,iteam+nthreads) + rho_grid(0   ,:,iteam+nthreads)
      rho_grid(: ,1,iteam+nthreads) = rho_grid(: ,1,iteam+nthreads) + rho_grid(:,ng+1,iteam+nthreads)
      rho_grid(:,ng,iteam+nthreads) = rho_grid(:,ng,iteam+nthreads) + rho_grid(:   ,0,iteam+nthreads)

      rhoe(1:ng,1:ng)  = rhoe(1:ng,1:ng) + rho_grid(1:ng,1:ng,iteam)
      delta(1:ng,1:ng) = delta(1:ng,1:ng) + rho_grid(1:ng,1:ng,iteam+nthreads)
    enddo
    deallocate(rho_grid)


    print*,''
    print*,'    delta : ',minval(delta-1),maxval(delta-1),sum((delta-1)*1d0)
    print*,'      save:',output_name('delta_c'//dsp_name)
    open(16,file=output_name('delta_c'//dsp_name),status='replace',access='stream')
    write(16) delta-1
    close(16)

    if (sm) then 
      call gassim(rhoe)
      call gassim(delta)
      out_name='s'//rho_name
    else
      out_name=rho_name
    endif

    print*,'    rho_c: ',minval(rhoe),maxval(rhoe),sum(rhoe*1d0)/ng/ng
    print*,'      save:',output_name(trim(adjustl(rho_name))//dsp_name//'_x')
    open(16,file=output_name(trim(adjustl(rho_name))//dsp_name//'_x'),status='replace',access='stream')
    write(16) rhoe
    close(16)
    count_bad = 0
    do i = 1, ng
    do j = 1, ng
      if (delta(i,j) < 1e-4) then
        count_bad = count_bad + 1
        rhoe(i,j) = 2e10
      else
        rhoe(i,j) = rhoe(i,j)/delta(i,j)
      endif
    enddo
    enddo   
    print*,''
    print*,'    rho_c: ',minval(rhoe),maxval(rhoe),sum(rhoe*1d0)/ng/ng
    print*,'    compute and write into'
    deallocate(delta)

    print*,''
    print*,'    org: ',minval(rhoe),maxval(rhoe), sum(rhoe*1d0)/ng/ng
    print*,'    bad  : ',count_bad  
    count_bad = -1
    c_count = 0

    do  
      c_count = c_count + 1
      count_bad_prev = count_bad
      count_bad = 0

      do ilayer = 1, 4
      !!! 
      !$omp paralleldo default(shared) &
      !!! 
      !$omp& private(i,j,di,dj,ii,jj,n,sum_d) reduction(+:count_bad)
      ! do i = 1, ng
      do i = ilayer, ng, 4
      do j = 1, ng
        if (is_bad(rhoe(i,j))) then
          n = 0
          sum_d = 0.0
          do di = -1, 1
            ii = modulo(i+di-1, ng)+1
          do dj = -1, 1
            if (dj == 0 .and. di == 0)  cycle
            jj = modulo(j+dj-1, ng)+1
            if (.not. is_bad(rhoe(ii,jj))) then
                sum_d = sum_d + rhoe(ii,jj)
                n = n + 1
            endif
          enddo
          enddo

          if (n >= n_min) then
            rhoe(i,j) = sum_d / real(n)
          else
            count_bad = count_bad + 1
          endif
        endif
      enddo
      enddo
      !!! 
      !$omp endparalleldo
      enddo !end ilayer

      ! 停止条件
      if (count_bad == 0) then
        print*,'    new: ',minval(rhoe),maxval(rhoe),sum(rhoe*1d0)/ng/ng
        print*,'        ',c_count,' cycles'
        state = 1
        exit
      else if (count_bad == count_bad_prev) then
        print *, 'Error:full '//rho_name//dsp_name//'_x stagnated, still', count_bad, " zeros left in ",c_count,'cycles'
        state = 0
        stop 'Error in full_rhoe'
      end if
    end do


    print*,''
    print*,'    compute and write into'
    print*,'      save:',output_name(trim(adjustl(rho_name))//dsp_name//'_xf')
    open(16,file=output_name(trim(adjustl(rho_name))//dsp_name//'_xf'),status='replace',access='stream')
    write(16) rhoe
    close(16)
  endsubroutine

  subroutine decompose_Mesh_D(dsp_name)
    implicit none
    character(len=*), intent(in) :: dsp_name
    integer, parameter,dimension(4) :: DD = [-2,-1,1,2]
    real,allocatable :: dsp(:,:,:),dsp_t(:,:,:)
    real(8),allocatable :: A_mesh(:,:,:,:),trace_A(:,:),det_A(:,:)
    real(4),allocatable :: kappa(:,:),gamma1(:,:),gamma2(:,:),omega(:,:),lambda1(:,:),lambda2(:,:),mu(:,:)

    allocate(dsp(2,-1:ng+2,-1:ng+2),dsp_t(2,ng,ng), A_mesh(2,2,ng,ng))
    dsp = 0
    print*,'  read:'
    print*,'    ',output_name('dsp_'//dsp_name)
    open(15,file=output_name('dsp_'//dsp_name),status='old',access='stream')
    read(15) dsp_t
    close(15)

    dsp(1,1:ng,1:ng) = dsp_t(1,1:ng,1:ng)
    dsp(2,1:ng,1:ng) = dsp_t(2,1:ng,1:ng)

    dsp(:,   -1: 0  , :) = dsp(:, ng-1:ng, :)
    dsp(:, ng+1:ng+2, :) = dsp(:,  1:2   , :)
    dsp(:, :,   -1: 0  ) = dsp(:, :, ng-1:ng)
    dsp(:, :, ng+1:ng+2) = dsp(:, :,  1:2   )
    deallocate(dsp_t)

    A_mesh(1,1,:,:) = 1
    A_mesh(2,1,:,:) = 0
    do i = 1,4
      A_mesh(1,1,:,:) = A_mesh(1,1,:,:) + dsp(1,1+DD(i):ng+DD(i),1:ng)*weight(i)
      A_mesh(2,1,:,:) = A_mesh(2,1,:,:) + dsp(2,1+DD(i):ng+DD(i),1:ng)*weight(i)
    enddo

    A_mesh(1,2,:,:) = 0
    A_mesh(2,2,:,:) = 1
    do j = 1,4
      A_mesh(1,2,:,:) = A_mesh(1,2,:,:) + dsp(1,1:ng,1+DD(j):ng+DD(j))*weight(j)   
      A_mesh(2,2,:,:) = A_mesh(2,2,:,:) + dsp(2,1:ng,1+DD(j):ng+DD(j))*weight(j)
    enddo
    deallocate(dsp)

    allocate(trace_A(ng,ng),det_A(ng,ng))
    trace_A = A_mesh(1,1,:,:) + A_mesh(2,2,:,:)
    det_A   = A_mesh(1,1,:,:)*A_mesh(2,2,:,:) - A_mesh(1,2,:,:)*A_mesh(2,1,:,:)


    allocate(mu(ng,ng))
    mu = abs(1/det_A)
    print*,'    u     : ',minval(mu),maxval(mu),sum(mu*1d0)/ng/ng
    open(21,file=output_name('u_'//dsp_name//'_q'),status='replace',access='stream')
    write(21) mu
    close(21)
    deallocate(det_A,mu)


    allocate(kappa(ng, ng))
    kappa = 1.0D0 - 0.5D0 * trace_A
    print*,'    kappa  : ',minval(kappa),maxval(kappa),sum(kappa*1d0)/ng/ng
    open(16,file=output_name('k_'//dsp_name//'_q'),status='replace',access='stream')
    write(16) kappa
    close(16)
    
    allocate(gamma1(ng, ng),gamma2(ng, ng))
    gamma1 = 0.5D0 * (A_mesh(1,1,:,:) - A_mesh(2,2,:,:))
    gamma2 = 0.5D0 * (A_mesh(1,2,:,:) + A_mesh(2,1,:,:))
    print*,'    gamma1 : ',minval(gamma1),maxval(gamma1),sum(gamma1*1d0)/ng/ng
    print*,'    gamma2 : ',minval(gamma2),maxval(gamma2),sum(gamma2*1d0)/ng/ng
    open(17,file=output_name('g1_'//dsp_name//'_q'),status='replace',access='stream')
    write(17) gamma1
    close(17)
    open(18,file=output_name('g2_'//dsp_name//'_q'),status='replace',access='stream')
    write(18) gamma2
    close(18)


    allocate(lambda1(ng, ng),lambda2(ng, ng))
    trace_A = sqrt(gamma1**2 + gamma2**2)
    print*,'    gamma  : ',minval(trace_A),maxval(trace_A),sum(trace_A*1d0)/ng/ng
    open(18,file=output_name('G_'//dsp_name//'_q'),status='replace',access='stream')
    write(18) trace_A
    close(18)
    lambda1 = -kappa + trace_A
    lambda2 = -kappa - trace_A
    print*,'    lambda1: ',minval(lambda1),maxval(lambda1),sum(lambda1*1d0)/ng/ng
    print*,'    lambda2: ',minval(lambda2),maxval(lambda2),sum(lambda2*1d0)/ng/ng
    open(20,file=output_name('l1_'//dsp_name//'_q'),status='replace',access='stream')
    write(20) lambda1
    close(20)
    open(20,file=output_name('l2_'//dsp_name//'_q'),status='replace',access='stream')
    write(20) lambda2
    close(20)
    gamma1 = 0
    where (lambda1 < 0 .and. lambda2 < 0) gamma1 = -1
    where (lambda1 > 0 .and. lambda2 > 0) gamma1 = 1

    open(20,file=output_name('A_'//dsp_name//'_q'),status='replace',access='stream')
    write(20) gamma1
    close(20)
    deallocate(gamma1,gamma2)
    deallocate(kappa)
    deallocate(lambda1,lambda2,trace_A)

    allocate(omega(ng, ng))
    omega = 0.5D0 * (A_mesh(2,1,:,:) - A_mesh(1,2,:,:))
    print*,'    omega  : ',minval(omega),maxval(omega),sum(omega*1d0)/ng/ng
    open(19,file=output_name('j_'//dsp_name//'_q'),status='replace',access='stream')
    write(19) omega
    close(19)
    deallocate(omega,A_mesh)
  endsubroutine decompose_Mesh_D

  subroutine decompose_Mesh_FFT(dsp_name)
    implicit none
    character(len=*), intent(in) :: dsp_name
    real,allocatable ::  dsp(:,:,:)
    complex,allocatable :: cdiv(:,:)
    integer i_dim,j_dim,i,j
    real(8),allocatable :: A_mesh(:,:,:,:),trace_A(:,:),det_A(:,:)
    real(4),allocatable :: kappa(:,:),gamma1(:,:),gamma2(:,:),omega(:,:),lambda1(:,:),lambda2(:,:),mu(:,:)

    allocate(dsp(2,ng,ng))
    dsp = 0
    print*,'  read:'
    print*,'    ',output_name('dsp_'//dsp_name)
    open(15,file=output_name('dsp_'//dsp_name),status='old',access='stream')
    read(15) dsp
    close(15)
    ! allocate(cdiv(ng/2+1,ng))
    ! cdiv=0
    allocate(A_mesh(2,2,ng,ng))
    A_mesh = 0
    do i_dim=1,2
    do j_dim=1,2
      print*,'     working on dim',int(i_dim,1),int(j_dim,1)
      rho1(1:ng,1:ng)=dsp(i_dim,1:ng,1:ng)
      call sfftw_execute( plan) ! Fourier transform
      
      !$omp paralleldo default(shared) &
      !$omp& private(i,j,kx,pdim)
      do j=1,ng
      do i=1,ng/2+1
        if (j == 1 .and. i == 1) cycle
        kx=modulo([i,j]+ng/2-1,ng)-ng/2 !k
        pdim=sin(2*pi*kx/ng)
        rho1k(i,j) = (0,1)*rho1k(i,j)*pdim(j_dim) !c means complex 
        ! if (j_dim == i_dim) cdiv(i,j)=cdiv(i,j) + rho1k(i,j)
      enddo
      enddo
      !$omp endparalleldo
      rho1k(1,1)=0
      call sfftw_execute(iplan) ! Fourier transform
      A_mesh(i_dim,j_dim,1:ng,1:ng) = real(rho1(1:ng,1:ng))/real(ng**2)
      ! if (j_dim == i_dim) cdiv(1:ng,1:ng)=cdiv(1:ng,1:ng) + A_mesh(i_dim,j_dim,1:ng,1:ng)
      print*,'    A_mesh : ',i_dim,j_dim,minval(A_mesh(i_dim,j_dim,1:ng,1:ng)),maxval(A_mesh(i_dim,j_dim,1:ng,1:ng)),sum((A_mesh(i_dim,j_dim,1:ng,1:ng))*1d0)
    enddo
    enddo ! i_dim
    deallocate(dsp)

    A_mesh(1,1,:,:) = A_mesh(1,1,:,:) + 1
    A_mesh(2,2,:,:) = A_mesh(2,2,:,:) + 1


    allocate(trace_A(ng,ng),det_A(ng,ng))
    trace_A = A_mesh(1,1,:,:) + A_mesh(2,2,:,:)
    det_A   = A_mesh(1,1,:,:)*A_mesh(2,2,:,:) - A_mesh(1,2,:,:)*A_mesh(2,1,:,:)


    allocate(mu(ng,ng))
    mu = abs(1/det_A)
    print*,'    u     : ',minval(mu),maxval(mu),sum(mu*1d0)/ng/ng
    open(21,file=output_name('u_'//dsp_name//'_q'),status='replace',access='stream')
    write(21) mu
    close(21)
    deallocate(det_A,mu)


    allocate(kappa(ng, ng))
    kappa = 1.0D0 - 0.5D0 * trace_A
    print*,'    kappa  : ',minval(kappa),maxval(kappa),sum(kappa*1d0)/ng/ng
    open(16,file=output_name('k_'//dsp_name//'_q'),status='replace',access='stream')
    write(16) kappa
    close(16)
    
    allocate(gamma1(ng, ng),gamma2(ng, ng))
    gamma1 = 0.5D0 * (A_mesh(1,1,:,:) - A_mesh(2,2,:,:))
    gamma2 = 0.5D0 * (A_mesh(1,2,:,:) + A_mesh(2,1,:,:))
    print*,'    gamma1 : ',minval(gamma1),maxval(gamma1),sum(gamma1*1d0)/ng/ng
    print*,'    gamma2 : ',minval(gamma2),maxval(gamma2),sum(gamma2*1d0)/ng/ng
    open(17,file=output_name('g1_'//dsp_name//'_q'),status='replace',access='stream')
    write(17) gamma1
    close(17)
    open(18,file=output_name('g2_'//dsp_name//'_q'),status='replace',access='stream')
    write(18) gamma2
    close(18)

    allocate(omega(ng, ng))
    omega = 0.5D0 * (A_mesh(2,1,:,:) - A_mesh(1,2,:,:))
    print*,'    omega  : ',minval(omega),maxval(omega),sum(omega*1d0)/ng/ng
    open(19,file=output_name('j_'//dsp_name//'_q'),status='replace',access='stream')
    write(19) omega
    close(19)
    deallocate(omega,A_mesh)


    allocate(lambda1(ng, ng),lambda2(ng, ng))
    trace_A = sqrt(gamma1**2 + gamma2**2)
    lambda1 = -kappa + trace_A
    lambda2 = -kappa - trace_A
    print*,'    lambda1: ',minval(lambda1),maxval(lambda1),sum(lambda1*1d0)/ng/ng
    print*,'    lambda2: ',minval(lambda2),maxval(lambda2),sum(lambda2*1d0)/ng/ng
    open(20,file=output_name('l1_'//dsp_name//'_q'),status='replace',access='stream')
    write(20) lambda1
    close(20)
    open(20,file=output_name('l2_'//dsp_name//'_q'),status='replace',access='stream')
    write(20) lambda2
    close(20)
    gamma1 = 0
    where (lambda1 < 0 .and. lambda2 < 0) gamma1 = -1
    where (lambda1 > 0 .and. lambda2 > 0) gamma1 = 1

    open(20,file=output_name('A_'//dsp_name//'_q'),status='replace',access='stream')
    write(20) gamma1
    close(20)
    deallocate(gamma1,gamma2)
    deallocate(kappa)
    deallocate(lambda1,lambda2,trace_A)
  endsubroutine decompose_Mesh_FFT

  subroutine gassim(rho0)
    implicit none
    real, intent(inout) :: rho0(ng,ng)
    ! 注意：rho1, rho1k, plan 等必须在 module 级定义或作为参数传入
    rho1(1:ng,1:ng) = rho0
    call sfftw_execute(plan) 
    rho1k = rho1k * grid_sk
    call sfftw_execute(iplan)
    rho0 = rho1(1:ng,1:ng)/ng/ng
  endsubroutine


  subroutine dsp_halo()
    implicit none
    integer(4) i,j,idx1,idx2
    real,allocatable :: dsp(:,:,:)
    integer np_head
    real,allocatable :: xp(:,:),qp(:,:)

    print*,''
    print*,'    Computing halo displacement field'

    open(11,file=output_name('halo'),status='old',access='stream')
    read(11) np_head ! skip-blink
    read(11) np_head ! skip-np_min
    read(11) np_head
    close(11)
    allocate(xp(2,np_head),qp(2,np_head),dsp(2,ng,ng))
    open(12,file=output_name('halo_xp_mean_only'),status='old',access='stream')
    read(12) xp(1:2,1:np_head)
    close(12)
    open(12,file=output_name('halo_qp_mean_only'),status='old',access='stream')
    read(12) qp(1:2,1:np_head)
    close(12)

    print*,'    np_head: ',np_head
    dsp = 0
    do i=1,np_head
      idx1 = mod(floor(xp(1,i)),ng)+1
      idx2 = mod(floor(xp(2,i)),ng)+1
      dsp(:,idx1,idx2) = modulo(xp(:,i)-qp(:,i) + ng/2, real(ng)) - ng/2
    enddo

    do i_dim=1,2
      print*, '      dsp_H: dimension',int(i_dim,1),'min,max values ='
      print*, '      ',minval(dsp(i_dim,:,:)), maxval(dsp(i_dim,:,:))
    enddo

    print*,'      save:',output_name('dsp_H')
    open(15,file=output_name('dsp_H'),status='replace',access='stream')
    write(15) dsp
    close(15)

    print*,''
    print*,'    Smoothing halo displacement field'
    do i_dim=1,2
      call gassim(dsp(i_dim,1:ng,1:ng))
      print*, '      dsp_sH: dimension',int(i_dim,1),'min,max values ='
      print*, '      ',minval(dsp(i_dim,:,:)), maxval(dsp(i_dim,:,:))
    enddo

    print*,'      save:',output_name('dsp_sH')
    open(15,file=output_name('dsp_sH'),status='replace',access='stream')
    write(15) dsp
    close(15)
    deallocate(dsp,xp,qp)
  endsubroutine

  subroutine dsp_void()
    implicit none
    integer(4) i,j,idx1,idx2
    real dsp(2,ng,ng)
    integer np_head
    real,allocatable :: xp(:,:),qp(:,:)

    
    open(11,file=output_name('void'),status='old',access='stream')
    read(11) np_head ! skip-blink
    read(11) np_head ! skip-np_min
    read(11) np_head
    close(11)
    allocate(xp(2,np_head),qp(2,np_head))
    open(12,file=output_name('void_xp_mean_only'),status='old',access='stream')
    read(12) xp(1:2,1:np_head)
    close(12)
    open(12,file=output_name('void_qp_mean_only'),status='old',access='stream')
    read(12) qp(1:2,1:np_head)
    close(12)

    dsp = 0
    do i=1,np_head
      idx1 = mod(floor(xp(1,i)),ng)+1
      idx2 = mod(floor(xp(2,i)),ng)+1
      dsp(:,idx1,idx2) = modulo(xp(:,i)-qp(:,i) + ng/2, real(ng)) - ng/2
    enddo

    do i_dim=1,2
      print*, '      dsp: dimension',int(i_dim,1),'min,max values ='
      print*, '      ',minval(dsp(i_dim,:,:)), maxval(dsp(i_dim,:,:))
      call gassim(dsp(i_dim,1:ng,1:ng))
      print*, '      ',minval(dsp(i_dim,:,:)), maxval(dsp(i_dim,:,:))
    enddo

    print*,'    Write Smoothed dsp into file:'
    print*,'      save:',output_name('dsp_sV')
    open(15,file=output_name('dsp_sV'),status='replace',access='stream')
    write(15) dsp
    close(15)
  endsubroutine

  subroutine gassim_file(rho_name)
    implicit none
    character(len=*), intent(in) :: rho_name

    open(15,file=output_name(rho_name),status='old',access='stream')
    read(15) rho1(1:ng,1:ng)
    close(15)
    call sfftw_execute(plan) 
    rho1k = rho1k * grid_sk
    call sfftw_execute(iplan)
    rho1(1:ng,1:ng) = rho1(1:ng,1:ng)/ng/ng

    open(15,file=output_name('s'//rho_name),status='replace',access='stream')
    write(15) rho1(1:ng,1:ng)
    close(15)
  endsubroutine

  logical function is_bad(val)
    use, intrinsic :: ieee_arithmetic
    real, intent(in) :: val
    is_bad = (abs(val)>1e10  .or. (val == 0.0) .or. (ieee_is_nan(val)) .or. (.not. ieee_is_finite(val)))
  endfunction
end