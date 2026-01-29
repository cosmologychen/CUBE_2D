! #define merge_projection
! #define density_nu
! #define density_matter
! #define power_matter

program density_power
  use omp_lib
  use parameters
  use powerspectrum
  implicit none
  save
  include 'fftw3.f'

  logical,parameter :: DTFE = .false.

  integer(8) plan,iplan

  real(8) rho8
  real(4),allocatable :: xp(:,:),rho_grid(:,:,:)
  
  integer i,j,l,iteam,idx(2),idx1(2),idx2(2),cur_checkpoint
  integer(8) np
  real pos1(2),dx1(2),dx2(2),xi(10,0:nbin)

  allocate(rho_grid(0:nw+1,0:nw+1,ncore))

  print*, 'cicpower on resolution: nw=',nw
  print*, 'checkpoint at:'
  open(16,file='./z_checkpoint.txt',status='old')
  do i=1,nmax_redshift
    read(16,end=71,fmt='(f8.4)') z_checkpoint(i)!; print*, z_checkpoint(i)
  enddo
  71 n_checkpoint=i-1
  close(16); print*,''

  call omp_set_num_threads(ncore)

  call sfftw_init_threads(l)
  print*, '    sfftw_init_threads status',l
  call sfftw_plan_with_nthreads(ncore)
  call sfftw_plan_dft_r2c_2d( plan,nw,nw,rho1k,rho1k,FFTW_MEASURE)
  call sfftw_plan_dft_c2r_2d(iplan,nw,nw,rho1k,rho1k,FFTW_MEASURE)


  ! do cur_checkpoint= 49,49
  ! do cur_checkpoint= n_checkpoint,n_checkpoint
  do cur_checkpoint= 1,n_checkpoint
    sim%cur_checkpoint=cur_checkpoint
    print*, ''
    print*, '==========================================='
    print*, '==========================================='
    print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
    !print*,output_name('info')
    open(11,file=output_name('info'),access='stream'); read(11) sim; close(11)
    np=sim%np
    print*, 'np =',np
    allocate(xp(2,np))
    open(11,file=output_name('xp'),access='stream'); read(11) xp; close(11)
    rho_grid=0
    
    !$omp paralleldo default(shared) schedule(dynamic)&
    !$omp& private(l,pos1,iteam,idx1,idx2,dx1,dx2)
    do l=1,np
      iteam=omp_get_thread_num()+1
      pos1=xp(:,l)*nic-0.5
      idx1=floor(pos1)+1; idx2=idx1+1
      dx1=idx1-pos1;      dx2=1-dx1
      rho_grid(idx1(1),idx1(2),iteam)=rho_grid(idx1(1),idx1(2),iteam)+dx1(1)*dx1(2)
      rho_grid(idx1(1),idx2(2),iteam)=rho_grid(idx1(1),idx2(2),iteam)+dx1(1)*dx2(2)
      rho_grid(idx2(1),idx1(2),iteam)=rho_grid(idx2(1),idx1(2),iteam)+dx2(1)*dx1(2)
      rho_grid(idx2(1),idx2(2),iteam)=rho_grid(idx2(1),idx2(2),iteam)+dx2(1)*dx2(2)
    enddo ! ip
    !$omp endparalleldo
              

    ! do iteam=1,ncore
    ! !$omp paralleldo default(shared) schedule(dynamic)&
    ! !$omp& private(i,j)
    ! do j=1,nw
    ! do i=1,nw
    !   rho1(i,j)=rho_grid(i,j,iteam)
    ! enddo
    ! enddo  
    ! !$omp endparalleldo
    ! enddo
    rho1 = 0
    do iteam=1,ncore
      rho1(:nw,:nw) = rho1(:nw,:nw) + rho_grid(1:nw,1:nw,iteam)
    enddo

    print*, 'check: min,max,sum of rho_grid = '
    ! 
    print*,minval(rho1),maxval(rho1),sum(rho1*1d0)
    rho8=sum(rho1*1d0)/real(nw*nw)
    print*,rho8
    do i=1,nw
      rho1(:,i)=rho1(:,i)/(rho8)-1
    enddo


    print*,'min',minval(rho1),'max',maxval(rho1),'mean',sum(rho1*1d0)/nw/nw;

    print*,'Write delta_c into',output_name('delta_c')
    open(11,file=output_name('delta_c'),status='replace',access='stream')
    write(11) rho1(:nw,:nw)
    close(11)

    ! open(101,file=output_name('delta_ic'),access='stream')
    ! do i=1,ngic
    !   read(101) rho1(1:ngic,i) ! write layer by layer to avoid bug
    ! enddo

    print*,'auto_power'
    call sfftw_execute( plan)
    rho1k=rho1k/real(nw*nw)
    call auto_power(xi,np,2)

    open(11,file=output_name('power'),status='replace',access='stream')
    write(11) xi
    close(11)
    
    deallocate(xp)
  enddo
  deallocate(rho_grid)
  call sfftw_destroy_plan( plan)
  call sfftw_destroy_plan(iplan)
  call sfftw_cleanup_threads()
  print*,'cicpower done'
endprogram