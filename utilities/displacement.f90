#define Cpower
program displacement
  use omp_lib
  use parameters
  use powerspectrum
  implicit none
  save
  include 'fftw3.f'




  integer :: ip,iq(2),i_dim,idx1(2),idx2(2)
  real :: dx1(2),dx2(2)
  integer(8) np,istat,nthreads,plan,iplan

  integer :: i,j,iteam,cur_checkpoint
  real :: pos0(2),pos1(2),dpos(2),kx(2),pdim(2),xi(10,0:nbin)
  real,allocatable :: rho_grid(:,:,:)
  real,allocatable :: dsp(:,:,:),xp(:,:)
  complex,allocatable :: cdiv(:,:),cphi(:,:),rhok_L(:,:),rhok_R(:,:),rhok_N(:,:)

  ! print*,modulo(512-1,512)+1,modulo(128-1,512)+1,modulo(0-1,512)+1
  ! stop

  print*, 'Displacement field analysis on resolution:'
  print*, 'ng=',ng
  print*, 'checkpoint at:'
  open(16,file='../z_checkpoint.txt',status='old')
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


  do cur_checkpoint= 1,6!n_checkpoint,n_checkpoint !2,-1
    print*, ''
    print*,'==========================================================='



    print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))

    sim%cur_checkpoint=cur_checkpoint
    open(11,file=output_name('info'),access='stream'); read(11) sim; close(11)
    np = sim%np
#ifdef Cpower
    sim%cur_checkpoint = 1
    open(11,file=output_name('delta_L'),status='old',access='stream')
    read(11) rho1(1:ng,1:ng)
    call sfftw_execute( plan)
    allocate(rhok_L(ng/2+1,ng))
    rhok_L=rho1k/real(ng*ng)
    sim%cur_checkpoint=cur_checkpoint

    open(11,file=output_name('delta_c'),status='old',access='stream')
    read(11) rho1(1:ng,1:ng)
    call sfftw_execute( plan)
    allocate(rhok_N(ng/2+1,ng))
    rhok_N=rho1k/real(ng*ng)
    call cross_power(xi,rhok_L,rhok_N,np,2)
    print*,'   save: ',output_name('Cpower_LN')
    open(11,file=output_name('Cpower_LN'),status='replace',access='stream')
    write(11) xi
    close(11)
    deallocate(rhok_N,rhok_L)
    ! cycle
#endif
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
    print*,'    Start computing delta_E'
    allocate(cdiv(ng/2+1,ng),cphi(ng/2+1,ng))
    cphi=0
    cdiv=0
    do i_dim=1,2
      ! print*,'     working on dim',int(i_dim,1)
      rho1(1:ng,1:ng)=dsp(i_dim,1:ng,1:ng)
      call sfftw_execute( plan) ! Fourier transform

      !$omp paralleldo default(shared) &
      !$omp& private(i,j,kx,pdim)
      do j=1,ng
      do i=1,ng/2+1
        if (j == 1 .and. i == 1) cycle
        kx=modulo([i,j]+ng/2-1,ng)-ng/2 !k
        pdim=sin(2*pi*kx/ng)
        cphi(i,j)=cphi(i,j)+(0,1)*rho1k(i,j)*pdim(i_dim)/(-sum(pdim**2)) !phik 
        cdiv(i,j)=cdiv(i,j)+(0,1)*rho1k(i,j)*pdim(i_dim) !c means complex 
      enddo
      enddo
      !$omp endparalleldo
    enddo ! i_dim

    cphi(1,1)=0
    cdiv(1,1)=0

    rho1=0
    rho1k=cdiv

    call sfftw_execute(iplan)
    rho1 = -rho1/(ng*ng)
    print*,''
    print*,'    write delta_E'
    print*,'    ',minval(rho1(1:ng,1:ng)),maxval(rho1(1:ng,1:ng)),sum(rho1(1:ng,1:ng)*1d0)
    open(15,file=output_name('E_q'),status='replace',access='stream')
    write(15) rho1(1:ng,1:ng)
    close(15)

#ifdef Cpower
    call sfftw_execute( plan)
    allocate(rhok_R(ng/2+1,ng))
    rhok_R=rho1k/real(ng*ng)
    sim%cur_checkpoint = 1
    open(11,file=output_name('delta_L'),status='old',access='stream')
    read(11) rho1(1:ng,1:ng)
    call sfftw_execute( plan)
    allocate(rhok_L(ng/2+1,ng))
    rhok_L=rho1k/real(ng*ng)
    call cross_power(xi,rhok_L,rhok_R,np,2)
    sim%cur_checkpoint = cur_checkpoint
    open(11,file=output_name('Cpower_LR'),status='replace',access='stream')
    write(11) xi
    close(11)
    deallocate(rhok_R,rhok_L)
#endif

    open(11,file=output_name('phik_E'),status='replace',access='stream')
    write(11) cphi
    close(11)
    rho1k=cphi
    deallocate(cphi)


    call sfftw_execute(iplan)
    rho1 = rho1/(ng*ng)
    open(15,file=output_name('phi_E'),status='replace',access='stream')
    write(15) rho1(1:ng,1:ng)
    close(15)

    do i_dim=1,2
      !$omp paralleldo default(shared) &
      !$omp& private(i,j,kx)
      do j=1,ng
      do i=1,ng/2+1
        if (j == 1 .and. i == 1) cycle
        kx=modulo([i,j]+ng/2-1,ng)-ng/2 !k
        rho1k(i,j)=(0,1)*cdiv(i,j)*sin(2*pi*kx(i_dim)/ng)/(-sum(sin(2*pi*kx/ng)**2))
      enddo
      enddo
      !$omp endparalleldo
      rho1k(1,1) = 0
      call sfftw_execute(iplan)
      rho1 = rho1/(ng*ng)
      dsp(i_dim,1:ng,1:ng) = rho1(1:ng,1:ng)
    enddo
    deallocate(cdiv)


    print*,''
    do i_dim=1,2
      print*, '   dsp_E: dimension',int(i_dim,1),'min,max values ='
      print*, '   ',minval(dsp(i_dim,:,:)), maxval(dsp(i_dim,:,:))
    enddo

    print*,'    Write dsp_E into file:'
    print*,'      save:',output_name('dsp_E')
    open(15,file=output_name('dsp_E'),status='replace',access='stream')
    write(15) dsp
    close(15)
    deallocate(dsp)
    

    call decompose_Mesh_FFT('D')
    ! call decompose_Mesh_D('D')


    ! ! call dep2delta_e('E','E',4,istat)
    ! call dep2delta_e('D','E',4,istat)
    ! call dep2delta_e('D','uD',4,istat)
    ! call dep2delta_e('D','jD',4,istat)
    ! call dep2delta_e('D','kD',4,istat)
#ifdef Cpower
    sim%cur_checkpoint = 1
    open(11,file=output_name('delta_L'),status='old',access='stream')
    read(11) rho1(1:ng,1:ng)
    call sfftw_execute( plan)
    allocate(rhok_L(ng/2+1,ng))
    rhok_L=rho1k/real(ng*ng)
    sim%cur_checkpoint=cur_checkpoint

    open(11,file=output_name('kD_q'),status='old',access='stream')
    read(11) rho1(1:ng,1:ng)
    call sfftw_execute( plan)
    allocate(rhok_N(ng/2+1,ng))
    rhok_N=rho1k/real(ng*ng)
    call cross_power(xi,rhok_L,rhok_N,np,2)
    print*,'   save: ',output_name('Cpower_Lkq')
    open(11,file=output_name('Cpower_Lkq'),status='replace',access='stream')
    write(11) xi
    close(11)
    deallocate(rhok_N,rhok_L)
    ! cycle
#endif
#ifdef Cpower
    sim%cur_checkpoint = 1
    open(11,file=output_name('delta_L'),status='old',access='stream')
    read(11) rho1(1:ng,1:ng)
    call sfftw_execute( plan)
    allocate(rhok_L(ng/2+1,ng))
    rhok_L=rho1k/real(ng*ng)
    sim%cur_checkpoint=cur_checkpoint

    open(11,file=output_name('uD_q'),status='old',access='stream')
    read(11) rho1(1:ng,1:ng)
    call sfftw_execute( plan)
    allocate(rhok_N(ng/2+1,ng))
    rhok_N=rho1k/real(ng*ng)
    call cross_power(xi,rhok_L,rhok_N,np,2)
    print*,'   save: ',output_name('Cpower_Luq')
    open(11,file=output_name('Cpower_Luq'),status='replace',access='stream')
    write(11) xi
    close(11)
    deallocate(rhok_N,rhok_L)
    ! cycle
#endif
! #ifdef Cpower
!   sim%cur_checkpoint=cur_checkpoint
!   open(11,file=output_name('delta_c'),status='old',access='stream')
!   read(11) rho1(1:ng,1:ng)
!   call sfftw_execute( plan)
!   allocate(rhok_L(ng/2+1,ng))
!   rhok_L=rho1k/real(ng*ng)

!   open(11,file=output_name('uD_q'),status='old',access='stream')
!   read(11) rho1(1:ng,1:ng)
!   call sfftw_execute( plan)
!   allocate(rhok_N(ng/2+1,ng))
!   rhok_N=rho1k/real(ng*ng)
!   call cross_power(xi,rhok_L,rhok_N,np,2)
!   print*,'   save: ',output_name('Cpower_Ruq')
!   open(11,file=output_name('Cpower_Ruq'),status='replace',access='stream')
!   write(11) xi
!   close(11)
!   deallocate(rhok_N,rhok_L)
!   ! cycle
! #endif
! #ifdef Cpower
!   sim%cur_checkpoint = 1
!   open(11,file=output_name('delta_L'),status='old',access='stream')
!   read(11) rho1(1:ng,1:ng)
!   call sfftw_execute( plan)
!   allocate(rhok_L(ng/2+1,ng))
!   rhok_L=rho1k/real(ng*ng)
!   sim%cur_checkpoint=cur_checkpoint

!   open(11,file=output_name('uDD_xf'),status='old',access='stream')
!   read(11) rho1(1:ng,1:ng)
!   call sfftw_execute( plan)
!   allocate(rhok_N(ng/2+1,ng))
!   rhok_N=rho1k/real(ng*ng)
!   call cross_power(xi,rhok_L,rhok_N,np,2)
!   print*,'   save: ',output_name('Cpower_Lux')
!   open(11,file=output_name('Cpower_Lux'),status='replace',access='stream')
!   write(11) xi
!   close(11)
!   deallocate(rhok_N,rhok_L)
!   ! cycle
! #endif
! #ifdef Cpower
!   sim%cur_checkpoint=cur_checkpoint
!   open(11,file=output_name('delta_c'),status='old',access='stream')
!   read(11) rho1(1:ng,1:ng)
!   call sfftw_execute( plan)
!   allocate(rhok_L(ng/2+1,ng))
!   rhok_L=rho1k/real(ng*ng)

!   open(11,file=output_name('uDD_xf'),status='old',access='stream')
!   read(11) rho1(1:ng,1:ng)
!   call sfftw_execute( plan)
!   allocate(rhok_N(ng/2+1,ng))
!   rhok_N=rho1k/real(ng*ng)
!   call cross_power(xi,rhok_L,rhok_N,np,2)
!   print*,'   save: ',output_name('Cpower_Rux')
!   open(11,file=output_name('Cpower_Rux'),status='replace',access='stream')
!   write(11) xi
!   close(11)
!   deallocate(rhok_N,rhok_L)
!   ! cycle
! #endif


! call decompose_Mesh_FFT('E')
! ! call decompose_Mesh_D('E')

! call dep2delta_e('E','E',4,istat)
! call dep2delta_e('E','uE',4,istat)
! call dep2delta_e('E','jE',4,istat)
! call dep2delta_e('E','kE',4,istat)
  enddo
  print*,'displacement done'

  contains

  ! subroutine full_delta1(namespace,rhoe,n_min,state)
  !   implicit none
  !   integer, intent(in)  :: n_min     ! 至少需要 n 个非零邻居
  !   real,intent(inout) :: rhoe(ng,ng)
  !   character(len=*), intent(in) :: namespace
  !   integer(8) , intent(inout) :: state

  !   integer :: i, j, di, dj, ii, jj, n, count_0, count_0_prev,c_count
  !   real :: sum_d
  !   real delta(ng, ng)

  !   if (n_min > 8) then
  !       print *, "n_min too large, exiting subroutine."
  !       return
  !   end if

  !   delta = 0
  !   open(15,file=output_name(namespace),status='old',access='stream')
  !   read(15) delta
  !   close(15)

  !   delta = delta+1
  !   print*,'    org: ',minval(delta),maxval(delta), sum(delta*1d0)
  !   count_0 = -1
  !   c_count = 0

  !   do  
  !       c_count = c_count + 1
  !       count_0_prev = count_0
  !       count_0 = 0

  !       !!! 
  !       !$omp paralleldo default(shared) &
  !       !!! 
  !       !$omp& private(i,j,di,dj,ii,jj,n,sum_d) reduction(+:count_0)
  !       do i = 1, ng
  !       do j = 1, ng
  !           if (delta(i,j) == 0.0) then
  !               n = 0
  !               sum_d = 0.0
  !               do di = -1, 1
  !                   ii = modulo(i+di-1, ng)+1
  !               do dj = -1, 1
  !                   if (dj == 0 .and. di == 0)  cycle
  !                   jj = modulo(j+dj-1, ng)+1
  !                   if (delta(ii,jj) /= 0.0) then
  !                       sum_d = sum_d + delta(ii,jj)
  !                       n = n + 1
  !                   endif
  !               enddo
  !               enddo

  !               if (n >= n_min) then
  !                   delta(i,j) = sum_d / real(n)
  !               else
  !                   count_0 = count_0 + 1
  !                   ! if (c_count == 24) then
  !                   !     print*,''
  !                   !     print*,'+++++++++++++++++++++++'
  !                   !     print*,'i,j',i,j

  !                   !     do di = -1, 1
  !                   !         ii = modulo(i+di-1, ng)+1
  !                   !         print*,'    ii',ii,i+di
  !                   !     do dj = -1, 1
  !                   !         if (dj == 0 .and. di == 0)  cycle
  !                   !         jj = modulo(j+dj-1, ng)+1
  !                   !         print*,'        jj',jj,j+dj
  !                   !         print*,'            d',delta(ii,jj)
  !                   !         if (delta(ii,jj) /= 0.0) then
  !                   !             sum_d = sum_d + delta(ii,jj)
  !                   !             n = n + 1
  !                   !         end if
  !                   !     enddo
  !                   !     enddo
  !                   !     if (count_0 > 10) stop
  !                   ! endif
  !               endif
  !           endif
  !       enddo
  !       enddo
  !       !!! 
  !       !$omp endparalleldo

  !       ! 停止条件
  !       if (count_0 == 0) then
  !           state = 1
  !           print*,'    new: ',minval(delta),maxval(delta), sum(delta*1d0)
  !           print*,c_count,'cycles'
  !           print*,'    write full_rho '//namespace//' into'
  !           print*,'        ',output_name(namespace//'f')
  !           open(16,file=output_name(namespace//'f'),status='replace',access='stream')
  !           write(16) delta-1
  !           close(16)
  !           exit
  !       else if (count_0 == count_0_prev) then
  !           print *, 'Error:full_rho '//namespace//' stagnated, still', count_0, " zeros left in ",c_count,'cycles'
  !           state = 0
  !           stop 'Error in full_delta'
  !       end if
  !   enddo

  !   rhoe = rhoe/delta
  !   print*,''
  !   print*,'    compute and write into'
  !   print*,'        ',output_name(namespace//'e')
  !   print*,'    ',minval(rhoe),maxval(rhoe),sum(rhoe*1d0)
  !   open(16,file=output_name(namespace//'e'),status='replace',access='stream')
  !   write(16) rhoe
  !   close(16)
  ! endsubroutine

  ! subroutine full_rho_Conservation(namespace,rho,n_min,state)
  !   implicit none
  !   integer, intent(in)  :: n_min     ! 至少需要 n 个非零邻居
  !   character(len=*), intent(in) :: namespace
  !   real, intent(out) :: rho(ng, ng)
  !   integer(8) , intent(inout) :: state

  !   integer :: ilayer,i, j, di, dj, ii, jj, n, count_0, count_0_prev,c_count
  !   real :: sum_d

  !   if (n_min > 8 .or. n_min < 2) then
  !       print *, "bad n_min = ",n_min, "exiting subroutine."
  !       return
  !   end if

  !   print*,'    org: ',minval(rho),maxval(rho), sum(rho*1d0)
  !   count_0 = -1
  !   c_count = 0

  !   do  
  !       c_count = c_count + 1
  !       count_0_prev = count_0
  !       count_0 = 0

  !       do ilayer = 0, 3
  !       !!! 
  !       !$omp paralleldo default(shared) &
  !       !!! 
  !       !$omp& private(i,j,di,dj,ii,jj,n,sum_d) reduction(+:count_0)
  !       do i = 1+ilayer, ng, 4
  !       do j = 1, ng
  !           if (rho(i,j) == 0.0) then
  !               n = 0
  !               sum_d = 0.0
  !               do di = -1, 1
  !                   ii = modulo(i+di-1, ng)+1
  !               do dj = -1, 1
  !                   if (dj == 0 .and. di == 0)  cycle
  !                   jj = modulo(j+dj-1, ng)+1
  !                   if (rho(ii,jj) /= 0.0) then
  !                       sum_d = sum_d + rho(ii,jj)
  !                       n = n + 1
  !                   endif
  !               enddo
  !               enddo

  !               if (n >= n_min) then

  !                   rho(i,j) = sum_d / real(n)
  !                   do di = -1, 1
  !                       ii = modulo(i+di-1, ng)+1
  !                   do dj = -1, 1
  !                       if (dj == 0 .and. di == 0)  cycle
  !                       jj = modulo(j+dj-1, ng)+1
  !                       if (rho(ii,jj) /= 0.0) then
  !                           rho(ii,jj) = rho(ii,jj)*(real(n-1)/ real(n))
  !                       endif
  !                   enddo
  !                   enddo
  !               else
  !                   count_0 = count_0 + 1
  !                   ! if (c_count == 24) then
  !                   !     print*,''
  !                   !     print*,'+++++++++++++++++++++++'
  !                   !     print*,'i,j',i,j

  !                   !     do di = -1, 1
  !                   !         ii = modulo(i+di-1, ng)+1
  !                   !         print*,'    ii',ii,i+di
  !                   !     do dj = -1, 1
  !                   !         if (dj == 0 .and. di == 0)  cycle
  !                   !         jj = modulo(j+dj-1, ng)+1
  !                   !         print*,'        jj',jj,j+dj
  !                   !         print*,'            d',rho(ii,jj)
  !                   !         if (rho(ii,jj) /= 0.0) then
  !                   !             sum_d = sum_d + rho(ii,jj)
  !                   !             n = n + 1
  !                   !         end if
  !                   !     enddo
  !                   !     enddo
  !                   !     if (count_0 > 10) stop
  !                   ! endif
  !               endif
  !           endif
  !       enddo
  !       enddo
  !       !!! 
  !       !$omp endparalleldo
  !       enddo

  !       ! 停止条件
  !       if (count_0 == 0) then
  !           state = 1
  !           print*,'    new: ',minval(rho),maxval(rho), sum(rho*1d0)
  !           print*,c_count,'cycles'
  !           print*,'    write full_rho '//namespace//' into'
  !           print*,'        ',output_name(namespace//'fC')
  !           open(16,file=output_name(namespace//'fC'),status='replace',access='stream')
  !           write(16) rho-1
  !           close(16)
  !           exit
  !       else if (count_0 == count_0_prev) then
  !           print *, 'Error:full_rho '//namespace//' stagnated, still', count_0, " zeros left in ",c_count,'cycles'
  !           state = 0
  !           stop 'Error in full_rho'
  !       end if
  !   end do
  ! endsubroutine

  ! subroutine full_delta_Conservation(namespace,rhoe,n_min,state)
  !   implicit none
  !   integer, intent(in)  :: n_min     ! 至少需要 n 个非零邻居
  !   character(len=*), intent(in) :: namespace
  !   real, intent(inout) :: rhoe(ng, ng)
  !   integer(8) , intent(inout) :: state

  !   integer :: i, j, di, dj, ii, jj, n, count_0, count_0_prev,c_count
  !   real :: sum_d
  !   real :: delta(ng, ng)

  !   if (n_min > 8 .or. n_min < 2) then
  !       print *, "bad n_min = ",n_min, "exiting subroutine."
  !       return
  !   end if

  !   delta = 0.0
  !   open(15,file=output_name(namespace),status='old',access='stream')
  !   read(15) delta
  !   close(15)

  !   delta = delta+1
  !   call full_rho_Conservation(namespace,delta,n_min,state)
  !   call full_rho_Conservation('delta_E',rhoe,n_min,state)
  !   rhoe = rhoe/delta

  !   print*,''
  !   print*,'    compute and write into'
  !   print*,'        ',output_name(namespace//'eC')
  !   print*,'    ',minval(rhoe),maxval(rhoe),sum(rhoe*1d0)
  !   open(16,file=output_name(namespace//'eC'),status='replace',access='stream')
  !   write(16) rhoe
  !   close(16)
  ! endsubroutine

  subroutine dep2delta_e(name1,name2,n_min,state)
    implicit none
    integer, intent(in)  :: n_min     ! 至少需要 n 个非零邻居
    character(len=*), intent(in) :: name1,name2
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
    print*, 'CIC interpolation'//name2//' by dsp_'//name1

    allocate(dsp(2,ng,ng))
    print*,size(rho1k(:,1)),size(rho1k(1,:))
    dsp = 0
    print*,'  read:'
    print*,'    ',output_name('dsp_'//name1)
    print*,'    ',output_name(name2//'_q')
    open(15,file=output_name('dsp_'//name1),status='old',access='stream')
    read(15) dsp
    close(15)

    rho1=0
    open(15,file=output_name(name2//'_q'),status='old',access='stream')
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
    print*,'      save:',output_name('delta_c'//name1)
    open(16,file=output_name('delta_c'//name1),status='replace',access='stream')
    write(16) delta-1
    close(16)

    print*,''
    print*,'    compute and write into'
    print*,'      save:',output_name(name2//name1//'_x')
    open(16,file=output_name(name2//name1//'_x'),status='replace',access='stream')
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
    print*,'    delta: ',minval(delta),maxval(delta), sum(delta*1d0)
    print*,'    rho_c: ',minval(rhoe),maxval(rhoe),sum(rhoe*1d0)
    print*,'    bad  : ',count_bad  
    deallocate(delta)

    print*,''
    print*,'    org: ',minval(rhoe),maxval(rhoe), sum(rhoe*1d0)
    count_bad = -1
    c_count = 0

    do  
      c_count = c_count + 1
      count_bad_prev = count_bad
      count_bad = 0

      ! do ilayer = 0, 3
      !!! 
      !$omp paralleldo default(shared) &
      !!! 
      !$omp& private(i,j,di,dj,ii,jj,n,sum_d) reduction(+:count_bad)
      do i = 1, ng
      ! do i = 1+ilayer, ng, 4
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
            ! if (c_count == 24) then
            !     print*,''
            !     print*,'+++++++++++++++++++++++'
            !     print*,'i,j',i,j

            !     do di = -1, 1
            !         ii = modulo(i+di-1, ng)+1
            !         print*,'    ii',ii,i+di
            !     do dj = -1, 1
            !         if (dj == 0 .and. di == 0)  cycle
            !         jj = modulo(j+dj-1, ng)+1
            !         print*,'        jj',jj,j+dj
            !         print*,'            d',rhoe(ii,jj)
            !         if (rhoe(ii,jj) /= 0.0) then
            !             sum_d = sum_d + rhoe(ii,jj)
            !             n = n + 1
            !         end if
            !     enddo
            !     enddo
            !     if (count_bad > 10) stop
            ! endif
          endif
        endif
      enddo
      enddo
      !!! 
      !$omp endparalleldo
      ! enddo !end ilayer

      ! 停止条件
      if (count_bad == 0) then
        print*,'    new: ',minval(rhoe),maxval(rhoe),sum(rhoe*1d0)
        print*,'        ',c_count,' cycles'
        state = 1
        exit
      else if (count_bad == count_bad_prev) then
        print *, 'Error:full '//name2//name1//'_x stagnated, still', count_bad, " zeros left in ",c_count,'cycles'
        state = 0
        stop 'Error in full_rhoe'
      end if
    end do


    print*,''
    print*,'    compute and write into'
    print*,'      save:',output_name(name2//name1//'_xf')
    open(16,file=output_name(name2//name1//'_xf'),status='replace',access='stream')
    write(16) rhoe
    close(16)
  endsubroutine

  subroutine decompose_Mesh_D(namespace)
    implicit none
    character(len=*), intent(in) :: namespace
    integer, parameter,dimension(4) :: DD = [-2,-1,1,2]
    real,allocatable :: dsp(:,:,:),dsp_t(:,:,:)
    real(8),allocatable :: A_mesh(:,:,:,:),trace_A(:,:),det_A(:,:)
    real(4),allocatable :: kappa(:,:),gamma1(:,:),gamma2(:,:),omega(:,:),lambda1(:,:),lambda2(:,:),mu(:,:)

    allocate(dsp(2,-1:ng+2,-1:ng+2),dsp_t(2,ng,ng), A_mesh(2,2,ng,ng))
    dsp = 0
    print*,'  read:'
    print*,'    ',output_name('dsp_'//namespace)
    open(15,file=output_name('dsp_'//namespace),status='old',access='stream')
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

    allocate(kappa(ng, ng))
    kappa = 1.0D0 - 0.5D0 * trace_A
    print*,'    kappa  : ',minval(kappa),maxval(kappa),sum(kappa*1d0)
    open(16,file=output_name('k'//namespace//'_q'),status='replace',access='stream')
    write(16) kappa
    close(16)
    deallocate(kappa)
    
    allocate(lambda1(ng, ng),lambda2(ng, ng))
    lambda1 = 0.5D0 * (trace_A + SQRT(trace_A**2 - 4.0D0 * det_A))
    lambda2 = 0.5D0 * (trace_A - SQRT(trace_A**2 - 4.0D0 * det_A))
    print*,'    lambda1: ',minval(lambda1),maxval(lambda1),sum(lambda1*1d0)
    print*,'    lambda2: ',minval(lambda2),maxval(lambda2),sum(lambda2*1d0)
    open(20,file=output_name('l1'//namespace//'_q'),status='replace',access='stream')
    write(20) lambda1
    close(20)
    open(20,file=output_name('l2'//namespace//'_q'),status='replace',access='stream')
    write(20) lambda2
    close(20)
    deallocate(lambda1,lambda2,trace_A)

    allocate(mu(ng,ng))
    mu = abs(1/det_A)
    print*,'    u     : ',minval(mu),maxval(mu),sum(mu*1d0)
    open(21,file=output_name('u'//namespace//'_q'),status='replace',access='stream')
    write(21) mu
    close(21)
    deallocate(det_A,mu)


    allocate(gamma1(ng, ng),gamma2(ng, ng))
    gamma1 = 0.5D0 * (A_mesh(1,1,:,:) - A_mesh(2,2,:,:))
    gamma2 = 0.5D0 * (A_mesh(1,2,:,:) + A_mesh(2,1,:,:))
    print*,'    gamma1 : ',minval(gamma1),maxval(gamma1),sum(gamma1*1d0)
    print*,'    gamma2 : ',minval(gamma2),maxval(gamma2),sum(gamma2*1d0)
    open(17,file=output_name('g1'//namespace//'_q'),status='replace',access='stream')
    write(17) gamma1
    close(18)
    open(17,file=output_name('g2'//namespace//'_q'),status='replace',access='stream')
    write(18) gamma2
    close(18)
    deallocate(gamma1,gamma2)

    allocate(omega(ng, ng))
    omega = 0.5D0 * (A_mesh(2,1,:,:) - A_mesh(1,2,:,:))
    print*,'    omega  : ',minval(omega),maxval(omega),sum(omega*1d0)
    open(19,file=output_name('j'//namespace//'_q'),status='replace',access='stream')
    write(19) omega
    close(19)
    deallocate(omega,A_mesh)
  endsubroutine decompose_Mesh_D

  subroutine decompose_Mesh_FFT(namespace)
    implicit none
    character(len=*), intent(in) :: namespace
    real,allocatable :: dsp(:,:,:)
    integer i_dim,j_dim,i,j
    real(8),allocatable :: A_mesh(:,:,:,:),trace_A(:,:),det_A(:,:)
    real(4),allocatable :: kappa(:,:),gamma1(:,:),gamma2(:,:),omega(:,:),lambda1(:,:),lambda2(:,:),mu(:,:)

    allocate(dsp(2,ng,ng))
    dsp = 0
    print*,'  read:'
    print*,'    ',output_name('dsp_'//namespace)
    open(15,file=output_name('dsp_'//namespace),status='old',access='stream')
    read(15) dsp
    close(15)
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
      enddo
      enddo
      !$omp endparalleldo
      call sfftw_execute(iplan) ! Fourier transform
      A_mesh(i_dim,j_dim,1:ng,1:ng) = real(rho1(1:ng,1:ng))/real(ng**2)
      print*,'    A_mesh : ',i_dim,j_dim,minval(A_mesh(i_dim,j_dim,1:ng,1:ng)),maxval(A_mesh(i_dim,j_dim,1:ng,1:ng)),sum((A_mesh(i_dim,j_dim,1:ng,1:ng))*1d0)
    enddo
    enddo ! i_dim
    deallocate(dsp)

    A_mesh(1,1,:,:) = A_mesh(1,1,:,:) + 1
    A_mesh(2,2,:,:) = A_mesh(2,2,:,:) + 1


    allocate(trace_A(ng,ng),det_A(ng,ng))
    trace_A = A_mesh(1,1,:,:) + A_mesh(2,2,:,:)
    det_A   = A_mesh(1,1,:,:)*A_mesh(2,2,:,:) - A_mesh(1,2,:,:)*A_mesh(2,1,:,:)

    allocate(kappa(ng, ng))
    kappa = 1.0D0 - 0.5D0 * trace_A
    print*,'    kappa  : ',minval(kappa),maxval(kappa),sum(kappa*1d0)
    open(16,file=output_name('k'//namespace//'_q'),status='replace',access='stream')
    write(16) kappa
    close(16)
    deallocate(kappa)
    
    allocate(lambda1(ng, ng),lambda2(ng, ng))
    lambda1 = 0.5D0 * (trace_A + SQRT(trace_A**2 - 4.0D0 * det_A))
    lambda2 = 0.5D0 * (trace_A - SQRT(trace_A**2 - 4.0D0 * det_A))
    print*,'    lambda1: ',minval(lambda1),maxval(lambda1),sum(lambda1*1d0)
    print*,'    lambda2: ',minval(lambda2),maxval(lambda2),sum(lambda2*1d0)
    open(20,file=output_name('l1'//namespace//'_q'),status='replace',access='stream')
    write(20) lambda1
    close(20)
    open(20,file=output_name('l2'//namespace//'_q'),status='replace',access='stream')
    write(20) lambda2
    close(20)
    deallocate(lambda1,lambda2,trace_A)

    allocate(mu(ng,ng))
    mu = abs(det_A)
    open(21,file=output_name('det'//namespace//'_q'),status='replace',access='stream')
    write(21) mu
    close(21)
    mu = abs(1/det_A)
    print*,'    u     : ',minval(mu),maxval(mu),sum(mu*1d0)
  open(21,file=output_name('u'//namespace//'_q'),status='replace',access='stream')
    write(21) mu
    close(21)
  open(21,file=output_name('u-1'//namespace//'_q'),status='replace',access='stream')
    write(21) mu-1
    close(21)
    deallocate(det_A,mu)


    allocate(gamma1(ng, ng),gamma2(ng, ng))
    gamma1 = 0.5D0 * (A_mesh(1,1,:,:) - A_mesh(2,2,:,:))
    gamma2 = 0.5D0 * (A_mesh(1,2,:,:) + A_mesh(2,1,:,:))
    print*,'    gamma1 : ',minval(gamma1),maxval(gamma1),sum(gamma1*1d0)
    print*,'    gamma2 : ',minval(gamma2),maxval(gamma2),sum(gamma2*1d0)
    open(17,file=output_name('g1'//namespace//'_q'),status='replace',access='stream')
    write(17) gamma1
    close(17)
    open(18,file=output_name('g2'//namespace//'_q'),status='replace',access='stream')
    write(18) gamma2
    close(18)
    deallocate(gamma1,gamma2)

    allocate(omega(ng, ng))
    omega = 0.5D0 * (A_mesh(2,1,:,:) - A_mesh(1,2,:,:))
    print*,'    omega  : ',minval(omega),maxval(omega),sum(omega*1d0)
    open(19,file=output_name('j'//namespace//'_q'),status='replace',access='stream')
    write(19) omega
    close(19)
    deallocate(omega,A_mesh)
  endsubroutine decompose_Mesh_FFT

  ! subroutine decompose_A(A, kappa, gamma1, gamma2, omega, l1, l2, det_A)
  !     IMPLICIT NONE
  !     real(8), intent(in)  :: A(2, 2)     ! 输入 2x2 矩阵
  !     real(8), intent(out) :: kappa       ! 扩张/收缩参数
  !     real(8), intent(out) :: gamma1      ! 剪切分量 1
  !     real(8), intent(out) :: gamma2      ! 剪切分量 2
  !     real(8), intent(out) :: omega       ! 旋转分量
  !     real(8), intent(out) :: l1, l2      ! 特征值
  !     real(8), intent(out) :: det_A          ! 给出行列式
    
  !     real(8) :: trace_A
  !     real(8) :: a11, a12, a21, a22
    
  !     ! 提取矩阵元素
  !     a11 = A(1, 1)
  !     a12 = A(1, 2)
  !     a21 = A(2, 1)
  !     a22 = A(2, 2)
    
  !     ! 计算迹 (trace) 和行列式 (determinant)
  !     trace_A = a11 + a22
  !     det_A = a11 * a22 - a12 * a21
  !     ! print*, det_A
    
  !     ! 计算 kappa (扩张/收缩)
  !     kappa = 1.0D0 - 0.5D0 * trace_A
    
  !     ! 计算 gamma1 和 gamma2 (剪切分量)
  !     gamma1 = 0.5D0 * (a22 - a11)
  !     gamma2 = 0.5D0 * (a12 + a21)
    
  !     ! 计算 omega (旋转分量)
  !     omega = 0.5D0 * (a21 - a12)
    
  !     ! 计算特征值 l1 和 l2
  !     l1 = 0.5D0 * (trace_A + SQRT(trace_A**2 - 4.0D0 * det_A))
  !     l2 = 0.5D0 * (trace_A - SQRT(trace_A**2 - 4.0D0 * det_A))
  ! endsubroutine decompose_A

  logical function is_bad(val)
    use, intrinsic :: ieee_arithmetic
    real, intent(in) :: val
    is_bad = (abs(val)>1e10  .or. (val == 0.0) .or. (ieee_is_nan(val)) .or. (.not. ieee_is_finite(val)))
  endfunction
end