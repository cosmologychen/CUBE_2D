
subroutine Green_2D(Gk,nglobal,n1,n2,apm_soft,apm_range,gridsize)
  use omp_lib
  use parameters
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer(8) nglobal,n1,n2
  real apm_soft,apm_range,gridsize
  real Gk(n1,n2)
  integer k1,k2,i1,i2
  real(8) fftspace(nglobal),kvec(2),u2n_2,Dvec(2),u2n,ur_cn(2),kvec_n(2)
  real(8) Uk2_n,kmag_n,ka,S_soft,S_range,Rvec_n(2),j1,j0

  call system_clock(t1,t_rate)
  print*,'  Initialize 2D_Green''s function'
  print*,'    mesh =',int(nglobal,kind=2)
  print*,'    n1,n2 =',int(n1,kind=2),int(n2,kind=2)
  print*,'    apm_soft, apm_range =',apm_soft,apm_range
  print*,'    gridsize =',gridsize
  print*,'    n_int =',int(n_int,kind=2)

  fftspace=2*pi*(1./nglobal)*(mod([(k1,k1=1,nglobal)]+nglobal/2-1,nglobal)-nglobal/2)
  
  !$omp paralleldo default(shared) schedule(dynamic)&
  !$omp& private(k2,k1,kvec,u2n_2,Dvec,u2n,ur_cn,i2,i1,kvec_n,kmag_n,Uk2_n)&
  !$omp& private(ka,S_soft,S_range,Rvec_n,j1,j0)
  do k2=1,n2
  do k1=1,n1
    kvec=[fftspace(k1),fftspace(k2)]
    ! print*,'========================================='
    ! print*, '    k1,k2 =',k1,k2
    ! print*, '    kvec =',kvec/2/pi
    u2n_2=product(1-sin(kvec/2)**2+(2./15)*sin(kvec/2)**4)
    Dvec=alpha*sin(kvec)+(1-alpha)*sin(2*kvec)/2
    
    ur_cn=0
    do i2=-n_int,n_int
    do i1=-n_int,n_int
      kvec_n=kvec+2*pi*[i1,i2]
      kmag_n=norm2(kvec_n)
      Uk2_n=product(merge(1d0,sin(kvec_n/2)/(kvec_n/2),kvec_n==0))**(2*p+2)
      !u2n=u2n+Uk2_n
      ka=apm_soft*kmag_n
      if (kmag_n==0) then
        S_soft=1
      else
        ! call Bessel(0,ka/2.,j0)
        ! call Bessel(1,ka/2.,j1)
        j0 = bessel_j0(ka/2)
        j1 = bessel_j1(ka/2)
        S_soft=(128*j1/ka-32*j0)*(ka**-2)
      endif
      if (apm_range==0) then
        S_range=0
      elseif (kmag_n==0) then
        S_range=1
      else
        ka=apm_range*kmag_n
        ! call Bessel(0,ka/2,j0)
        ! call Bessel(1,ka/2,j1)
        j0 = bessel_j0(ka/2)
        j1 = bessel_j1(ka/2)
        S_range=(128*j1/ka-32*j0)*(ka**-2)
      endif

      Rvec_n=kvec_n*(S_soft**2-S_range**2)/kmag_n**2
      ur_cn=ur_cn+Uk2_n*Rvec_n
      
    enddo
    enddo
    u2n=product(1-sin(kvec/2)**2+(2./15)*sin(kvec/2)**4);
    Gk(k1,k2)=2*pi*sum(Dvec*ur_cn)/sum(Dvec**2)/u2n**2/gridsize**(2-1)
    if (k1 > 1 .and. k2 > 1 .and. (.not. ieee_is_finite(Gk(k1,k2)))) then
      print*,'      Gk(',k1,',',k2,') =',Gk(k1,k2)
      print*,'      u2n =',u2n
      print*,'      Uk2_n =',Uk2_n
      print*,'      u2n_2 =',u2n_2
      print*,'      Dvec =',Dvec
      print*,'      ur_cn =',ur_cn
      print*,'      kvec =',kvec/2/pi
      print*,'      kvec_n =',kvec_n/2/pi
      print*,'      kmag_n =',kmag_n
      stop
    endif
  enddo
  enddo
  !$omp endparalleldo
  Gk(::nglobal/2,::nglobal/2)=0;
  ! print*,maxval(Gk),minval(Gk)
  ! stop
  

  call system_clock(t2,t_rate)
  print*,'    elapsed time =',real(t2-t1)/t_rate,'secs';
  print*,''
endsubroutine

! subroutine Bessel(a,x,j)
!   use parameters
!   use, intrinsic :: ieee_arithmetic
!   implicit none
!   real(8), intent(in)  :: x
!   integer, intent(in)  :: a
!   real(8), intent(out) :: j
!   real(16) f,ti
!   integer(8) mi,mn

!   mn = 1000

    
!   f = 2*pi/mn
!   j=0
!   ti=0
!   do mi = 1, mn
!     ti = ti+f
!     j = j + cos(a*ti-x*sin(ti))*f
!   end do
!   j = j/2/pi
  
!   if (.not. ieee_is_finite(j)) then
!     print *, 'j 是 Inf 或 NaN'
!     print *, 'a = ', a
!     print *, 'x = ', x
!     stop
!   endif
! endsubroutine
