
module powerspectrum
    use omp_lib
    use parameters
    integer(8),parameter :: nbin=nint(nyquist*sqrt(3.))

    real        rho1(nw+2,nw)
    complex     rho1k(nw/2+1,nw)
    equivalence(rho1,rho1k)

contains

subroutine cross_power(xi,rhok1,rhok2,n_particle,n_interp)
  use omp_lib
  implicit none
  integer i,j,ig,jg,ibin,n_interp
  integer(8) n_particle
  real kr,kx(2),rbin,C1k(2),Dk,amp11,amp12,amp21,amp22,xi(10,0:nbin)
  complex rhok1(nw/2+1,nw),rhok2(nw/2+1,nw)

  xi=0
  do j=1,nw
    do i=1,nw/2+1
        kx(2)=mod(j+nw/2-1,nw)-nw/2
        kx(1)=i-1
        if (j==1 .and. i==1) cycle ! zero frequency
        kr=sqrt(kx(1)**2+kx(2)**2)
        ibin=nint(kr)
        xi(1,ibin)=xi(1,ibin)+1 ! number count
        xi(2,ibin)=xi(2,ibin)+kr ! k count
        amp11=real(rhok1(i,j)*conjg(rhok1(i,j)))
        amp22=real(rhok2(i,j)*conjg(rhok2(i,j)))
        amp12=real(rhok1(i,j)*conjg(rhok2(i,j)))
        ! amp21=real(rhok2(i,j)*conjg(rhok1(i,j)))
        ! print *,'k=',kr,amp11
        if (n_interp==1) then ! NGP
        C1k=1
        elseif (n_interp==2) then ! CIC
        C1k=1-(2./3.)*sin(pi*kx/nw)**2
        elseif (n_interp==3) then ! TSC
        C1k=1-sin(pi*kx/nw)**2+(2./15.)*sin(pi*kx/nw)**4
        endif
        Dk=(C1k(1)*C1k(2))/n_particle
        ! print*,kr,amp11,amp22,amp12,Dk
        xi(3 ,ibin)=xi(3 ,ibin)+amp11 ! raw power
        xi(4 ,ibin)=xi(4 ,ibin)+(amp11-Dk) ! P_r(k)
        xi(5 ,ibin)=xi(5 ,ibin)+Dk ! P_r(k)
        xi(6 ,ibin)=xi(6 ,ibin)+amp22
        xi(7 ,ibin)=xi(7 ,ibin)+(amp22-Dk)
        xi(9 ,ibin)=xi(9 ,ibin)+abs(amp12)
        xi(10,ibin)=xi(10,ibin)+(amp12-Dk)
    enddo
    ! stop
  enddo

    xi(2,:)=xi(2,:)/xi(1,:)

    xi(3,:)=xi(3,:)/xi(1,:) ! raw power11
    xi(4,:)=xi(4,:)/xi(1,:) ! raw power11 - Dk

    xi(6,:)=xi(6,:)/xi(1,:) ! raw power22
    xi(7,:)=xi(7,:)/xi(1,:) ! raw power22 - Dk
    xi(8,:)=xi(7,:)
    xi(9 ,:)=xi(9 ,:)/xi(1,:)
    ! xi(10,:)=xi(10,:)/xi(1,:)

    call pk_correction(xi,n_interp,5,3)
    call pk_correction(xi,n_interp,5,3)
    call pk_correction(xi,n_interp,5,3)
    call pk_correction(xi,n_interp,5,3)
    call pk_correction(xi,n_interp,8,3)
    call pk_correction(xi,n_interp,8,3)
    call pk_correction(xi,n_interp,8,3)
    call pk_correction(xi,n_interp,8,3)
    ! divide and normalize
    xi(2 ,:)=xi(2 ,:)*(2 * pi)/box ! k_phys  
    xi(3 ,:)=xi(3 ,:)*(box**2) ! power11_phys
    xi(4 ,:)=xi(4 ,:)*(box**2) ! power11_phys
    xi(5 ,:)=xi(5 ,:)*(box**2) ! power11_phys
    xi(6 ,:)=xi(6 ,:)*(box**2)
    xi(7 ,:)=xi(7 ,:)*(box**2)
    xi(8 ,:)=xi(8 ,:)*(box**2)
    xi(9 ,:)=xi(9 ,:)*(box**2)
    xi(10,:)=xi(10,:)*(box**2)
    ! print*,xi(9,:)
  
endsubroutine cross_power


subroutine auto_power(xi,n_particle,n_interp)
  use omp_lib
  implicit none
  integer i,j,ig,jg,ibin,n_interp
  integer(8) n_particle
  real kr,kx(2),sincx,sincy,sinc,rbin,C1k(2),Dk,amp11,xi(10,0:nbin)

  xi=0
  do j=1,ngic
    do i=1,nyquist+1
        kx(2)=mod(j+nyquist-1,ngic)-nyquist
        kx(1)=i-1
        if (j==1 .and. i==1) cycle ! zero frequency
        kr=sqrt(kx(1)**2+kx(2)**2)
        ibin=nint(kr)
        xi(1,ibin)=xi(1,ibin)+1 ! number count
        xi(2,ibin)=xi(2,ibin)+kr ! k count
        amp11=real(rho1k(i,j)*conjg(rho1k(i,j)))
        ! print *,'k=',kr,amp11
        if (n_interp==1) then ! NGP
        C1k=1
        elseif (n_interp==2) then ! CIC
        C1k=1-(2./3.)*sin(pi*kx/nw)**2
        elseif (n_interp==3) then ! TSC
        C1k=1-sin(pi*kx/nw)**2+(2./15.)*sin(pi*kx/nw)**4
        endif
        Dk=(C1k(1)*C1k(2))/n_particle
        xi(3,ibin)=xi(3,ibin)+amp11 ! raw power
        xi(4,ibin)=xi(4,ibin)+(amp11-Dk) ! P_r(k)
    enddo
  enddo

    xi(2,:)=xi(2,:)/xi(1,:)
    xi(3,:)=xi(3,:)/xi(1,:) ! raw power
    xi(4,:)=xi(4,:)/xi(1,:) ! raw power - Dk
    xi(5,:)=xi(4,:)
    call pk_correction(xi,n_interp,5,3)
    call pk_correction(xi,n_interp,5,3)
    call pk_correction(xi,n_interp,5,3)
    call pk_correction(xi,n_interp,5,3)
    ! divide and normalize
    xi(2,:)=xi(2,:)*(2*pi)/box ! k_phys  
    xi(3,:)=xi(3,:)*(box**2) ! power_phys
    xi(4,:)=xi(4,:)*(box**2) ! power_phys
    xi(5,:)=xi(5,:)*(box**2) ! power_phys
    ! print*, 'xi(3,:)',xi(3,:)
  
endsubroutine auto_power

subroutine pk_correction(xi,p,n_xi,n_int)
  use omp_lib
  implicit none
  integer i,j,n_xi,n_int,in,jn,ibin,nplocal,icore,p
  real alpha,kvec(2),kmag,kmagn,kvecn(2),ks(2),Wk2Pk,Pk,cdata(0:nbin,0:ncore,3),xi(10,0:nbin)
  call omp_set_num_threads(ncore)
  alpha=(log(interp1(xi(2,:),xi(n_xi,:),real(nyquist)))-log(interp1(xi(2,:),xi(n_xi,:),real(nyquist)/2)))/log(2.)
  print*,'pk_correction: p,n_int =',p,n_int
  print*,'  P(k_N),P(k_N/2) =',interp1(xi(2,:),xi(n_xi-1,:),real(nyquist)),interp1(xi(2,:),xi(n_xi-1,:),real(nyquist)/2)
  print*,'  alpha =',alpha
  cdata=0
  call system_clock(t1,t_rate)
  !$omp paralleldo default(shared) schedule(dynamic)&
  !$omp& private(i,icore,j,kvec,kmag,ibin,Wk2Pk,Pk,in,jn,kvecn,kmagn,ks)
  do i=1,nyquist+1
    icore=omp_get_thread_num()+1
    do j=1,i
        kvec=[i,j]-1.0
        kmag=norm2(kvec)
        ibin=nint(kmag)
        Wk2Pk=0
        Pk=kmag**alpha
        do in=-n_int,n_int
        do jn=-n_int,n_int
          kvecn=kvec+[in,jn]*nw
          kmagn=norm2(kvecn)
          ks=pi*kvecn/nw
          Wk2Pk=Wk2Pk+(product(merge(1.,sin(ks)/ks,ks==0))**(2*p)) * (kmagn**alpha)
        enddo
        enddo
        cdata(ibin,icore,:)=cdata(ibin,icore,:)+[1.,kmag,Wk2Pk/Pk]
    enddo
  enddo
  !$omp endparalleldo
  call system_clock(t2,t_rate)
  print*, '  integration time =',real(t2-t1)/t_rate,'secs'
  cdata(1:nbin,0,:)=sum(cdata(1:nbin,1:ncore,:),dim=2)
  cdata(1:nbin,0,2)=cdata(1:nbin,0,2)/cdata(1:nbin,0,1)
  cdata(1:nbin,0,3)=cdata(1:nbin,0,3)/cdata(1:nbin,0,1)
  xi(n_xi,1:nbin)=xi(n_xi-1,1:nbin)/cdata(1:nbin,0,3)
endsubroutine


real function interp1(xdata,ydata,xq)
    implicit none
    integer(4) i_mid,i1,i2
    real xdata(nbin),ydata(nbin),xq
    i1=1; i2=nbin
    do while (i2-i1>1)
        i_mid=(i1+i2)/2
        if (xq>xdata(i_mid)) then
        i1=i_mid
        else
        i2=i_mid
        endif
    enddo
    interp1=ydata(i1)+(xq-xdata(i1))/(xdata(i2)-xdata(i1))*(ydata(i2)-ydata(i1))
endfunction

endmodule