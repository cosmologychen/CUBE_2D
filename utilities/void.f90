program new_void
  use parameters
  implicit none

  real,parameter :: b_link=0.20 ! linking length for FoF
  real :: min_r,  b_link_r,np_void_min
  integer :: void_np,  nvoids, nvoids_r

  integer(8) i,j,ii,jj,l,cur_checkpoint,ip,ip1,ip2,np_void,idl(2),cm(2),rm,iq(2)
  real,allocatable :: xp(:,:),triangles(:,:,:),down_r(:),rho(:,:)
  integer(8),allocatable :: hoc(:,:),ll(:),counts(:,:),ll_void(:),llgp(:),hcgp(:),ecgp(:),iph_void_all(:),iph_void(:)
  integer :: iunit, iostat
  character(len=400) :: cmd
  integer :: unit, filesize, num_reals, ntri ,t ,n_void,np_low

  integer np_head,np_iso,np_mem,halo_np
  real :: xp_mean(2),center(2),radius, a(2), b(2), c(2), r2,x2,u2, d2,xpos(2),qp_mean(2)
  integer :: tmax(2),tmin(2),current_pos

  min_r = 1
  b_link_r = 2
  np_void_min = 50

  open(16,file='../z_checkpoint.txt',status='old')
  do i=1,nmax_redshift-1
    read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
  enddo
  71 n_checkpoint=i-1
  close(16)
  if (n_checkpoint==0) stop 'z_checkpoint.txt empty'
  
  sim%cur_checkpoint = 5

  open(11,file=output_name('info'),access='stream'); read(11) sim; close(11)
  print*, 'np =',sim%np
  allocate(xp(2,sim%np))
  open(11,file=output_name('xp'),access='stream'); read(11) xp; close(11)

  allocate(ll(sim%np),hoc(ng,ng),counts(ng,ng),ll_void(sim%np/4),llgp(sim%np),hcgp(sim%np),ecgp(sim%np),down_r(sim%np))
  ll = 0; hoc = 0; counts = 0; ll_void=0; llgp=0; hcgp=0; ecgp=0


  allocate(rho(0:ng+1,0:ng+1))
  open(15,file=output_name('uD_q'),status='old',access='stream')
  read(15) rho(1:ng,1:ng)
  close(15)

  open(11,file=output_name('np_low_xp'),status='replace',access='stream')
  np_low = 0
  rm = min_r
  do ip=1,sim%np
    xpos=xp(:,ip)
    idl=xpos2mesh(xpos,ng)
    ll(ip)=hoc(idl(1),idl(2))
    hoc(idl(1),idl(2))=ip
    counts(idl(1),idl(2))=counts(idl(1),idl(2))+1
  enddo ! ip
  close(11)

!   do ip=1,sim%np
!     xpos=xp(:,ip)
!     idl=xpos2mesh(xpos,ng)
!     ll(ip)=hoc(idl(1),idl(2))
!     hoc(idl(1),idl(2))=ip
!     counts(idl(1),idl(2))=counts(idl(1),idl(2))+1
!   enddo ! ip


!   open(11,file=output_name('np_low_xp'),status='replace',access='stream')
!   np_low = 0
!   rm = min_r
!   do i = 1, ng
!   do j = 1, ng
!     if (counts(i,j) == 1) then
!       ip1 = hoc(i,j)
!       do ii = -rm,rm
!       do jj = -rm,rm
!         if (ii == 0 .and. jj == 0) cycle
!         idl(1) = mod(i+ii+ng-1,ng)+1
!         idl(2) = mod(j+jj+ng-1,ng)+1

!         if (counts(idl(1),idl(2)) == 0 ) cycle

!         if (counts(idl(1),idl(2)) > 2 ) goto 73
!         ip2 = hoc(idl(1),idl(2))
!         do while (ip2 /= 0)
!           ! print*, ip2
!           if (sum(pbc_vec(xp(:,ip1)-xp(:,ip2))**2) < min_r) goto 73
!           ip2 = ll(ip2)
!         enddo
!       enddo
!       enddo
!       np_low = np_low + 1
!       ll_void(np_low) = ip1
!       hcgp(ip1) = ip1
!       counts(i,j) = -1
!       write(11) xp(:,ip1)-0.5
! 73  continue
!     endif

!   enddo
!   enddo
!   close(11)
     
  print*, 'np_low =', np_low,real(np_low)/real(sim%np)
  ! stop

  open(11,file=output_name('void_test'),status='replace',access='stream')
  down_r = min_r**2
  d2 = min_r**2
  rm = b_link_r
  u2 = b_link_r**2
  ip2 = 0
  do l = 1, np_low
    ip1 = ll_void(l)
    ! if (ip2 > 0) ip1 = ip2
    ! if (ip2 == 0) print*,'new ',ip1 
    center = xp(:,ip1)
    cm = xpos2mesh(center,ng)
    ! print*, 'ip cm', ip1,cm,center
    ip2 = 0
    ! u2 = b_link_r**2
    do i=-rm-1,rm+1
    do j=-rm-1,rm+1
      if (i ==0 .and. j == 0) cycle
      idl(1) = mod(cm(1)+i+ng-1,ng)+1
      idl(2) = mod(cm(2)+j+ng-1,ng)+1
      if (counts(idl(1),idl(2)) == -1) then
        ip = hoc(idl(1),idl(2))
        ! d2 = min(down_r(ip1),down_r(ip))+1e-5
        r2 = sum(pbc_vec(center-xp(:,ip))**2)
        ! print*,ip, d2, r2, u2
        if ((r2 < u2) .and. (r2 > d2)) then
          ip2 = ip
          ! u2  = r2
          call merge_chain(ip1,ip2,d2)
        endif
      endif
    enddo
    enddo
    ! if (ip2 > 0) then
    !   write(11) xp(:,ip1) - 0.5
    !   ! print*, '  closest', ip1,ip2, u2
    !   call merge_chain(ip1,ip2,d2)
    ! endif
    ! print*,''
    ! if (ip2 > 0) call merge_chain(ip1,ip2,d2)
  enddo


  if (ip2 > 0) write(11) xp(:,ip2) - 0.5
  close(11)
  print*, 'linking done',rm
  ! stop


  allocate(iph_void_all(sim%np),iph_void(sim%np))
  np_iso = 0; np_head = 0; iph_void=0; iph_void_all=0
  do i=1,sim%np
    if (hcgp(i)==i) then
      np_iso=np_iso+1
    elseif (hcgp(i)/=0) then

      ip1 = hcgp(i)
      void_np = 0
      do while (ip1 /= 0)
        void_np = void_np + 1
        ip1 = llgp(ip1)
      enddo
      if (void_np > np_void_min) then
        np_head=np_head+1
        iph_void_all(np_head) = hcgp(i)
        iph_void(np_head) = void_np
      endif
    endif
  enddo
  ! stop

  if (np_head <= 1) stop 'np_head <= 1'
  print*, 'np_iso, np_head', np_iso, np_head,maxval(iph_void(:np_head)),minval(iph_void(:np_head)),sum(iph_void)+np_iso

  call indexed_sort(np_head,-iph_void(:np_head),ecgp(:np_head))

  open(11,file=output_name('void'),status='replace',access='stream')
  open(12,file=output_name('void_xp_mean_only'),status='replace',access='stream')
  open(13,file=output_name('void_qp_mean_only'),status='replace',access='stream')
  write(11) b_link,np_void_min,np_head
  do j=1,np_head
    i = ecgp(j)
    ip = iph_void_all(i)
    void_np = iph_void(i)
    print*, j,'void', i,void_np
    xp_mean = 0
    ip1 = ip
    write(11) void_np
    do while (ip1 /= 0)
      xp_mean = xp_mean + pbc_vec(xp(:,ip1)-xp(:,ip))
      write(11) xp(:,ip1)-0.5
      ip1 = llgp(ip1)
    enddo  
    xp_mean = wrap_position2(xp_mean / void_np+xp(:,ip))
    ! print*,'  ',xp_mean,void_np

    ip1 = iph_void_all(i)
    qp_mean = 0
    do while (ip1 /= 0)
      iq(1)=(ip1-1)/ng
      iq(2)=modulo(ip1-1,int(ng,4))
      xpos=iq+0.5
      qp_mean = qp_mean + pbc_vec(xpos-xp_mean)
      ! write(11) xpos
      ip1 = llgp(ip1)
    enddo
    qp_mean = wrap_position2(qp_mean / void_np + xp_mean)
    write(11) xp_mean
    write(12) xp_mean
    write(13) qp_mean

    ! print*,iph_void_all(i), 'np', void_np, 'xp_mean', xp_mean
  enddo
  close(11)
  close(12)
  close(13)

  deallocate(xp,ll,llgp,hcgp,ecgp,hoc)
  print*,output_name('void_xp_mean')


  contains





  subroutine merge_chain(ii,jj,d2)
    integer(8) ii,jj,ihead,jhead,iend,jend,ipart
    real d2

    down_r(ii) = d2 
    ! down_r(jj) = d2
    jend=merge(jj,ecgp(jj),ecgp(jj)==0)
    iend=merge(ii,ecgp(ii),ecgp(ii)==0)
    if (iend==jend) return ! same chain
    ihead=max(hcgp(ii),hcgp(iend))
    jhead=max(hcgp(jj),hcgp(jend))
    ipart=jhead ! change eoc of the chain-j
    do while (ipart/=0)
      ecgp(ipart)=iend ! set chain-j's eoc to be iend
      ipart=llgp(ipart) ! next particle
    enddo
    llgp(jend)=ihead ! link j group to i group
    ecgp(ii)=iend
    hcgp(iend)=jhead ! change hoc
    hcgp(jend)=0 ! set jend as a member
  endsubroutine merge_chain

  subroutine indexed_sort(N,ARRIN,INDX)
    implicit none
    integer(4) N ! number of voids to sort
    integer(4) IR
    integer(8) INDX(:),INDXT,J,L,I
    integer(8) ARRIN(N),Q

    IF (N.LE.1) RETURN

    DO 11 J=1,N
      INDX(J)=J
    11    CONTINUE
    L=N/2+1
    IR=N
    10    CONTINUE
      IF(L.GT.1)THEN
        L=L-1
        INDXT=INDX(L)
        Q=ARRIN(INDXT)
      ELSE
        INDXT=INDX(IR)
        Q=ARRIN(INDXT)
        INDX(IR)=INDX(1)
        IR=IR-1
        IF(IR.EQ.1)THEN
          INDX(1)=INDXT
          RETURN
        ENDIF
      ENDIF
      I=L
      J=L+L
    20      IF(J.LE.IR)THEN
        IF(J.LT.IR)THEN
          IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
        ENDIF
        IF(Q.LT.ARRIN(INDX(J)))THEN
          INDX(I)=INDX(J)
          I=J
          J=J+J
        ELSE
          J=IR+1
        ENDIF
      GO TO 20
      ENDIF
      INDX(I)=INDXT
    GO TO 10
  endsubroutine indexed_sort




end program new_void