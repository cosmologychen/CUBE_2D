program CUBE_FoF
  use parameters
  implicit none

  real,parameter :: density = 2.5
  integer(4) ij(2,n_neighbor)
  integer(8) i,j,l,cur_checkpoint,ip,ip1,ip2,i_neighbor,idl(2),iq(2)
  integer(8),allocatable :: hoc(:,:),ll(:),llgp(:),hcgp(:),ecgp(:),iph_halo_all(:),iph_halo(:)

  integer(8) np,number_halo
  real(4),allocatable :: xv(:,:)

  integer np_head,np_iso,np_mem,halo_np
  real(4) xp_mean(2),xpos(2),qp_mean(2)

  real r2,link2

  open(16,file='../z_checkpoint.txt',status='old')
  do i=1,nmax_redshift-1
    read(16,end=71,fmt='(f8.4)') z_checkpoint(i)
  enddo
  71 n_checkpoint=i-1
  close(16)
  if (n_checkpoint==0) stop 'z_checkpoint.txt empty'

  print*,'  initialize FoF cell neighbors'
  l=0
  do j=-nrange,-1
  do i=-nrange,nrange
    l=l+1 
    ij(:,l)=[i,j]
  enddo
  enddo
  j=0
  do i=-nrange,-1
    l=l+1
    ij(:,l)=[i,j]
  enddo

  link2 = b_link**2
  ! print*, 'link2 =',link2
  ! do i_neighbor = 1,n_neighbor
  !   print*,ij(:,i_neighbor)
  ! enddo
  ! stop


  do cur_checkpoint= 5,5!1,n_checkpoint
    sim%cur_checkpoint=cur_checkpoint
    print*, ''
    print*, '==========================================='
    print*, '==========================================='
    print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
    !print*,output_name('info')
    open(11,file=output_name('info'),access='stream'); read(11) sim; close(11)
    np=sim%np
    print*, 'np =',np
    allocate(xv(4,np))
    open(11,file=output_name('xp'),access='stream'); read(11) xv(1:2,:); close(11)
    open(12,file=output_name('vp'),access='stream'); read(12) xv(3:4,:); close(12)

    allocate(ll(sim%np),llgp(sim%np),hcgp(sim%np),ecgp(sim%np),hoc(0:ng+1,0:ng+1))
    ll = 0 
    do ip=1,sim%np
        xpos=xv(1:2,ip)
        idl=xpos2mesh(xpos,ng)
        ll(ip)=hoc(idl(1),idl(2))
        hoc(idl(1),idl(2))=ip
        hcgp(ip)=ip
    enddo ! ip
    
    hoc(0,:)     = hoc(ng,:)
    hoc(ng+1,:) = hoc(1,:)
    hoc(:,0)     = hoc(:,ng)
    hoc(:,ng+1) = hoc(:,1)

    print*, '  hoc done'


    llgp=0; ecgp=0;

    do i = 1,ng
    do j = 1,ng
      ip1 = hoc(i,j)
      do while (ip1 /= 0)
        ip2 = ll(ip1)
        do while (ip2 /= 0)
          r2 = sum(pbc_vec(xv(1:2,ip1)-xv(1:2,ip2))**2)
          if (r2 <= link2) call merge_chain(ip1,ip2)
          ip2 = ll(ip2)
        enddo
        do i_neighbor = 1,n_neighbor
          ip2 = hoc(i+ij(1,i_neighbor),j+ij(2,i_neighbor))
          do while (ip2 /= 0)
            r2 = sum(pbc_vec(xv(1:2,ip1)-xv(1:2,ip2))**2)
            if (r2 <= link2) call merge_chain(ip1,ip2)
            ip2 = ll(ip2)
          enddo
        enddo
        ip1 = ll(ip1)
      enddo
    enddo
    enddo

    print*,'fof done'


    ! np_iso=0; np_mem=0; np_head=0;
    ! do i=1,sim%np
    !   if (hcgp(i)==i) then
    !     np_iso=np_iso+1
    !   elseif (hcgp(i)==0) then
    !     np_mem=np_mem+1
    !   else
    !     np_head=np_head+1
    !   endif
    ! enddo
    ! print*,'N_iso,mem,head =',np_iso,np_mem,np_head

    ! ! stop

    allocate(iph_halo_all(sim%np),iph_halo(sim%np))
    np_iso = 0; np_head = 0
    do i=1,sim%np
      if (hcgp(i)==i) then
        np_iso=np_iso+1
      elseif (hcgp(i)/=0) then

        ip1 = hcgp(i)
        halo_np = 0
        xp_mean = 0
        do while (ip1 /= 0)
          halo_np = halo_np + 1
          ip1 = llgp(ip1)
        enddo
        if (halo_np > np_halo_min) then
          np_head=np_head+1
          iph_halo_all(np_head) = hcgp(i)
          iph_halo(np_head) = halo_np
        endif
      endif
    enddo
     
    print*, 'np_iso, np_head', np_iso, np_head,maxval(iph_halo(:np_head)),minval(iph_halo(:np_head))

    call indexed_sort(np_head,-iph_halo(:np_head),ecgp(:np_head))

    number_halo = min((density*1e-4)**(2/3)*box*box,real(np_head))
    print*, 'number_halo',number_halo,iph_halo(ecgp(1))

    open(11,file=output_name('halo'),status='replace',access='stream')
    open(12,file=output_name('halo_xp_mean_only'),status='replace',access='stream')
    open(13,file=output_name('halo_qp_mean_only'),status='replace',access='stream')
    write(11) b_link,np_halo_min,np_head
    do j=1,number_halo
      i = ecgp(j)
      print*, j,'halo', i
      ip1 = iph_halo_all(i)
      halo_np = iph_halo(i)
      xp_mean = 0
      do while (ip1 /= 0)
        xp_mean = xp_mean + xv(1:2,ip1)
        write(11) xv(1:2,ip1)
        ip1 = llgp(ip1)
      enddo
      xp_mean = xp_mean / halo_np

      ip1 = iph_halo_all(i)
      qp_mean = 0
      do while (ip1 /= 0)
        iq(1)=(ip1-1)/ng
        iq(2)=modulo(ip1-1,int(ng,4))
        xpos=iq+0.5
        qp_mean = qp_mean + pbc_vec(xpos-xp_mean)
        write(11) xv(1:2,ip1)
        ip1 = llgp(ip1)
      enddo
      qp_mean = wrap_position2(qp_mean / halo_np + xp_mean)
      write(11) halo_np
      write(11) xp_mean
      write(12) xp_mean
      write(13) qp_mean

      print*,iph_halo_all(i), 'np', halo_np, 'xp_mean', xp_mean
    enddo
    close(11)
    close(12)
    close(13)

    deallocate(xv,ll,llgp,hcgp,ecgp,hoc)
    print*,output_name('halo_xp_mean')

  enddo

  contains





  subroutine merge_chain(ii,jj)
    integer(8) ii,jj,ihead,jhead,iend,jend,ipart
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
    integer(4) N ! number of halos to sort
    integer(4) IR
    integer(8) INDX(:),INDXT,J,L,I
    integer(8) ARRIN(N),Q

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



end program CUBE_FoF
