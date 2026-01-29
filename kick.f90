subroutine kick
    use omp_lib
    use variables
    use, intrinsic :: ieee_arithmetic
    implicit none

    logical,parameter :: PP=.true.

    integer itile,nl(2),nh(2),i,j,it,nptile,ix,iy,ip,ip1,ip2,i1,i2,idl(2),idx(2,3)
    integer ipeny,ipenx,ny,nx,i_iter,n_iter,i_neighbor,ncount,tshift(2)!,Pm1count(ncore)
    real,allocatable :: rho_th(:,:,:),phi_th(:,:),force_th(:,:,:,:)
    real,allocatable :: phi(:,:),force(:,:,:)
    real,allocatable :: af(:,:)
    real(8) pt,fpp(2)
    real rmag,rvec(2)
    real dx(2,p+1),l2(2),dv(2),xpos(2)
    integer,allocatable :: ll(:),hoc(:,:)!,ll2(:)

    real :: x1(2)
    integer :: x_wrapped1(2)

    ! allocate(rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt),rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt,ncore))
    ! allocate(rhoc(nt,nt,nnt,nnt),rhoc_th(nt,nt,nnt,nnt,ncore))
    ! rhoc_th = 0
    ! rhoc  = 0
    ! Pm1count = 0
    print*,'vmax',maxval(abs(vp(1,:))),maxval(abs(vp(2,:)))
    

    allocate(ll(sim%np),hoc(1-ncb:nc+ncb,1-ncb:nc+ncb))!,ll2(sim%np))
    ! ll2=0
    ll = 0
    hoc = 0
    
    do ip=1,sim%np
        ! xp(:,ip) = [255.995865,180.817261]
        xpos=xp(:,ip)/ratio_cs
        idl=xpos2mesh(xpos,nc)
        if (maxval(idl) > nc  .or. minval(idl) < 1) then
            print*, 'xp out of range in kick ll',ip
            print*, 'xp(1,ip) =',xp(1,ip),'  xp(2,ip) =',xp(2,ip)
            print*, 'xpos =',xpos
            print*, 'fxpos =',floor(xpos)
            print*, 'idl =',idl

            x1 = xpos
            print*, ''
            print*, ''
            print*, ''
            do i=1,2
                print*,i,'x(i)=',x1(i),nc
                x_wrapped1(i) = floor(wrap_position(x1(i)))
                print*,'  x_wrapped1(i)=',x_wrapped1(i),x_wrapped1(i) >= nc,x_wrapped1(i) <  0.0
                if (x_wrapped1(i) >= nc ) x_wrapped1(i) = x_wrapped1(i) - nc  
                if (x_wrapped1(i) <  0.0) x_wrapped1(i) = x_wrapped1(i) + nc
                print*,'  x_wrapped1d(i)=',x_wrapped1(i)
                x_wrapped1(i) = x_wrapped1(i)+1
            enddo
            stop 
        endif
        ll(ip)=hoc(idl(1),idl(2))
        hoc(idl(1),idl(2))=ip
    enddo ! ip
    
    hoc(1-ncb:0,:) = hoc(nc-ncb+1:nc,:)
    hoc(nc+1:,:)   = hoc(1:ncb,:)
    hoc(:,1-ncb:0) = hoc(:,nc-ncb+1:nc)
    hoc(:,nc+1:)   = hoc(:,1:ncb)
    if (maxval(ll) > ng**2 .or. minval(ll) < 0 .or. maxval(hoc) > ng**2 .or. minval(hoc) < 0) then
        stop 'll or hoc out of range in kick ll'
    endif

    ! ip = 0
    ! do i = 1,nc
    ! do j = 1,nc
    ! ip1=hoc(i,j)
    ! do while(ip1/=0) ! particle A
    !     ip = ip + 1
    !     ip1=ll(ip1)
    ! enddo
    ! enddo
    ! enddo
    ! print*,'  np =',ip
    ! stop



    print*,'PM1' ! =====================================================
    call tic(11)
        allocate(rho_th(1-p:nw+p,1-p:nw+p,ncore))
        rho_th=0
        !$omp paralleldo default(shared) schedule(dynamic)&
        !$omp& private(ip,xpos,iteam,idx,dx,l2,i1,i2)
        do ip=1,sim%np
            iteam=omp_get_thread_num()+1
            xpos=xp(:,ip)/ratio_cs-0.5
            idx(:,2)=floor(xpos)+1
            idx(:,1)=idx(:,2)-1
            idx(:,3)=idx(:,2)+1
            l2=xpos-floor(xpos)
            dx(:,1)=(1-l2)**2/2
            dx(:,3)=l2**2/2
            dx(:,2)=1-dx(:,1)-dx(:,3)
            if (maxval(idx) > nw+p  .or. minval(idx) < 1-p ) then
                print*, 'xp out of range in PM1 rho_th',ip
                print*, 'xp(1,ip) =',xp(1,ip),'  xp(2,ip) =',xp(2,ip)
                print*, 'idx =',idx
                print*, 'idx range',1-p,nw+p
                stop 
            endif
            do i2=1,p+1
            do i1=1,p+1
                rho_th(idx(1,i1),idx(2,i2),iteam)=rho_th(idx(1,i1),idx(2,i2),iteam)+dx(1,i1)*dx(2,i2)
            enddo
            enddo
            ! Pm1count(iteam) = Pm1count(iteam) + 1
        enddo ! ip
        !$omp endparalleldo
        
        rho1 = 0
        do iteam=1,ncore
            ! print*,'PM1 rho_th  ',iteam,sum(rho_th(:,:,iteam)),Pm1count(iteam),sum(rho_th(:,:,iteam))-Pm1count(iteam)
            ! print*,'PM1 rho_th0  ',iteam,sum(rho_th(:nw,:nw,iteam)),Pm1count(iteam),sum(rho_th(:nw,:nw,iteam))-Pm1count(iteam)
            rho_th(nw-p+1:nw,:,iteam) = rho_th(nw-p+1:nw,:,iteam) + rho_th(:0,:,iteam)
            rho_th(:,nw-p+1:nw,iteam) = rho_th(:,nw-p+1:nw,iteam) + rho_th(:,:0,iteam)
            rho_th(1:p,:,iteam)       = rho_th(1:p,:,iteam)       + rho_th(nw+1:nw+p,:,iteam)
            rho_th(:,1:p,iteam)       = rho_th(:,1:p,iteam)       + rho_th(:,nw+1:nw+p,iteam)
            ! print*,'PM1 rho_th  ',iteam,sum(rho_th(1:nw,1:nw,iteam)),Pm1count(iteam),sum(rho_th(1:nw,1:nw,iteam))-Pm1count(iteam)

            rho1(:nw,:nw) = rho1(:nw,:nw) + rho_th(1:nw,1:nw,iteam)
        enddo
        ! stop
        deallocate(rho_th)!,rhoc_th)
        ! print*,'PM1 rho ',sum(rho1(:nw,:nw)),sum(Pm1count),maxval(rho1(:nw,:nw)),minval(rho1(:nw,:nw))
        ! print*, 'Write rho1 into',output_name_step('rho1')
        ! open(11,file=output_name_step('rho1'),status='replace',access='stream')
        ! write(11) rho1(:nw,:nw)!-16
        ! close(11)
        ! stop
        call sfftw_execute( plan)
        rho1k=rho1k*Gk1
        call sfftw_execute(iplan)
        rho1 = rho1 /(nc*nc)

        allocate(phi(-1-p:nw+p+2,-1-p:nw+p+2))
        phi(1:nc,1:nc)=rho1(1:nc,1:nc)
        phi(:0,:)=phi(nc-p-1:nc,:)
        phi(nc+1:,:)=phi(1:p+2,:)
        phi(:,:0)=phi(:,nc-p-1:nc)
        phi(:,nc+1:)=phi(:,1:p+2)
        ! print*, 'Write phi1 into',output_name_step('phi1')
        ! open(11,file=output_name_step('phi1'),status='replace',access='stream')
        ! write(11) phi
        ! close(11)
        ! stop     

        allocate(force(2,1-p:nw+p,1-p:nw+p))
        call phiential_to_force(force,phi,-1-p,nc+p+2,1-p,nw+p)
        deallocate(phi)
        pm%f2max=maxval(sum(force(:,1:nc,1:nc)**2,1))
        !$omp paralleldo default(shared) schedule(dynamic)&
        !$omp& private(ip,xpos,idx,dx,l2,i1,i2)
        do ip=1,sim%np
            ! iteam=omp_get_thread_num()+1
            xpos=xp(:,ip)/ratio_cs - 0.5
            idx(:,2)=floor(xpos)+1
            idx(:,1)=idx(:,2)-1
            idx(:,3)=idx(:,2)+1
            l2=xpos-floor(xpos)
            dx(:,1)=(1-l2)**2/2
            dx(:,3)=l2**2/2
            dx(:,2)=1-dx(:,1)-dx(:,3)
            dv=0
            do i2=1,p+1
            do i1=1,p+1
                dv=dv+force(:,idx(1,i1),idx(2,i2))*dx(1,i1)*dx(2,i2)
            enddo
            enddo
            if (maxval(idx) > nw+p  .or. minval(idx) < 1-p ) then
                print*, 'xp out of range in PM1 vp',ip
                print*, 'xp(1,ip) =',xp(1,ip),'  xp(2,ip) =',xp(2,ip)
                print*, 'idx =',idx
                print*, 'idx range',1-p,nw+p
                stop 
            endif
            vp(:,ip)=vp(:,ip)+dv*a_mid*dt*G_grid
        enddo ! ip
        !$omp endparalleldo
        deallocate(force)
    call toc(11)
        
    vmax = [maxval(abs(vp(1,:))),maxval(abs(vp(2,:)))]
    sim%dt_pm1=0.5*sqrt( 1. / (sqrt(pm%f2max)*a_mid*G_grid) )
    print*,'  pm%f2max =',pm%f2max; print*,'  vmax =',vmax
    print*,'  real time =',tcat(11,istep),'secs'; print*,''



    print*,'PM2' ! =====================================================
    call tic(12)
        vmax=0; f2max=0; pt=0;f2max_team=0;vmax_team=0
        ! allocate(rho_th(1-ngb-p:ngp+ngb+p,1-ngb-p:ngp+ngb+p,ncore))
        allocate(force_th(2,1-p:ngp+p,1-p:ngp+p,ncore))


        !$omp paralleldo schedule(dynamic) default(shared) &
        !$omp& private(itile,iteam,ix,iy,i,j,ip,xpos,idx,dx,i1,i2,l2,dv,tshift)
        do itile=0,nnt**2-1
            iteam=omp_get_thread_num()+1
            itx = itile/nnt
            ity = mod(itile,nnt)
            ix = itx*nt
            iy = ity*nt
            ! rho_th(:,:,iteam) = 0
            rho2(:,:,iteam) = 0
            ! print*, 'tile', itile, 'team', iteam, 'itx', itx, 'ity', ity
            do j = iy-ncb+1,iy+nt+ncb
            do i = ix-ncb+1,ix+nt+ncb
                tshift = -[ix,iy]*ratio_cs
                if (i < 1  ) tshift(1) = tshift(1)-ng
                if (i > nc ) tshift(1) = tshift(1)+ng

                if (j < 1  ) tshift(2) = tshift(2)-ng
                if (j > nc ) tshift(2) = tshift(2)+ng
                ! print*, 'i,j,tshift', i,j,tshift

                ip = hoc(i,j)
                ! print*,'ip',ip,i,j
                do while(ip/=0)
                    xpos=xp(:,ip)+tshift-0.5
                    idx(:,2)=floor(xpos)+1
                    idx(:,1)=idx(:,2)-1
                    idx(:,3)=idx(:,2)+1
                    l2=xpos-floor(xpos)
                    dx(:,1)=(1-l2)**2/2
                    dx(:,3)=l2**2/2
                    dx(:,2)=1-dx(:,1)-dx(:,3)
                    if (maxval(idx) > ngp+ngb+p  .or. minval(idx) < 1-ngb-p ) then
                        print*, 'itile,iteam',itile+1,iteam,tshift
                        print*, 'itx,ity,ix,iy',itx,ity,ix,iy,i,j
                        print*, 'ip',ip,-ncb+1,nt+ncb
                        print*, 'xp_c =',xp(:,ip)/ratio_cs,floor(xp(:,ip)/ratio_cs)+1
                        print*, 'xp_ct =',xpos
                        print*, 'idx =',idx
                        print*, 'idx range',1-ngb-p,ngp+ngb+p
                        stop 'xp out of range in PM2 rho_th'
                    endif
                    do i2=1,p+1
                    do i1=1,p+1
                        ! rho_th(idx(1,i1),idx(2,i2),iteam)=rho_th(idx(1,i1),idx(2,i2),iteam) + dx(1,i1)*dx(2,i2)
                        if (maxval(idx) > ngp+ngb  .or. minval(idx) < 1-ngb ) cycle
                        rho2(idx(1,i1),idx(2,i2),iteam)=rho2(idx(1,i1),idx(2,i2),iteam) + dx(1,i1)*dx(2,i2)
                    enddo
                    enddo
                    ip=ll(ip)
                enddo
            enddo
            enddo
            ! ! print *,'rho_th done'
            ! do i = 1-ngb,ngp+ngb
            !     rho2(1-ngb:ngp+ngb,i,iteam) = rho_th(1-ngb:ngp+ngb,i,iteam)
            ! enddo

            ! if (itile == 1) then
            !     ! print*, 'Write rho_th into',output_name_step('rho_th')
            !     ! open(11,file=output_name_step('rho_th'),status='replace',access='stream')
            !     ! write(11) rho_th(:,:,iteam)
            !     ! close(11)
            !     print*, 'Write rho2 into',output_name_step('rho2')
            !     open(11,file=output_name_step('rho2'),status='replace',access='stream')
            !     write(11) rho2(1-ngb:ngp+ngb,1-ngb:ngp+ngb,iteam)
            !     close(11)
            ! endif
            
            
            call sfftw_execute( plan2(iteam))
            ! print*,itile,maxval(rho2(:,:,iteam)),minval(rho2(:,:,iteam))
            rho2k(:,:,iteam)=rho2k(:,:,iteam)*Gk2
            ! print*,itile,maxval(rho2(:,:,iteam)),minval(rho2(:,:,iteam))
            call sfftw_execute(iplan2(iteam))
            rho2(:,:,iteam) = rho2(:,:,iteam)/(ngt*ngt)
            ! print*,'fft done'

            ! if (itile == 1) then
            !     print*, 'Write phi2 into',output_name_step('phi2')
            !     open(11,file=output_name_step('phi2'),status='replace',access='stream')
            !     write(11) rho2(1-ngb:ngp+ngb,1-ngb:ngp+ngb,iteam)
            !     close(11)
            ! endif

            call phiential_to_force(force_th(:,:,:,iteam),rho2(1-ngb:ngp+ngb,1-ngb:ngp+ngb,iteam),1-ngb,ngp+ngb,1-p,ngp+p)
            ! print*,'force done',itile,iteam
            f2max_team(iteam)=max(f2max_team(iteam),maxval(sum(force_th(:,1:ngp,1:ngp,iteam)**2,1)))

            tshift = -[ix,iy]*ratio_cs
            do j = iy+1,iy+nt
            do i = ix+1,ix+nt
                ! if (i < 1  ) tshift(1) = tshift(1) - ng
                ! if (j < 1  ) tshift(2) = tshift(2) - ng
                ! if (i > nc ) tshift(1) = tshift(1) + ng
                ! if (j > nc ) tshift(2) = tshift(2) + ng

                ip = hoc(i,j)
                do while(ip/=0)
                    xpos=xp(:,ip)+tshift-0.5
                    idx(:,2)=floor(xpos)+1
                    idx(:,1)=idx(:,2)-1
                    idx(:,3)=idx(:,2)+1
                    l2=xpos-floor(xpos)
                    dx(:,1)=(1-l2)**2/2
                    dx(:,3)=l2**2/2
                    dx(:,2)=1-dx(:,1)-dx(:,3)
                    if (maxval(idx) > ngp+p  .or. minval(idx) < 1-p ) then
                        print*, 'itile,iteam',itile+1,iteam,tshift
                        print*, 'itx,ity,ix,iy',itx,ity,ix,iy,i,j
                        print*, 'ip',ip,-ncb+1,nt+ncb
                        print*, 'xp_c =',xp(:,ip)/ratio_cs,floor(xp(:,ip)/ratio_cs)+1
                        print*, 'xp_ct =',xpos
                        print*, 'idx =',idx
                        print*, 'idx range',1-p,ngp+p
                        stop 'xp out of range in PM2 rho_th'
                    endif
                    dv = 0
                    do i2=1,p+1
                    do i1=1,p+1
                        dv=dv+force_th(:,idx(1,i1),idx(2,i2),iteam)*dx(1,i1)*dx(2,i2)
                    enddo
                    enddo
                    !!! !omp atomic
                    vp(:,ip)=vp(:,ip)+dv*a_mid*dt*G_grid
                    ip=ll(ip)
                enddo
            enddo
            enddo

            ! stop
        enddo
        !$omp endparalleldo
        deallocate(force_th)
    call toc(12)

    f2max=maxval(f2max_team); vmax = [maxval(abs(vp(1,:))),maxval(abs(vp(2,:)))]
    sigma_vi=sigma_vi_new;
    sim%vsim2phys=(1.5/sim%a)*box*h0*100.*sqrt(omega_m)/ng
    sim%dt_pm2=0.5*sqrt( 1. / (sqrt(f2max)*a_mid*G_grid) )
    print*,'  vmax =',vmax; print*,'  f2max =',f2max
    print*,'  real time =',tcat(12,istep),'secs'; print*,''



    if (PP) then ! =====================================================
        print*,'PP'
        call tic(14)
            vmax=0; f2max=0; pt=0;f2max_team=0;vmax_team=0
            allocate(af(2,sim%np));af=0
            ! self pp
            !$omp paralleldo default(shared) schedule(dynamic)&
            !$omp& private(j,i,ip1,ip2,rvec,rmag,fpp)
            do j=1,nc
            do i=1,nc
            ip1=hoc(i,j)
            do while(ip1/=0) ! particle A
                ip2=ll(ip1)
                do while (ip2/=0) ! particle B in the same cell
                rvec=xp(:,ip2)-xp(:,ip1); rmag=norm2(rvec)
                if (0.<rmag .and. rmag<apm2) then
                    fpp=rvec/rmag*(F_ra(rmag,app)-F_ra(rmag,apm2))
                    ! print*,fpp
                    if (any(.not. ieee_is_finite(fpp))) then
                        print*,'fpp is NaN'
                        print*,xp(:,ip1),xp(:,ip2)
                        print*,rvec,rmag
                        print*,fpp
                        print*,ip1,ip2
                        print*,i,j,ipenx,ipeny
                        print*,itile
                        print*,iteam
                        print*,hoc(i,j)
                        print*,ll(ip1)
                        stop
                    endif
                    af(:,ip1)=af(:,ip1)+fpp
                    af(:,ip2)=af(:,ip2)-fpp
                endif
                ip2=ll(ip2)
                enddo ! particle B
                ip1=ll(ip1)
            enddo ! particle A
            enddo
            enddo
            !$omp endparalleldo

            !print*,'neigh pp'
            ny=(nc)/2
            nx=(nc)/4
            n_iter=ny*nx ! number of pencils
            do ipeny=1,2
            do ipenx=1,4
            !$omp paralleldo default(shared) schedule(dynamic)&
            !$omp& private(i_iter,iy,ix,j,i,ip1,i_neighbor,idl,ip2,rvec,rmag,fpp,tshift)
            do i_iter=0,n_iter-1
                iteam=omp_get_thread_num()+1
                iy=i_iter/nx
                ix=mod(i_iter,nx)
                j=ipeny+2*iy
                i=ipenx+4*ix
                ! print*,'loc',i,j
                if (i>nc .or. j>nc .or. i<1 .or. j<1) then
                    print*,'i,j',i,j
                    print*,'ipenx,ipenx',ipenx,ipenx
                    print*,'ix,iy',ix,iy
                    stop 'particle outside the box'
                endif
                ip1=hoc(i,j)
                do while(ip1/=0) ! particle A
                    ! print*,'ip1',ip1
                    ! Pm1count(iteam) = Pm1count(iteam) + 1
                    ! if (ip1 == 2) then
                    !     ip1=ll(ip1)
                    !     cycle
                    ! endif
                    ! ll2(ip1) = ip1
                    do i_neighbor=1,n_neighbor ! neighbor cells
                        tshift = 0
                        idl=[i,j]+ij(:,i_neighbor)
                        if (minval(idl)<0 .or. maxval(idl)>nc+1) then
                            print*,'ix,iy',ix,iy
                            print*,'i,j',i,j
                            print*,'i_neighbor',i_neighbor,ij(:,i_neighbor)
                            print*,'idl',idl
                            stop  'particle outside the box'
                        endif
                        ip2=hoc(idl(1),idl(2))

                        if (idl(1) < 1  )  tshift(1)  = -ng
                        if (idl(2) < 1  )  tshift(2)  = -ng
                        if (idl(1) > nc  ) tshift(1)  = ng
                        if (idl(2) > nc  ) tshift(2)  = ng
                        do while (ip2/=0) ! particle B
                            rvec=xp(:,ip2)+tshift-xp(:,ip1); rmag=norm2(rvec)
                            if (rmag<apm2) then
                                fpp=rvec/rmag*(F_ra(rmag,app)-F_ra(rmag,apm2))
                                if (any(.not. ieee_is_finite(fpp))) then
                                    print*,'fpp is NaN'
                                    print*,xp(:,ip1),xp(:,ip2)
                                    print*,rvec,rmag
                                    print*,fpp
                                    print*,ip1,ip2
                                    print*,i,j,ipenx,ipeny
                                    print*,itile
                                    print*,iteam
                                    print*,hoc(i,j)
                                    print*,ll(ip1)
                                    stop
                                endif
                                af(:,ip1)=af(:,ip1)+fpp
                                af(:,ip2)=af(:,ip2)-fpp
                            endif
                            ip2=ll(ip2)
                        enddo
                    enddo
                    ip1=ll(ip1)
                enddo
            enddo ! i_iter
            !$omp endparalleldo

            enddo ! ipenx
            enddo ! ipeny

            vp=vp+af*a_mid*dt*G_grid*sim%mass_p_cdm
        call toc(14)
        f2max=maxval(((af(1,:)**2+af(2,:)**2))); vmax = [maxval(abs(vp(1,:))),maxval(abs(vp(2,:)))]
        sim%dt_pp=0.1*sqrt( 1. / (sqrt(f2max)*a_mid*G_grid) )  
        print*,'  vmax =',vmax; print*,'  f2max =',f2max
        print*,'  real time =',tcat(14,istep),'secs'; print*,''
        ! print*,sum(Pm1count)
        ! if (any(ll2==0))then 
        !     print*,'  WARNING: zero-length links'
        !     stop
        ! endif
        deallocate(af)
    endif ! PP

    deallocate(hoc,ll)

    sim%dt_vmax=vbuf*14/maxval(vmax)
    sim%vz_max=vmax(2)
    ! stop

    contains


    subroutine phiential_to_force(force,phi,nps,npe,nfs,nfe)
        use omp_lib
        implicit none
        integer i_0,j_0,i_n(4),j_n(4)
        integer(8) nps,npe,nfs,nfe
        real force(2,nfs:nfe,nfs:nfe)
        real phi(nps:npe,nps:npe)
        ! print*,'phiential_to_force'

        force=0
        
        do j_0=nfs,nfe
        do i_0=nfs,nfe
            ! print*,i_0,j_0
            i_n=i_0+[-2,-1,1,2]
            j_n=j_0+[-2,-1,1,2]
            force(1,i_0,j_0)=sum(phi(i_n,j_0)*weight)
            force(2,i_0,j_0)=sum(phi(i_0,j_n)*weight)
        enddo
        enddo
    endsubroutine
endsubroutine