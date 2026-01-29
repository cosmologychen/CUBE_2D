! #define Cpower
program displacement
    use omp_lib
    use parameters
    use powerspectrum
    implicit none
    save
    include 'fftw3.f'


    integer :: ip,iq(2),i_dim,idx1(2),idx2(2),dx1(2),dx2(2)
    integer(8) np,istat,nthreads,plan,iplan

    integer :: i,j,iteam,cur_checkpoint
    real :: pos0(2),pos1(2),dpos(2),kx(2),pdim(2),xi(10,0:nbin)
    real,allocatable :: rho_grid(:,:,:)
    real,allocatable :: dsp(:,:,:),xp(:,:)
    complex,allocatable :: cdiv(:,:),cphi(:,:),rhok_L(:,:),rhok_R(:,:),rhok_N(:,:)

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


    do cur_checkpoint= n_checkpoint,2,-1
        print*, ''
        print*,'===========================================================',cur_checkpoint
        
#ifdef Cpower
        ! 计算LN的互功率谱
        sim%cur_checkpoint = 1
        open(11,file=output_name('delta_L'),status='old',access='stream')
        read(11) rho1(1:ng,1:ng)
        call sfftw_execute( plan)
        allocate(rhok_L(ng/2+1,ng))
        rhok_L=rho1k/real(ng*ng)
        sim%cur_checkpoint=cur_checkpoint

        rho1(1:ng,1:ng) = DTFE_py('xp')
        call sfftw_execute( plan)
        allocate(rhok_N(ng/2+1,ng))
        rhok_N=rho1k/real(ng*ng)
        call cross_power(xi,rhok_L,rhok_N,np,2)
        open(11,file=output_name('Cpower_LN'),status='replace',access='stream')
        write(11) xi
        close(11)
        deallocate(rhok_N,rhok_L)
#endif

        print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
        sim%cur_checkpoint=cur_checkpoint
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
        ! rho1 = 1
        ! rho1(1:ng,1:ng) = DTFE_py_e('xp',dsp,rho1(1:ng,1:ng))
        ! stop

        do i_dim=1,2
            print*, '   dsp: dimension',int(i_dim,1),'min,max values ='
            print*, '   ',minval(dsp(i_dim,:,:)), maxval(dsp(i_dim,:,:))
        enddo

        print*,'    Write dsp into file:'
        print*,'        ',output_name('dsp')
        open(15,file=output_name('dsp'),status='replace',access='stream')
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
        open(15,file=output_name('delta_E'),status='replace',access='stream')
        write(15) rho1(1:ng,1:ng)
        close(15)
        rho1(1:ng,1:ng) = DTFE_py_e('xp',dsp,rho1(1:ng,1:ng))









    enddo
    print*, 'Finished'

    contains

    function DTFE_py(zipname) result(DTFE) !!! 没有考虑周期性边界条件
        implicit none
        character(*), intent(in) ::  zipname
        character(len=400) :: cmd
        real  rho(0:ng+1,0:ng+1),DTFE(ng,ng)
        real, allocatable :: triangles(:,:,:)
        integer :: unit, filesize, num_reals, ntri

        integer :: t, i, j, xmin_t, xmax_t, ymin_t,ymax_t
        real(8) p(2), a(2), b(2), c(2), v0(2), v1(2), v2(2),area
        real(8) dot00, dot01, dot02, dot11, dot12, invDenom, u, v

        print*, 'Reading DTFE...'

        ! ! 构建调用命令
        ! write(cmd, '(A,A)') "~/anaconda3/envs/camb/bin/python3 DTFE.py ", trim(output_name(zipname))

        ! ! 调用 Python 脚本
        ! print*, "Calling Python script..."
        ! print*, "Command: ", trim(cmd)
        ! call system(trim(cmd))

        inquire(file="triangles.bin", size=filesize)
        num_reals = filesize / 4
        ntri = num_reals / (2 * 3)  ! 每个三角形有 2*3 个坐标

        allocate(triangles(2,3,ntri))
        open(unit=20, file="triangles.bin", form="unformatted", access="stream")
        read(20) triangles
        close(20)

        print*, "ntri: ", ntri
        print*,maxval(triangles),minval(triangles)

        rho = 0

        do t = 1, ntri
            a = triangles(:,1,t)
            b = triangles(:,2,t)
            c = triangles(:,3,t)  
            area = 1/(0.5d0 * abs((b(1) - a(1)) * (c(2) - a(2)) - (c(1) - a(1)) * (b(2) - a(2))))
            if (area > 1d10) cycle  
            ! print*, "t: ", t
            ! print*, "   a: ", a
            ! print*, "   b: ", b
            ! print*, "   c: ", c
            ! print*, "   area: ", area

            xmin_t = floor(  min(a(1), b(1), c(1)))
            xmax_t = ceiling(max(a(1), b(1), c(1)))
            ymin_t = floor(  min(a(2), b(2), c(2)))
            ymax_t = ceiling(max(a(2), b(2), c(2)))
            ! print*, "   xmin, xmax, ymin, ymax: ", xmin_t, xmax_t, ymin_t, ymax_t

            do j = ymin_t, ymax_t
            do i = xmin_t, xmax_t
                p = [i, j]
                v0 = c - a
                v1 = b - a
                v2 = p - a

                dot00 = dot_product(v0, v0)
                dot01 = dot_product(v0, v1)
                dot02 = dot_product(v0, v2)
                dot11 = dot_product(v1, v1)
                dot12 = dot_product(v1, v2)

                invDenom = 1.0d0 / (dot00 * dot11 - dot01 * dot01)
                u = (dot11 * dot02 - dot01 * dot12) * invDenom
                v = (dot00 * dot12 - dot01 * dot02) * invDenom

                if ((u >= 0.d0) .and. (v >= 0.d0) .and. (u + v <= 1.d0)) then
                    rho(i,j) = area
                    ! print*,t,i,j,rho(i,j)
                endif
            enddo
            enddo
        enddo

        ! !$omp paralleldo default(shared) &
        ! !$omp& private(t,a,b,c,area, xmin_t, xmax_t, ymin_t,ymax_t,i, j,v0, v1, v2,dot00, dot01, dot02, dot11, dot12, invDenom, u, v)
        ! do t = 1, ntri
        !     a = triangles(:,1,t)
        !     b = triangles(:,2,t)
        !     c = triangles(:,3,t)
        !     area = 1/(0.5d0 * abs((b(1) - a(1)) * (c(2) - a(2)) - (c(1) - a(1)) * (b(2) - a(2))))
        !     if (area > 1d10) cycle

        !     xmin_t = floor(  min(a(1), b(1), c(1)))
        !     xmax_t = ceiling(max(a(1), b(1), c(1)))
        !     ymin_t = floor(  min(a(2), b(2), c(2)))
        !     ymax_t = ceiling(max(a(2), b(2), c(2)))

        !     do j = ymin_t, ymax_t
        !     do i = xmin_t, xmax_t
        !         p = [i, j]
        !         v0 = c - a
        !         v1 = b - a
        !         v2 = p - a

        !         dot00 = dot_product(v0, v0)
        !         dot01 = dot_product(v0, v1)
        !         dot02 = dot_product(v0, v2)
        !         dot11 = dot_product(v1, v1)
        !         dot12 = dot_product(v1, v2)

        !         invDenom = 1.0d0 / (dot00 * dot11 - dot01 * dot01)
        !         u = (dot11 * dot02 - dot01 * dot12) * invDenom
        !         v = (dot00 * dot12 - dot01 * dot02) * invDenom

        !         if ((u >= 0.d0) .and. (v >= 0.d0) .and. (u + v <= 1.d0)) then
        !             rho(i,j) = area
        !         endif
        !     enddo
        !     enddo
        ! enddo
        ! !$omp endparalleldo


        rho(1 ,:) = rho(1 ,:) + rho(ng+1,:)
        rho(ng,:) = rho(ng,:) + rho(0   ,:)
        rho(: ,1) = rho(: ,1) + rho(:,ng+1)
        rho(:,ng) = rho(:,ng) + rho(:   ,0)
        print*,'    ',minval(rho(1:ng,1:ng)),maxval(rho(1:ng,1:ng)),sum(rho(1:ng,1:ng)*1d0),ng**2

        DTFE(1:ng,1:ng) = rho(1:ng,1:ng)/(sum(rho(1:ng,1:ng)*1d0)/(ng**2))-1
        print*,'    ',minval(DTFE(1:ng,1:ng)),maxval(DTFE(1:ng,1:ng)),sum(DTFE(1:ng,1:ng)*1d0),ng**2
        print*,' DTFE done'
        print*,''
        print*,'    write delta_dD'
        open(15,file=output_name('delta_dD'),status='replace',access='stream')
        write(15) DTFE(1:ng,1:ng)
        close(15)
        deallocate(triangles)
    endfunction


    function DTFE_py_e(zipname,dsp,rho0) result(DTFE) !!! 没有考虑周期性边界条件
        implicit none
        character(*), intent(in) ::  zipname
        character(len=400) :: cmd
        real  rho(0:ng+1,0:ng+1),DTFE(ng,ng),dsp(2,ng,ng),rho0(ng,ng)

        integer ip,ti, iq(2), idl(2)
        real pos(2),dist2,min_dist2

        real, allocatable :: triangles(:,:,:),La(:)
        integer :: unit, filesize, num_reals, ntri

        integer :: t, i, j, xmin_t, xmax_t, ymin_t,ymax_t
        real(8) p(2), a(2,3),rho_E(3), v0(2), v1(2), v2(2),area
        real(8) dot00, dot01, dot02, dot11, dot12, invDenom, u, v
        ! real,allocatable :: xp(:,:)
        integer,allocatable :: ll(:),hoc(:,:),Lc(:),tip(:,:)


        print*, 'Reading DTFE...'
        print*,'    ',minval(rho0(1:ng,1:ng)),maxval(rho0(1:ng,1:ng)),sum(rho0(1:ng,1:ng)*1d0)
        ! stop
        allocate(ll(np),hoc(0:nc+1,0:nc+1))
        
        ll = 0
        hoc = 0
        
        ! print*,np,ratio_cs
        do ip=1,np
            iq(1)=(ip-1)/ng
            iq(2)=modulo(ip-1,int(ng,4))
            pos0=wrap_position2(iq+0.5+dsp(:,iq(1)+1,iq(2)+1))/ratio_cs
            idl=floor(pos0)+1
            ! print*,ip
            ! print*,'  qi :',iq
            ! print*,'  pos:',pos0
            ! print*,'  idl:',idl
            if (maxval(idl) > nc  .or. minval(idl) < 1) then
                print*, 'xp out of range in kick ll',ip
                print*, 'xp(:,ip) =',pos0
                print*, 'fxpos =',floor(pos0)
                print*, 'idl =',idl
                stop 
            endif
            ll(ip)=hoc(idl(1),idl(2))
            hoc(idl(1),idl(2))=ip
        enddo ! ip

        hoc(0   ,:) = hoc(nc,:)
        hoc(nc+1,:) = hoc(1 ,:)
        hoc(: ,0  ) = hoc(:,nc)
        hoc(:,nc+1) = hoc(: ,1)

        print*,'hoc done'
        if (maxval(ll) > ng**2 .or. minval(ll) < 0 .or. maxval(hoc) > ng**2 .or. minval(hoc) < 0) then
            print*,'ll or hoc out of range in kick ll'
            stop 
        endif
        ! stop

        ! ! 构建调用命令
        ! write(cmd, '(A,A)') "~/anaconda3/envs/camb/bin/python3 DTFE.py ", trim(output_name(zipname))

        ! ! 调用 Python 脚本
        ! print*, "Calling Python script..."
        ! print*, "Command: ", trim(cmd)
        ! call system(trim(cmd))

        inquire(file="triangles.bin", size=filesize)
        num_reals = filesize / 4
        ntri = num_reals / (2 * 3)  ! 每个三角形有 2*3 个坐标

        allocate(triangles(2,3,ntri))
        open(unit=20, file="triangles.bin", form="unformatted", access="stream")
        read(20) triangles
        close(20)

        print*, "ntri: ", ntri
        print*,maxval(triangles),minval(triangles)


        allocate(Lc(np),La(np),tip(3,ntri))
        Lc = 0
        tip = 0
        !$omp paralleldo default(shared) &
        !$omp& private(t, a, area, ti, pos0, idl, min_dist2, i, j, ip, iq, dist2)
        do t = 1, ntri
            a = triangles(:,:,t)   
            area = (0.5d0 * abs((a(1,2) - a(1,1)) * (a(2,3) - a(2,1)) - (a(1,3) - a(1,1)) * (a(2,2) - a(2,1))))
            if (area > 1d10) cycle
            ! print*, "t: ", t
            ! print*, "   a: ", a(:,1)
            ! print*, "   b: ", a(:,2)
            ! print*, "   c: ", a(:,3)
            ! print*, "   area: ", area
            do ti = 1, 3
                pos0 = wrap_position2(real(a(:,ti),4))/ratio_cs
                idl=floor(pos0)+1

                min_dist2 = 0.01  ! 初始最小距离设为一个很大的数

                do i = -1,1
                do j = -1,1
                    ip = hoc(idl(1)+i, idl(2)+j)
                    do while (ip /= 0)
                        iq(1)=(ip-1)/ng+1
                        iq(2)=modulo(ip-1,int(ng,4))+1
                        pos0=wrap_position2(iq-0.5+dsp(:,iq(1),iq(2)))
                        dist2 = sum((pos0 - a(:,ti))**2)

                        if (dist2 < min_dist2) then
                            min_dist2 = dist2
                            tip(ti,t) = ip
                        endif

                        ip = ll(ip)
                    enddo
                enddo
                enddo

                if (min_dist2 > 0.001) then
                    print*, '       WARNING: min_dist2 > 0.01'
                    print*, ' ti',ti,t
                    print*, ' a',a(:,ti)
                    print*, ' wa',wrap_position2(real(a(:,ti),4))/ratio_cs
                    print*, '   idl:',idl
                    print*, '       min_dist2',min_dist2

                    min_dist2 = 16  ! 初始最小距离设为一个很大的数
                    rho_E(ti) = -1e6     ! 默认密度值

                    do i = -1,1
                    do j = -1,1
                        ip = hoc(idl(1), idl(2))
                        print*,'============================='
                        print*, ' a',a(:,ti)
                        print*, '   idl:',idl
                        print*, '   ip',ip
                        do while (ip /= 0)
                            iq(1)=(ip-1)/ng
                            iq(2)=modulo(ip-1,int(ng,4))
                            pos0=wrap_position2(iq+0.5+dsp(:,iq(1)+1,iq(2)+1))
                            dist2 = sum((pos0 - a(:,ti))**2)
                            print*, '       '
                            print*, '       ip',ip
                            print*, '       iq',iq
                            print*, '       pos0',pos0
                            print*, '       dis ',dist2

                            if (dist2 < min_dist2) then
                                min_dist2 = dist2
                                tip(i,t) = ip
                            endif

                            ip = ll(ip)
                        enddo
                    enddo
                    enddo

                    stop
                endif
                Lc(tip(ti,t)) = Lc(tip(ti,t)) + 1
                La(tip(ti,t)) = La(tip(ti,t)) + area
            enddo
        enddo
        !$omp endparalleldo

        rho = 0


        ! do t = 1, ntri
        !     a = triangles(:,:,t)   
        !     area = 1/(0.5d0 * abs((a(1,2) - a(1,1)) * (a(2,3) - a(2,1)) - (a(1,3) - a(1,1)) * (a(2,2) - a(2,1))))
        !     if (area > 1d10) cycle
        !     ! print*, "t: ", t
        !     ! print*, "   a: ", a(:,1)
        !     ! print*, "   b: ", a(:,2)
        !     ! print*, "   c: ", a(:,3)
        !     ! print*, "   area: ", area

        !     do ti = 1, 3
        !         pos0 = wrap_position2(real(a(:,ti),4))/ratio_cs
        !         idl=floor(pos0)+1

        !         min_dist2 = 16  ! 初始最小距离设为一个很大的数
        !         rho_E(ti) = -1e6     ! 默认密度值

        !         do i = -1,1
        !         do j = -1,1
        !             ip = hoc(idl(1)+i, idl(2)+j)
        !             do while (ip /= 0)
        !                 iq(1)=(ip-1)/ng
        !                 iq(2)=modulo(ip-1,int(ng,4))
        !                 pos0=wrap_position2(iq+0.5+dsp(:,iq(1)+1,iq(2)+1))
        !                 dist2 = sum((pos0 - a(:,ti))**2)

        !                 if (dist2 < min_dist2) then
        !                     min_dist2 = dist2
        !                     rho_E(ti) = rho0(iq(1)+1, iq(2)+1)
        !                 endif

        !                 ip = ll(ip)
        !             enddo
        !         enddo
        !         enddo
        !         ! print*, '       min_dist2',min_dist2
        !         ! print*, '       rho_E(ti)',rho_E(ti)
        !         ! print*,'++++++++++++++++++++++++++++++++++'
        !         ! print*,''
        !         if (min_dist2 > 0.01) then
        !             print*, '       WARNING: min_dist2 > 0.01'
        !             print*, ' ti',ti,t
        !             print*, ' a',a(:,ti)
        !             print*, ' wa',wrap_position2(real(a(:,ti),4))/ratio_cs
        !             print*, '   idl:',idl
        !             print*, '       min_dist2',min_dist2
        !             print*, '       rho_E(ti)',rho_E(ti)

        !             min_dist2 = 16  ! 初始最小距离设为一个很大的数
        !             rho_E(ti) = -1e6     ! 默认密度值

        !             do i = -1,1
        !             do j = -1,1
        !                 ip = hoc(idl(1), idl(2))
        !                 print*,'============================='
        !                 print*, ' a',a(:,ti)
        !                 print*, '   idl:',idl
        !                 print*, '   ip',ip
        !                 do while (ip /= 0)
        !                     iq(1)=(ip-1)/ng
        !                     iq(2)=modulo(ip-1,int(ng,4))
        !                     pos0=wrap_position2(iq+0.5+dsp(:,iq(1)+1,iq(2)+1))
        !                     dist2 = sum((pos0 - a(:,ti))**2)
        !                     print*, '       '
        !                     print*, '       ip',ip
        !                     print*, '       iq',iq
        !                     print*, '       pos0',pos0
        !                     print*, '       dis ',dist2

        !                     if (dist2 < min_dist2) then
        !                         min_dist2 = dist2
        !                         rho_E(ti) = rho0(iq(1)+1, iq(2)+1)
        !                     endif

        !                     ip = ll(ip)
        !                 enddo
        !             enddo
        !             enddo

        !             stop
        !         endif
        !     enddo
        !     area = sum(rho_E)/3/(0.5d0 * abs((a(1,2) - a(1,1)) * (a(2,3) - a(2,1)) - (a(1,3) - a(1,1)) * (a(2,2) - a(2,1))))

        !     xmin_t = floor(  minval(a(1,:)))
        !     xmax_t = ceiling(maxval(a(1,:)))
        !     ymin_t = floor(  minval(a(2,:)))
        !     ymax_t = ceiling(maxval(a(2,:)))
        !     ! print*, "   xmin, xmax, ymin, ymax: ", xmin_t, xmax_t, ymin_t, ymax_t

        !     do j = ymin_t, ymax_t
        !     do i = xmin_t, xmax_t
        !         v0 = a(:,3) - a(:,1)
        !         v1 = a(:,2) - a(:,1)
        !         v2 = [i, j] - a(:,1)

        !         dot00 = dot_product(v0, v0)
        !         dot01 = dot_product(v0, v1)
        !         dot02 = dot_product(v0, v2)
        !         dot11 = dot_product(v1, v1)
        !         dot12 = dot_product(v1, v2)

        !         invDenom = 1.0d0 / (dot00 * dot11 - dot01 * dot01)
        !         u = (dot11 * dot02 - dot01 * dot12) * invDenom
        !         v = (dot00 * dot12 - dot01 * dot02) * invDenom

        !         if ((u >= 0.d0) .and. (v >= 0.d0) .and. (u + v <= 1.d0)) then
        !             ! rho(i,j) = area*(u*rho_E(1)+v*rho_E(2)+(1.d0-u-v)*rho_E(3))
        !             rho(i,j) = area
        !             ! print*,t,i,j,rho(i,j)
        !         endif
        !     enddo
        !     enddo
        ! enddo

        !$omp paralleldo default(shared) &
        !$omp& private(t, a, area, rho_E, ti, idl, min_dist2, ip, iq, pos0, dist2, xmin_t, xmax_t, ymin_t,ymax_t,i, j,v0, v1, v2,dot00, dot01, dot02, dot11, dot12, invDenom, u, v)
        do t = 1, ntri
            a = triangles(:,:,t)   
            area = 1/(0.5d0 * abs((a(1,2) - a(1,1)) * (a(2,3) - a(2,1)) - (a(1,3) - a(1,1)) * (a(2,2) - a(2,1))))
            if (area > 1d10) cycle
            ! print*, "t: ", t
            ! print*, "   a: ", a(:,1)
            ! print*, "   b: ", a(:,2)
            ! print*, "   c: ", a(:,3)
            ! print*, "   area: ", area

            ! rho_E = 1
            ! do ti = 1, 3
            !     pos0 = wrap_position2(real(a(:,ti),4))/ratio_cs
            !     idl=floor(pos0)+1

            !     min_dist2 = 0.01  ! 初始最小距离设为一个很大的数
            !     rho_E(ti) = -1e6     ! 默认密度值

            !     do i = -1,1
            !     do j = -1,1
            !         ip = hoc(idl(1)+i, idl(2)+j)
            !         do while (ip /= 0)
            !             iq(1)=(ip-1)/ng+1
            !             iq(2)=modulo(ip-1,int(ng,4))+1
            !             pos0=wrap_position2(iq-0.5+dsp(:,iq(1),iq(2)))
            !             dist2 = sum((pos0 - a(:,ti))**2)

            !             if (dist2 < min_dist2) then
            !                 min_dist2 = dist2
            !                 rho_E(ti) = rho0(iq(1), iq(2))
            !             endif

            !             ip = ll(ip)
            !         enddo
            !     enddo
            !     enddo
            !     ! print*, '       min_dist2',min_dist2
            !     ! print*, '       rho_E(ti)',rho_E(ti)
            !     ! print*,'++++++++++++++++++++++++++++++++++'
            !     ! print*,''
            !     if ((min_dist2 > 0.001) .or. (rho_E(ti) == -1e6)) then
            !         if (min_dist2 > 0.001) print*, '       WARNING: min_dist2 > 0.01'
            !         if (rho_E(ti) == -1e6) print*, '       WARNING: rho_E(ti) == -1e6'
            !         print*, ' ti',ti,t
            !         print*, ' a',a(:,ti)
            !         print*, ' wa',wrap_position2(real(a(:,ti),4))/ratio_cs
            !         print*, '   idl:',idl
            !         print*, '       min_dist2',min_dist2
            !         print*, '       rho_E(ti)',rho_E(ti)

            !         min_dist2 = 16  ! 初始最小距离设为一个很大的数
            !         rho_E(ti) = -1e6     ! 默认密度值

            !         do i = -1,1
            !         do j = -1,1
            !             ip = hoc(idl(1), idl(2))
            !             print*,'============================='
            !             print*, ' a',a(:,ti)
            !             print*, '   idl:',idl
            !             print*, '   ip',ip
            !             do while (ip /= 0)
            !                 iq(1)=(ip-1)/ng
            !                 iq(2)=modulo(ip-1,int(ng,4))
            !                 pos0=wrap_position2(iq+0.5+dsp(:,iq(1)+1,iq(2)+1))
            !                 dist2 = sum((pos0 - a(:,ti))**2)
            !                 print*, '       '
            !                 print*, '       ip',ip
            !                 print*, '       iq',iq
            !                 print*, '       pos0',pos0
            !                 print*, '       dis ',dist2

            !                 if (dist2 < min_dist2) then
            !                     min_dist2 = dist2
            !                     rho_E(ti) = rho0(iq(1)+1, iq(2)+1)
            !                 endif

            !                 ip = ll(ip)
            !             enddo
            !         enddo
            !         enddo

            !         stop
            !     endif
            ! enddo
            do ti = 1, 3
                ip = tip(ti,t)
                iq(1)=(ip-1)/ng+1
                iq(2)=modulo(ip-1,int(ng,4))+1
                rho_E(ti) = rho0(iq(1), iq(2))
            enddo
            area = sum(rho_E)/3/(0.5d0 * abs((a(1,2) - a(1,1)) * (a(2,3) - a(2,1)) - (a(1,3) - a(1,1)) * (a(2,2) - a(2,1))))

            xmin_t = floor(  minval(a(1,:)))
            xmax_t = ceiling(maxval(a(1,:)))
            ymin_t = floor(  minval(a(2,:)))
            ymax_t = ceiling(maxval(a(2,:)))
            ! print*, "   xmin, xmax, ymin, ymax: ", xmin_t, xmax_t, ymin_t, ymax_t

            do j = ymin_t, ymax_t
            do i = xmin_t, xmax_t
                v0 = a(:,3) - a(:,1)
                v1 = a(:,2) - a(:,1)
                v2 = [i, j] - a(:,1)

                dot00 = dot_product(v0, v0)
                dot01 = dot_product(v0, v1)
                dot02 = dot_product(v0, v2)
                dot11 = dot_product(v1, v1)
                dot12 = dot_product(v1, v2)

                invDenom = 1.0d0 / (dot00 * dot11 - dot01 * dot01)
                u = (dot11 * dot02 - dot01 * dot12) * invDenom
                v = (dot00 * dot12 - dot01 * dot02) * invDenom

                if ((u >= 0.d0) .and. (v >= 0.d0) .and. (u + v <= 1.d0)) then
                    ! rho(i,j) = area*(u*rho_E(1)+v*rho_E(2)+(1.d0-u-v)*rho_E(3))
                    rho(i,j) = area
                    ! print*,t,i,j,rho(i,j)
                endif
            enddo
            enddo
        enddo
        !$omp endparalleldo


        rho(1 ,:) = rho(1 ,:) + rho(ng+1,:)
        rho(ng,:) = rho(ng,:) + rho(0   ,:)
        rho(: ,1) = rho(: ,1) + rho(:,ng+1)
        rho(:,ng) = rho(:,ng) + rho(:   ,0)
        print*,'    ',minval(rho(1:ng,1:ng)),maxval(rho(1:ng,1:ng)),sum(rho(1:ng,1:ng)*1d0),ng**2

        DTFE(1:ng,1:ng) = rho(1:ng,1:ng)/(sum(rho(1:ng,1:ng)*1d0)/(ng**2))-1
        print*,'    ',minval(DTFE(1:ng,1:ng)),maxval(DTFE(1:ng,1:ng)),sum(DTFE(1:ng,1:ng)*1d0),ng**2
        print*,' DTFE done'
        print*,''
        print*,'    write delta_dDE'
        open(15,file=output_name('delta_dDE'),status='replace',access='stream')
        write(15) DTFE(1:ng,1:ng)
        close(15)
        ! deallocate(ll,hoc,triangles)
    endfunction


end