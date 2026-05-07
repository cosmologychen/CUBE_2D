module variables
   use omp_lib
   use parameters
   implicit none


   real stime(istep_max),s2a(istep_max),s2tau(istep_max),s2chi(istep_max),dtau
   integer(8),parameter :: np_image=(nc*np_nc)**2*merge(2,1,body_centered_cubic) ! average number of particles per image
   integer(8),parameter :: np_image_max=np_image*(nte*1./nt)**2
   integer(8),parameter :: np_tile_max=np_image/nnt**2*(nte*1./nt)**2*tile_buffer
   real,parameter :: dt_max=1

   integer itx,ity,ifx,ify
   integer(8) plan,plan2(ncore),iplan,iplan2(ncore) ! FFT plans
   integer(4) ij(2,n_neighbor),iteam,ipm2,ixy2(2,nnt**2)

   real dt,dt_old,dt_mid,dt_e,da,a_mid,f2_max_pp
   real vmax(2),vmax_team(2,ncore),f2max,f2max_team(ncore),overhead_tile,overhead_image,sigma_vi,sigma_vi_new,svz(500,2),svr(100,2)
   real(8) testrho,std_vsim_c,std_vsim_res,std_vsim
   integer ia

   real,allocatable :: Gk1(:,:),Gk2(:,:) ! Green's functions
   real(4),allocatable :: xp(:,:),xp_new(:,:)
   real(4),allocatable :: vp(:,:),vp_new(:,:)
   integer(izipi),allocatable :: pid(:),pid_new(:)!,id_new(:)

   ! one_run_lightcone: 全局格点对应的尺度因子
   ! a_grid(i,j) 表示格点(i,j)实际对应的尺度因子
   ! 初始值来自ic.f90，随模拟演化更新
   real(4),allocatable :: a_grid(:,:),D_grid(:,:)


   real        rho1(nw+2,nw)
   complex     rho1k(nw/2+1,nw)
   equivalence(rho1,rho1k)

   real rho2(1-ngb:ngp+ngb+2,1-ngb:ngp+ngb,ncore)
   complex rho2k(ngt/2+1,ngt,ncore)
   equivalence(rho2,rho2k)

   integer(8),dimension(1-ncb:nt+ncb,nnt,nnt) :: idx_b_l,idx_b_r
   integer(8),dimension(nt,nnt,nnt) :: ppl0,pplr,pprl,ppr0,ppl,ppr

   type type_pm
      integer pm_layer ! 1=coarse, 2=standard, 3=fine
      ! integer iapm
      integer nwork ! grid number to assign density. = nt,ngt,nft(iapm)
      integer nstart(2),nend(2) ! density index range; for cubefft, do FFT padding
      integer tile1(2),tile2(2) ! tiles to work on
      integer nloop ! additional layer of coarse cells to loop over particles
      integer nex ! additional layers of grids for density assignment
      integer nphy ! physical density assignment per tile. =nt,ngp,nfp(iapm)
      real gridsize ! determines gravity strength
      real f2max ! get maximum force squared
      real sigv1,sigv2 ! in/out-come velocity conversion coefficient
      integer tile_shift,utile_shift(2) ! tile and subtile offset
      integer m1,m2 ! density assignment per tile. = [1,nt],[1-ngb,ngp+ngb],[1-nfb(iapm),nfp(iapm)+nfb(iapm)]
      integer nforce ! grid number to update velocity. = nc,ngp,nfp(iapm)
      integer m1phi(2),m2phi(2) ! potential index range. [-2,nc+3],
      integer nc1(2),nc2(2) ! coarse grid iteration
   endtype

   type(type_pm) pm

contains

   subroutine spine_tile(rhoce,idx_ex_r,pp_l,pp_r,ppe_l,ppe_r)
      !! make a particle index (cumulative sumation) on tile
      !! used in update_particle, initial_conditions
      !! input:
      !! rhoce -- particle number density on tile, with 2x buffer depth
      !! output:
      !! idx_ex_r -- last extended index on extended right boundary
      !! ppe_r -- last extended index on physical right boundary
      !! pp_l -- first physical index on physical left boundary
      !! pp_r -- last physical index on physical right boundary
      !! ppe_l -- first extended index on physical left boundary
      implicit none
      integer(4),intent(in) :: rhoce(1-2*ncb:nt+2*ncb,1-2*ncb:nt+2*ncb)
      integer(8),dimension(1-2*ncb:nt+2*ncb),intent(out) :: idx_ex_r
      integer(8),dimension(nt),intent(out) :: pp_l,pp_r,ppe_l,ppe_r
      integer(8) nsum,np_phy
      integer igy
      ! spine := yz-plane to record cumulative sum
      nsum=0
      do igy=1-2*ncb,nt+2*ncb
         nsum=nsum+sum(rhoce(:,igy))
         idx_ex_r(igy)=nsum ! right index
      enddo

      do igy=1,nt
         ppe_r(igy)=idx_ex_r(igy)-sum(rhoce(nt+1:,igy))
      enddo

      nsum=0
      do igy=1,nt
         pp_l(igy)=nsum+1
         np_phy=sum(rhoce(1:nt,igy))
         nsum=nsum+np_phy
         pp_r(igy)=nsum
         ppe_l(igy)=ppe_r(igy)-np_phy+1
      enddo
   endsubroutine

   subroutine spine_image(rhoc,idx_b_l,idx_b_r,ppe_l0,ppe_lr,ppe_rl,ppe_r0,ppl,ppr)
      !! make a particle index (cumulative sumation) on image
      !! used in buffer_grid
      !! input:
      !! rhoc -- particle number density on image, with 1x buffer depth
      !! output:
      !! idx_b_r -- last extended index on extended right boundary
      !! idx_b_l -- zeroth extended index on extended left boundary
      !! ppr     -- last physical index on physical right boundary
      !! ppl     -- zeroth physical index on physical left boundary
      !! ppe_r0  -- last extended index on physical right boundary
      !! ppe_l0  -- zeroth extended index on physical left boundary
      !! ppe_rl  -- last extended index on inner right boundary
      !! ppe_lr  -- zeroth extended index on inner left boundary
      implicit none
      integer(4),intent(in) :: rhoc(1-ncb:nt+ncb,1-ncb:nt+ncb,nnt,nnt)
      integer(8),dimension(1-ncb:nt+ncb,nnt,nnt),intent(out) :: idx_b_l,idx_b_r
      integer(8),dimension(nt,nnt,nnt),intent(out) :: ppe_l0,ppe_lr,ppe_rl,ppe_r0,ppl,ppr
      integer(8) nsum,nsum_p,np_phy,ihy,ihx,igy,ctile_mass(nnt,nnt),ctile_mass_p(nnt,nnt)
      ! spine := yz-plane to record cumulative sum
      nsum=0;nsum_p=0
      do ihy=1,nnt ! sum cumulative tile mass first
         do ihx=1,nnt
            ctile_mass(ihx,ihy)=nsum
            ctile_mass_p(ihx,ihy)=nsum_p
            nsum=nsum+sum(rhoc(:,:,ihx,ihy))
            nsum_p=nsum_p+sum(rhoc(1:nt,1:nt,ihx,ihy))
         enddo
      enddo

      do ihy=1,nnt ! calculate extended spine cumulative index on both sides
         do ihx=1,nnt
            nsum=ctile_mass(ihx,ihy)
            do igy=1-ncb,nt+ncb
               idx_b_l(igy,ihx,ihy)=nsum
               nsum=nsum+sum(rhoc(:,igy,ihx,ihy))
               idx_b_r(igy,ihx,ihy)=nsum
            enddo
         enddo
      enddo

      do ihy=1,nnt ! calculate physical spine
         do ihx=1,nnt
            nsum_p=ctile_mass_p(ihx,ihy)
            do igy=1,nt
               ppl(igy,ihx,ihy)=nsum_p
               np_phy=sum(rhoc(1:nt,igy,ihx,ihy))
               nsum_p=nsum_p+np_phy
               ppr(igy,ihx,ihy)=nsum_p

               ppe_r0(igy,ihx,ihy)=idx_b_r(igy,ihx,ihy)-sum(rhoc(nt+1:,igy,ihx,ihy))
               ppe_rl(igy,ihx,ihy)= ppe_r0(igy,ihx,ihy)-sum(rhoc(nt-ncb+1:nt,igy,ihx,ihy))

               ppe_l0(igy,ihx,ihy)=idx_b_l(igy,ihx,ihy)+sum(rhoc(:0,igy,ihx,ihy))
               ppe_lr(igy,ihx,ihy)=ppe_l0(igy,ihx,ihy) +sum(rhoc(1:ncb,igy,ihx,ihy))
            enddo
         enddo
      enddo

      if (nsum_p  /= sim%np ) then
         print*,'nsum_p,sim%np',nsum_p,sim%np
         stop
      end if
   endsubroutine

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
      print*,term_z,term_r
   endfunction

endmodule
