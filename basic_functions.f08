subroutine tic(i)
  integer i
  call system_clock(tictoc(1,i),t_rate)
endsubroutine

subroutine toc(i)
  integer i
  
  call system_clock(tictoc(2,i),t_rate)
  tcat(i,istep)=real(tictoc(2,i)-tictoc(1,i))/t_rate

endsubroutine

pure function z2str(z)
  character(:),allocatable :: z2str
  character(20) :: str
  real,intent(in) :: z
  write(str,'(f7.3)') z
  z2str=trim(adjustl(str))
endfunction

function output_name(zipname)
  character(*) ::  zipname
  character(:),allocatable :: output_name
  character(20) :: str_z
  write(str_z,'(f7.3)') z_checkpoint(sim%cur_checkpoint)
  output_name=opath//trim(adjustl(str_z))//'_'//zipname//'.bin'
endfunction

function output_name_step(zipname)
  character(*) ::  zipname
  character(:),allocatable :: output_name_step
  character(20) :: str_z
  write(str_z,'(i4)') istep
  output_name_step=opath//'/runtime/'//zipname//'_'//trim(adjustl(str_z))//'.bin'
endfunction

function output_name_ng(zipname)
  character(*) ::  zipname
  character(:),allocatable :: output_name_ng
  character(20) :: str_z
  write(str_z,'(i4)') ng
  output_name_ng=opath//zipname//'_'//trim(adjustl(str_z))//'.bin'
endfunction


! function output_name_halo(zipname)
!   character(*) ::  zipname
!   character(:),allocatable :: output_name_halo
!   character(20) :: str_z
!   write(str_z,'(f7.3)') z_halofind(sim%cur_halofind)
!   output_name_halo=opath//trim(adjustl(str_z))//'_'//zipname//'.bin'
! endfunction

pure function wrap_position(x) result(x_wrapped)
    real, intent(in) :: x
    real :: x_wrapped
    x_wrapped = x - floor(x / rng) * rng
    if (x_wrapped >= ng ) x_wrapped = x_wrapped - rng  
    if (x_wrapped <  0.0) x_wrapped = x_wrapped + rng
end function wrap_position

pure function wrap_position2(x) result(x_wrapped)
    real, intent(in) :: x(2)
    real :: x_wrapped(2)
    integer :: i
    do i=1,2
      x_wrapped(i) = wrap_position(x(i))
    enddo
end function wrap_position2

pure function pbc_vec(x) result(wrapped)
    real, intent(in) :: x(:)
    real :: wrapped(size(x))
    
    wrapped = x - rng * nint(x / rng)
end function pbc_vec


pure function xpos2mesh(x,nmesh) result(x_wrapped)
    real, intent(in) :: x(2)
    integer(8), intent(in) :: nmesh
    integer :: x_wrapped(2)
    integer :: i
    do i=1,2
      x_wrapped(i) = floor(wrap_position(x(i)))
      if (x_wrapped(i) >= nmesh) x_wrapped(i) = x_wrapped(i) - nmesh
      if (x_wrapped(i) <  0.0) x_wrapped(i) = x_wrapped(i) + nmesh
      x_wrapped(i) = x_wrapped(i) + 1
    enddo
end function xpos2mesh

elemental real function F_ra(r,apm)
  real,intent(in) :: r,apm
  real ep
  ep=2*r/apm
  if (apm==0 .or. ep>2) then
    F_ra=r**(-1)
  elseif (ep>=1) then
    F_ra=(1./apm)*(1.539967/ep - 6.823192 + 15.107021*ep - 11.856245*ep**2 + 4.081230*ep**3 - 0.524104*ep**4)
  elseif (ep>=0) then
    F_ra=(1./apm)*(-0.000082 + 2.672707*ep - 0.057116*ep**2 - 1.832992*ep**3 + 0.743081*ep**4)
  else
    F_ra=0
  endif
  ! F_ra = 2 * F_ra
endfunction
