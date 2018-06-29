real function userff (npar, data, myid)

  ! Uses

  use em_lib
  use ISO_FORTRAN_ENV

  ! No implicit typing

  implicit none

  ! Variables

  integer :: n
  integer :: ns=0
  integer :: ll
  integer :: l0=0
  integer :: l1=0
  integer :: l2=0
  integer :: l3=0
  integer :: id
  integer :: i
  integer :: j
  integer :: t_code

  integer, intent(in)   :: npar
  integer, intent(in)   :: myid
  integer, dimension(5) :: spec
  integer, dimension(:), allocatable :: nl0
  integer, dimension(:), allocatable :: nl1
  integer, dimension(:), allocatable :: nl2
  integer, dimension(:), allocatable :: nl3

  character(len=1) :: lchr

  real(dp) :: f
  real(dp) :: e
  real(dp) :: Teff_o
  real(dp) :: Teff_e
  real(dp) :: logg_o
  real(dp) :: logg_e
  real(dp) :: FeH_o
  real(dp) :: FeH_e
  real(dp) :: R_o
  real(dp) :: R_e
  real(dp) :: L_o
  real(dp) :: L_e
  real(dp) :: Teff
  real(dp) :: logg
  real(dp) :: FeH
  real(dp) :: R
  real(dp) :: L
  real(dp) :: M_mod
  real(dp) :: Z_mod
  real(dp) :: Y_mod
  real(dp) :: alpha_mod
  real(dp) :: age
  real(dp) :: resid
  real(dp) :: chisq_r=0.d0

  real(dp), dimension(npar), intent(in) :: data
  real(dp), dimension(:), allocatable :: nu0
  real(dp), dimension(:), allocatable :: dnu0
  real(dp), dimension(:), allocatable :: Enorm0
  real(dp), dimension(:), allocatable :: nu1
  real(dp), dimension(:), allocatable :: dnu1
  real(dp), dimension(:), allocatable :: Enorm1
  real(dp), dimension(:), allocatable :: nu2
  real(dp), dimension(:), allocatable :: dnu2
  real(dp), dimension(:), allocatable :: Enorm2
  real(dp), dimension(:), allocatable :: nu3
  real(dp), dimension(:), allocatable :: dnu3
  real(dp), dimension(:), allocatable :: Enorm3

  type(freq_t) :: fr_obs(0:3)
  type(freq_t) :: fr_mod(0:3)
  type(freq_t) :: fr_cor(0:3)

  ! First pass to set array sizes

  open(unit=100, file='obs.dat', status='old', action='read')
  read(100,*) n
  do i=1,n
     read(100,*) ll
     if(ll.eq.0) l0=l0+1
     if(ll.eq.1) l1=l1+1
     if(ll.eq.2) l2=l2+1
     if(ll.eq.3) l3=l3+1
  enddo
  ns=0
  do i=1,5
     read(100,*,end=50) lchr
     ns=ns+1
  enddo
50 close(100)

  allocate ( nl0(l0) )
  allocate ( nl1(l1) )
  allocate ( nl2(l2) )
  allocate ( nl3(l3) )
  allocate ( nu0(l0) )
  allocate ( nu1(l1) )
  allocate ( nu2(l2) )
  allocate ( nu3(l3) )
  allocate ( dnu0(l0) )
  allocate ( dnu1(l1) )
  allocate ( dnu2(l2) )
  allocate ( dnu3(l3) )
  allocate ( Enorm0(l0) )
  allocate ( Enorm1(l1) )
  allocate ( Enorm2(l2) )
  allocate ( Enorm3(l3) )

  ! Read observational constraints

  open(unit=100, file='obs.dat', status='old', action='read')
  read(100,*) n
  l0=0
  l1=0
  l2=0
  l3=0
  spec(1)=0
  spec(2)=0
  spec(3)=0
  spec(4)=0
  spec(5)=0
  do i=1,n+ns
     ll=-1
     read(100,*) lchr,f,e
     if (lchr.eq."T") then
        Teff_o=f
        Teff_e=e
        spec(1)=1
     elseif (lchr.eq."G") then
        logg_o=f
        logg_e=e
        spec(2)=1
     elseif (lchr.eq."M") then
        FeH_o=f
        FeH_e=e
        spec(3)=1
     elseif (lchr.eq."R") then
        R_o=f
        R_e=e
        spec(4)=1
     elseif (lchr.eq."L") then
        L_o=f
        L_e=e
        spec(5)=1
     else
        read(lchr,*) ll
        if(ll.eq.0) then
           l0=l0+1
           nl0(l0)=l0
           nu0(l0)=f
           dnu0(l0)=e
           Enorm0(l0)=1.d0
        elseif(ll.eq.1) then
           l1=l1+1
           nl1(l1)=l1
           nu1(l1)=f
           dnu1(l1)=e
           Enorm1(l1)=1.d0
        elseif(ll.eq.2) then
           l2=l2+1
           nl2(l2)=l2
           nu2(l2)=f
           dnu2(l2)=e
           Enorm2(l2)=1.d0
        elseif(ll.eq.3) then
           l3=l3+1
           nl3(l3)=l3
           nu3(l3)=f
           dnu3(l3)=e
           Enorm3(l3)=1.d0
        endif
     endif
  enddo
  close(100)

  fr_obs = freq_t()
  fr_obs(0) = freq_t(nu0,dnu0,Enorm0,nl0)
  fr_obs(1) = freq_t(nu1,dnu1,Enorm1,nl1)
  fr_obs(2) = freq_t(nu2,dnu2,Enorm2,nl2)
  fr_obs(3) = freq_t(nu3,dnu3,Enorm3,nl3)

  call set_obs_freqs(0, fr_obs(0))
  call set_obs_freqs(1, fr_obs(1))
  call set_obs_freqs(2, fr_obs(2))
  call set_obs_freqs(3, fr_obs(3))

  ! Create a star

  M_mod = data(1)+0.75
  Z_mod = 10.**(1.4*data(2)-2.7)
  Y_mod = 0.10*data(3)+0.22
  alpha_mod = 2.0*data(4)+1.0

  id = create_star( &
       M = M_mod, &
       Z = Z_mod, &
       Y = Y_mod, &
       alpha = alpha_mod, &
       f_ov=0.0d0, &
       max_age=1D11)

  ! Evolve it to the ZAMS

  call evolve_star_to_zams(id, t_code)

  if (t_code == t_ok) then
     print *,'Evolve to ZAMS: OK'
  else
     print *,'Evolve to ZAMS: Failed, termination code =', t_code
     stop
  end if

  ! Evolve it until seismic constraints are met

  print *,' calling evolve_star_seismic'
  call evolve_star_seismic(id, t_code)
  print *,' done'

  if (t_code == t_ok) then
     print *,'Evolve to seismic: OK'
  elseif (t_code == t_max_age) then
     print *,'Evolve to seismic: Reached maximum age'
  else
     print *,'Evolve to seismic: Failed, termination code =', t_code
     stop
  endif

  ! Analyze results if the seismic run went OK

  if (t_code == t_ok) then

     ! Get model data

     call get_mod_data(Teff, logg, FeH, R, L, age)

     ! Calculate spectroscopic chisq

     if (spec(1).eq.1) chisq_r = chisq_r + (Teff_o-Teff)*(Teff_o-Teff)/(Teff_e*Teff_e)
     if (spec(2).eq.1) chisq_r = chisq_r + (logg_o-logg)*(logg_o-logg)/(logg_e*logg_e)
     if (spec(3).eq.1) chisq_r = chisq_r + (FeH_o-FeH)*(FeH_o-FeH)/(FeH_e*FeH_e)
     if (spec(4).eq.1) chisq_r = chisq_r + (R_o-R)*(R_o-R)/(R_e*R_e)
     if (spec(5).eq.1) chisq_r = chisq_r + (L_o-L)*(L_o-L)/(L_e*L_e)

     ! Get model frequencies

     fr_mod = freq_t()
     fr_mod(0) = get_mod_freqs(0)
     fr_mod(1) = get_mod_freqs(1)
     fr_mod(2) = get_mod_freqs(2)
     fr_mod(3) = get_mod_freqs(3)

     ! Get corrected frequencies

     fr_cor = freq_t()
     call apply_combined_correction(fr_mod, fr_obs, fr_cor)

     ! Calculate seismic chisq

     do i = 0, 3
        do j = 1, fr_obs(i)%n
           resid = (fr_obs(i)%nu(j)-fr_cor(i)%nu(j))/fr_obs(i)%dnu(j)
           chisq_r = chisq_r + resid*resid
        end do
     end do
  
     userff = float(n+ns-5)/chisq_r

  end if

  ! Finish

end function userff

