!
! Copyright (C) 2017
! Cyrille Barreteau <mailto:cyrille.barreteau@cea.fr>,
! Mathieu Cesar <mailto:mathieu.cesar@cea.fr>,
! Pascal Thibaudeau <mailto:pascal.thibaudeau@cea.fr>.
!
! This software is a computer program whose purpose is TBKOSTER.
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and inRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.
!
! atom.f90
! TBKOSTER
module atom_mod
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use element_mod
   use lattice_mod
   use math_mod, only: deg2rad, epsilon, i_unit, pi, rad2deg, sph2cart
   use precision_mod, only: rp
   use string_mod, only: sl, int2str, lower, real2str, &
                         fixedreal2str, unique_str, TBKOSTER_flush
   use units_mod
   use constant_mod, only: a_0
   implicit none
   private

   !> Derived type properties for i/o methods
   character(len=sl), dimension(*), parameter :: property_list = &
      [character(len=sl) :: 'ns', 'na', 'ntag', 'stag', 'tag', &
      'ia2ie', 'pbc', 'k_spiral', 'r_coord', 'm_coord', 'r', &
      'p', 'm', 'lambda_pen']

   type, public :: atom
      !> Units
      class(units), pointer :: u
      !> Elements
      class(element), pointer :: e
      !> Real space lattice
      class(lattice), pointer :: l_r
      !> Reciprocal space lattice
      class(lattice), pointer :: l_k

      !> @defgroup Spin Spin-related variables
      !> @{
      !> Spin polarization ; options:
      !>  'unpolarized' : spin-unpolarized (default)
      !>  'collinear    : spin-polarized collinear
      !>  'noncollinear': spin-polarized noncollinear
      !character(len=12) :: spin_polarization
      !> Number of spins, options:
      !>  1: spin-unpolarized
      !>  2: spin-polarized collinear
      !>  4: spin-polarized noncollinear
      integer :: ns
      !> Number of spin loops (ns=1,4->1, ns=2->2)
      integer :: nsl
      !> Spin polarization (ns=1->1, ns=2,4->2)
      integer :: nsp
      !> Spin degeneracy (ns=1->2, ns=2,4->1)
      integer :: g_s
      !> 2spin/2spin-to-4spin index
      integer, dimension(2, 2) :: iss2is
      !> @}

      !> @defgroup Atom Atom-related variables
      !> @{
      !> Number of atoms
      integer :: na
      !> Number of tags
      integer :: ntag
      !> Strides between tags
      integer, dimension(:), allocatable :: stag
      !> Tags
      character(len=sl), dimension(:), allocatable :: tag
      !> Atom-to-element index
      integer, dimension(:), allocatable :: ia2ie
      !> Number of electrons
      real(rp) :: nel
      !> @}

      !> Periodic boundary conditions
      integer, dimension(3) :: pbc
      !> Spin spiral vector
      real(rp), dimension(3) :: k_spiral

      !> @defgroup Dynamics Dynamics-related variables
      !> @{
      !> Position coordinates; options:
      !>  'cartesian'
      !>  'direct' (default)
      character(len=9) :: r_coord
      !> Positions, (na,3)
      real(rp), dimension(:, :), allocatable :: r
      !> Momentums, (na,3)
      real(rp), dimension(:, :), allocatable :: p
      !> Expression of the magnetization; options:
      !> 'cartesian'
      !> 'spherical' (default)
      character(len=9) :: m_coord
      !> Magnetizations in spherical coordinates (r,theta,phi) (radial, polar and
      !> azimuthal component), (na,3)
      real(rp), dimension(:, :), allocatable :: m, m_pen
      !> @}

      !> @defgroup Neighbours Neighbours-related variables
      !> @{

      !> Number of neighbours
      integer, dimension(:), allocatable :: nn
      !> Maximum number of neighbours
      integer :: nn_max
      !> Atom/neighbour-to-atom index
      integer, dimension(:, :), allocatable :: ian2ia
      !> Atom/periodic-boundary-to-neighbour index
      integer, dimension(:, :, :, :, :), allocatable :: iapbc2in
      !> Neighbour position
      real(rp), dimension(:, :, :), allocatable :: rn
      !> @}

      !> @defgroup Magnetic_penalization Magnetic penalization related variables
      !> @{

      !> Magnetic penalization.\n
      !> See \cite Gebauer1999 page 45-48.\n
      !> To impose a magnetic constraint on a system, an extra term is added to the
      !> total energy functional \f$ E_{\mathrm{tot}}[N,\mathbf{M}] \f$ in the form
      !> of a penalization energy functional.\n
      !> The penalization is chosen such that it is zero if the constraint is
      !> fulfilled and large and positive otherwise. In this way, the constraint is
      !> imposed approximately by energetically penalizing any deviation from a
      !> certain condition.\n
      !> When trying to impose a given magnetization \f$ \mathbf{M}_{\mathrm{pen}}
      !> \f$, the penalization can be chosen to be proportional to the mean square
      !> deviation of the actual magnetization \f$ \mathbf{M} \f$ from its target
      !> \f$ \mathbf{M}_{\mathrm{pen}} \f$.\n
      !> The new energy functional reads:
      !> \f{equation}{
      !>   E_{\lambda}[N,\mathbf{M}] = E_{\mathrm{tot}}[N,\mathbf{M}] + \lambda \int
      !>   \mathrm{d}\mathbf{r} (\mathbf{M}(\mathbf{r}) - \mathbf{M}_{\mathrm{pen}}
      !>   (\mathbf{r}))^2
      !> \f}
      !> where \f$ \lambda \f$ is a large positive multiplier.\n
      !> For a finite \f$ \lambda \f$, the constraint is not imposed exactly as the
      !> difference \f$ \mathbf{M}(\mathbf{r}) - \mathbf{M}_{\mathrm{pen}}
      !> (\mathbf{r}) \f$ remains finite too.\n
      !> In the limit \f$ \lambda \rightarrow \infty \f$, the difference \f$
      !> \mathbf{M}(\mathbf{r}) - \mathbf{M}_{\mathrm{pen}} (\mathbf{r}) \f$ tends
      !> to zero. This however comes at the cost of a slower convergence.\n
      !> In practical calculations, \f$ \lambda \f$ must therefore be chosen
      !> sufficiently large such that \f$ \mathbf{M} \f$ is close enough to \f$
      !> \mathbf{M}_{\mathrm{pen}} \f$, while still being as small as possible in
      !> order to allow fast convergence of the self-consistent cycle.
      !> Multiplier for the magnetic penalization
      real(rp), dimension(:), allocatable   :: lambda_pen
      !> Effective magnetic induction from the magnetic penalization, 'na,3)
      real(rp), dimension(:, :), allocatable :: b_pen
      !> @}

   contains
      final :: destructor
      procedure :: build_c_k
      procedure :: build_dk_c_k
      procedure :: calculate_neighbours
      procedure :: calculate_nel
      procedure :: calculate_spin
      procedure :: read_txt
      procedure :: write_lammps
      procedure :: write_txt
      procedure :: write_txt_formatted
      procedure :: write_xyz
      procedure :: read_namelist => read_txt
      procedure :: write_namelist_formatted => write_txt_formatted
   end type atom

   ! Constructor
   interface atom
      procedure :: constructor
   end interface atom

contains
   function constructor(e, l_r, l_k) result(obj)
      class(element), target, intent(in) :: e
      class(lattice), target, intent(in) :: l_r, l_k
      type(atom) :: obj
      obj%u => e%u
      obj%e => e
      obj%l_r => l_r
      obj%l_k => l_k
      obj%iss2is = reshape((/1, 4, 3, 2/), (/2, 2/))
   end function constructor

   subroutine destructor(obj)
      type(atom) :: obj
      if (allocated(obj%stag)) deallocate (obj%stag)
      if (allocated(obj%tag)) deallocate (obj%tag)
      if (allocated(obj%ia2ie)) deallocate (obj%ia2ie)
      if (allocated(obj%r)) deallocate (obj%r)
      if (allocated(obj%p)) deallocate (obj%p)
      if (allocated(obj%m)) deallocate (obj%m)
      if (allocated(obj%m_pen)) deallocate (obj%m_pen)
      if (allocated(obj%nn)) deallocate (obj%nn)
      if (allocated(obj%ian2ia)) deallocate (obj%ian2ia)
      if (allocated(obj%iapbc2in)) deallocate (obj%iapbc2in)
      if (allocated(obj%rn)) deallocate (obj%rn)
      if (allocated(obj%lambda_pen)) deallocate (obj%lambda_pen)
      if (allocated(obj%b_pen)) deallocate (obj%b_pen)
   end subroutine destructor

   function build_ia2ie(tag_in) result(ia2ie)
      character(len=*), dimension(:), intent(in) :: tag_in
      integer, dimension(size(tag_in)) :: ia2ie
      character(len=len(tag_in)), dimension(:), allocatable :: cdummy
      integer, dimension(:), allocatable :: idummy
      character(len=2), dimension(:), allocatable :: tug
      allocate (tug(size(tag_in)))
      tug = tag_in(:)(1:2)
      call unique_str(tug, cdummy, idummy, ia2ie)
      deallocate (tug)
   end function build_ia2ie

   !> calculate the factor e^{i k.(R+tau_j-tau_i)
   function build_c_k(obj, k) result(c_k)
      class(atom), intent(in) :: obj
      real(rp), dimension(3), intent(in) :: k
      complex(rp), dimension(obj%na, obj%nn_max, obj%nsp) :: c_k
      real(rp) :: k_t(3), RR_c(3), RR_u(3), P_R
      integer :: ia, in_idx, is
      c_k = cmplx(0.0_rp, 0.0_rp, kind=rp)
      do ia = 1, obj%na
         do in_idx = 1, obj%nn(ia)
            RR_c = obj%rn(ia, in_idx, :)
            RR_u = matmul(obj%l_k%v, RR_c)
            if (obj%ns <= 2) then
               P_R = dot_product(k, RR_u)
               c_k(ia, in_idx, 1:obj%ns) = exp(2*pi*i_unit*P_R)
            else
               do is = 1, 2
                  k_t = k + (2*is - 3)*0.5_rp*obj%k_spiral
                  P_R = dot_product(k_t, RR_u)
                  c_k(ia, in_idx, is) = exp(2*pi*i_unit*P_R)
               end do
            end if
         end do
      end do
   end function build_c_k

   !> calculate the derivative of the factor e^{i k.(R+tau_j-tau_i)
   function build_dk_c_k(obj, k) result(dk_c_k)
      class(atom), intent(in) :: obj
      real(rp), dimension(3), intent(in) :: k
      complex(rp), dimension(3, obj%na, obj%nn_max, obj%nsp) :: dk_c_k
      real(rp) :: k_t(3), RR_c(3), RR_u(3), P_R
      integer :: ia, in_idx, is
      dk_c_k = cmplx(0.0_rp, 0.0_rp, kind=rp)
      do ia = 1, obj%na
         do in_idx = 1, obj%nn(ia)
            RR_c = obj%rn(ia, in_idx, :)
            RR_u = matmul(obj%l_k%v, RR_c)
            if (obj%ns <= 2) then
               P_R = dot_product(k, RR_u)
               dk_c_k(:, ia, in_idx, 1:obj%ns) = &
                  spread(i_unit*RR_c*exp(2*pi*i_unit*P_R), 2, obj%ns)
            else
               do is = 1, 2
                  k_t = k + (2*is - 3)*0.5_rp*obj%k_spiral
                  P_R = dot_product(k_t, RR_u)
                  dk_c_k(:, ia, in_idx, is) = &
                     i_unit*RR_c*exp(2*pi*i_unit*P_R)
               end do
            end if
         end do
      end do
   end function build_dk_c_k

   subroutine calculate_neighbours(obj, r_max, type_in)
      class(atom), intent(inout) :: obj
      real(rp), intent(inout) :: r_max
      character(len=3) :: type_in
      integer, dimension(:), allocatable :: nn_l
      integer, dimension(:, :), allocatable :: ian_l
      integer, dimension(:, :, :, :, :), allocatable :: iapbc_l
      real(rp), dimension(:, :, :), allocatable :: rn_l
      integer :: ia1, ia2, in_l, ip1, ip2, ip3, d1, d2, d3
      real(rp) :: rm2, r2, rv(3)
      real(rp), dimension(:,:,:,:), allocatable :: s
      d1 = 2*obj%pbc(1)+1
      d2 = 2*obj%pbc(2)+1
      d3 = 2*obj%pbc(3)+1
      allocate(s(3, d1, d2, d3))
      do ip1 = 1, d1
         do ip2 = 1, d2
            do ip3 = 1, d3
               s(:, ip1, ip2, ip3) = matmul([real(ip1-obj%pbc(1)-1,rp), &
                  real(ip2-obj%pbc(2)-1,rp), real(ip3-obj%pbc(3)-1,rp)], &
                  obj%l_r%v)
            end do
         end do
      end do
      if (type_in == 'mod' .or. type_in == 'wan') then
         rm2 = 0.0_rp
         do ia1 = 1, obj%na
            do ia2 = 1, obj%na
               do ip1 = 1, d1
                  do ip2 = 1, d2
                     do ip3 = 1, d3
                        rv = s(:, ip1, ip2, ip3) + obj%r(ia2, :) - obj%r(ia1, :)
                        r2 = dot_product(rv, rv)
                        if (r2 > rm2) rm2 = r2
                     end do
                  end do
               end do
            end do
         end do
         rm2 = rm2 + 4*epsilon
         r_max = sqrt(rm2)
      else
         rm2 = r_max*r_max
      end if
      allocate (nn_l(obj%na))
      nn_l = 0
      do ia1 = 1, obj%na
         do ia2 = 1, obj%na
            do ip1 = 1, d1
               do ip2 = 1, d2
                  do ip3 = 1, d3
                     rv = s(:, ip1, ip2, ip3) + obj%r(ia2, :) - obj%r(ia1, :)
                     r2 = dot_product(rv, rv)
                     if (r2 < rm2 .and. r2 > epsilon) nn_l(ia1) = nn_l(ia1) + 1
                  end do
               end do
            end do
         end do
      end do
      obj%nn_max = maxval(nn_l)
      allocate (ian_l(obj%na, obj%nn_max), &
                iapbc_l(obj%na, obj%na, d1, d2, d3), &
                rn_l(obj%na, obj%nn_max, 3))
      iapbc_l = 0
      do ia1 = 1, obj%na
         in_l = 0
         do ia2 = 1, obj%na
            do ip1 = 1, d1
               do ip2 = 1, d2
                  do ip3 = 1, d3
                     rv = s(:, ip1, ip2, ip3) + obj%r(ia2, :) - obj%r(ia1, :)
                     r2 = dot_product(rv, rv)
                     if (r2 < rm2 .and. r2 > epsilon) then
                        in_l = in_l + 1
                        ian_l(ia1, in_l) = ia2
                        iapbc_l(ia1, ia2, ip1, ip2, ip3) = in_l
                        rn_l(ia1, in_l, :) = rv
                     end if
                  end do
               end do
            end do
         end do
      end do
      call move_alloc(nn_l, obj%nn)
      call move_alloc(ian_l, obj%ian2ia)
      call move_alloc(iapbc_l, obj%iapbc2in)
      call move_alloc(rn_l, obj%rn)
      if (allocated(s)) deallocate(s)
   end subroutine calculate_neighbours

   subroutine calculate_nel(obj)
      class(atom), intent(inout) :: obj
      integer :: ia
      obj%nel = 0.0_rp
      do ia = 1, obj%na
         obj%nel = obj%nel + obj%e%q(obj%ia2ie(ia))
      end do
   end subroutine calculate_nel

   subroutine calculate_spin(obj)
      class(atom), intent(inout) :: obj
      select case (obj%ns)
      case (1)
         obj%nsl = 1
         obj%nsp = 1
         obj%g_s = 2
      case (2, 4)
         obj%nsl = merge(2, 1, obj%ns == 2)
         obj%nsp = 2
         obj%g_s = 1
      end select
   end subroutine calculate_spin

   !> Check the validity of the listing
   subroutine check_listing(l)
      character(len=*), intent(in) :: l
      if (l /= 'by_atom' .and. l /= 'by_tag') error stop 'listing error'
   end subroutine check_listing

   !> Check the validity of the number of tags
   subroutine check_ntag(na, ntag)
      integer, intent(in) :: na, ntag
      if (ntag < 1 .or. ntag > na) error stop 'ntag error'
   end subroutine check_ntag

   !> Check the validity of the position coordinates
   subroutine check_r_coord(rc)
      character(len=*), intent(in) :: rc
      if (rc /= 'cartesian' .and. rc /= 'direct') error stop 'r_coord error'
   end subroutine check_r_coord

   subroutine check_direction(d)
      character(len=*), intent(in) :: d
      if (d /= 'from' .and. d /= 'to') error stop 'direction error'
   end subroutine check_direction
   
   !> Check the validity of the magnetization coordinates
   subroutine check_m_coord(mc)
      character(len=*), intent(in) :: mc
      if (mc /= 'cartesian' .and. mc /= 'spherical') error stop 'm_coord error'
   end subroutine check_m_coord

   subroutine check_stag(na, stag)
      integer, intent(in) :: na
      integer, dimension(:), intent(in) :: stag
      if (sum(stag) /= na) error stop 'stag error'
   end subroutine check_stag

   subroutine convert_coordinates(obj, direction, m_coord)
      class(atom), intent(inout) :: obj
      character(len=*), intent(in) :: direction, m_coord
      character(len=9) :: c1, c2
      real(rp) :: mx, my, mz, mr, mt, mp
      integer :: ia
      call check_direction(direction)
      if (direction == 'from') then
         c1 = m_coord
         c2 = obj%m_coord
      else
         c1 = obj%m_coord
         c2 = m_coord
      end if
      if (trim(c1) == 'cartesian' .and. trim(c2) == 'spherical') then
         do ia = 1, obj%na
            mx = obj%m(ia, 1)
            my = obj%m(ia, 2)
            mz = obj%m(ia, 3)
            mr = sqrt(mx*mx + my*my + mz*mz)
            mt = acos(mz/mr)
            mp = atan(my/mx)
            obj%m(ia, :) = (/mr, mt, mp/)
         end do
      else if (trim(c1) == 'spherical' .and. trim(c2) == 'cartesian') then
         do ia = 1, obj%na
            mr = obj%m(ia, 1)
            mt = obj%m(ia, 2)
            mp = obj%m(ia, 3)
            mx = mr*sin(mt)*cos(mp)
            my = mr*sin(mt)*sin(mp)
            mz = mr*cos(mt)
            obj%m(ia, :) = (/mx, my, mz/)
         end do
      end if
   end subroutine convert_coordinates

   subroutine initialize_dyn(na, r, p, m)
      integer, intent(in) :: na
      real(rp), dimension(:, :), allocatable, intent(out) :: r, p, m
      allocate (r(na, 3), p(na, 3), m(na, 3))
      r = 0.0_rp
      p = 0.0_rp
      m = 0.0_rp
   end subroutine initialize_dyn

   subroutine initialize_ia2ie(ntag, stag, tag_loc, ia2ie)
      integer, intent(in) :: ntag
      integer, dimension(ntag), intent(in) :: stag
      character(len=*), dimension(ntag), intent(in) :: tag_loc
      integer, dimension(:), allocatable, intent(out) :: ia2ie
      integer, dimension(ntag) :: it2e; integer :: i1, i2, it
      allocate (ia2ie(sum(stag)))
      it2e = build_ia2ie(tag_loc)
      i1 = 1
      do it = 1, ntag
         i2 = i1 + stag(it) - 1
         ia2ie(i1:i2) = spread(it2e(it), 1, stag(it))
         i1 = i2 + 1
      end do
   end subroutine initialize_ia2ie

   subroutine initialize_pbc(pbc_v, ks_v)
      integer, dimension(3), intent(out) :: pbc_v
      real(rp), dimension(3), intent(out) :: ks_v
      pbc_v = 5; ks_v = 0.0_rp
   end subroutine initialize_pbc

   subroutine initialize_pen(na, m_pen, lambda_pen, b_pen)
      integer, intent(in) :: na
      real(rp), dimension(:), allocatable, intent(out) :: lambda_pen
      real(rp), dimension(:, :), allocatable, intent(out) :: b_pen, m_pen
      allocate (lambda_pen(na), b_pen(na, 3), m_pen(na, 3))
      lambda_pen = 0.0_rp; b_pen = 0.0_rp; m_pen = 0.0_rp
   end subroutine initialize_pen

   subroutine initialize_tag(na, ntag, stag, tag_loc)
      integer, intent(in) :: na, ntag
      integer, dimension(:), allocatable, intent(out) :: stag
      character(len=sl), dimension(:), allocatable, intent(out) :: tag_loc
      allocate (stag(ntag), tag_loc(ntag))
      if (ntag == 1) then; stag = na; elseif (ntag == na) then; stag = 1; end if
   end subroutine initialize_tag

   !> Read object in text format from file (default: 'in_atom.txt')
   subroutine read_txt(obj, file)
      class(atom), intent(inout) :: obj
      character(len=*), intent(in), optional :: file
      character(len=:), allocatable :: file_rt
      integer :: iostatus
      ! Namelist variables
      integer :: ns, na, ntag
      integer, dimension(:), allocatable :: stag
      character(len=sl), dimension(:), allocatable :: tag
      integer, dimension(:), allocatable :: ia2ie
      integer, dimension(3) :: pbc
      real(rp), dimension(3) :: k_spiral
      character(len=9) :: r_coord
      real(rp), dimension(:, :), allocatable :: r
      real(rp), dimension(:, :), allocatable :: p
      character(len=7) :: m_listing
      character(len=9) :: m_coord
      real(rp), dimension(:, :), allocatable :: m, m_pen
      character(len=7) :: lambda_pen_listing
      real(rp), dimension(:), allocatable :: lambda_pen
      real(rp), dimension(:, :), allocatable :: b_pen
      logical :: isopen
      ! Namelist
      namelist /atom/ ns, na, ntag, stag, tag, ia2ie, pbc, k_spiral, r_coord, &
         m_coord, r, p, m_listing, m, lambda_pen_listing, lambda_pen
      ! Local variables
      integer :: ia1, ia2, itag

      if (present(file)) then
         file_rt = trim(file)
      else
         file_rt = 'in_atom.txt'
      end if

      inquire (unit=10, opened=isopen)
      if (isopen) then
         write (error_unit, '(a)') 'atom%read_txt() : Unit 10 is already open'
         error stop
      else
         open (unit=10, file=file_rt, action='read', iostat=iostatus, status='old')
      end if
      if (iostatus /= 0) then
         write (error_unit, *) 'atom%read_txt(): file ', file_rt, ' not found'
         error stop
      end if

      call initialize_pbc(pbc, k_spiral)
      r_coord = 'direct'
      m_coord = 'spherical'
      m_listing = 'by_atom'
      lambda_pen_listing = 'by_atom'
      allocate (stag(0))
      read (10, nml=atom, iostat=iostatus)
      deallocate (stag)
      call check_ntag(na, ntag)
      call initialize_tag(na, ntag, stag, tag)
      call initialize_dyn(na, r, p, m)
      call initialize_pen(na, m_pen, lambda_pen, b_pen)
      allocate (ia2ie(0))
      rewind (10)
      read (10, nml=atom, iostat=iostatus)
      deallocate (ia2ie)
      call check_stag(na, stag)
      call initialize_ia2ie(ntag, stag, tag, ia2ie)
      rewind (10)
      read (10, nml=atom)
      r_coord = trim(lower(r_coord))
      call check_r_coord(r_coord)
      if (r_coord == 'cartesian') then
         r = r*obj%u%convert_length('to', 'hau')
      else !if(trim(r_coord) == 'direct') then
         r_coord = 'cartesian'
         r = obj%l_r%dir2cart(r)
      end if
      m_coord = trim(lower(m_coord))
      call check_m_coord(m_coord)
      m_listing = trim(lower(m_listing))
      call check_listing(m_listing)
      if (m_listing == 'by_tag') then
         ia2 = na
         do itag = ntag, 1, -1
            ia1 = ia2 - stag(itag) + 1
            m(ia1:ia2, :) = spread(m(itag, :), 1, stag(itag))
            ia2 = ia1 - 1
         end do
      end if
      lambda_pen_listing = trim(lower(lambda_pen_listing))
      call check_listing(lambda_pen_listing)
      if (lambda_pen_listing == 'by_tag') then
         ia2 = na
         do itag = ntag, 1, -1
            ia1 = ia2 - stag(itag) + 1
            lambda_pen(ia1:ia2) = spread(lambda_pen(itag), 1, stag(itag))
            ia2 = ia1 - 1
         end do
      end if

      obj%ns = ns
      obj%na = na
      obj%ntag = ntag
      call move_alloc(stag, obj%stag)
      call move_alloc(tag, obj%tag)
      call move_alloc(ia2ie, obj%ia2ie)
      obj%pbc = pbc
      obj%k_spiral = k_spiral
      obj%r_coord = r_coord

      call move_alloc(r, obj%r)
      call move_alloc(p, obj%p)
      call move_alloc(m, obj%m)
      call move_alloc(m_pen, obj%m_pen)
      if ((ns == 4) .and. (m_coord == 'spherical')) then
         obj%m(:, 2:3) = obj%m(:, 2:3)*deg2rad
      end if

      obj%m_coord = 'spherical' ! internal representation is spherical
      if ((ns == 4) .and. (m_coord == 'cartesian')) then
         call convert_coordinates(obj, 'from', 'cartesian')
      end if
      obj%m_pen = obj%m
      lambda_pen = lambda_pen*obj%u%convert_energy('to', 'hau')
      call move_alloc(lambda_pen, obj%lambda_pen)
      call move_alloc(b_pen, obj%b_pen)

      call obj%calculate_nel()
      call obj%calculate_spin()

      close (unit=10)
      !deallocate(file_rt)
   end subroutine read_txt

   !> For timestep t, write lammps file to unit (default: 10), if it's a file
   !> its name is set to file (default: 'out_lammps.lammpstrj')
   subroutine write_lammps(obj, t, file, unit)
      class(atom), intent(in) :: obj; real(rp), intent(in) :: t
      character(len=*), intent(in), optional :: file
      character(len=:), allocatable :: f_rt; integer, intent(in), optional :: unit
      integer :: u_rt, ia; real(rp) :: lx, xy, xz, ly, yz, lz, xlo, xhi, ylo, &
         yhi, zlo, zhi, xb, xh, yb, yh, zb, zh, m_v(3)
      f_rt = 'out_lammps.lammpstrj'
      if (present(file)) f_rt = file
      u_rt = merge(unit, 10, present(unit))
      if (.not. present(unit)) &
         open (unit=u_rt, file=f_rt, action='write', position='append')
      write (u_rt, '(a)') 'ITEM: TIMESTEP'
      write (u_rt, '(a)') real2str(t*obj%u%convert_time('from', 'hau'))
      write (u_rt, '(a)') 'ITEM: NUMBER OF ATOMS'
      write (u_rt, '(a)') int2str(obj%na)
      write (u_rt, '(a)') 'ITEM: BOX BOUNDS xy xz yz'
      lx = obj%l_r%v(1,1); xy = obj%l_r%v(1,2); xz = obj%l_r%v(1,3)
      ly = obj%l_r%v(2,2); yz = obj%l_r%v(2,3); lz = obj%l_r%v(3,3)
      xlo = 0.0_rp; ylo = 0.0_rp; zlo = 0.0_rp
      do ia = 1, obj%na
         xlo = min(xlo, obj%r(ia, 1))
         ylo = min(ylo, obj%r(ia, 2))
         zlo = min(zlo, obj%r(ia, 3))
      end do
      xhi = xlo+lx; yhi = ylo+ly; zhi = zlo+lz
      xb = xlo+min(0.0_rp, min(xy, min(xz, xy+xz)))
      xh = xhi+max(0.0_rp, max(xy, max(xz, xy+xz)))
      yb = ylo+min(0.0_rp, yz); yh = yhi+max(0.0_rp, yz); zb = zlo; zh = zhi
      xb = xb*obj%u%convert_length('from', 'hau')
      xh = xh*obj%u%convert_length('from', 'hau')
      yb = yb*obj%u%convert_length('from', 'hau')
      yh = yh*obj%u%convert_length('from', 'hau')
      zb = zb*obj%u%convert_length('from', 'hau')
      zh = zh*obj%u%convert_length('from', 'hau')
      write (u_rt, '(a)') real2str(xb)//' '//real2str(xh)//' '//real2str(xy)
      write (u_rt, '(a)') real2str(yb)//' '//real2str(yh)//' '//real2str(xz)
      write (u_rt, '(a)') real2str(zb)//' '//real2str(zh)//' '//real2str(yz)
      write (u_rt, '(a)') 'ITEM: ATOMS type x y z vx vy vz'
      do ia = 1, obj%na; m_v = sph2cart(obj%m(ia, :))
         write (u_rt, '(a)') int2str(obj%ia2ie(ia))//' '// &
            real2str(obj%r(ia,1)*obj%u%convert_length('from','hau'))//' '// &
            real2str(obj%r(ia,2)*obj%u%convert_length('from','hau'))//' '// &
            real2str(obj%r(ia,3)*obj%u%convert_length('from','hau'))//' '// &
            real2str(m_v(1))//' '//real2str(m_v(2))//' '//real2str(m_v(3))
      end do
      call TBKOSTER_flush(u_rt); if (.not. present(unit)) close (u_rt)
   end subroutine write_lammps

   !> Write object in text format to unit (default: 10), if it's a file
   !> its name is set to file (default: 'out_atom.txt')
   subroutine write_txt(obj, file, unit)
     class(atom), intent(in) :: obj
      character(len=*), intent(in), optional :: file
      character(len=:), allocatable         :: file_rt
      integer, intent(in), optional :: unit
      integer                     :: unit_rt
      ! Namelist variables
      integer :: ns, na, ntag
      integer, dimension(obj%ntag) :: stag
      character(len=sl), dimension(obj%ntag) :: tag
      integer, dimension(obj%na) :: ia2ie
      integer, dimension(3) :: pbc
      real(rp), dimension(3) :: k_spiral
      character(len=len(obj%r_coord)) :: r_coord
      real(rp), dimension(obj%na, 3) :: r, p, m
      character(len=len(obj%m_coord)) :: m_coord
      real(rp), dimension(obj%na)   :: lambda_pen
      ! Namelist
      namelist /atom/ ns, na, ntag, stag, tag, ia2ie, pbc, k_spiral, r_coord, &
         m_coord, r, p, m, lambda_pen

      if (present(file)) then
         file_rt = file
      else
         file_rt = 'out_atom.txt'
      end if
      if (present(unit)) then
         unit_rt = unit
      else
         unit_rt = 10
      end if

      if (.not. present(unit)) then
         open (unit=unit_rt, file=file_rt, action='write')
      end if

      ns = obj%ns
      na = obj%na
      ntag = obj%ntag
      stag = obj%stag
      tag = obj%tag
      ia2ie = obj%ia2ie
      pbc = obj%pbc
      r_coord = obj%r_coord
      if (obj%r_coord == 'cartesian') then
         r = obj%r*obj%u%convert_length('from', 'hau')
      else
         r = obj%r
      end if
      p = obj%p
      m_coord = obj%m_coord
      m = obj%m
      if ((obj%ns == 4) .and. (obj%m_coord == 'spherical')) then
         m(:, 2:3) = m(:, 2:3)*rad2deg
      end if

      lambda_pen = obj%lambda_pen*obj%u%convert_energy('from', 'hau')

      write (unit_rt, nml=atom)

      call TBKOSTER_flush(unit_rt); if (.not. present(unit)) close (unit_rt)
   end subroutine write_txt

   !> Write property (default: property_list) in text format to unit
   !> (default: 10), if it's a file its name is set to file (default:
   !> 'out_atom.txt'), if tag (default: .true.) the namelist opening and closing
   !> tags are written
   subroutine write_txt_formatted(obj, file, property, tag, unit, access)
      class(atom), intent(in) :: obj
      character(len=*), intent(in), optional :: file, access
      character(len=*), dimension(:), intent(in), optional :: property
      character(len=:), allocatable :: f_rt, p_rt(:)
      logical, intent(in), optional :: tag; logical :: t_rt
      integer, intent(in), optional :: unit
      integer :: u_rt, ia_f, in_f, ip_f, it_f
      real(rp), allocatable :: r_f(:, :), rn_f(:, :, :), lp_f(:)
      f_rt = 'out_atom.txt'
      if (present(file)) f_rt = file
      p_rt = property_list
      if (present(property)) p_rt = property
      t_rt = .true.
      if (present(tag)) t_rt = tag
      u_rt = merge(unit, 10, present(unit))
      if (.not. present(unit)) then
         if (present(access)) then
            open (unit=u_rt, file=f_rt, action='write', access=access)
         else
            open (unit=u_rt, file=f_rt, action='write')
         end if
      end if
      if (t_rt) write (u_rt, '(a)') '&atom'
      do ip_f = 1, size(p_rt)
         select case (lower(trim(p_rt(ip_f))))
         case ('ns'); call wi('ns', obj%ns)
         case ('na'); call wi('na', obj%na)
         case ('ntag'); call wi('ntag', obj%ntag)
         case ('stag')
            do it_f = 1, obj%ntag
               call wi('stag('//int2str(it_f)//')', obj%stag(it_f))
            end do
         case ('tag')
            do it_f = 1, obj%ntag
               call ws('tag('//int2str(it_f)//')', &
                  ''''//trim(obj%tag(it_f))//'''')
            end do
         case ('ia2ie')
            do ia_f = 1, obj%na
               call wi('ia2ie('//int2str(ia_f)//')', obj%ia2ie(ia_f))
            end do
         case ('nel'); call wr('nel', obj%nel)
         case ('pbc'); call wv('pbc', [real(obj%pbc,rp)])
         case ('k_spiral'); call wv('k_spiral', obj%k_spiral)
         case ('r_coord'); call ws('r_coord', ''''//obj%r_coord//'''')
         case ('m_coord'); call ws('m_coord', ''''//obj%m_coord//'''')
         case ('r'); allocate(r_f(obj%na,3)); r_f = obj%r
            if (obj%r_coord == 'cartesian') &
               r_f = r_f*obj%u%convert_length('from', 'hau')
            do ia_f = 1, obj%na
               call wfx('r('//int2str(ia_f)//',:)', r_f(ia_f, :))
            end do; deallocate(r_f)
         case ('p'); do ia_f = 1, obj%na
               call wv('p('//int2str(ia_f)//',:)', obj%p(ia_f, :))
            end do
         case ('m'); do ia_f = 1, obj%na
               call wv('m('//int2str(ia_f)//',:)', &
                  [obj%m(ia_f, 1), obj%m(ia_f, 2:3)*rad2deg])
            end do
         case ('nn'); do ia_f = 1, obj%na
               call wi('nn('//int2str(ia_f)//')', obj%nn(ia_f))
            end do
         case ('nn_max'); call wi('nn_max', obj%nn_max)
         case ('ian2ia')
            do ia_f = 1, obj%na
               do in_f = 1, obj%nn(ia_f)
                  call wi('ian2ia('//int2str(ia_f)//','//int2str(in_f)//')', &
                     obj%ian2ia(ia_f, in_f))
               end do
            end do
         case ('rn'); allocate(rn_f(obj%na,obj%nn_max,3))
            rn_f = obj%rn*obj%u%convert_length('from', 'hau')
            do ia_f = 1, obj%na
               do in_f = 1, obj%nn(ia_f)
                  call wv('rn('//int2str(ia_f)//','//int2str(in_f)//',:)', &
                     rn_f(ia_f, in_f, :))
               end do
            end do; deallocate(rn_f)
         case ('lambda_pen'); allocate(lp_f(obj%na))
            lp_f = obj%lambda_pen*obj%u%convert_energy('from', 'hau')
            do ia_f = 1, obj%na
               call wr('lambda_pen('//int2str(ia_f)//')', lp_f(ia_f))
            end do; deallocate(lp_f)
         case ('b_pen')
            do ia_f = 1, obj%na
               call wv('b_pen('//int2str(ia_f)//',:)', obj%b_pen(ia_f, :))
            end do
         end select
      end do
      if (t_rt) write (u_rt, '(a)') ' /'
      call TBKOSTER_flush(u_rt); if (.not. present(unit)) close (u_rt)
   contains
      subroutine ws(n, v); character(*), intent(in) :: n, v
         write (u_rt, '(a)') ' '//n//' = '//v
      end subroutine
      subroutine wi(n, v); character(*), intent(in) :: n; integer, intent(in) :: v
         call ws(n, int2str(v))
      end subroutine
      subroutine wr(n, v); character(*), intent(in) :: n; real(rp), intent(in) :: v
         call ws(n, real2str(v))
      end subroutine
      subroutine wv(n, v); character(*), intent(in) :: n; real(rp), intent(in) :: v(:)
         integer :: k; character(len=256) :: s
         s = real2str(v(1))
         do k = 2, size(v); s = trim(s)//', '//real2str(v(k)); end do
         call ws(n, s)
      end subroutine
      subroutine wfx(n, v); character(*), intent(in) :: n; real(rp), intent(in) :: v(:)
         integer :: k; character(len=256) :: s
         s = fixedreal2str(v(1))
         do k = 2, size(v); s = trim(s)//', '//fixedreal2str(v(k)); end do
         call ws(n, s)
      end subroutine
   end subroutine write_txt_formatted

   !> Write object in extended .xyz format to unit (default: 10), if it's a file
   !> its name is set to file (default: 'out_atom.xyz')
   subroutine write_xyz(obj, file, unit)
      class(atom), intent(in) :: obj; character(len=*), intent(in), optional :: file
      character(len=:), allocatable :: f_rt; integer, intent(in), optional :: unit
      integer :: u_rt, ia; real(rp) :: mc(3), lm(3)
      f_rt = 'out_atom.xyz'
      if (present(file)) f_rt = file
      u_rt = merge(unit, 10, present(unit))
      if (.not. present(unit)) open (unit=u_rt, file=f_rt, action='write')
      write (u_rt, '(a)') int2str(obj%na)
      write (u_rt, '(a)') 'Lattice="'// &
         real2str(obj%l_r%v(1,1)*obj%u%convert_length('from','hau'))//' '// &
         real2str(obj%l_r%v(1,2)*obj%u%convert_length('from','hau'))//' '// &
         real2str(obj%l_r%v(1,3)*obj%u%convert_length('from','hau'))//' '// &
         real2str(obj%l_r%v(2,1)*obj%u%convert_length('from','hau'))//' '// &
         real2str(obj%l_r%v(2,2)*obj%u%convert_length('from','hau'))//' '// &
         real2str(obj%l_r%v(2,3)*obj%u%convert_length('from','hau'))//' '// &
         real2str(obj%l_r%v(3,1)*obj%u%convert_length('from','hau'))//' '// &
         real2str(obj%l_r%v(3,2)*obj%u%convert_length('from','hau'))//' '// &
         real2str(obj%l_r%v(3,3)*obj%u%convert_length('from','hau'))//'" '// &
         'Properties=species:S:1:pos:R:3:dipoles:R:3'
      do ia = 1, obj%na; lm = obj%m(ia, :); mc = sph2cart(lm)
         write (u_rt, '(a)') obj%e%symbol(obj%ia2ie(ia))//' '// &
            real2str(obj%r(ia,1)*obj%u%convert_length('from','hau'))//' '// &
            real2str(obj%r(ia,2)*obj%u%convert_length('from','hau'))//' '// &
            real2str(obj%r(ia,3)*obj%u%convert_length('from','hau'))//' '// &
            real2str(mc(1))//' '//real2str(mc(2))//' '//real2str(mc(3))
      end do
      call TBKOSTER_flush(u_rt); if (.not. present(unit)) close (u_rt)
   end subroutine write_xyz
end module atom_mod
