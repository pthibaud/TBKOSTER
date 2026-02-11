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
!  atom.f90 - TBKOSTER
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
      integer :: ns
      integer :: nsl
      integer :: nsp
      integer :: g_s
      integer, dimension(2, 2) :: iss2is
      !> @}

      !> @defgroup Atom Atom-related variables
      !> @{
      integer :: na
      integer :: ntag
      integer, dimension(:), allocatable :: stag
      character(len=sl), dimension(:), allocatable :: tag
      integer, dimension(:), allocatable :: ia2ie
      real(rp) :: nel
      !> @}

      integer, dimension(3) :: pbc
      real(rp), dimension(3) :: k_spiral

      !> @defgroup Dynamics Dynamics-related variables
      !> @{
      character(len=9) :: r_coord
      real(rp), dimension(:, :), allocatable :: r
      real(rp), dimension(:, :), allocatable :: p
      character(len=9) :: m_coord
      real(rp), dimension(:, :), allocatable :: m, m_pen
      !> @}

      !> @defgroup Neighbours Neighbours-related variables
      !> @{
      integer, dimension(:), allocatable :: nn
      integer :: nn_max
      integer, dimension(:, :), allocatable :: ian2ia
      integer, dimension(:, :, :, :, :), allocatable :: iapbc2in
      real(rp), dimension(:, :, :), allocatable :: rn
      !> @}

      !> @defgroup Magnetic_penalization Magnetic penalization related variables
      !> @{
      real(rp), dimension(:), allocatable   :: lambda_pen
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

   subroutine check_listing(l)
      character(len=*), intent(in) :: l
      if (l /= 'by_atom' .and. l /= 'by_tag') error stop 'listing error'
   end subroutine check_listing

   subroutine check_ntag(na, ntag)
      integer, intent(in) :: na, ntag
      if (ntag < 1 .or. ntag > na) error stop 'ntag error'
   end subroutine check_ntag

   subroutine check_r_coord(rc)
      character(len=*), intent(in) :: rc
      if (rc /= 'cartesian' .and. rc /= 'direct') error stop 'r_coord error'
   end subroutine check_r_coord

   subroutine check_direction(d)
      character(len=*), intent(in) :: d
      if (d /= 'from' .and. d /= 'to') error stop 'direction error'
   end subroutine check_direction

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

   subroutine read_txt(obj, file)
      class(atom), intent(inout) :: obj
      character(len=*), intent(in), optional :: file
      character(len=:), allocatable :: f_rt
      integer :: ios, i1, i2, it, ns, na, ntag, pbc_l(3)
      integer, dimension(:), allocatable :: stag, ia2e
      character(len=sl), dimension(:), allocatable :: tag_l
      real(rp), dimension(3) :: ks
      character(len=9) :: rc, mc
      real(rp), dimension(:, :), allocatable :: r_l, p_l, m_l, mp, bp
      character(len=7) :: ml, lpl
      real(rp), dimension(:), allocatable :: lp_l
      logical :: op
      namelist /atom/ ns, na, ntag, stag, tag_l, ia2e, pbc_l, ks, rc, mc, &
                      r_l, p_l, ml, m_l, lpl, lp_l
      f_rt = 'in_atom.txt'
      if (present(file)) f_rt = trim(file)
      inquire (unit=10, opened=op); if (op) error stop 'Unit 10 open'
      open (unit=10, file=f_rt, action='read', iostat=ios, status='old')
      if (ios /= 0) error stop 'file not found'
      call initialize_pbc(pbc_l, ks); rc = 'direct'; mc = 'spherical'
      ml = 'by_atom'; lpl = 'by_atom'
      allocate (stag(0)); read (10, nml=atom, iostat=ios); deallocate (stag)
      call check_ntag(na, ntag); call initialize_tag(na, ntag, stag, tag_l)
      call initialize_dyn(na, r_l, p_l, m_l)
      call initialize_pen(na, mp, lp_l, bp)
      allocate (ia2e(0)); rewind (10); read (10, nml=atom, iostat=ios)
      deallocate (ia2e)
      call check_stag(na, stag); call initialize_ia2ie(ntag, stag, tag_l, ia2e)
      rewind (10); read (10, nml=atom)
      rc = trim(lower(rc)); call check_r_coord(rc)
      if (rc == 'cartesian') then
         r_l = r_l*obj%u%convert_length('to', 'hau')
      else
         rc = 'cartesian'
         r_l = obj%l_r%dir2cart(r_l)
      end if
      mc = trim(lower(mc)); call check_m_coord(mc); ml = trim(lower(ml))
      call check_listing(ml)
      if (ml == 'by_tag') then
         i2 = na
         do it = ntag, 1, -1
            i1 = i2-stag(it)+1
            m_l(i1:i2,:) = spread(m_l(it,:), 1, stag(it))
            i2 = i1-1
         end do
      end if
      lpl = trim(lower(lpl)); call check_listing(lpl)
      if (lpl == 'by_tag') then
         i2 = na
         do it = ntag, 1, -1
            i1 = i2-stag(it)+1
            lp_l(i1:i2) = spread(lp_l(it), 1, stag(it))
            i2 = i1-1
         end do
      end if
      obj%ns = ns; obj%na = na; obj%ntag = ntag
      call move_alloc(stag, obj%stag)
      call move_alloc(tag_l, obj%tag)
      call move_alloc(ia2e, obj%ia2ie)
      obj%pbc = pbc_l; obj%k_spiral = ks; obj%r_coord = rc
      call move_alloc(r_l, obj%r); call move_alloc(p_l, obj%p)
      call move_alloc(m_l, obj%m); call move_alloc(mp, obj%m_pen)
      if (ns == 4 .and. mc == 'spherical') obj%m(:, 2:3) = obj%m(:, 2:3)*deg2rad
      obj%m_coord = 'spherical'
      if (ns == 4 .and. mc == 'cartesian') &
         call convert_coordinates(obj, 'from', 'cartesian')
      obj%m_pen = obj%m; lp_l = lp_l*obj%u%convert_energy('to', 'hau')
      call move_alloc(lp_l, obj%lambda_pen); call move_alloc(bp, obj%b_pen)
      call obj%calculate_nel(); call obj%calculate_spin()
      close (10)
   end subroutine read_txt

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

   subroutine write_txt(obj, file, unit)
      class(atom), intent(in) :: obj; character(len=*), intent(in), optional :: file
      character(len=:), allocatable :: f_rt; integer, intent(in), optional :: unit
      integer :: u_rt, ns_o, na_o, ntag_o
      integer, allocatable :: s_o(:), i2e_o(:); integer :: pbc_o(3)
      character(len=sl), allocatable :: tag_o(:); real(rp) :: k_s_o(3)
      real(rp), allocatable :: r_o(:, :), p_o(:, :), m_o(:, :), lp_o(:)
      character(len=9) :: rc_o, mc_o
      namelist /atom_o/ ns_o, na_o, ntag_o, s_o, tag_o, i2e_o, pbc_o, &
                        k_s_o, rc_o, mc_o, r_o, p_o, m_o, lp_o
      f_rt = 'out_atom.txt'; if (present(file)) f_rt = file
      u_rt = merge(unit, 10, present(unit))
      if (.not. present(unit)) open (unit=u_rt, file=f_rt, action='write')
      ns_o = obj%ns; na_o = obj%na; ntag_o = obj%ntag; s_o = obj%stag
      tag_o = obj%tag; i2e_o = obj%ia2ie; pbc_o = obj%pbc
      rc_o = obj%r_coord; r_o = obj%r
      if (rc_o == 'cartesian') r_o = r_o*obj%u%convert_length('from', 'hau')
      p_o = obj%p; mc_o = obj%m_coord; m_o = obj%m
      if (obj%ns == 4 .and. mc_o == 'spherical') &
         m_o(:, 2:3) = m_o(:, 2:3)*rad2deg
      lp_o = obj%lambda_pen*obj%u%convert_energy('from', 'hau')
      write (u_rt, nml=atom_o)
      call TBKOSTER_flush(u_rt); if (.not. present(unit)) close (u_rt)
   end subroutine write_txt

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
