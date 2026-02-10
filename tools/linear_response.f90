!
! Copyright (C) 2019
! Cyrille Barreteau <mailto:cyrille.barreteau@cea.fr>,
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
!  bands.f90
!  TBKOSTER
program linear_response
   use, intrinsic :: iso_fortran_env, only: output_unit
   use precision_mod
   use string_mod
   implicit none
   integer, parameter :: unit_mesh = 11
   integer, parameter :: unit_mesh_out = 111
   integer, parameter :: unit_hamiltonian_tb = 12
   integer, parameter :: unit_LR_in = 13
   integer, parameter :: unit_band_out = 14
   integer, parameter :: unit_log_out = 15
   integer :: iostatus, nx,ia,na, ih_min,ih_max,ih, nh, ns, nsl, no_max, n_band,na_band
   logical :: isopen, close_file
   character(len=*), parameter :: dir = 'band/'
   character(len=*), parameter :: file_LR_in = 'in_LR.txt'
   character(len=*), parameter :: file_band_out = 'out_band.txt'
   character(len=*), parameter :: file_mesh_in = 'in_mesh.txt'
   character(len=*), parameter :: file_mesh_out = 'out_mesh.txt'
   character(len=*), parameter :: file_hamiltonian_tb = 'out_hamiltonian_tb.txt'
   character(len=*), parameter :: file_log_out = 'out_log.txt'
   character(len=9) :: x_coord
   character(len=4) :: type, type_case
   character(len=80) :: line
   character(len=20) :: proj
   integer, dimension(:), allocatable :: ih_band,ia_band, iband2io
   real(rp), dimension(:, :, :), allocatable :: en_k
   real(rp), dimension(:), allocatable :: w
   real(rp), dimension(:, :), allocatable :: x
   real(rp) :: v_factor
   real(rp), dimension(3, 3) :: v, vrec
   real(rp), dimension(:, :, :, :), allocatable :: w_band_spin
   real(rp), dimension(:, :, :, :), allocatable :: w_band_orb
   real(rp), dimension(:, :, :, :), allocatable :: band_velocity
   real(rp) :: en_min, en_max, en_f, Eref, E_min, E_max
   namelist /lattice/ v_factor, v, vrec
   namelist /atom/ na
   namelist /mesh/ type
   namelist /mesh_out/ type, nx, x_coord, x, w
   namelist /hamiltonian_tb/ nh, ns
   namelist /LR/ proj,Eref, E_min,E_max, ih_min,ih_max, na_band, ia_band
   namelist /band_out/ en_k, iband2io, w_band_spin, w_band_orb, band_velocity

   ! inquire(unit=unit_energy_in,opened=isopen)
   ! if (isopen) then
   !   write(*,'(a)') 'energy%read_txt() : unit 14 is already open'
   !   close(unit_energy_in)
   ! else
   !   open(unit=unit_energy_in,file=dir//file_energy_in,action='read',iostat=iostatus,status='old')
   ! end if
   ! if(iostatus /= 0) then
   !   write(*,*) 'energy%read_txt(): file ', file_energy_in, ' not found'
   !   error stop
   !  end if
   open (unit_LR_in, file=dir//file_LR_in, action='read', iostat=iostatus, status='old')
   open (unit_band_out, file=dir//file_band_out, action='read', iostat=iostatus, status='old')
   open (unit_mesh, file=dir//file_mesh_in, action='read', iostat=iostatus, status='old')
   open (unit_mesh_out, file=dir//file_mesh_out, action='read', iostat=iostatus, status='old')
   open (unit_hamiltonian_tb, file=dir//file_hamiltonian_tb, action='read', iostat=iostatus, status='old')
   open (unit_log_out, file=dir//file_log_out, action='read', iostat=iostatus, status='old')

   read (unit_log_out, nml=lattice, iostat=iostatus)
   rewind(unit_log_out)
   read (unit_log_out, nml=atom, iostat=iostatus)
   close (unit_log_out)
   ! read k mesh in (test type of mesh)
   read (unit_mesh, nml=mesh, iostat=iostatus)
   type_case = lower(type)
   close (unit_mesh)
   ! read k mesh out
   allocate (x(0, 0), w(0))
   read (unit_mesh_out, nml=mesh_out, iostat=iostatus)
   deallocate (x, w)
   select case (type_case)
   case ('path')
      write (*, *) 'path is not suited for linear response analysis'
      stop
   case ('list')
      write (*, *) 'list is  suited for linear response analysis but beware to have a good sampling'
   case ('mp')
      write (*, *) 'monkhorst pack is suited for linear response analysis'
   end select

   allocate (x(nx, 3))
   allocate (w(nx))
   x = 0.0_rp
   w = 1.0_rp/nx
   rewind (unit_mesh_out)
   read (unit_mesh_out, nml=mesh_out)
   close (unit_mesh_out)
   x_coord = lower(x_coord)
   if (x_coord == 'direct') then
      write (*, *) 'x_coord must be cartesian'
      stop
   end if

   close (unit_mesh_out)

   ! read information about the hamiltonian size
   read (unit_hamiltonian_tb, nml=hamiltonian_tb, iostat=iostatus)
   close (unit_hamiltonian_tb)

   select case (ns)
   case (1)
      nsl = 1
    write (*, *) 'ns=1 not suited to linear-response it should be equal to 4 (non-collinear)'
      stop   
   case (2)
      nsl = 2
     write (*, *) 'ns=2 not suited to linear-response it should be equal to 4 (non-collinear)'
      stop    
   case (4)
      nsl = 1
      write (*, *) 'ns=4 is suited to linear-response'
   end select
   
   E_min=-1.0_rp
   E_max=1.0_rp

   na_band = 0
   allocate (ia_band(0))
   read (unit_LR_in, nml=LR, iostat=iostatus)
   deallocate (ia_band)
   allocate (ia_band(na_band))
   rewind (unit_LR_in)
   read (unit_LR_in, nml=LR, iostat=iostatus)
       if(na_band > 0.and.na_band.ne.na) then
    rewind(unit_LR_in)
    read(unit_LR_in,nml=LR,iostat=iostatus)
    rewind(unit_LR_in)
    elseif(na_band==na) then
      do ia=1,na
      ia_band(ia)=ia
      end do
    end if
   read (unit_LR_in, nml=LR, iostat=iostatus)
   if(TRIM(proj).ne.'linear-response') then ! 
     stop
   endif
   ih_min=0
   ih_max=0
   rewind(unit_LR_in)
    read(unit_LR_in,nml=LR,iostat=iostatus)
    if(ih_min/=0.and.ih_max/=0) then
    allocate(ih_band(ih_max-ih_min+1))
       do ih=1,ih_max-ih_min+1
         ih_band(ih)=ih+ih_min-1
        end do
    elseif(ih_min==0.and.ih_max==0) then
       ih_min=1
       ih_max=nh
     allocate(ih_band(ih_max-ih_min+1))
       do ih=1,ih_max-ih_min+1
         ih_band(ih)=ih
        end do
    endif 
    n_band=ih_max-ih_min+1
   close (unit_LR_in)

   ! read eigenvalues and different types of projections: sites or spin or orbital
   allocate (en_k(nh, nx, nsl))
   read (unit_band_out, nml=band_out, iostat=iostatus)

   rewind (unit_band_out)
      allocate (w_band_spin(0, 0, 0, 0))
      allocate (w_band_orb(0, 0, 0, 0))
      read (unit_band_out, nml=band_out, iostat=iostatus)
      deallocate (w_band_spin)
      deallocate (w_band_orb)
      rewind (unit_band_out)
      allocate (w_band_spin(n_band, nx, 0:na_band, 4))
      allocate (w_band_orb(n_band, nx, 0:na_band, 4))
      allocate (band_velocity(n_band, nx, nsl, 4))
      read (unit_band_out, nml=band_out, iostat=iostatus)
   close (unit_band_out)
   ! set the Fermi energy to zero
   call get_fermi_scf(en_f)
   en_k = en_k - en_f

   call build_linear_response(x, w, en_k, Eref, E_min,E_max, w_band_spin, w_band_orb, band_velocity,&
                               ih_band, n_band,na_band, nh,  nx, nsl)
contains

    subroutine build_linear_response(x, w, en_k, Eref, E_min, E_max,w_band_spin,w_band_orb, band_velocity,&
                                    ih_band, n_band,na_band, nh,  nx, nsl)
      use precision_mod
      use string_mod
      implicit none
      integer, intent(in) :: n_band, na_band, nsl, nh, nx
      integer, intent(in) :: ih_band(n_band)
      real(rp), intent(in) :: x(nx, 3), w(nx), en_k(nh, nx, nsl)
      real(rp), intent(in) :: w_band_spin(n_band, nx, 0:na_band, 4),w_band_orb(n_band, nx, 0:na_band, 4)
      real(rp), intent(in) :: band_velocity(n_band, nx, nsl, 4)
      real(rp), intent(in) :: Eref, E_min, E_max
      logical :: isopen
      integer :: ih, ix, ie, ir, ia_band,ne_max
      logical :: close_file
      integer, parameter :: unit_spin_E = 10, unit_orb_E= 11, unit_spin_site = 12, unit_orb_site= 13
      real(rp) :: energ
      real(rp), dimension(3) :: LR_spin,LR_orb
      character(len=*), parameter :: dir = 'band/'
      character(len=*), parameter :: file_spin_E = 'LR_spin_E.dat'
      character(len=*), parameter :: file_orb_E  = 'LR_orb_E.dat'
      character(len=*), parameter :: file_spin_site = 'LR_spin_site.dat'
      character(len=*), parameter :: file_orb_site  = 'LR_orb_site.dat'
      character(len=80) :: fmt

      inquire (unit=unit_spin_E, opened=isopen)
      if (.not. isopen) then
         open (unit=unit_spin_E, file=dir//file_spin_E, action='write')
      end if
      inquire (unit=unit_orb_E, opened=isopen)
      if (.not. isopen) then
         open (unit=unit_orb_E, file=dir//file_orb_E, action='write')
      end if
      inquire (unit=unit_spin_site, opened=isopen)
      if (.not. isopen) then
         open (unit=unit_spin_site, file=dir//file_spin_site, action='write')
      end if
      inquire (unit=unit_orb_site, opened=isopen)
      if (.not. isopen) then
         open (unit=unit_orb_site, file=dir//file_orb_site, action='write')
      end if

      write(*,*) 'n_band= ',n_band
      write (unit_spin_E, *) '@# spin linear response vs energy'
      write (unit_orb_E, *)  '@# orbit linear response vs energy'
      write (unit_spin_site, *) '@# spin linear response vs site'
      write (unit_orb_site,*)  '@# orbit linear response vs site'

      fmt = trim('(4f12.7)')
      ne_max=151
      do ie=1,ne_max
         energ=E_min+(E_max-E_min)*(ie-1)/(ne_max-1)
         LR_spin=0.0_rp
         LR_orb=0.0_rp
         do ix = 1, nx
            do ih = 1, n_band
              do ir=1,3
                 LR_spin(ir)=LR_spin(ir) +w(ix)* w_band_spin(ih, ix, 0, ir+1)*band_velocity(ih, ix, 1, 2)&
                                          *delta_function(en_k(ih_band(ih), ix, 1) - energ)
                 LR_orb(ir)=LR_orb(ir) +w(ix)* w_band_orb(ih, ix, 0, ir+1)*band_velocity(ih, ix, 1, 2)&
                                          *delta_function(en_k(ih_band(ih), ix, 1) - energ)
              end do
            end do
         end do
         write (unit_spin_E, fmt) energ,LR_spin(:)
         write (unit_orb_E, fmt)  energ,LR_orb(:)
      end do

     close (unit_spin_E)
     close (unit_orb_E)

     fmt = trim('(I2,3f12.7)')
      do ia_band = 1, na_band
         LR_spin=0.0_rp
         LR_orb=0.0_rp
         do ix = 1, nx
            do ih = 1, n_band
              do ir=1,3
                 LR_spin(ir)=LR_spin(ir) +w(ix)* w_band_spin(ih, ix, ia_band, ir+1)*band_velocity(ih, ix, 1, 2)&
                                          *delta_function(en_k(ih_band(ih), ix, 1) - Eref)
                 LR_orb(ir)=LR_orb(ir) +w(ix)* w_band_orb(ih, ix, ia_band, ir+1)*band_velocity(ih, ix, 1, 2)&
                                          *delta_function(en_k(ih_band(ih), ix, 1) - Eref)
              end do
            end do
         end do
         write (unit_spin_site, fmt) ia_band,LR_spin(:)
         write (unit_orb_site, fmt)  ia_band,LR_orb(:)
      end do     

     close (unit_spin_site)
     close (unit_orb_site)
   end subroutine build_linear_response

   subroutine get_fermi_scf(en_f)
      use precision_mod
      implicit none
      integer, parameter :: unit_energy_scf = 10
      integer :: iostatus
      character(len=*), parameter :: file_energy_scf = 'out_energy.txt'
      real(rp) :: en_f
      namelist /energy/ en_f

      open (unit_energy_scf, file=file_energy_scf, action='read', iostat=iostatus, status='old')
      read (unit_energy_scf, nml=energy, iostat=iostatus)

      close (unit_energy_scf)
   end subroutine get_fermi_scf

   real(rp) function delta_function(x)
      use precision_mod
    !
    ! --> 'fd': derivative of the Fermi-Dirac smearing. 0.5/(1.0+cosh(x))
    !
      real(rp), intent(in) :: x
    ! output: the value of the function
    ! input: the point where to compute the function
    ! local variable
      real(rp), parameter :: smearing = 0.005_rp  !smearing in eV
    ! Fermi-Dirac smearing
      if (abs(x) <= 36.0_rp) then
         delta_function = 1.0_rp/(smearing*(2.0_rp + exp(-x/smearing) + exp(+x/smearing)))
         ! in order to avoid problems for large values of x in the e
      else
         delta_function = 0.0_rp
      end if
   end function delta_function

end program linear_response
