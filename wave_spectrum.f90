!!
!!  Copyright (C) 2011-2017  Johns Hopkins University
!!
!!  This file is part of lesgo.
!!
!!  lesgo is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  lesgo is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with lesgo.  If not, see <http://www.gnu.org/licenses/>.
!!

!*******************************************************************************
module wave_spectrum
!*******************************************************************************
use types, only : rprec

save
private rprec
public

integer :: jx, jy
real(rprec) :: grav, delta_kx, delta_ky, k_filter, omega_filter, sigma_spec
real(rprec), dimension(:,:),allocatable :: Theta_wave, rand
contains

!*******************************************************************************
subroutine spectrum_calc ()
!*******************************************************************************
! This subroutine will initialize all variables for the ocean spectrum, this will
! only be called once at the beginning of the simulation in initialize.f90.

use param, only : L_x, L_y, pi, nx, ny, dx, dy, alpha_spec, wp_spec, kp_spec, &
                  filt_spec, ld
use sim_param, only: K_wave, D_spread, S_omega_2D, S_k_2D, phi, kx_wave, ky_wave, omega_wave
implicit none

allocate (rand(ld,ny)); rand = 0.0_rprec  
allocate (Theta_wave(ld,ny)); Theta_wave = 0.0_rprec

grav = 9.81_rprec
call init_random_seed
call random_number(rand)

do jy = 1, ny
 do jx = 1, ld

        kx_wave(jx,jy) = (2_rprec*pi/L_x)*(jx-1)
        ky_wave(jx,jy) = (2_rprec*pi/L_y)*(jy-1)

        K_wave(jx,jy) = sqrt(kx_wave(jx,jy)**2 + ky_wave(jx,jy)**2)
        omega_wave(jx,jy) = sqrt(grav*K_wave(jx,jy))

        Theta_wave(jx,jy) =  atan(ky_wave(jx,jy)/kx_wave(jx,jy))

        phi(jx,jy) = rand(jx,jy)*2*pi
 end do
end do

Theta_wave(1,1) = 0._rprec

k_filter = pi/sqrt(dx*dy)
omega_filter = sqrt(grav*k_filter)

! Creating the 2D (x,y) spectrum using the JONSWAP specturm with directional spreading 
! (Same as Yang, et al 2014). This was tested and worked on 04/07/2024 MA

do jy = 1, ny
 do jx = 1, ld

        if (omega_wave(jx,jy) <= wp_spec) then
            sigma_spec = 0.07_rprec
        else 
            sigma_spec = 0.09_rprec
        end if

        if (abs(Theta_wave(jx,jy)) <= pi/2.0_rprec) then
            D_spread(jx,jy) = (2.0_rprec/pi)*cos(Theta_wave(jx,jy))**2
        else
            D_spread(jx,jy) = 0.0
        end if

        S_omega_2D(jx,jy) = (alpha_spec*grav**2) / omega_wave(jx,jy)**5 * exp(-1.25_rprec*(wp_spec/omega_wave(jx,jy))**4) * &
                         (3.3_rprec)**exp(-((omega_wave(jx,jy) - wp_spec)**2 / (2.0_rprec*sigma_spec**2*wp_spec**2)))
        S_k_2D(jx,jy) = S_omega_2D(jx,jy) * D_spread(jx,jy) * 0.5_rprec * sqrt(grav/K_wave(jx,jy)) *         &
                      (1.0_rprec/K_wave(jx,jy))

        S_omega_2D(1,1) = 0.0_rprec
        S_k_2D(1,1) = 0.0_rprec

        if ( filt_spec .and. omega_wave(jx,jy) > omega_filter) then
            S_omega_2D(jx,jy) = 0.0_rprec
            S_k_2D(jx,jy) = 0.0_rprec
        end if


  end do
end do

end subroutine spectrum_calc


end module wave_spectrum
