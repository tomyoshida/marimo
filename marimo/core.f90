!
! f2py -m core --fcompiler=gfortran --f90flags='-fopenmp' -lgomp -c core.f90
!

module core
        use iso_fortran_env
        !$ use omp_lib
        implicit none


        real(8), parameter :: cc = 29979245800.0
        real(8), parameter :: hp = 6.6260755e-27
        real(8), parameter :: kb = 1.380658e-16
        real(8), parameter :: mp = 1.672621924e-24
        real(8), parameter :: GG = 6.67430e-8
        real(8), parameter :: atm = 1013250.0

        real(8), parameter :: pi = 3.14159265358979323846264338

        
        integer(4), save :: progress = 0
        integer(4), save :: max_prog = 1

        !integer(kind=int32) :: progress
        !integer(kind=int32) :: max_prog

contains


subroutine integrate(Varr, Xarr, Yarr, dX, rmin, subpix_max, subpix_pow, PA, incl, z2r, mu, nu0, Jl, Aul, nzmax, &
        sigma_g0, gamma_g, r0, rexp, rin, rout, rcav, dcav, d2g, k_d, &
        tmid0, pmid, tatm0, patm, mu_gas, mstar, xco, vturb, &
        gamma_h2, r_h2, gamma_he, r_he, tref, use_dust, chi_d, contsub, &
        use_voigt, kind, D, nthreads,  I, nVarr, nXarr, nYarr)
        !f2py threadsafe

        real(8), intent(in) :: PA, incl, z2r, mu, nu0, Jl, Aul, &
        sigma_g0, gamma_g, r0, rexp, rin, rout, rcav, dcav, d2g, k_d, &
        tmid0, pmid, tatm0, patm, mu_gas, mstar, xco, vturb, &
        gamma_h2, r_h2, gamma_he, r_he, tref, chi_d

        integer, intent(in) :: nXarr, nYarr, nVarr, nzmax, use_voigt, kind, nthreads, use_dust, contsub
        real(8), intent(in) :: Xarr(nXarr), Yarr(nYarr), Varr(nVarr), dX, rmin, subpix_max, subpix_pow
        real(8), intent(out) :: I(nYarr, nXarr, nVarr)
        real(8) :: X, Y, dZ, s, I_subpix !, X_pos, Y_pos
        real(8) :: dX_sub, R_tmp
        real(8) :: rr, tt, zz, temp, nh2, nco, sig, gam, dt, dnu, alpha, phi, xv, yv, continuum, dt_dust, sigma_d, sigma_g
        real(8) :: temp2, tau_l, Z, vnow, vkep, Zs, Ze, zz_old, dr_pos, R_pos_old, dt_pos, D
        integer :: iz, iv, ix, iy, flag_cont, ix_sub, iy_sub, Nsubpix
        integer :: proc_num, nproc, Nsubpix_count

        real(8) :: sqrtln2 = 0.8325546111576977d0

        
        call omp_set_num_threads(nthreads)

        progress = 0

        !$omp parallel default(firstprivate), shared(I, progress, max_prog)

        max_prog = nVarr * nXarr * nYarr
        proc_num = omp_get_thread_num()

        nproc = omp_get_num_procs()
        
        !$omp do
        do iv=1, nVarr
                vnow = Varr(iv)

                !R_pos_old = 0.0d0

                do ix=1, nXarr
                
                        !X = Xarr(ix)

                        !dR_pos = R_pos - R_pos_old

                        do iy=1, nYarr
                                !Y = Yarr(iy)

                                !if (proc_num == 0) then
                                progress = progress + 1
                                !end if

                                !! sub-pixelize !!
                                ! call xyz2rtz( Xarr(ix), Yarr(iy), 0.0d0, PA, incl, rr, tt, zz)
                                R_tmp = sqrt( Xarr(ix)**2 +  Yarr(iy)**2 )

                                !if (R_tmp > rmin) then
                                
                                Nsubpix = ceiling( subpix_max * ( R_tmp / rmin )**(-subpix_pow) )
                                ! 1, 2, 3, 4, 5 ...

                                if (Nsubpix > 1) then
                                        dX_sub = dX / (2 * (Nsubpix-1) )
                                        ! dX/2, dX/4, dX/6...
                                else
                                        dX_sub = 0
                                end if

                                I_subpix = 0.0d0

                                Nsubpix_count = 0

                                do ix_sub = 1, 2 * Nsubpix - 1

                                        X = Xarr(ix) + dX_sub * (ix_sub - Nsubpix)
                                        ! Nsubpix = 1 & ix_sub = 1 -> X = Xarr(ix)
                                        ! Nsubpix = 2 & ix_sub = 1 -> X = Xarr(ix) - dX / 2
                                        ! Nsubpix = 2 & ix_sub = 2 -> X = Xarr(ix)
                                        ! Nsubpix = 2 & ix_sub = 3 -> X = Xarr(ix) + dX / 2

                                        ! Nsubpix = 3 & ix_sub = 1 -> X = Xarr(ix) - 2 * dX / 4
                                        ! Nsubpix = 3 & ix_sub = 2 -> X = Xarr(ix) - 1 * dX / 4
                                        ! Nsubpix = 3 & ix_sub = 3 -> X = Xarr(ix)
                                        ! Nsubpix = 3 & ix_sub = 4 -> X = Xarr(ix) + 1 * dX / 4
                                        ! Nsubpix = 3 & ix_sub = 5 -> X = Xarr(ix) + 2 * dX / 4

                                do iy_sub = 1, 2 * Nsubpix - 1

                                        Y = Yarr(iy) + dX_sub * (iy_sub - Nsubpix)

                                        R_tmp = sqrt( X**2 +  Y**2 )

                                        if (R_tmp > rmin) then
                                       
                                        call zse(X, Y, PA, incl, z2r, Zs, Ze)
                                        dZ = (Ze-Zs) / nzmax

                                        s = 0.0d0
                                        Z = Zs
                                        tau_l = 0.0d0

                                        call xyz2rtz(X, Y, Z, PA, incl, rr, tt, zz)

                                        temp2 = f_temp(rr, zz, tmid0, pmid, tatm0, patm, r0, mu_gas, mstar)    

                                        flag_cont = 1 
                                        continuum = 0.0d0

                                        zz_old = 100*Ze

                                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                        !!! integration along the line of sight !!!
                                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                        do iz=1, nzmax
                                                        
                                        call xyz2rtz(X, Y, Z, PA, incl, rr, tt, zz)

                                        temp = f_temp(rr, zz, tmid0, pmid, tatm0, patm, r0, mu_gas, mstar)
                                                        
                                        nh2 = f_nh2(rr, zz, sigma_g0, gamma_g, r0, rexp, rin, rout, &
                                        rcav, dcav, tmid0, pmid, mu_gas, mstar)

                                        nco = nh2 * xco

                                        sig = sqrt(2*kb*temp/(mu*mp) * ( 1.0d0 + vturb**2 * (mu/2.37) ) )

                                        gam = f_voigt_c(temp, gamma_h2, r_h2, gamma_he, r_he, tref) * nh2

                                        vkep = f_vkep(rr, tt, zz, incl, mstar)
                                        dnu = (vnow-vkep)/cc * nu0

                                        alpha = cc / sig / nu0
                                        xv = alpha * dnu
                                        yv = alpha * gam

                                        !xv = sqrtln2 * abs(vnow-vkep)/sig / sqrt(2.0)
                                        !yv = sqrtln2 * alpha * gam
                                                                
                                        if (use_voigt == 1) then

                                                if (kind == 1) then
                                                        !dt = af(nu0, jl, aul, temp) * nco &
                                                        !* humlicek_w4(xv, yv) * (alpha / sqrt(pi) ) * dZ
                                                        dt = af(nu0, jl, aul, temp) * nco &
                                                        * voigt_humlicek(xv, yv) * (alpha / sqrt(pi) ) * dZ
                                                else if (kind == 2) then
                                                        dt = af(nu0, jl, aul, temp) * nco &
                                                        * cpf12(xv, yv) * (alpha/ sqrt(pi) ) * dZ
                                                end if
                                        else
                                                dt = af(nu0, jl, aul, temp) * nco &
                                                * exp(-(vnow-vkep)**2/sig**2) / (nu0/cc * sig * sqrt(pi)) * dZ
                                                !dt = af_old(nu0, jl, aul, mu, temp) * nco &
                                                !* exp(-(vnow-vkep)**2/sig**2) * dZ

                                        end if

                                        if ( dt < 1.0e-40 ) then
                                                dt = 1.0e-40
                                        end if
                                                
                                        if (temp > 2000) then
                                               temp = 2000
                                        end if

                                        s = s * exp(-dt) &
                                        + (1.0-dt*exp(-dt)-exp(-dt)) / dt * planck(nu0, temp) &
                                        + (dt-1.0+exp(-dt))/dt * planck(nu0, temp2)

                                        ! when the runner reaches the midplane for the first time,
                                        if ( zz >= 0.0d0 .and. flag_cont == 1 .and. rr > rin .and. use_dust == 1) then

                                                ! at the midplane
                                                if (rr > rin .and. rr < rout) then
                                                        sigma_g = sigma_g0 * (rr/r0)**(-gamma_g) * exp( -( rr/rexp )**(2.0-gamma_g))
                
                                                        if ( rr < rcav ) then
                                                                sigma_g = sigma_g * dcav
                                                        end if
                                                else
                                                        sigma_g = 0.0
                                                end if

                                                sigma_d = sigma_g * d2g
                                                dt_dust = sigma_d * k_d

                                                continuum = chi_d * planck(nu0, temp2) * ( 1 - exp(-dt_dust) )
                                                s = s * exp(-dt_dust) + continuum
                                                
                                                flag_cont = 0
                                        end if

                                        zz_old = zz
                                        Z = Z + dZ
                                        tau_l = tau_l + dt
                                        temp2 = temp

                                        end do


                                        if (contsub == 1 .and. use_dust == 1) then
                                                I_subpix = I_subpix + ( s - continuum )
                                        else
                                                I_subpix = I_subpix + s
                                        end if  

                                        Nsubpix_count = Nsubpix_count + 1

                                        end if
                                                
                                end do
                                        
                                ! sub-pix finished
                                end do 


                                if (Nsubpix_count > 0) then                   
                                        I(iy, ix, iv) = I_subpix * 1.0d0 / Nsubpix_count
                                        ! I(iy, ix, iv) = tau_l / Nsubpix_count
                                else
                                        I(iy, ix, iv) = 0.0d0
                                end if

                                !else
                                
                                !I(iy, ix, iv) = 0.0d0

                                !end if

                        end do

                end do

        end do
        !$omp end do
        !$omp end parallel

end subroutine integrate


subroutine xyz2rtz(X, Y, Z, PA, incl, rr, tt, zz)

        real(8), intent(in):: X, Y, Z, PA, incl
        real(8) xi, yi, zi, xd, yd, zd
        real(8), intent(out):: rr, tt, zz
        
        ! 1. image (X, Y, Z) to image (xi, yi, zi)
        xi = cos(PA) * X - sin(PA) * Y
        yi = sin(PA) * X + cos(PA) * Y
        zi = Z

        ! 2. image (xi, yi, zi) to disk (xd, yd, zd)
        ! considering rotation around x-axis
        xd = xi
        yd = cos(-incl) * yi - sin(-incl) * zi
        zd = sin(-incl) * yi + cos(-incl) * zi

        ! disk (xd, yd, zd) to disk (r, t, z)
        rr = sqrt(xd**2 + yd**2)
        tt = atan2(yd, xd)
        zz = zd

end subroutine xyz2rtz

subroutine zse(x, y, pa, incl, z2r, zs, ze)

        real(8), intent(in) :: x, y, pa, incl, z2r
        real(8), intent(out) :: zs, ze
        real(8) :: xim, yim, b, a1, a2, a3, sol1, sol2

        ! find starting/ending points
        ! image (X, Y) to image (xim, yim)
        xim = cos(pa) * x - sin(pa) * y
        yim = sin(pa) * x + cos(pa) * y

        ! the disk boundary corn z=f(r, t)
        b = z2r / sqrt( 1.0 + z2r**2 )
        a1 = cos(incl)**2 - b**2
        a2 = - sin(2.0*incl) * yim
        a3 = yim**2 * sin(incl)**2 - (xim**2 + yim**2)*b**2

        sol1 = ( -a2 + sqrt( a2**2 - 4.0*a1*a3 ) )/2.0/a1
        sol2 = ( -a2 - sqrt( a2**2 - 4.0*a1*a3 ) )/2.0/a1

        if (sol1 < sol2) then
                zs = sol1
                ze = sol2
        else
                zs = sol2
                ze = sol1
        end if

end subroutine zse


real(8) function au2rad(l, D)

        real(8), intent(in) :: l, D
        real(8) :: au = 14959787070000.0d0

        au2rad = l / au / D / 3600.0d0 * 0.0174533

        return
end function au2rad

real(8) function planck(nu, temp)

        real(8), intent(in) :: nu, temp

        planck = 2.0 * hp * nu**3 / cc**2 / ( exp(hp*nu/kb/temp) - 1 )

        return
end function planck

real(8) function af_old(nu0, jl, aul, mu, t)

        real(8), intent(in) :: nu0, jl, aul, mu, t
        real(8) :: m, a1, a2, a3, a4
        
        m = mu * mp

        a1 = cc**3 / ( 8.0 * pi * sqrt(2.0*pi) * nu0**3 ) * (jl + 3.0/2.0) / (jl+1.0)
        a2 = aul * hp * nu0 / (kb*t) * sqrt( m/(kb*t) )
        a3 = ( 1.0 - exp(-hp*nu0/kb/t) )
        a4 = exp(-hp*nu0/2.0/kb/t * ( jl+1.0/(3*(jl+1)) ) )

        af_old = a1 * a2 * a3 *a4
        return 
end function af_old

real(8) function af(nu0, jl, aul, t)

        real(8), intent(in) :: nu0, jl, aul, t
        real(8) :: ju, gu, gl, nl, El, B0, U

        ju = jl + 1

        gl = 2 * jl + 1
        gu = 2 * ju + 1

        B0 = nu0 / 2.0d0 / ju
        El = hp * B0 * jl * ( jl + 1 )

        U = kb * t / hp / B0 * exp( hp * B0 / 3.0 / kb/ t )
        
        nl = gl * exp( -El/kb/t ) / U

        af = cc**2 / ( 8.0 * pi * nu0**2) * gu/gl * aul * nl * (1 - exp(-hp*nu0/kb/t))

        return 
end function af

real(8) function f_vkep(r, t, z, incl, mstar)

        real(8), intent(in) :: r, t, z, incl, mstar

        f_vkep = sqrt( GG*mstar * r**2/( (r**2 + z**2)**(1.5d0)) ) * sin(incl) * cos(t)

end function f_vkep

real(8) function f_voigt_c(t, gamma_h2, r_h2, gamma_he, r_he, tref)

        real(8), intent(in) :: t, gamma_h2, r_h2, gamma_he, r_he, tref

        !tref = 296
        !gamma_h2 = 0.0711
        !r_h2 = 0.66
        !gamma_he = 0.0472
        !r_he = 0.51
        !atm = 1013250

        f_voigt_c = ( gamma_h2 * (tref/t)**r_H2 + 0.19 * gamma_He * (tref/t)**r_He ) * cc*kb*t/atm

end function f_voigt_c


real(8) function f_hg(r, tmid, mu_gas, mstar)

        real(8), intent(in) :: r, tmid, mu_gas, mstar
        real(8) :: cs, omega_k

        cs = sqrt(kb*tmid/(mu_gas*mp))
        omega_k = sqrt( GG * mstar / r**3 )

        f_hg = cs / omega_k

end function f_hg

real(8) function f_temp(r, z, tmid0, pmid, tatm0, patm, r0, mu_gas, mstar)

        real(8), intent(in) :: r, z, tmid0, pmid, tatm0, patm, r0, mu_gas, mstar
        real(8) :: tmid, tatm, hg, l1, l2

        tmid = tmid0 * ( r/r0 ) ** (-pmid)
        tatm = tatm0 * ( r/r0 ) ** (-patm)

        hg = f_hg( r, tmid, mu_gas, mstar)

        f_temp = 3.0d0

        l1 = 1.0d0
        l2 = 4.0d0

        if (abs(z) < l1 *hg) then
                f_temp = tmid 
        else if ( abs(z) < l2 * hg ) then
                f_temp = tmid + ( tatm - tmid ) * (abs(z) - l1*hg) / (l2-l1)/hg
        else 
                f_temp = tatm
        end if

        return  
end function f_temp
                

real(8) function f_nh2(r, z, sigma_g0, gamma_g, r0, rexp, rin, rout, rcav, dcav, tmid0, pmid, mu_gas, mstar)

        real(8), intent(in) :: r, z, sigma_g0, gamma_g, r0, rexp, rin, rout, rcav, dcav, tmid0, pmid, mu_gas, mstar
        real(8) :: tmid, sigma_g, hg

        tmid = tmid0 * ( r/r0 ) ** (-pmid)
        hg =  f_hg( r, tmid, mu_gas, mstar)
        
        if (r > rin .and. r < rout) then
                sigma_g = sigma_g0 * (r/r0)**(-gamma_g) * exp( -( r/rexp )**(2.0-gamma_g))
                
                if ( r < rcav ) then
                        sigma_g = sigma_g * dcav
                end if
        else
                sigma_g = 0.0
        end if

        f_nh2 = sigma_g / ( mu_gas * mp * sqrt(2*pi) * hg ) * exp( -0.5*(z/hg)**2 )

end function f_nh2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                                                                !!
real(8) function cpf12(x, y)
!!                                                                                                                                !!
!!    complex probability function for complex argument Z=X+iY                                                                    !!
!!    real part = voigt function K(x,y)                                                                                           !!
!!                                                                                                                                !!
!!    Source:   J. Humlicek, JQSRT 21, pp. 309-313, 1979                                                                          !!
!!              An efficient method for evaluation of the complex probability function: the Voigt function and its derivatives    !!
!!                                                                                                                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT none
	
    DOUBLE PRECISION, PARAMETER ::  zero  = 0.d0
    DOUBLE PRECISION, PARAMETER ::  one   = 1.d0
    DOUBLE PRECISION, PARAMETER ::  two   = 2.d0
    

    INTEGER, PARAMETER ::   kdp=KIND(zero)

    DOUBLE PRECISION, INTENT(IN)  ::   x, y
    COMPLEX(KIND=kdp) ::  prbFct


!   local constants
    DOUBLE PRECISION, PARAMETER :: y0=1.5d0, y0sq=2.25d0, bigExp=170.d0
    DOUBLE PRECISION, PARAMETER :: t(6)     = [ .3142403762544,  .9477883912402, 1.5976826351526, &
                                                 2.2795070805011, 3.0206370251209, 3.88972489786978 ]
    DOUBLE PRECISION, PARAMETER :: alpha(6) = [ -1.393236997981977,    -0.2311524061886763, +0.1553514656420944, &
                                                -0.006218366236965554, -9.190829861057117E-5, +6.275259577E-7 ]
    DOUBLE PRECISION, PARAMETER :: beta (6) = [ 1.011728045548831,    -0.7519714696746353, 0.01255772699323164, &
                                                0.01002200814515897,  -2.420681348155727E-4,  5.008480613664576E-7 ]
!   local scratch arrays and some further variables
    DOUBLE PRECISION ::            wr, wi
    DOUBLE PRECISION ::            yy0,yy0sq, y0yy0, y2y0

!   ================================================================================================================================

!   initialize
    
    prbFct = zero
    
    yy0   = y+y0
    yy0sq = yy0*yy0
    y0yy0 = y0 * (y+y0)
    y2y0  = y + two*y0

!   --------------------------------------------------------------------------------------------------------------------------------
    IF (y>0.85) THEN

!       region I only
        
        wr = + (alpha(1)*(x-t(1)) +  beta(1)*yy0) / ((x-t(1))**2+yy0sq) &
                 - (alpha(1)*(x+t(1)) -  beta(1)*yy0) / ((x+t(1))**2+yy0sq) &
                 + (alpha(2)*(x-t(2)) +  beta(2)*yy0) / ((x-t(2))**2+yy0sq) &
                 - (alpha(2)*(x+t(2)) -  beta(2)*yy0) / ((x+t(2))**2+yy0sq) &
                 + (alpha(3)*(x-t(3)) +  beta(3)*yy0) / ((x-t(3))**2+yy0sq) &
                 - (alpha(3)*(x+t(3)) -  beta(3)*yy0) / ((x+t(3))**2+yy0sq) &
                 + (alpha(4)*(x-t(4)) +  beta(4)*yy0) / ((x-t(4))**2+yy0sq) &
                 - (alpha(4)*(x+t(4)) -  beta(4)*yy0) / ((x+t(4))**2+yy0sq) &
                 + (alpha(5)*(x-t(5)) +  beta(5)*yy0) / ((x-t(5))**2+yy0sq) &
                 - (alpha(5)*(x+t(5)) -  beta(5)*yy0) / ((x+t(5))**2+yy0sq) &
                 + (alpha(6)*(x-t(6)) +  beta(6)*yy0) / ((x-t(6))**2+yy0sq) &
                 - (alpha(6)*(x+t(6)) -  beta(6)*yy0) / ((x+t(6))**2+yy0sq)
        wi = +  (beta(1)*(x-t(1)) - alpha(1)*yy0) / ((x-t(1))**2+yy0sq) &
                 +  (beta(1)*(x+t(1)) + alpha(1)*yy0) / ((x+t(1))**2+yy0sq) &
                 +  (beta(2)*(x-t(2)) - alpha(2)*yy0) / ((x-t(2))**2+yy0sq) &
                 +  (beta(2)*(x+t(2)) + alpha(2)*yy0) / ((x+t(2))**2+yy0sq) &
                 +  (beta(3)*(x-t(3)) - alpha(3)*yy0) / ((x-t(3))**2+yy0sq) &
                 +  (beta(3)*(x+t(3)) + alpha(3)*yy0) / ((x+t(3))**2+yy0sq) &
                 +  (beta(4)*(x-t(4)) - alpha(4)*yy0) / ((x-t(4))**2+yy0sq) &
                 +  (beta(4)*(x+t(4)) + alpha(4)*yy0) / ((x+t(4))**2+yy0sq) &
                 +  (beta(5)*(x-t(5)) - alpha(5)*yy0) / ((x-t(5))**2+yy0sq) &
                 +  (beta(5)*(x+t(5)) + alpha(5)*yy0) / ((x+t(5))**2+yy0sq) &
                 +  (beta(5)*(x-t(5)) - alpha(5)*yy0) / ((x-t(5))**2+yy0sq) &
                 +  (beta(5)*(x+t(5)) + alpha(5)*yy0) / ((x+t(5))**2+yy0sq)
        prbFct = CMPLX(wr, wi, kdp)

!       ---------------------------------------------------------------------------------------------------------------------------
    ELSE
!       ---------------------------------------------------------------------------------------------------------------------------
        
        IF (abs(x)<18.1*y+1.65) THEN
!             - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!             region I
              wr = + (alpha(1)*(x-t(1)) +  beta(1)*yy0) / ((x-t(1))**2+yy0sq) &
                   - (alpha(1)*(x+t(1)) -  beta(1)*yy0) / ((x+t(1))**2+yy0sq) &
                   + (alpha(2)*(x-t(2)) +  beta(2)*yy0) / ((x-t(2))**2+yy0sq) &
                   - (alpha(2)*(x+t(2)) -  beta(2)*yy0) / ((x+t(2))**2+yy0sq) &
                   + (alpha(3)*(x-t(3)) +  beta(3)*yy0) / ((x-t(3))**2+yy0sq) &
                   - (alpha(3)*(x+t(3)) -  beta(3)*yy0) / ((x+t(3))**2+yy0sq) &
                   + (alpha(4)*(x-t(4)) +  beta(4)*yy0) / ((x-t(4))**2+yy0sq) &
                   - (alpha(4)*(x+t(4)) -  beta(4)*yy0) / ((x+t(4))**2+yy0sq) &
                   + (alpha(5)*(x-t(5)) +  beta(5)*yy0) / ((x-t(5))**2+yy0sq) &
                   - (alpha(5)*(x+t(5)) -  beta(5)*yy0) / ((x+t(5))**2+yy0sq) &
                   + (alpha(6)*(x-t(6)) +  beta(6)*yy0) / ((x-t(6))**2+yy0sq) &
                   - (alpha(6)*(x+t(6)) -  beta(6)*yy0) / ((x+t(6))**2+yy0sq)
              wi = +  (beta(1)*(x-t(1)) - alpha(1)*yy0) / ((x-t(1))**2+yy0sq) &
                   +  (beta(1)*(x+t(1)) + alpha(1)*yy0) / ((x+t(1))**2+yy0sq) &
                   +  (beta(2)*(x-t(2)) - alpha(2)*yy0) / ((x-t(2))**2+yy0sq) &
                   +  (beta(2)*(x+t(2)) + alpha(2)*yy0) / ((x+t(2))**2+yy0sq) &
                   +  (beta(3)*(x-t(3)) - alpha(3)*yy0) / ((x-t(3))**2+yy0sq) &
                   +  (beta(3)*(x+t(3)) + alpha(3)*yy0) / ((x+t(3))**2+yy0sq) &
                   +  (beta(4)*(x-t(4)) - alpha(4)*yy0) / ((x-t(4))**2+yy0sq) &
                   +  (beta(4)*(x+t(4)) + alpha(4)*yy0) / ((x+t(4))**2+yy0sq) &
                   +  (beta(5)*(x-t(5)) - alpha(5)*yy0) / ((x-t(5))**2+yy0sq) &
                   +  (beta(5)*(x+t(5)) + alpha(5)*yy0) / ((x+t(5))**2+yy0sq) &
                   +  (beta(5)*(x-t(5)) - alpha(5)*yy0) / ((x-t(5))**2+yy0sq) &
                   +  (beta(5)*(x+t(5)) + alpha(5)*yy0) / ((x+t(5))**2+yy0sq)
              prbFct = CMPLX(wr, wi, kdp)
!             - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ELSE
!             - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!             region II
              wr = +( beta(1)*((x-t(1))**2-y0yy0) -alpha(1)*(x-t(1))*y2y0 ) / (((x-t(1))**2+y0sq)*((x-t(1))**2+yy0sq)) &
                   +( beta(1)*((x+t(1))**2-y0yy0) +alpha(1)*(x+t(1))*y2y0 ) / (((x+t(1))**2+y0sq)*((x+t(1))**2+yy0sq)) &
                   +( beta(2)*((x-t(2))**2-y0yy0) -alpha(2)*(x-t(1))*y2y0 ) / (((x-t(2))**2+y0sq)*((x-t(2))**2+yy0sq)) &
                   +( beta(2)*((x+t(2))**2-y0yy0) +alpha(2)*(x+t(1))*y2y0 ) / (((x+t(2))**2+y0sq)*((x+t(2))**2+yy0sq)) &
                   +( beta(3)*((x-t(3))**2-y0yy0) -alpha(3)*(x-t(1))*y2y0 ) / (((x-t(3))**2+y0sq)*((x-t(3))**2+yy0sq)) &
                   +( beta(3)*((x+t(3))**2-y0yy0) +alpha(3)*(x+t(1))*y2y0 ) / (((x+t(3))**2+y0sq)*((x+t(3))**2+yy0sq)) &
                   +( beta(4)*((x-t(4))**2-y0yy0) -alpha(4)*(x-t(1))*y2y0 ) / (((x-t(4))**2+y0sq)*((x-t(4))**2+yy0sq)) &
                   +( beta(4)*((x+t(4))**2-y0yy0) +alpha(4)*(x+t(1))*y2y0 ) / (((x+t(4))**2+y0sq)*((x+t(4))**2+yy0sq)) &
                   +( beta(5)*((x-t(5))**2-y0yy0) -alpha(5)*(x-t(1))*y2y0 ) / (((x-t(5))**2+y0sq)*((x-t(5))**2+yy0sq)) &
                   +( beta(5)*((x+t(5))**2-y0yy0) +alpha(5)*(x+t(1))*y2y0 ) / (((x+t(5))**2+y0sq)*((x+t(5))**2+yy0sq)) &
                   +( beta(6)*((x-t(6))**2-y0yy0) -alpha(6)*(x-t(1))*y2y0 ) / (((x-t(6))**2+y0sq)*((x-t(6))**2+yy0sq)) &
                   +( beta(6)*((x+t(6))**2-y0yy0) +alpha(6)*(x+t(1))*y2y0 ) / (((x+t(6))**2+y0sq)*((x+t(6))**2+yy0sq))
              wi = + (beta(1)*(x-t(1)) - alpha(1)*yy0) / ((x-t(1))**2+yy0sq) &
                   + (beta(1)*(x+t(1)) + alpha(1)*yy0) / ((x+t(1))**2+yy0sq) &
                   + (beta(2)*(x-t(2)) - alpha(2)*yy0) / ((x-t(2))**2+yy0sq) &
                   + (beta(2)*(x+t(2)) + alpha(2)*yy0) / ((x+t(2))**2+yy0sq) &
                   + (beta(3)*(x-t(3)) - alpha(3)*yy0) / ((x-t(3))**2+yy0sq) &
                   + (beta(3)*(x+t(3)) + alpha(3)*yy0) / ((x+t(3))**2+yy0sq) &
                   + (beta(4)*(x-t(4)) - alpha(4)*yy0) / ((x-t(4))**2+yy0sq) &
                   + (beta(4)*(x+t(4)) + alpha(4)*yy0) / ((x+t(4))**2+yy0sq) &
                   + (beta(5)*(x-t(5)) - alpha(5)*yy0) / ((x-t(5))**2+yy0sq) &
                   + (beta(5)*(x+t(5)) + alpha(5)*yy0) / ((x+t(5))**2+yy0sq) &
                   + (beta(6)*(x-t(6)) - alpha(6)*yy0) / ((x-t(6))**2+yy0sq) &
                   + (beta(6)*(x+t(6)) + alpha(6)*yy0) / ((x+t(6))**2+yy0sq)
              prbFct = CMPLX(y*wr+exp(-DMIN1(x*x,bigExp)), wi, kdp)
!             - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          END IF
!       ---------------------------------------------------------------------------------------------------------------------------
    END IF

    cpf12 = real(prbFct)

    return
!   ================================================================================================================================

END function cpf12



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                                                                !!
real(8) function zpf16  (x, y)
!!                                                                                                                                !!
!!    complex probability function for complex argument Z=X+iY                                                                    !!
!!    real part = voigt function K(x,y)                                                                                           !!
!!                                                                                                                                !!
!!  large x+y     J  Humlicek, JQSRT 27, 437-444, 1982   asymptotic region I approximation                                        !!
!!  small x+y:    J  Humlicek, JQSRT 21, 309-313, 1979   16-term rational approximation Eq. (6)                                   !!
!!                                                                                                                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IMPLICIT none

    INTEGER, PARAMETER ::   kdp=KIND(1.D0)

    DOUBLE PRECISION, INTENT(IN)  ::   x, y
    COMPLEX(KIND=kdp) ::  w

!   Humlicek cpf16 constants
    COMPLEX(KIND=kdp), PARAMETER ::  a(16) =          &
                [(41445.0374210222, 0.0),   &
                 (0.0, -136631.072925829),  &
                 (-191726.143960199, 0.0),  &
                 (0.0, 268628.568621291),   &
                 (173247.907201704,  0.0),  &
                 (0.0, -179862.56759178),   &
                 (-63310.0020563537, 0.0),  &
                 (0.0, 56893.7798630723),   &
                 (11256.4939105413,  0.0),  &
                 (0.0, -9362.62673144278),  &
                 (-1018.67334277366, 0.0),  &
                 (0.0, 810.629101627698),   &
                 (44.5707404545965, 0.0),   &
                 (0.0, -34.5401929182016),  &
                 (-0.740120821385939, 0.0), &
                 (0.0, 0.564189583547714)]
    DOUBLE PRECISION, PARAMETER ::  b(8) =  &
               [ 7918.06640624997,         &
                 -126689.0625,        &
                 295607.8125,        &
                 -236486.25,        &
                 84459.375,        &
                 -15015.0,        &
                 1365.0,        &
                 -60.0]

!   further variables
    COMPLEX(KIND=kdp) ::              z, zz, numer, denom

!   ================================================================================================================================
    
    z     = CMPLX(x, y+1.31183,kdp)
    zz    = z*z
    numer = ((((((a(16)*z+a(15))*z+a(14))*z+a(13))*z+a(12))*z+a(11))*z+a(10))*z+a(9)
    numer = ((((((((numer*z+a(8))*z+a(7))*z+a(6))*z+a(5))*z+a(4))*z+a(3))*z+a(2))*z+a(1)) 
    denom = b(1)+(b(2)+(b(3)+(b(4)+(b(5)+(b(6)+(b(7)+b(8)*zz)*zz)*zz)*zz)*zz)*zz)*zz
    w  = numer/denom
    
    zpf16 = real(w)

    return
!   ================================================================================================================================

END function zpf16



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                              
!!                                                                   
!!
real(8) function humlicek_w4(x, y)
!!                                                                                                                                !!
!!    complex probability function for complex argument Z=X+iY                                                                    !!
!!    real part = voigt function K(x,y)                                                                                           !!
!!                                                                                                                                !!
!!    source:   J. Humlicek, JQSRT 27, 437, 1982                                                                                  !!
!!                                                                                                                                !!
!!    the stated accuracy is claimed to be 1.0E-04 by the author.                                                                 !!
!!    R.H.Norton has checked the accuracy by comparing values computed using a program written by B.H.Armstrong,                  !!
!!    and the accuracy claim seems to be warranted.                                                                               !!
!!                                                                                             12/91, converted to f90 may 2009   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    IMPLICIT none

    DOUBLE PRECISION, PARAMETER ::  half  = 0.5d0
    DOUBLE PRECISION, PARAMETER ::  recSqrtPi   = 1.0d0 / 1.7724538509055159d0

    INTEGER, PARAMETER ::   kdp=KIND(half)

    DOUBLE PRECISION, INTENT(IN)  ::   x, y
    COMPLEX(KIND=kdp) ::  prbFct !?! vgtK(0:NX),vgtL(0:NX)

    DOUBLE PRECISION ::     s, ax
    COMPLEX(kdp)::   t, u
!   humlicek prbFct region bounds
    DOUBLE PRECISION, PARAMETER ::    s15=15.d0,  s55=5.5d0

!   statement functions
    COMPLEX(kdp) ::   approx1, approx2, approx3, approx4
    approx1(t)   = t*recSqrtPi / (half + t*t)
    approx2(t,u) = (t * (1.410474 + u*recSqrtPi))/ (.75 + (u *(3.+u)))
    approx3(t)   = ( 16.4955 + t * (20.20933 + t * (11.96482 + t * (3.778987 + 0.5642236*t)))) &
                   / ( 16.4955 + t * (38.82363 + t * (39.27121 + t * (21.69274 + t * (6.699398 + t)))))
    approx4(t,u) = (t * (36183.31 - u * (3321.99 - u * (1540.787 - u *(219.031 - u *(35.7668 - u *(1.320522 - u *recSqrtPi)))))) &
                   / (32066.6 - u * (24322.8 - u * (9022.23 - u * (2186.18 - u * (364.219 - u * (61.5704 - u * (1.84144 - u))))))))

!   ================================================================================================================================

    
    t  = CMPLX(y,-x,kdp)
    ax = ABS(x)
    s  = ax + y

    IF (s>=s15) THEN
!          region I
           prbFct= approx1(t)
    ELSE IF (s<s15 .AND. s>=s55) THEN
!          region II
           u = t * t
           prbFct= approx2(t,t*t)
    ELSE IF (s<s55 .AND. y>=(0.195*ax-0.176)) THEN
!         region III
          prbFct= approx3(t)
    ELSE
!         region IV
          u  = t * t
          prbFct= EXP(u) - approx4(t,u)
    END IF
    
    
    humlicek_w4 = real(prbFct)

    return
!   ================================================================================================================================

END function humlicek_w4


!--------------------------------------------------------------------------
!              VOIGT FUNCTION APPROXIMATION BY HUMLICEK
!
! This function is a modified version of a subroutine taken from the 
! following web site: 
!
!    http://www.op.dlr.de/oe/ir/voigt.html
!
! and written by F. Schreier (DLR - Institute for Optoelectronics
! Oberpfaffenhofen, 82230 Wessling) from the publication
!
!  F. Schreier, J. Quant. Spectros. Radiat. Transfer 48, 743-762 (1992)
!
! which bases the algorithm on the paper by Humlicek:
!
!  J. Humlicek, J. Quant. Spectros. Radiat. Transfer 27, 437 (1982)
!
! It calculates the real part of the complex probability function for
! complex argument Z=X+iY. The real part = voigt function K(x,y).
!
! Parameters:                                                      
!   X      X argument                                           in  
!   Y      Voigt function parameter, ratio of lorentz/doppler   in  
!
! Result:
!   function value = Voigt function
!  
! The stated accuracy is claimed to be 1.0E-04 by the author.  R. H. Norton
! has checked the accuracy by comparing values computed using a program
! written by B. H. Armstrong, and the accuracy claim seems to be warranted.
!
! NOTE 1: The normalization of this function is such that 
!
!  /
!  | voigt(x) dx = sqrt(pi)
!  /
!
! NOTE 2: The x coordinate is defined such that the Gaussian component
!         is always exp(-x^2). The parameter y simply changes the 
!         width of the Lorenz component compared to the Gaussian.
!
! fgs 12/91
! Modified: CPD 2011
!--------------------------------------------------------------------------

real(8) function voigt_humlicek(x, y)

  real(8), intent(in) :: x, y
  real(8) :: s, ax
  COMPLEX*16 :: t, u, z
  COMPLEX*16 :: APPROX1, APPROX2, APPROX3, APPROX4
  !
  ! Four different approximations, as inline formulae
  !
  APPROX1(T)   = (T * .5641896) / (.5 + (T * T))
  APPROX2(T,U) = (T * (1.410474 + U*.5641896))/ (.75 + (U *(3.+U)))
  APPROX3(T)   = ( 16.4955 + T * (20.20933 + T * (11.96482 + &
                 T * (3.778987 + 0.5642236*T)))) &
                 / ( 16.4955 + T * (38.82363 + T * &
                 (39.27121 + T * (21.69274 + T * (6.699398 + T)))))
  APPROX4(T,U) = (T * (36183.31 - U * (3321.99 - U * (1540.787 - U &
                 *(219.031 - U *(35.7668 - U *(1.320522 - U * .56419)))))) &
                 / (32066.6 - U * (24322.8 - U * (9022.23 - U * (2186.18 &
                 - U * (364.219 - U * (61.5704 - U * (1.84144 - U))))))))
  !
  ! If y=0, then we have a Gaussian profile, so we are done quickly
  !
  if(y.eq.0.d0) then
     voigt_humlicek = exp(-x**2)
     return
  endif
  !
  ! Voigt Profile...
  !
  ! Distinguish different y-regimes
  !
  if(y.gt.15.d0) then
     !
     ! All points are in region I
     !
     t = cmplx(y,-x)
     z = APPROX1(t)
     !
  elseif(y.ge.5.5d0) then
     !
     ! Points are in region I or region II
     !
     t  = cmplx(y,-x)
     s  = abs(x) + y
     if(s.ge.15.d0) then
        z = APPROX1(t)
     else
        u = t * t
        z = APPROX2(t,u)
     endif
     !
  elseif(y.gt.0.75d0) then
     t  = cmplx(y,-x)
     s  = abs(x) + y
     if(s.ge.15.d0) then
        z = APPROX1(t)
     elseif(s.lt.5.5d0) then
        z = APPROX3(t)
     else
        u = t * t
        z = APPROX2(t,u)
     endif
  else
     t  = cmplx(y,-x)
     ax = abs(x)
     s  = ax + y
     if(s.ge.15.d0) then
        !
        ! region I
        !
        z = APPROX1(t)
     elseif(s.ge.5.5d0) then
        !
        ! region II
        !
        u = t * t
        z = APPROX2(t,u)
     elseif(y.ge.(0.195d0*ax-0.176d0)) then
        !
        ! region III
        !
        z = APPROX3(t)
     else
        !
        ! region IV
        !
        u = t * t
        z = cdexp(u) - APPROX4(t,u)
     endif
  endif
  !
  voigt_humlicek = real(z)

  return
end function voigt_humlicek



end module core



