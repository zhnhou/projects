module expansion

    use constants
    use precision
    !use settings
    use Recombination
    use RecData
    use model_density
    use massive_neutrinos

    !use cmbtypes

    implicit none

    private
    
    real(8),    parameter   :: PI = const_pi
    real(8),    parameter   :: Tcmb = 2.725d0
    real(8),    parameter   :: OmegakFlat = 5.00d-7

    !! parameters relevant to thermo-expansion *AFTER* BBN
    type ThermExp_Params
		real(8)     :: ombh2, omb
		real(8)     :: omch2, omc
        real(8)     :: omnuh2, omnu
        real(8)     :: ommh2, omm
		real(8)     :: omk
		real(8)     :: omv, w
		real(8)     :: H0
		real(8)     :: neff
		real(8)     :: Yp
        real(8)     :: nu_massive
        real(8)     :: nu_massless
        real(8)     :: m_nu

        logical     :: flat, closed, open
        real(8)     :: curv,r, Ksign !CP%r = 1/sqrt(|CP%curv|), CP%Ksign = 1,0 or -1
        real(8)     :: tau0,chi0 !time today and rofChi(CP%tau0/CP%r)
    end	type ThermExp_Params

    type(ThermExp_Params)               :: Params
    type(RecombinationParams),  private :: Recomb
    
    !public SetTEParams
    public ThermExp_Params, Params, z_equality, zeq_simple, z_recomb, Theta_Sound, Theta_Damp, &
           ThetaSToH0, ThetaDToYp, ThetaDToOmbh2, ZeqToOmch2, SetTEParams, z_baryondrag, &
           DvtoRs, DvtoRs_zdrag, D_v, Rs_drag, sum_mnu_ev, ThetaSToOmnuh2

    contains

    subroutine SetTEParams(H0, first)
        real(8),    intent(in)      :: H0
        logical,    intent(in)      :: first

        Params%H0    = H0

        w_lam = Params%w

        Params%ommh2 = Params%ombh2 + Params%omch2 + Params%omnuh2
        Params%omm   = Params%ommh2/(Params%H0/100.d0)**2.0d0
        Params%omb   = Params%ombh2/(Params%H0/100.d0)**2.0d0
        Params%omc   = Params%omch2/(Params%H0/100.d0)**2.0d0
        Params%omnu  = Params%omnuh2/(Params%H0/100.d0)**2.0d0
        Params%omv   = 1.00d0 - Params%omb - Params%omc - Params%omnu - Params%omk

        if (Params%omnu .eq. 0.0_dl) then
            Params%nu_massless = Params%neff
            Params%nu_massive  = 0.0_dl
        else
            Params%nu_massive  = Params%neff
            Params%nu_massless = 0.0_dl
        endif

        Params%flat = (abs(Params%omk) <= OmegakFlat)
        Params%closed = Params%omk < -OmegakFlat

        Params%open = .not. Params%flat .and. .not. Params%closed
        if (Params%flat) then
            Params%curv=0
            Params%Ksign=0
            Params%r=1._dl !so we can use tau/CP%r, etc, where CP%r's cancel
        else   
            Params%curv = -Params%omk/((c/1000)/Params%H0)**2
            Params%Ksign = sign(1._dl,Params%curv)
            Params%r = 1._dl/sqrt(abs(Params%curv))
        end if


        grhom = 3*Params%H0**2/c**2*1000**2
        grhog = kappa/c**2*4*sigma_boltz/c**3*Tcmb**4*Mpc**2
        grhor = 7.d0/8.d0*(4.d0/11.d0)**(4.d0/3.d0)*grhog
        grhoc = grhom*Params%omc
        grhob = grhom*Params%omb
        grhon = grhom*Params%omnu
        grhov = grhom*Params%omv
        grhok = grhom*Params%omk

        grhonu_massless = Params%nu_massless * grhor
        grhonu_massive  = Params%nu_massive * grhor

        mass_nu = 0.0_dl
        if (grhonu_massive .ne. 0) then
            call massive_nu_init(Params%omnu, Params%nu_massive, Params%H0)
        endif
        Params%m_nu = 1.68e-4*mass_nu  !! eV mass of each species
        
        if (first) then
            call Recombination_SetDefParams(Recomb)
            N0 = Params%omb * (1.00d0 - Params%Yp) * grhom*c**2/kappa/m_H/Mpc**2
            athom = sigma_thomson*N0*Mpc 
        endif

        return
    end subroutine SetTEParams


    function rofChi(Chi) !sinh(chi) for open, sin(chi) for closed.
        real(dl) Chi,rofChi

        if (Params%closed) then
            rofChi=sin(chi)
        else if (Params%open) then
            rofChi=sinh(chi)
        else
            rofChi=chi
        endif

        return
    end function rofChi


    function sum_mnu_ev(nu_massive)
        
        real(8),    intent(in)  :: nu_massive
        real(8)                 :: sum_mnu_ev

        sum_mnu_ev = nu_massive * Params%m_nu
        
        return
    end function sum_mnu_ev


    subroutine ZeqToOmch2(z_eq, ombh2, neff, omch2)
        
        real(8),    intent(in)  :: z_eq
        real(8),    intent(in)  :: ombh2
        real(8),    intent(in)  :: neff
        real(8),    intent(out) :: omch2

        real(8)     :: f_nu

        f_nu = neff*(7.d0/8.d0)*(4.d0/11.d0)**(4.d0/3.d0)
        omch2 = (z_eq+1.0d0)*(8.d0*Pi*G*a_rad*(1.d0+f_nu)*Tcmb**4)/(3.d0*(c*bigH)**2) - ombh2

        return
    end subroutine ZeqToOmch2


    subroutine ThetaDToYp(ThetaD, ThetaS, zstar)

        real(8),    intent(in)  :: ThetaD
        real(8),    intent(in)  :: ThetaS
        real(8),    intent(inout)   :: zstar
        
        integer     :: i
        real(8)     :: DA
        real(8)     :: D_b, D_t, D_try, try_b, try_t, lasttry
        
        DA = ThetaD
        try_b = 0.01d0
        Params%Yp = try_b
        call ThetaSToH0(ThetaS, zstar)
        D_b = Theta_Damp(zstar)
        
        try_t = 0.70d0
        Params%Yp = try_t
        call ThetaSToH0(ThetaS, zstar)
        D_t = Theta_Damp(zstar)
        
        if (DA < D_b .or. DA > D_t) then
            write(*,*) DA, D_b, D_t
            write(*,*) "CANNOT find Yp, STOP!" 
            stop
        else
            lasttry = -1
            i = 0
            do
                i = i + 1
                if (i .eq. 100) stop 'Yp finder did not converge'
                Params%Yp = (try_b+try_t)/2.d0
                call ThetaSToH0(ThetaS, zstar)
                call z_recomb(zstar)
                D_try = Theta_Damp(zstar)
        
                if (D_try < DA) then
                    try_b = (try_b+try_t)/2
                else
                    try_t = (try_b+try_t)/2
                end if

                if (abs(D_try - lasttry) .lt. 1e-8) then
                    exit
                endif
                lasttry = D_try
            end do
        endif
        
        return
    end subroutine ThetaDToYp


    subroutine ThetaDToOmbh2(ThetaD, ThetaS, zstar, ommh2)

        use bbn

        real(8),    intent(in)  :: ThetaD
        real(8),    intent(in)  :: ThetaS
        real(8),    intent(inout)   :: zstar
        real(8),    intent(in), optional :: ommh2
        
        integer     :: i
        real(8)     :: DA
        real(8)     :: D_b, D_t, D_try, try_b, try_t, lasttry
        
        DA = ThetaD
        try_b = 1.5d-2
        Params%ombh2 = try_b
        if (present(ommh2)) Params%omch2 = ommh2 - Params%ombh2
        Params%Yp    = yp_bbn(Params%ombh2, Params%neff-3.0460d0)

        call ThetaSToH0(ThetaS, zstar)
        D_b = Theta_Damp(zstar)
        
        try_t = 3.0d-2
        Params%ombh2 = try_t
        if (present(ommh2)) Params%omch2 = ommh2 - Params%ombh2
        Params%Yp    = yp_bbn(Params%ombh2, Params%neff-3.0460d0)
        call ThetaSToH0(ThetaS, zstar)
        D_t = Theta_Damp(zstar)
        
        if (DA > D_b .or. DA < D_t) then
            write(*,*) DA, D_b, D_t
            write(*,*) "CANNOT find Omegabh2, STOP!" 
            stop
        else
            lasttry = -1
            i = 0
            do
                i = i + 1
                if (i .eq. 100) stop 'Omegabh2 finder did not converge'
                Params%ombh2 = (try_b+try_t)/2.d0
                if (present(ommh2)) Params%omch2 = ommh2 - Params%ombh2
                Params%Yp    = yp_bbn(Params%ombh2, Params%neff-3.0460d0)

                call ThetaSToH0(ThetaS, zstar)
                call z_recomb(zstar)
                D_try = Theta_Damp(zstar)
        
                if (D_try > DA) then
                    try_b = (try_b+try_t)/2
                else
                    try_t = (try_b+try_t)/2
                end if

                if (abs(D_try - lasttry) .lt. 1e-8) then
                    exit
                endif
                lasttry = D_try
            end do
        endif
        
        return
    end subroutine ThetaDToOmbh2


    subroutine ThetaSToH0(ThetaS, zstar)

		real(8),    intent(in)  :: ThetaS
		real(8),    intent(inout)   :: zstar
!		real(8),	intent(out)	:: H0

        integer     :: i
		real(8)     :: DA
     	real(8)     :: D_b, D_t, D_try, try_b, try_t, lasttry

        DA = ThetaS
        try_b = 30.0d0
        call SetTEParams(try_b, .true.)
        call z_recomb(zstar)
        D_b = Theta_Sound(zstar)
        
        try_t = 150.d0
        call SetTEParams(try_t, .false.)
        call z_recomb(zstar)
        D_t = Theta_Sound(zstar)
        
        if (DA < D_b .or. DA > D_t) then
            write(*,*) DA, D_b, D_t
            write(*,*) "CANNOT find H_0, STOP!" 
            stop
        else
            lasttry = -1
            i = 0
            do
                i = i + 1
                if (i .eq. 100) stop 'H_0 finder did not converge'
                call SetTEParams((try_b+try_t)/2.0d0, .false.)
                call z_recomb(zstar)
                D_try = Theta_Sound(zstar)
        
                if (D_try < DA) then
                    try_b = (try_b+try_t)/2
                else
                    try_t = (try_b+try_t)/2
                end if
        
                if (abs(D_try - lasttry) .lt. 1e-8) then
                    exit
                endif
                lasttry = D_try
            end do
        endif
        
        return
    end subroutine ThetaSToH0


    subroutine ThetaSToOmnuh2(ThetaS, zstar, H0)

		real(8),    intent(in)  :: ThetaS, H0
		real(8),    intent(inout)   :: zstar
!		real(8),	intent(out)	:: H0

        integer     :: i
		real(8)     :: DA
     	real(8)     :: D_b, D_t, D_try, try_b, try_t, lasttry

        DA = ThetaS
        try_b = 0.0d0
        Params%omnuh2 = try_b
        call SetTEParams(H0, .true.)
        call z_recomb(zstar)
        D_b = Theta_Sound(zstar)
        
        try_t = 0.1d0
        Params%omnuh2 = try_t
        call SetTEParams(H0, .false.)
        call z_recomb(zstar)
        D_t = Theta_Sound(zstar)
        
        if (DA < D_b .or. DA > D_t) then
            write(*,*) DA, D_b, D_t
            write(*,*) "CANNOT find omnuh2, STOP!" 
            stop
        else
            lasttry = -1
            i = 0
            do
                i = i + 1
                if (i .eq. 100) stop 'omnuh2 finder did not converge'
                Params%omnuh2 = (try_b+try_t)/2.0d0
                call SetTEParams(H0, .false.)
                call z_recomb(zstar)
                D_try = Theta_Sound(zstar)
        
                if (D_try < DA) then
                    try_b = (try_b+try_t)/2
                else
                    try_t = (try_b+try_t)/2
                end if
        
                if (abs(D_try - lasttry) .lt. 1e-8) then
                    exit
                endif
                lasttry = D_try
            end do
        endif
        
        return
    end subroutine ThetaSToOmnuh2


    !subroutine z_equality(zEQ)

	!	real(8),    intent(out) :: zEQ

	!	real(8)     :: bigH, H0p, fnu
	!	real(8)     :: foe
	!	real(8)     :: omega_radh2, ommh2
	!	real(8)     :: aEQ

    !    bigH= 100.0D3/Mpc
    !    H0p = Params%H0/100._dl*bigH
    !    fnu = Params%neff*(7.d0/8.d0)*(4.d0/11.d0)**(4.d0/3.d0)

    !    omega_radh2 = (8.d0*Pi*G*a_rad*(1.d0+fnu)*Tcmb**4)/(3.d0*(c*bigH)**2)
    !    ommh2 = Params%ombh2 + Params%omch2
    !    aEQ = omega_radh2/ommh2

    !    zEQ = 1.00d0/aEQ - 1.00d0

    !    return
    !end subroutine z_equality

    subroutine z_equality(zEQ)

		real(8),    intent(out) :: zEQ
        
        real(8)     :: omnuh2
        real(8)     :: grhorad, grhomat
		real(8)     :: aEQ
        real(8)     :: z1, z2, avg, diff, ratio
        integer     :: i

        omnuh2 = Params%omnuh2

        if (omnuh2 .eq. 0) then
            grhorad = grhog + grhonu_massless
            aEQ = grhorad / (grhob + grhoc)
            zEQ = 1.00d0/aEQ - 1.00d0
        else
            z1 = 0.0d0
            z2 = 10000.d0

            i = 0
            diff = 10.0d0
            do while (diff .gt. 1d-8)
                i = i + 1
                if (i .eq. 200) stop 'redshift finder did not converge'

                diff = ratio_radmat_massivenu(z2) - ratio_radmat_massivenu(z1)
                avg = 0.5d0*(z2 + z1)
                ratio = ratio_radmat_massivenu(avg)
                if (ratio .gt. 1.d0) then
                    z2 = avg
                else
                    z1 = avg
                end if
            enddo

            zEQ = avg
        endif

        return
    end subroutine z_equality


    subroutine zeq_simple(zEQ)

		real(8),    intent(out) :: zEQ
        
        real(8)     :: grhorad, grhomat
		real(8)     :: aEQ

        grhorad = grhog + Params%neff * grhor
        aEQ = grhorad / (grhob + grhoc)
        zEQ = 1.00d0/aEQ - 1.00d0

        return
    end subroutine zeq_simple


    function ratio_radmat_massivenu(z)
        real(8),    intent(in)     :: z
        real(8)     :: ratio_radmat_massivenu

        real(8)     :: a, rhonu
        real(8)     :: grhomat, grhorad

        a = 1.0d0/(1.0d0+z)
        grhomat = (grhob + grhoc) * a
        call nu_density(a*mass_nu, rhonu)
        grhorad = rhonu*grhonu_massive + grhog + grhonu_massless

        ratio_radmat_massivenu = grhorad / grhomat

        return
    end function

     
    subroutine z_recomb(zstar)
     
	 	real(8),    intent(out) :: zstar
	 	real(8)     :: try1,try2,diff,avg
         integer     :: i
         
         zstar = 0.0d0
         
         call Recombination_init(Recomb, Params%omc, Params%omb, 0.0d0, Params%omv, &
         & Params%H0, Tcmb, Params%Yp, Params%neff)
!	 	call Recombination_init(Recomb, Pa%OmegaC, Pa%OmegaB, 0.0d0, Pa%OmegaV, Pa%H0, &
!	 	& Tcmb, Pa%YHe)
     
         try1 = 0.d0
         try2 = 10000.d0
     
         i = 0
         diff = 10.d0
         do while (diff .gt. 1d-8)
            i=i+1
            if (i .eq. 100) stop 'optical depth redshift finder did not converge'

            diff = tau(try2) - tau(try1)
            avg = 0.5d0*(try2+try1)
            if (tau(avg) .gt. 1.d0) then
                try2 = avg
            else
                try1 = avg
            end if
        end do

        zstar = avg

        return
    end subroutine z_recomb


    subroutine z_baryondrag(zdrag)

        real(8),    intent(out) :: zdrag
		real(8)     :: try1,try2,diff,avg
        integer     :: i

        real(8)     :: tau_drag, a_drag

        zdrag = 0.0d0
        call Recombination_init(Recomb, Params%omc, Params%omb, 0.0d0, Params%omv, &
        & Params%H0, Tcmb, Params%Yp, Params%neff)

        try1 = 0.d0
        try2 = 5000.d0

        i = 0
        diff = 10.d0
        do while (diff .gt. 1d-6)
            i=i+1
            if (i .eq. 100) stop 'optical depth redshift finder did not converge'

            diff = tau(try2) - tau(try1)
            avg = 0.5d0*(try2+try1)
            
            a_drag = 1.00d0/(1.00d0+avg)
            tau_drag = a_drag * 3.d0/4.d0*Params%omb*grhom/grhog

            if (tau(avg) .gt. tau_drag) then
                try2 = avg
            else
                try1 = avg
            end if
        end do

        zdrag = avg

        return
    end subroutine z_baryondrag


    function dtau_dz(z)

		real(dl) :: dtau_dz
		real(dl), intent(in) :: z
		real(dl) :: a

        a = 1._dl/(1._dl+z)
        !ignoring reionisation, not relevant for distance measures
        dtau_dz = Recombination_xe(a) * athom * deta_da(a)

        return
    end function dtau_dz


    function tau(z)
		real(dl) :: rombint2
		real(dl) tau
		real(dl),intent(in) :: z

        tau = rombint2(dtau_dz, 0.d0, z, 1d-6, 20, 100)

        return
    end function tau


    function dRs_da(a)

		real(8),    intent(in)  :: a
		real(8)     :: dRs_da
		real(8)     :: R, cs

        R = 3*grhob*a / (4*grhog)
        cs = 1.0d0/sqrt(3.0d0*(1.0d0+R))
        dRs_da = deta_da(a) * cs
        
        return
    end function dRs_da

    
    function Theta_Sound(zstar, Rs_out, Da_out)

		real(8),    intent(in)  :: zstar
		real(8),    intent(out),optional    :: Rs_out, Da_out

		real(8)     :: astar, Theta_Sound
		real(8)     :: Rs, Da_raw, Da, rombint

        astar = 1.0d0/(zstar+1.0d0)
        Rs = rombint(dRs_da,1d-8,astar,1d-7)
        Da_raw = rombint(deta_da,astar,1.0d0,1.0d-7)
        Da = Params%r * rofchi(Da_raw / Params%r)

        if (present(Rs_out)) Rs_out = Rs
        if (present(Da_out)) Da_out = Da

        Theta_Sound = Rs/Da

        return
    end function Theta_Sound


    function Theta_Damp(zstar, Rd_out, Da_out)

		real(8),    intent(in)  :: zstar
		real(8),    intent(out),optional    :: Rd_out, Da_out

		real(8)     :: Theta_Damp
		real(8)     :: rombint
		real(8)     :: Rd2, Rd, astar, Da_raw, Da

        astar = 1.0d0/(1.0d0+zstar)

        call Recombination_init(Recomb, Params%omc, Params%omb, 0.0d0, Params%omv, &
        & Params%H0, Tcmb, Params%Yp, Params%neff)

        Rd2 = rombint(dRd2_da, 1d-8, astar,1d-7)
        Rd = Pi*sqrt(Rd2)
        Da_raw = rombint(deta_da,astar,1.0d0,1.0d-7)
        Da = Params%r * rofchi(Da_raw / Params%r)

        if (present(Rd_out)) Rd_out = Rd
        if (present(Da_out)) Da_out = Da

        Theta_Damp = Rd/Da

        return
    end function Theta_Damp


    function dRd2_da(a)

        implicit none

		real(8)     :: a
		real(8)     :: xe
		real(8)     :: R
		real(8)     :: dRd2_da

        R = 3*grhob*a / (4*grhog)

        xe = Recombination_xe(a)
        dRd2_da = deta_da(a)/(xe*athom/a**2.0d0)
        dRd2_da = dRd2_da/(6.d0*(1.d0+R)**2.d0)*(R**2.d0+16.d0/15.d0*(1.d0+R))

        return    
    end function dRd2_da


    function Rs_drag(zdrag)
        real(8)     :: zdrag
        real(8)     :: adrag, Rs_drag
        real(8)     :: rombint

        adrag = 1.0d0/(1.0d0+zdrag)
        Rs_drag = rombint(dRs_da,1.0d-8, adrag, 1d-7)

        return
    end function


    function D_v(z)
        real(8),    intent(in)  :: z

        real(8)     :: D_v
        real(8)     :: a, Da_raw, Da
        real(8)     :: grhoam2, hz
        real(8)     :: am2, am3, am4, rhonu
        real(8)     :: rombint
        
        a = 1.00d0/(1.00d0+z)
        Da_raw = rombint(deta_da,a,1.0d0,1.0d-7)
        Da = Params%r * rofchi(Da_raw / Params%r)

        !grhoa2 = grhok*a2 + (grhoc+grhob)*a + grhog + grhonu_massless
        !if (w_lam .eq. -1.00d0) then
        !    grhoa2 = grhoa2 + grhov*a2**2
        !else
        !    grhoa2 = grhoa2 + grhov*a**(1.0d0-3.0d0*w_lam)
        !end if

        !if (grhonu_massive .ne. 0) then
        !    call nu_density(a*mass_nu, rhonu)
        !    grhoa2 = grhoa2 + rhonu*grhonu_massive
        !endif

        am2 = a**(-2.0d0)
        am3 = a**(-3.0d0)
        am4 = a**(-4.0d0)

        grhoam2 = grhok*am2 + (grhoc+grhob+grhon)*am3 + (grhog + grhonu_massless)*am4
        grhoam2 = grhoam2 + grhov*am3**(1.0d0+w_lam)
        if (grhonu_massive .ne. 0) then
            call nu_density(a*mass_nu, rhonu)
            grhoam2 = grhoam2 + rhonu*grhonu_massive*am4
        endif

        hz = sqrt(grhoam2/3.0d0)

        D_v = (Da**2.d0*z/Hz)**(1.d0/3.d0)
        
        return
    end function D_v


    function DvtoRs_zdrag(z, zdrag)
        real(8)     :: z, zdrag
        real(8)     :: DvtoRs_zdrag

        DvtoRs_zdrag = D_v(z)/Rs_drag(zdrag)

        return
    end function


    function DvtoRs(z, rs)
        real(8)     :: z, rs
        real(8)     :: DvtoRs

        DvtoRs = D_v(z)/rs

        return
    end function


    function Acoustic(z)

        real(8),    intent(in)  :: z
        real(8)     :: Acoustic
        real(8)     :: ckm

        ckm = c/1e3_dl
        Acoustic = 100*D_v(z)*sqrt(Params%ommh2)/(ckm*z)
        
        return
    end function


end module expansion
