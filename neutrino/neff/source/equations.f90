module massive_neutrinos

    use precision
    use constants
    
    implicit none
    private 
    
    real(dl), parameter  :: pi = const_pi 
    real(dl), parameter  :: const  = 7._dl/120*pi**4 ! 5.68219698_dl
    !const = int q^3 F(q) dq = 7/120*pi^4
    real(dl), parameter  :: const2 = 5._dl/7/pi**2   !0.072372274_dl
    real(dl), parameter  :: zeta3  = 1.2020569031595942853997_dl
    real(dl), parameter  :: zeta5  = 1.0369277551433699263313_dl
    real(dl), parameter  :: zeta7  = 1.0083492773819228268397_dl

    integer, parameter  :: nrhopn=2000  
    real(dl), parameter :: am_min = 0.01_dl  !0.02_dl
    !smallest a*m_nu to integrate distribution function rather than using series
    real(dl), parameter :: am_max = 600._dl 
    !max a*m_nu to integrate
          
    real(dl),parameter  :: am_minp=am_min*1.1
    real(dl), parameter :: am_maxp=am_max*0.9
   
    real(dl) dlnam, mass_nu

    real(dl), dimension(:), allocatable ::  r1,p1,dr1,dp1,ddr1

    !Sample for massive neutrino momentum
    !These settings appear to be OK for P_k accuate at 1e-3 level
    integer, parameter :: nqmax0=80 !maximum array size of q momentum samples 
    real(dl) :: nu_q(nqmax0), nu_int_kernel(nqmax0)
 
    integer nqmax !actual number of q modes evolves
 
    public massive_nu_init, mass_nu, nu_density

    contains

    subroutine massive_nu_init(omnu, nu_massive, H0)
      
!  Initialize interpolation tables for massive neutrinos.
!  Use cubic splines interpolation of log rhonu and pnu vs. log a*m.
        
        real(dl), intent(in) :: omnu, nu_massive, H0
        integer i
        real(dl) dq,dlfdlq, q, am, rhonu,pnu
        real(dl) grhom_in, grhog_in, grhor_in
        real(dl) spline_data(nrhopn)
     
!  nu_masses=m_nu(i)*c**2/(k_B*T_nu0).
!  Get number density n of neutrinos from
!  rho_massless/n = int q^3/(1+e^q) / int q^2/(1+e^q)=7/180 pi^4/Zeta(3)
!  then m = Omega_nu/N_nu rho_crit /n
!  Error due to velocity < 1e-5
        
        !do i=1, CP%Nu_mass_eigenstates 
        !    nu_masses(i) = const/(1.5d0*zeta3)*grhom/grhor*CP%omegan*CP%Nu_mass_fractions(i) &
        !                   /CP%Nu_mass_degeneracies(i)
        !end do
        
        grhom_in = 3*H0**2/c**2*1000**2
        grhog_in = kappa/c**2*4*sigma_boltz/c**3*2.725_dl**4*Mpc**2
        grhor_in = 7.0_dl/8.0_dl*(4.0_dl/11.0_dl)**(4.0_dl/3.0_dl)*grhog_in
        mass_nu = const/(1.5d0*zeta3)*grhom_in/grhor_in * omnu / nu_massive

        if (allocated(r1)) return
        allocate(r1(nrhopn),p1(nrhopn),dr1(nrhopn),dp1(nrhopn),ddr1(nrhopn))

        nqmax=3
        !if (AccuracyLevel >1) nqmax=4
        !if (AccuracyLevel >2) nqmax=5
        !if (AccuracyLevel >3) nqmax=nint(AccuracyLevel*10) 
          !note this may well be worse than the 5 optimized points

        if (nqmax > nqmax0) call MpiStop('massive_nu_init: qmax > nqmax0')

        !We evolve evolve 4F_l/dlfdlq(i), so kernel includes dlfdlnq factor
        !Integration scheme gets (Fermi-Dirac thing)*q^n exact,for n=-4, -2..2
        !see CAMB notes
        if (nqmax==3) then
          !Accurate at 2e-4 level
          nu_q(1:3) = (/0.913201, 3.37517, 7.79184/)
          nu_int_kernel(1:3) = (/0.0687359, 3.31435, 2.29911/)
          
        else if (nqmax==4) then
          !This seems to be very accurate (limited by other numerics)
           nu_q(1:4) = (/0.7, 2.62814, 5.90428, 12.0/)
           nu_int_kernel(1:4) = (/0.0200251, 1.84539, 3.52736, 0.289427/)
      
        else if (nqmax==5) then
        !exact for n=-4,-2..3 
        !This seems to be very accurate (limited by other numerics)
         nu_q(1:5) = (/0.583165, 2.0, 4.0, 7.26582, 13.0/)  
         nu_int_kernel(1:5) = (/0.0081201, 0.689407, 2.8063, 2.05156, 0.126817/) 
  
        else
         dq = (12 + nqmax/5)/real(nqmax)
         do i=1,nqmax
            q=(i-0.5d0)*dq
            nu_q(i) = q 
            dlfdlq=-q/(1._dl+exp(-q))
            nu_int_kernel(i)=dq*q**3/(exp(q)+1._dl) * (-0.25_dl*dlfdlq) !now evolve 4F_l/dlfdlq(i)
            
         end do
        end if
        nu_int_kernel=nu_int_kernel/const
        
        dlnam=-(log(am_min/am_max))/(nrhopn-1)
 

        !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(STATIC) &
        !$OMP & PRIVATE(am, rhonu,pnu) 
        do i=1,nrhopn
          am=am_min*exp((i-1)*dlnam)
          call nuRhoPres(am,rhonu,pnu)
          r1(i)=log(rhonu)
          p1(i)=log(pnu)
        end do
        !$OMP END PARALLEL DO


        call splini(spline_data,nrhopn)
        call splder(r1,dr1,nrhopn,spline_data)
        call splder(p1,dp1,nrhopn,spline_data)
        call splder(dr1,ddr1,nrhopn,spline_data)

    end subroutine massive_nu_init

    subroutine nuRhoPres(am,rhonu,pnu)
!  Compute the density and pressure of one eigenstate of massive neutrinos,
!  in units of the mean density of one flavor of massless neutrinos.

        real(dl),  parameter :: qmax=30._dl
        integer, parameter :: nq=100
        real(dl) dum1(nq+1),dum2(nq+1)
        real(dl), intent(in) :: am
        real(dl), intent(out) ::  rhonu,pnu
        integer i
        real(dl) q,aq,v,aqdn,adq
      

!  q is the comoving momentum in units of k_B*T_nu0/c.
!  Integrate up to qmax and then use asymptotic expansion for remainder.
        adq=qmax/nq
        dum1(1)=0._dl
        dum2(1)=0._dl
        do  i=1,nq
          q=i*adq
          aq=am/q
          v=1._dl/sqrt(1._dl+aq*aq)
          aqdn=adq*q*q*q/(exp(q)+1._dl)
          dum1(i+1)=aqdn/v
          dum2(i+1)=aqdn*v
        end do
        call splint(dum1,rhonu,nq+1)
        call splint(dum2,pnu,nq+1)
!  Apply asymptotic corrrection for q>qmax and normalize by relativistic
!  energy density.
        rhonu=(rhonu+dum1(nq+1)/adq)/const
        pnu=(pnu+dum2(nq+1)/adq)/const/3._dl   
    end subroutine nuRhoPres

    subroutine nu_density(am,rhonu)
        use precision
        !use ModelParams
        real(dl), intent(in) :: am
        real(dl), intent(out) :: rhonu

!  Compute massive neutrino density in units of the mean
!  density of one eigenstate of massless neutrinos.  Use cubic splines to
!  interpolate from a table.

        real(dl) d
        integer i
      
        if (am <= am_minp) then
          rhonu=1._dl + const2*am**2  
          return
        else if (am >= am_maxp) then
          rhonu = 3/(2*const)*(zeta3*am + (15*zeta5)/2/am)
          return
        end if
        
        d=log(am/am_min)/dlnam+1._dl
        i=int(d)
        d=d-i
       
!  Cubic spline interpolation.
        rhonu=r1(i)+d*(dr1(i)+d*(3._dl*(r1(i+1)-r1(i))-2._dl*dr1(i) &
               -dr1(i+1)+d*(dr1(i)+dr1(i+1)+2._dl*(r1(i)-r1(i+1)))))
        rhonu=exp(rhonu)
    end subroutine nu_density

end module massive_neutrinos


module model_density

    real(8)     :: w_lam = -1.00d0
    real(8)     :: grhom, grhog, grhor, grhoc, grhon, grhob, grhov, grhok
    real(8)     :: grhonu_massless, grhonu_massive
	real(8)     :: N0, athom

    public grhom, grhog, grhor, grhoc, grhon, grhob, grhov, grhok, grhonu_massless, &
    & grhonu_massive, N0, athom, w_lam

    contains

    function deta_da(a)

        use massive_neutrinos
    
    	real(8),    intent(in)  :: a
    	real(8)     :: deta_da
    	real(8)     :: a2, grhoa2
        real(8)     :: rhonu
    
        a2 = a**2.d0
        grhoa2 = grhok*a2 + (grhoc+grhob)*a + grhog + grhonu_massless
        if (w_lam .eq. -1.00d0) then
            grhoa2 = grhoa2 + grhov*a2**2
        else
            grhoa2 = grhoa2 + grhov*a**(1.0d0-3.0d0*w_lam)
        end if
    
        if (grhonu_massive .ne. 0) then
            call nu_density(a*mass_nu, rhonu)
            grhoa2 = grhoa2 + rhonu*grhonu_massive
        endif
        
        deta_da = sqrt(3.d0/grhoa2)
        
        return
    end function deta_da

end module model_density

