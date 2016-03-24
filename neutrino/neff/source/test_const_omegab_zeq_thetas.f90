PROGRAM get_samples
    
    use bbn
    use expansion
    use model_density
    use ModelParams

    implicit none
    
    character(*),   parameter   :: prefix = 'samples/params_base_camspec_lowl_lowlike_const_omegab_zeq_thetas'
    character(*),   parameter   :: params_file = prefix//'.txt'
    character(*),   parameter   :: pnames_file  = prefix//'.paramnames'
    
    real(8),        parameter   :: ombh2_best  = 0.022032d0
    real(8),        parameter   :: omch2_best  = 0.12038d0
    real(8),        parameter   :: neff_best   = 3.046d0
    real(8),        parameter   :: H0_best     = 67.04d0
    real(8),        parameter   :: omnuh2_best = 0.000d0
    !real(8),        parameter   :: omnuh2_best = 0.060d0/94.30d0

    real(8),        parameter   :: neff_min = 1.0d0
    real(8),        parameter   :: neff_max = 7.0d0
    real(8),        parameter   :: delta_neff = 1.0d0
    
    real(8),        parameter   :: tau_best    = 0.09250d0
    real(8),        parameter   :: omk_best    = 0.0000000E+00
    real(8),        parameter   :: ns_best     = 0.96190d0
    real(8),        parameter   :: nt_best     = 0.0000000E+00
    real(8),        parameter   :: nrun_best   = 0.0000000E+00
    real(8),        parameter   :: logA_best   = 3.0980d0
    real(8),        parameter   :: r_best      = 0.0000000E+00
    real(8),        parameter   :: w_best      = -1.000000E+00

    character(20)   :: pnames(10)
    
    integer(4)  :: i, nsamples, iz
    real(8)     :: thetas_best, thetad_best, As_best, zeq_best, Yp_best
    real(8)     :: omch2, neff, H0, Yp, zstar, theta_s, theta_d, zeq, omnuh2, rs, Da
    real(8)     :: params_list(10), weight, loglike
    
    call paramnames(pnames)
    open(unit=10, file=pnames_file, status='unknown')
    do i=1, 10
        write(10,*) trim(pnames(i))//'    '//trim(pnames(i))
    enddo 
    close(10)
    
    weight = 1.00d0
    loglike = 0.00d0

    As_best = exp(logA_best)/10.0d0

    Params%ombh2 = ombh2_best
    Params%omch2 = omch2_best
    Params%omnuh2= omnuh2_best
    Params%omk   = omk_best
    Params%neff  = neff_best
    Params%Yp    = yp_bbn(ombh2_best, neff_best-3.0460d0)
    Params%H0    = H0_best
    Params%w     = w_best

    call SetTEParams(H0_best, .true.)
    call z_recomb(zstar)

    thetas_best = Theta_Sound(zstar, rs, Da)
    thetad_best = Theta_Damp(zstar)
    call zeq_simple(zeq_best)
    
    nsamples = (neff_max - neff_min)/delta_neff + 1

    params_list(1) = Params%ombh2
    params_list(2) = Params%omch2
    params_list(3) = Params%omnuh2
    params_list(4) = Params%omk
    params_list(5) = Params%Neff
    params_list(6) = Params%Yp
    params_list(7) = Params%H0
    params_list(8) = tau_best
    params_list(9) = ns_best
    params_list(10)= As_best

    write(*,100), neff_best, omch2_best, zeq_best, thetas_best, thetad_best

    open(unit=10, file=params_file, status='unknown')
    !write(10,102) weight, loglike, params_list

    do i=0, nsamples-1
        neff = neff_min + i*delta_neff

        call zeqtoomch2(zeq_best, ombh2_best, neff, omch2)

        Params%ombh2 = ombh2_best
        Params%omch2 = omch2
        Params%omnuh2= omnuh2_best
        Params%omk   = omk_best
        Params%neff  = neff
        Params%Yp    = yp_bbn(ombh2_best, neff-3.0460d0)
        
        call ThetaSToH0(thetas_best, zstar)
        !call ThetaDToYp(thetad_best, thetas_best, zstar)
    
        theta_s = Theta_Sound(zstar)
        theta_d = Theta_Damp(zstar)

        params_list(1) = Params%ombh2
        params_list(2) = Params%omch2
        params_list(3) = Params%omnuh2
        params_list(4) = Params%omk
        params_list(5) = Params%Neff
        params_list(6) = Params%Yp
        params_list(7) = Params%H0
        params_list(8) = tau_best
        params_list(9) = ns_best
        params_list(10)= As_best

        call zeq_simple(zeq)
        
        write(10,102) weight, loglike, params_list
        write(*,100), neff, Params%omch2, zeq, theta_s, theta_d
    enddo

    close(10)

100 format (100E16.7)
102 format (E16.7,E16.7,100E16.7)


END PROGRAM get_samples


subroutine paramnames(pnames)

    character(20)   :: pnames(10)

    pnames(1) = 'omegabh2'
    pnames(2) = 'omegach2'
    pnames(3) = 'omeganuh2*'
    pnames(4) = 'omegak'
    pnames(5) = 'nnu'
    pnames(6) = 'yheused*'
    pnames(7) = 'H0*'
    pnames(8) = 'tau'
    pnames(9) = 'ns'
    pnames(10)= 'A*'

    return

end subroutine paramnames
