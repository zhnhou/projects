program thetastar

    use bbn
    use expansion
    use ModelParams
    !use extension,  only: getArgument, nArguments

    implicit none

    character(*),   parameter   :: path = '/home/zhenhou/data/projects/sptxhfi/params/thetastar/chain'
    character(*),   parameter   :: prefix1 = 'test'
    character(*),   parameter   :: prefix2 = 'test_thetastar'

    !character(256)  :: arg, path, prefix1, prefix2
    !integer(4)      :: nArgs

    character(256)  :: pname_file, chain_file_in, chain_file_out

    character(1)    :: cidx
    character(2)    :: pidx
    character(10)   :: pname

    integer(4)      :: num_params
    integer(4)      :: file_unit, unit1, unit2
    real(8)         :: zstar, rs, Da, thetas

    real(8),    allocatable :: param_in_line(:), param_out_line(:)

    !nArgs = nArguments()

    !if (nArgs .ne. 3) then
    !    write(*,*) 'Usage: thetastar [path prefix1 prefix2]'
    !    stop
    !endif

    !call getArgument(1, arg)
    !path = trim(adjustl(arg))

    !call getArgument(1, arg)
    !prefix1 = trim(adjustl(arg))

    !call getArgument(1, arg)
    !prefix2 = trim(adjustl(arg))

    num_params = 0
    pname_file = trim(path)//'/'//trim(prefix1)//'.paramnames'
    open(newunit=file_unit, file=trim(pname_file), status='old', action='read')
    do while (.not. eof(file_unit))
        read(file_unit,'(A10)') pname
        num_params = num_params+1
    enddo
    close(file_unit)

    allocate(param_in_line(0:num_params-1))
    allocate(param_out_line(0:num_params))

    chain_file_in = trim(path)//'/'//trim(prefix1)//'.txt'
    chain_file_out = trim(path)//'/'//trim(prefix2)//'.txt'
    open(newunit=unit1, file=trim(chain_file_in), status='old', action='read')
    open(newunit=unit2, file=trim(chain_file_out), status='unknown', action='write')
    do while (.not. eof(unit1))
        read(unit1,'(100E18.9)') param_in_line

        params%ombh2  = param_in_line(2)
        params%omch2  = param_in_line(3)
        params%omnuh2 = param_in_line(4)
        params%omk    = 0.00d0
        params%neff   = 3.04600d0
        params%Yp     = yp_bbn(params%ombh2, params%neff-3.0460d0)
        params%H0     = param_in_line(7)

        call SetTEParams(params%H0, .true.)
        call z_recomb(zstar)

        thetas = Theta_Sound(zstar, rs, Da)

        param_out_line(0:num_params-1) = param_in_line(0:num_params-1)
        param_out_line(num_params) = thetas

        write(unit2, '(100E18.9)') param_out_line
    enddo

    close(unit1)
    close(unit2)

    deallocate(param_in_line, param_out_line)
end program
