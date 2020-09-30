module migration_mod

    ! Variables controlling migration
    integer :: type1 = 0 ! 0=off, 1=on
    real*8  :: type1factor = 10.0d0 ! tmigI = C*tmigI
    integer :: type2 = 0 ! 0=off, 1=on
    integer :: randomMig = 0 ! 0=off, 1=on
    integer :: sigorb = 2 ! Used for random migration calculations
    integer :: type1dir
    real*8  :: type1dt, type1time, type2dt, type2time
    real*8  :: drbdt, tauMig, taudisc, mu, alpha
    real*8  :: betaMig, type1drbdt, type1max30, dlogSdlogr
    real*8  :: vel,dist,norb(100),norbit,ranpm,rannos(100000)

end module migration_mod
