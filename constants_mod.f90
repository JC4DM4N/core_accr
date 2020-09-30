module constants_mod

    ! Define constants
    real*8 :: k = 1.3807d-16
    real*8 :: pi = DACos(-1.0d0)
    real*8 :: AU = 1.496d13 ! Astronomical unit, cm
    real*8 :: G = 6.672d-8 ! Gravitational constant, cm^3 g^-1 s^-2
    real*8 :: stefan = 5.670d-5 ! boltzmann constant
    real*8 :: Msun = 1.989d33 ! Solar mass, g
    real*8 :: Me = 5.979d27 ! Earth mass, g
    real*8 :: Mj = 317.8*5.979d27 ! Jupiter mass, g
    real*8 :: MuH = 2.0d0 ! Mean hydrogren mass, g
    real*8 :: muHe = 4.0d0 ! Mean helium mass, g
    real*8 :: mH = 1.673d-24 ! Hydrogen mass
    real*8 :: beta = 1.0d0 ! Planetesimal surface density exponent
    real*8 :: yr = 365.25D0*24.0D0*60.0D0*60.0D0 ! s
    real*8 :: kappa = 1.0d0

    COMMON /mol/ muH, muHe, mH, XX, YY
    COMMON /const/ AU, G, pi, k, Msun, stefan

end module constants_mod
