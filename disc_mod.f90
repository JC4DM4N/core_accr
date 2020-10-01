module disc_mod

    implicit none

    integer numpts
    Parameter (numpts = 500000)

    integer nok,nob,nml,nex
    integer Tary, Rary
    real*8 sigmagas(500000),sigmadust(500000)
    real*8 AddUpGas(100), AddUpEatenDust(100), initialcoresR(100)
    real*8 marker(100,2), rad(500000), temper(500000)
    real*8 RADarray(500), Tarray(1000), SIGarray(1000,500), OSA(1000,500)
    real*8 initialMD(100)
    real*8 Mgas(100)
    real*8 sigma5au, signorm
    real*8 sigmadust_temp
    real*8 Mstar
    real*8 Tfrac, SIGfrac, gasSig
    real*8 aIce, etaIce
    real*8 rIn, rOut, delR

end module disc_mod
