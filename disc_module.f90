module disc_module

    implicit none

    integer numpts
    Parameter (numpts = 500000)
    real*8 sigmagas(500000),sigmadust(500000)
    real*8 AddUpGas(100), AddUpEatenDust(100), initialcoresR(100)
    real*8 marker(100,2), rad(500000), temper(500000)
    real*8 RADarray(500), SIGarray(1000,500), OSA(1000,500)
    real*8 initialMD(100)
    real*8 Mgas(100)

end module disc_module
