module planet_mod

    integer Tindex, sigindex
    integer index, index1, index2
    integer ocean
    real*8 Mdotcore, dMcore, Mdotgas
    real*8 totmass, TauKH, Mcrit, Mgiso
    real*8 Mcore, Mpl, massLimit, Mdust
    real*8 Tm, csm, H, rhom
    real*8 rhocore
    real*8 RH, Rc, af
    real*8 adddust
    real*8 finalcoresR(100)
    real*8 finalcoreM(100), finalpM(100)
    real*8 check(100)
    real*8 Fg, omega
    real*8 tolG, tolD, temptolD, temptolG
    real*8 Rcore, rm

end module planet_mod
