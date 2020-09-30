module everything_mod
    integer l, i, j, Tary, Rary, done
    integer Tindex, sigindex
    integer numstar,numlhr,numok
    integer index, index1, index2, numpl, numdisc
    integer ocean, nok,nob,nml,nex
    integer iprint
    real*8 Fg, omega
    real*8 Tm, csm, H, rhom
    real*8 rhocore
    real*8 RH, Rc, af
    real*8 tolG, tolD, temptolD, temptolG
    real*8 aIce, etaIce
    real*8 finalcoresR(100)
    real*8 finalMD(100), finalcoreM(100), finalpM(100)
    real*8 check(100)
    real*8 Rcore, rIn, rOut, delR, rm
    real*8 time, dt, origdt
    real*8 disclife
    real*8 Mstar, Mcore, Mpl, massLimit, Mdust
    real*8 Mdotcore, dMcore, Mdotgas
    real*8 fgas, adddust
    real*8 totmass
    real*8 TauKH, Mcrit, Mgiso
    real*8 sigmadust_temp
    real*8 Tarray(1000)
    real*8 Tfrac, SIGfrac, gasSig, scaleM, scaleF
    real*8 sigma5au, signorm
    real rand
end module everything_mod
