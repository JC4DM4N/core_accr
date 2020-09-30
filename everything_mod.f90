module everything_mod
    integer l, i, j, Tary, Rary, done
    integer Tindex, sigindex
    integer numstar,numlhr,numok
    integer index, index1, index2, numpl, numdisc
    integer ocean, nok,nob,nml,nex
    integer type1, type2, randomMig, sigorb, type1dir
    integer iprint
    real*8 AU, G, pi, k, stefan
    real*8 Fg, omega, beta
    real*8 Tm, csm, H, rhom
    real*8 rhocore, kappa
    real*8 RH, Rc, af
    real*8 drbdt, tauMig, taudisc, mu, alpha
    real*8 betaMig, type1drbdt, type1max30, dlogSdlogr
    real*8 tolG, tolD, temptolD, temptolG
    real*8 aIce, etaIce
    real*8 finalcoresR(100)
    real*8 finalMD(100), finalcoreM(100), finalpM(100)
    real*8 check(100)
    real*8 Rcore, rIn, rOut, delR, rm
    real*8 time, dt, type1dt, type1time, type2dt, type2time, origdt,yr
    real*8 type1factor
    real*8 disclife
    real*8 Me, Msun, Mstar, Mcore, Mpl, massLimit, Mj, Mdust
    real*8 Mdotcore, dMcore, Mdotgas
    real*8 fgas, adddust
    real*8 totmass
    real*8 TauKH, Mcrit, Mgiso
    real*8 sigmadust_temp
    real*8 muH, muHe, mH, XX, YY
    real*8 vel,dist,norb(100),norbit,ranpm,rannos(100000)
    real*8 Tarray(1000)
    real*8 Tfrac, SIGfrac, gasSig, scaleM, scaleF
    real*8 sigma5au, signorm
    real rand
    COMMON /mol/ muH, muHe, mH, XX, YY
    COMMON /const/ AU, G, pi, k, Msun, stefan
end module everything_mod
