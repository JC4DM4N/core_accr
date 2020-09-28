program p1

implicit none
integer numpts, l, m, i, j, n, Tary, Rary, done
integer Tindex, sigindex, sigIn, sigOut
integer numstar,numlhr,numok
integer index, index1, index2, numpl, numdisc
integer ocean, nok,nob,nml,nex
integer type1, type2, randomMig, sigorb, type1dir, obt
integer iprint
Parameter (numpts = 500000)
real*8 AU, G, pi, k, stefan
real*8 Fg, omega, beta
real*8 Tm, csm, H, rhom
real*8 rho, rhocore, tauL, kappa
real*8 rL, RH, Rc, Ra, af
real*8 drbdt, tauMig, taudisc, mu, alpha
real*8 betaMig, type1drbdt, type1max30, dlogSdlogr
real*8 tolG, tolD, temptolD, temptolG, marker(100,2)
real*8 rad(numpts),sigmagas(numpts),sigmadust(numpts)
real*8 temper(numpts), aIce, etaIce
real*8 initialcoresR(100), finalcoresR(100)
real*8 initialMD(100), finalMD(100), finalcoreM(100), finalpM(100)
real*8 AddUpGas(100), AddUpEatenDust(100), check(100), Mgas(100)
real*8 r, Rcore, rIn, rOut, delR, rm
real*8 time, dt, type1dt, type1time, type2dt, type2time, origdt,yr
real*8 type1factor
real*8 disclife
real*8 Me, Msun, Mstar, Mcore, Mpl, massLimit, Mj, Mdisc, Mdust
real*8 Mdotcore, dMcore, Mdotgas
real*8 fdust, fgas, fd, maxfd, minfd, range, adddust
real*8 totmass
real*8 TauKH, Mcrit, Mgiso
real*8 sigmadust_temp
real*8 muH, muHe, mH, XX, YY
real*8 vel,dist,orb,norb(100),norbit,ranpm,pm,rannos(100000),tdt2
real*8 Tarray(1000),RADarray(500),SIGarray(1000,500),OSA(1000,500)
real*8 Tfrac, SIGfrac, gasSig, scaleM, scaleF
real*8 sigma5au, signorm
real rand, x, y
COMMON /mol/ muH, muHe, mH, XX, YY
COMMON /const/ AU, G, pi, k, Msun, stefan

! open random numbers file
open(unit=25,file='randnums.dat',status='old')
read(25,*) rannos
close(25)

! Define constants
k = 1.3807d-16
pi = DACos(-1.0d0)
AU = 1.496d13 ! Astronomical unit, cm
G = 6.672d-8 ! Gravitational constant, cm^3 g^-1 s^-2
stefan = 5.670d-5 ! boltzmann constant
Msun = 1.989d33 ! Solar mass, g
Me = 5.979d27 ! Earth mass, g
Mj = 317.8*Me ! Jupiter mass, g
MuH = 2.0d0 ! Mean hydrogren mass, g
muHe = 3.971d0 ! Mean helium mass, g
muHe = 4.0d0 ! Mean helium mass, g
mH = 1.673d-24 ! Hydrogen mass
beta = 1.0d0 ! Planetesimal surface density exponent
yr = 365.25D0*24.0D0*60.0D0*60.0D0 ! s

Mstar = 2.4d0*Msun ! Star mass

! Read params
open(unit=19,file='core_accr.params',status='old')
read(19,*) rm ! semi-major axis of planet
read(19,*) Mcore ! Initial core mass of planet
read(19,*) sigma5au ! sigma_p @ 5au, used as scaling factor
close(19)

rm = rm*AU
Mcore = Mcore*Me

! Read in evolution of surface density profile of gas disc
if(Mstar.eq.Msun)then
 open(unit=26,file='discMstar1.dat',status='old')
 print*,'Mstar = 1.0Msun'
else if(Mstar.eq.2.4d0*Msun)then
 open(unit=26,file='discABAur.dat',status='old')
 print*,'Mstar = 2.4Msun'
else if(Mstar.eq.0.5d0*Msun)then
 open(unit=26,file='discMstar05.dat',status='old')
 print*,'Mstar = 0.5Msun'
else if(Mstar.ge.0.29d0*Msun)then
 open(unit=26,file='discMstar03.dat',status='old')
 print*,'Mstar = 0.3Msun'
else if(Mstar.ge.0.09d0*Msun)then
 open(unit=26,file='discMstar01.dat',status='old')
 print*,'Mstar = 0.1Msun'
endif
do Tary=1,1000
 read (26,*) Tarray(Tary)
 do Rary=1,500
    read (26,*) RADarray(Rary), OSA(Tary,Rary)
 enddo
enddo
close(26)

origdt = 10000.0D0*yr
disclife = Tarray(1000)/yr
massLimit = 15.0 !max mass is 15 jupiter masses
scaleM = 0.01d0 !Mgas * this = Mdust
scaleF = 0.41d0

nok = 0 ! total number of planets
nob = 0 ! total number out of bounds
nml = 0 ! total munber over the mass limit

! Variables controlling migration
type1 = 0 ! 0=off, 1=on
type1factor = 10.0d0 ! tmigI = C*tmigI
type2 = 0 ! 0=off, 1=on
randomMig = 0 ! 0=off, 1=on
sigorb = 2 ! Used for random migration calculations

tolD = 0.02 ! limits the dust accretion at each timestep, dmcore (*Me)
tolG = 0.02 ! limits the gas accretion at each timestep, dmgas (*Mj)
numpl = 1 ! number of different radii to try for each disc
numdisc = 1 ! number of stars to simulate

maxfd = 2.492533*scaleM*scaleF*(Mstar/Mj) ! fdust scales with star mass
minfd = 0.1*(Mstar/Msun) ! minimum fdust is 0.1(Mstar/Msun)
fd = maxfd - 0.1*(Mstar/Msun) ! fd = fd/(numdisc-1)

do l = 1, numdisc ! loop over discs
  nex = 0 ! counts number of exits due to negligible mass gain
  ! a loop to set all planets as being initially in the simulation
  do i=1,numpl
    check(i) = 1 !1=ok, 2=rh>disc, 3=mass limit exceeded
  enddo

!*************************************************************************
!                         setup procedures
!*************************************************************************

  ! calculate fdust for the current disc
  fdust = minfd + fd*(l-1)
  if(numdisc.eq.1) fdust = maxfd
  print*,' '
  print*,'Running disc number ',l,'/',numdisc,', dustfrac =',fdust

  kappa = 1.0d0

  ! loop over number of simulated planets
  do m = 1,numpl
    done = 1 ! ensure only 1 print per planet(line 640)
    AddUpGas(m) = 0.0 ! total gas mass of planet
    AddUpEatenDust(m) = 0.0 ! total dust mass of planet

    print*, 'Initial planet radius, ', rm/AU, ' AU'

    ! stuff for migration
    initialcoresR(m)=rm/Au
    marker(m,1) = 0.0d0
    marker(m,2) = 0.0d0

    ! limit on maximum mass acreation allowed at each timestep
    temptolG = tolG*(alog(real(1+rm/Au)))**3.0d0 ! gas
    temptolD = tolD*(alog(real(1+rm/Au)))**3.0d0 ! dust

    rIn =  RADarray(1) ! inner edge of grid
    rOut = RADarray(500) ! outer edge of grid
    delR = (rOut-rIn)/numpts ! divide grid into N=numpts

    rad(1) = rIn

    ! create array of N=numpts radius values
    do i = 2, numpts
       rad(i) = rad(i-1) + delR
    enddo

    ! set the temperature at each radius interval
    do i = 1, numpts-1
       temper(i) = 280.0d0*(1.0d0*AU/rad(i))**0.5d0*(Mstar/Msun)
    enddo

    ! set the snow line
    aIce = 2.7d0*(Mstar/Msun)**2.0d0*AU
    ocean = 0 ! would the planet have a liquid water ocean (0=yes, 1=no)
    if (rm.gt.aIce) then
      ocean = 1
    endif

    ! tracks the total disc mass
    totmass = 0.0d0

    ! scales the gas surface density (not scaled here)
    do i=1,1000
      do j=1,500
        ! SIGarray(i,j) = OSA(i,j)*fdust/maxfd*scaleF
        ! SIGarray(i,j) = OSA(i,j)*fdust*0.09242d0*(Msun/Mstar)
        SIGarray(i,j) = OSA(i,j)
      enddo
    enddo

    ! calculate sigmagas, extrapolated from values in SIGarray
    i = 1 ! iradius for extrapolated values
    do j=1,500 ! iradius for known values
      do while(RADarray(j).ge.rad(i))
        if(i.eq.1)then !do the 1st point seperately (no 0 element)
          sigmagas(i) = SIGarray(1,j)
        else
          SIGfrac = SIGarray(1,j-1) - SIGarray(1,j)
          SIGfrac = SIGfrac/(RADarray(j) - RADarray(j-1))
          sigmagas(i) = SIGarray(1,j-1) - SIGfrac*(rad(i) - RADarray(j-1))
        endif
        if (i.eq.numpts) goto 42
        i = i+1
      enddo
    enddo
42  continue

    ! calculate sigmadust profile
    iprint = 0 ! only want to print dust profile once
    do i = 1, numpts-1
      ! scale factor accounting for ice line
      if (rad(i).lt.aIce) then
        etaIce = 1.0d0
      else
        etaIce = 4.2d0
      endif
      sigmadust(i) = fdust*etaIce*10.0d0*(1.0d0*AU/rad(i))**(beta)
      totmass = totmass + &
                (pi*rad(i+1)**2.0D0-pi*rad(i)**2.0D0)*sigmagas(i) + &
                (pi*rad(i+1)**2.0D0-pi*rad(i)**2.0D0)*sigmadust(i)

      if (rad(i).le.rm .and. rad(i+1).gt.rm .and. iprint.eq.0) then
	    signorm = sigmadust(i)/sigma5au ! scale factor for sigmadust
        print*, ''
        print*, 'SIGMA AT 5AU, GAS:', sigmagas(i), 'DUST:', sigmadust(i)/signorm, 'gcm-2'
        print*, ''
        iprint = 1 ! don't print again
        index = i ! array index where the planet is located
      endif
    enddo
    sigmadust = sigmadust/signorm ! normalise sigmadust array

    ! write gas and dust profiles to a file
    open(12345, file='densities.dat', status='unknown')
    do i=1,numpts-1
      write(12345, *) rad(i), sigmagas(i), sigmadust(i)
    enddo
    close(12345)

    initialMD(m) = 0.0 ! initial dust mass
    do i=1,numpts-1
      initialMD(m)=initialMD(m)+(sigmadust(i)*(pi*rad(i+1)**2-pi*rad(i)**2))
    enddo

    ! determine index in RADarray where planet is located
    do i = 1, 499
      if((RADarray(i).le.rm).and.(RADarray(i+1).gt.rm)) then
        sigindex = i
      endif
    enddo

    Mpl = Mcore
    Mgiso = 50.0d0*(sigmagas(index)/2.4d3)**1.5d0*(rm/AU)**3.0*(Mstar/Msun)**(-0.5d0) ! gas isolation mass

    time = 0.0d0
    dt = origdt  ! keeps track of what dt is origionally set as

    drbdt = 0.0d0  !initial core migration is zero

    Mgas(m) = 0.0d0
    do i=1,numpts-1
       Mgas(m)=Mgas(m)+(sigmagas(i)*(pi*rad(i+1)**2-pi*rad(i)**2))
    enddo

    print*,'Mgas=',Mgas(m)/Msun,' Msol, Mdust=',initialMD(m)/Msun,' Msol'
    print*, 'Current dust-to-gas ratio=', initialMD(m)/Mgas(m)
    print*, 'Initial dust-to-gas ratio=', 2.0d0*initialMD(m)/Mstar ! Assumes star-disc intially had q=0.5
    print*,' '

!*************************************************************************
!      disc setup complete...evolve over the lifetime of the disc
!*************************************************************************
    do while (time .lt. disclife)

      Tm = temper(index) ! sets the planets T to that of the surounding dust
      csm = DSqrt(7.0d0/5.0d0*k*Tm/muH/mH) ! local sound speed
      Omega = DSqrt(G*Mstar/rm**3.0d0) ! local Keplerian frequency
      H = csm/omega ! local scale height
      ! get the current time array element
      do i=1,1000
        if(Tarray(i+1).gt.time*yr) then
          Tindex = i
          exit
        endif
      enddo
      ! calculates how far along to the next element the current time is
      Tfrac = time*yr - Tarray(Tindex)
      Tfrac = Tfrac/(Tarray(Tindex+1) - Tarray(Tindex))
      ! extrapolate the sigma gas array to get the density at this partial timestep
      i = 1
      do j=1,500
        do while(RADarray(j).ge.rad(i))
          if(i.eq.1)then !do  the 1st point seperately (no 0 element)
            sigmagas(i) = (SIGarray(Tindex,j) - SIGarray(Tindex+1,j))*Tfrac
            sigmagas(i) = SIGarray(Tindex,j) - sigmagas(i)
          else
            SIGfrac = SIGarray(Tindex,j-1) - SIGarray(Tindex,j)
            SIGfrac = SIGfrac/(RADarray(j) - RADarray(j-1)) ! last timestep
            gasSig = SIGarray(Tindex+1,j-1) - SIGarray(Tindex+1,j)
            gasSig = SIGfrac/(RADarray(j) - RADarray(j-1)) ! next timestep
            sigmagas(i) = (SIGfrac - gasSig)*Tfrac
            sigmagas(i) = SIGarray(Tindex,j-1) - sigmagas(i)
          endif
          if(sigmagas(i).lt.0.0d0) then
            sigmagas(i) = 0.0d0 ! prevents unphysical negative sigmagas
          endif
          if(i.eq.numpts) goto 43
          i = i+1
        enddo
      enddo
43    continue

      ! checks the planet hasnt exceeded its isolation mass or mass limit
      if ((Mpl/Me.lt.Mgiso).and.(Mpl/Mj.lt.massLimit)) then
        rhom = gasSig/H ! gas density at planet location
        rhocore = 3.2d0 ! density of core, gcm^-3
        Rcore = (3.0d0*Mcore/(4.0d0*pi*rhocore))**(1.0d0/3.0d0)
        Rc = Rcore

        RH = rm*(Mpl/3.0d0/Mstar)**(1.0d0/3.0d0) ! calculates the hill radius

        call calcFg(Fg,Rc,RH,af,omega) !calculates the grav enhancement factor

        ! set the inner and outer indices of the planet's hill sphere
        do j = 1, numpts-1
          if ((rad(j).le.(rm-af)).and.(rad(j+1).gt.(rm-af))) then
            index1 = j
          endif
          if ((rad(j).lt.(rm+af)).and.(rad(j+1).ge.(rm+af))) then
            index2 = j+1
          endif
        enddo

        ! condition to end simulation and discard planet for hill radius outwith disc
        if((rm-af).lt.rIn.or.(rm+af).gt.rOut) then
          check(m) = 2
          goto 17
        endif

        ! sum dust mass within the Hill sphere
        Mdust = 0.0d0
        adddust = 0.0d0
        ! inner edge of the Hill sphere
        adddust =(pi*rad(index1+1)**2.0d0-pi*(rm-af)**2.0d0)*sigmadust(index1)
        Mdust = Mdust + adddust
        ! loop over hill sphere
        do j = index1+1, index2-2
          adddust = 0.0d0
          adddust =(pi*rad(j+1)**2.0D0-pi*rad(j)**2.0D0)*sigmadust(j)
          Mdust = Mdust + adddust
        enddo
        ! outer edge of hill sphere
        adddust = 0.0d0
        adddust =(pi*(rm+af)**2.0D0-pi*rad(index2-1)**2.0D0)*sigmadust(index2-1)
        Mdust = Mdust + adddust

        if (Mdust.lt.0.0d0) then
          Mdust = 0.0d0 !prevents unphisical Mdust < 0
        endif

        ! calculate average sigma within the Hill sphere
        sigmadust_temp = Mdust/(pi*(rm+af)**2.0d0-pi*(rm-af)**2.0d0)
        if (sigmadust_temp.lt.0.0d0) then
          sigmadust_temp = 0.0d0 ! stops unphysical sigma<0
        endif

        ! calculate Safronov (1969) core mass accretion
        Mdotcore = pi*Rc**2.0D0*sigmadust_temp*omega*Fg
        dMcore = Mdotcore*dt

        ! stop core accretion if Hill radius is too large
        if (RH.gt.0.84d0*H) then
          Mdotcore = 0.0d0
        endif

        ! reduce dt if dMcore is greater than the amount of dust available
        do while (dMcore.gt.Mdust)
          dt = dt/2.0d0
          dMcore = 0.0d0
          dMcore = Mdotcore*dt
        enddo
        ! also reduce dt if dMcore is greater than the user-defined tolerance
        do while (dMcore.gt.temptolD*Me)
          dt = dt/2.0d0
          dMcore = 0.0d0
          dMcore = Mdotcore*dt
        enddo

        TauKH = 1.0d9*((Mpl/Me)**(-3.0d0))*kappa ! Kelvin-Helmholtz timescale for gas envelope accretion
        Mcrit = 10.0d0*(Mdotcore*yr/Me/1.0d-6)**(1.0d0/4.0d0)*kappa**(1.0d0/4.0d0) ! critical mass for gas accreation
        ! gas is accreated if the mass exceeds the critical threshold
        if (Mcore/Me.gt.Mcrit) then
          Mdotgas = Mpl/TauKH/yr
          if (type2.ne.0) then
            if(rH.ge.0.84d0*H) then
              mu = 3.0d-7*rm**1.5d0
              Mdotgas = 3.0d0*3.14159*mu*sigmagas(index)*0.9
            endif
          endif
          else
            Mdotgas = 0.0d0 ! otherwise gas can escape
          endif
          ! calculates the gas isolation mass at the current timestep
          Mgiso = 50.0d0*(sigmagas(index)/2.4d3)**1.5d0*(rm/AU)**3.0*(Mstar/Msun)**(-0.5d0) ! includes decaying sigmagas
          ! if isolation mass is exceeded, no gas is accreated
          if (Mpl/Me.gt.Mgiso) then
            Mdotgas = 0.0d0
          endif

          ! reduce timestep if accreted gas exceeds user-defined tolerance
          do while (Mdotgas*dt.gt.temptolG*Mj)
            dt = dt/2.0d0
          enddo

!*************************************************************************
!  dt is now fixed for this time step, so calculations with dt can be made
!*************************************************************************

          ! recalculate temporary sigma for remaining dust, after accretion
          sigmadust_temp = (Mdust-dMcore)/(pi*(rm+af)**2.0D0-pi*(rm-af)**2.0D0)
          if(sigmadust_temp.le.0.0)then
            sigmadust_temp = 0.0d0
          endif
          ! recalculates the inner edge sigmadust
          sigmadust(index1)=(sigmadust(index1)*                             &
                             (pi*(rm-af)**2.0D0-pi*rad(index1)**2.0D0)+     &
                             sigmadust_temp*                                &
                             (pi*rad(index1+1)**2.0D0-pi*(rm-af)**2.0D0))/  &
                             (pi*rad(index1+1)**2.0D0-pi*rad(index1)**2.0D0)
         ! all the sigmadust inside hill sphere of planet are set equal
          do j = index1+1, index2-2
             sigmadust(j) = sigmadust_temp
          enddo
          ! recalculates the outside edge sigmadust
          sigmadust(index2-1)=(sigmadust(index2-1)*                     &
                 (pi*rad(index2)**2.0D0-                               &
                 pi*(rm+af)**2.0D0)+sigmadust_temp*                    &
                 (pi*(rm+af)**2.0D0-pi*rad(index2-1)**2.0D0))/         &
                 (pi*rad(index2)**2.0D0-pi*rad(index2-1)**2.0D0)

          do j=1,numpts-1 ! prevent unphysical sigmadust < 0
            if(sigmadust(j).lt.0.0) then
              sigmadust(j)=0.0d0
            endif
          enddo

          ! increment core mass
          Mcore = Mcore + dt*Mdotcore
          dMcore = dt*Mdotcore
          ! increment the total planet mass
          Mpl = Mpl + dMcore + dt*Mdotgas
          ! add up gas and dust which has been accreted
          AddUpGas(m) = AddUpGas(m) + dt*Mdotgas
          AddUpEatenDust(m) = AddUpEatenDust(m) + dMcore

        endif ! ends planet accretion loop

        finalMD(m) = 0.0d0
        do i=1,numpts-1
          finalMD(m)=finalMD(m)+(sigmadust(i)*((pi*rad(i+1))**2-(pi*rad(i))**2))
        enddo


!*************************************************************************
!       now mass accretion is complete we can calculate migration
!*************************************************************************
       ! type II migration
       if(type2.ne.0) then
         if(rH.ge.0.84d0*H)then
           mu = 3.0d-7*rm**1.5d0
           taudisc = rm**2.0d0/mu
           alpha = mu/0.05d0**2.0d0/rm**2.0d0/omega
           tauMig = (8.0d5*(Mpl/Mj)*(Msun/Mstar)*               &
                    (1.0d-4/alpha)*(rm/au)**0.5d0*yr)/          &
                    (sigmagas(index)*(rm/au)**1.5d0/2.4d3)

           if (tauMig.lt.taudisc) then
             tauMig = taudisc
           endif
           type2time = 0.0d0
           type2dt = dt
           drbdt = rm/tauMig
           do while((drbdt*type2dt)/Au.gt.0.001)
             type2dt = type2dt/2.0d0
           enddo
           do while(type2time.lt.dt)
             rm = rm - drbdt*type2dt
             marker(m,2) = marker(m,2) + drbdt*type2dt
             mu = 3.0d-7*rm**1.5d0
             taudisc = rm**2.0d0/mu
             tauMig = (8.0d5*(Mpl/Mj)*(Msun/Mstar)*              &
                      (1.0d-4/alpha)*(rm/au)**0.5d0*yr)/         &
                      (sigmagas(index)*(rm/au)**1.5d0/2.4d3)
             if (tauMig.lt.taudisc) then
               tauMig = taudisc
             endif
             drbdt = rm/tauMig
             type2time = type2time + type2dt
           enddo

           ! recalculate conditions for new radius
           RH = rm*(Mpl/3.0d0/Mstar)**(1.0d0/3.0d0)
           Tm = temper(index)
           csm = DSqrt(7.0d0/5.0d0*k*Tm/muH/mH)
           Omega = DSqrt(G*Mstar/rm**3.0d0)
           H = csm/omega
         endif
       endif

       ! type I migration
       if (type1.ne.0) then
         if (rH.lt.0.84d0*H) then
           if (SIGarray(Tindex,sigindex).gt.0.0d0) then
             dlogSdlogr = Dlog(SIGarray(Tindex,sigindex+1)) -            &
                          Dlog(SIGarray(Tindex,sigindex-1))
             dlogSdlogr = dlogSdlogr/(DLog(RADarray(sigindex+            &
                          1)) - Dlog(RADarray(sigindex-1)))
           else
             dlogSdlogr = 0.0d0
           endif
           betaMig = -1.0d0*dlogSdlogr
           type1drbdt = (2.7d0+1.1d0*betaMig)*(Mpl/Mstar)                &
                        *(sigmagas(index)*(rm**2.0d0)/Mstar)*(rm*Omega   &
                        /csm)**2.0d0*Omega*rm
           type1drbdt = type1drbdt/type1factor
           type1max30 = (2.7d0+1.1d0*betaMig)*50.0d0*Me/Mstar*           &
                        (sigmagas(index)*(rm**2.0d0)/Mstar)*             &
                        (rm*Omega/csm)**2.0d0*Omega*rm
           type1dir = 1.0d0 !the default direction is decreasing radius
           if (randomMig.ne.0) then !code for turbulent type 1
             type1time = 0.0d0
             do while (type1time.lt.dt)
               ranpm = rand(0)*(sigorb*2.0d0)
               norbit = (10.0d0-sigorb) + ranpm !~10 orbits before random dir chang
               vel = sqrt(Mstar*G/rm)
               dist = 2*pi*rm*norbit
               type1dt = dist/vel
               norb(m) = norb(m) + norbit
               ranpm = rannos(int(rand(0)*100000))
               if (ranpm.gt.0.5) then
                 type1dir = -1.0d0 !50/50 chance for each direction
               else
                 type1dir = 1.0d0
               endif
               if (type1drbdt.lt.type1max30) then
                 rm = rm - type1dir*type1drbdt*type1dt
               else
                 rm = rm - type1dir*type1max30*type1dt -(type1drbdt-type1max30)*type1dt
               endif
               marker(m,1) = marker(m,1) + sqrt((type1drbdt*type1dt)**2.0d0)
               type1time = type1time + type1dt
               if (type1time.gt.dt) then !corrects too much orbiting
                 rm = rm + type1dir*type1drbdt*(type1time-dt)
                 type1time = type1time - (type1time-dt)
                 exit
               endif
             enddo
             !type 1 with no turbulance
             else
               rm = rm - type1dir*type1drbdt*dt
               marker(m,1) = marker(m,1) + type1drbdt*dt
             endif

             ! recalculate conditions for new radius
             RH = rm*(Mpl/3.0D0/Mstar)**(1.0D0/3.0D0) !recalculates the hill radius
             Tm = temper(index)
             csm = DSqrt(7.0d0/5.0d0*k*Tm/muH/mH) !
             Omega = DSqrt(G*Mstar/rm**3.0d0) !a Kepler frequency
             H = csm/omega !scale height

             if (rH.ge.0.84d0*H) goto 55

          endif
       endif

       ! remove planet out of bounds
       if (rm.gt.rOut) then
         check(m) = 2
         exit
       elseif(rm.lt.rIn) then
         check(m) = 2
         exit
       endif
55     continue

       ! recalculate rad array index where the planet is located
       do i = 1, numpts-1
         if ((rad(i).le.rm).and.(rad(i+1).gt.rm)) then
           index = i
         endif
       enddo
       do i = 1, 499
         if((RADarray(i).le.rm).and.(RADarray(i+1).gt.rm)) then
            sigindex = i
         endif
       enddo

       if(numpl.eq.1.and.numdisc.eq.1)then
          open(unit=26,file='migration.dat',status='unknown')
          write(26,*) time,' ',rm/Au
       endif

       if(numpl.eq.1.and.numdisc.eq.1)then
          open(unit=28,file='mvt.dat',status='unknown')
          write(28,*) time,' ',mcore/me,' ',mpl/me
       endif

       if(numpl.eq.1.and.numdisc.eq.1)then
          open(unit=29,file='mvr.dat',status='unknown')
          write(29,*) rm/Au,' ',mpl/me
       endif

       if(numpl.eq.1.and.numdisc.eq.1)then
          open(unit=30,file='hvr.dat',status='unknown')
          write(30,*) rm/Au,' ',H/Au
       endif

       if(numpl.eq.1.and.numdisc.eq.1)then
          open(unit=32,file='hvrh.dat',status='unknown')
          write(32,*) rm/Au,' ',rH/Au
       endif

       if(rH.ge.0.84d0*H.and.done.eq.1)then
          open(unit=31,file='t2start.dat',status='unknown')
          write(31,*) rm/Au,' ',mpl/me
          done = 0
       endif

       time = time + dt/yr !counts the time

!     exit if acreation virtually stopped
       if(dMcore/me.lt.1d-30.and.dt*Mdotgas/me.lt.1d-6)then
          print*,'exit, growth',m,' negligible'
          nex = nex + 1
          goto 17
       endif
       ! resets origioal timestep
       dt = origdt
    end do !ends time loop

!*************************************************************************
!             disc evolution complete, now write to files
!*************************************************************************

! this continue is for a simulation exit due to too big RH
17  continue

    ! store final values
    finalcoresR(m)=rm/au
    finalcoreM(m)=Mcore
    finalpM(m)=Mpl

    open(unit=20,file='final.dat')
    write(19,*) time, fdust, fgas, rm/au, Mgiso, Mpl/Me, Mcore/Me, ocean, &
                check(m), rH/Au, af/Au, marker(m,1)/Au, marker(m,2)/Au
    close(20)

    ! sum final dust mass
    finalMD(m) = 0.0
      do i=1,numpts-1
        finalMD(m)=finalMD(m)+(sigmadust(i)*(pi*rad(i+1)**2-pi*rad(i)**2))
      enddo

      if ((Mpl/Mj).gt.masslimit) then
        check(m)= 3
      endif

      if(check(m).eq.1) then
        print*,'type1 for',marker(m,1)/Au,' Au'
        print*,'type2 for',marker(m,2)/Au,' Au'
        print*,'initial core ',m,'''s orbit rad =',initialcoresR(m)
        print*,'final core ',m,'''s orbital radius =',finalcoresR(m)
        print*,'Mdust added to core = ',AddUpEatenDust(m)/Me,' Me'
        print*,'Mgas added to planet =',AddUpGas(m)/Mj,' Mj'
        print*,'fianl mass of core',m,' =',finalcoreM(m)/Me,' Me'
        print*,'final planet mass =',finalpM(m)/Me,' Me'
        print*,' '
      else
        print*,'planet failed'
        print*,'type1 for',marker(m,1)/Au,' Au'
        print*,'type2 for',marker(m,2)/Au,' Au'
        print*,'initial core ',m,'''s orbit rad =',initialcoresR(m)
        print*,'final core ',m,'''s orbital radius =',finalcoresR(m)
        print*,'final planet mass =',finalpM(m)/Me,' Me'
        print*,' '
      endif

      ! total  type 1 and type 2 migration
      open(unit=33,file='totalt1.dat',status='unknown') !total amaount of t1 migration
      open(unit=34,file='totalt2.dat',status='unknown')
      write(33,*) marker(m,1)/Au,' ',Mpl/Me
      write(34,*) marker(m,2)/Au,' ',Mpl/Me

    ! end num planets loop
    enddo

    ! count the number of planets exceeding the mass limit for each disc
    numstar = 0
    numlhr = 0
    numok = 0
    do i=1,numpl
      if (check(i).eq.3) then
        numstar = numstar + 1
      elseif(check(i).eq.2) then
        numlhr = numlhr + 1
      elseif(check(i).eq.1) then
        numok = numok + 1
      endif
   enddo
   print*,' '
   print*,'planets exceeding mass limit for this fdust=',numstar
   print*,'number stopped due to negligible mass gain',nex

   open(unit=12,file='coreoutput.dat',status='unknown') !coremass vs radius
   open(unit=13,file='ploutput.dat',status='unknown') !planetmass vs radius
   open(unit=14,file='allpl.dat',status='unknown') !all planets (over mass limit + out of bounds)
   open(unit=15,file='allcore.dat',status='unknown') !ditto for coremass vs radius
   open(unit=16,file='simrep.dat',status='unknown') !summing up of each disc type
   open(unit=17,file='discmass.dat',status='unknown') !total mass of dust at start
   open(unit=18,file='gasmass.dat',status='unknown') !total mass of gas at start
   open(unit=27,file='init.dat',status='unknown') !initial positions of planets vs end mass

   write(16,*) 'disc',l
   write(16,*) 'num exceeding mass limit',numstar
   write(16,*) 'num planets out of bounds',numlhr
   write(16,*) 'num of planets produced',numok
   write(16,*) ' '
   nok = nok + numok
   nob = nob + numlhr
   nml = nml + numstar

   write(17,*) fdust,' ',initialMD(1)/Mstar
   write(18,*) fdust,' ',Mgas(1)/Mstar

   open(unit=20,file='fd1.dat',status='unknown')
   open(unit=21,file='fd2.dat',status='unknown')
   open(unit=22,file='fd3.dat',status='unknown')
   open(unit=23,file='fd4.dat',status='unknown')
   open(unit=24,file='fd5.dat',status='unknown')
   open(unit=25,file='fd6.dat',status='unknown')
   range = maxfd/6.0d0
   do i=1,numpl
     if (check(i).ne.2.and.check(i).ne.3) then
       if(fdust.le.range)then
         write(20,*) finalcoresR(i),(finalpM(i)/Me)
       elseif(fdust.le.2*range)then
         write(21,*) finalcoresR(i),(finalpM(i)/Me)
       elseif(fdust.le.3*range)then
         write(22,*) finalcoresR(i),(finalpM(i)/Me)
       elseif(fdust.le.4*range)then
         write(23,*) finalcoresR(i),(finalpM(i)/Me)
       elseif(fdust.le.5*range)then
         write(24,*) finalcoresR(i),(finalpM(i)/Me)
       elseif(fdust.gt.5*range)then
         write(25,*) finalcoresR(i),(finalpM(i)/Me)
       endif
     endif
   enddo

   do i=1,numpl
     if(check(i).eq.1)then
       write(12,*) finalcoresR(i), (finalcoreM(i)/Me)
       write(13,*) finalcoresR(i), (finalpM(i)/Me)
       write(27,*) initialcoresR(i), finalpM(i)/Me
     elseif(check(i).ne.1)then
       write(14,*) finalcoresR(i), (finalpM(i)/Me)
       write(15,*) finalcoresR(i), (finalcoreM(i)/Me)
     endif
   enddo

enddo ! ends loop over discs, line 128

print*,' '
print*,'total number of planets',nok
print*,'total number out of bounds',nob
print*,'total munber over the mass limit',nml

close(12)
close(13)
close(14)
close(15)
close(16)
close(17)
close(18)
close(20)
close(21)
close(22)
close(23)
close(24)
close(25)
close(26)
close(27)
close(28)
close(29)
close(30)
close(31)
close(32)
close(33)
close(34)

! end the main program
end

subroutine calcFg(Fg,Rc,RH,af,omega)

  real*8 Fg, Rc, RH, af, omega
  real*8 erf
  real*8 rhopl, Mpl, rpl, vescpl
  real*8 i, iH, e, eH, dc
  real*8 S1, S2, S3, beta, B, eps1, eps2
  real*8 AU, G, pi, k, Msun, stefan
  external erf
  common /const/ AU, G, pi, k, Msun, stefan
  rpl = 1.0d7
  rhopl = 3.2d0
  Mpl = rhopl*4.0d0/3.0d0*pi*rpl**3.0d0
  vescpl = DSqrt(2.0d0*G*Mpl/rpl)
  iH = vescpl/Sqrt(3.0d0)/omega/RH
  eH = max(2.0d0*iH,2.0d0)
  af = DSqrt(12.0d0 + eH**2.0d0)*RH
  dc = Rc/RH
  eps1 = 0.2d0*iH

  if (dc.lt.0.0284d0) then
    eps2 = Sqrt(2.2d0*dc)/iH + eps1
  else
    eps2 = 0.25d0/iH + eps1
  endif

  S1 = eps1*Sqrt(pi)*(erf(eps1)-erf(eps2))
  S1 = Exp(-1.0*eps1**2.0)-Exp(-1.0*eps2**2.0)+S1
  S1 = 4.0*Exp(eps1**2.0)/dc**1.5*S1

  if (dc.lt.0.0284) then
    S2 = erf(0.25/iH)-erf(Sqrt(2.2*dc)/iH)
    S2 = 6.69/dc/iH*S2
  else
    S2 = 0.0
  endif
  S3 = 1.0 + 3.31/dc/iH**2.0

  beta = 1.0/6.0*((log10(iH)-0.131)/0.343)**3.0
  beta = (Dlog10(iH)-0.131d0)/0.343d0 + beta
  beta = 1.0d0 + Exp(beta)
  beta = pi/2.0d0/beta

  B = 1.0D0/24.0D0*((log10(iH)+0.109D0)/0.231D0)**4.0D0
  B = 1.0D0 + 0.422D0*Exp(-1.0D0*B)

  Fg = B*((S1+S2)*Sin(beta)+S3*cos(beta))

  return
end

function erf(x)

  real*8 erf, x
  real*8 gammp

  if (x.lt.0.0d0) then
    erf = -gammp(0.5d0,x**2.0d0)
  else
    erf = gammp(0.5d0,x**2.0d0)
  endif

  return
end

function gammp(a,x)
  real*8 a,gammp,x
  real*8 gammcf, gamser, gln

  if(x.lt.0.0d0.or.a.le.0.0d0)stop
  if(x.lt.a+1.0d0)then
    call gser(gamser,a,x,gln)
    gammp=gamser
  else
    call gcf(gammcf,a,x,gln)
    gammp=1.-gammcf
  endif
  return
end

subroutine gser(gamser,a,x,gln)

  integer itmax, n
  real*8 a,gamser,gln,x,eps
  real*8 ap,del,sun,gammln
  parameter (itmax=100,eps=3.d-7)
  external gammln
  gln=gammln(a)

  if(x.le.0.0d0)then
    if(x.lt.0.0d0)stop
    gamser=0.0d0
    return
  endif
  ap=a
  sum=1.0d0/a
  del=sum
  do n=1,itmax
    ap=ap+1.0d0
    del=del*x/ap
    sum=sum+del
    if(abs(del).lt.abs(sum)*eps)go to 1
  enddo
  stop 'a too large, itmax too small'
  1    gamser=sum*exp(-x+a*dlog(x)-gln)

  return
end

subroutine gcf(gammcf,a,x,gln)

  integer itmax
  real*8 a,gammcf,gln,x,eps,fpmin
  real*8 an,b,c,d,del,h,gammln
  parameter (itmax=100,eps=3.e-7)
  external gammln
  gln=gammln(a)
  gold=0.0d0
  a0=1.0d0
  a1=x
  b0=0.0d0
  b1=1.0d0
  fac=1.0d0
  do n=1,itmax
    an=float(n)
    ana=an-a
    a0=(a1+a0*ana)*fac
    b0=(b1+b0*ana)*fac
    anf=an*fac
    a1=x*a0+anf*a1
    b1=x*b0+anf*b1
    if(a1.ne.0.0d0)then
      fac=1./a1
      g=b1*fac
      if(abs((g-gold)/g).lt.eps)go to 1
      gold=g
    endif
  enddo
  stop 'a too large, itmax too small'
  1    gammcf=dexp(-x+a*dlog(x)-gln)*g

  return
end

function gammln(xx)
  real*8 gammln,xx
  real*8 cof(6),stp,half,one,fpf,x,tmp,ser
  data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,-1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
  data half,one,fpf/0.5d0,1.0d0,5.5d0/
  x=xx-one
  tmp=x+fpf
  tmp=(x+half)*dlog(tmp)-tmp
  ser=one

  do j=1,6
    x=x+one
    ser=ser+cof(j)/x
  enddo
  gammln=tmp+dlog(stp*ser)

  return
end
