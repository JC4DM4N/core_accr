program p1

implicit none
integer numpts, l, m, i, j, n, Tary, Rary, done
integer Tindex, sigindex, sigIn, sigOut
integer numstar,numlhr,numok
integer index, index1, index2, numpl, numdisc
integer ix, iseed, iz, ocean, nok,nob,nml,nex
integer type1, type2, randomMig, numorb, sigorb, type1dir, obt
integer*4 timeArray(3)    ! Holds the hour, minute, and second
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
real*8 temper(numpts), a_ice, eta_ice
real*8 initialcoresR(100), finalcoresR(100)
real*8 initialMD(100), finalMD(100), finalcoreM(100), finalpM(100)
real*8 AddUpGas(100), AddUpEatenDust(100), check(100), Mgas(100)
real*8 r, Rcore, r_in, r_out, del_r, ro, rm, rm1, rm2, origrm
real*8 time, dt, type1dt, type1time, type2dt, type2time, origdt,yr
real*8 type1factor
real*8 disclife, decayTime
real*8 Me, Msun, Mstar, Mcore, Mpl, massLimit, Mj, Mdisc, Mdust
real*8 Mdotcore, dMcore, Mdotgas
real*8 fdust, fgas, fd, maxfd, minfd, range, adddust
real*8 totmass, totmass_old
real*8 TauKH, Mcrit, Mgiso
real*8 sigmadust_temp
real*8 muH, muHe, mH, XX, YY
real*8 rand_no
real*8 vel,dist,orb,norb(100),norbit,ranpm,pm,rannos(100000),tdt2
real*8 Tarray(1000),RADarray(500),SIGarray(1000,500),OSA(1000,500)
real*8 Tfrac, SIGfrac, gasSig, scaleM, scaleF
real*8 sigma5au, signorm
real rand, x, y
COMMON /mol/ muH, muHe, mH, XX, YY
COMMON /const/ AU, G, pi, k, Msun, stefan

open(unit=25,file='randnums.dat',status='old')
read(25,*) rannos
close(25)

! Define constants
k = 1.3807D-16
pi = DACos(-1.0D0)

AU = 1.496D13 !cm
G = 6.672D-8 !cm^3 g^-1 s^-2
stefan = 5.670D-5 !boltzmann constant

Msun = 1.989D33
Me = 5.979D27
Mj = 317.8*Me

MuH = 2.0D0
muHe = 3.971D0
muHe = 4.0D0
mH = 1.673D-24

Mstar = 2.4D0*Msun

beta = 1.0d0

!read params
open(unit=19,file='core_accr.params',status='old')
read(19,*) rm
read(19,*) Mcore
read(19,*) sigma5au
close(19)

rm = rm*AU
Mcore = Mcore*Me

!    read in gas disc values for 500,000 years for disc from 0.0506Au to 19.854Au
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

yr = 365.25D0*24.0D0*60.0D0*60.0D0
origdt = 10000.0D0*yr
disclife = Tarray(1000)/yr
massLimit = 15.0 !max mass is 15 jupiter masses
scaleM = 0.01d0 !Mgas * this = Mdust
scaleF = 0.41d0

!     In order to get a different sequence each time, we initialize the
!     seed of the random number function with the sum of the current
!     hour, minute, and second.

call itime(timeArray)     ! Get the current time
iz = rand (timeArray(1)+timeArray(2)+timeArray(3))

!     sets random number interval [0.25, 12] for selecting rm values
rm1 = 0.25d0
rm2 = 15.0d0
nok = 0
nob = 0
nml = 0

type1 = 0
type1factor = 10.0d0
type2 = 0
randomMig = 0
numorb = 10
sigorb = 2

tolD = 0.02 !tolerance for mdotcore (*Me)
tolG = 0.02 !tolerance for mdotgas (*Mj)
numpl = 1 ! number of different radii to try for each disc
numdisc = 1 !also used for calculating fdust, line 105
!     fdust scales with star mass
!     maximum mass varies with solar mass
      maxfd = 2.492533*scaleM*scaleF*(Mstar/Mj)
!     minimum fdust is 0.1(Mstar/Msun)
fd = maxfd - 0.1*(Mstar/Msun)
!fd = fd/(numdisc-1)
minfd = 0.1*(Mstar/Msun)

do l = 1, numdisc !loop over different discs
nex = 0!counts number of exits(by small Mdot) in each disc
!     a loop to set all planets as being initially in the simulation
 do 94 i=1,numpl
    check(i) = 1 !1=ok, 2=rh>disc, 3=mass limit exceeded
94      continue

!     calculate fdust for the current disc
fdust = minfd + fd*(l-1)
if(numdisc.eq.1) fdust = maxfd
print*,' '
print*,'l =',l,'  fdust =',fdust, numdisc

kappa = 1.0D0

ro = 1.0D0*AU          ! Normalisation radius for densities
origrm = 0.25D0*AU          ! Initial planet radius - usually 0.25

!     a loop to place the core at different radii in the disc from inner to outer disc rad
do m = 1,numpl
    done = 1 !a parameter to ensure only 1 print per planet(line 640)
    AddUpGas(m) = 0.0
    AddUpEatenDust(m) = 0.0

    !rm = ((rand(0)*(rm2-rm1))+rm1)*Au  !randomly place the core
    !rm = 10*AU

    print*, 'Initial planet radius, ', rm/AU, ' AU'

    initialcoresR(m)=rm/Au !notes initial radius
    marker(m,1) = 0.0d0
    marker(m,2) = 0.0d0
!     Mass acreation tolerance adjusted according to radius
    temptolD = tolD*(alog(real(1+rm/Au)))**3.0d0
    temptolG = tolG*(alog(real(1+rm/Au)))**3.0d0

    r_in =  RADarray(1)   ! Inner edge of grid
    r_out = RADarray(500)   ! Outer edge of grid
    del_r = (r_out-r_in)/numpts  !divides up disc rad for even spacing

    rad(1) = r_in

!  spaces out rad evenly from the inner radius
    do i = 2, numpts
       rad(i) = rad(i-1) + del_r
    enddo

!   sets the temperature at depending on radius
    do i = 1, numpts-1
       temper(i) = 280.0d0*(ro/rad(i))**0.5d0*(Mstar/Msun)
    enddo

    a_ice = 2.7d0*(Mstar/Msun)**2.0d0*AU  !sets the snow line
    ocean = 0
    If (rm.gt.a_ice) ocean=1
    totmass = 0.0D0  !starts at 0 and then added up each again in each loop

!     this loop scales the gas surface density
    do i=1,1000
       do j=1,500
!                  SIGarray(i,j) = OSA(i,j)*fdust/maxfd*scaleF
!                  SIGarray(i,j) = OSA(i,j)*fdust*0.09242d0*(Msun/Mstar)
           SIGarray(i,j) = OSA(i,j)
       enddo
    enddo

!     this loop calculates sigma gas based on a realistic model
    i = 1               !start rad() off at 1
    do j=1,500
       do while(RADarray(j).ge.rad(i))
          if(i.eq.1)then !do the 1st point seperately (no 0 element)
             sigmagas(i) = SIGarray(1,j)
          else
             SIGfrac = SIGarray(1,j-1) - SIGarray(1,j)
             SIGfrac = SIGfrac/(RADarray(j) - RADarray(j-1))
             sigmagas(i) = SIGarray(1,j-1) - SIGfrac*(rad(i) - RADarray(j-1))
          endif
          if(i.eq.numpts)goto 42
          i = i+1
       enddo
    enddo
42         continue

!     this loop calculates sigmadust based on a simple model, then adds up the total mass
    iprint = 0
    do i = 1, numpts-1
       If (rad(i).lt.a_ice) Then
          eta_ice = 1.0d0
       Else
          eta_ice = 4.2d0
       EndIf

       sigmadust(i) = fdust*eta_ice*10.0d0*(ro/rad(i))**(beta)

       !sigmadust(i) = sigmadust(i)/1.53d0    ! 15gcm^-2
       !sigmadust(i) = sigmadust(i)/2.295d0   ! 10gcm^-2
       !sigmadust(i) = sigmadust(i)/3.06d0    ! 7.5gcm^-2

       totmass = totmass + &
                 (pi*rad(i+1)**2.0D0-pi*rad(i)**2.0D0)*sigmagas(i) + &
                 (pi*rad(i+1)**2.0D0-pi*rad(i)**2.0D0)*sigmadust(i)

       if (rad(i).lt.5.001*au .and. rad(i).gt.4.999*au .and. iprint.eq.0) then
			  signorm = sigmadust(i)/sigma5au
          print*, ''
          print*, 'SIGMA AT 5AU, GAS:', sigmagas(i), 'DUST:', sigmadust(i)/signorm, 'gcm-2'
          print*, ''
          iprint = 1
       endif

    enddo !end ln 124
    ! normalise sigmadust
    sigmadust = sigmadust/signorm

    open(12345, file='densities.dat', status='unknown')
    do i=1,numpts-1
      write(12345, *) rad(i), sigmagas(i), sigmadust(i)
    enddo
    close(12345)

    initialMD(m) = 0.0
    do 2 i=1,numpts-1
       initialMD(m)=initialMD(m)+(sigmadust(i)*(pi*rad(i+1)**2-pi*rad(i)**2))
2          continue

    totmass_old = totmass !keeps track of the origional mass

    !     this loop sets "index" to be the number in the radius array where the planet is
    do i = 1, numpts-1
       If ((rad(i).le.rm).and.(rad(i+1).gt.rm)) Then
          index = i
       EndIf
    enddo
    do i = 1, 499
       if((RADarray(i).le.rm).and.(RADarray(i+1).gt.rm)) then
          sigindex = i
       endif
    enddo

    !Mcore = 0.01d0*Me  !sets the core mass to 0.1 earth mass
    Mpl = Mcore
    Mgiso = 50.0d0*(sigmagas(index)/2.4d3)**1.5d0*(rm/AU)**3.0*(Mstar/Msun)**(-0.5d0)

    time = 0.0D0
    dt = origdt  !keeps track of what dt is origionally set as

    ix = -1  !a seed for random number generators
    iseed = -1 !a seed for random number generators

    call ranset(ix,iseed) !generate 1st random number

    drbdt = 0.0d0  !initial core migration is zero

    Mgas(m) = 0.0d0
    do 65 i=1,numpts-1
       Mgas(m)=Mgas(m)+(sigmagas(i)*(pi*rad(i+1)**2-pi*rad(i)**2))
65         continue

    print*,'Mgas=',Mgas(m)/Msun,' Msol, Mdust=',initialMD(m)/Msun,' Msol'
    print*, 'Current dust-to-gas ratio=', initialMD(m)/Mgas(m)
    print*, 'Initial dust-to-gas ratio=', 2.0d0*initialMD(m)/Mstar
    print*,' '


!     now that all the initial conditions have been set a time loop is used to itterate
!     over the lifetime of the disc

    Do while (time .lt. disclife)

       Tm = temper(index) !sets the planets T to that of the surounding dust
       csm = DSqrt(7.0d0/5.0d0*k*Tm/muH/mH) !
       Omega = DSqrt(G*Mstar/rm**3.0d0) !a Kepler frequency
       H = csm/omega    !scale height
       decayTime = DExp(-time/disclife)
!     finds out the current array element for time
       do i=1,1000
          Tindex = i
          if(Tarray(i+1).gt.time*yr) exit
       enddo
!     calculates how far along to the next element the current time is
       Tfrac = time*yr - Tarray(Tindex)
       Tfrac = Tfrac/(Tarray(Tindex+1) - Tarray(Tindex))
!     now the sigmagas array is updated for this timestep
        i = 1               !start rad() off at 1
        do j=1,500
           do while(RADarray(j).ge.rad(i))
              if(i.eq.1)then !do the 1st point seperately (no 0 element)
                 sigmagas(i) = (SIGarray(Tindex,j) - SIGarray(Tindex+1,j))*Tfrac
                 sigmagas(i) = SIGarray(Tindex,j) - sigmagas(i)
              else
                 SIGfrac = SIGarray(Tindex,j-1) - SIGarray(Tindex,j)
                 SIGfrac = SIGfrac/(RADarray(j) - RADarray(j-1)) !last timestep
                 gasSig = SIGarray(Tindex+1,j-1) - SIGarray(Tindex+1,j)
                 gasSig = SIGfrac/(RADarray(j) - RADarray(j-1))  !next timestep
                 sigmagas(i) = (SIGfrac - gasSig)*Tfrac
                 sigmagas(i) = SIGarray(Tindex,j-1) - sigmagas(i)
              endif
              if(sigmagas(i).lt.0.0d0)then
                 sigmagas(i) = 0.0d0 !prevents unphysical sigmagas
              endif
              if(i.eq.numpts)goto 43
              i = i+1
           enddo
        enddo
    43         continue

!  checks the planet hasnt exceeded its isolation mass or mass limit
       If ((Mpl/Me.lt.Mgiso).and.(Mpl/Mj.lt.massLimit)) Then
          rhom = gasSig/H
          rhocore = 3.2D0 ! Density of core
          Rcore = (3.0D0*Mcore/(4.0D0*pi*rhocore))**(1.0D0/3.0D0)
          Rc = Rcore

          call Dran(ix,rand_no)  !create random number

          RH = rm*(Mpl/3.0D0/Mstar)**(1.0D0/3.0D0) !calculates the hill radius

          call calcFg(Fg,Rc,RH,af,omega) !calculates the G enhancement

!     this sets the index1(inside) & index2(outside) values to the areas around the planets
!     area by calculating weather they are within its hill sphere
          do j = 1, numpts-1
             If ((rad(j).le.(rm-af)).and.(rad(j+1).gt.(rm-af))) Then
                index1 = j
             EndIf
             If ((rad(j).lt.(rm+af)).and.(rad(j+1).ge.(rm+af))) Then
                index2 = j+1
             EndIf
          enddo

!     condition to end simulation and discard planet for hill radius outwith disc
          if((rm-af).lt.r_in.or.(rm+af).gt.r_out) then
             check(m) = 2
             goto 17
          endif

          Mdust = 0.0D0  !firstly dust mass reset
!     1st addition is from the inner edge of the hill sphere
          adddust = 0.0d0
          adddust =(pi*rad(index1+1)**2.0D0-pi*(rm-af)**2.0D0)*sigmadust(index1)
          Mdust = Mdust + adddust

!     this loop then covers all the areas within the planets hill sphere and adds
!     the dust mass in them to the value Mdust, the planets potential feeding mass
          do j = index1+1, index2-2
             adddust = 0.0d0
             adddust =(pi*rad(j+1)**2.0D0-pi*rad(j)**2.0D0)*sigmadust(j)
             Mdust = Mdust + adddust
          enddo

!     final contribution from outside edge of hill sphere
          adddust = 0.0d0
          adddust =(pi*(rm+af)**2.0D0-pi*rad(index2-1)**2.0D0)*sigmadust(index2-1)
          Mdust = Mdust + adddust

          If (Mdust.lt.0.0D0) Then
             Mdust = 0.0D0 !prevents unphisical Mdust < 0
          EndIf

!     sets a temporart sigma over feeding zone to the average of outer & inner edges
          sigmadust_temp = Mdust/(pi*(rm+af)**2.0D0-pi*(rm-af)**2.0D0)

          If (sigmadust_temp.lt.0.0D0) Then
             sigmadust_temp = 0.0D0 !stops unphysical sigma<0
          EndIf

!     Equations to increase the core mass by the appropriate ammount
          Mdotcore = pi*Rc**2.0D0*sigmadust_temp*omega*Fg
          dMcore = Mdotcore*dt

          If (RH .gt. 0.84d0*H) Then  ! Stop core accretion..
            Mdotcore = 0.0d0
          EndIf

          do while(dMcore.gt.Mdust)
             dt = dt/2.0d0
             dMcore = 0.0d0
             dMcore = Mdotcore*dt
          enddo

          do while(dMcore.gt.temptolD*Me)
             dt = dt/2.0d0
             dMcore = 0.0d0
             dMcore = Mdotcore*dt !dont acreate more than 0.05earth mass at a time
          enddo

          TauKH = 1.0d9*((Mpl/Me)**(-3.0d0))*kappa !Kelvin Helmholtz timescale
          Mcrit = 10.0d0*(Mdotcore*yr/Me/1.0d-6)**(1.0d0/4.0d0)*kappa**(1.0d0/4.0d0) !critical mass for gas accreation


!     gas is accreated if the mass exceeds the critical threshold
          If (Mcore/Me.gt.Mcrit) Then
             Mdotgas = Mpl/TauKH/yr

             if(type2.ne.0) then
               if(rH.ge.0.84d0*H)then
                 mu = 3.0d-7*rm**1.5d0

                 Mdotgas = 3.0d0*3.14159*mu*sigmagas(index)*0.9
               EndIf
             EndIf
          Else
             Mdotgas = 0.0d0 !otherwise gas can escape
          EndIf

!     calculates the gas isolation mass corresponding to gap opening
          Mgiso = 50.0d0*(sigmagas(index)/2.4d3)**1.5d0*(rm/AU)**3.0*(Mstar/Msun)**(-0.5d0) !includes decaying sigmagas

          If (Mpl/Me.gt.Mgiso) then
             Mdotgas = 0.0d0 !if isolation is achieved no gas accreated
          EndIf

!     This bit reduces the time step so a smaller amount of gas is
!     added each run, so that Mgiso isn't greatly exceeded
          Do while (Mdotgas*dt.gt.temptolG*Mj)
             dt = dt/2.0d0
          EndDo

!*************************************************************************
!c dt is now fixed for this time step, so calculations with dt can be made
!*************************************************************************

!     recalculate temporary sigma for remaining dust
          sigmadust_temp = (Mdust-dMcore)/(pi*(rm+af)**2.0D0-pi*(rm-af)**2.0D0)
          if(sigmadust_temp.le.0.0)then
             sigmadust_temp = 0.0d0
          endif

!     recalculates the inside edge sigmadust
          sigmadust(index1)=(sigmadust(index1)*                             &
                             (pi*(rm-af)**2.0D0-pi*rad(index1)**2.0D0)+     &
                             sigmadust_temp*                                &
                             (pi*rad(index1+1)**2.0D0-pi*(rm-af)**2.0D0))/  &
                             (pi*rad(index1+1)**2.0D0-pi*rad(index1)**2.0D0)

!     all the sigmadust inside hill sphere of planet are set equal
          do j = index1+1, index2-2
             sigmadust(j) = sigmadust_temp
          enddo

!     recalculates the outside edge sigmadust
          sigmadust(index2-1)=(sigmadust(index2-1)*                     &
                 (pi*rad(index2)**2.0D0-                               &
                 pi*(rm+af)**2.0D0)+sigmadust_temp*                    &
                 (pi*(rm+af)**2.0D0-pi*rad(index2-1)**2.0D0))/         &
                 (pi*rad(index2)**2.0D0-pi*rad(index2-1)**2.0D0)

          do 15 j=1,numpts-1!prevent unphysical sigmadust < 0
             if(sigmadust(j).lt.0.0)then
                sigmadust(j)=0.0d0
             endif
15              continue

          Mcore = Mcore + dt*Mdotcore !changes the cores mass due to accretion
          dMcore = dt*Mdotcore !makes a note of the magnitude of the change

!     the cores mass is altered dependint on total accreation
          Mpl = Mpl + dMcore + dt*Mdotgas
          AddUpGas(m) = AddUpGas(m) + dt*Mdotgas
          AddUpEatenDust(m) = AddUpEatenDust(m) + dMcore

!                  print*, Mcore/Me, Mpl/Me, Mcrit, rm/AU

       end if !end mass limit $ Mgiso

       finalMD(m) = 0.0d0
       do 66 i=1,numpts-1
       finalMD(m)=finalMD(m)+(sigmadust(i)*((pi*rad(i+1))**2-(pi*rad(i))**2))
66            continue


!**************** Type II migration *****************

       if(type2.ne.0) then
          if(rH.ge.0.84d0*H)then
             mu = 3.0d-7*rm**1.5d0
             taudisc = rm**2.0d0/mu

             alpha = mu/0.05d0**2.0d0/rm**2.0d0/omega
!                     alpha = 1.0d-4

             tauMig = (8.0d5*(Mpl/Mj)*(Msun/Mstar)*                 &
                        (1.0d-4/alpha)*(rm/au)**0.5d0*yr)/          &
                        (sigmagas(index)*(rm/au)**1.5d0/2.4d3)

             If (tauMig.lt.taudisc) Then
                tauMig = taudisc
             EndIf

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
                         (1.0d-4/alpha)*(rm/au)**0.5d0*yr)/        &
                         (sigmagas(index)*(rm/au)**1.5d0/2.4d3)

                If (tauMig.lt.taudisc) Then
                   tauMig = taudisc
                EndIf
                drbdt = rm/tauMig
                type2time = type2time + type2dt
             enddo

!     recalculate conditions for new radius
             RH = rm*(Mpl/3.0D0/Mstar)**(1.0D0/3.0D0) !recalculates the hill radius
             Tm = temper(index)
             csm = DSqrt(7.0d0/5.0d0*k*Tm/muH/mH) !
             Omega = DSqrt(G*Mstar/rm**3.0d0) !a Kepler frequency
             H = csm/omega !scale height

          endif
       endif

!**************** Type I migration ******************

       if(type1.ne.0)then
          if(rH.lt.0.84d0*H)then
             If (SIGarray(Tindex,sigindex).gt.0.0d0) Then
              dlogSdlogr = Dlog(SIGarray(Tindex,sigindex+1)) -          &
                     Dlog(SIGarray(Tindex,sigindex-1))
              dlogSdlogr = dlogSdlogr/(DLog(RADarray(sigindex+          &
                     1)) - Dlog(RADarray(sigindex-1)))
             Else
                dlogSdlogr = 0.0d0
             EndIf
             betaMig = -1.0d0*dlogSdlogr

             type1drbdt = (2.7d0+1.1d0*betaMig)*(Mpl/Mstar)             &
                 *(sigmagas(index)*(rm**2.0d0)/Mstar)*(rm*Omega         &
                  /csm)**2.0d0*Omega*rm

             type1drbdt = type1drbdt/type1factor

             type1max30 = (2.7d0+1.1d0*betaMig)*50.0d0*Me/Mstar*       &
                   (sigmagas(index)*(rm**2.0d0)/Mstar)*(rm*Omega/csm)**2.0d0*Omega*rm

             type1dir = 1.0d0 !the default direction is decreasing radius

             if(randomMig.ne.0) then !code for turbulent type 1
                type1time = 0.0d0
                do while (type1time.lt.dt)
                   ranpm = rand(0)*(sigorb*2.0d0)
                   norbit = (10.0d0-sigorb) + ranpm !~10 orbits before random dir chang
                   vel = sqrt(Mstar*G/rm)
                   dist = 2*pi*rm*norbit
                   type1dt = dist/vel
                   norb(m) = norb(m) + norbit
                   ranpm = rannos(int(rand(0)*100000))
                   if(ranpm.gt.0.5) then
                      type1dir = -1.0d0 !50/50 chance for each direction
                   else
                      type1dir = 1.0d0
                   endif

                   If (type1drbdt.lt.type1max30) Then
                     rm = rm - type1dir*type1drbdt*type1dt
                   Else
                     rm = rm - type1dir*type1max30*type1dt -(type1drbdt-type1max30)*type1dt
                   EndIf
                   marker(m,1) = marker(m,1) + sqrt((type1drbdt*type1dt)**2.0d0)
                   type1time = type1time + type1dt
                   if(type1time.gt.dt)then !corrects too much orbiting
                      rm = rm + type1dir*type1drbdt*(type1time-dt)
                      type1time = type1time - (type1time-dt)
                      exit
                   endif
                enddo
             else       !type 1 with no turbulance
                rm = rm - type1dir*type1drbdt*dt
                marker(m,1) = marker(m,1) + type1drbdt*dt
             endif

!     recalculate conditions for new radius
             RH = rm*(Mpl/3.0D0/Mstar)**(1.0D0/3.0D0) !recalculates the hill radius
             Tm = temper(index)
             csm = DSqrt(7.0d0/5.0d0*k*Tm/muH/mH) !
             Omega = DSqrt(G*Mstar/rm**3.0d0) !a Kepler frequency
             H = csm/omega !scale height

             if(rH.ge.0.84d0*H)goto 55

          endif !line416
       endif !line415

       if(rm.gt.r_out) then !removes planets out of bounds
          check(m) = 2
          exit
       elseif(rm.lt.r_in) then
          check(m) = 2
          exit
       endif
55            continue

!     this loop sets "index" to be the number in the radius array where the planet is
       do i = 1, numpts-1
          If ((rad(i).le.rm).and.(rad(i+1).gt.rm)) Then
             index = i
          EndIf
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

       dt = origdt !resets dt to origional value
    end do !ends time loop

17         continue !this continue is for a simulation exit due to too big RH

    finalcoresR(m)=rm/au
    finalcoreM(m)=Mcore
    finalpM(m)=Mpl
!     creates the final output, t,fd,fg,rad,Mgiso,Mplanet,Mcore,ocean(if snow line was crossed),
!     check(the status of planet),hill rad,feeding area,t1 magnitude, t2 magnitude
    write(19,*) time, fdust, fgas, rm/au, Mgiso, Mpl/Me, Mcore/Me, ocean, &
                check(m), rH/Au, af/Au, marker(m,1)/Au, marker(m,2)/Au

       finalMD(m) = 0.0
       do 3 i=1,numpts-1
          finalMD(m)=finalMD(m)+(sigmadust(i)*(pi*rad(i+1)**2-pi*rad(i)**2))
3             continue

       if((Mpl/Mj).gt.masslimit)then
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

       open(unit=33,file='totalt1.dat',status='unknown') !total amaount of t1 migration
       open(unit=34,file='totalt2.dat',status='unknown') !total amaount of t2 migration
       write(33,*) marker(m,1)/Au,' ',Mpl/Me
       write(34,*) marker(m,2)/Au,' ',Mpl/Me

 end do  !ends loop over placing the core at different radii in disc, line 125

!     a loop to count the number of planets exceeding the mass limit for each disc
 numstar = 0
 numlhr = 0
 numok = 0
 do 39 i=1,numpl
    if(check(i).eq.3)then
       numstar = numstar + 1
    elseif(check(i).eq.2)then
       numlhr = numlhr + 1
    elseif(check(i).eq.1)then
       numok = numok + 1
    endif
39      continue
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
 do 74 i=1,numpl
    if(check(i).ne.2.and.check(i).ne.3)then
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
74      continue

 do 22 i=1,numpl
    if(check(i).eq.1)then
       write(12,*) finalcoresR(i), (finalcoreM(i)/Me)
       write(13,*) finalcoresR(i), (finalpM(i)/Me)
       write(27,*) initialcoresR(i), finalpM(i)/Me
    elseif(check(i).ne.1)then
       write(14,*) finalcoresR(i), (finalpM(i)/Me)
       write(15,*) finalcoresR(i), (finalcoreM(i)/Me)
    endif
22      continue

end do !ends loop over discs, line 104

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


end !ends the main program

Subroutine calcFg(Fg,Rc,RH,af,omega)

Real*8 Fg, Rc, RH, af, omega
Real*8 erf
Real*8 rhopl, Mpl, rpl, vescpl
Real*8 i, iH, e, eH, dc
Real*8 S1, S2, S3, beta, B, eps1, eps2
Real*8 AU, G, pi, k, Msun, stefan
External erf
COMMON /const/ AU, G, pi, k, Msun, stefan

rpl = 1.0D7
rhopl = 3.2D0
Mpl = rhopl*4.0D0/3.0D0*pi*rpl**3.0D0

vescpl = DSqrt(2.0D0*G*Mpl/rpl)

iH = vescpl/Sqrt(3.0D0)/omega/RH

eH = max(2.0D0*iH,2.0D0)

af = DSqrt(12.0D0 + eH**2.0D0)*RH

dc = Rc/RH

eps1 = 0.2D0*iH


if (dc.lt.0.0284D0) then
 eps2 = Sqrt(2.2D0*dc)/iH + eps1
else
 eps2 = 0.25D0/iH + eps1
endif


S1 = eps1*Sqrt(pi)*(erf(eps1)-erf(eps2))
S1 = Exp(-1.0*eps1**2.0)-Exp(-1.0*eps2**2.0)+S1
S1 = 4.0*Exp(eps1**2.0)/dc**1.5*S1

If (dc.lt.0.0284) then
 S2 = erf(0.25/iH)-erf(Sqrt(2.2*dc)/iH)
 S2 = 6.69/dc/iH*S2
Else
 S2 = 0.0
EndIf


S3 = 1.0 + 3.31/dc/iH**2.0

beta = 1.0/6.0*((log10(iH)-0.131)/0.343)**3.0
beta = (Dlog10(iH)-0.131D0)/0.343D0 + beta
beta = 1.0D0 + Exp(beta)
beta = pi/2.0D0/beta

B = 1.0D0/24.0D0*((log10(iH)+0.109D0)/0.231D0)**4.0D0
B = 1.0D0 + 0.422D0*Exp(-1.0D0*B)

Fg = B*((S1+S2)*Sin(beta)+S3*cos(beta))

Return
End



Function erf(x)

Real*8 erf, x
Real*8 gammp

if (x.lt.0.0D0) then
 erf = -gammp(0.5D0,x**2.0D0)
else
 erf = gammp(0.5D0,x**2.0D0)
endif

return
end



FUNCTION GAMMP(A,X)

REAL*8 A,GAMMP,X
REAL*8 GAMMCF, GAMSER, GLN

IF(X.LT.0.0D0.OR.A.LE.0.0D0)STOP

IF(X.LT.A+1.0D0)THEN
 CALL GSER(GAMSER,A,X,GLN)
 GAMMP=GAMSER
ELSE
 CALL GCF(GAMMCF,A,X,GLN)
 GAMMP=1.-GAMMCF
ENDIF

RETURN
END



SUBROUTINE GSER(GAMSER,A,X,GLN)

INTEGER ITMAX, N
REAL*8 A,GAMSER,GLN,X,EPS
REAL*8 AP,DEL,SUN,GAMMLN
PARAMETER (ITMAX=100,EPS=3.D-7)
EXTERNAL GAMMLN
GLN=GAMMLN(A)

IF(X.LE.0.0D0)THEN
 IF(X.LT.0.0D0)STOP
 GAMSER=0.0D0
 RETURN
ENDIF

AP=A
SUM=1.0D0/A
DEL=SUM

DO 11 N=1,ITMAX
 AP=AP+1.0D0
 DEL=DEL*X/AP
 SUM=SUM+DEL
 IF(ABS(DEL).LT.ABS(SUM)*EPS)GO TO 1
11   CONTINUE
STOP 'A too large, ITMAX too small'
1    GAMSER=SUM*EXP(-X+A*DLOG(X)-GLN)

RETURN
END




SUBROUTINE GCF(GAMMCF,A,X,GLN)

INTEGER ITMAX
REAL*8 A,GAMMCF,GLN,X,EPS,FPMIN
REAL*8 AN,B,C,D,DEL,H,GAMMLN
PARAMETER (ITMAX=100,EPS=3.E-7)
EXTERNAL GAMMLN
GLN=GAMMLN(A)
GOLD=0.0D0
A0=1.0D0
A1=X
B0=0.0D0
B1=1.0D0
FAC=1.0D0

DO 11 N=1,ITMAX
 AN=FLOAT(N)
 ANA=AN-A
 A0=(A1+A0*ANA)*FAC
 B0=(B1+B0*ANA)*FAC
 ANF=AN*FAC
 A1=X*A0+ANF*A1
 B1=X*B0+ANF*B1

 IF(A1.NE.0.0D0)THEN
    FAC=1./A1
    G=B1*FAC
    IF(ABS((G-GOLD)/G).LT.EPS)GO TO 1
    GOLD=G
 ENDIF
11   CONTINUE
STOP 'A too large, ITMAX too small'
1    GAMMCF=DEXP(-X+A*DLog(X)-GLN)*G

RETURN
END




FUNCTION GAMMLN(XX)

REAL*8 GAMMLN,XX
REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,-1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/
X=XX-ONE
TMP=X+FPF
TMP=(X+HALF)*DLOG(TMP)-TMP
SER=ONE

DO 11 J=1,6
 X=X+ONE
 SER=SER+COF(J)/X
11   CONTINUE
GAMMLN=TMP+DLOG(STP*SER)

RETURN
END


!
!     *******************************************************************************
!  Black box random number generator
! *******************************************************************************
!
SUBROUTINE RANSET(ir,iseed)
!
!     This routine is used to initialize all of the
!     random number generators in RANPAK.
!
!     The seeding is performed as determined by ir
!     and iseed as follows:
!
!        ir <= 0, initialize all seeds
!        ir > 0, initialize only seed ir
!        iseed < 1, the default seed is used
!        iseed >= 1, use the value passed in iseed
!
!     After initialization, the value of the jth
!     random number seed should be accessed with
!     the RANGET routine.
!
!     CODE DEPENDENCIES: (internal RANPAK COMMON block)
!
!     DATE: DEC. 1, 1994
!     AUTHOR: R.D. STEWART
!
INTEGER ir,iseed,iset,i
INTEGER nrg
PARAMETER (nrg=100)
INTEGER isd(nrg)
COMMON /RanPak/ isd

IF (iseed.GT.0) THEN
 iset = iseed
ELSE
!     Use the "demonic" default seed
 iset = 666
ENDIF

IF (ir.LT.1) THEN
!     Initialize all random number generator seeds
!     to iseed
 DO i=1,nrg
    isd(i) = iset
 ENDDO
ELSEIF (ir.LE.nrg) THEN
 isd(ir) = iset
ENDIF

ir = 1

RETURN
END





SUBROUTINE DRAN(ir,rand_no)
!
!     PREFERRED DOUBLE PRECISION RANDOM NUMBER GENERATOR
!
!     The random number seed is set initially by a call to the
!     RANSET routine.  ir is a pointer to an internal random
!     number seed and NOT a random number seed.
!
!     COMMENTS: This routine will work correctly on any computer
!               with a maximum integer greater than or equal
!               to 2**31 - 1.  NOTE: The proper function of
!               the routine can verified by checking to see
!               that the random number seed after the generation
!               of 10000 random numbers is 1,043,618,065
!
!               For more information refer to the classic paper
!
!               Park, S.K. and Miller, K.W., "Random number generators:
!               good ones are hard to find."  Communications of the
!               ACM, Vol 31, No 10 (Oct. 1988).
!
!     CODE DEPENDENCIES: (internal RANPAK COMMON block)
!
!     DATE: DEC. 1, 1994
!     AUTHOR: R.D. STEWART
!
INTEGER nrg,ir
PARAMETER (nrg=100)
INTEGER isd(nrg)
COMMON /RanPak/ isd

DOUBLE PRECISION tmp,a,seed,rand_no
DOUBLE PRECISION m,Minv
PARAMETER (a=16807.0D+00,m=2147483647.0D+00)
PARAMETER (Minv=1.0D+00/m)

seed = isd(ir)

tmp = A*seed
seed = tmp - M*DINT(tmp*Minv)
rand_no = seed*Minv

isd(ir) = seed

RETURN
END
