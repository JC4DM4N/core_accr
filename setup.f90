subroutine setup(m)

    use disc_mod
    use constants_mod
    use everything_mod

    integer, intent(in) :: m

    print*, '-------------------'
    print*, 'begin disc setup...'
    print*, ' '

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
      SIGarray(i,j) = OSA(i,j)
    enddo
    enddo

    ! calculate sigmagas, extrapolated between values in SIGarray
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
    sigmadust(i) = etaIce*(1.0d0*AU/rad(i))**(beta)
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
    ! define initial gas isolation mass
    Mgiso = 50.0d0*(sigmagas(index)/2.4d3)**1.5d0*(rm/AU)**3.0*(Mstar/Msun)**(-0.5d0)

    time = 0.0d0
    ! keeps track of what dt is origionally set as
    dt = origdt

    ! initial core migration is zero
    drbdt = 0.0d0
    ! initial gas mass of planet is zero
    Mgas(m) = 0.0d0
    do i=1,numpts-1
      Mgas(m)=Mgas(m)+(sigmagas(i)*(pi*rad(i+1)**2-pi*rad(i)**2))
    enddo

    print*, 'Mgas=',Mgas(m)/Msun,' Msol, Mdust=',initialMD(m)/Msun,' Msol'
    print*, 'Current dust-to-gas ratio=', initialMD(m)/Mgas(m)
    ! Assumes star-disc intially had q=0.5
    print*, 'Initial dust-to-gas ratio=', 2.0d0*initialMD(m)/Mstar
    print*, ' '
    print*, '...setup complete'
    print*, '-----------------'
    print*, ' '

    return
end subroutine setup
