program core_accr

    use constants_mod
    use disc_mod
    use everything_mod
    use migration_mod

    implicit none

    integer m

    ! open random numbers file
    open(unit=25,file='randnums.dat',status='old')
    read(25,*) rannos
    close(25)

    ! Read params
    open(unit=19,file='core_accr.params',status='old')
    read(19,*) rm ! semi-major axis of planet
    read(19,*) Mcore ! Initial core mass of planet
    read(19,*) sigma5au ! sigma_p @ 5au, used as scaling factor
    close(19)

    rm = rm*AU         ! initial core semi-major axis
    Mcore = Mcore*Me   ! initial core mass
    Mstar = 2.4d0*Msun ! star mass

    ! Read in evolution of surface density profile of gas disc
    if (Mstar.eq.Msun) then
        open(unit=26,file='discMstar1.dat',status='old')
        print*,'Mstar = 1.0Msun'
    else if(Mstar.eq.real(2.4d0*Msun))then
        open(unit=26,file='discABAur.dat',status='old')
        print*,'Mstar = 2.4Msun'
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

    tolD = 0.02 ! limits the dust accretion at each timestep, dmcore (*Me)
    tolG = 0.02 ! limits the gas accretion at each timestep, dmgas (*Mj)
    numpl = 1 ! number of different radii to try for each disc
    numdisc = 1 ! number of stars to simulate

    do l = 1, numdisc ! loop over discs
        nex = 0 ! counts number of exits due to negligible mass gain
        ! a loop to set all planets as being initially in the simulation
        do i=1,numpl
            check(i) = 1 !1=ok, 2=rh>disc, 3=mass limit exceeded
        enddo

!*************************************************************************
!                         setup procedures
!*************************************************************************
        print*,' '
        print*,'Running disc number ',l,'/',numdisc

        ! loop over number of simulated planets
        do m = 1,numpl
            call setup(m)

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
43              continue

!*************************************************************************
!         now calculate planet accretion for this timestep
!*************************************************************************
                ! checks the planet hasnt exceeded its isolation mass or mass limit
                if ((Mpl/Me.lt.Mgiso).and.(Mpl/Mj.lt.massLimit)) then
                    rhom = gasSig/H ! gas density at planet location
                    rhocore = 3.2d0 ! density of core, gcm^-3
                    Rcore = (3.0d0*Mcore/(4.0d0*pi*rhocore))**(1.0d0/3.0d0)
                    Rc = Rcore
                    RH = rm*(Mpl/3.0d0/Mstar)**(1.0d0/3.0d0) ! the hill radius

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
                    !prevents unphisical Mdust < 0
                    if (Mdust.lt.0.0d0) then
                        Mdust = 0.0d0
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
                    ! Kelvin-Helmholtz timescale for gas envelope accretion
                    TauKH = 1.0d9*((Mpl/Me)**(-3.0d0))*kappa
                    ! critical mass for gas accreation
                    Mcrit = 10.0d0*(Mdotcore*yr/Me/1.0d-6)**(1.0d0/4.0d0)*kappa**(1.0d0/4.0d0)
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
                        ! otherwise gas can escape
                        Mdotgas = 0.0d0
                    endif
                    ! calculates the gas isolation mass at the current timestep, with decaying sigmagas
                    Mgiso = 50.0d0*(sigmagas(index)/2.4d3)**1.5d0*(rm/AU)**3.0*(Mstar/Msun)**(-0.5d0)
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
                    sigmadust(index1)=(sigmadust(index1)*                              &
                                        (pi*(rm-af)**2.0D0-pi*rad(index1)**2.0D0)+     &
                                        sigmadust_temp*                                &
                                        (pi*rad(index1+1)**2.0D0-pi*(rm-af)**2.0D0))/  &
                                        (pi*rad(index1+1)**2.0D0-pi*rad(index1)**2.0D0)
                    ! all the sigmadust inside hill sphere of planet are set equal
                    do j = index1+1, index2-2
                         sigmadust(j) = sigmadust_temp
                    enddo
                    ! recalculates the outside edge sigmadust
                    sigmadust(index2-1)=(sigmadust(index2-1)*                      &
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
                call calc_migration(m)

                ! remove planet out of bounds
                if (rH.ge.0.84d0*H) goto 55
                if (rm.gt.rOut) then
                    check(m) = 2
                    exit
                elseif (rm.lt.rIn) then
                    check(m) = 2
                    exit
                endif
55              continue

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
17          continue

            ! store final values
            finalcoresR(m)=rm/au
            finalcoreM(m)=Mcore
            finalpM(m)=Mpl

            open(unit=19,file='final.dat')
            write(19,*) time, fgas, rm/au, Mgiso, Mpl/Me, Mcore/Me, ocean, &
                        check(m), rH/Au, af/Au, marker(m,1)/Au, marker(m,2)/Au
            close(19)

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

        enddo ! end num planets loop

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
        print*,'planets exceeding mass limit=',numstar
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

        write(17,*) initialMD(1)/Mstar
        write(18,*) Mgas(1)/Mstar

        do i=1,numpl
            if (check(i).eq.1) then
                write(12,*) finalcoresR(i), (finalcoreM(i)/Me)
                write(13,*) finalcoresR(i), (finalpM(i)/Me)
                write(27,*) initialcoresR(i), finalpM(i)/Me
            elseif (check(i).ne.1) then
                write(14,*) finalcoresR(i), (finalpM(i)/Me)
                write(15,*) finalcoresR(i), (finalcoreM(i)/Me)
            endif
        enddo

    enddo ! ends loop over discs

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
    close(27)
    close(28)
    close(29)
    close(30)
    close(31)
    close(32)
    close(33)
    close(34)

! end the main program
end program core_accr
