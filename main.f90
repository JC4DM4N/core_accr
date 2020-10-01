program core_accr

    use constants_mod
    use disc_mod
    use migration_mod
    use planet_mod
    use setup_mod

    implicit none

    integer m
    integer :: ierr = 0

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

    call read_discfile()

    origdt = 10000.0d0*yr
    disclife = Tarray(1000)/yr
    massLimit = 15.0 !max mass is 15 jupiter masses

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
                    call calc_growth(m,ierr)
                    if (ierr.eq.1) goto 17
                endif

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
                ! write planet data to files for this timestep
                if(numpl.eq.1.and.numdisc.eq.1)then
                    open(unit=26,file='migration.dat',status='unknown')
                    write(26,*) time,' ',rm/Au
                    open(unit=28,file='mvt.dat',status='unknown')
                    write(28,*) time,' ',mcore/me,' ',mpl/me
                    open(unit=29,file='mvr.dat',status='unknown')
                    write(29,*) rm/Au,' ',mpl/me
                    open(unit=30,file='hvr.dat',status='unknown')
                    write(30,*) rm/Au,' ',H/Au
                    open(unit=32,file='hvrh.dat',status='unknown')
                    write(32,*) rm/Au,' ',rH/Au
                endif

                ! evolve the time
                time = time + dt/yr

                ! exit if accretion has stopped
                if (dMcore/me.lt.1d-30.and.dt*Mdotgas/me.lt.1d-6) then
                    print*,'exit, growth',m,' negligible'
                    nex = nex + 1
                    goto 17
                endif
               ! resets dt to original timestep
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
            write(19,*) time, rm/au, Mgiso, Mpl/Me, Mcore/Me, ocean, &
                        check(m), rH/Au, af/Au, marker(m,1)/Au, marker(m,2)/Au
            close(19)

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

        nok = nok + numok
        nob = nob + numlhr
        nml = nml + numstar


    enddo ! ends loop over discs

    print*,' '
    print*,'total number of planets',nok
    print*,'total number out of bounds',nob
    print*,'total munber over the mass limit',nml


    ! files ive checked
    close(26)
    close(28)
    close(29)
    close(30)
    close(32)
    close(33)
    close(34)

! end the main program
end program core_accr
