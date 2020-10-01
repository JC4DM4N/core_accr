subroutine calc_growth(m,ierr)

    use disc_mod
    use constants_mod
    use planet_mod
    use setup_mod

    integer, intent(in)    :: m
    integer, intent(inout) :: ierr

    ierr = 0

    rhom = gasSig/H ! gas density at planet location
    rhocore = 3.2d0 ! density of core, gcm^-3
    Rcore = (3.0d0*Mcore/(4.0d0*pi*rhocore))**(1.0d0/3.0d0)
    Rc = Rcore
    RH = rm*(Mpl/3.0d0/Mstar)**(1.0d0/3.0d0) ! Hill radius
    ! calculate the grav enhancement factor
    call calcFg(Fg,Rc,RH,af,omega)
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
        ierr = 1
        return
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
end subroutine
