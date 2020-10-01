subroutine calc_migration(m)

    use constants_mod
    use disc_mod
    use migration_mod
    use planet_mod
    use setup_mod

    integer, intent(in) :: m

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
        csm = sqrt(7.0d0/5.0d0*k*Tm/muH/mH)
        Omega = sqrt(G*Mstar/rm**3.0d0)
        H = csm/omega
      endif
    endif

    ! type I migration
    if (type1.ne.0) then
      if (rH.lt.0.84d0*H) then
        if (SIGarray(Tindex,sigindex).gt.0.0d0) then
          dlogSdlogr = log(SIGarray(Tindex,sigindex+1)) -            &
                       log(SIGarray(Tindex,sigindex-1))
          dlogSdlogr = dlogSdlogr/(Log(RADarray(sigindex+            &
                       1)) - log(RADarray(sigindex-1)))
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
            ranpm = rand*(sigorb*2.0d0)
            norbit = (10.0d0-sigorb) + ranpm !~10 orbits before random dir chang
            vel = sqrt(Mstar*G/rm)
            dist = 2*pi*rm*norbit
            type1dt = dist/vel
            norb(m) = norb(m) + norbit
            ranpm = rannos(int(rand*100000))
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

       endif
    endif
    return
end subroutine calc_migration
