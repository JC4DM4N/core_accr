subroutine calcFg(Fg,Rc,RH,af,omega)

  real*8, intent(inout) :: Fg, af
  real*8, intent(in) :: Rc, RH, omega
  real*8 erf
  real*8 rhopl,Mpl,rpl,vescpl
  real*8 iH,eH,dc
  real*8 S1,S2,S3,beta,B,eps1,eps2
  real*8 AU,G,pi,k,Msun,stefan
  external erf
  common /const/ AU,G,pi,k,Msun,stefan

  rpl = 1.0d7
  rhopl = 3.2d0
  Mpl = rhopl*4.0d0/3.0d0*pi*rpl**3.0d0
  vescpl = sqrt(2.0d0*G*Mpl/rpl)
  iH = vescpl/sqrt(3.0d0)/omega/RH
  eH = max(2.0d0*iH,2.0d0)
  af = sqrt(12.0d0 + eH**2.0d0)*RH
  dc = Rc/RH
  eps1 = 0.2d0*iH

  if (dc.lt.0.0284d0) then
    eps2 = sqrt(2.2d0*dc)/iH + eps1
  else
    eps2 = 0.25d0/iH + eps1
  endif

  S1 = eps1*sqrt(pi)*(erf(eps1)-erf(eps2))
  S1 = exp(-1.0*eps1**2.0)-exp(-1.0*eps2**2.0)+S1
  S1 = 4.0*exp(eps1**2.0)/dc**1.5*S1

  if (dc.lt.0.0284) then
    S2 = erf(0.25/iH)-erf(sqrt(2.2*dc)/iH)
    S2 = 6.69/dc/iH*S2
  else
    S2 = 0.0
  endif
  S3 = 1.0 + 3.31/dc/iH**2.0

  beta = 1.0/6.0*((log10(iH)-0.131)/0.343)**3.0
  beta = (log10(iH)-0.131d0)/0.343d0 + beta
  beta = 1.0d0 + exp(beta)
  beta = pi/2.0d0/beta

  B = 1.0D0/24.0D0*((log10(iH)+0.109D0)/0.231D0)**4.0D0
  B = 1.0D0 + 0.422D0*exp(-1.0D0*B)

  Fg = B*((S1+S2)*sin(beta)+S3*cos(beta))

  return
end

function erf(x)

  real*8 erf, x
  !real*8 gammp

  if (x.lt.0.0d0) then
    erf = -gammp(0.5,x**2.0)
  else
    erf = gammp(0.5,x**2.0)
  endif

  return
end

function gammp(a,x)
  real*8 a,gammp,x
  real*8 gammcf, gamser, gln

  if(x.lt.0.0d0.or.a.le.0.0d0)stop
  if(x.lt.a+1.0d0) then
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
  real*8 ap,del,gammln
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
  1    gamser=sum*exp(-x+a*log(x)-gln)

  return
end

subroutine gcf(gammcf,a,x,gln)

  integer itmax
  real*8 a,gammcf,gln,x,eps
  real*8 an,gammln
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
  1    gammcf=exp(-x+a*log(x)-gln)*g

  return
end

function gammln(xx)
  real*8 gammln,xx
  real*8 cof(6),stp,half,one,fpf,x,tmp,ser
  data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,-1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
  data half,one,fpf/0.5d0,1.0d0,5.5d0/
  x=xx-one
  tmp=x+fpf
  tmp=(x+half)*log(tmp)-tmp
  ser=one

  do j=1,6
    x=x+one
    ser=ser+cof(j)/x
  enddo
  gammln=tmp+log(stp*ser)

  return
end
