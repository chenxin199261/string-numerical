      subroutine intpol(md,n0,xi0,x0,n1,xi1,x1)
      implicit double precision (a-h,o-z)
C
C third order polynomial interpolation
C by W. Ren
      parameter (mi=4)
C
C Reparameterization of the string
C Input: xi0, x0: normalized arclength and coordinates of 
C                 the images along the string before reparameterzation
C        xi1: normalized arclength at which interpolation is done
C
C output: x1: images along the string after interpolation
C
      double precision xi0(0:n0),x0(md,0:n0),xi1(0:n1),x1(md,0:n1)
      double precision xx(mi),dy
C
      do 1 i=1,n1-1
        call hunt(xi0,n0+1,xi1(i),jlo)
        kk=min(max(jlo-(mi-1)/2,1),n0+2-mi)-1
        do 2 k=1,md
          do 3 j=1,mi
            xx(j)=x0(k,kk+j-1)
  3       continue
          call polint(xi0(kk),xx(1),mi,xi1(i),x1(k,i),dy)
  2     continue
  1   continue
C
C the initial and final states are fixed
      do 4 k=1,md
        x1(k,0)=x0(k,0)
        x1(k,n1)=x0(k,n0)
  4   continue
C
      return
      end
      
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      implicit double precision (a-h,o-z)
      INTEGER n,NMAX
      double precision dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      double precision den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END

      SUBROUTINE hunt(xx,n,x,jlo)
      implicit double precision (a-h,o-z)
      INTEGER jlo,n
      double precision x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)then
        if(x.eq.xx(n))jlo=n-1
        if(x.eq.xx(1))jlo=1
        return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END
