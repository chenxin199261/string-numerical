      subroutine tangent(x,t)
      implicit double precision (a-h,o-z)
C
C calculate the unit tangent vectors along the string
C by W. Ren
C
C Input x: coordinates of the images along the string
C Output t: unit tangent vectors 
C
      parameter (md=2,n=20)
      double precision x(md,0:n),t(md,n-1)
      double precision e(0:n),e1,e2,tmp1,tmp2
C
C upwind scheme based on the potential energy
C
      call energy(x,e)
C
      do 1 k=1,n-1
        if(e(k+1).ge.e(k).and.e(k).ge.e(k-1)) then
          do 2 i=1,md
            t(i,k)=x(i,k+1)-x(i,k)
  2       continue
        else if(e(k+1).le.e(k).and.e(k).le.e(k-1)) then
          do 3 i=1,md
            t(i,k)=x(i,k)-x(i,k-1)
  3       continue
        else
          tmp1=dabs(e(k-1)-e(k))
          tmp2=dabs(e(k)-e(k+1))
          e1=dmax1(tmp1,tmp2)
          e2=dmin1(tmp1,tmp2)
          if(e(k+1).ge.e(k-1)) then
            do 4 i=1,md
              t(i,k)=e1*(x(i,k+1)-x(i,k))+
     .             e2*(x(i,k)-x(i,k-1))
  4         continue
          else
            do 5 i=1,md
              t(i,k)=e2*(x(i,k+1)-x(i,k))
     .             +e1*(x(i,k)-x(i,k-1))
  5         continue
          endif
        endif
C nomalization
        tmp1=0.d0
        do 6 i=1,md
          tmp1=tmp1+t(i,k)**2.d0
  6     continue
        tmp1=dsqrt(tmp1)
        do 7 i=1,md
          t(i,k)=t(i,k)/tmp1
  7     continue
C
  1   continue
C
      return
      end
