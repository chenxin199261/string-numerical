      subroutine mvtan(x)
      implicit double precision (a-h,o-z)
      parameter (md=2,n=20)
C
C Reparameterization of the string by polynomial interpolation
C By W. Ren
C
      double precision x(md,0:n),s(0:n),s1(0:n),x1(md,0:n)
C
C calculate the distance between neighboring images along the string
      do 1 k=1,n
        s(k)=0.d0
        do 2 i=1,md
          s(k)=s(k)+(x(i,k)-x(i,k-1))**2.d0
  2     continue
        s(k)=dsqrt(s(k))
  1   continue
C
C normaized arclength starting from the initial state
      s(0)=0.d0
      do 3 k=1,n
        s(k)=s(k-1)+s(k)
  3   continue
      do 31 k=1,n
        s(k)=s(k)/s(n)
  31  continue
C
C Equal-arclength: the locations at which new images are calculated
C by interpolation 
      do 4 k=0,n
        s1(k)=dble(k)/dble(n)
  4   continue
C
C interpolation to get new images
C
      call intpol(md,n,s,x,n,s1,x1)
C
      call dcopy(md*(n+1),x1,1,x,1)
C
      return
      end
