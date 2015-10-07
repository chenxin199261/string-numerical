      program main
C
C generate initial path for Mueller potential
C by W. Ren
C
      parameter (n=20)
      double precision x(2,0:n),t
C
C initial state
      x(1,0)= -1
      x(2,0)= 0
C
C final state
      x(1,n)= 1
      x(2,n)= 0
C
C linear interpolation between the initial and final state
      do 2 i=1,n-1
        t=dble(i)/dble(n)
        x(1,i)=x(1,0)*(1.d0-t)+x(1,n)*t
        x(2,i)=x(2,0)*(1.d0-t)+x(2,n)*t
  2   continue
C
C output the path
      call output('data0',2,n,x)
C
      end
