      subroutine force(md,x,b)
      implicit double precision (a-h,o-z)
C calculate force subroutione by Chen Xin
C df/dx=(f(x+h)-f(x-h))/(2h) 
C set h=dx/10.0 
C Input x: coordinate of the image
C output b: potential force
C 
      parameter (nx=100,ny=100)
      double precision x(md),b(md)
      double precision hx,hy
      double precision e1,e2
      double precision x1_PES(nx),x2_PES(ny),y_PES(nx,ny)
      common /PES/x1_PES,x2_PES,y_PES
C
C
      b(1)=0.d0
      b(2)=0.d0
      hx=(x1_PES(2)-x1_PES(1))/10.0
      hy=(x2_PES(2)-x2_PES(1))/10.0
C      print *,x(1),x(2)
      call SPLINE_2d(x1_PES,x2_PES,y_PES,nx,ny,x(1)+hx,x(2),e1)
      call SPLINE_2d(x1_PES,x2_PES,y_PES,nx,ny,x(1)-hx,x(2),e2)
      b(1)=-(e1-e2)/(2*hx)
      call SPLINE_2d(x1_PES,x2_PES,y_PES,nx,ny,x(1),x(2)+hy,e1)
      call SPLINE_2d(x1_PES,x2_PES,y_PES,nx,ny,x(1),x(2)-hy,e2)
      b(2)=-(e1-e2)/(2*hy)
C      print *,x(1),x(2),b(1),b(2),e1


C        
      return
      end
