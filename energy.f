      subroutine energy(x,e)
      implicit double precision (a-h,o-z)
C
C by Chen Xin @JLU
C Calculate energy or image by spline interplorate
C Input : x(2,0:n), the coordinates of each image along the string
C Output : e(0:n), the potential energy of each image
C
      parameter (md=2,n=20)
      parameter (nx=100, ny=100)
C
      double precision x(md,0:n),e(0:n)
      double precision x1_PES(nx),x2_PES(ny),y_PES(nx,ny)
      common /PES/x1_PES,x2_PES,y_PES
C---------------------------------------------------
C    Energy grid input zone 
C    make sure that axix input with as ascend array
C    
C--------------------------------------------------- 
      do i=0,n
        call SPLINE_2d(x1_PES,x2_PES,y_PES,nx,ny,x(1,i),x(2,i),e(i))
C        print *, x(1,i),x(2,i),e(i)
      end do
C
      end
