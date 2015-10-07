      subroutine input(fname,md,n,x)
      integer md,n
      double precision x(md,0:n)
      character*20,fname
      open(11,file=fname,status='unknown',form='formatted')
        do 1 k=0,n
          read(11,*)x(1,k),x(2,k)
C          print *,x(1,k),x(2,k)
  1     continue
      close(11)
      return
      end

      subroutine output(fname,md,n,x)
      integer md,n
      double precision x(md,0:n)
      character*20,fname
      open(11,file=fname,status='unknown',form='formatted')
        do 1 k=0,n
          write(11,*)x(1,k),x(2,k)
  1     continue
      close(11)
      return
      end
      
      subroutine readPES(fname)
C     read PES in subroutine by ChenXin
C
C
C     nx: number of x grid and ny:number in y grid
      parameter(nx=100,ny=100)
      integer i,j
      character*20,fname
      double precision x1_PES(nx),x2_PES(ny),y_PES(nx,ny)

      common /PES/x1_PES,x2_PES,y_PES

      open(11,file=fname,status='unknown',form='formatted')
      do i=1,nx
         do j=1,ny
            read(11,*) x1_PES(i),x2_PES(j),y_PES(i,j)
         end do
      end do
      close(11)
      end
        
