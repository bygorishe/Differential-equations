      real function f(t,y)
      implicit none
      real t,y 
      f=2*t*y
      end
      
      real function g(t)   !analiticheski y'=2ty  y(t)=Ce^(t*t)
      implicit none
      real t 
      g=exp(t*t)                  !y(0)=1   y(t)=e^t^2
      end
      
      real function Fx(t,y,yy,h)
      implicit none
      real t,y,yy,h,f
      Fx=y-yy-h*f(t,y)
      end
      
      real function Fxx(t,y,yy,h)
      implicit none
      real t,y,yy,h,f
      Fxx=y-yy-h/2*(f(t,y)+f(t-h,yy))
      end
      
      real function Gx(t,h)
      implicit none
      real t,h
      Gx=1-h*2*t     !Fx/dx
      end
      
      real function Gxx(t,h)
      implicit none
      real t,h
      Gxx=1-h*t         !Fxx/dx
      end
      
      subroutine newton(t,h,m1)
      implicit none
      common/var/y,ypr
      real t,y,yy,ypr,h,Fx,Fxx,Gx,Gxx  
      integer m1
      yy=0
      if(m1.eq.1)then
       do while(abs(y-yy).gt.0.00001)  
        yy = y 
        y=y-Fx(t,y,ypr,h)/Gx(t,h)
        !write(*,*) y
       end do
      endif
      if(m1.eq.2)then
       do while(abs(y-yy).gt.0.00001)  
        yy = y 
        y=y-Fxx(t,y,ypr,h)/Gxx(t,h)
        !write(*,*) y
       end do
      endif
      end

      subroutine EL(h,t)   !Method Eulera
      implicit none
      common/var/y,ypr 
      real h,t,y,ypr
      real f
       y=ypr+h*f(t,y)    
      end
      
      subroutine ELM(h,t)  !Method Eulera mod
      implicit none
      common /var/ y,ypr 
      real h,t,y,ypr
      real f
       y=ypr+h/2*(f(t-h,ypr)+f(t,y))
      end

      program main
      implicit none
      common/var/y,ypr 
      integer m1,m2,n,i
      real t,y,ypr,h,g
      y=1.0
      t=0.
      open(1,FILE='out.txt')
      write(1,70)
      write(1,80)
      write(1,70)
2     print*, 'Select:'
      print*, '1) h=0.1'
      print*, '2) h=0.05'
      print*, '3) h=0.025'
      read(*,*) m2
      if(m2.eq.1)then
       h=0.1
       goto 11
      endif
      if(m2.eq.2)then
       h=0.05
       goto 11      
      endif   
      if(m2.eq.3)then
       h=0.025
       goto 11      
      endif 
11    print*, 'Select method:'
      print*, '<1> el'
      print*, '<2> elmod'   
      read(*,*) m1  
      n=1/h+1 
      if(m1.eq.1)then
       do i=0,n
        write(*,300) t,y,g(t),abs(y-g(t))
        write(1,100) t,y,g(t),abs(y-g(t))
        write(1,70)
        t=t+h
        ypr=y
        call newton(t,h,m1)
        call EL(h,t)
       enddo 
      endif
      if(m1.eq.2)then
       do i=0,n
        write(*,300) t,y,g(t),abs(y-g(t))
        write(1,100) t,y,g(t),abs(y-g(t))
        write(1,70)
        t=t+h
        ypr=y
        call newton(t,h,m1)
        call ELM(h,t)
       enddo
      endif 
      print*, 'Completed.'
      close(1)
      pause  
        
   70 format(52('-'))   
   80 format(' t          | y          | y*         | |y-y*|     |')  
  300 format('t = 'E11.3 ' y = 'E16.8 ' y* = 'E16.8 ' |y-y*| = 'E16.8)
  100 format(E11.3' |'E11.4' |'E11.4' |'E11.4' |')
      end
