      real*8 function f(t,y)
      implicit none
      real*8 t,y 
      f=2*t*y
      end
      
      real*8 function g(t)   !analiticheski y'=2ty  y(t)=Ce^(t*t)
      implicit none
      real*8 t 
      g=exp(t*t)                  !y(0)=1   y(t)=e^t^2
      end
      
      real*8 function ff(t,y)
      implicit none
      real*8 t,y 
      ff=-25*y+cos(t)+25*sin(t)
      end
      
      real*8 function gg(t)   !analiticheski y'=-25*y+cos(t)+25*sin(t)  y(t)=Ce^(-25t)+sin(t)
      implicit none
      real*8 t 
      gg=exp(-25*t)+sin(t)                  !y(0)=1   y(t)=e^(-25t)+sin(t)
      end
      
      real*8 function Fx(t,y,yy,h)
      implicit none
      real*8 t,y,yy,h,ff
      Fx=y-yy-h*ff(t,y)
      end
      
      real*8 function Fxx(t,y,yy,h)
      implicit none
      real*8 t,y,yy,h,ff
      Fxx=y-yy-h/2*(ff(t,y)+ff(t-h,yy))
      end
      
      real*8 function Gx(h)
      implicit none
      real*8 h
      Gx=1+h*25     !Fx/dx
      end
      
      real*8 function Gxx(h)
      implicit none
      real*8 h
      Gxx=1+h*25/2         !Fxx/dx
      end
      
      subroutine newton(t,h,m1)
      implicit none
      common/var/y,ypr
      real*8 t,y,yy,ypr,h,Fx,Fxx,Gx,Gxx  
      integer m1
      yy=0.0D0
      if(m1.eq.2)then
       do while(abs(y-yy).gt.0.00001)  
        yy = y 
        y=y-Fx(t,y,ypr,h)/Gx(h)
        write(*,*) y
       end do
      endif
      if(m1.eq.3)then
       do while(abs(y-yy).gt.0.00001)  
        yy = y 
        y=y-Fxx(t,y,ypr,h)/Gxx(h)
        write(*,*) y
       end do
      endif
      end

      subroutine ELun(h,t)   !Method Eulera un
      implicit none
      common/var/y,ypr 
      real*8 h,t,y,ypr,ff
       y=ypr+h*ff(t,y)    
      end
      
      subroutine TRAPEZE(h,t)  !Method tr
      implicit none
      common /var/ y,ypr 
      real*8 h,t,y,ypr,ff
       y=ypr+h/2*(ff(t-h,ypr)+ff(t,y))
      end
      
      
      real*8 function EL(h,t,y)   !Method Eulera
      implicit none
      real*8 h,t,y,ff
       EL=y+h*ff(t,y)
      end
      
      real*8 function ELM(h,t,y)  !Method Eulera mod
      implicit none
      real*8 h,t,y,ff
       ELM=y+h/2*(ff(t,y)+ff((t+h),(y+h*ff(t,y))))
      end

      subroutine RungeKutta(h,t)
      implicit none
      common /var/ y 
      real*8 h,t,y,f,k1,k2,k3,k4
      k1=f(t,y)
      k2=f(t+h/2,y+h/2*k1)
      k3=f(t+h/2,y+h/2*k2)
      k4=f(t+h,y+h*k3)
      y=y+h*(k1+2*k2+2*k3+k4)/6
      end
      
      subroutine RungeKutta2(h,t) !for -25*y+cos(t)+25*sin(t)
      implicit none
      common /var/ y 
      real*8 h,t,y,ff,k1,k2,k3,k4
      k1=ff(t,y)
      k2=ff(t+h/2,y+h/2*k1)
      k3=ff(t+h/2,y+h/2*k2)
      k4=ff(t+h,y+h*k3)
      y=y+h*(k1+2*k2+2*k3+k4)/6
      end

      program main
      implicit none
      common/var/y,ypr 
      integer m1,m2,n,i
      real*8 t,y,ypr,h,g,gg,EL,ELM
      y=1.0D0
      t=0.0D0
      open(1,FILE='out.txt')
      write(1,70)
      write(1,80)
      write(1,70)
2     print*, 'Select:'
      print*, '1) h=0.1    (rungekutta, trapeze, mod euler, euler)'
      print*, '2) h=0.05   (rungekutta(2ty), trapeze, mod euler, euler)'
      print*, '3) h=0.025  (rungekutta(2ty), euler)'
      print*, '4) h=0.2    (trapeze, mod euler)'
      print*, '5) h=0.5    (rungekutta)'
      print*, '6) h=0.25   (rungekutta)'
      print*, '7) h=0.125  (rungekutta)'
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
      if(m2.eq.4)then
       h=0.2
       goto 11      
      endif 
      if(m2.eq.5)then
       h=0.5
       goto 11      
      endif 
      if(m2.eq.6)then
       h=0.25
       goto 11      
      endif 
      if(m2.eq.7)then
       h=0.125
       goto 11      
      endif 
11    print*, 'Select method:'
      print*, '<1> rungekutta (2ty)'
      print*, '<2> Euler unmanifest'  
      print*, '<3> TRAPEZE' 
      print*, '<4> Euler' 
      print*, '<5> Euler mod' 
      print*, '<6> rungekutta '
      read(*,*) m1  
      n=1/h+1 
      if(m1.eq.1)then
       do i=0,n
        write(*,300) t,y,g(t),abs(y-g(t))
        write(1,100) t,y,g(t),abs(y-g(t))
        write(1,70)
        call RungeKutta(h,t)
        t=t+h
       enddo 
      endif
      if(m1.eq.2)then
       do i=0,2*n
        write(*,300) t,y,gg(t),abs(y-gg(t))
        write(1,100) t,y,gg(t),abs(y-gg(t))
        write(1,70)
        ypr=y
        t=t+h
        call newton(t,h,m1)
        call Elun(h,t)
       enddo 
      endif 
      if(m1.eq.3)then
       do i=0,2*n
        write(*,300) t,y,gg(t),abs(y-gg(t))
        write(1,100) t,y,gg(t),abs(y-gg(t))
        write(1,70)
        ypr=y
        t=t+h
        call newton(t,h,m1)
        call TRAPEZE(h,t)
       enddo 
      endif
      if(m1.eq.4)then
       do i=0,2*n
        write(*,300) t,y,gg(t),abs(y-gg(t))
        write(1,100) t,y,gg(t),abs(y-gg(t))
        write(1,70)
        y = EL(h,t,y)
        t=t+h
       enddo 
      endif 
      if(m1.eq.5)then
       do i=0,2*n
        write(*,300) t,y,gg(t),abs(y-gg(t))
        write(1,100) t,y,gg(t),abs(y-gg(t))
        write(1,70)
        y = ELM(h,t,y)
        t=t+h
       enddo 
      endif
      if(m1.eq.6)then
       do i=0,2*n-2
        write(*,300) t,y,gg(t),abs(y-gg(t))
        write(1,100) t,y,gg(t),abs(y-gg(t))
        write(1,70)
        call RungeKutta2(h,t)
        t=t+h
       enddo 
      endif
      print*, 'Completed.'
      close(1)
      pause  
        
   70 format(52('-'))   
   80 format(' t          | y          | y*         | |y-y*|     |')  
  300 format('t = 'E11.4 ' y = 'E16.8 ' y* = 'E16.8 ' |y-y*| = 'E16.8)
  100 format(E11.4' |'E11.4' |'E11.4' |'E11.4' |')
      end
