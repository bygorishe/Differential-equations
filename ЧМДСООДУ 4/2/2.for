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
      
      real*8 function Fx(t,h)
      implicit none
      common/var/y,ypr1,ypr2,ypr3
      real*8 t,y,h,f,ypr1,ypr2,ypr3
      Fx=y-ypr1-h/12*(5*f(t,y)+8*f(t-h,ypr1)-1*f(t-h*2,ypr2))
      end
      
      real*8 function Gx(h,t)
      implicit none
      real*8 h,t
      Gx=1-5*t*h/6     !Fx/dx
      end
      
      real*8 function Fxx(t,h)
      implicit none
      common/var/y,ypr1,ypr2,ypr3
      real*8 t,y,h,f,ypr1,ypr2,ypr3
      Fxx=y-ypr1-h/24*(9*f(t,y)+19*f(t-h,ypr1)-5*f(t-h*2,ypr2)
     *+1*f(t-h*3,ypr3))
      end
      
      real*8 function Gxx(h,t)
      implicit none
      real*8 h,t
      Gxx=1-9*t*h/12    !Fx/dx
      end
      
      subroutine newton(t,h,m1)
      implicit none
      common/var/y,ypr1,ypr2,ypr3
      real*8 t,y,yy,h,Fx,Gx,Fxx,Gxx,ypr1,ypr2,ypr3
      integer m1
      yy=0.0D0
      if(m1.eq.3)then
       do while(abs(y-yy).gt.0.00000001)  
        yy = y 
        y=y-Fx(t,h)/Gx(h,t)
        !write(*,*) y
       end do
      endif
      if(m1.eq.4)then
       do while(abs(y-yy).gt.0.00000001)  
        yy = y 
        y=y-Fxx(t,h)/Gxx(h,t)
        !write(*,*) y
       end do
      endif
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
     
      real*8 function AD3(h,t)  
      implicit none
      common/var/y,ypr1,ypr2,ypr3
      real*8 h,t,y,f,ypr1,ypr2,ypr3
       AD3=y+h/12*(23*f(t,y)-16*f(t-h,ypr1)+5*f(t-h*2,ypr2))
      end
     
      real*8 function AD4(h,t)  
      implicit none
      common/var/y,ypr1,ypr2,ypr3
      real*8 h,t,y,f,ypr1,ypr2,ypr3
       AD4=y+h/24*(55*f(t,y)-59*f(t-h,ypr1)+37*f(t-h*2,ypr2)
     *-9*f(t-h*3,ypr3))
      end
      
      
      real*8 function ADU3(h,t)  
      implicit none
      common/var/y,ypr1,ypr2,ypr3
      real*8 h,t,y,f,ypr1,ypr2,ypr3
       ADU3=ypr1+h/12*(5*f(t,y)+8*f(t-h,ypr1)-1*f(t-h*2,ypr2))
      end
     
      real*8 function ADU4(h,t)  
      implicit none
      common/var/y,ypr1,ypr2,ypr3
      real*8 h,t,y,f,ypr1,ypr2,ypr3
       ADU4=ypr1+h/24*(9*f(t,y)+19*f(t-h,ypr1)-5*f(t-h*2,ypr2)
     *+1*f(t-h*3,ypr3))
      end
    

     

      program main
      implicit none
      common/var/y,ypr1,ypr2,ypr3
      integer m1,m2,n,i
      real*8 t,y,ypr1,q,ypr2,ypr3,h,g,AD3,AD4,ADU3,ADU4
      y=1.0D0
      ypr1=0.0D0
      ypr2=0.0D0
      ypr3=0.0D0
      t=0.0D0
      open(1,FILE='out.txt')
      write(1,70)
      write(1,80)
      write(1,70)
2     print*, 'Select:'
      print*, '1) h=0.1  '
      print*, '2) h=0.05 '
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
      print*, '<1> see'
      print*, '<2> un'
      print*, '<3> 3' 
      print*, '<4> 4' 
      read(*,*) m1 
      
      n=1/h+1 
      if(m1.eq.1)then                             !nahodim pervie k tochki 
       do i=0,1
        write(*,300) t,y,g(t),abs(y-g(t))
        write(1,100) t,y,g(t),abs(y-g(t))
        write(1,70)
        ypr2=ypr1
        ypr1=y
        call RungeKutta(h,t)
        t=t+h  
       enddo
       do i=2,n
        q=y
        write(*,300) t,y,g(t),abs(y-g(t))
        write(1,100) t,y,g(t),abs(y-g(t))
        write(1,70)
        y = AD3(h,t)
        t=t+h
        ypr2=ypr1
        ypr1=q
       enddo 
      endif
      if(m1.eq.2)then
       do i=0,2
        write(*,300) t,y,g(t),abs(y-g(t))
        write(1,100) t,y,g(t),abs(y-g(t))
        write(1,70)
        ypr3=ypr2
        ypr2=ypr1
        ypr1=y
        call RungeKutta(h,t)
        t=t+h  
       enddo
       do i=3,n
        q=y
        write(*,300) t,y,g(t),abs(y-g(t))
        write(1,100) t,y,g(t),abs(y-g(t))
        write(1,70)
        y = AD4(h,t)
        t=t+h
        ypr3=ypr2
        ypr2=ypr1
        ypr1=q
       enddo 
      endif 
      if(m1.eq.3)then

        write(*,300) t,y,g(t),abs(y-g(t))
        write(1,100) t,y,g(t),abs(y-g(t))
        write(1,70)
        ypr2=ypr1
        ypr1=y
        call RungeKutta(h,t)
        t=t+h  

       do i=1,n
        write(*,300) t,y,g(t),abs(y-g(t))
        write(1,100) t,y,g(t),abs(y-g(t))
        write(1,70)
        t=t+h
        ypr2=ypr1
        ypr1=y
        call newton(t,h,m1)
        y = ADU3(h,t)
       enddo 
      endif
      if(m1.eq.4)then
       do i=0,1
        write(*,300) t,y,g(t),abs(y-g(t))
        write(1,100) t,y,g(t),abs(y-g(t))
        write(1,70)
        ypr3=ypr2
        ypr2=ypr1
        ypr1=y
        call RungeKutta(h,t)
        t=t+h  
       enddo
       do i=2,n
        write(*,300) t,y,g(t),abs(y-g(t))
        write(1,100) t,y,g(t),abs(y-g(t))
        write(1,70)
        t=t+h
        ypr3=ypr2
        ypr2=ypr1
        ypr1=y
        call newton(t,h,m1)
        y = ADU4(h,t)
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
