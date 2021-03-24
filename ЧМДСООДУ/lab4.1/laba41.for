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

      real function EL(h,t,y)   !Method Eulera
      implicit none
      real h,t,y
      real f
       EL=y+h*f(t,y)
      end
      
      real function ELM(h,t,y)  !Method Eulera mod
      implicit none
      real h,t,y
      real f
       ELM=y+h/2*(f(t,y)+f((t+h),(y+h*f(t,y))))
      end
      
      real function ELMM(h,t,y)   !Method Eulera usoversh
      implicit none
      real h,t,y
      real f
       ELMM=y+h*f((t+h/2),(y+h/2*f(t,y)))
      end

     
 
      program main
      implicit none
      integer m1,m2,n,i
      real EL,ELM,ELMM,t,y,h,g
      y=1.
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
      print*, '<3> elmodmod'    
      read(*,*) m1  
      n=1/h+1 
      if(m1.eq.1)then
       do i=0,n
        !write(*,300) t,y,g(t),abs(y-g(t))
        write(1,100) t,y,g(t),abs(y-g(t))
        write(1,70)
        y = EL(h,t,y)
        t=t+h
       enddo 
      endif
      if(m1.eq.2)then
       do i=0,n
        !write(*,300) t,y,g(t),abs(y-g(t))
        write(1,100) t,y,g(t),abs(y-g(t))
        write(1,70)
        y = ELM(h,t,y)
        t=t+h
       enddo
      endif 
      if(m1.eq.3)then
       do i=0,n
        !write(*,300) t,y,g(t),abs(y-g(t))
        write(1,100) t,y,g(t),abs(y-g(t))
        write(1,70)
        y = ELMM(h,t,y)
        t=t+h
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


