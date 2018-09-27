        program fghm 
        implicit none
        integer :: nx,i,j,k,np,Ts,NE
        INTEGER :: LWORK,INFO
        real*8 ::dp,dx,xmin,xmax,pi,mass,D,beta,Xe,Kb,Z,P,F
        real*8 ::Z1,Z2,Z3,Z4
        real*8, allocatable :: x(:),E(:),F1(:),F2(:),F3(:),A(:)
        REAL*8, ALLOCATABLE :: WORK(:)
        real*8, allocatable :: H(:,:),S(:,:),C(:,:),T(:,:)
        real*8, allocatable :: V(:)
        pi=2.0*asin(1.0)
        mass=12861.93
        xmin=0.0
        xmax=10.0
        D=0.3640
        beta=1.4285
        Xe=2.075
        Kb=3.16685d-6
        NE=67
        Z=0.0 
         write(6,*)'inserte el numero de puntos (nx), temperatura(K)'
        read(5,*)nx,Ts
        
        np=(nx-1)/2
        dx=(xmax-xmin)/(nx-1)
        dp=(2.0*pi)/((nx-1)*dx)
        
        

        allocate (x(nx))
        do i=1,nx
        x(i)=xmin+(i-1)*dx
        enddo
        allocate(H(nx,nx),S(nx,nx),T(nx,nx))
        allocate (V(nx))
        do i=1,nx
        do j=i,nx
          T(i,j)=0.0
        do k=1,np
          T(i,j)=T(i,j)+(((k*dp)**2)*cos(k*dp*(x(j)-x(i))))
        
        enddo
        T(i,j)=T(i,j)/((nx-1)*mass)
        S(i,j)=0.0


        H(i,j)=T(i,j)
        if(i.eq.j)then
          
        V(i)=D*(1.0-dexp(-beta*(x(i)-Xe)))**2
        H(i,j)=H(i,j)+V(i)
        S(i,j)=1.0
               
        endif
        H(j,i)=H(i,j)
        enddo
        enddo
        !como ahorrar tiempo de calculo?? Se hace aprovechar la
        !hermiticidad del hamiltoniano
        !write(6,*)H(1,2),H(2,1),S(1,1),S(1,2),T(1,2),T(2,1),T(3,3)
        
        allocate(C(nx,nx),E(nx))
        LWORK=3*nx
        ALLOCATE(WORK(LWORK))
        CALL DSYGV(1,'V','U',nx,H,nx,S,nx,E,WORK,LWORK,INFO)
        DEALLOCATE(WORK)
        C=H/(dsqrt(dx))
        do i=1,nx
        !if(E(i).lt.D) write(6,*)i,E(i)
        enddo

      !  stop 'ok'
!!!!!!!!!!!!!!!! sirve para parar el programa
        do i=1,67
        Z=Z+(dexp(-E(i)/(Kb*Ts)))
        enddo
        write(6,*) 'Funcion de particion Z=',Z
       
        F=(dexp(-E(i)/(Kb*Ts)))/Z
        

        F=(dexp(-E(i)/(kb*300)))/Z
       !vamos a hacer una funcion de particion para cada temperatura
        
        Z1=0.0
        Z2=0.0
        Z3=0.0
        do i=1,67
        
        Z1=Z1+(dexp(-E(i)/(Kb*300)))
        Z2=Z2+(dexp(-E(i)/(Kb*1000)))
        Z3=Z3+(dexp(-E(i)/(Kb*4000)))
        enddo
        
        
        allocate(F1(i),F2(i),F3(i))
        do i=1,67
        F1(i)=(dexp(-E(i)/(kb*300)))/Z1
        
        F2(i)=(dexp(-E(i)/(kb*1000)))/Z2
        
        F3(i)=(dexp(-E(i)/(kb*4000)))/Z3
        enddo
              
       
        open(1,file='n2.dat')   
        do i=1,67
         write(1,*)i,F1(i),F2(i),F3(i)
         enddo
        close(1)
        
        allocate(A(nx))
        A(nx)=0.d0
        do j=1,nx
        do i=1,nx
        A(j)=A(j)+((C(i,j)*dexp(-(x(i)-4.0)**2)))/(dsqrt((2*pi)*(0.1)))
        enddo
        enddo    
        A=A*dx
        write(6,*)'El area es=', sum(A(:)*A(:)), 'valor de dx=', dx  
        write(6,*)'Valor de raiz de dx=', dsqrt(dx)
       
        !A(j)=A(j)+(C(i,j)*((dexp(-(x(i)-4.0)**2))/dsqrt((2*pi)*(0.1))))
         deallocate(H,T,E,S,C)
      !  write(6,*)'info', INFO
        

        end program
