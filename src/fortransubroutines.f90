!***********************************************************************************
!************************* Modelagem ***********************************************
!***********************************************************************************

SUBROUTINE modelagem(r,nt,Fx,Fz)

	IMPLICIT NONE
    
   CHARACTER                    :: r
   INTEGER                      :: k, Fx, Fz, zi           ! counter, Source Position
   INTEGER                      :: i,j					      !
   INTEGER,INTENT(in)           :: Nx,Nz                   ! Grid Elements
   REAL                         :: t, aux, fc              ! delay time and auxiliar
   REAL                         :: fonte, amort            ! initial time and wavelet
   REAL,INTENT(in)              :: nt,dh,dt
   REAL, DIMENSION(Nz,Nx)       :: aux_vel
   REAL, DIMENSION(Nz,Nx)       :: vel                     ! model
   REAL, DIMENSION(Nz,Nx)       :: P                       !Pressure Matrix
   REAL, DIMENSION(Nz,Nx)       :: Pf

   zi = 3 ! Profundidade dos receptores
   
   open(10, file = 'wavelet_ricker.dat')
   
   fonte = read(10,*)
   
   
              
do k = 1,nt	

	if k .le. size(fonte) - 1
               
		P(Fz,Fx) = P(Fz,Fx) + fonte(int(k))
		
	end if

!Chamando a subrotina operador_quarta_ordem

	call operador_quarta_ordem(Nz,Nx,dh,dt,vel,P,Pf)
    
!Chamando a Subrotina bordas_cerjan    

  ! call bordas_cerjan( ) 
	
   !Troca de Matrizes
   P = P 
   P = Pf
   
   if r == 'y'
   
       call snapshots(Nsnap,nt)

! Sismogramas
    
       
   do j = 0,Nx
                              
        Sismo(k,j) = Pf(zi,j)
        
   end do
   
! Salvando Sismogramas
   
   open(2, file='Sismograma.dat',&
         status='unknown',form='formatted')
   
   write(2,*) Sismo
   
end do   
        
   
	
END SUBROUTINE modelagem	




!***********************************************************************************
!************************* 4nd ORDER OPERATOR IN SPACE****************************
!***********************************************************************************  
SUBROUTINE operador_quarta_ordem(Nz,Nx,dh,dt,vel,P,Pf)
  IMPLICIT NONE
  INTEGER                                       :: i,j
  INTEGER,INTENT(in)                            :: Nx,Nz          !Grid Elements
  REAL,DIMENSION(Nz,Nx)                         :: aux_vel
  REAL,INTENT(in)                               :: dh,dt
  REAL, DIMENSION(Nz,Nx)                        :: vel            ! model
  REAL, DIMENSION(Nz,Nx)                        :: P              !Pressure Matrix
  REAL, DIMENSION(Nz,Nx)                        :: Pf

  aux_vel = (vel*vel)*(dt*dt)/(12*(dh*dh))  ! Remove it. Wasting cpu time

  do i=3,Nx-2
     do j=3,Nz-2
        !4th order in space and 2nd order in time
        Pf(j,i)=2*P(j,i)-Pf(j,i) + aux_vel(j,i)*&
             &(-(P(j,i-2) + P(j-2,i) + P(j+2,i) + P(j,i+2)) + & 
             &16*(P(j,i-1) + P(j-1,i) + P(j+1,i) + P(j,i+1))- &
             &60*P(j,i))
     end do
  end do

  !Bounary Condition - 2nd order in space
  !Left
  i=2
  do j=2,Nz-1
     Pf(j,i)=2*P(j,i)-Pf(j,i)+aux_vel(j,i)*(P(j,i-1)-2*P(j,i)+P(j,i+1) + P(j-1,i)-2*P(j,i)+P(j+1,i))*12
  end do
  !Right
  i=Nx-1
  do j=2,Nz-1
     Pf(j,i)=2*P(j,i)-Pf(j,i)+aux_vel(j,i)*(P(j,i-1)-2*P(j,i)+P(j,i+1) + P(j-1,i)-2*P(j,i)+P(j+1,i))*12
  end do
  !Top
  j=2
  do i=3,Nx-2
     Pf(j,i)=2*P(j,i)-Pf(j,i)+aux_vel(j,i)*(P(j,i-1)-2*P(j,i)+P(j,i+1) + P(j-1,i)-2*P(j,i)+P(j+1,i))*12
  end do
  !Bottom
  j=Nz-1
  do i=3,Nx-2
     Pf(j,i)=2*P(j,i)-Pf(j,i)+aux_vel(j,i)*(P(j,i-1)-2*P(j,i)+P(j,i+1) + P(j-1,i)-2*P(j,i)+P(j+1,i))*12
  end do
  RETURN 
END SUBROUTINE operador_quarta_ordem



!***********************************************************************************
!************************* WAVELET CALCULATION *************************************
!***********************************************************************************

SUBROUTINE wavelet(switch,dtime,MaxAmp,freqcut)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Wavele Ricker Calculation.          
    ! If you're using NSG operator the Ricker is defined by 2nd    
    ! derivative of gaussian function. And if you're using
    ! SSG operator the Ricker is defined by 1st derivative of
    ! gaussian function.
    ! INPUT:  
    ! dtime         = Time increment
    ! MaxAmp        = Max Amplitude of Source. 
    ! switch        = Select beewten NSG(1) and SSG(2) and ESG(3) operator 
    ! freqcut       =  Cut Frequency of wavelet
    ! 
    ! OUTPUT: ../analysis_files/wavelet_ricker.dat
    ! 
    ! 
    ! Code Written by Felipe Timoteo
    !                 Last update: May 15th, 2017
    !
    ! Copyright (C) 2017 Grupo de Imageamento Sísmico e Inversão Sísmica (GISIS)
    !                    Departamento de Geologia e Geofísica
    !                    Universidade Federal Fluminense
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IMPLICIT NONE

    INTEGER                                    :: k, Ntsource           ! counter, Source elements
    REAL                                       :: t, aux, fc            ! delay time and auxiliar
    REAL                                       :: time_0,vector_source  ! initial time and wavelet

    INTEGER,INTENT(in)                         :: switch
    REAL, INTENT(in)                           :: dtime,freqcut,MaxAmp
    REAL, parameter :: pi = 4.0 * atan(1.0)

    fc       = freqcut/(3.*sqrt(pi))        ! Ajust to cut of gaussian function
    time_0   = 2*sqrt(pi)/freqcut           ! Initial time source
    Ntsource = nint(2*time_0/dtime) + 1        ! Number of elements of the source

    open(77, file='wavelet_ricker.dat',&
         status='unknown',form='formatted')

    select case(switch)    !NSG or SSG modelling
    case(1)

       do k=1,Ntsource                          !Nts=nint(tm/dt)+1
          t=(k-1)*dtime-time_0                    !Delay Time
          aux=pi*(pi*fc*t)*(pi*fc*t)
          vector_source = MaxAmp*(2*aux-1)*exp(-aux)
          write(77,*) t,vector_source
       end do

    case(2)
       do k=1,Ntsource                          !Nts=nint(tm/dt)+1
          t=(k-1)*dtime-time_0                    !Delay Time
          aux=-pi*(pi*fc*t)*(pi*fc*t)
          vector_source = -MaxAmp*t*exp(aux)       
          write(77,*) t,vector_source
       end do

    case(3)

       do k=1,Ntsource                          !Nts=nint(tm/dt)+1
          t=(k-1)*dtime-time_0                    !Delay Time
          aux=pi*(pi*fc*t)*(pi*fc*t)
          vector_source = MaxAmp*(2*aux-1)*exp(-aux)
          write(77,*) t,vector_source
       end do

    end select
    close(77)
END SUBROUTINE wavelet


!***********************************************************************************
!************************* Bordas ***********************************************
!***********************************************************************************


!***********************************************************************************
!************************* Snapshots ***********************************************
!***********************************************************************************

SUBROUTINE snapshots(Nsnap,nt)
   
   
   
    if k % int(nt/Nsnap) == 0
      
        character(len=90) :: filename
                 
              fileloop: do j = 1, Nsnap
              
                write(filename, 'snapshot',I1,'.dat') j 
                      
                open(unit = j, status = 'old', file = filename)
                
                close j
                
              end do fileloop
                
    end if
    
END SUBROUTINE snapshots