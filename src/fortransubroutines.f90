!***********************************************************************************
!************************* Modelagem ***********************************************
!***********************************************************************************

SUBROUTINE nucleomodelagem(Nz,Nx,Nt,dh,dt,NpCA,shot,shotshow,NSx,NSz,fonte,Nfonte,Nsnap,regTTM,caminho_modelo)


  ! SOCORRO: Valores de Nsnap e Nfonte estao trocados mas funcionando mesmo assim :o 
  ! Esse problema esta na linha 151 do codigo em python

  ! PROBLEMA POSSIVELMENTE RESOLVIDO: O Nfonte e tratado pela f2py como um argumento opcional portanto nao precisamos chama-lo na funcao do python, somente no fortran. Isso porque o Nfonte ja indentificado como tamanho do vetor

  !Veja linha 29

  !Na duvida leia: https://stackoverflow.com/questions/35528927/f2py-order-of-function-arguments-messed-up


  IMPLICIT NONE  

  INTEGER                        :: k,aux
  INTEGER                        :: Nzz,Nxx                     !Expanded dimensions
  CHARACTER(len=256)             :: caminho_modelo

  INTEGER,INTENT(in)             :: Nsnap!,Nshot
  INTEGER                        :: count_snap
  INTEGER,INTENT(in)             :: shot,shotshow,NSx,NSz,Nfonte     ! Related source
  INTEGER,INTENT(in)             :: Nx,Nz,Nt,NpCA                    ! Grid Elements

  INTEGER,INTENT(in)             :: regTTM                         ! Condition Transit Time Matrix

  REAL,INTENT(in)                :: dh,dt                            
  REAL,DIMENSION(Nfonte)         :: fonte                            ! Source  
  REAL,DIMENSION(NpCA)           :: func_Am  
  REAL,DIMENSION(Nt,Nx)          :: Seism                             
  REAL,DIMENSION(Nz,Nx)          :: TTM, ATTM                        !Related Transit Time Matrix
  REAL,ALLOCATABLE,DIMENSION(:,:):: P,Pf,vel                          

  Nxx = NpCA + Nx + NpCA
  Nzz = Nz + NpCA
  
  ALLOCATE(P(Nzz,Nxx))
  ALLOCATE(Pf(Nzz,Nxx))
  ALLOCATE(vel(Nzz,Nxx))

  ! Load Damping Function
  open(20,file='f_amort.dat',&
       status='unknown',form='formatted')
  do k=1,NpCA
     read(20,*)func_Am(k)
  end do
  close(20)

  aux = Nsnap
  aux = Nt/aux              ! evaluate number of snapshots

  count_snap = 0
  TTM = 0.0
  ATTM = 0.0 

 ! revisar nome de entrada do modelo 

  CALL   LoadVelocityModelExpanded(Nz,Nzz,Nx,Nxx,NpCA,trim(caminho_modelo),vel)
 
  
  P    = 0.0                   !Pressure field
  Pf   = 0.0                   !Pressure field in future  

  do k=1,Nt

     if (k <= Nfonte) then        !Nts=nint(tm/dt)+1!number of source time elements
        P(Nsz,Nsx + NpCA)= P(Nsz,Nsx+NpCA) - fonte(k) ! NpCA para o modelo expandido
     end if

     CALL operador_quarta_ordem(Nzz,Nxx,dh,dt,vel,P,Pf)

     CALL Updatefield(Nzz,Nxx,P,Pf)

     CALL CerjanNSG(Nzz,Nxx,NpCA,func_Am,P,Pf)            

     CALL ReynoldsEngquistNSG(Nzz,Nxx,dh,dt,vel,P,Pf)

     ! Revisar posicionamento dos receptores

     if (regTTM == 0) then
       Seism(k,:) = P(10,NpCA+1:NpCA+Nx)
     end if 
     
     if ( (mod(k,aux)==0)  .and. shotshow >0 .and. shotshow == shot) then 
        ! print *, "k=",k , "time=", (k-1)*dt

        CALL snap(Nzz,Nxx,count_snap,shotshow,"Marmousi","../snapshot/",P(1:Nz,NpCA+1:NpCA+Nx))
     end if
     
     if (regTTM == 1) then 
        CALL TransitTimeMatrix(Nz,Nx,k,P(1:Nz,NpCA+1:NpCA+Nx),TTM,ATTM,shot)
     end if

  end do

  if (regTTM == 0) then

     if (caminho_modelo  == "marmousi_vp_383x141.bin") then
        CALL Seismogram(Nt,Nx,shot,"Marmousi","../sismograma/",Seism)
     end if


     if (caminho_modelo == '../modelos_utilizados/velocitymodel_Hmgns_wtrly.bin') then   
        CALL Seismogram(Nt,Nx,shot,"Homogeneo","../sismograma_modelo_camada_de_agua/",Seism)
     end if
end if

     if (regTTM == 1) then
        CALL writematrix(Nz,Nx,shot,TTM, "Marmousi","../matriz_tempo_transito/")
     end if
  
END SUBROUTINE nucleomodelagem


!***********************************************************************************
!************************* Migracao ************************************************
!***********************************************************************************  



SUBROUTINE migracao(Nz,Nx,Nt,dh,dt,NpCA,zr,shot,shotshow,Nsnap,caminho_modelo)

  IMPLICIT NONE  

  INTEGER                        :: k,aux
  INTEGER                        :: Nzz,Nxx                     !Expanded dimensions
  CHARACTER(len=256)             :: caminho_modelo

  INTEGER,INTENT(in)             :: Nsnap
  INTEGER                        :: count_snap
  INTEGER,INTENT(in)             :: shot,shotshow               ! Related source
  INTEGER,INTENT(in)             :: Nx,Nz,Nt,NpCA                    ! Grid Elements

  REAL,DIMENSION(NpCA)           :: func_Am 
  INTEGER,INTENT(in)             :: zr                         ! Condition Transit Time Matrix


  REAL,INTENT(in)                 :: dh,dt                                                         
  REAL,DIMENSION(Nt,Nx)           :: Seism                             
  REAL,DIMENSION(Nz,Nx)           :: TTM                        !Related Transit Time Matrix
  REAL,ALLOCATABLE,DIMENSION(:,:) :: P,Pf,vel,Imagem                        

  Nxx = NpCA + Nx + NpCA
  Nzz = Nz + NpCA
  
  ALLOCATE(P(Nzz,Nxx))
  ALLOCATE(Pf(Nzz,Nxx))
  ALLOCATE(vel(Nzz,Nxx))
  ALLOCATE(Imagem (Nz,Nx))

  ! Load Damping Function
  open(20,file='f_amort.dat',&
       status='unknown',form='formatted')
  do k=1,NpCA
     read(20,*)func_Am(k)
  end do

!  caminho_modelo = '../modelo_suavizado/Suave_v15_marmousi_vp_383x141.bin'

! Abrindo o Modelo Suavizado, Sismograma sem a onda direta e a Matriz de Tempo de Transito

  
  CALL  LoadVelocityModelExpanded(Nz,Nzz,Nx,Nxx,NpCA,trim(caminho_modelo),vel)
  CALL  LoadSeismogram(Nt,Nx,shot,"Marmousi","../sismograma_sem_onda_direta/",Seism)
  CALL  LoadVelocityModel(Nz,Nx,'../matriz_tempo_transito/Marmousi_shot001.bin',TTM)

  P    = 0.0                   !Pressure field
  Pf   = 0.0                   !Pressure field in future  
  Imagem = 0.0
  count_snap = 0


  aux = Nsnap
  aux = Nt/aux              ! evaluate number of snapshots


  do k = Nt,1,-1

       P(zr,1:Nx) =  Seism(k,1:Nx) + P(zr,1:Nx)

       CALL operador_quarta_ordem(Nzz,Nxx,dh,dt,vel,P,Pf)

       CALL Updatefield(Nzz,Nxx,P,Pf)

       CALL CerjanNSG(Nzz,Nxx,NpCA,func_Am,P,Pf)            

       CALL ReynoldsEngquistNSG(Nzz,Nxx,dh,dt,vel,P,Pf)
       
       CALL ImagingConditionMaxAmP(k,Nz,Nx,P(1:Nz,NpCA+1:NpCA+Nx),TTM,Imagem)
       
       
     if ( (mod(k,aux)==0)  .and. shotshow > 0 .and. shotshow == shot) then 

        CALL snap(Nzz,Nxx,count_snap,shotshow,"Marmousi","../snapshot_migracao_rtm/",P(1:Nz,NpCA+1:NpCA+Nx))

     end if

  end do
  
       CALL writematrix(Nz,Nx,shot,Imagem,"Imagem_Marmousi","../Imagem/")

END SUBROUTINE migracao

!***********************************************************************************
!************************* 4nd ORDER OPERATOR IN SPACE******************************
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
!*************************LOAD VELOCITY MODEL***************************************
!***********************************************************************************

SUBROUTINE LoadVelocityModel(Nz,Nx,modelpath,vel)
  IMPLICIT NONE
  INTEGER                                         :: i,j                     !Counters
  LOGICAL                                         :: filevel                 !Check if file exists

  CHARACTER(*),INTENT(in)                         :: modelpath
  INTEGER,INTENT(in)                              :: Nz,Nx                   !Grid parameter
  REAL, DIMENSION(Nz,Nx),INTENT(out)              :: vel                     !Rock physics

  INQUIRE(file=modelpath,exist=filevel) !verify if parameters file exist

  if (filevel) then

     OPEN(23, FILE=modelpath, STATUS='unknown',&
          &FORM='unformatted',ACCESS='direct', RECL=(Nx*Nz*4))

     read(23,rec=1) ((vel(j,i),j=1,nz),i=1,nx)       ! read velocity matrix

     close(23)
  else
     print*, ''
     print*,'============================================================================='
     print*, 'Velocity model input File doesnt exist!'
     print*, 'Please create a velocity model and '
     print*, 'put it in the Models folder.'
     print*,'============================================================================='
     print*, ''
     print*, 'PRESS RETURN TO EXIT...   '
     read(*,*)
     stop

  end if

  RETURN
END SUBROUTINE LoadVelocityModel

!***********************************************************************************
!************************ LOAD EXPANDED VELOCITY MODEL   ***************************
!***********************************************************************************

SUBROUTINE LoadVelocityModelExpanded(Nz,Nzz,Nx,Nxx,NpCA,modelpath,vel)
  ! Load a Velocity model and expand it to consider the absorbing boundaries
  ! outside of valid model.
  ! 
  ! INPUT: 
  ! Nz            = Number of grid points in z direction
  ! Nzz           = Nz + NpCA
  ! Nx            = Number of grid points in z direction
  ! Nxx           = Nx + NpCA
  ! NpCA          = Number of point in absorving boundary
  ! modelpath     = Model path
  ! 
  ! OUTPUT: 
  ! vel           = Expanded velocity model matrix
  !
  ! Code Written by Felipe Timoteo
  !                 Last update: 31th Jan, 2018
  !
  ! Copyright (C) 2018 Grupo de Imageamento Sísmico e Inversão Sísmica (GISIS)
  !                    Departamento de Geologia e Geofísica
  !                    Universidade Federal Fluminense
  
  IMPLICIT NONE
  INTEGER                                         :: i,j           !Counters    

  LOGICAL                                         :: filevel       !Check if file exists
  CHARACTER(*),INTENT(in)                         :: modelpath
  INTEGER,INTENT(in)                              :: NpCA
  INTEGER,INTENT(in)                              :: Nz,Nx         !Grid parameter
  INTEGER,INTENT(in)                              :: Nzz,Nxx       !Expanded Grid parameter 
  REAL, DIMENSION(Nzz,Nxx),INTENT(out)            :: vel           !Rock physics


  INQUIRE(file=modelpath,exist=filevel) !verify if parameters file exist

  if (filevel) then

     OPEN(23, FILE=modelpath, STATUS='unknown',&
          &FORM='unformatted',ACCESS='direct', RECL=(Nx*Nz*4))

     read(23,rec=1) ((vel(j,i),j=1,Nz),i=NpCA+1,NpCA+Nx)       ! read velocity matrix

     ! Left Boundary
     do i = 1,NpCA
        do j = 1,Nz
           vel(j,i) = vel(j,NpCA+1)
        end do
     end do

     ! Right Boundary
     do i = NpCA+Nx+1,Nxx
        do j = 1,Nz
           vel(j,i) = vel(j,NpCA+Nx)
        end do
     end do

     !Bottom Boundary
     do i = 1,Nxx
        do j = Nz+1,Nzz
           vel(j,i) = vel(Nz,i)
        end do
     end do

     close(23)

  else
     print*, ''
     print*,'============================================================================='
     print*, 'Velocity model input File doesnt exist!'
     print*, 'Please create a velocity model and '
     print*, 'put it in the Models folder.'
     print*,'============================================================================='
     print*, ''
     print*, 'PRESS RETURN TO EXIT...   '
     read(*,*)
     stop

  end if
  RETURN
END SUBROUTINE LoadVelocityModelExpanded

!***********************************************************************************
!************************* Updatefield *********************************************
!***********************************************************************************

SUBROUTINE Updatefield(Nz,Nx,P,Pf)
  IMPLICIT NONE

  INTEGER                              :: i,j
  REAL                                 :: Update                 !Updating

  INTEGER,INTENT(in)                   :: Nx,Nz
  REAL,DIMENSION(Nz,Nx),INTENT(inout)  :: P
  REAL,DIMENSION(Nz,Nx),INTENT(inout)  :: Pf
  do i=1,Nx
     do j=1,Nz

        update = P(j,i) !Updating
        P(j,i)=Pf(j,i)  !Updating
        Pf(j,i)=update  !Updating

     end do
  end do
END SUBROUTINE Updatefield


!***********************************************************************************
!************************* WRITTING SEISMOGRAM *************************************
!***********************************************************************************

SUBROUTINE Seismogram(Ntime,Nxspace,Nshot,outfilename,select_folder,Seismmatrix)
  ! Writting Seismogram in a binary file
  ! 
  ! INPUT: 
  ! 
  ! Ntime         = Total Number of Samples in Time
  ! Nxspace       = Total Number of Grid Points in X direction
  ! Nshot         = Shot Number
  ! outfilename   = Prefix in Seismogram filename
  ! 
  ! select_folder = Folder of Seismogram file
  ! myID          = Number of identification of process (MPI Parameter)
  ! proc_name     = Name of processor (MPI Parameter)
  ! 
  ! OUTPUT: ../select_folder/outfile_SeismogramShot.bin
  ! SeismMatrix   = Matrix (Ntime,Nxspace) that will writting Seismogram

  ! Code Written by Felipe Timoteo
  !                 Last update: May 23th, 2016
  !
  ! Copyright (C) 2017 Grupo de Imageamento Sísmico e Inversão Sísmica (GISIS)
  !                    Departamento de Geologia e Geofísica
  !                    Universidade Federal Fluminense
  IMPLICIT NONE

  CHARACTER(len=3)                               :: num_shot           ! write differents files
  CHARACTER(LEN=*),INTENT(in)                    :: select_folder      ! folder
  INTEGER,INTENT(in)                             :: Nxspace,Ntime      ! Total Number of Samples in Space and Time
  INTEGER,INTENT(in)                             :: Nshot              ! Number of Shot and ID processor
  CHARACTER(LEN=*),INTENT(in)                    :: outfilename        ! output filename pattern
  REAL, DIMENSION(Ntime,Nxspace),INTENT(in)      :: Seismmatrix        ! Seismogram

  write(num_shot,"(i3.3)")Nshot ! write shot counter in string to write differentes Seismograms

  OPEN(11, FILE=trim(select_folder)//trim(outfilename)//'_sismograma'//num_shot//'.bin', STATUS='unknown',&
       &FORM='unformatted',ACCESS='direct', RECL=(Ntime*Nxspace*4))


  write(11,rec=1) Seismmatrix
  CLOSE(11)

  RETURN
END SUBROUTINE Seismogram


!************************************************************************************
!************************* LOADING SEISMOGRAM  *************************************
!***********************************************************************************
SUBROUTINE LoadSeismogram(Ntime,Nxspace,Nshot,infilename,select_folder,SeismMatrix)
  ! Load a Seismogram from a binary file
  ! 
  ! INPUT:  ../select_folder/outfile_SeismogramShot.bin
  ! 
  ! Ntime         = Total Number of Samples in Time
  ! Nxspace       = Total Number of Grid Points in X direction
  ! Nshot         = Shot Number
  ! infilename   = Prefix in Seismogram filename
  ! 
  ! select_folder = Folder of Seismogram file
  ! myID          = Number of identification of process (MPI Parameter)
  ! proc_name     = Name of processor (MPI Parameter)
  ! SeismMatrix   = Matrix (Nt,Nx) that will receive Seismogram
  ! 
  ! OUTPUT: None
  ! 
  ! Code Written by Felipe Timoteo
  !                 Last update: May 23th, 2016
  !
  ! Copyright (C) 2017 Grupo de Imageamento Sísmico e Inversão Sísmica (GISIS)
  !                    Departamento de Geologia e Geofísica
  !                    Universidade Federal Fluminense


  IMPLICIT NONE
  CHARACTER(len=3)                               :: num_shot           !write differents files
  INTEGER                                        :: kk,ii              !Counter   
  LOGICAL                                        :: fileSeis           !Check if file exists

  CHARACTER(LEN= *),INTENT(in)                   :: select_folder      !folder
  CHARACTER(LEN= *),INTENT(in)                   :: infilename         !output filename pattern
  INTEGER,INTENT(in)                             :: Nxspace,Ntime,Nshot

  REAL, DIMENSION(Ntime,Nxspace),INTENT(out)     :: SeismMatrix       !Seismogram

  ! print*,'...............................................'
  ! write(*,"(A11,A10,A1,i3,A22,i3)"), 'Processor:',proc_name,'-',myID,'Loading Seismogram',Nshot
  ! write(*,"(A14,A20,A3)"),'in the folder ',select_folder ,'...'
  ! print*,'...............................................'

  write(num_shot,"(i3.3)")Nshot ! write shot counter in string to write differentes Seismograms

  INQUIRE(file=trim(select_folder)//trim(infilename)//'_sismograma'//num_shot//'.bin',&
       exist=fileSeis) !verify if parameters file exist

  if (fileSeis) then

     OPEN(11, FILE=trim(select_folder)//trim(infilename)//'_sismograma'//num_shot//'.bin', STATUS='unknown',&
          &FORM='unformatted',ACCESS='direct', RECL=(Ntime*Nxspace*4))
     read(11,rec=1) ((SeismMatrix(kk,ii),kk=1,Ntime),ii=1,Nxspace)
     close(11)

  else


     print*, ''
     print*,'============================================================================='
     print*, 'Seismogram ',infilename,num_shot, ' NOT FOUND. Do you have sure that this Seismogram'
     print*, ' is in the folder:', select_folder, '? Please, if you not sure'
     print*, 'check the folder and try again.'
     print*,'============================================================================='
     print*, ''
     print*, 'PRESS RETURN TO EXIT...   '
     read(*,*)
     stop

  end if

  RETURN
END SUBROUTINE LoadSeismogram



!***********************************************************************************  
!**************************** SNAPSHOT *********************************************
!***********************************************************************************

SUBROUTINE snap(Nz,Nx,count_snap,shot,outfile,folder,Field)
  IMPLICIT NONE
  CHARACTER(len=3)                               :: num_shot,num_snap  !write differents files

  CHARACTER(LEN=*),INTENT(in)                    :: outfile,folder            !output filename pattern
  INTEGER,INTENT(inout)                          :: count_snap
  INTEGER,INTENT(in)                             :: Nx,Nz,shot
  REAL, DIMENSION(Nz,Nx),INTENT(in)              :: Field

  count_snap = count_snap + 1          !snap couting

  write(num_shot,"(i3.3)")shot         !write shot counter in string 
  write(num_snap,"(i3.3)")count_snap   !change in string


  OPEN(10, FILE=trim(folder)//trim(outfile)//'_shot'//num_shot//'snap'//num_snap //'.bin', STATUS='unknown',&
       &FORM='unformatted',ACCESS='direct', RECL=(Nz*Nx*4))

  write(10,rec=1) Field !write Pressure matrix    
  close(10)

  RETURN
END SUBROUTINE snap

!***********************************************************************
!*******************CERJAN - DAMPING LAYER - NSG ***********************
!***********************************************************************

SUBROUTINE CerjanNSG(Nz,Nx,NpCA,func_Am,P,Pf)
  IMPLICIT NONE
  INTEGER                           :: i,j,k
  REAL                              :: aux

  INTEGER,INTENT(in)                :: Nx,Nz,NpCA
  REAL, DIMENSION(npca),INTENT(in)  :: func_Am
  REAL,DIMENSION(Nz,Nx),INTENT(inout) :: P, Pf

  !Amortecimento a esquerda
  DO i=3, NpCA !i=2, NpCA
     aux = func_Am(i)
     DO j=2, Nz-1
        P(j,i)  = P(j,i)*aux
        Pf(j,i) = Pf(j,i)*aux
     ENDDO
  ENDDO

  !Amortecimento a direita
  k = NpCA
  DO i=Nx-NpCA+1, Nx-2 !i=Nx-NpCA+1, Nx-1
     aux = func_Am(k)
     DO j=2, Nz-1
        P(j,i)  = P(j,i)*aux
        Pf(j,i) = Pf(j,i)*aux
     ENDDO
     k= k -1
  ENDDO

  !    Amortecimento fundo
  DO i=2, Nx-1
     k = NpCA
     DO j=Nz-NpCA+1, Nz-2  !j=Nz-NpCA+1, Nz-1
        aux=func_Am(k)
        P(j,i)  = P(j,i)*aux
        Pf(j,i) = Pf(j,i)*aux
        k = k - 1
     ENDDO
  ENDDO

  ! !    Amortecimento cima
  ! DO i=2, Nx-1
  !    DO j=2, NpCA
  !       aux = func_Am(j)
  !       P(j,i)  = P(j,i)*aux
  !       Pf(j,i) = Pf(j,i)*aux          
  !    ENDDO
  ! ENDDO
  RETURN
END SUBROUTINE CerjanNSG

!***********************************************************************
!**********Reynolds/Engquist - DAMPING LAYER - NSG *********************
!***********************************************************************

SUBROUTINE ReynoldsEngquistNSG(Nz,Nx,dh,dt,vel,P,Pf)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ReynoldsEngquistNSG calculates Non Reflective Bounary Condition.          
  !
  ! INPUT:  
  ! dtime         = Time increment
  ! time_0        = Initial time source
  ! Ntsource      = Total Number of time elements of Source
  ! MaxAmp        = Max Amplitude of Source. 
  ! switch        = Select beewten NSG(1) and SSG(2) operatior 
  ! freqcut       =  Cut Frequency of wavelet

  ! OUTPUT: 
  ! vectorsource  = Source amplitude vector
  ! 
  ! Code Written by Felipe Timoteo
  !                 Last update: May 15th, 2017
  !
  ! Copyright (C) 2017 Grupo de Imageamento Sísmico e Inversão Sísmica (GISIS)
  !                    Departamento de Geologia e Geofísica
  !                    Universidade Federal Fluminense
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IMPLICIT NONE
  INTEGER                                       :: i,j
  INTEGER,INTENT(in)                            :: Nx,Nz          !Grid Elements
  REAL,DIMENSION(Nz,Nx)                         :: aux_vel
  REAL,INTENT(in)                               :: dh,dt
  REAL, DIMENSION(Nz,Nx),INTENT(in)             :: vel            ! model
  REAL, DIMENSION(Nz,Nx),INTENT(inout)          :: Pf,P           !Pressure Matrix

  aux_vel = vel*dt/dh
  !      Inferior
  DO i=2,Nx-1
     Pf(Nz,i) = P(Nz,i) - aux_vel(Nz,i)*(P(Nz,i)-P(Nz-1,i))
  ENDDO

  !      Direita
  DO j=1,Nz
     Pf(j,Nx) = P(j,Nx) - aux_vel(j,Nx)*(P(j,Nx)-P(j,Nx-1))
  ENDDO

  !      Esquerda
  DO j=1,Nz
     Pf(j,1)  = P(j,1)  + aux_vel(j,1)*(P(j,2)-P(j,1))
  ENDDO

  ! !      Superior
  ! DO i=2,Nx-1
  !    Pf(1,i) = P(1,i) + aux_vel(1,i)*(P(2,i)-P(1,i))
  ! ENDDO

END SUBROUTINE ReynoldsEngquistNSG

!***********************************************************************************
!************************* Trannsit Time Matrix ************************************
!***********************************************************************************

SUBROUTINE TransitTimeMatrix(Nz,Nx,k,P,TTM,ATTM)
  ! Calculation of Transit Time Matrix
  ! The Criteria choose to evaluate a Transit Time Matrix was 
  ! maximum amplitude of pressure field
  ! INPUT:  
  ! Nz            = Total Number of Grid Points in Z direction   
  ! Nx            = Total Number of Grid Points in X direction   
  ! P             = Current Pressure Field
  ! TTM           = update Transit Time Matrix
  ! ATTM          = previous Transit Time Matrix

  ! OUTPUT: 
  ! TTM           = updated Transit Time Matrix
  ! ATTM          = previous Transit Time Matrix
  ! 
  ! steplength    = Step Length
  ! 
  ! Code Written by Felipe Timoteo
  !                 Last update: May 1st, 2017
  !
  ! Copyright (C) 2017 Grupo de Imageamento Sísmico e Inversão Sísmica (GISIS)
  !                    Departamento de Geologia e Geofísica
  !                    Universidade Federal Fluminense
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IMPLICIT NONE

  INTEGER                                :: i,j
  INTEGER,INTENT(in)                     :: Nx,Nz,k

  REAL,DIMENSION(Nz,Nx),INTENT(inout)    :: P,TTM,ATTM

  do i = 1,Nx
     do j = 1,Nz
        if (abs(P(j,i)) > abs(ATTM(j,i))) then
           ATTM(j,i) = P(j,i)
           TTM(j,i)  = k 
        end if
     end do
  end do
  RETURN
END SUBROUTINE TransitTimeMatrix

!***********************************************************************************
!***************************** IMAGING CONDITION ***********************************
!***********************************************************************************

SUBROUTINE ImagingConditionMaxAmP(k,Nz,Nx,P,TTM,Image)

  IMPLICIT NONE

  INTEGER                              :: i,j

  INTEGER,INTENT(in)                   :: Nx,Nz,k
  REAL,DIMENSION(Nz,Nx),INTENT(in)     :: P,TTM
  REAL,DIMENSION(Nz,Nx),INTENT(inout)  :: Image


  do i = 1,Nx
     do j = 1,Nz

        if (k == TTM(j,i)) then
           Image(j,i) = P(j,i)
        end if

     end do
  end do

  RETURN
END SUBROUTINE ImagingConditionMaxAmP


!***********************************************************************************
!***************************** WRITING MATRIX **************************************
!***********************************************************************************

SUBROUTINE writematrix(Nz,Nx,shot,Matrix,outfile,folder)
  IMPLICIT NONE

  CHARACTER(len=3)                    :: num_shot
  CHARACTER(LEN=*),INTENT(in)         :: outfile,folder

  INTEGER,INTENT(in)                  :: Nz
  INTEGER,INTENT(in)                  :: Nx
  INTEGER,INTENT(in)                 :: shot
  REAL,DIMENSION(Nz,Nx), INTENT(in)   :: Matrix

  write(num_shot,"(i3.3)")shot         !write shot counter in string   

  OPEN(12, FILE=trim(folder)//trim(outfile)//'_shot'//num_shot //'.bin', STATUS='unknown',&
       &FORM='unformatted',ACCESS='direct', RECL=(Nz*Nx*4))


  write(12,rec=1) Matrix !write  matrix    
  close(12)

END SUBROUTINE writematrix


!*********************************************************************************
!************************ REMOVE ONDA DIRETA *************************************
!*********************************************************************************

SUBROUTINE removeondadireta(Nt,Nx,shot)

 IMPLICIT NONE

 INTEGER,INTENT(in)      :: Nx,Nt,shot

 REAL, DIMENSION(Nt,Nx)  :: Sismograma_Homogeneo,Sismograma_Real,Sismograma_sem_onda_direta !Seismogram


!Ler Sismograma Homogeneo

  CALL LoadSeismogram(Nt,Nx,shot,"Homogeneo","../sismograma_modelo_camada_de_agua/",Sismograma_Homogeneo)

!Ler Sismograma Real

  CALL LoadSeismogram(Nt,Nx,shot,"Marmousi","../sismograma/",Sismograma_Real)

! Retira a onda direta

Sismograma_sem_onda_direta = Sismograma_Real - Sismograma_Homogeneo

!grava sismograma sem onda direta

CALL Seismogram(Nt,Nx,shot,"Marmousi","../sismograma_sem_onda_direta/",Sismograma_sem_onda_direta)

END SUBROUTINE removeondadireta
 
