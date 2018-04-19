def readbinaryfile(dim1,dim2,filename):
      """
      readbinaryfile - Functions that read a binary file.
      Usage
      Input:
      dim1     = Number of sample of 1st Dimension
      dim2     = Number of sample of 2nd Dimension
      filename = path of binary file
      Output:
      
      """
      from numpy import fromfile, reshape,float32

      with open(filename, 'rb') as f:    
            data   = fromfile(f, dtype=float32, count= dim1*dim2)
            matrix = reshape(data, [dim1,dim2], order='F')
      return matrix

 
def estabilidade(C,h,beta,dt):
    
    from numpy import max

    if dt > h / beta * max(max(C)):

       raise ValueError ("Erro de Estabilidade")

def dispersao(C,h,alfa,f_corte):
    
    from numpy import min   

    if h > min(min(C)) / (alfa * f_corte):

        raise ValueError ("Erro de Dispersao Numerica")   

   
def plotgraphics(ID,filename,color):
        
        """
    
            plotgraphics - Function that read a file and plot a graphic.
    
            REMEMBER: The function 'np.loadtxt' only can be used with files written
              with dot 
              
              ID = number of columns of the file
    
         """
        from numpy import loadtxt
        from matplotlib.pylab import figure, plot, draw

        if ID == 2:
             
            X,Y = loadtxt(filename, unpack = True)

            figure()
            plot(X,Y, c = color)
            draw()
        
        if ID == 1:
            
            X = loadtxt(filename, unpack = True)

            figure()
            plot(X, c = color)
            draw()
            
        return 
	
def plotmodel(matrix,colormap):

      from matplotlib.pylab import figure, imshow, colorbar, draw,cm

      figure()
      imshow(matrix,cmap= colormap)

      colorbar()
      draw() # drawing figure to be plotted later
      
def plotseism(Sismograma,Nt,Nx):
      
      from matplotlib.pylab import figure, imshow, draw, cm, xlabel, ylabel

      fig = figure()
      ax = fig.add_subplot(1, 1, 1)
      ax.xaxis.set_ticks_position("top")    
      imshow(Sismograma,aspect="auto",cmap = cm.gray, extent = [1,parametro.Nx,parametro.T,0])
      draw()


def plotsnaps(dim1,dim2,Nsnap,snap_path):

      from numpy import arange, fromfile, reshape, float32
      from matplotlib.pylab import figure, imshow, draw, cm, show
      import matplotlib.animation as animation

      fig2 = figure()
      shot = 1
      movie = []

      for snap in arange(1,Nsnap):
            inputfilename= snap_path + "shot" + '%03d'%(shot) + "snap" + '%03d'%(snap) + ".bin"
            data   = fromfile(inputfilename, dtype=float32, count= dim1*dim2)
            base   = reshape(data, [dim1,dim2], order='F')
            snapshot=imshow(base, cmap = cm.gray, animated=True)
            movie.append([snapshot])
      im_ani = animation.ArtistAnimation(fig2, movie, interval=250, repeat_delay=3000,
                                               blit=True)
  
      show()


def plotimage(imagepath,N_shot):


      from numpy import arange
  
      for shot in arange(1,N_shot+1):
            filename_imagem = imagepath + '%03d'%(shot) + ".bin"
            Imagem  =  readbinaryfile(parametro.Nz,parametro.Nx,filename_imagem)     
            StackImage = Imagem + StackImage
            # plotmodel(Imagem,cm.gray)
            plotmodel(StackImage,'jet')    

def amort(fat_amort,n_grid):

	"""
	amort - Cria a funcao de amortecimento que vai ser usada nas bordas do modelo e a salva em um arquivo de texto
	
	"""
	from numpy import zeros, arange, exp, savetxt

	w = zeros(n_grid)
  
	for i in arange(0,(n_grid)):
		w[i] = exp(-(fat_amort * (n_grid-i)) ** 2)
  
	  
	savetxt('f_amort.dat',w, delimiter='.')
        return w

def posicao_fonte(Nz,Nx,N_shot,Fx0,Fz0,SpaFonte):

      from numpy import zeros,arange,savetxt

      posicao = zeros((N_shot,2))

      Fz = arange(0,N_shot)*0 + Fz0
      Fx = arange(0,N_shot)*SpaFonte + Fx0


      # if (N_shot*SpaFonte + Fx0) > (Nx -1):
      #       raise ValueError ("Fonte fora do modelo valido (horizontal). Modifique ou o espacamento da fonte ou a posicao inicial Fx0")

      # if (N_shot*SpaFonte + Fz0) > (Nz -1):
      #       raise ValueError ("Fonte fora do modelo valido (vertical). Modifique a posicao inicial Fz0")
     
      posicao[0:N_shot,0] = Fx
      posicao[0:N_shot,1] = Fz

      savetxt("posicoes_fonte.dat",posicao,fmt = '%i')

def modelagem_acustica(regTTM,modelo_modelagem,ID_modelo):      
      

      from numpy import loadtxt,size,arange
      from matplotlib.pylab import cm
      from fortransubroutines import wavelet
      from fortransubroutines import nucleomodelagem
      from parametro import N_shot

  
      # Modelo de Velocidade Usado

      C = readbinaryfile(parametro.Nz,parametro.Nx,parametro.modeloreal)

      #plotmodel(C,'jet')    
      
      # Condicao de estabilidade e dispersao
      
      estabilidade(C,parametro.h,parametro.beta,parametro.dt)
      dispersao(C,parametro.h,parametro.alfa,parametro.f_corte)
      
      # Fonte Sismica
      
      wavelet(1,parametro.dt,1,parametro.f_corte) 
      
      #plotgraphics(2,'wavelet_ricker.dat', 'k')
      
      lixo, fonte = loadtxt('wavelet_ricker.dat', unpack = True)
      Nfonte      = size(fonte)
      
     # Funcao Amortecedora
      
      func_amort = amort(parametro.fat,parametro.nat)
      
      #plotgraphics(1,'f_amort.dat', 'b')
            
      # Modelagem

      # Funcao Posicao das Fontes

      if parametro.gera_pos_fonte == 1: 
            posicao_fonte(parametro.Nz,parametro.Nx,parametro.N_shot,parametro.Fx0,parametro.Fz0,parametro.SpaFonte)
      
      Fx, Fz = loadtxt('posicoes_fonte.dat',dtype = 'int',unpack = True)
      N_shot = size(Fx)

      if N_shot == 1:
            print "Fx =", Fx, "Fz =", Fz, "shot",N_shot
            nucleomodelagem(parametro.Nz,parametro.Nx,parametro.Nt,\
                                  parametro.h,parametro.dt,parametro.nat,\
                                  N_shot,parametro.shotshow,\
                                  Fx,Fz,fonte,parametro.Nsnap,regTTM,modelo_modelagem,parametro.zr,ID_modelo)
            # m1 = multiprocessing.Process(target=nucleomodelagem, args=(parametro.Nz,parametro.Nx,parametro.Nt,\
            #                       parametro.h,parametro.dt,parametro.nat,\
            #                       N_shot,parametro.shotshow,\
            #                       Fx,Fz,fonte,parametro.Nsnap,regTTM,modelo_modelagem,parametro.zr,ID_modelo,))
            # m1.start()
            # m1.join()
      else:
            for shot in arange(0,N_shot):
                  print "Fx =", Fx[shot], "Fz =", Fz[shot], "shot", shot+1

                  # m1 =multiprocessing.Process(target=nucleomodelagem, args = (parametro.Nz,parametro.Nx,parametro.Nt,\
                  #                 parametro.h,parametro.dt,parametro.nat,\
                  #                 shot+1,parametro.shotshow,\
                  #                 Fx[shot],Fz[shot],fonte,parametro.Nsnap,\
                  #                 regTTM,modelo_modelagem,parametro.zr,ID_modelo,))

                  nucleomodelagem(parametro.Nz,parametro.Nx,parametro.Nt,\
                                  parametro.h,parametro.dt,parametro.nat,\
                                  shot+1,parametro.shotshow,\
                                  Fx[shot],Fz[shot],fonte,parametro.Nsnap,\
                                  regTTM,modelo_modelagem,parametro.zr,ID_modelo)
                  
                  # m1.start()
                  # m1.join()

      # SOCORRO: Valores de Nsnap e Nfonte estao trocados mas funcionando mesmo assim :o
      # Esse problema esta na linha 5 do codigo em fortran

      #Problema Resolvido: Olhar o codigo em fortran: da linha 10 a linha 14!


      # if regTTM == 0:

      #       for shot in arange(1,N_shot+1):
      #             filename_sismograma_real = "../sismograma/Marmousi_sismograma" + '%03d'%(shot) + ".bin"
            
      #             Sismograma_Real = readbinaryfile(parametro.Nt,parametro.Nx,filename_sismograma_real)
      #             plotseism(Sismograma_Real,parametro.T,parametro.Nx)
      
      
      # if regTTM == 1:

      #       for shot in arange(1,N_shot+1):

      #             filename_matriz_tempo_transito = "../matriz_tempo_transito/Marmousi_" + 'shot' + '%03d'%(shot) + ".bin"
      #             matriz_tempo_transito = readbinaryfile(parametro.Nz,parametro.Nx, filename_matriz_tempo_transito)
      
      #             plotmodel(matriz_tempo_transito,'jet')

      # if parametro.shotshow > 0:
      #       plotsnaps(parametro.Nz,parametro.Nx,parametro.N_snap)
  
def remove_onda_direta():
      
      from fortransubroutines import removeondadireta
      from numpy import arange,loadtxt,size
      from parametro import N_shot
      
      Fx, Fz = loadtxt('posicoes_fonte.dat',dtype = 'int',unpack = True)
      N_shot = size(Fx)

      if N_shot == 1:
            print "Fx =", Fx, "Fz =", Fz, "shot",N_shot
            removeondadireta(parametro.Nt,parametro.Nx,N_shot)
      else:
            for shot in arange(0,N_shot):             
                  print "Fx =", Fx[shot], "Fz =", Fz[shot], "shot", shot+1
                  removeondadireta(parametro.Nt,parametro.Nx,shot+1)


      # for shot in arange(1,N_shot + 1): 
      #       filename_sismograma_camada_agua = "../sismograma_modelo_camada_de_agua/"+'Homogeneo_sismograma'+'%03d'%(shot) + '.bin'
      #       filename_sismograma_com_onda_direta = "../sismograma/Marmousi_sismograma" + '%03d'%(shot) + ".bin"
      #       filename_sismograma_sem_onda_direta = "../sismograma_sem_onda_direta/Marmousi_sismograma" + '%03d'%(shot) + ".bin"           

            # Sismograma_Real = readbinaryfile(parametro.Nt,parametro.Nx,filename_sismograma_com_onda_direta)
            # plotseism(Sismograma_Real,parametro.T,parametro.Nx)
            
            # Sismograma_Camada_Agua = readbinaryfile(parametro.Nt,parametro.Nx,filename_sismograma_camada_agua)
            # plotseism(Sismograma_Camada_Agua,parametro.T,parametro.Nx)
            
            # Sismograma =  readbinaryfile(parametro.Nt,parametro.Nx,filename_sismograma_sem_onda_direta)
            # plotseism(Sismograma,parametro.T,parametro.Nx)
  

def migracao_rtm(modelo_migracao):

      from matplotlib.pylab import cm
      from fortransubroutines import migracao
      from numpy import arange, loadtxt,size,zeros
      from parametro import N_shot

      StackImage = zeros((parametro.Nz,parametro.Nx))

      Fx, Fz = loadtxt('posicoes_fonte.dat',dtype = 'int',unpack = True)
      N_shot = size(Fx)

      if N_shot == 1:
            print "Fx =", Fx, "Fz =", Fz, "shot",N_shot
            migracao(parametro.Nz,parametro.Nx,parametro.Nt,parametro.h,parametro.dt,parametro.nat,\
                           parametro.zr,N_shot,parametro.shotshow,parametro.Nsnap,modelo_migracao)

      
      else:      
            for shot in arange(0,N_shot):
                  print "Fx =", Fx[shot], "Fz =", Fz[shot], "shot", shot+1
                  migracao(parametro.Nz,parametro.Nx,parametro.Nt,parametro.h,parametro.dt,parametro.nat,\
                           parametro.zr,shot+1,parametro.shotshow,parametro.Nsnap,modelo_migracao)
   
      for shot in arange(1,N_shot+1):
            filename_imagem = "../Imagem/Imagem_Marmousi_shot" + '%03d'%(shot) + ".bin"
            Imagem  =  readbinaryfile(parametro.Nz,parametro.Nx,filename_imagem)     
            StackImage = Imagem + StackImage
            #plotmodel(Imagem,cm.gray)
      
      plotmodel(StackImage,'jet')

      # if parametro.shotshow > 0:
      #       plotsnaps(parametro.Nz,parametro.Nx )
 

if __name__ == '__main__':
    
      """
      Structure prepared for object-oriented programming
      """
    
      from matplotlib.pylab import show
      import time
      import parametro
      from numpy import arange
      import multiprocessing 


      regTTM = 0

      start_time = time.time()
      
      ID_modelo = 1
      print "Modelagem_Sismogramas_Modelo_Real"
      m1 = multiprocessing.Process(target=modelagem_acustica, args = (regTTM,parametro.modeloreal,ID_modelo,))
      #modelagem_acustica(regTTM,parametro.modeloreal,ID_modelo)
     
      ID_modelo = 2

      print "Modelagem_Sismogramas_Camada_de_Agua"
      #modelagem_acustica(regTTM,parametro.modelocamadadeagua,ID_modelo) 
      m2 = multiprocessing.Process(target=modelagem_acustica, args = (regTTM,parametro.modelocamadadeagua,ID_modelo,))

      regTTM =1
      ID_modelo = 0

      print "Modelagem_Matriz_de_Tempo_de_Transito" 
      #modelagem_acustica(regTTM,parametro.modelosuavizado,ID_modelo)
      m3 = multiprocessing.Process(target=modelagem_acustica, args = (regTTM,parametro.modelosuavizado,ID_modelo,))

      m1.start()
      m2.start()
      m3.start()

      m1.join()
      m2.join()
      m3.join()
    
      print "Removendo Onda Direta"
      remove_onda_direta()

      print "Migracao"
      migracao_rtm(parametro.modelosuavizado)

      elapsed_time_python = time.time() - start_time

      print ("Tempo de processamento python = ", elapsed_time_python, "s")
    
      show() # Showing all figures draw
