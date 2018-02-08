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
	
def plotmodel(matrix):

      from matplotlib.pylab import figure, imshow, colorbar, draw

      figure()
      imshow(matrix,cmap='jet')
      colorbar()
      draw() # drawing figure to be plotted later
      
def plotseism(Sismograma,Nt,Nx):
      
      from matplotlib.pylab import figure, imshow, draw, cm, xlabel, ylabel

      fig = figure()
      ax = fig.add_subplot(1, 1, 1)
      ax.xaxis.set_ticks_position("top")    
      imshow(Sismograma,aspect="auto",cmap = cm.gray, extent = [1,parametro.Nx,parametro.T,0])
      xlabel('Canal') 
      ylabel('Tempo (segundos)')
      draw()


def plotsnaps(dim1,dim2):

      from numpy import arange, fromfile, reshape, float32
      from matplotlib.pylab import figure, imshow, draw, cm, show
      import matplotlib.animation as animation

      fig2 = figure()
      shot = 1
      movie = []

      for snap in arange(1,21):
            inputfilename="../snapshot/Marmousi_" + "shot" + '%03d'%(shot) + "snap" + '%03d'%(snap) + ".bin"
            data   = fromfile(inputfilename, dtype=float32, count= dim1*dim2)
            base   = reshape(data, [dim1,dim2], order='F')
            snapshot=imshow(base, cmap = cm.gray, animated=True)
            movie.append([snapshot])
      im_ani = animation.ArtistAnimation(fig2, movie, interval=250, repeat_delay=3000,
                                               blit=True)
  
      show()

def amort(fat_amort,n_grid):

	"""
	amort - Cria a funcao de amortecimento que vai ser usada nas bordas do modelo e a salva em umarquivo de texto
	
	"""
	from numpy import zeros, arange, exp, savetxt

	w = zeros(n_grid)
  
	for i in arange(0,(n_grid)):
		w[i] = exp(-(fat_amort * (n_grid-i)) ** 2)
  
	  
	savetxt('f_amort.dat',w, delimiter='.')
        return w
        
def modelagemacustica(regTTM):      
      '''
      Main program is here      
      '''
      from numpy import loadtxt,size
      from fortransubroutines import wavelet
      from fortransubroutines import nucleomodelagem
      
      
      # Modelo de Velocidade

     # C = readbinaryfile(parametro.Nz,parametro.Nx,parametro.modeloreal)
      C = readbinaryfile(parametro.Nz,parametro.Nx,parametro.modelohomogeneo)
      plotmodel(C)    
      
      # Condicao de estabilidade e dispersao
      
      estabilidade(C,parametro.h,parametro.beta,parametro.dt)
      dispersao(C,parametro.h,parametro.alfa,parametro.f_corte)
      
      # Fonte Sismica
      
      wavelet(1,parametro.dt,1,parametro.f_corte) 
      
      plotgraphics(2,'wavelet_ricker.dat', 'k')
      
      lixo, fonte = loadtxt('wavelet_ricker.dat', unpack = True)
     # Nfonte      = size(fonte)
      Nfonte = 2364

      # Funcao Amortecedora
      
      func_amort = amort(parametro.fat,parametro.nat)
      
      plotgraphics(1,'f_amort.dat', 'b')
            
      # Modelagem
      
      Fx = int(parametro.Nx/2)               # Posicao da Fonte (x)
      Fz = 5  #int(parametro.Nz/2) #10       # Posicao da Fonte (z)


      # print 'parametro.Nz       ' ,parametro.Nz       
      # print 'parametro.Nx       ' ,parametro.Nx       
      # print 'parametro.Nt       ' ,parametro.Nt       
      # print 'parametro.h        ' ,parametro.h
      # print 'parametro.dt       ' ,parametro.dt       
      # print 'parametro.nat      ' ,parametro.nat     
      # print 'parametro.shot     ' ,parametro.shot     
      # print 'parametro.shotshow ' ,parametro.shotshow
      # print 'Fx                 ' ,Fx
      # print 'Fz                 ' ,Fz
      # print 'Nfonte             ' ,Nfonte
      # print 'parametro.Nsnap    ' ,parametro.Nsnap                          
      # print np.size(fonte),np.size(func_amort)

      nucleomodelagem(parametro.Nz,parametro.Nx,parametro.Nt,\
                      parametro.h,parametro.dt,parametro.nat,\
                      parametro.shot,parametro.shotshow,\
                      Fx,Fz,fonte,parametro.Nsnap,regTTM)
      

      # SOCORRO: Valores de Nsnap e Nfonte estao trocados mas funcionando mesmo assim :o
      # Esse problema esta na linha 5 do codigo em fortran

      #Problema Resolvido: Olhar o codigo em fortran: linha 10 a linha 14!

      if regTTM == 0:
            Sismograma_Real = readbinaryfile(parametro.Nt,parametro.Nx,"../sismograma/Marmousi_sismograma001.bin")

            plotseism(Sismograma_Real,parametro.T,parametro.Nx)

            Sismograma_H = readbinaryfile(parametro.Nt,parametro.Nx,"../sismograma_modelo_camada_de_agua/Homogeneo_sismograma001.bin")
      
            plotseism(Sismograma_H,parametro.T,parametro.Nx)

      
      if regTTM == 1:
            matriz_tempo_transito = readbinaryfile(parametro.Nz,parametro.Nx,"../matriz_tempo_transito/Marmousi_shot001.bin")
      
            plotmodel(matriz_tempo_transito)

      if parametro.shotshow > 0:
            plotsnaps(parametro.Nz,parametro.Nx) 

  
if __name__ == '__main__':
    
      """
      Structure prepared for object-oriented programming
      """
    
      from matplotlib.pylab import show
      from fortransubroutines import removeondadireta
      import time
      import parametro
      

      regTTM = 1

      start_time = time.time()


      modelagemacustica(regTTM) # Call main function

      removeondadireta(parametro.Nt,parametro.Nx,parametro.shot)
      Sismograma =  readbinaryfile(parametro.Nt,parametro.Nx,"../sismograma_sem_onda_direta/Marmousi_sismograma001.bin")
      plotseism(Sismograma,parametro.T,parametro.Nx)
      
      #migracao()

      elapsed_time_python = time.time() - start_time
      print ("Tempo de processamento python = ", elapsed_time_python, "s")
    
      show() # Showing all figures draw

