
# RTM_python-fortran [![Waffle.io - Columns and their card count](https://badge.waffle.io/GISIS-UFF/RTM_python-fortran.svg?columns=all&style=flat-square)](https://waffle.io/GISIS-UFF/RTM_python-fortran) 

Neste projeto foram abordadas duas técnicas fundamentais do processamento sísmico: A Modelagem
Sísmica e a Migração Reversa no Tempo (Reverse Time Migration - RTM). O objetivo principal será utilizar um código em Fortran para a execução dos cálculos que exigem grande poder de processamento enquanto que o código em Python será responsável pelos elementos visuais e interação com o usuário. 

Esse projeto visa a criação de um ambiente de desenvolvimento para a implementação da técnica Revere Time Migration que utiliza dados sísmicos pré Stack para a criação de uma imagem sísmica em profundidade. 

## Instruções para executar o programa

fonte: @Caicau

OBS: A migração utiliza a condição de imagem por tempo de excitação.

### Pré-requisitos 

* Python

* Compilador Fortran

* Antes de executar o seu programa, verifique/modifique o script "parametro.py" de acordo com a sua necessidade.

* Outra pasta que talvez você queira modificar é a "modelos_utilizados", lá você poderá adicionar o seu próprio arquivo com seu modelo em binário.

* Lembre-se sempre que para realizar a modelagem e a migração você precisará de um modelo real, um modelo suavizado e um modelo homogêneo (apresentando, normalmente, a velocidade da camada de água), todos em binário. Certifique-se de que você possui os três antes de começar. 

**Na pasta src:**

1º Abra o Terminal

2º Execute o script "prepare_fortran_subroutines.py". Esse script é o responsável pela criação do arquivo .so que será lido pelo Python.

$ python prepare_fortran_subroutines.py 

3º Verifique se o arquivo "fortransubroutines.so" foi criado.

4º Antes de realizar a modelagem e/ou a migração, verifique o script "parametro.py". Nele encontram-se todos os parametros necessários para a realização dos demais scripts. Normalmente, você irá mudar o tamanho e os nomes dos modelos.

5º No seu terminal (ainda na pasta src) execute o script "seismic_modeling.py". Ele é o responsável pela modelagem acústica.

$ python seismic_modeling.py

6º Os scripts remove_direct_wave.py e transit_matrix.py são os responsáveis pela remoção da onda direta e pela crição das matrizes de tempo de trânsito, respectivamente. 
 OBS: Eles devem ser feitos depois que já existir sismogramas do modelo real. E o transit_matrix.py deve ser rodado depois do remove_direct_wave.py
 
$ python remove_direct_wave.py
$ python transit_matrix.py

7º No seu terminal (ainda na pasta src) execute o script "migration.py". Ele é o responsável pela migração sísmica:

 $ python migration.py

## Instruções para manter repositório forked atualizado.
fonte: https://stackoverflow.com/questions/7244321/how-do-i-update-a-github-forked-repository

No seu repositório local:

1º Adicione a url ao remote name "upstream"

$ git remote add upstream https://github.com/GISIS-UFF/RTM_python-fortran

2º Atualize as informações do repositório remoto

$ git fetch upstream 

3º Verifique se está no branch master

$ git checkout master

4º Reescreva seu branch master - Atenção seus commits serão sobrescritos

$ git rebase upstream/master

Se não quiser sobreescrever os commits faça:

$ git merge upstream/master

Se você rebased seu branch em upstream/master você precisa forçar o push
para seu repositório forked

$ git push -f origin master




 
