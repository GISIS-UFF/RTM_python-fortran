# RTM_python-fortran

Esse projeto visa a criação de um ambiente de desenvolvimento para a implementação da técnica Revere Time Migration que utiliza dados sísmicos pré Stack para a criação de uma imagem sísmica em profundidade.

Incluir pré requisitos para execução do programa

Incluir instruções para execução do programa.

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

## Instruções para executar o programa

OBS0: Até o momento a migração utiliza a condição de imagem por tempo de excitação. 

fonte: Cintia

OBS1: Antes de executar o seu programa verifique o script parametro.py e modifique os parâmetros apresentados de acordo com a sua necessidade.

OBS2: Outra pasta que talvez você queira modificar é "modelos_utilizados", lá você poderá adicionar o seu próprio arquivo com seu modelo em binário.

OBS3: Lembre-se sempre que para realizar a modelagem e a migração você precisará de um modelo normal, um modelo suavizado e um modelo homogêneo(apresentando a velocidade da camada de água). Certifique-se de que você possui os três antes de começar 

Na pasta src:

1º Abra o Terminal

2º Execute o script "prepare_fortran_subroutines.py". Esse script é o responsável pela criação do arquivo .so que será lido pelo Python.

$ python prepare_fortran_subroutines.py 

3º Verifique se o arquivo fortransubroutines.so foi criado

4º No seu terminal (ainda na pasta src) execute o script "programa_imageamento.py"

OBS 4: Antes de executar este comando abra o script e verifique se você quer realizar tanto a modelagem quanto a migração. Caso o seu interesse seja apenas em um dos dois métodos, certifique-se de comentar o outro no "main" do programa antes de executá-lo.

$ python programa_imageamento.py



 