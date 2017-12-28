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
