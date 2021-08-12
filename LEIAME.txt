## Instruções de execução e compilação ##

1) Linguagem e bibliotecas
O código foi desenvolvido em Python 3.7.
Para executar o algoritmo é necessário ter a biblioteca NumPy.

2) Instruções para execução
2.1) Ao rodar o código, será perguntado qual a questão de interesse.
	Deve-se responder A, B, C para selecionar a questão.
2.2) Caso seja escolhido A: 
	O código responderá a tarefa Testes do cálculo de autovetores e autovalores da matriz simétrica contida no item a) do enunciado
2.3) Caso seja escolhido B: 
	O código responderá a tarefa Testes do cálculo de autovetores e autovalores da matriz simétrica contida no item b) do enunciado
	
	Para os testes A e B serão printados: 
	 • A matriz em análise
	 • Os autovalores encontrados
	 • Os autovetores encontrados
	 • Rotina de verificação Av = lv
	 • Verificação da ortogonalidade da matriz de autovetores

2.4) Caso seja escolhido C:
	O código irá responder com as 5 menores frequências e seus respectivos modos de vibração para a estrutura apresentada no enunciado

	TAREFA BÔNUS:
	-> Após a execução da tarefa C, o usuário será questionado se deseja verificar o algoritmo usado para a criação da animação.
	-> Caso responda S, o usuário será questionado sobre o modo de vibração de interesse, o número de figuras
	e o intervalo de tempo entre elas e serão criados alguns arquivos na pasta ./outputs:
	 • Um arquivo contendo de forma legível os valores de x(t) para a verificação do monitor
	 • Arquivos que serão usados no input da aplicação criada.
	-> Para acessar e conseguir visualizar as animações, basta abrir o arquivo animation.html em seu navegador, 
	de preferência o Chrome, mas versões novas de todos devem suportá-lo.
	-> Ao acessar a aplicação será possível escolher o modo de vibração usando o slider superior
	-> Para controlar a amplitude dos movimentos, basta usar o slider inferior

3) Após executar um exercício, o programa irá perguntar se deseja executar novamente.
	Deve-se responder S/s ou N/n. Em caso S, realizará as operações explicadas acima. Em caso N, o algoritmo será finalizado.
	
