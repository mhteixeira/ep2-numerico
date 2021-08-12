const loadTrussesData = (arquivo) => {
  arquivo = splitLines(arquivo)
  dados = arquivo[0].split(' ')
  nosFixos = int(dados[0]);
  nosLivres = int(dados[1]);
  barras = int(dados[2])
  
  posicoes = {
    1: [0, 0]
  }
  matrizDeBarras = []
  for(i = 2; i <= barras+1;i++) {
    infosDaBarra = arquivo[i].split(' ');
    primeiroNo = infosDaBarra[0];
    segundoNo = infosDaBarra[1];
    matrizDeBarras[i-2] = [primeiroNo, segundoNo]
    posicao_x = [0, 25,15,35,25,15,5,35,25,15,5,25,15,30,10]
    posicao_y = [0, 45,45,35,35,35,35,25,25,25,25,15,15,5,5]
  }

  return [nosFixos, nosLivres, barras, posicao_x, posicao_y, matrizDeBarras];
}