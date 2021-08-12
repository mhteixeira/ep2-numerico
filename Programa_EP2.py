# EP2 - Algoritmo QR para matrizes simétricas
# Dupla: 
# Lucas Domingues Boccia, NUSP 11262320
# Murillo Hierocles Alves de Sá Teixeira, NUSP 11325264

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#---------------------------------------------------------------#
#-------------------- CÓDIGOS TIRADOS DO EP1 -------------------#
#---------------------------------------------------------------#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#-----------------#
#---Importações---#
#-----------------#
 
import numpy as np

#-------------------------#
#    Rotação de Givens    #
#-------------------------#
 
def calculaQ(A, k):
  ak = A[k-1][k-1]
  bk = A[k][k-1]
  if abs(ak) > abs(bk):
      tau = -bk/ak
      ck = 1/(1+tau**2)**.5
      sk = ck*tau
  else:
      tau = -ak/bk
      sk = 1/(1+tau**2)**.5
      ck = sk*tau
  Q = np.identity(len(A))
  Q[k-1][k-1] = ck
  Q[k][k-1] = sk 
  Q[k][k] = ck
  Q[k-1][k] = -sk
  return Q
 
#--------------------#
#    Fatoração QR    #
#--------------------#
 
def fatoracaoQR(A):
    Qant = np.zeros(len(A))
    R = A
    Q = np.identity(len(A))
    for i in range(len(A)-1):
        Qant = calculaQ(R, i + 1)
        R = np.matmul(Qant, R)
        Q = np.matmul(Q, Qant.T)
    return[Q, R]
 
#-------------------------------#
#    Heurística de Wilkinson    #
#-------------------------------#
 
def calculaMi(A, m):
  ak_ant = A[m-2][m-2]
  ak = A[m-1][m-1]
  dk = (ak_ant - ak)/2
  sgn_dk = 1 if dk >=0 else -1
  bk = A[m-1][m-2]
  uk = ak + dk - sgn_dk*(dk**2 + bk**2)**.5
  return uk
 
#----------------------#
#    O Algoritmo QR    #
#----------------------#
 
# Ele se encontra modificado em relação à versão do EP1
# para que consiga recever uma matriz inicial para V diferente
# da identidade

def algoritmoQR(A, eps, H, desloc):
    Ld = A
    n = len(A)
    V = H
    k = 0
    uk = 0
    AutoValores = np.zeros(n)
    AutoVetores = np.zeros(n**2).reshape((n, n))
    for i in range(n-1):
      m = n - i
      while abs(Ld[m-1][m-2]) >= eps: 
        if (k > 0 and desloc==True):
          uk = calculaMi(Ld, m)
        Desloc = uk*np.identity(m)
        [Q, R] = fatoracaoQR(Ld - Desloc)
        Ld = np.matmul(R, Q) + Desloc
        k += 1
        V[:, :m] = np.matmul(V[:, :m], Q)
      AutoValores[i] = Ld[-1][-1]
      Ld = Ld[:len(Ld)-1,:len(Ld)-1]
    AutoValores[i+1] = Ld[0]
    AutoValores = np.flip(AutoValores)
    AutoVetores = V
 
    return [AutoValores, AutoVetores, k]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#---------------------------------------------------------------#
#------------------------ INÍCIO DO EP2 ------------------------#
#---------------------------------------------------------------#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#-------------------------------#
#    Cálculo dos vetores w_i    #
#-------------------------------#

def create_wi(A):
    col1 = A[:,0].copy() 
    signal = -1 if col1[1] >=0 else 1
    factor = signal*np.sqrt(col1[1:].dot(col1[1:]))
    w = col1; w[0] = 0; w[1] = col1[1] - factor
    return w

#--------------------------------------#
#    Cálculo da matriz H transposta    #
#--------------------------------------#

def createHT(X, w):
    IH = np.zeros(X.T.shape)
    for i in range(0, len(X.T)):
        col = X[:, i].copy()
        IH[i] = col - (2*w.dot(col)/w.dot(w))*w
    return IH
    
#------------------------------------#
#    Transformação de Householder    #
#------------------------------------#

def householder(A, w):
    col1 = A[:,0].copy()
    col1 = col1 - (2*(w.dot(col1))/w.dot(w))*w
    col1[2:] = np.zeros(len(col1)-2)
    
    AH = np.zeros(A.shape)
    AH[0] = col1
    for i in range(len(A)-1):
        col = A[:, i+1]
        AH[i+1] = col - (2*(w.dot(col))/w.dot(w))*w
    return AH

#--------------------------------------#
#    Algoritmo de tridiagonalização    #
#--------------------------------------#

def tridiagonalize(A):
    n = len(A)
    T = A.copy()
    HT = np.eye(n)

    for i in range(n-2):
        w = create_wi(T[i:, i:])
        tempT = householder(T[i:, i:], w)
        T[i:, i:] = householder(tempT, w)
        HT[:, i:] = createHT(HT[:, i:].T, w)

    return T, HT

#-----------------------------------#
#    Carga dos dados da tarefa 1    #
#-----------------------------------#

def readMatrix(fileName):
    f = open(fileName, 'r')
    n = int(f.readline().strip())
    M = []
    for i in range(n):
        line = f.readline().strip().replace('  ', ' ').split(' ')
        M.append([float(element) for element in line])
    f.close()
    M = np.array(M)
    return M
    
#---------------------------#
#    Solução da tarefa 1    #
#---------------------------#

def task1(question):
    if(question == 'A'):
        M = readMatrix('./inputs/input-a')

    if(question == 'B'):
        M = readMatrix('./inputs/input-b')

    print('Matriz em análise = ')
    print(M)
    T, H = tridiagonalize(M)
    L, V, _ = algoritmoQR(T, 1e-6, H, True)
    print('\nAutovalores calculados = ')
    print(np.around(L, 5))
    print('\nAutovetores calculados = ')
    print(np.around(V, 5))
    for i in range(len(M)):
        print(f'\nVerificação Av = lv para l = {np.around(L[i], 5)}')
        print(f'\n     A x {i+1}° autovetor  = ')
        print('     ', np.around(M@V[:, i], 5))
        print(f'\n     {i+1}° autovalor x {i+1}° autovetor  = ')
        print('     ', np.around(L[i]*V[:, i], 5))
    print('\nVerificação da ortogonalidade: V x Vt = I => Vt = V^(-1)')
    print('V x Vt = ')
    print(np.around((np.matmul(V, V.T)), 5))

#---------------------------------------------#
#    Carga dos dados de input das treliças    #
#---------------------------------------------#

def loadTrussesData(fileName):
    f = open(fileName, "r")
    totalJoints, freeJoints, numberOfBars = f.readline().split(" ")[:3]
    density, A, E = f.readline().split(" ")[:3]
    params = [int(totalJoints), 
                int(freeJoints), 
                int(numberOfBars),
                float(density), 
                float(A), 
                float(E)*1e9]
    
    lengths = []
    angles = []
    bars = []
    
    for i in range(int(numberOfBars)):
        infos = f.readline().rstrip().split(" ")
        bars.append(
            [int(infos[0]), int(infos[1])]
        )
        angles.append(float(infos[2]))
        lengths.append(float(infos[3]))
    f.close()

    return params, lengths, angles, bars

#-----------------------------------#
#    Criação da matriz de Massas    #
#-----------------------------------#

def createMassMatrix(numberOfFreeJoints, numberOfBars, density, A, bars, lengths):
    M = np.zeros(numberOfBars)
    
    for i in range(numberOfBars):
        M[2*(bars[i][0]-1)] += 0.5*A*density*lengths[i]
        M[2*(bars[i][0]-1)+1] += 0.5*A*density*lengths[i]
        M[2*(bars[i][1]-1)] += 0.5*A*density*lengths[i]
        M[2*(bars[i][1]-1)+1] += 0.5*A*density*lengths[i]
    return M[:-4]

#------------------------------------------------#
#    Criação da matriz de rididez para barras    #
#------------------------------------------------#

def createBarStiffnessMatrix(A, E, angle, length):
    C = np.cos(np.deg2rad(angle))
    S = np.sin(np.deg2rad(angle))
    vectorCS = np.array([-C, -S, C, S]).reshape(4, 1)
    
    K = (A*E/length)*vectorCS@vectorCS.T
    return K

#-----------------------------------------------------#
#    Criação da matriz de rididez para a estrutura    #
#-----------------------------------------------------#

def createTrussesStiffnessMatrix(params, bars, lengths, angles):
    numberOfFreeJoints = params[1]
    numberOfBars = params[2]
    A = params[4]
    E = params[5]
    
    K = np.zeros((2*numberOfFreeJoints, 2*numberOfFreeJoints))
    for i in range(numberOfBars):
        firstJoint = bars[i][0];
        secondJoint = bars[i][1];
        length, angle = [lengths[i], angles[i]]
        Kij = createBarStiffnessMatrix(A, E, angle, length)
    
        K[2*(firstJoint-1), 2*(firstJoint-1)] += Kij[0, 0]
        K[2*(firstJoint-1)+1, 2*(firstJoint-1)] += Kij[1, 0]
        K[2*(firstJoint-1), 2*(firstJoint-1)+1] += Kij[0, 1]
        K[2*(firstJoint-1)+1, 2*(firstJoint-1)+1] += Kij[1, 1]
    
        if(secondJoint <= numberOfFreeJoints):
            K[2*(firstJoint-1), 2*(secondJoint-1)] += Kij[0, 2]
            K[2*(firstJoint-1)+1, 2*(secondJoint-1)] += Kij[1, 2]
            K[2*(firstJoint-1), 2*(secondJoint-1)+1] += Kij[0, 3]
            K[2*(firstJoint-1)+1, 2*(secondJoint-1)+1] += Kij[1, 3]
        
            K[2*(secondJoint-1), 2*(secondJoint-1)] += Kij[2, 2]
            K[2*(secondJoint-1)+1, 2*(secondJoint-1)] += Kij[3, 2]
            K[2*(secondJoint-1), 2*(secondJoint-1)+1] += Kij[2, 3]
            K[2*(secondJoint-1)+1, 2*(secondJoint-1)+1] += Kij[3, 3]
            
            K[2*(secondJoint-1), 2*(firstJoint-1)] += Kij[2, 0]
            K[2*(secondJoint-1)+1, 2*(firstJoint-1)] += Kij[3, 0]
            K[2*(secondJoint-1), 2*(firstJoint-1)+1] += Kij[2, 1]
            K[2*(secondJoint-1)+1, 2*(firstJoint-1)+1] += Kij[3, 1]
    return K

# ------------------------------#
#     Criação da matriz Ktil    #
# ------------------------------#

def createKtilmatrix(K, M):
    Ktil = np.zeros(K.shape, dtype='float')
    for i in range(len(K)):
        Ktil[i] = np.sqrt(1/M[i])*K[i]
    for i in range(len(K)):
        Ktil[:, i] = np.sqrt(1/M[i])*Ktil[:, i]
    return Ktil

# ---------------------------#
#     Solução da tarefa 2    #
# ---------------------------#

def task2():
    # Solução da equação da Tarefa 2
    # Import dos parâmetros dos arquivos
    params, lengths, angles, bars = loadTrussesData("./inputs/input-c")
    numberOfJoints = params[0]
    numberOfFreeJoints = params[1]
    numberOfBars = params[2]
    density = params[3]
    A = params[4]
    
    # Parâmetros da solução da equação:
    # MX" + KX = 0
    M = createMassMatrix(numberOfFreeJoints, numberOfBars,  density, A, bars, lengths)
    K = createTrussesStiffnessMatrix(params, bars, lengths, angles)

    Ktil = createKtilmatrix(K, M)

    T, H = tridiagonalize(Ktil)

    L, Y, _ = algoritmoQR(T, 1e-6, H, desloc=True)
    w = np.sqrt(L)

    Z = np.zeros(Y.shape)
    for i in range(len(K)):
        Z[i] = np.sqrt(1/M[i])*Y[i]

    n = 5
    freqs = np.sort(w)[0:n]
    modosDeVibracao = np.array([Z[:, list(w).index(freqs[0])],
                                Z[:, list(w).index(freqs[1])],
                                Z[:, list(w).index(freqs[2])],
                                Z[:, list(w).index(freqs[3])],
                                Z[:, list(w).index(freqs[4])]])

    for i in range(n):
        print(f'\nFrequência {i+1} =', np.around(freqs[i], 3), 'rad/s')
        print(f'Modo de vibração {i+1} = ')
        print(modosDeVibracao[i])

    resposta = input('\nDeseja conferir os dados para animação (S/N)? ')
    if resposta.upper().strip()=='S':
        createOutputFilesForAnimation(freqs, modosDeVibracao)


#-------------------------------------------#
#    Criação de arquivos para a animação    #
#-------------------------------------------#

def createOutputFilesForAnimation(freqs, modosDeVibracao):
        modo = input('\nEscolha um modo de vibração (1, 2, 3, 4 ou 5): ')
        numeroDeFiguras = input('\nQual o número de figuras que se deseja criar? ')
        tempoEntreFiguras = input('\nQual o tempo entre figuras (em segundos)? ')

        modo = int(modo)
        numeroDeFiguras = int(numeroDeFiguras) + 1
        tempoEntreFiguras = float(tempoEntreFiguras)
        

        tempo = np.linspace(0, numeroDeFiguras*tempoEntreFiguras, numeroDeFiguras)
        x = np.zeros(24*(len(tempo)+1)).reshape((24, len(tempo)+1))
        
        headerStr = "       i"
        x[:, 0] = range(1, 25)
        for t in range(len(tempo)):
            headerStr += f'   X({t+1:05d})'
            print()
            x[:, t+1] = (modosDeVibracao[modo-1]*np.cos(freqs[modo-1]*tempo[t])).T

        # Criação dos arquivos que serão utilizados pela aplicação web para a execução da animação
        
        # Arquivo para a verificação do algoritmo utilizado
        np.savetxt('./outputs/x(t)', x, delimiter=" ", fmt=("%10.7f"), header=headerStr)

        #Arquivo JavaScript com os modos de vibração (para evitar o uso de bibliotecas adicionais!)
        np.savetxt('./outputs/modos_de_vibracao.js', modosDeVibracao, delimiter=" ")
        with open('./outputs/modos_de_vibracao.js','r') as contents:
            save = contents.read()
        with open('./outputs/modos_de_vibracao.js','w') as contents:
            contents.write("let modosDeVibracao = `")
        with open('./outputs/modos_de_vibracao.js','a') as contents:
            contents.write(save)
            contents.write('`')
        
        #Arquivo JavaScript com as frequências de vibração (para evitar o uso de bibliotecas adicionais!)
        np.savetxt('./outputs/frequencias.js', freqs, delimiter=" ")
        with open('./outputs/frequencias.js','r') as contents:
            save = contents.read()
        with open('./outputs/frequencias.js','w') as contents:
            contents.write("let frequencias = `")
        with open('./outputs/frequencias.js','a') as contents:
            contents.write(save)
            contents.write('`')
        
        #Arquivo JavaScript com o input-c (para evitar o uso de bibliotecas adicionais!)
        file = open("./outputs/input-c.js", "w") 
        with open('./inputs/input-c','r') as contents:
            save = contents.read()
        with open("./outputs/input-c.js", "w") as contents:
            contents.write("let inputC = `")
        with open("./outputs/input-c.js",'a') as contents:
            contents.write(save)
            contents.write('`')

#-------------------------------#
#    Interface com o usuário    #
#-------------------------------#

def userInterface():
  print('EP2 - Algoritmo QR para matrizes simétricas')
  print('Dupla:')
  print(' Lucas Domingues Boccia, NUSP 11262320')
  print(' Murillo Hierocles Alves de Sá Teixeira, NUSP 11325264')
  print('\n-------------')
  querAnalisar = True

  while(querAnalisar):
    print('\nQuestões:')
    print('A. Teste da tridiagonalização, item A')
    print('B. Teste da tridiagonalização, item B')
    print('C. Aplicação: Treliças Planas')
    print('\n-------------')
    
    questaoInvalida = True
    while(questaoInvalida):
        questao = input('\nQual a questão de interesse (A/B/C)? ')
        print('\n')

        questaoInvalida = False
        if (questao.upper().strip() == 'A' or questao.upper().strip() == 'B'):
            task1(questao.upper().strip())

        elif (questao.upper().strip() == 'C'): 
            task2()

        else:
            questaoInvalida = True
            print('Questão inválida!')

    print('\n-------------')
    resposta = input('\nDeseja testar outra questão (S/N)? ')
    querAnalisar = True if resposta.upper().strip()=='S' else False
  print('\n-------------')
  print('\nFim da execução!')

userInterface()