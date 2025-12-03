import tkinter as tk
from tkinter import *
from tkinter import ttk
import numpy as np
import time
import math

def lerArquivo(arquivo):
    try:
        with open(arquivo, 'r') as f:
            linhas = f.readlines()
            linhas = [linha.strip() for linha in linhas if linha.strip()]
            return linhas
    except Exception as e:
        print(f'Erro ao abrir o arquivo: {e}')
        return None

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.geometry('500x500')
        self.title('Aplicativo')

        self.formulas = ('Eliminação de gauss(sem pivoteamento)', 
                         'Eliminação de gauss(com pivoteamento)', 
                         'Pivoteamento completo', 
                         'Fatoração LU', 
                         'Fatoração Cholensky', 
                         'Gauss-Jacobi', 
                         'Gauss-Seidel',
                         'Bisseccao',
                         'MIL',
                         'Newton',
                         'Secante',
                         'Regula-falsi')

        self.opcao_var = tk.StringVar(self)
        self.resposta_var = tk.StringVar(self)

        self.inicializaWidgets()
    
    def inicializaWidgets(self):

        paddings = {'padx': 10, 'pady': 5}

        # Frame superior
        frame_top = tk.Frame(self)
        frame_top.pack(fill="x")

        # Label
        ttk.Label(frame_top, text='Selecione o método: ').grid(column=0, row=0, sticky=tk.W, **paddings)

        # OptionMenu
        ttk.OptionMenu(
            frame_top,
            self.opcao_var,
            self.formulas[0],
            *self.formulas,
        ).grid(column=1, row=0, sticky=tk.W, **paddings)

        # Botão resolve
        ttk.Button(frame_top, text='Resolve', command=self.resolve).grid(column=0, row=1, sticky=tk.W, **paddings)

        # Label de saída
        ttk.Label(frame_top, textvariable=self.resposta_var, foreground="red").grid(
            column=0, row=2, columnspan=2, sticky=tk.W, **paddings
        )

        # Tabela
        self.aplicativoFrame = tk.Frame(self)
        self.aplicativoFrame.pack(anchor="w", fill="x")

        scroll = tk.Scrollbar(self.aplicativoFrame, orient='vertical')
        self.tabela = ttk.Treeview(self.aplicativoFrame, yscrollcommand=scroll.set)
        self.tabela.pack(side=tk.LEFT, padx=(10,0))
        scroll.config(command=self.tabela.yview)

        
        
        scroll.pack(side=tk.LEFT, fill=tk.Y)


    def resolve(self):
        linhas = lerArquivo("variaveis.txt")
        for row in self.tabela.get_children():
            self.tabela.delete(row)
        resultado=''
        funcao_str = linhas[0]
        m = int(linhas[1])
        A = []
        for i in range(m):
            A.append([float(x) for x in linhas[2 + i].split()])
        A = np.array(A, dtype=float)
        k = int(linhas[2 + m])
        B = []
        for j in range(k):
            B.append([float(x) for x in linhas[3 + m + j].split()])
        B = np.array(B, dtype=float).reshape(-1, 1) 
        start = 3 + m + k
        a = float(linhas[start])
        b = float(linhas[start+1])
        delta = float(linhas[start+2])
        n = int(linhas[start+3])
        phi_str = linhas[start+4]
        derivada = linhas[start+5]
        eps = float(linhas[start+6])
        x0Secante = 1
        x1Secante = 2
        x0 = 1.5
        

        ponto_teste = (a + b) / 2
        resultado_teste = Funcao(ponto_teste, funcao_str)
        metodo = self.opcao_var.get()
    
        if A.shape[0] != A.shape[1]:
            resultado = "ERRO: MATRIZ NAO QUADRADA"
        elif A.shape[0] != B.shape[0]:
            resultado = "ERRO: MATRIZ B DE TAMANHO DIFERENTE"
        else:
            start_time = time.time()
            try:
                if(metodo == 'Eliminação de gauss(sem pivoteamento)'):
                    resultado = resolver_sistema_gauss(A.copy(), B.copy(), 'sem')
            except ValueError as e:
                resultado = (f"ERRO SEM PIVOTEAMENTO: {e}")
            try:
                if(metodo == 'Eliminação de gauss(com pivoteamento)'):
                    resultado = resolver_sistema_gauss(A.copy(), B.copy(), 'parcial')
            except ValueError as e:
                resultado = (f"ERRO PARCIAL: {e}")
            try:
                if(metodo == 'Pivoteamento completo'):
                    resultado = resolver_sistema_gauss(A.copy(), B.copy(),'completo')
            except ValueError as e:
                resultado = (f"ERRO COMPLETO: {e}")
            try:
                if(metodo == 'Fatoração Cholensky'):
                    resultado = cholensky(A, B)
            except ValueError as e:
                resultado = (f"ERRO COMPLETO: {e}")
            try:
                if(metodo == 'Fatoração LU'):
                    resultado = resolver_sistema_com_pivoteamento(A, B)
            except ValueError as e:
                resultado = (f"ERRO COMPLETO: {e}")
            try: 
                if(metodo == 'Gauss-Jacobi'):
                    x = jacobi(A, B, eps)
                    self.criaTabela(x)
            except ValueError as e:
                resultado = (f"ERRO COMPLETO: {e}")
            try:
                if(metodo == 'Gauss-Seidel'):
                    x = seidel(A, B, eps)
                    self.criaTabela(x)
            except ValueError as e:
                resultado = (f"ERRO COMPLETO: {e}")
            try: 
                if(metodo == 'Bisseccao'):
                    x = Bisseccao(funcao_str, a, b, delta, n)
                    self.criaTabela(x['resultados'])
            except ValueError as e:
                resultado = (f"ERRO COMPLETO: {e}")
            try: 
                if(metodo == 'MIL'):
                    x = MIL(funcao_str, phi_str, x0, delta, n)
                    self.criaTabela(x['resultados'])
            except ValueError as e:
                resultado = (f"ERRO COMPLETO: {e}")
            try:
                if(metodo == 'Newton'):
                    x = Newton(funcao_str, derivada, x0, delta, n)
                    self.criaTabela(x['resultados'])
            except ValueError as e:
                resultado = (f"ERRO COMPLETO: {e}")
            try: 
                if(metodo == 'Secante'):
                    x = Secante(funcao_str, x0Secante, x1Secante, delta, n)
                    self.criaTabela(x['resultados'])
            except ValueError as e:
                resultado = (f"ERRO COMPLETO: {e}")
            try:
                if(metodo == 'Regula-falsi'):
                    x = RegulaFalsi(funcao_str, a, b, delta, delta, n)
                    self.criaTabela(x['resultados'])
            except ValueError as e:
                resultado = (f"ERRO COMPLETO: {e}")
            end_time = time.time()  
            tempo_execucao = end_time - start_time
            exeecucao = f'\nTempo decorrido(ms): {tempo_execucao}'
            
        
        self.resposta_var.set(str(resultado) + exeecucao)
    def criaTabela(self, lista):

        if isinstance(lista[0], dict):
            cols = ["Iteração"] + list(lista[0].keys())
            self.tabela.config(columns=cols, show="headings")
            for col in cols:
                self.tabela.heading(col, text=col)
                self.tabela.column(col, width=80, anchor="center")
            for i, linha in enumerate(lista):
                valores = [i] + [linha[k] for k in linha] 
                self.tabela.insert("", "end", values=valores)
        else:  
            num_vars = len(lista[0])
            cols = ["Iteração"] + ["X" + str(i+1) for i in range(num_vars)]
            self.tabela.config(columns=cols, show="headings")
            for col in cols:
                self.tabela.heading(col, text=col)
                self.tabela.column(col, width=80, anchor="center")
            for i, linha in enumerate(lista):
                valores = [i] + linha
                self.tabela.insert("", "end", values=valores)


#PIVVVVOOOO

def pivot_parcial(Aa, j):
    n = Aa.shape[0]
    idx_relativo = np.argmax(np.abs(Aa[j:n, j]))
    idx_pivo = j + idx_relativo

    if idx_pivo != j:
        Aa[[j, idx_pivo]] = Aa[[idx_pivo, j]]
        return True
    return False


def pivot_completo(Aa, j, P):
    n = Aa.shape[0]
    submatrix = Aa[j:n, j:n]


    idx_max_flat = np.argmax(np.abs(submatrix))
    idx_row_rel, idx_col_rel = np.unravel_index(idx_max_flat, submatrix.shape)

    idx_row_abs = j + idx_row_rel
    idx_col_abs = j + idx_col_rel

    troca_ocorreu = False


    if idx_row_abs != j:
        Aa[[j, idx_row_abs]] = Aa[[idx_row_abs, j]]
        troca_ocorreu = True


    if idx_col_abs != j:

        Aa[:, [j, idx_col_abs]] = Aa[:, [idx_col_abs, j]]


        P[[j, idx_col_abs]] = P[[idx_col_abs, j]]
        troca_ocorreu = True

    return troca_ocorreu


def sistemaTriangularSuperior(U, y):
    n = U.shape[0]
    x = np.zeros((n, 1))

    for i in range(n - 1, -1, -1):
        soma_termos = np.dot(U[i, i + 1:n], x[i + 1:n, 0])

        if U[i, i] == 0:
            raise ValueError("Erro: Pivô zero. O sistema é singular.")


        x[i, 0] = (y[i, 0] - soma_termos) / U[i, i]

    return x




def resolver_sistema_gauss(A, b, pivoteamento='sem'):
    n = A.shape[0]

    Aa = np.concatenate((A, b), 1)


    P = np.arange(n)


    for j in range(n - 1):  # Coluna de eliminação


        if pivoteamento == 'parcial':
            pivot_parcial(Aa, j)
        elif pivoteamento == 'completo':
            pivot_completo(Aa, j, P)


        pivo = Aa[j, j]

        if np.isclose(pivo, 0):
            raise ValueError(f"Sistema Singular: Pivô zero encontrado (j={j}). \nNão há solução única.")


        for i in range(j + 1, n):
            fator = Aa[i, j] / pivo
            Aa[i, j:] = Aa[i, j:] - fator * Aa[j, j:]

    U = Aa[:, 0:n]
    y = Aa[:, n:]

    x_permut = sistemaTriangularSuperior(U, y)

    if pivoteamento == 'completo':
        x_solucao = np.zeros_like(x_permut)
        x_solucao[P] = x_permut
        return x_solucao

    return x_permut

def verificarQuadrada(matriz):
    return matriz.shape[0] == matriz.shape[1]
def verificarPositiva(matriz):
    if not np.allclose(matriz, matriz.T):
        return False
    autovalores = np.linalg.eigvals(matriz)
    return np.all(autovalores > 0)

def calcularTriangularInferior(matriz):
    n = matriz.shape[0]
    G = np.zeros_like(matriz)

    for i in range(n):
        for j in range(i + 1):
            if i == j:
                G[i, j] = np.sqrt(matriz[i, i] - np.sum(G[i, :j] ** 2))
            else:
                G[i, j] = (matriz[i, j] - np.sum(G[i, :j] * G[j, :j])) / G[j, j]

    return G

def transposta(matriz):
    return [[matriz[j][i] for j in range(len(matriz))] for i in range(len(matriz[0]))]

def multiplicacaoMatriz(matriz1, matriz2):
    return np.dot(matriz1, matriz2)

def verificarIgualdade(matrizResultado, matrizComparacao):
    return matrizComparacao.all() == matrizResultado.all()

def calcularMatrizColuna(matriz1, matriz2):

    if np.linalg.det(matriz1) == 0:
        return "\nA matriz G não é invertível."

    Y = np.linalg.solve(matriz1, matriz2)
    return Y

def cholensky(A, b):

    if verificarQuadrada(A):

        if np.allclose(A, A.T):
            if verificarPositiva(A):
                G = calcularTriangularInferior(A)
                G_T = transposta(G)
                resultado = multiplicacaoMatriz(G, G_T)
                if verificarIgualdade(resultado, A):
                    matriz_Y = calcularMatrizColuna(G, b)

                    if isinstance(matriz_Y, str):
                        raise ValueError("A matriz G não é invertível.")
                    else:
                        matriz_X = calcularMatrizColuna(G_T, matriz_Y)

                        if isinstance(matriz_X, str):
                            raise ValueError("A matriz G não é invertível.")
                        else:
                            return matriz_X
                else:
                    raise ValueError("O resultado da multiplicação G * G_T não é igual à matriz A.")

            else:
                raise ValueError("A matriz 'A' é simétrica, mas não é positiva.")
        else:
            raise ValueError("A matriz 'A' não é simétrica nem positiva.")
    else:
        raise ValueError("A matriz 'A' não é quadrada.")

def norma(V, X):
    n = len(V)
    maxnum = 0
    maxden = 0

    for i in range(0, n):
        num = abs(V[i] - X[i])
        if num > maxnum:
            maxnum = num
        den = abs(V[i])
        if den > maxden:
            maxden = den
    if den == 0:
        return None

    return maxnum/maxden

def jacobi(A, b, epsilon, iterMax=50):
    n = len(A)
    lista = list()
    x = n * [0]
    v = n * [0]

    for i in range(0, n):
        for j in range(0, n):
            if i!=j:
                A[i][j] = A[i][j]/A[i][i]
        b[i] = b[i]/A[i][i]
        x[i] = b[i]
    for k in range(1, iterMax+1):
        for i in range(0, n):
            S = 0
            for j in range(0, n):
                if i != j:
                    S = S + A[i][j] * x[j]
            v[i] = b[i] - S


        lista.append(v[:])
        
        d = norma(v, x)
        if d <= epsilon:
            return lista

        x = v[:]
    return lista

def seidel(A, b, epsilon, iterMax=50):
    n = len(A)
    lista = list()
    x = n * [0]
    v = n * [0]

    for i in range(0, n):
        for j in range(0, n):
            if i != j:
                A[i][j] = A[i][j] / A[i][i]
        b[i] = b[i] / A[i][i]
    for k in range(1, iterMax + 1):
        for i in range(0, n):
            S = 0
            for j in range(0, n):
                if i != j:
                    S = S + A[i][j] * x[j]
            x[i] = b[i] - S
        d = norma(x, v)
        lista.append(x[:])
        if d <=epsilon:
            return lista

        v = x[:]
    raise ValueError("Erro: Número máximo de iterações atingido")
    return x
def substituicoes_sucessivas(L, b):
    n = len(L)
    y = [0] * n
    for i in range(n):
        soma = 0
        for j in range(i):
            soma += L[i][j] * y[j]
        y[i] = (b[i] - soma) / L[i][i]
    return y


def substituicoes_retroativas(U, y):
    n = len(U)
    x = [0] * n
    for i in range(n - 1, -1, -1):
        soma = 0
        for j in range(i + 1, n):
            soma += U[i][j] * x[j]
        x[i] = (y[i] - soma) / U[i][i]
    return x


def identidade(n):
    return [[1 if i == j else 0 for j in range(n)] for i in range(n)]


def lu_com_pivoteamento(A):
    n = len(A)
    L = identidade(n)
    U = [linha[:] for linha in A]
    P = identidade(n)

    for k in range(n - 1):

        max_index = k
        for i in range(k + 1, n):
            if abs(U[i][k]) > abs(U[max_index][k]):
                max_index = i

        if max_index != k:
            U[k], U[max_index] = U[max_index], U[k]
            P[k], P[max_index] = P[max_index], P[k]
            for j in range(k):
                L[k][j], L[max_index][j] = L[max_index][j], L[k][j]


        for i in range(k + 1, n):
            if U[k][k] == 0:
                raise ValueError("Matriz singular")
            L[i][k] = U[i][k] / U[k][k]
            for j in range(k, n):
                U[i][j] -= L[i][k] * U[k][j]

    return P, L, U


def resolver_sistema_com_pivoteamento(A, b):

    P, L, U = lu_com_pivoteamento(A)


    n = len(b)
    b_perm = [0] * n
    for i in range(n):
        for j in range(n):
            if P[i][j] == 1:
                b_perm[i] = b[j]


    y = substituicoes_sucessivas(L, b_perm)


    x = substituicoes_retroativas(U, y)
    x = np.array(x).reshape(-1, 1)
    return x


def Funcao(x, expr):
    try:
        ambiente = {
            'x': x,
            'sin': math.sin,
            'cos': math.cos,
            'tan': math.tan,
            'log': math.log,
            'log10': math.log10,
            'exp': math.exp,
            'sqrt': math.sqrt,
            'pi': math.pi,
            'e': math.e,
            'abs': abs,
            'pow': pow
        }
        return eval(expr, {"__builtins__": {}}, ambiente)

    except Exception as e:
        print("Erro")
        return None


def Bisseccao(funcao_str, a, b, delta, n, arquivo_saida=None):
    k = 0
    meio = 0

    fa = Funcao(a, funcao_str)
    fb = Funcao(b, funcao_str)




    resultados = []

    while (abs(b - a) > delta and k < n):
        k += 1
        meio = (a + b) / 2
        fmeio = Funcao(meio, funcao_str)



        fa = Funcao(a, funcao_str)

        resultado_iteracao = {
            'a': a,
            'b': b,
            'meio': meio,
            'f_a': fa,
            'f_meio': fmeio,
            'erro': abs(b - a)
        }
        resultados.append(resultado_iteracao)


        if fa * fmeio < 0:
            b = meio
        else:
            a = meio
        print("\n")

    print("\n")
    if k == n:
        print(f"Número máximo de iterações ({n}) atingido")


    return {
        'raiz': meio,
        'f_raiz': Funcao(meio, funcao_str),
        'erro_final': abs(b - a),
        'resultados': resultados
    }


def MIL(funcao, phi, x0, delta, n, arquivo_saida=None):
    k = 0
    x_atual = x0
    resultados = []

    fx0 = Funcao(x0, funcao)


    if abs(fx0) < delta:
        return {
            'raiz': x0,
            'f_raiz': fx0,
            'erro_final': 0,
            'resultados': []
        }

    k = 1

    while k <= n:
        xk = Funcao(x_atual, phi)


        fx_novo = Funcao(xk, funcao)


        diferenca = abs(xk - x_atual)

        resultado_iteracao = {
            'x_anterior': x_atual,
            'x_atual': xk,
            'f_x_atual': fx_novo,
            'diferenca': diferenca
        }
        resultados.append(resultado_iteracao)

        if abs(fx_novo) < delta or diferenca < delta:
            break

        x_atual = xk
        k += 1

    return {
        'raiz': xk,
        'f_raiz': fx_novo,
        'erro_final': diferenca,
        'resultados': resultados
    }

def Newton(funcao, derivada, x0, delta, n, arquivo_saida=None):
    resultados = []
    fx0 = Funcao(x0, funcao)
    if abs(fx0) > delta:
        k = 1
        xant = x0
        flinha = Funcao(xant, derivada)


        xi = xant - (fx0/flinha)

        fx0 = Funcao(xi, funcao)


        resultado_iteracao = {
            'iteracao': k,
            'x_anterior': xant,
            'x_atual': xi,
            'f_x_anterior' : Funcao(xant, funcao),
            'f(x)_atual' : fx0,
            'derivada' : flinha,
            'diferenca' : abs(xi - xant)
        }

        resultados.append(resultado_iteracao)

        while (abs(fx0) > delta and abs(xi - xant) > delta and k < n):
            k += 1
            xant = xi
            flinha = Funcao(xant, derivada)
            xi = xant - (fx0/flinha)
            fx0 = Funcao(xi, funcao)

            resultado_iteracao = {
                'iteracao': k,
                'x_anterior': xant,
                'x_atual': xi,
                'f_x_anterior': Funcao(xant, funcao),
                'f(x)_atual': fx0,
                'derivada': flinha,
                'diferenca': abs(xi - xant)
            }

            resultados.append(resultado_iteracao)
        raiz = xi
    else:
        raiz = x0
        fx0 = Funcao(x0, funcao)
    return {
        'raiz': raiz,
        'f_raiz': fx0,
        'erro_final': abs(xi - xant) if 'x_i' in locals() else 0,
        'resultados': resultados
    }

def Secante(funcao, x0, x1, delta, n, arquivo_saida=None):
    resultados = []
    fx0 = Funcao(x0, funcao)
    fx1 = Funcao(x1, funcao)

    if abs(fx0) < delta:
        return {
            'raiz': x0,
            'f_raiz': fx0,
            'erro_final': 0,

            'resultados': []
        }
    if abs(fx1) < delta or abs(x1 - x0) < delta:
        return {
            'raiz': x1,
            'f_raiz': fx1,
            'erro_final': 0,
            'resultados': []
        }
    k = 1

    xant = x0
    xatual = x1
    fant = fx0
    fatual = fx1

    while k < n:
        if fatual - fant == 0:
            break

        x2 = xatual - (fatual * (xatual - xant)) / (fatual - fant)

        fx2 = Funcao(x2, funcao)

        resultado_iteracao = {
            'x0': xant,
            'x1': xatual,
            'x2': x2,
            'f(x0)': fant,
            'f(x1)': fatual,
            'f(x2)': fx2,
            'diferenca': abs(x2 - xatual)
        }

        resultados.append(resultado_iteracao)


        if abs(fx2) < delta or abs(x2 - xatual) < delta or k > n:
            break

        xant = xatual
        xatual = x2
        fant = fatual
        fatual = fx2

        k += 1
    raiz = x2

    return {
        'raiz': raiz,
        'f_raiz': fx2,
        'erro_final': abs(x2 - xatual),
        'resultados': resultados
    }


def RegulaFalsi(funcao_str, a, b, delta1, delta2, n, arquivo_saida=None):
    resultados = []
    fa = Funcao(a, funcao_str)
    fb = Funcao(b, funcao_str)

    if fa is None or fb is None:
        return None


    if fa * fb >= 0:
        return None

    if abs(b - a) < delta1:
        raiz = (a + b) / 2
        return {
            'raiz': raiz,
            'f_raiz': Funcao(raiz, funcao_str),
            'erro_final': abs(b - a),
            'resultados': []
        }

    if abs(fa) < delta2:
        return {
            'raiz': a,
            'f_raiz': fa,
            'erro_final': 0,
            'resultados': []
        }
    if abs(fb) < delta2:
        return {
            'raiz': b,
            'f_raiz': fb,
            'erro_final': 0,
            'resultados': []
        }

    k = 1
    a_atual = a
    b_atual = b
    fa_atual = fa
    fb_atual = fb

    while k <= n:

        x = (a_atual * fb_atual - b_atual * fa_atual) / (fb_atual - fa_atual)

        fx = Funcao(x, funcao_str)

        if fx is None:
            print(f"Erro ao calcular f(x) na iteração {k}")
            break

        resultado_iteracao = {
            'a': a_atual,
            'b': b_atual,
            'x': x,
            'f(a)': fa_atual,
            'f(b)': fb_atual,
            'f(x)': fx,
            'erro': abs(b_atual - a_atual)
        }

        resultados.append(resultado_iteracao)


        if abs(fx) < delta2:
            raiz = x
            break


        if fa_atual * fx < 0:
            b_atual = x
            fb_atual = fx
        else:
            a_atual = x
            fa_atual = fx

        if abs(b_atual - a_atual) < delta1:
            raiz = (a_atual + b_atual) / 2
            break

        k += 1

    if k > n:
        raiz = x
        print(f"Número máximo de iterações ({n}) atingido")

    return {
        'raiz': raiz,
        'f_raiz': fx,
        'erro_final': abs(b_atual - a_atual),
        'resultados': resultados
    }




if __name__ == "__main__":
    app = App()
    app.mainloop()