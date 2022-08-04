from tkinter import messagebox
from PyQt5 import QtCore, QtGui, QtWidgets
import numpy as np
from scipy.spatial import distance
import math
import sys


# Trabalho 2 - PO II
# PNL - Multivariavel irrestrito
# Luciano Tanaka 
# Caio Regal 

# Executavel nao funcionando, entao preferi enviar apenas o codigo fonte para correcao dos metodos :)

n=2
# f = funcao,df = derivada da funcao, x0 = valor inicial, eps = delta, itmax = max de it
# Metodo de Newton para realizar busca na reta (minf(x) com x -> R)
def newton_mono(f, e):
    df = np.diff(exec(f))
    
    for n in range(0,20):
        f_xn = exec(f(xn))
        if abs(f_xn) < e:
            return xn
        df_xn = exec(df(xn))
        if df_xn == 0:
            return None
        xn = xn - f_xn/df_xn
    return None

def grad(f, x1, x2, e):

    x_i = np.array([])
    d_i = np.array([])
    l_i = np.array([])

    k=1 
    grad = np.gradient(f)
    
    #teste = list(grad)
    #print(teste)
    #arrumar gradient
    #mapped = list(map(np.gradient(), lambda x: f(x)))
    #print(mapped)

    while(distance.euclidean(grad[0], grad[1]) <= e):
        d_i[k] = -grad
        l = newton_mono(f, e) 
        x_i[k+1] =  x_i[k] + l*d_i[k]
        k+=1
    x_otimo = x_i[k]
    return x_otimo

# criteiro de parada
# grad = np.gradient(f[K+1])
# distance.euclidean(grad)

def err(string):
    print(string)
    input('Pressione tecla para sair')
    sys.exit()

# troca de linhas
def swapRows(v,i,j):
    if len(v.shape) == 1:
        v[i],v[j] = v[j],v[i]
    else:
        v[[i,j],:] = v[[j,i],:]
        
def swapCols(v,i,j):
    v[:,[i,j]] = v[:,[j,i]]

# Pivotamento de Gauss
def gaussPivot(a, b, e=1.0e-12):
    n = len(b)
    
    s = np.zeros(n)
    for i in range(n):
        s[i] = max(np.abs(a[i,:]))
        
    for k in range(0, n-1):
        
    #   Troca de linhas
        p = np.argmax(np.abs(a[k:n,k])/s[k:n]) + k
        if abs(a[p,k]) < e: 
            err('Matriz nao tem inversa')
        if p != k:
            swapRows(b,k,p)
            swapRows(s,k,p)
            swapRows(a,k,p)
            
    #   Eliminacao
        for i in range(k+1,n):
            if a[i,k] != 0.0:
                lam = a[i,k]/a[k,k]
                a[i,k+1:n] = a[i,k+1:n] - lam*a[k,k+1:n]
                b[i] = b[i] - lam*b[k]
    if abs(a[n-1, n-1]) < e: 
        err('Matriz nao tem inversa')
        
#   Colocando na inversa
    b[n-1] = b[n-1]/a[n-1,n-1]
    for k in range(n-2,-1,-1):
        b[k] = (b[k] - np.dot(a[k,k+1:n],b[k+1:n]))/a[k,k]
    return b

# Metodo de Newton multivariavel
def newton_multi(f, x1, x2, e):

    # Construcao da matriz jacobiana
    def jacobiano(f,x):
        h = 1.0e-4
        n = len(x)
        jac = np.zeros((n,n))
        f0 = exec(f(x))
        for i in range(x):
            aux = x[i]
            x[i] = aux + h
            f1 = f(x)
            x[i] = aux
            jac[:,i] = (f1 - f0)/h
        return [jac, f0]

    for i in range(25):
        jac, f0 = jacobiano(f, x1)
        if math.sqrt(np.dot(f0, f0)/len(x1) < e):
            return x1
        dx = gaussPivot(jac, -f0)
        x1+= dx
        if math.sqrt(np.dot(dx, dx) < e):
            return x1


def fletcher_reeves(f, x1, x2, e, n):
    x_i = np.array([])
    l_i = np.array([])
    d_i = np.array([])

    # calcular
    g0 = np.gradient(exec(f(x1)))
    d0 = -g0

    # passo 1
    grad = np.gradient(exec(f(x1)))
    if(distance.euclidean(grad[0], grad[1]) < e):
        return x1

    # passo 2
    for k in range(0, n-1):
        x_i[k+1] = x_i[k] + l_i[k]*d_i[k]
    return

def davidon_fletcher_powell(f, x1, x2, e, n):
    s = np.array(1,0,1,0)
    grad = np.gradient(exec(f(x1)))
    g  = np.array([])

    g[0] = grad
    k=0
    i=0

    q = np.array([])
    p = np.array([])
    d = np.array([])
    l = np.array([])
    x = np.array([])

    while(distance.euclidean(grad[0], grad[1]) >= e):
        d[k] = -s[k]*g[k]
        l[k] = newton_mono(f)
        x[k+1] = x[k] + l[k]*d[k]

        if (k < n-1):
            g[k+1] = np.gradient(f(x[k+1]))
            q[k] = g[k+1] - g[k]
            p[k] = l[k] * d[k]
            s[k+1] = s[k] + (p[k] * np.transpose(p[k])) / (np.transpose(p[k]) * q[k]) - (s[k] * q[k]*np.transpose(q[k]) * s[k]) / (np.transpose(q[k]) * s[k] * q[k])
            k +=1
        else:
            x[i+1] = x[n]
            i+=1
            x[0] = x[n]
            g[0] = np.gradient(exec(f(x[0])))
            k = 0
    return

class Ui_mainWindow(object):
    def setupUi(self, mainWindow):
        mainWindow.setObjectName("mainWindow")
        mainWindow.setWindowModality(QtCore.Qt.WindowModal)
        mainWindow.resize(750, 600)
        mainWindow.setStyleSheet("background-color: rgb(255, 255, 255);")
        self.centralwidget = QtWidgets.QWidget(mainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.botao_calcular = QtWidgets.QPushButton(self.centralwidget)
        self.botao_calcular.setGeometry(QtCore.QRect(240, 450, 111, 41))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.botao_calcular.setFont(font)
        self.botao_calcular.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.botao_calcular.setStyleSheet("background-color: rgb(232, 232, 232);")
        self.botao_calcular.setAutoRepeat(True)
        self.botao_calcular.setAutoDefault(True)
        self.botao_calcular.setDefault(False)
        self.botao_calcular.setFlat(False)
        self.botao_calcular.setObjectName("botao_calcular")
        self.metodos = QtWidgets.QComboBox(self.centralwidget)
        self.metodos.setGeometry(QtCore.QRect(260, 210, 301, 51))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        font.setBold(True)
        font.setWeight(75)
        self.metodos.setFont(font)
        self.metodos.setCursor(QtGui.QCursor(QtCore.Qt.PointingHandCursor))
        self.metodos.setStyleSheet("background-color: rgb(232, 232, 232);")
        self.metodos.setObjectName("metodos")
        self.metodos.addItem("")
        self.metodos.addItem("")
        self.metodos.addItem("")
        self.metodos.addItem("")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(30, 380, 111, 31))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(30, 220, 211, 31))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(-20, -10, 621, 211))
        font = QtGui.QFont()
        font.setFamily("Arial Black")
        font.setPointSize(16)
        font.setStrikeOut(False)
        font.setKerning(True)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(30, 280, 131, 31))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.form_e = QtWidgets.QLineEdit(self.centralwidget)
        self.form_e.setGeometry(QtCore.QRect(140, 280, 121, 31))
        self.form_e.setObjectName("form_e")
        self.form_func = QtWidgets.QLineEdit(self.centralwidget)
        self.form_func.setGeometry(QtCore.QRect(140, 380, 401, 31))
        self.form_func.setObjectName("form_func")
        self.out_x = QtWidgets.QLineEdit(self.centralwidget)
        self.out_x.setGeometry(QtCore.QRect(90, 570, 131, 41))
        self.out_x.setObjectName("out_x")
        self.line = QtWidgets.QFrame(self.centralwidget)
        self.line.setGeometry(QtCore.QRect(0, 450, 591, 31))
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        self.label_7.setGeometry(QtCore.QRect(30, 600, 51, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        self.label_8.setGeometry(QtCore.QRect(310, 600, 151, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.label_8.setFont(font)
        self.label_8.setObjectName("label_8")
        self.out_fx = QtWidgets.QLineEdit(self.centralwidget)
        self.out_fx.setGeometry(QtCore.QRect(390, 590, 121, 41))
        self.out_fx.setObjectName("out_fx")
        self.label_10 = QtWidgets.QLabel(self.centralwidget)
        self.label_10.setGeometry(QtCore.QRect(140, 510, 311, 51))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.label_9 = QtWidgets.QLabel(self.centralwidget)
        self.label_9.setGeometry(QtCore.QRect(280, 280, 71, 31))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.form_x1 = QtWidgets.QLineEdit(self.centralwidget)
        self.form_x1.setGeometry(QtCore.QRect(360, 280, 51, 31))
        self.form_x1.setObjectName("form_x1")
        self.form_x2 = QtWidgets.QLineEdit(self.centralwidget)
        self.form_x2.setGeometry(QtCore.QRect(420, 280, 51, 31))
        self.form_x2.setObjectName("form_x2")
        self.label_11 = QtWidgets.QLabel(self.centralwidget)
        self.label_11.setGeometry(QtCore.QRect(480, 280, 71, 31))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")
        self.out_x_2 = QtWidgets.QLineEdit(self.centralwidget)
        self.out_x_2.setGeometry(QtCore.QRect(90, 620, 131, 41))
        self.out_x_2.setObjectName("out_x_2")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(140, 340, 311, 31))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(14)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.metodos.raise_()
        self.label.raise_()
        self.label_2.raise_()
        self.label_3.raise_()
        self.label_6.raise_()
        self.form_e.raise_()
        self.form_func.raise_()
        self.out_x.raise_()
        self.line.raise_()
        self.label_7.raise_()
        self.label_8.raise_()
        self.out_fx.raise_()
        self.label_10.raise_()
        self.label_9.raise_()
        self.form_x1.raise_()
        self.botao_calcular.raise_()
        self.form_x2.raise_()
        self.label_11.raise_()
        self.out_x_2.raise_()
        self.label_5.raise_()
        mainWindow.setCentralWidget(self.centralwidget)
        self.statusbar = QtWidgets.QStatusBar(mainWindow)
        self.statusbar.setObjectName("statusbar")
        mainWindow.setStatusBar(self.statusbar)

        # Tamanho da tela fixo
        mainWindow.setFixedSize(mainWindow.height(), mainWindow.width())

        # Desabilitar tela se necessario
        #self.metodos.currentIndexChanged.connect(self.modifica_campos)

         # rotina caso botao calcular seja clickado
        self.botao_calcular.clicked.connect(self.calcular)


        self.retranslateUi(mainWindow)
        QtCore.QMetaObject.connectSlotsByName(mainWindow)

    # funcao do botao calcular
    def calcular(self):
        ops = self.metodos.currentText()

        if(ops == "Gradiente"):
            try:
                f = self.form_func.text()
                x1 = float(self.form_x1.text())
                x2 = float(self.form_x2.text())
                e = float(self.form_e.text())

                output = grad(f, x1, x2, e)
                self.out_x.setText(self.out_x.text() + output)
                self.out_fx.setText(exec(f(output)))

            except ValueError:
                messagebox.showwarning('Erro ao calcular!')

        elif ops  == "Newton":
            try:
                f = self.form_func.text()
                x1 = float(self.form_x1.text())
                x2 = float(self.form_x2.text())
                e = float(self.form_e.text())

                output = newton_multi(f, x1, x2, e)
                self.out_x.setText(self.out_x.text() + output)
                self.out_fx.setText(exec(f(output)))

            except ValueError:
                messagebox.showwarning('Erro ao calcular!!')
                
        elif ops  == "Fletcher & Reeves":
            try:
                f = self.form_func.text()
                x1 = float(self.form_x1.text())
                x2 = float(self.form_x2.text())
                e = float(self.form_e.text())

                output = fletcher_reeves(f, x1, x2, e, n)
                self.out_x.setText(self.out_x.text() + output)
                self.out_fx.setText(exec(f(output)))

            except ValueError:
                messagebox.showwarning('Erro ao calcular!!!')

        elif ops  == "Davidon-Fletcher-Powell":
            try:
                f = self.form_func.text()
                x1 = float(self.form_x1.text())
                x2 = float(self.form_x2.text())
                e = float(self.form_e.text())

                output = davidon_fletcher_powell(f, x1, x2, e, n)
                self.out_x.setText(self.out_x.text() + output)
                self.out_fx.setText(exec(f(output)))

            except ValueError:
                messagebox.showwarning('Erro ao calcular!!!!')


    def retranslateUi(self, mainWindow):
        _translate = QtCore.QCoreApplication.translate
        mainWindow.setWindowTitle(_translate("mainWindow", "PNL : MONOVARIÁVEL"))
        self.botao_calcular.setText(_translate("mainWindow", "Calcular"))
        self.metodos.setItemText(0, _translate("mainWindow", "Gradiente"))
        self.metodos.setItemText(1, _translate("mainWindow", "Newton"))
        self.metodos.setItemText(2, _translate("mainWindow", "Fletcher & Reeves"))
        self.metodos.setItemText(3, _translate("mainWindow", "Davidon-Fletcher-Powell"))
        self.label.setText(_translate("mainWindow", "<html><head/><body><p>Min f(x) = </p></body></html>"))
        self.label_2.setText(_translate("mainWindow", "<html><head/><body><p><span style=\" font-size:14pt;\">Selecione o método:</span></p></body></html>"))
        self.label_3.setText(_translate("mainWindow", "<html><head/><body><p align=\"center\"><span style=\" font-size:18pt; font-weight:600;\">PO II - T2</span></p><p align=\"center\">PROGRAMAÇÃO NÃO LINEAR:</p><p align=\"center\">MULTIVARIÁVEL</p></body></html>"))
        self.label_6.setText(_translate("mainWindow", "<html><head/><body><p>Com  E =</p></body></html>"))
        self.form_e.setToolTip(_translate("mainWindow", "Tolerância"))
        self.form_func.setToolTip(_translate("mainWindow", "Função Objetivo"))
        self.out_x.setToolTip(_translate("mainWindow", "Ponto mínimo"))
        self.label_7.setText(_translate("mainWindow", "x* ="))
        self.label_8.setText(_translate("mainWindow", "f(x*) = "))
        self.out_fx.setToolTip(_translate("mainWindow", "Número de Iterações"))
        self.label_10.setText(_translate("mainWindow", "<html><head/><body><p align=\"center\"><span style=\" font-size:14pt; font-weight:600;\">Ponto minímo encontrado:</span></p></body></html>"))
        self.label_9.setText(_translate("mainWindow", "<html><head/><body><p>X_0 = (</p></body></html>"))
        self.form_x1.setToolTip(_translate("mainWindow", "Tolerância"))
        self.form_x2.setToolTip(_translate("mainWindow", "Tolerância"))
        self.label_11.setText(_translate("mainWindow", "<html><head/><body><p>)</p></body></html>"))
        self.out_x_2.setToolTip(_translate("mainWindow", "Ponto mínimo"))
        self.label_5.setText(_translate("mainWindow", "<html><head/><body><p>ex: (x1 - 2)**4 + (x1 - 3*x2)**2</p></body></html>"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    mainWindow = QtWidgets.QMainWindow()
    ui = Ui_mainWindow()
    ui.setupUi(mainWindow)
    mainWindow.show()
    sys.exit(app.exec_())

