# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 18:05:44 2018

@author: casenave
"""

# ************ CODE POUR L'ALGORITHME de NEWTON *******************************

import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as scint
from matplotlib import patches as pat
import scipy.optimize as scop# pour utiliser des algorithmes d'optimisation déjà codés sous python

plt.close('all') # ferme toutes les figures


# Illustration du principe de l'algorithme de Newton pour trouver les zéros d'une fonction
# ----------------------------------------------------------------------------------------
def fonction_f(S):
    # parametres du modele de Monod
    coeffk = 1.34 # coefficient k
    coeffa = 1.57 # coefficient a
    QsurV = 1 # valeur de Q/V
    
    # calcul de mu(S)
    mu = coeffk*S/(coeffa+S)
    # calcul de la dérivée de mu
    muprime = coeffk*coeffa/((coeffa+S)**2)
    # calcul de f=mu(S)-QsurV
    f = mu-QsurV
    # calcul de la dérivée de f
    fprime = muprime
    return [f, fprime]

# Calcul de la fonction f=mu(S)-Q/V pour le tracé
N = 200 # nombre de points pour le tracé
S = np.linspace(0,7,num=N) # mesures du substrat
f = fonction_f(S) # fonction f

def plot_iter(i):
    global xval
    # ième itération de l'algorithme de Newton
    x0 = xval[i] # valeur initiale du zéro
    temp = fonction_f(x0) # calcul de f=mu(S)-Q/V et sa dérivée
    fx0 = temp[0] # fonction f
    fprimex0 = temp[1] # dérivée de f
    x1 = x0 - fx0/fprimex0 # mise à jour de la valeur du zéro à l'itération i
    
    # tracé
    plt.figure(1)
    plt.plot(S,f[0],'b')
    plt.plot(x0,fx0,'ko')
    plt.text(x0*0.93,-0.98,'$x_{'+str(i)+'}$',fontsize=14)
    plt.plot(x1,0,'ko')
    plt.text(x1*1.02,-0.98,'$x_{'+str(i+1)+'}$',fontsize=14)
    plt.plot(np.array([0, 3.7]), fprimex0*(np.array([0, 3.7])-x0)+fx0,'r',
             label='tangente à $f$ en $x_{'+str(i)+'}$ d\'équation $y=f^\prime(x_{'+str(i)+'})(x-x_{'+str(i)+'})+f(x_{'+str(i)+'})$')
    plt.text(3.1,0.3,'tangente à $f$ en $x_{'+str(i)+'}$' ,fontsize=12)
    plt.text(3.1,0.2,'$y=f^\prime(x_{'+str(i)+'})(x-x_{'+str(i)+'})+f(x_{'+str(i)+'})$',fontsize=12)
    plt.plot(np.array([0, 7]), np.array([0, 0]),'k--')
    plt.plot(np.array([x0, x0]), np.array([-1, fx0]),'k--')
    plt.plot(np.array([x1, x1]), np.array([-1, 0]),'k--')
    plt.xlabel('Substrat $S$')
    plt.ylabel('$f(S)=\mu(S)-1$')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylim(-1,0.4)
    plt.show()

# Algorithme complet de Newton
Nbiter = 8 # nombre d'itérations
x0 = 0.5 # valeur initiale du zéro
xval = np.zeros((Nbiter+1,1)) # initialisation
xval[0]=x0 # stockage de la valeur initiale
for i in np.arange(1,Nbiter,1): # boucle d'itération
    temp = fonction_f(x0) # calcul de f=mu(S)-Q/V et de sa dérivée
    fx0 = temp[0] # fonction f
    fprimex0 = temp[1] # dérivée de f
    xval[i] = x0 - fx0/fprimex0 # mise à jour de la valeur du zéro à l'itération i
    x0 = xval[i] #stockage de la nouvelle estimation du zéro à l'itération i
    # Illustration du principe de l'algorithme de Newton pour trouver les zéros d'une fonction
# ----------------------------------------------------------------------------------------
def fonction_f(S):
    # parametres du modele de Monod
    coeffk = 1.34 # coefficient k
    coeffa = 1.57 # coefficient a
    QsurV = 1 # valeur de Q/V
    
    # calcul de mu(S)
    mu = coeffk*S/(coeffa+S)
    # calcul de la dérivée de mu
    muprime = coeffk*coeffa/((coeffa+S)**2)
    # calcul de f=mu(S)-QsurV
    f = mu-QsurV
    # calcul de la dérivée de f
    fprime = muprime
    return [f, fprime]

# Calcul de la fonction f=mu(S)-Q/V pour le tracé
N = 200 # nombre de points pour le tracé
S = np.linspace(0,7,num=N) # mesures du substrat
f = fonction_f(S) # fonction f

def plot_iter(i):
    global xval
    # ième itération de l'algorithme de Newton
    x0 = xval[i] # valeur initiale du zéro
    temp = fonction_f(x0) # calcul de f=mu(S)-Q/V et sa dérivée
    fx0 = temp[0] # fonction f
    fprimex0 = temp[1] # dérivée de f
    x1 = x0 - fx0/fprimex0 # mise à jour de la valeur du zéro à l'itération i
    
    # tracé
    plt.figure(1)
    plt.plot(S,f[0],'b')
    plt.plot(x0,fx0,'ko')
    plt.text(x0*0.93,-0.98,'$x_{'+str(i)+'}$',fontsize=14)
    plt.plot(x1,0,'ko')
    plt.text(x1*1.02,-0.98,'$x_{'+str(i+1)+'}$',fontsize=14)
    plt.plot(np.array([0, 7]), fprimex0*(np.array([0, 7])-x0)+fx0,'r')
    plt.text(3.1,0.45,'itération '+str(i),fontsize=12)
    plt.text(3.1,0.3,'tangente à $f$ en $x_{'+str(i)+'}$' ,fontsize=12)
    plt.text(3.1,0.2,'$y=f^\prime(x_{'+str(i)+'})(x-x_{'+str(i)+'})+f(x_{'+str(i)+'})$',fontsize=12)
    plt.plot(np.array([0, 7]), np.array([0, 0]),'k--')
    plt.plot(np.array([x0, x0]), np.array([-1, fx0]),'k--')
    plt.plot(np.array([x1, x1]), np.array([-1, 0]),'k--')
    plt.xlabel('Substrat $S$')
    plt.ylabel('$f(S)=\mu(S)-Q/V$')
    plt.ylim(-1,0.4)
    plt.show()

# Algorithme complet de Newton
Nbiter = 8 # nombre d'itérations
x0 = 0.5 # valeur initiale du zéro
xval = np.zeros((Nbiter+1,1)) # initialisation
xval[0]=x0 # stockage de la valeur initiale
for i in np.arange(1,Nbiter+1,1): # boucle d'itération
    temp = fonction_f(x0) # calcul de f=mu(S)-Q/V et de sa dérivée
    fx0 = temp[0] # fonction f
    fprimex0 = temp[1] # dérivée de f
    xval[i] = x0 - fx0/fprimex0 # mise à jour de la valeur du zéro à l'itération i
    x0 = xval[i] #stockage de la nouvelle estimation du zéro à l'itération i
    plot_iter(i)
