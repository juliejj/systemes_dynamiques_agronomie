# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 18:05:44 2018

@author: casenave
"""

# ************ CODE POUR L'IDENTIFICATION DES PARAMETRES a ET k ***************
# *********************** D'UNE FONCTION DE MONOD *****************************

import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as scint
from matplotlib import patches as pat
import scipy.optimize as scop# pour utiliser des algorithmes d'optimisation déjà codés sous python

plt.close('all') # ferme toutes les figures

# 1. IDENTIFICATION LINEAIRE avec l'estimateur des moindres carrés linéaire
# -------------------------------------------------------------------------
def identifMonod(N,sigma1,sigma2):
    # N : nombre d'observation
    # sigma1 : ecart type bruit de mesure sur mu
    # sigma2 : ecart type bruit de mesure sur S
    
    # Generation de donnees d'observation
    # -----------------------------------
    # parametres du modele de Monod
    # mu(S)=k*S/(S+a)
    coeffk = 1.34 # coefficient k
    coeffa = 1.57 # coefficient a
    
    # calcul des mesures bruitées de S et de mu
    S = np.linspace(0,15,num=N) # mesures du substrat
    mu = coeffk*S/(coeffa+S) # sorties du modèle correspondantes
    mub = mu+sigma1*(np.random.rand(N)-0.5) # bruitage des sorties
    Sb = S + sigma2*(np.random.rand(N)-0.5) # bruitage des entrées
    
    # calcul du modèle exact
    Sabs = np.linspace(0,15,num=200) # mesures du substrat
    muexact = coeffk*Sabs/(coeffa+Sabs) # sorties du modèle correspondantes
    
    # trace des observations
    plt.figure(1)
    plt.plot(Sabs,muexact,'r',label='modèle exact')
    plt.plot(Sb,mub,'ro',label='mesures bruitées')
    plt.xlabel('Substrat')
    plt.ylabel('Taux de croissance')

    # Identification des parametres k et a du modele de Monod
    # -------------------------------------------------------
    beta = np.ones(N) # coefficients de pondération
    
    # calcul de l'estimateur de moindres carrés
    matM = np.zeros((2,2))
    vecb = np.zeros((2,1))
    for i in np.arange(0,N,1):
        phii = np.array([[S[i]], [-mub[i]]]) # regresseur
        matM = matM+beta[i]*(phii*np.transpose(phii)) # matrice à inverser
        vecb = vecb + beta[i]*phii*S[i]*mub[i] 
    theta_hat = np.linalg.solve(matM, vecb) #estimateur des moindres carrés
    
    # calcul du modèle identifié
    muhat = theta_hat[0]*Sabs/(theta_hat[1]+Sabs)  # calcul du mu correspondant
    
    # tracé du modèle identifié et des valeurs des paramètres k et a identifiés
    plt.plot(Sabs,muhat,'b',label='modèle identifié')
    plt.text(16.5,0.8,'coefficient k exact ='+str(coeffk)+'; estimé ='+str(theta_hat[0][0]))
    plt.text(16.5,0.7,'coefficient a exact ='+str(coeffa)+'; estimé ='+str(theta_hat[1][0]))
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
    plt.show()
  
N=74
sigma1=0.2
sigma2=0.8
identifMonod(N,sigma1,sigma2)

# 2. IDENTIFICATION NON LINEAIRE à l'aide de l'algorithme du gradient
# ------------------------------------------------------------------
# les algorithmes de gradient à pas optimal et à pas constant sont testés

# -----------------------------------------------------------------------------
#                      PARAMETRES POUVANT ETRE MODIFIES
# -----------------------------------------------------------------------------
# *** Paramètres des données
sigma1=0.2 #  ecart type bruit de mesure sur mu
sigma2=0.9 # ecart type bruit de mesure sur S
N =100 # nombre d'observation

# *** Paramètres des algorithmes de gradient
x0 = np.array([0.01,0.01]) # valeur initiale pour l'algorithme à pas optimal
x0c = np.array([0.01,0.01]) # valeur initiale pour l'algorithme à pas constant

x0 = np.array([1.9,4.5]) # valeur initiale pour l'algorithme à pas optimal
x0c = np.array([1.9,4.5]) # valeur initiale pour l'algorithme à pas constant

alpha_constant = 2 # valeur de alpha pour l'algorithme de gradient à pas constant

Nbiter = 35 # nombre d'itérations

precision = 0.001 # seuil de précision pour le critère d'arrêt des algorithmes

# *** Parametres du modèle de Monod
coeffk = 1.34 # coefficient k
coeffa = 1.57 # coefficient a

# ------------------------------------------------------------------------------

print('Critères d\'arret pour les algorithmes')
print('1. Nombre d\'itération maximal = '+str(Nbiter))
print('2. Seuil de précision = '+str(precision))

# *** Generation de donnees d'observation
S = np.linspace(0,15,num=N) # mesures du substrat
mu = coeffk*S/(coeffa+S) # sorties du modèle: calcul du mu
mub = mu+sigma1*(np.random.rand(N)-0.5) # bruitage des sorties mesurées
Sb = S + sigma2*(np.random.rand(N)-0.5) # bruitage des entrées mesurées
  
# *** fonction qui calcule la valeur de la fonction f à minimiser, ainsi que son gradient
def fonction_f(paramk,parama):
    global mub, Sb, N
    # f : fonction f
    # partialf_k : dérivée partielle de f par rapport à k
    # partialf_a : dérivée partielle de f par rapport à a
    # gradf : gradient de f
    
    # initialisation
    if len(paramk.shape)>0: # si les entrées sont vectorielles
        f = np.zeros((np.size(paramk,0),np.size(paramk,1)))
        partialf_k = np.zeros((np.size(paramk,0),np.size(paramk,1)))
        partialf_a = np.zeros((np.size(paramk,0),np.size(paramk,1)))
    else: # si les entrées sont scalaires
        f = 0
        partialf_k = 0
        partialf_a = 0
        
    # Boucle sur les données
    for i in np.arange(0,len(Sb)-1,1):
        mu = paramk*Sb[i]/(parama+Sb[i])
        # calcul de f qui est définie comme la somme sur les données (i=1:N) de
        # (k*S_i/(a+S_i)-mu_i)^2 divisée par N
        f = f + (mu-mub[i])**2
        # calcul de la dérivée partielle de f par rapport à k qui est définie comme
        # la somme sur les données (i=1:N) de 2*S_i/(a+S_i)*(k*S_i/(a+S_i)-mu_i) divisée par N
        partialf_k = partialf_k + 2*(mu-mub[i])*Sb[i]/(parama+Sb[i])
        # calcul de la dérivée partielle de f par rapport à a qui est définie comme
        # la somme sur les données (i=1:N) de -2*k*S_i/(a+S_i)^2*(k*S_i/(a+S_i)-mu_i) divisée par N
        partialf_a = partialf_a - 2*(mu-mub[i])*paramk*Sb[i]/((parama+Sb[i])**2)
    # division par N
    f = f/N
    partialf_k = partialf_k/N
    partialf_a = partialf_a/N
    # stockage du gradient
    gradf = np.array([partialf_k, partialf_a])
    return [f, gradf]

# *** Calcul de la fonction f pour le tracé
# on fait un maillage de l'espace des valeurs des paramètres k et a
# On se limite ici aux valeurs de k comprises entre 0.5 et 2
# et les valeurs de a comprises entre 0.5 et 5: ce choix est arbitraire
x , y = np.meshgrid(np.linspace(0.5,2,201),np.linspace(0.5,5,200))
z = fonction_f(x,y)
z = z[0]

# tracé des courbes de niveau de f, c'est à dire des courbes sur lesquelles f est constante
plt.figure(1)
graphe = plt.contour(x,y,z,400)
plt.xlabel('coefficient k')
plt.ylabel('coefficient a')

# fonction à minimiser dans l'algorithme de gradient à pas optimal pour obtenir la valeur de alpha
def fonction_falpha(alpha,x,gradfx):
    # dans l'algorithme de gradient à pas optimal, on cherche la valeur de alpha qui
    # minimise f(x_k-alpha gradf(x_k)) pour un x_k donné à chaque itération
    x1 = x - alpha*gradfx
    f = fonction_f(x1[0],x1[1])
    f = f[0]
    return f

# *** Algorithme du gradient

# 1. Gradient à pas optimal
# initialisation des vecteurs où seront stockées les valeurs intermédiaires obtenues avec les algorithmes itératifs
xval = np.zeros((Nbiter+1,2))

# stockage de la valeur initiale: le choix de cette valeur est effectué en début de code
xval[0]=x0
# le choix du nombre d'itérations ainsi que le seuil de précision du critère d'arrêt sont
# également effectués en début de code

# Tracé de la valeur initiale dans l'espace des paramètres
plt.plot(x0[0],x0[1],'ro')

# initialisation de la valeur du gradient de f
temp = fonction_f(np.asarray(x0[0]),np.asarray(x0[1]))
fx0 = temp[0] # valeur de f en x_k
gradfx0 = temp[1] # gradient de f en x_k
    
# Boucle itérative de l'algorithme de gradient à pas optimal
noiter = 1 # numéro de l'itération courante
while (noiter < Nbiter+1) and (np.sqrt(gradfx0[0]**2+gradfx0[1]**2)>precision):
     
    # calcul de la valeur de alpha optimale en utilisant ici une méthode de minimisation codée sous
    # python appelée algorithme de Nelder-Mead
    alphaopt = scop.minimize(fonction_falpha,0,args=(x0,gradfx0),method='Nelder-Mead')
    if alphaopt.success == True:
        alpha_k = alphaopt.x
    else:
        alpha_k = alpha_constant    
        
    # stockage et mise à jour des valeurs de paramètres
    xval[noiter] = x0 - alpha_k*gradfx0 # algorithme de gradient optimal
    
    # Tracé de la direction de descente dans l'espace des paramètres
    plt.plot(np.array([x0[0],xval[noiter][0]]),np.array([x0[1],xval[noiter][1]]),'r')
    
    # mise à jour pour l'itération d'après
    x0 = xval[noiter]
    noiter = noiter + 1
    
    # calcul de f et du gradient de f en x_k
    # pour le gradient à pas optimal
    temp = fonction_f(np.asarray(x0[0]),np.asarray(x0[1]))
    fx0 = temp[0] # valeur de f en x_k
    gradfx0 = temp[1] # gradient de f en x_k
    
    # Tracé des points dans l'espace des paramètres
    plt.plot(x0[0],x0[1],'r.')

# Affichage des légendes et des résultats
plt.plot(np.array([x0[0],xval[noiter-2][0]]),np.array([x0[1],xval[noiter-2][1]]),'r',label='algorithme du gradient à pas optimal')
plt.text(2.2,3.3,'** Gradient à pas optimal **')
plt.text(2.2,3.0,'nombre d\'itération ='+str(noiter-1))
plt.text(2.2,2.7,'norme du gradient ='+str(np.sqrt(gradfx0[0]**2+gradfx0[1]**2)))
plt.text(2.2,2.4,'coefficient k exact ='+str(coeffk)+'; estimé ='+str(x0[0]))
plt.text(2.2,2.1,'coefficient a exact ='+str(coeffa)+'; estimé ='+str(x0[1]))
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

# 2. Gradient à pas constant
# initialisation des vecteurs où seront stockées les valeurs intermédiaires obtenues avec les algorithmes itératifs
xvalc = np.zeros((Nbiter+1,2))

# stockage de la valeur initiale: le choix de cette valeur est effectué en début de code
xvalc[0]=x0c
# le choix du nombre d'itérations ainsi que le seuil de précision du critère d'arrêt sont
# également effectués en début de code

# Tracé de la valeur initiale dans l'espace des paramètres
plt.plot(x0c[0],x0c[1],'ro')

# initialisation de la valeur du gradient de f
temp = fonction_f(np.asarray(x0c[0]),np.asarray(x0c[1]))
fx0c = temp[0] # valeur de f en x_k
gradfx0c = temp[1] # gradient de f en x_k

# Boucle itérative de l'algorithme de gradient à pas constant
noiter = 1 # numéro de l'itération courante
while (noiter < Nbiter+1) and (np.sqrt(gradfx0c[0]**2+gradfx0c[1]**2)>precision):
        
    # stockage et mise à jour des valeurs de paramètres
    xvalc[noiter] = x0c - alpha_constant*gradfx0c # algorithme de gradient à pas constant
    
    # Tracé de la direction de descente dans l'espace des paramètres
    plt.plot(np.array([x0c[0],xvalc[noiter][0]]),np.array([x0c[1],xvalc[noiter][1]]),'g')
    
    # mise à jour pour l'itération d'après
    x0c = xvalc[noiter]
    noiter = noiter + 1
    
    # calcul de f et du gradient de f en x_k
    # pour le gradient à pas constant
    temp = fonction_f(np.asarray(x0c[0]),np.asarray(x0c[1]))
    fx0c = temp[0] # valeur de f en x_k
    gradfx0c = temp[1] # gradient de f en x_k  
    
    # Tracé des points dans l'espace des paramètres
    plt.plot(x0c[0],x0c[1],'g.')
    
# Affichage des légendes et des résultats
plt.plot(np.array([x0c[0],xvalc[noiter-2][0]]),np.array([x0c[1],xvalc[noiter-2][1]]),'g',label='algorithme du gradient à pas constant')
plt.plot(coeffk,coeffa,'ko',label='solution exacte')
plt.text(2.2,1.7,'**Gradient à pas constant**')
plt.text(2.2,1.4,'nombre d\'itération ='+str(noiter-1))
plt.text(2.2,1.1,'norme du gradient ='+str(np.sqrt(gradfx0c[0]**2+gradfx0c[1]**2)))
plt.text(2.2,0.8,'coefficient k exact ='+str(coeffk)+'; estimé ='+str(x0c[0]))
plt.text(2.2,0.5,'coefficient a exact ='+str(coeffa)+'; estimé ='+str(x0c[1]))
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.show()