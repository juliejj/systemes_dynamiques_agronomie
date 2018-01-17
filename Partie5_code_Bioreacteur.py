# -*- coding: utf-8 -*-
"""
Created on Sat Jan 13 18:05:44 2018

@author: casenave
"""

# ************ CODE POUR LA SIMULATION ET LE CONTROLE DU BIOREACTEUR **********

import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as scint
from matplotlib import patches as pat
plt.close('all') # ferme toutes les figures


# modèle de réacteur continu
def reacteur(x,t,k,muast,KS,KI,Qin,V,S0,Sast,control_type,coeffcontrol,disturb):
    # x : variables d'état du modèle c'est à dire B et S dans notre cas + 
    #       une variable supplémentaire dans le cas où il y a un terme intégral
    #       dans la loi de commande
    # t : temps
    # muast, KS, KI : paramètres de la fonction de Haldane utilisée pour le taux de croissance
    # Qin : débit d'entrée dans le réacteur
    # S0 : concentration en sucre dans le milieu qui alimente le réacteur
    # Sast : consigne en concentration en sucre pour la commande (= valeur que l'on veut atteindre)
    # control_type : type de loi de commande à appliquer. Pour l'instant un type est possible control_type='BO'
    #                pour la boucle ouverte
    # coeffcontrol : paramètres utilisés dans la loi de commande
    # disturb : perturbation sur la commande c'est à dire valeur telle que Qréel = Qcalc*(1+disturb)
    
    # récupération des valeurs des variables d'état
    B = x[0] # biomasse
    S = x[1] # substrat
    # Calcul de la commande: dans fonction_u la loi de commande est calculée, disturb permet de prendre
    # en compte d'éventuelles perturbation sur la valeur réellement appliquée de la commande
    Q = fonction_u(t,x,Sast,k,muast,KI,KS,V,S0,control_type,coeffcontrol)*(1+disturb)
    # initialisation de dx, second membre du modèle correspondant à la dérivée des variables d'état dB/dt et dS/dt
    if control_type in ['PI','PID']: # si il y a un terme intégrale dans la commande, on rajoute une équation pour
        # calculer l'intégrale Sast-S au cours du temps
        dx = np.zeros(3)
        # second membre de l'equation qui calcule l'intégrale de Sast-S qui sera stockée dans x[2]
        dx[2] = Sast-S
    else: # si il n'y a pas de terme intégrale dans la commande, on ne rajoute pas d'équation
        dx = np.zeros(2)
    # taux de croissance (fonction de Haldane)
    mu = muast*S/(KS+S+S**2/KI)
    # second membre de l'équation en B
    dx[0] = mu*B-Q/V*B
    # second membre de l'équatio en S
    dx[1] = -k*mu*B+Q/V*(S0-S)
    return dx


# ** fonction qui simule le modèle avec la loi de commande demandée et qui trace ensuite la solution
def culture_cont(Sast,control_type,coeffcontrol,disturb):
    
    # vecteur de temps pour la simulation
    tmax = 150
    temps = np.linspace(0,tmax,2000) 

    # paramètres du modèle
    k = 0.6; muast = 2.3; KS = 10; KI = 0.1; Qin = 0.01; V = 0.5; 
    
    # conditions initiales du modèle (valeurs initiales de la biomasse et de la concentration en sucre)
    B0 = 9; S0 = 3.2; 
    
    # si on utilise un terme intégrale, il faudra rajouter une équation dans le modèle et du coup rajouter
    # la condition initiale correspondante qui est égale à 0 (voir plus loin dans le paragraphe sur le terme
    # intégral)
    if control_type in ['I','PI','PID']: 
        x0 = np.array([B0,S0,0])
    else: 
        x0 = np.array([B0,S0])
        
    # intégration numérique de l'EDO
    x = scint.odeint(reacteur,x0,temps,args=(k,muast,KS,KI,Qin,V,S0,Sast,control_type,coeffcontrol,disturb))
    
    # re-calcul de la commande appliquée
    u = fonction_u(temps,x,Sast,k,muast,KI,KS,V,S0,control_type,coeffcontrol)

    # tracé des solutions
    plt.figure(figsize = (10, 3))
    plt.subplots_adjust(hspace=0.4,wspace=0.4)
    plt.subplot2grid((1,2),(0,0))
    plt.plot(temps,x[:,0],'r',label='Biomasse')
    plt.plot(temps,x[:,1],'g',label='Substrat')
    plt.plot(np.array([0,temps[-1]]),np.array([Sast,Sast]),'g--',label='S*')
    plt.legend(); plt.xlabel('time (h)')
    
    plt.subplot2grid((1,2),(0,1))
    plt.plot(temps,u,'r',label='Debit')
    plt.legend(); plt.xlabel('time (h)')
    plt.show()
    
# Loi de commande
def fonction_u(t,x,Sast,k,muast,KI,KS,V,S0,control_type,coeffcontrol):
    if control_type == 'BO': # BOUCLE OUVERTE
        # loi donnée par Qast=mu(Sast)*V
        Qast = muast*Sast/(KS+Sast+Sast**2/KI)*V
        if type(t)==float:  # cas où t est scalaire 
            valu = Qast
        else: # cas où t est un vecteur
            valu = np.ones(len(t))*Qast
    elif control_type == 'P': # BOUCLE FERMEE action PROPORTIONNELLE
        # récupération des paramètres de la loi de commande
        kprop = coeffcontrol
        # et de la valeur de S
        if type(t)==float: # cas où t est scalaire 
            valS = x[1]
        else: # cas où t est un vecteur
            valS = x[:,1]
        # loi donnée par mu(Sast)*V+kprop*(Sast-S)
        valu = muast*Sast/(KS+Sast+Sast**2/KI)*V+kprop*(Sast-valS)  
    elif control_type == 'PI': # BOUCLE FERMEE action PROPORTIONNELLE INTEGRALE
        # récupération des paramètres de la loi de commande
        kprop = coeffcontrol[0]
        kint = coeffcontrol[1]
        # et de la valeur de valint, qui est l'intégrale entre 0 et t de Sast-S
        if type(t)==float: # cas où t est scalaire 
            valint = x[2]; valS = x[1]
        else: # cas où t est un vecteur
            valint = x[:,2]; valS = x[:,1]
        # loi donnée par mu(Sast)*V+kprop*(Sast-S) + kint*valint
        # où valint est l'intégrale entre 0 et t de Sast-S
        Qast = muast*Sast/(KS+Sast+Sast**2/KI)*V
        valu = Qast+kprop*(Sast-valS)+kint*valint    
    elif control_type == 'PID': # BOUCLE FERMEE action PROPORTIONNELLE INTEGRALE DERIVEE
        # récupération des paramètres de la loi de commande
        kprop = coeffcontrol[0]
        kint = coeffcontrol[1]
        kderiv = coeffcontrol[2]
        # et des valeurs de valint (intégrale de Sast-S), de S et de B
        if type(t)==float: # cas où t est scalaire 
            valint = x[2]; valS = x[1]; valB=x[0]
        else: # cas où t est un vecteur
            valint = x[:,2]; valS = x[:,1]; valB = x[:,0]
        # loi donnée par mu(Sast)*V+kprop*(Sast-S) + kint*valint + kderiv*(dSast/dt-dS/dt)
        mu = muast*valS/(KS+valS+valS**2/KI)
        Qast = muast*Sast/(KS+Sast+Sast**2/KI)*V
        valu = (Qast+kprop*(Sast-valS)+kint*valint+kderiv*k*mu*valB)/(1+kderiv/V*(S0-valS))
    return valu

# Essai de la commande PID
Sast= 1.2  
kprop = 0.05
kint = 0.01
kderiv = 0.4
disturb = 0.1  
culture_cont(Sast,'PID',np.array([kprop,kint,kderiv]),disturb)
