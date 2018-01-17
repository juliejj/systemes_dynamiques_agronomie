# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 21:49:11 2018

@author: casenave
"""

# ************ CODE POUR LA SIMULATION ET LE CONTROLE DU RESERVOIR **********

import numpy as np
from matplotlib import pyplot as plt
import scipy.integrate as scint
from matplotlib import patches as pat
plt.close('all') # ferme toutes les figures

# t est le temps courant
# T est la durée pendant laquelle on va pomper dans la nappe fréatique pour remplir le réservoir
# dmax est la valeur maximale du débit
# base est l'aire de la base carrée du réservoir
# cote est la longueur du côté de la base
# precip et conso sont des variables qui contiennent des données de précipitations et de consommation par les
#        les agriculteurs: cest variables ne sont pour l'instant pas utilisées
# x est la variable d'état du système, c'est à dire ici la hauteur d'eau dans le réservoir
# Vc est le volume que l'on veut atteindre dans le réservoir

# modèle du réservoir
def reservoir(x,t,T,dmax,base,precip,conso,Vc):
    # le modèle du réservoir est donné par dx/dt = u/base
    if np.fix(t)>len(precip)-1: 
        dx = fonction_u(t,T,dmax,base,Vc,precip,conso,x)/base
    else : 
        dx = (fonction_u(t,T,dmax,base,Vc,precip,conso,x[0])+precip[int(np.fix(t))]+conso[int(np.fix(t))])/base
    return dx

def remplissage(pluie,arrosage):
    # pluie = 1 si on prend en compte les précipitations, 0 sinon
    # arrosage = 1 si on prend en compte la consommation des agriculteurs, 0 sinon
    
    # paramètres du modèles
    dmax  = 2.; cote  = 3.; base  = cote**2; Vmax  = 100; Vc = 90; T = Vc/dmax 
    
    # Prise en compte des précipitations pendant 1 an
    if pluie == 0: # si on ne considère pas les précipitations
        precip = np.zeros(365)
        tmax = min(T*1.5,365)
        nblinefig = 2
    else: # si on considère les précipitations, on va charger un fichier de données
        precip0 = np.loadtxt('data/precipitations.txt') # données en mm par jour (par m² de terrain)
        precip = 1e-3*precip0*base # conversion des données en m³ par jour
    
    # Prise en compte des arrosages
    conso = np.zeros(len(precip))
    if arrosage == 1: 
        # ajout de periodes dans l'année pendant lesquelles les agriculteurs utilisent de l'eau du réservoir
        conso[10:15] = -1
        conso[40:45] = -1
        conso[70:75] = -1
        conso[100:105] = -1
        conso[130:135] = -1
    
    # integration numerique de l'EDO
    if pluie + arrosage >0: 
        tmax=365 # temps maximal de simulation
        nblinefig = 4 # nombre de figure pour l'affichage
    temps = np.linspace(0,tmax,2000) # vecteur temps
    x0 = 0 # condition initiale
    h = scint.odeint(reservoir,x0,temps,args=(T,dmax,base,precip,conso,Vc)) # integration de l'équation
    u = fonction_u(temps,T,dmax,base,Vc,precip,conso,h) # re-calcul de la commande qui a été appliquée

    # tracé des solutions
    plt.figure(figsize = (10, 7))
    plt.subplots_adjust(hspace=0.4,wspace=0.4)
    plt.subplot2grid((nblinefig,3),(0,0),rowspan=nblinefig)
    axes=plt.gca()
    axes.add_artist(pat.Rectangle((0, 0), cote, h[-1], color = 'blue'))
    plt.ylim([0,Vmax/base])
    plt.xlim([0,cote])
    plt.subplot2grid((nblinefig,3),(0,1),colspan=2)
    plt.plot(temps, h, color="red", linewidth="1")
    plt.plot(np.array([0,tmax]), np.array([Vc/base,Vc/base]), color="black", linewidth="1")
    plt.xlim([0,tmax])
    plt.ylim([0,Vmax/base])
    plt.title("Hauteur d'eau (m)")
    plt.subplot2grid((nblinefig,3),(1,1),colspan=2)
    plt.xlim([0,tmax])
    plt.plot(temps, u, color="red", linewidth="1")
    plt.ylim([-0.2,dmax*1.1])
    plt.title("débit d'entrée ($m^3$ par jour)")
    if pluie + arrosage == 0:
        plt.xlabel("Temps (jour)")
        plt.xlim([0,tmax])
    else:
        plt.subplot2grid((nblinefig,3),(2,1),colspan=2)
        plt.plot(np.arange(365),precip0,color="red", linewidth="1")
        plt.xlim([0,tmax])
        plt.ylim([0,60])
        plt.title("Précipitations (mm) par jour")
        plt.subplot2grid((nblinefig,3),(3,1),colspan=2)
        plt.plot(np.arange(365),conso,color="red", linewidth="1")
        plt.xlim([0,tmax])
        plt.ylim([-1.2,0])
        plt.xlabel("Temps (jour)")
        plt.title("Consommation en $m^3$ par jour")
    plt.show()

# Loi de commande boucle ouverte
def fonction_u(t,T,dmax,base,Vc,precip,conso,x):
   
    # Si t<T, la commande vaut dmax, sinon elle vaut 0
    if type(t)==float: valu = (t<=T)*dmax # cas où l'entrée t de la fonction est un scalaire
    else: # cas où l'entrée t de la fonction est un vecteur
        valu = np.zeros(len(t))
        valu[t<=T] = dmax  
    return valu
    
# Remplissage d'un réservoir par une loi de commande boucle ouverte, 
# sans prise en compte ni des précipitations ni de l'utilisation de l'eau par les agriculteurs
remplissage(0,0)

# Remplissage d'un réservoir par une loi de commande boucle ouverte, 
# AVEC prise en compte des précipitations et de l'utilisation de l'eau par les agriculteurs
remplissage(1,1)

# Remplissage d'un réservoir par une loi de commande boucle fermée
# AVEC prise en compte ni des précipitations ni de l'utilisation de l'eau par les agriculteurs
def fonction_u(t,T,dmax,base,Vc,precip,conso,x):
    coeffalpha=0.1
    valu = base*coeffalpha*(Vc/base-x)
    return valu

remplissage(1,1)