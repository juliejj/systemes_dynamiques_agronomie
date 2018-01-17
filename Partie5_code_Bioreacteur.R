# ************ CODE POUR LA SIMULATION ET LE CONTROLE DU BIOREACTEUR **********

rm(list=ls(all = TRUE)) 
set.seed(2)
setwd(dir="/home/casenave/Documents/Hubic/Boulot/2011_INRA/Enseignement/Cours/Modelisation_et_Controle_numerique_systemes_dynamiques_agronomie/Version_2018/TP_Cas_etude/Code_R")

## Packages a installer/charger
library(deSolve);

# modèle de réacteur continu
reacteur <- function(t, x, parms) {
  with(parms, {

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
    B = x[1] # biomasse
    S = x[2] # substrat

    # Calcul de la commande: dans fonction_u la loi de commande est calculée, disturb permet de prendre
    # en compte d'éventuelles perturbation sur la valeur réellement appliquée de la commande
    Q = fonction_u(t,x,Sast,k,muast,KI,KS,V,S0,control_type,coeffcontrol)*(1+disturb)

    # taux de croissance (fonction de Haldane)
    mu = muast*S/(KS+S+S**2/KI)
    # second membre de l'équation en B
    dB = mu*B-Q/V*B
    # second membre de l'équatio en S
    dS = -k*mu*B+Q/V*(S0-S)

    # initialisation de dx, second membre du modèle correspondant à la dérivée des variables d'état dB/dt et dS/dt
    if ((control_type =='PI')||(control_type == 'PID')){
        # si il y a un terme intégral dans la commande, on rajoute une équation pour
        # calculer l'intégrale Sast-S au cours du temps
        # second membre de l'equation qui calcule l'intégrale de Sast-S
        dSint = Sast-S
        dx = list(c(dB,dS,dSint))}
    else{ 
        # si il n'y a pas de terme intégral dans la commande, on ne rajoute pas d'équation
        dx = list(c(dB,dS))}

    return(dx)
  })
}


# ** fonction qui simule le modèle avec la loi de commande demandée et qui trace ensuite la solution
culture_cont <- function(Sastval,control_typeval,coeffcontrolval,disturbval){
    
    # vecteur de temps pour la simulation
    times <- seq(from = 0, to = 150, by = 0.075)

    # paramètres du modèle
    kval = 0.6
    muastval = 2.3
    KSval = 10
    KIval = 0.1
    Qinval = 0.01
    Vval = 0.5
    S0val=3.2
    
    parms <- list(
    k=kval,
    muast=muastval,
    KS=KSval,
    KI=KIval,
    Qin=Qinval,
    V=Vval,
    S0=S0val,
    Sast=Sastval,
    control_type=control_typeval,
    coeffcontrol=coeffcontrolval,
    disturb=disturbval
    ) 
    
    # conditions initiales du modèle (valeurs initiales de la biomasse et de la concentration en sucre)
    B0 = 9; S0 = 3.2; 
    
    # si on utilise un terme intégrale, il faudra rajouter une équation dans le modèle et du coup rajouter
    # la condition initiale correspondante qui est égale à 0 (voir plus loin dans le paragraphe sur le terme
    # intégral)
    if ((control_typeval =='PI')||(control_typeval == 'PID')){
        x0  <- c(B0, S0, 0)}
    else{
        x0  <- c(B0, S0)}
        
    # intégration numérique de l'EDO
    out <- ode(x0, times, reacteur, parms)
    # on renomme les colonnes de la matrice de sortie
    if ((control_typeval =='PI')||(control_typeval == 'PID')){
        colnames(out) <- c("time","B", "S","Sint")}
    else{
        colnames(out) <- c("time","B", "S")}
    
    # re-calcul de la commande appliquée
    u = fonction_u(times,out,Sastval,kval,muastval,KIval,KSval,Vval,S0,control_typeval,coeffcontrolval)

    # tracé des solutions

    par(mfrow=c(1,2))
    plot(out[,"time"],out[,"B"],xlab="time",ylab="",type = 'l',col="red",ylim=c(0, 10))
    lines(out[,"time"],out[,"S"],xlab="time",type = 'l',col="green")
    lines(c(0,150),c(Sastval,Sastval),type='l',lty=2,col='green')
    legend(60,10, legend=c("Biomasse", "Sucre","Sast"),
           col=c("red", "green","green"), pch=c(NA,NA,NA),lty=c(1,1,2), cex=0.6)
    
    plot(out[,"time"],u,xlab="time",ylab="",type = 'l',col="red")
    legend(60,0.08, legend=c("Q"),
           col=c("red"), pch=c(NA),lty=c(1), cex=0.6)
    
}
    
# Loi de commande
fonction_u <- function(t,x,Sast,k,muast,KI,KS,V,S0,control_type,coeffcontrol){
    if(control_type == 'BO'){ # BOUCLE OUVERTE
        # loi donnée par Qast=mu(Sast)*V
        Qast = muast*Sast/(KS+Sast+Sast**2/KI)*V
        valu = matrix(1,ncol(t),1)*Qast
    }
    else if(control_type == 'P'){ # BOUCLE FERMEE action PROPORTIONNELLE
        # récupération des paramètres de la loi de commande
        kprop = coeffcontrol
        # et de la valeur de S
        if(length(t)==1){ # cas où t est scalaire 
          valS = x[2]
        }
        else{ # cas où t est un vecteur
          valS = x[,"S"]
        }
        # loi donnée par mu(Sast)*V+kprop*(Sast-S)
        valu = muast*Sast/(KS+Sast+Sast**2/KI)*V+kprop*(Sast-valS)
    }
    else if(control_type == 'PI'){ # BOUCLE FERMEE action PROPORTIONNELLE INTEGRALE
        # récupération des paramètres de la loi de commande
        kprop = coeffcontrol[1]
        kint = coeffcontrol[2]
        # et de la valeur de valint, qui est l'intégrale entre 0 et t de Sast-S
        if(length(t)==1){ # cas où t est scalaire 
            valint = x[3]; valS = x[2]
        }
        else{ # cas où t est un vecteur
            valint = x[,"Sint"]; valS = x[,"S"]
        }
        # loi donnée par mu(Sast)*V+kprop*(Sast-S) + kint*valint
        # où valint est l'intégrale entre 0 et t de Sast-S
        Qast = muast*Sast/(KS+Sast+Sast**2/KI)*V
        valu = Qast+kprop*(Sast-valS)+kint*valint
    } 
    else if(control_type == 'PID'){ # BOUCLE FERMEE action PROPORTIONNELLE INTEGRALE DERIVEE
        # récupération des paramètres de la loi de commande
        kprop = coeffcontrol[1]
        kint = coeffcontrol[2]
        kderiv = coeffcontrol[3]
        # et des valeurs de valint (intégrale de Sast-S), de S et de B
        if(length(t)==1){ # cas où t est scalaire 
            valint = x[3]; valS = x[2]; valB=x[1]
        }
        else{ # cas où t est un vecteur
            valint = x[,"Sint"]; valS = x[,"S"]; valB = x[,"B"]
        }
        # loi donnée par mu(Sast)*V+kprop*(Sast-S) + kint*valint + kderiv*(dSast/dt-dS/dt)
        mu = muast*valS/(KS+valS+valS**2/KI)
        Qast = muast*Sast/(KS+Sast+Sast**2/KI)*V
        valu = (Qast+kprop*(Sast-valS)+kint*valint+kderiv*k*mu*valB)/(1+kderiv/V*(S0-valS))
    }
    else{
        valu=0
    }
    return(valu)
}

# Essai de la commande PID
Sast= 1.2  
kprop = 0.05
kint = 0.01
kderiv = 0.4
disturb = 0.1  
culture_cont(Sast,'PID',c(kprop,kint,kderiv),disturb)

