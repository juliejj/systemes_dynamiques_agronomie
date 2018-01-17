# ********************************************************************
#      TP3: IDENTIFICATION des PARAMETRES DU MODELE DE FERMENTATION  *
# ********************************************************************

rm(list=ls(all = TRUE)) 
set.seed(2)
setwd(dir="/home/casenave/Documents/Hubic/Boulot/2011_INRA/Enseignement/Cours/Modelisation_et_Controle_numerique_systemes_dynamiques_agronomie/Version_2018/TP_Cas_etude/Code_R")

## Packages a installer/charger
library(deSolve)
library(stats4)

## Chargement et affichage des donnees
data <- read.csv("batch_R.txt")
data[,"N"]=data[,"N"]/1000
data[,"dCO2dt"]=data[,"dCO2dt"]/100
Nbdata <- length(data[,"time"]) # nombre de données

## Tracé de N en fonction de B et de S en fonction de E
## -----------------------------------------------------
par(mfrow=c(1,2))
plot(data[,"X"],data[,"N"],xlab="X",ylab = "N",col="red")
plot(data[,"E"],data[,"S"],xlab="E",ylab = "S",col="red")

# Identification des paramètres k1 et k2 du modèle de fermentation en réacteur batch
# ----------------------------------------------------------------------------------

# identification de k1
# --------------------
fit <- lm(-data[,"N"] ~ data[,"X"])
c1_est <- fit$coefficients[1]
k1_est <- fit$coefficients[2]

par(mfrow=c(1,2))
plot(data[,"X"],data[,"N"],xlab="X",ylab = "N",col="red")
lines(data[,"X"],-k1_est*data[,"X"]-c1_est,col="red")
sprintf('estimation de k1 = %e',as.numeric(k1_est))

# identification de k2
# --------------------
fit <- lm(-data[,"S"] ~ data[,"E"])
c2_est <- fit$coefficients[1]
k2_est <- fit$coefficients[2]

plot(data[,"E"],data[,"S"],xlab="time",ylab = "X",col="red")
lines(data[,"E"],-k2_est*data[,"E"]-c2_est,col="red")
sprintf('estimation de k2 = %e',as.numeric(k2_est))

## Modele de fermenteur en réaction batch
# ---------------------------------------
fermenteur <- function(t, x, parms) {
  with(parms, {
    X=x[1];N=x[2];E=x[3];S=x[4]
    # calcul de mu1(S)
    mu1 = mu1max*N/(KN+N)
    # calcul de mu2(E,S)
    mu2 = mu2max*S/(KS+S)*KE/(KE+E)
    
    # second membre de l'équation de B
    dX = mu1*X
    # second membre de l'équation de N
    dN = -k1*mu1*X
    # second membre de l'équation de E
    dE = mu2*X
    # second membre de l'équation de S
    dS = -k2*mu2*X
    
    return(list(c(dX,dN,dE,dS)))
  })
}

# Fonction cout à minimiser
cout <- function(parms,cas_identif,beta,data,k2_est) {
  
  # affichage des paramètres qu'on cherche à identifier
  print(parms)
  
  # récupération des valeurs des paramètres
  k1val <- as.numeric(parms["k1"])
  if (cas_identif == 1){
    k2val <- as.numeric(parms["k2"])}
  else
  {k2val <- k2_est}
  mu1maxval <- as.numeric(parms["mu1max"])
  mu2maxval <- as.numeric(parms["mu2max"])
  KNval <- as.numeric(parms["KN"])
  KEval <- as.numeric(parms["KE"])
  KSval <- as.numeric(parms["KS"])
  X0val <- as.numeric(parms["X0"])
  N0val <- as.numeric(parms["N0"])
  E0val <- as.numeric(parms["E0"])
  S0val <- as.numeric(parms["S0"])
  
  # stockage des valeurs des paramètres nécessaires pour la simulation du modèle
  param <- list(k1 = k1val, k2 = k2val, mu1max = mu1maxval,
                mu2max = mu2maxval, KN = KNval, KE = KEval,
                KS = KSval)
  
  times <- seq(from = 0, to = 90, by = 0.1)
  # condition initiale
  x0  <- c(X0val,N0val,E0val,S0val)
  # intégration numérique du modèle
  out <- ode(x0, data[,"time"], fermenteur, param)
  # on donne des noms aux colonnes
  colnames(out) <- c("time","X", "N", "E","S")
  
  # tracé des sorties simulées et comparaison avec les données
  # par(mfrow=c(2,2))
  # plot(data[,"time"],data[,"X"],ylab="X",xlab="time")
  # lines(out[,"time"],out[,"X"])
  # plot(data[,"time"],data[,"N"],ylab="N",xlab="time")
  # lines(out[,"time"],out[,"N"])
  # plot(data[,"time"],data[,"E"],ylab="E",xlab="time")
  # lines(out[,"time"],out[,"E"])
  # plot(data[,"time"],data[,"S"],ylab="S",xlab="time")
  # lines(out[,"time"],out[,"S"])

  # Calcul de la fonction coût à minimiser donnée par la somme
  # des carrés des écarts
  SumSquares <- as.numeric(beta["betaX"])*sum(((data[1:Nbdata-1,"X"]-out[1:Nbdata-1,"X"])**2
                           +(data[2:Nbdata,"X"]-out[2:Nbdata,"X"])**2)/2
                          *(data[2:Nbdata,"time"]-data[1:Nbdata-1,"time"]))
  SumSquares <- SumSquares+as.numeric(beta["betaN"])*sum(((data[1:Nbdata-1,"N"]-out[1:Nbdata-1,"N"])**2
                           +(data[2:Nbdata,"N"]-out[2:Nbdata,"N"])**2)/2
                          *(data[2:Nbdata,"time"]-data[1:Nbdata-1,"time"]))
  SumSquares <- SumSquares+as.numeric(beta["betaE"])*sum(((data[1:Nbdata-1,"E"]-out[1:Nbdata-1,"E"])**2
                           +(data[2:Nbdata,"E"]-out[2:Nbdata,"E"])**2)/2
                          *(data[2:Nbdata,"time"]-data[1:Nbdata-1,"time"]))
  SumSquares <- SumSquares+as.numeric(beta["betaS"])*sum(((data[1:Nbdata-1,"S"]-out[1:Nbdata-1,"S"])**2
                           +(data[2:Nbdata,"S"]-out[2:Nbdata,"S"])**2)/2
                          *(data[2:Nbdata,"time"]-data[1:Nbdata-1,"time"]))
  
  return(SumSquares)
}

# coefficients de pondération beta
betaXval <- 1/sum(((data[1:Nbdata-1,"X"])**2+(data[2:Nbdata,"X"])**2)/2
               *(data[2:Nbdata,"time"]-data[1:Nbdata-1,"time"]))
betaNval <- 1/sum(((data[1:Nbdata-1,"N"])**2+(data[2:Nbdata,"N"])**2)/2
               *(data[2:Nbdata,"time"]-data[1:Nbdata-1,"time"]))
betaEval <- 1/sum(((data[1:Nbdata-1,"E"])**2+(data[2:Nbdata,"E"])**2)/2
               *(data[2:Nbdata,"time"]-data[1:Nbdata-1,"time"]))
betaSval <- 1/sum(((data[1:Nbdata-1,"S"])**2+(data[2:Nbdata,"S"])**2)/2
               *(data[2:Nbdata,"time"]-data[1:Nbdata-1,"time"]))
beta <- c(betaX=betaXval, betaN=betaNval,betaE=betaEval,betaS=betaSval)

# Cas 1: on identifie tous les paramètres du modèle + les conditions initiales
cas_identif <-1
# Jeu de parametres de départ
StartList <- c(k1 = 0.01, k2 = 2.0, mu1max = 1.2, mu2max = 1.2,
               KN = 1.6, KE = 12., KS = 0.03, X0=data[1,"X"], N0=data[1,"N"],
               E0=data[1,"E"], S0=data[1,"S"])

Estim1 <- optim(StartList, ## Param de départ
               cout, ## fonction à optimiser
               gr=NULL,
               cas_identif,beta,data,k2_est,
               method= "L-BFGS-B", lower=c(0,0,0,0,0,0,0,0,0,0,0)
               ## algo choisi (CG ) gradient conjugué algo possible : 
               #method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
)

# Cas 2: on identifie tous les paramètres du modèle sauf k2 + les conditions initiales
cas_identif <-2
# Jeu de parametres de départ
StartList <- c(k1 = 0.01, mu1max = 1.2, mu2max = 1.2,
               KN = 1.6, KE = 12., KS = 0.03, X0=data[1,"X"], N0=data[1,"N"],
               E0=data[1,"E"], S0=data[1,"S"])

Estim2 <- optim(StartList, ## Param de départ
               cout, ## fonction à optimiser
               gr=NULL,
               cas_identif,beta,data,k2_est,
               method= "L-BFGS-B", lower=c(0,0,0,0,0,0,0,0,0,0)
               ## algo choisi (CG ) gradient conjugué algo possible : 
               #method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
)


# simulation du modèle identifié et comparaison aux données
# ---------------------------------------------------------
# vecteur de temps
times <- seq(from = 0, to = 90, by = 0.1)

# ** modéle identifié 1
# parametres du modele
parms1 <- list(
  k1 = Estim1$par["k1"],
  k2 = Estim1$par["k2"],
  mu1max = Estim1$par["mu1max"],
  mu2max = Estim1$par["mu2max"],
  KN = Estim1$par["KN"],
  KE = Estim1$par["KE"],
  KS = Estim1$par["KS"]
)

# condition initiale
x0  <- c(Estim1$par["X0"],Estim1$par["N0"],Estim1$par["E0"],Estim1$par["S0"])

# intégration numérique
out1 <- ode(x0, times, fermenteur, parms1)
# on renomme les colonnes de la matrice de sortie
colnames(out1) <- c("time","X", "N", "E", "S")

# ** modéle identifié 2
# parametres du modele
parms2 <- list(
  k1 = Estim2$par["k1"],
  k2 = k2_est,
  mu1max = Estim2$par["mu1max"],
  mu2max = Estim2$par["mu2max"],
  KN = Estim2$par["KN"],
  KE = Estim2$par["KE"],
  KS = Estim2$par["KS"]
)

# condition initiale
x0  <- c(Estim2$par["X0"],Estim2$par["N0"],Estim2$par["E0"],Estim2$par["S0"])

# intégration numérique
out2 <- ode(x0, times, fermenteur, parms2)
# on renomme les colonnes de la matrice de sortie
colnames(out2) <- c("time","X", "N", "E", "S")

# trace des solutions simulees et des données
par(mfrow=c(3,2))
plot(out1[,"time"],out1[,"X"],xlab="time",ylab="biomasse",type = 'l',col="magenta")
lines(out2[,"time"],out2[,"X"],xlab="time",type = 'l',col="blue")
points(data[,"time"],data[,"X"],xlab="time",col="red")
legend(50, 5, legend=c("identif 1","identif 2", "données"),
       col=c("magenta", "blue","red"),pch=c(NA,NA,1),lty=c(1,1,NA), cex=0.8)

plot(out1[,"time"],out1[,"N"],xlab="time",ylab="azote",type = 'l',col="magenta")
lines(out2[,"time"],out2[,"N"],xlab="time",type = 'l',col="blue")
points(data[,"time"],data[,"N"],xlab="time",col="red")

plot(out1[,"time"],out1[,"E"],xlab="time",ylab="ethanol",type = 'l',col="magenta")
lines(out2[,"time"],out2[,"E"],xlab="time",type = 'l',col="blue")
points(data[,"time"],data[,"E"],xlab="time",col="red")

plot(out1[,"time"],out1[,"S"],xlab="time",ylab="sucre",type = 'l',col="magenta")
lines(out2[,"time"],out2[,"S"],xlab="time",type = 'l',col="blue")
points(data[,"time"],data[,"S"],xlab="time",col="red")

dCO2_1 <- as.numeric(parms1["mu2max"])*out1[,"S"]/(as.numeric(parms1["KS"])+out1[,"S"])*as.numeric(parms1["KE"])/(as.numeric(parms1["KE"])+out1[,"E"])*out1[,"X"]
dCO2_2 <- as.numeric(parms2["mu2max"])*out2[,"S"]/(as.numeric(parms2["KS"])+out2[,"S"])*as.numeric(parms2["KE"])/(as.numeric(parms2["KE"])+out2[,"E"])*out2[,"X"]

plot(out1[,"time"],dCO2_1,xlab="time",ylab="dCO2/dt",type = 'l',col="magenta")
lines(out2[,"time"],dCO2_2,xlab="time",type = 'l',col="blue")
points(data[,"time"],data[,"dCO2dt"],xlab="time",col="red")

