# **********************************************************************
#      TP1: SIMULATION du MODELE de FERMENTATION EN REACTEUR BATCH     *
# **********************************************************************

rm(list=ls(all = TRUE)) 
set.seed(2)
setwd(dir="/home/casenave/Documents/Hubic/Boulot/2011_INRA/Enseignement/Cours/Modelisation_et_Controle_numerique_systemes_dynamiques_agronomie/Version_2018/TP_Cas_etude/Code_R")

## Packages a installer/charger
library(deSolve);

# Modèle de fermentation en réacteur continu
# ------------------------------------------

# fonction second membre du systeme d'equations
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

# parametres du modele
parms <- list(
  k1 = 0.0606, # -
  k2 = 2.17, # -
  mu1max = 1.34, # 1/h
  mu2max = 1.45, # 1/h
  KN = 1.57, # g/L
  KE = 14.1, # g/L
  KS = 0.0154 # g/L
)

# simulation du modèle
# --------------------
# vecteur de temps
times <- seq(from = 0, to = 90, by = 0.1)

# condition initiale
x0  <- c(0.01, 0.425, 0, 200)

# intégration numérique
out <- ode(x0, times, fermenteur, parms)
# on renomme les colonnes de la matrice de sortie
colnames(out) <- c("time","X", "N", "E", "S")

## Chargement des données
data <- read.csv("batch_R.txt")
data[,"N"]=data[,"N"]/1000
data[,"dCO2dt"]=data[,"dCO2dt"]/100

# trace des solutions simulees et des données
par(mfrow=c(2,2))
plot(out[,"time"],out[,"X"],xlab="time",ylab="",type = 'l',col="red")
points(data[,"time"],data[,"X"],xlab="time",col="red")
lines(out[,"time"],out[,"N"]*10,xlab="time",type = 'l',col="green")
points(data[,"time"],data[,"N"]*10,xlab="time",col="green")
legend(50, 5, legend=c("Biomasse X (simulée)","Biomasse X (données)", "Azote N (x10) (simulée)", "Azote N (x10) (données)"),
      col=c("red", "red","green","green"),pch=c(NA,1,NA,1),lty=c(1,NA,1,NA), cex=0.8)

plot(out[,"time"],out[,"S"],xlab="time",ylab="",type = 'l',col="black")
points(data[,"time"],data[,"S"],xlab="time",col="black")
lines(out[,"time"],out[,"E"],xlab="time",type = 'l',col="blue")
points(data[,"time"],data[,"E"],xlab="time",col="blue")
legend(60,200, legend=c("Sucre S (simulée)","Sucre S (données)", "Ethanol E (simulée)","Ethanol E (données)"),
       col=c("black","black", "blue","blue"),pch=c(NA,1,NA,1),lty=c(1,NA,1,NA), cex=0.8)

dCO2 <- as.numeric(parms["mu2max"])*out[,"S"]/(as.numeric(parms["KS"])+out[,"S"])*as.numeric(parms["KE"])/(as.numeric(parms["KE"])+out[,"E"])*out[,"X"]
plot(out[,"time"],dCO2,xlab="time",ylab="",type = 'l',col="magenta")
points(data[,"time"],data[,"dCO2dt"],xlab="time",col="magenta")
legend(60,2, legend=c("dCO2/dt (simulée)", "dCOE/dt (données)"),
       col=c("magenta", "magenta"), pch=c(NA,1),lty=c(1,NA), cex=0.8)

