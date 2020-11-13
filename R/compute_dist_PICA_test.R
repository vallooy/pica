# Charger toutes les fonctions d abord :
Calcul_fourmis_distance <- function() {
  parcelle_ID <- "Estagnol"
  annee <- 2017
  selected <- c(76,66,39,44,48,79)
  
  parcelle <- read.csv("https://cloud.opencpu.org/ocpu/apps/vallooy/pica/www/Parcelle_simulated_Estagnol_2017.csv", sep = ",", dec = ".",header = T)  
  #La liste des points de la parcelle
  
  liste_sites <- read.csv("https://cloud.opencpu.org/ocpu/apps/vallooy/pica/www/data_to_python_matrix_estagnol.csv", sep = ";", dec = ",",header = T)  
  # La liste des individus associés aux ligne de la matrice de distance (un individu apparait deux fois, une fois pour chaque interang
    
  dist_mat <- read.csv("https://cloud.opencpu.org/ocpu/apps/vallooy/pica/www/test_dist_estagnol2.csv", sep = ";", dec = ",",header = T)  
  # La matrice de distance, assez long à charger
  
  #entree_file <- paste0(dossier,"Point_entree_",parcelle_ID,".txt")
  #entree <- read.table(entree_file, sep = ",", dec = ".",header = T)
  # Le point d'entrée dans la parcelle
  
  ID_start <- NA
  # Les bords de rang utilises comme point d'entree dans la parcelle
  if (parcelle_ID == "Arnel") { ID_start <- 2310 }
  if (parcelle_ID == "Estagnol") { ID_start <- 1021 }
  if (parcelle_ID == "Larzat") { ID_start <- 3946 }
  
  N <- length(selected) #nombre de site de mesure
  
  
  parcelle_selected <- parcelle[parcelle[,"Site"] %in% c(selected),]
  #sites de mesure seletionnes
  
  parcelle_selected <- Sites_duplicate(parcelle_selected = parcelle_selected,
                                       liste_sites = liste_sites)
  #On les duplique en considerant qu'il peuvent etre associes a deux inerangs differents
  ID_selected <- c(ID_start)
  
  #on fait la liste des identifiants dans la matrice de distance en commencant par le point de depart
  for (i in 1:length(parcelle_selected$Site)) {
    ID <- liste_sites$ID[(liste_sites$X == parcelle_selected$X[i]) & 
                           (liste_sites$Y == parcelle_selected$Y[i]) & 
                           (liste_sites$Rang == parcelle_selected$Rang[i])]
    ID_selected <- c(ID_selected, ID)
  }
  # On retrouve le numero de chaque site dans la matrice de distance
  
  Label2 <- c(LETTERS[1],paste(LETTERS[2:(N + 1)],1,sep = ""),paste(LETTERS[2:(N + 1)],2,sep = ""))
  
  #On fait des groupes pour s'assurer de bien selectionner une seule fois chaque site de mesure
  # ( ils corespondent a deux sites chacun, il serait possible de les compter deux fois sinon)
  groupe_selected <- c(LETTERS[1],LETTERS[2:(N + 1)],LETTERS[2:(N + 1)])
  
  # la matrice de distance
  mat_selected <- Real_dist_matrix(ID_selected = ID_selected,
                                   dist_mat = dist_mat)
  
  names(groupe_selected) <- rownames(mat_selected) <- colnames(mat_selected) <- Label2
  
  mat_selected <- as.dist(mat_selected)
  
  #On appelle l'algo des fourmis
  solution <- fColonie(MatDist = mat_selected, Groupe = groupe_selected, Nfourmi = 30, 
                   Niter = 20, Alpha = 1.5, Beta = 1, Gamma = 0, Taux = 0.1, 
                   Renforcement = FALSE, kPermute = 2, Qtot = 1, 
                   Rho = .5, Boucle = TRUE, txQmin = 0.0001, txQmax = Inf,
                   Plot = FALSE, Coord = NULL, xlim = c(0,1), ylim = c(0,1), Lwd = 5)
  
  distance <- solution$Best$Longueur
  
  #On rajoute la distance entre le point d'entrée et le bord de rang
  #distance <- distance + 2*sqrt( (entree$X - liste_sites$X[liste_sites$ID == ID_start])^2 + 
  #                                 (entree$Y - liste_sites$Y[liste_sites$ID == ID_start])^2 )

  chemin <- solution$Best$Chemin[1,]
  chemin <- chemin[ 1:(N + 1)]
  
  while (chemin[1] != "A") { # on reordonne pour que ça parte de A, le point de depart
    chemin <- c(chemin[2:(N + 1)],chemin[1])
  }
  
  
  parcours <- liste_sites[liste_sites$ID == ID_selected[names(groupe_selected) == "A"],]
  #on cree un tableau qui contient les coordonnees et rangs de tous les points dans l ordre de passage
  
  parcours_sites <- c(0)
  
  for (k in 2:length(chemin)) {
    parcours <- rbind(parcours,
                      liste_sites[liste_sites$ID == ID_selected[names(groupe_selected) == chemin[k]],])
    
    ID_site <- parcelle_selected$Site[(parcelle_selected$X == parcours$X[k]) & 
                                        (parcelle_selected$Y == parcours$Y[k]) & 
                                        (parcelle_selected$Rang == parcours$Rang[k])]
    
    parcours_sites <- c(parcours_sites, ID_site)
  }
  
  # Pour garder les memes identifiants qu'au debut
  parcours$Site <- parcours_sites
  colnames(parcours) <- c("X","Y","Interang","Type","ID_matrix","ID_site")
  
  sortie <- list(distance = distance, parcours = parcours)
  return(sortie)
  
  
}

Real_dist_matrix <- function(ID_selected,dist_mat){
  
  size <- length(ID_selected)
  
  mat_dist <- matrix(data = 0, nrow = size, ncol = size)
  for (i in 2:size) {
    for (j in 1:(i - 1)) {
      
      ind1 <- ID_selected[i]

      ind2 <- ID_selected[j]

      mat_dist[i,j] <- mat_dist[j,i] <- dist_mat[ind1,ind2]
    }
  }
  return(mat_dist)
}


Sites_duplicate <- function(parcelle_selected, liste_sites) {
  parcelle_selected2 <- parcelle_selected
  parcelle_selected2$Rang <- parcelle_selected2$Rang + 1
  
  parcelle_selected2$Rang[parcelle_selected2$Rang > max(liste_sites$Rang)] <- max(liste_sites$Rang)
  
  parcelle_selected <- rbind(parcelle_selected,parcelle_selected2)
  
  return(parcelle_selected)
}



# Tout ce qui est en dessous est un ensemble de fonctions crees avec Giles Le moguedec 
# au cours de la these qui me servent a calculer la distance via la methode des fourmis


# Probleme du voyageur de commerce
# Resolution par colonie de fourmis
# On ne garde que les fonctions

# Possibilite de creer des groupes de points pour chacun de ces points
# Il s'agira alors de passer une et une seule fois par chacun de ces groupes

# Creation des points par groupe 
# Le point original est au centre de gravite des points de son groupe
# Coord : matrice des coordonnees des points originaux
# NGroupe : vecteur des effectifs associes au points originaux
# Rayon : Distance des nouveaux points a leur centre de gravite
# Ces nouveaux points sont regulierement repartis autour du centre de gravite
createGroupe<-function(Coord=NULL,NGroupe=NULL,Rayon=.01){
  Sortie<-NULL
  if (!is.null(Coord)){
    if (is.null(NGroupe)){ CoordFin<-Coord}
    else {
      NGroupe<-rep(NGroupe,length=nrow(Coord))
      Indice<-1:nrow(Coord)
      ListSortie<-sapply(X=Indice,FUN=function(x){
        Label<-rownames(Coord)[x]
        if (is.null(Label)){ Label<-paste(1:nrow(Coord),".",sep = "")}
        n<-NGroupe[x]
        if (n<=1){
          Sortie<-matrix(data=Coord[x,],nrow=1)
        }
        else {
          Label<-paste(Label,1:n,sep="")
          Theta<- (1:n)*2*pi/n
          Sortie<-Coord[rep(x,n),]+matrix(data=Rayon*c(cos(Theta),sin(Theta)),ncol=2)
        }
        rownames(Sortie)<-Label
        return(Sortie)
      })
      CoordFin<-NULL
      for (k in 1:length(ListSortie)){ CoordFin<-rbind(CoordFin,ListSortie[[k]])}
      names(dimnames(CoordFin))<-names(dimnames(Coord))
    }
    Groupe<-rep(names(NGroupe),times=NGroupe)
    names(Groupe)<-rownames(CoordFin)
    Sortie<-list(Coord=CoordFin,Groupe=Groupe)
  }
  return(Sortie)
}

# Calcul des distances
# Distorsion des distances pour introduire eventuellement des asymetries
# On decide arbitrairement de ne distordre que selon l'axe Oy
# Si Sym = TRUE, la fonction retourne un objet de type "dist"
# Sinon, c'est une matrice asymetrique
fDistorsion<-function(Coord, Sym=FALSE, ylim=c(0,1), method = "euclidean", diag = FALSE, upper = FALSE, p = 2){
  Distance<-dist(x=Coord, method = method, diag = diag, upper = upper, p = p)
  if (!Sym){
    Distance<-as.matrix(Distance)
    
    Delta<-Distance
    for ( k in colnames(Delta)){
      Delta[,k]<-Coord[rownames(Delta),"y"] - Coord[k,"y"]
    }
    Delta<-1+Delta/(ylim[2]-ylim[1])
    Delta<-ifelse(Delta<=0,0,ifelse(Delta>=1,1,Delta))
    Distance<-Distance*Delta
    
  }
  return(Distance)
}


# Nombre de chemins possibles
# Si Groupe != NULL, la valeur de n prise en compte
# est le nombre d'elements uniques de Groupe
# Sinon on prend le n declare
# n : Nombre de points du chemin
# Boucle : le chemin fait-il une boucle ?
# Sym : la matrice de distances est-elle symetrique ?
# Groupe : vecteur nomme pour la composition des groupes
fnbChemin<-function(n=1,Boucle=FALSE,Sym=TRUE,Groupe=NULL){
  if (!is.null(Groupe)){
    n<-length(unique(Groupe))
    Prod<-prod(as.vector(table(Groupe)))
  }
  else { Prod<-1}
  Boucle<-rep(Boucle,length(n))
  Sym<-rep(Sym,length(n))
  n1<-ifelse(Boucle,n-1,n)
  nb<-ifelse(n>0,factorial(n1),0)
  nb<-ifelse(Sym & n1>1, nb/2,nb)
  nb<-Prod*nb
  return(nb)
}

# Nombre de combinaisons par permutations de k sommets successif d'un chemin de longueur n
fnbpermute<-function(n=1,k=1,Boucle=FALSE){
  n<-pmax(0,ifelse(Boucle,n-1,n))
  Nb<-ifelse(k<=1,1,ifelse(k>=n,factorial(n),1+(factorial(k)-1)*(n-k+1)))
  return(Nb)
}

# Cheminement pour une fourmi
# Start : point de depart si impose
# MatProb : Matrice de probabilite du choix d'aller d'un point vers un autre
# Boucle : le chemin fait-il une boucle ?
# Groupe : vecteur nomme pour la composition des groupes
# Le resultat est un vecteur donnant la succession des points d'un chemin
fCheminement<-function(Start=NULL,MatProb=NULL,Boucle=FALSE,Groupe=NULL){
  Chemin<-NULL
  if (!is.null(MatProb)){
    if (is.null(Groupe)){ 
      Groupe<-rownames(MatProb)
      if (is.null(Groupe)){ Groupe<-1:nrow(Matrice)}
      names(Groupe)<-Groupe
    }
    Label<-names(Groupe)
    Prob<-diag(MatProb)
    if (all(Prob<=0)) { Prob<-NULL }
    if (is.null(Start)){ Start=sample(x=Label,size=1,prob=Prob) }
    Chemin<-Start
    Courant<-Start
    Reste<-Label[!(Groupe %in% Groupe[Chemin])]
    while (length(Reste)>0) {
      Prob<-MatProb[Courant,Reste]
      if (sum(abs(Prob))==0) { Prob<-NULL }
      Prochain<-sample(x=Reste,prob=Prob,size=1)
      Chemin<-c(Chemin,Prochain)
      Courant<-Prochain
      Reste<-Label[!(Groupe %in% Groupe[Chemin])]
    }
    if (Boucle) { Chemin<-c(Chemin,Start) }
  }
  return(Chemin) 
}

# Longueur d'un chemin
# Chemin : vecteur des noms des points du chemin
# MatDist : matrice des distances
lchemin<-function(Chemin=NULL,MatDist=NULL){
  Longueur<-NULL
  if (!is.null(Chemin) & !is.null(MatDist)){
    MatDist<-as.matrix(MatDist)
    Longueur<-0
    if (length(Chemin)>0){
      for (k in 1:(length(Chemin)-1)){ Longueur<-Longueur+MatDist[Chemin[k],Chemin[k+1]] }
    }
  }  
  return(Longueur)
}

# Representation : intensite de la liaison entre chaque couple de points
# Coord : Coordonnees des points a representer
# Groupe : affectation de chacun de ces points a un groupe
# Si non null, on represente les points autour de leur centre de gravite
# et alors seul ce dernier est nomme
# MatPher : Matrice des Pheromones sur chaque arrete
# composante diagonale : intensite pour le point de depart (si le chemin n'est pas une boucle)
# Si diagonale est nulle, le chemin fait une boucle
# QPlafond : quantite maximale possible de pheromone par arrete (pour etalonner les couleurs)
# Lref <- longueur de reference pour standardiser QPlafond
# de maniere a etre homogene avec la quantite utilisee dans le programme :
# une quantite de pheromones par unite de longueur.
# Sym : la matrice des distances est-elle symetrique ?
# Titre : titre du graphique
# xlim, ylim : limites des axes pour la representation
# Lwd : epaisseur du trait entre deux points
# La sortie est une matrice des intensites sur une echelle de 0 a 1
fPlot<-function(Coord=NULL,Groupe=NULL,MatPher=NULL,QPlafond=0,Lref=1,Sym=TRUE,Titre=NULL,xlim=c(0,1),ylim=c(0,1),Lwd=5){
  Sortie<-NULL
  if (!is.null(Coord)){
    if (is.null(Groupe)){ Coord0<-Coord }
    else {
      Coord0<-aggregate(x=Coord,by=list(Groupe=Groupe),FUN=mean)
      rownames(Coord0)<-Coord0$Groupe
      Coord0<-as.matrix(Coord0[,-1])
    }
    if (Lref<=0){ Lref<-1}
    QPlafond<-QPlafond/Lref
    plot(x=Coord0,asp=1,xlim=xlim,ylim=ylim,pch=19,main=Titre)
    if (!is.null(Groupe)){ points(x=Coord,pch=21,col="grey")}
    if (!is.null(MatPher) & QPlafond>0){
      if (nrow(MatPher)>1){
        if (Sym) { MatPher<-MatPher+t(MatPher) }
        MatPher<-ifelse(MatPher>QPlafond,1,MatPher/QPlafond)
        if (any(diag(MatPher)>0)){points(x=Coord,pch=19,cex=3,col=rgb(red=1,blue=1,green=0,alpha=diag(MatPher)))}
        for (k1 in 1:(nrow(MatPher)-1)){
          for (k2 in (k1+1):nrow(MatPher)){
            if (Sym){ if (MatPher[k1,k2]>0){ lines(x=Coord[c(k1,k2),],col=rgb(red=1,blue=1,green=0,alpha=MatPher[k1,k2]),lwd=Lwd) } }
            else {
              if ( MatPher[k1,k2]==MatPher[k2,k1] ){
                if (MatPher[k1,k2]>0){ arrows(x0=Coord[k1,"x"],y0=Coord[k1,"y"],
                                              x1=Coord[k2,"x"],y1=Coord[k2,"y"],
                                              col=rgb(red=1,blue=1,green=0,alpha=MatPher[k1,k2]),lwd=Lwd,code=3) }
              }
              else {
                Lwd2<-ceiling(Lwd/2)
                if (MatPher[k1,k2]>MatPher[k2,k1]){
                  if (MatPher[k2,k1]>0){ arrows(x0=Coord[k1,"x"],y0=Coord[k1,"y"],
                                                x1=Coord[k2,"x"],y1=Coord[k2,"y"],
                                                col=rgb(red=1,blue=0,green=0,alpha=MatPher[k2,k1]),lwd=Lwd,code=1) }
                  arrows(x0=Coord[k1,"x"],y0=Coord[k1,"y"],
                         x1=Coord[k2,"x"],y1=Coord[k2,"y"],
                         col=rgb(red=0,blue=1,green=0,alpha=MatPher[k1,k2]),lwd=Lwd2,code=2)
                }
                else {
                  if (MatPher[k1,k2]>0){  arrows(x0=Coord[k1,"x"],y0=Coord[k1,"y"],
                                                 x1=Coord[k2,"x"],y1=Coord[k2,"y"],
                                                 col=rgb(red=1,blue=0,green=0,alpha=MatPher[k1,k2]),lwd=Lwd,code=2) }
                  arrows(x0=Coord[k1,"x"],y0=Coord[k1,"y"],
                         x1=Coord[k2,"x"],y1=Coord[k2,"y"],
                         col=rgb(red=0,blue=1,green=0,alpha=MatPher[k2,k1]),lwd=Lwd2,code=1)
                } 
              } 	
            }
          }
        }
      }
      Sortie<-MatPher
    }
    points(x=Coord0,pch=19)
    text(x=Coord0,labels=rownames(Coord0),asp=1,adj=c(1,1))
  }
  return(Sortie)
}


# Calcul de la quantite maximale de pheromone possible sur une arrete
# (Sans tenir compte du renforcement eventuel)
# Nfourmi : nombre de fourmis pour la recherche
# Qtot : Quantite de pheromone deposee par unite de longueur sur une arrete
# Rho : taux de conservation des pheromones (= 1- taux d'evaporation) par iteration
# Niter : Nombre d'iteration maximal
# Taux : Taux de fourmis selectionnees pour le depot des pheromones
fPlafond<-function(Nfourmi=1,Qtot=1,Rho=0,Niter=0,Taux=1){
  Rho<-pmin(1,pmax(0,Rho))
  Niter<-Niter+1
  Sortie<- ceiling(Nfourmi*Taux)*Qtot*ifelse(Rho<1,(1-Rho^Niter)/(1-Rho),Niter*Rho)
  return(Sortie)
}

# Filtrage des solutions equivalentes en cas de symetrie :
# Parcourir un chemin dans un sens ou dans l'autre est equivalent,
# On ne garde qu'un seul sens si les deux figurent dans la liste
# Tab chemin : matrice des chemins : une ligne par chemin
# Le resultat est une matrice de chemins
filtreSym<-function(TabChemin){
  if (!is.null(TabChemin)){
    nBest<-nrow(TabChemin)
    if (nBest>1){
      Keep<-rep(TRUE,nBest)
      for (k1 in 1:(nBest-1)){
        if (Keep[k1]) {
          for (k2 in (k1+1):nBest) { if (all(TabChemin[k1,]==rev(TabChemin[k2,]))) { Keep[k2]<-FALSE }}
        }
      }
      TabChemin<-TabChemin[Keep,]
      if (is.vector(TabChemin)) { TabChemin<-matrix(data=TabChemin,nrow=1) }
    }
  }
  return(TabChemin)
}

# Filtrage d'un cycle
# Si le chemin est une boucle, peu importe le point de depart
# Cette fonction elimine les chemins superflus dans ce cas
# Matrice : matrice des chemins, une ligne par chemin
# Le resultat est une matrice de chemins
filtreCycle<-function(Matrice){
  fextrait<-function(n,Chemin){
    ncol<-length(Chemin)
    if (length(n)==0){n<-0}
    if (n>=1 & n<=ncol){
      if (n>1){ Indice<-c(n:ncol,1:(n-1))}
      else { Indice<-1:ncol}
      Chemin<-Chemin[Indice]
    }
    return(Chemin)
  }
  Sortie<-Matrice
  if (is.matrix(Matrice)){
    if (all(Matrice[,1]==Matrice[,ncol(Matrice)])){
      if (ncol(Matrice)>1){
        Matrice<-Matrice[,-ncol(Matrice)]
        if (is.vector(Matrice)){ Matrice<-matrix(data = Matrice, nrow = 1)}
        Ref<-apply(X=Matrice,MARGIN=1,FUN=function(x){return(which(x==Matrice[1,1]))})
        Sortie<-t(mapply(n=Ref,Chemin=as.list(as.data.frame(t(Matrice))),FUN = fextrait))
        if (is.vector(Sortie)){ Sortie<-matrix(data = Sortie,nrow=1)}
        Sortie<-cbind(Sortie,Sortie[,1])
        Sortie<-unique(Sortie)
      }
    }
  }
  return(Sortie)
}

# Algorithme de recherche par colonie de fourmis
# MatDist : Matrice des distances entre points
# Groupe : repartition eventuelle de ces points en groupe
# Start : Point de depart de tous les chemins testes (NULL par defaut). Sans interet si Boucle=TRUE.
# Dans ce cas on cherche un chemin qui passe une et une seule fois dans chaque groupe
# MatPher : matrice initiale des pheromones. Par defaut on part d'une matrice nulle.
# Nfourmi : Nombre de fourmis se deplacant par iteration
# Niter : nombre d'iterations de l'algorithme
# Alpha, Beta, Gamma : parametres pour la conversion d'une matrice de pheromones en matrice de probabilites
# Alpha : puissance associee a la matrice des pheromones 
# Beta : puissance associee a la matrice de visibilite/desirabilite (= matrice des inverses des distances)
# Gamma : Constante positive ou nulle visant a garantir une probabilite minimale pour chaque arrete.
# Elle sert a eviter la convergence prematuree
# Qtot : quantite de pheromone deposee par unite de longueur sur chaque arrete
# Sans effet dans la version actuelle mais figure dans toutes les versions de l'algorithme que j'ai trouvees
# Je n'ai pas compris son interet, je l'ai conserve pour des evolutions ulterieures.
# Rho : taux de conservation des pheromones (= 1- taux d'evaporation) par iteration
# Boucle : le chemin fait-il une boucle ?
# Le programme calcule QPlafond, la quantite maximale de pheromones qu'il est theoriquement
# possible de deposer sur une arrete (sans tenir compte du renforcement eventuel)
# Il associe aux deux parametres suivants qui servent a limiter les risques de convergence prematuree
# txMin : fraction de QPlafond servant de valeur minimale garantie de pheromones sur chaque arrete ;
# txMax : fraction de QPlafond servant de valeur maximale garantie de pheromones sur chaque arrete 
# Taux : Taux de fourmis selectionnees pour le depot des pheromones a chaque iteration
# Renforcement : si TRUE, a chaque iteration une fourmi supplementaire depose des pheromones sur le meilleur chemin trouve jusqu'ici
# kPermute : pour chaque chemin parcouru, on teste les permutations de kPermutes sommets contigus
# Si kPermute < n (=longueur du chemin, longueur-1 si boucle), il y a 1+ (k!-1)*(n-k+1) combinaisons a tester
# => Ne pas depasser kPermute = 3 sous peine d'exploser les temps de calcul
# Ce qui suit ne sert qu'aux representations graphiques (voir fonction fPlot)
# Plot : si TRUE : a chaque iteration on produit un graphique de la quantite de pheromones deposees
# Coord : Coordonnees des points 
# xlim, ylim : limites des axes
# Lwd : epaisseur du trait entre deux points
# Lref : Longueur par laquelle on divise QPlafond pour standardiser l'echelle de couleurs
# Le resultat est une liste avec les composantes suivantes :
# - Parametres : Parametrage de l'algorithme + indicateurs divers pour le probleme
# (Nombre de chemins possibles, QPlafond, la matrice de distance est-elle symetrique ?, ...)
# - Pheromones : Matrice des Pheromones.
# - Matrice des probabilites de passage d'un point a l'autre
# Pour les deux matrices precedentes : si la diagonale est non nulle, le chemin n'est pas une boucle
# La diagonale est alors soit la quantite de pheromone pour le choix d'un point de depart,
# soit la probabilite de choisir ce point pour le depart
# - Best : liste pour le meilleur chemin rencontre avec les composantes suivantes :
#   - Iteration : iterations au cours desquelle au moins une fourmi a suivi le meilleur chemin rencontre
#   - Chemin : meilleur chemin rencontre
#   - Longueur : longueur du meilleur chemin rencontre
#   - Historique : Longueur du meilleur chemin rencontre depuis le debut a chaque iteration
fColonie<-function(MatDist=NULL,Groupe=NULL,Start=NULL,MatPher=NULL,Nfourmi=0,Niter=0,Alpha=0,Beta=0,Gamma=0,Qtot=1,Rho=0,Boucle=FALSE,
                   txQmin=0,txQmax=Inf,Taux=1,Renforcement=TRUE,kPermute=1,Plot=FALSE,Coord=NULL,xlim=c(0,1),ylim=c(0,1),Lwd=5,Lref=NA){
  Sortie<-NULL
  # Calcul de la matrice de probabilite
  fMatProba<-function(MatPher,MatVisi,Alpha,Gamma,Groupe=NULL){
    if (sum(abs(MatPher))==0) {
      MatProb<-Gamma+MatVisi
      DiagP<-diag(MatProb)
    }
    else {
      MatProb<-Gamma+(MatPher^Alpha)*MatVisi
      DiagP<-Gamma+diag(MatPher)^Alpha
    }
    diag(MatProb)<-0
    if (!is.null(Groupe)){
      ListeGroupe<-unique(Groupe)
      if (length(ListeGroupe)<length(Groupe)){
        for (kG in ListeGroupe){
          IdGroupe<-names(Groupe)[Groupe %in% kG]
          MatProb[IdGroupe,IdGroupe]<-0
        }
      }
    }
    MatProb<- MatProb / rowSums(MatProb)
    MatProb<-ifelse(is.nan(MatProb)|is.na(MatProb),0,MatProb)
    if (sum(DiagP)>0) {diag(MatProb)<-DiagP/sum(DiagP)}
    return(MatProb) 
  }
  ### Pour le calcul des permutations
  # Permutations des composantes d'un vecteur
  permutation<-function(x=NULL){
    Sortie<-NULL
    if (!is.null(x)){
      if (length(x)<=1){ Sortie<-x}
      else {
        Sortie<-matrix(data = x[1])
        for (k in 2:length(x)){
          Base<-Sortie
          Vec<-rep(x[k],nrow(Sortie))
          Sortie<-cbind(Vec,Base)
          if (ncol(Base)>1){
            for (Id in 1:(ncol(Base)-1)){
              Temp<-cbind(Base[,1:Id],Vec,Base[,(Id+1):ncol(Base)])
              Sortie<-rbind(Temp,Sortie)
            }
          }
          Temp<-cbind(Base,Vec)
          Sortie<-rbind(Temp,Sortie)
        }
        dimnames(Sortie)<-NULL
      }
    }
    return(Sortie)
  }
  # Permutation de k elements contigus dans un vecteur
  fPermuteContigu<-function(x=NULL,k=1){
    Sortie<-NULL
    if (!is.null(x)){
      if (k<=1){ Sortie<-x}
      else {
        if (k>=length(x)){ Sortie<-permutation(x)}
        else {
          Perm<-permutation(x = 0:(k-1))
          Perm<-Perm[-1,]
          if (is.vector(Perm)){ Perm<-matrix(data = Perm, nrow=1)}
          Sortie<-matrix(data = x, nrow=1)
          for (Id in 1:(length(x)-k+1)){
            Base<-matrix(data = x, byrow = TRUE, ncol=length(x),nrow=nrow(Perm))
            Temp<-t(apply(X=1+Perm,MARGIN=1,FUN=function(a){return((x[Id:(Id+k-1)])[a])}))
            Base[,Id:(Id+k-1)]<-Temp
            Sortie<-rbind(Sortie,Base)
          }
        }
      }
    }
    return(Sortie)
  }
  # Permutations entre k sommets contigus d'un chemin
  ftestPermute<-function(Chemin=NULL,MatDist=NULL,k=1,FixedStart=FALSE){
    Sortie<-NULL
    if (!is.null(Chemin)&!is.null(MatDist)){
      Boucle<- Chemin[1]==Chemin[length(Chemin)]
      if (Boucle){ Chemin<-Chemin[-length(Chemin)]}
      if (FixedStart){
        Start<-Chemin[1]
        Chemin<-Chemin[-1]
      }
      MatChemin<-fPermuteContigu(x=Chemin,k=k)
      if (is.vector(MatChemin)){ MatChemin<-matrix(data = MatChemin, nrow=1)}
      if (FixedStart){ MatChemin<-cbind(rep(Start,nrow(MatChemin)),MatChemin) }
      if (Boucle){ MatChemin<-cbind(MatChemin,MatChemin[,1])}
      LChemin<-apply(X=MatChemin,MARGIN=1,MatDist=MatDist,FUN=lchemin)
      Keep<-which.min(LChemin)
      Sortie<-list(Chemin=MatChemin[Keep,],Longueur=LChemin[Keep])
    }
    return(Sortie)
  }
  # Idem pour une matrice de chemins
  fMatTestPermute<-function(MatChemin=NULL,MatDist=NULL,k=1,FixedStart=FALSE){
    Sortie<-NULL
    if (!is.null(MatChemin)&!is.null(MatDist)){
      Liste<-apply(X=MatChemin,MARGIN=1,FUN=ftestPermute,MatDist=MatDist,k=k,FixedStart=FixedStart)
      MatChemin<-lapply(X=Liste,FUN=function(x){return(x$Chemin)})
      MatChemin<-matrix(data=unlist(MatChemin),byrow=TRUE,nrow=length(Liste))
      VecLongueur<-lapply(X=Liste, function(x){return(x$Longueur)})
      VecLongueur<-unlist(VecLongueur)
      Sortie<-list(MatChemin=MatChemin,Longueur=VecLongueur)
    }
    return(Sortie)
  }
  #### Debut de l'algorithme ####
  if (is.null(MatDist) & !is.null(Coord)){ MatDist<-dist(x=Cord, method = "euclidean", diag = FALSE, upper = FALSE, p = 2) } 
  if (!is.null(MatDist)){
    Rho<-min(1,max(0,Rho))
    Taux<-min(1,max(0,Taux))
    txQmin<-min(1,max(0,txQmin))
    txQmax<-max(0,txQmax)
    Label<-labels(MatDist)
    if (is.list(Label)){ Label<-Label[[1]] }
    ListLabel<-list("Depart"=Label,"Arrivee"=Label)
    Keep<-rep(FALSE,Nfourmi)
    Keep[1:ceiling(Taux*Nfourmi)]<-TRUE
    FixedStart<- !is.null(Start)
    # Matrice de visibilite
    MatVisi<- as.matrix(1/MatDist)^Beta
    Sym<-isSymmetric(matrix(data=MatVisi,nrow=nrow(MatVisi)))
    diag(MatVisi)<-0
    dimnames(MatVisi)<-ListLabel
    NbChemin<-fnbChemin(n=nrow(MatVisi)-as.numeric(FixedStart&!Boucle),Boucle=Boucle,Sym=Sym,Groupe=Groupe)
    if (is.null(Groupe)){ NbSommet<-nrow(MatVisi)}
    else { NbSommet<-length(unique(Groupe))}
    NbSommet<-NbSommet+as.numeric(Boucle)-as.numeric(FixedStart&!Boucle)
    NbPermute<-fnbpermute(n=NbSommet,k=kPermute,Boucle = Boucle)
    # Initialisation de la matrice de pheromones
    if (is.null(MatPher)){MatPher<-matrix(data=0,nrow=nrow(MatVisi),ncol=ncol(MatVisi),dimnames=ListLabel)}
    if (Nfourmi>0 & Niter>0){
      k<-0
      Continue<-TRUE
      BestLongueur<-Inf
      BestChemin<-NULL
      BestIter<-k
      MatProb<-fMatProba(MatPher=MatPher,MatVisi=MatVisi,Alpha=Alpha,Gamma=Gamma,Groupe=Groupe)
      QPlafond<-fPlafond(Nfourmi=Nfourmi,Qtot=Qtot,Rho=Rho,Niter=Niter,Taux=Taux)
      Historique<-c()
      # Par defaut on standardise QPlafond avec la longueur moyenne d'un chemin
      if (is.na(Lref)){
        Lref<-c(MatDist)
        Lref<-mean(Lref[Lref>0])*(length(unique(Groupe))+as.numeric(Boucle)-1)
      }
      if (Plot &!is.null(Coord)){
        ListPlot<-as.list(formals(fPlot))
        Titre<-"Situation initiale"
        for (Nom in names(ListPlot)){ ListPlot[[Nom]]<-get(Nom)}
        ListPlot$MatPher<-MatPher
        OutPlot<-do.call(what = fPlot, args = ListPlot)
      }
      while(Continue){
        k<-k+1
        # Cheminement des fourmis
        MatChemin<-replicate(n=Nfourmi,expr=fCheminement(Start=Start,Groupe=Groupe,MatProb=MatProb,Boucle=Boucle))
        # Evaluation
        if (kPermute<=1){ Longueur<-apply(X=MatChemin,MARGIN=2,FUN=lchemin,MatDist=MatDist) }
        else {
          Liste<-fMatTestPermute(MatChemin = t(MatChemin), MatDist = MatDist, k=kPermute, FixedStart=FixedStart)
          MatChemin<-t(Liste$MatChemin)
          Longueur<-Liste$Longueur
        }
        # On garde le meilleur chemin rencontre
        LongMin<-min(Longueur)
        Historique<-c(Historique,min(BestLongueur,LongMin))
        if (LongMin<=BestLongueur){
          BestInd<-which(Longueur==LongMin)
          BestTemp<-t(MatChemin[,BestInd])
          if (LongMin==BestLongueur){
            BestChemin<-unique(rbind(BestChemin,BestTemp))
            BestIter<-c(BestIter,k)
          }
          else {
            BestChemin<- BestTemp
            BestLongueur<-LongMin
            BestIter<-k
          }
          # Si chemin ferme, on filtre pour les solutions redondantes
          if (Boucle) { BestChemin<-filtreCycle(BestChemin)}
          # Si symetrie, parcourir le chemin dans un sens ou dans l'autre est equivalent. On elimine l'eventuel doublon
          if (Sym){ BestChemin<-filtreSym(BestChemin) }
          # Par securite
          BestChemin<-unique(BestChemin)
        }
        # Depot des pheromones
        Ordre<-order(Longueur)
        Garde<-(1:Nfourmi)[Ordre[Keep]]
        Depot<-matrix(data=0,nrow=length(Label),ncol=length(Label),dimnames=ListLabel)
        # Seules les meilleures fourmis contribuent
        for (k1 in Garde){
          Dep<-MatChemin[1,k1]
          # Si le chemin n'est pas une boucle, la diagonale renseigne sur le point de depart
          if (!Boucle) {
            Depot[Dep,Dep]<-Depot[Dep,Dep]+Qtot/Longueur[k1]
            if (Sym) {
              Last<-MatChemin[nrow(MatChemin),k1]
              Depot[Last,Last]<-Depot[Last,Last]+Qtot/Longueur[k1]
            }
          }
          for (k2 in 2:nrow(MatChemin)){
            Arr<-MatChemin[k2,k1]
            Depot[Dep,Arr]<-Depot[Dep,Arr]+Qtot/Longueur[k1]
            if (Sym) { Depot[Arr,Dep]<-Depot[Arr,Dep]+Qtot/Longueur[k1] }
            Dep<-Arr
          }		
        }
        # Si Renforcement, le meilleur chemin rencontre participe systematiquement
        if (Renforcement){
          for (k1 in 1:nrow(BestChemin)){
            Dep<-BestChemin[k1,1]
            if (!Boucle) {
              Depot[Dep,Dep]<-Depot[Dep,Dep]+Qtot/BestLongueur
              if (Sym) {
                Last<-MatChemin[k1,ncol(BestChemin)]
                Depot[Last,Last]<-Depot[Last,Last]+Qtot/BestLongueur
              }
            }
            for (k2 in 2:ncol(BestChemin)){
              Arr<-BestChemin[k1,k2]
              Depot[Dep,Arr]<-Depot[Dep,Arr]+Qtot/BestLongueur 
              if (Sym) { Depot[Arr,Dep]<-Depot[Arr,Dep]+Qtot/BestLongueur }
              Dep<-Arr
            }		
          }
        }
        # Mise a jour de la matrice de pheromones
        MatPher<- Rho*MatPher + Depot
        MatPher<-ifelse(MatPher>txQmax*QPlafond,txQmax*QPlafond,MatPher)
        MatPher<-ifelse(MatPher<txQmin*QPlafond,txQmin*QPlafond,MatPher)
        if (!is.null(Groupe)){
          ListeGroupe<-unique(Groupe)
          if (length(ListeGroupe)<length(Groupe)){
            DiagPher<-diag(MatPher)
            for (kG in ListeGroupe){
              IdGroupe<-names(Groupe)[Groupe %in% kG]
              MatPher[IdGroupe,IdGroupe]<-0
            }
            diag(MatPher)<-DiagPher
          }
        }
        
        MatProb<-fMatProba(MatPher=MatPher,MatVisi=MatVisi,Alpha=Alpha,Gamma=Gamma,Groupe = Groupe)
        if (Boucle) {
          diag(MatPher)<-0
          diag(MatProb)<-0
        }
        if (Plot &!is.null(Coord)){
          ListPlot$Titre<-paste("Iteration nÂ°",k," : meilleur chemin rencontre = ",format(BestLongueur,digits=5),sep="")
          ListPlot$MatPher<-MatPher
          OutPlot<-do.call(what = fPlot, args = ListPlot)
          VecCol<-"black"
          nBest<-nrow(BestChemin)
          if (nBest>1) {VecCol<-rainbow(n=nBest)}
          for (kb in 1:nBest){ lines(x=Coord[BestChemin[kb,],],col=VecCol[kb]) }
        }
        Continue<-(k<Niter)
      }
      if (k>0) {names(Historique)<-1:k}
    }
    Best<-list(Iteration=BestIter,Chemin=BestChemin,Longueur=BestLongueur,Historique=Historique)
    Parametres<-list(Nfourmi=Nfourmi,Niter=Niter,Alpha=Alpha,Beta=Beta,Gamma=Gamma,Qtot=Qtot,Rho=Rho,Boucle=Boucle,
                     txQmin=txQmin,txQmax=txQmax,Taux=Taux,Sym=Sym,QPlafond=QPlafond,Lref=Lref,NbChemin=NbChemin,
                     kPermute=kPermute,NbPermute=NbPermute)
    Sortie<-list(Parametres=Parametres,Pheromones=MatPher,Proba=MatProb,Best=Best)
  }
  return(Sortie)
}

# Reconstitution du meilleur chemin a partir de la matrice de sortie
# Matrice : Matrice a partir de laquelle contruire ce meilleur chemin.
# Soit le chemin parcouru en suivant le maximum de pheromones
# Soit le chemin parcouru en suivant le maximum des probabilites
# Dans tous les cas, le programme determine si le chemin fait une boucle
# et s'il est indifferent de suivre un chemin dans un sens ou dans l'autre
# a partir de la structure de la matrice (composante diagonale + symetrie)
# La sortie peut etre constituee de plusieurs chemins
# Dans ce dernier cas, il faudra verifier s'ils sont ex aequo ou non
reBuildChemin<-function(Matrice=NULL,Groupe=NULL){
  Chemin<-NULL
  nextStep<-function(Chemin,Matrice,Groupe){
    if (is.null(Groupe)){ 
      Groupe<-rownames(Matrice)
      if (is.null(Groupe)){ Groupe<-1:nrow(Matrice)}
      names(Groupe)<-Groupe
    }
    Label<-names(Groupe)
    if (is.null(Label)) {
      Label<-1:nrow(Matrice)
      dimnames(Matrice)<-list(Label,Label)
    }
    Reste<-Label[!(Groupe %in% Groupe[Chemin])]
    Last<-Chemin[length(Chemin)]
    Score<-Matrice[Last,Reste]
    Next<-Reste[which(Score == max(Score)) ] 
    Sortie<-cbind(matrix(data=rep(Chemin,length(Next)),byrow=TRUE,nrow=length(Next)),Next)
    return(t(Sortie))
  }
  if (!is.null(Matrice)) {
    if (is.null(Groupe)){ 
      Groupe<-rownames(Matrice)
      if (is.null(Groupe)){ Groupe<-1:nrow(Matrice)}
      names(Groupe)<-Groupe
    }
    Label<-names(Groupe)
    NbGroupe<-length(unique(Groupe))
    # Chemin sous forme de boucle ou non ?
    Boucle<-all(diag(Matrice)<=0) ;
    # Chemin oriente ?
    # C'est tordu mais il faut se debarasser du nom des dimensions
    Sym<-isSymmetric(matrix(data=Matrice,nrow=nrow(Matrice)))
    # Depart
    if (Boucle) { Chemin<-Label[1] }
    else { Chemin <- Label[ which( diag(Matrice) == max(diag(Matrice))) ] }
    Chemin<-matrix(data=Chemin,nrow=1)
    if (NbGroupe>1){
      for (k in 2:NbGroupe) {
        Chemin<-apply(X=Chemin,MARGIN=2,FUN=nextStep,Matrice=Matrice,Groupe=Groupe)
        if (is.list(Chemin)){ Chemin<-matrix(data = unlist(Chemin),nrow=k) }
        Chemin<-matrix(data = Chemin,nrow=k)
      }
    }
    if (Boucle) { Chemin<-rbind(Chemin,Chemin[1,]) }
    Chemin<-unique(t(Chemin))
    if (Boucle) { Chemin<-filtreCycle(Chemin)}
    if (Sym){ Chemin<-filtreSym(Chemin) }
  }
  return(Chemin)
}
