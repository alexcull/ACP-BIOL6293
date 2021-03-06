#Pr�sentation PCA
#Alex Cull
#G�n�ration des donn�es pour graphique de pr�sentation ----
ls()
rm(list = ls())
ls()

Var.1 <- rnorm(6, mean = 10, sd = 2)
Var.2 <- rnorm(6, mean = 8, sd = 3)
Var.3 <- rnorm(6, mean = 3, sd = 0.8)
Var.4 <- rnorm(6, mean = 0.5, sd = 0.25)
Var.1 <- round(Var.1, 1)
Var.2 <- round(Var.2, 1)
Var.3 <- round(Var.3, 1)
Var.4 <- round(Var.4, 1)
Var.1
#[1] 10.0 11.1  9.6 13.7 11.3 12.2
Var.2
#[1] 8.3 6.8 4.4 3.8 3.8 3.5
Var.3
#[1] 3.7 2.9 1.4 2.9 3.5 3.8
Var.4
#[1] 0.9 0.7 0.4 0.3 0.5 0.3
DF <- data.frame(Var.1, Var.2, Var.3, Var.4)
DF <-t(DF)
write.csv(x = DF, file = 'DFt.csv')
#PCA dans R avec prcomp ----
ls()
rm(list = ls())
ls()
donn�es <- matrix(nrow = 100, ncol=10)
colnames(donn�es) <- c(
  paste("ct",1:5, sep = ""),
  paste("tt", 1:5, sep=""))
rownames(donn�es) <- paste("g�ne", 1:100, sep= "")
for(i in 1:100) {
  valeurs.ct <- rpois(5, lambda = sample(x=10:1000, size =1))
  valeurs.tt <- rpois(5, lambda = sample(x=10:1000, size =1))
  
  donn�es[i,] <- c(valeurs.ct, valeurs.tt)
}
head(donn�es)
#prcomp veux des �chantillons dans les colonnes et les variables dans les rang�es pour donner x, sdev 
#et rotation. x = PC. Le 1er repr�sente le plus de variation, 2e est le 2e. 
pca <- prcomp(t(donn�es), scale = TRUE)
plot(pca$x[,1], pca$x[,2])

#Variance des PC
pca.var <- pca$sdev^2
#On transforme en pourcentage
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
#Scree Plot
barplot(pca.var.per, main = 'Scree Plot',
        xlab = "PC", 
        ylab = "Variation en %")

#ggplot
library("ggplot2", lib.loc="~/R/win-library/3.5")
#Formattage pour ggplot
pca.donn�es <- data.frame(Sample = rownames(pca$x),
                       X = pca$x[,1],
                       Y= pca$x[,2])
pca.donn�es


ggplot(data=pca.donn�es, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pca.var.per[1], "%", sep = " ")) + #titre d'axes
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep = " ")) + # % variation
  ggtitle("Tableau PCA") #titre

#Loading scores = vecteur propre = eigenvector = rotation
loading_scores <- pca$rotation[,1]
#scores (-) d�place les �chantillons vers la gauche, (+) vers la droite
g�ne_scores <-abs(loading_scores)
g�ne_score_en_ordre <- sort(g�ne_scores, decreasing = TRUE)
top_10_g�nes <- names(g�ne_score_ranked[1:10])
#Top 10 g�nes ayant le plus d'impact
top_10_g�nes
#Quels effets ont les top 10 g�nes
pca$rotation[top_10_g�nes,1] ## show the scores and + or -




