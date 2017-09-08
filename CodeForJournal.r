#####MANTEL TESTS FOR CORRELATION BETWEEN TAXONOMIC DISTANCE AND INTERACTION STRUCTURE

require(prabclus)
require(ape)
require(vegan)
require(ade4)

#TAXONOMIC DISTANCE MATRIX
#If two species belong to the same GENUS, their distance is 0
#If they belong to the same FAMILY but are from different GENERA, their distance is 1
#If they belong to the same ORDER but are from different FAMILIES, their distance is 2
#And so on...

##############################################################

int=read.table("interactionmatrix.txt",header=T)   #binary interaction matrix for plants (or insects)
mat=as.matrix(int)        #coerce to matrix
jac=jaccard(mat)          #use jaccard algorithm to make distance matrix

taxonomies <-read.table("taxonomies.txt",header=TRUE) #taxonomic hierarchy of plants (or insects)
attach(taxonomies)              #a three column dataframe of the family, subfamily, and genus each species belongs to
dO<-weight.taxo(Order)
dF <- weight.taxo(Family)
dG <- weight.taxo(Genus)
dS <- weight.taxo(Species)
Dist <- dO+ dF + dG + dS
taxo <- (max(Dist)-Dist)  #taxo is the matrix of taxonomic distances
detach(taxonomies)

mantel(taxo,jac)      #the mantel test prints a correlation statistic and significance
                      #default is 999 permutations

#This is a one way test, so it can be done to compare the relationship between
#the interaction structure of the insects and their taxonomies, or the plants and
#their taxonomies.



######FOURTH CORNER ALGORITHM FOR CORRELATIONS BETWEEN TAXONOMIES OF INTERACTING
######PARTNERS

tabL<-read.table('interactionmatrix.txt',header=T)
tabR<-read.table('orddist.txt',header=T)            #matrix of insect species in order
tabQ<-read.table('orddist_plants.txt',header=T)     #matrix of plant species in order

ordfourth=fourthcorner(tabR,tabL,tabQ,nrepet=999)

save(ordfourth,file="fourth_ord.RData")

sumfour=summary(ordfourth)

write(sumfour,file="sumfour_ord.txt")

tabL<-read.table('interactionmatrix.txt',header=T) 
tabR<-read.table('famdist.txt',header=T)            #matrix of insect species in families
tabQ<-read.table('famdist_plants.txt',header=T)     #matrix of plant species in families

famfourth=fourthcorner(tabR,tabL,tabQ,nrepet=999)

save(famfourth,file="fourth_fam.RData")

sumfour=summary(famfourth)

write(sumfour,file="sumfour_fam.txt")
