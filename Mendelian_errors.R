#load packages
library(data.table)
library(plyr)
library(dplyr)

#read in map file, genotype file with all chromosomes as matrix file and pedigree
map<-read.table('Input_files/W1WC_new.map')
genom<-read.table('Input_files/W1WC.geno',header=FALSE)
correct_animals<-read.table('Input_files/data_Lisette.txt',header=TRUE)
correct_snps<-read.table('Input_files/keep_WC.txt')
ped<-read.delim("Input_files/ped.txt",header=FALSE)
colnames(ped)<-c('self','father','mother')

colnames(genom)<-c('self', as.character(map[,2]))
genom<-as.matrix(genom)

#############################1####################################
#remove animals with pedigree errors etc. 
genom_offs<-subset(genom, genom[,1]>=8000000)
genom_par<-subset(genom, genom[,1]<8000000)
genom_offs<-subset(genom_offs, genom_offs[,1] %in% correct_animals[,1])
genom<-rbind(genom_par,genom_offs)

#remove snps after quality check
wrong_snps<-which(!colnames(genom) %in% as.character(correct_snps[,1]))
wrong_snps[1]<-0
genom<-genom[,-c(wrong_snps)]

##final files
ped<-subset(ped, self %in% genom[,1])
genom<-subset(genom, genom[,1] %in% ped$self)

############################2#######################################
#find all mismatches for each triplet of offspring, sire, dam
nIDs<-nrow(ped)
mismatch<-(matrix(as.integer(0),nrow=dim(genom)[1],ncol=dim(genom)[2]))            #store all single mismatches
n_evaluated<-0
n_notevaluated<-0
for (self in 1:nIDs) {                                                    
  geno_self<-genom[self,]
  sire_pos<-which(ped$self==ped$father[self])
  dam_pos<-which(ped$self==ped$mother[self])
  if((length(sire_pos)==1)&(length(dam_pos)==1)) {                                  #both sire and dam genotyped
    n_evaluated<-n_evaluated+1
    geno_sire<-genom[sire_pos,]
    geno_dam<-genom[dam_pos,]
    mismatch[self,which(geno_sire==0&geno_dam==0&geno_self!=0&geno_self!=9)]<-as.integer(1)                 #store mismatches per self, missing self never gives an error
    mismatch[self,which(geno_sire==0&geno_dam==1&(geno_self==2))]<-as.integer(1)  #as.integer saves half the memory             
    mismatch[self,which(geno_sire==0&geno_dam==2&geno_self!=1&geno_self!=9)]<-as.integer(1)               
    mismatch[self,which(geno_sire==1&geno_dam==0&(geno_self==2))]<-as.integer(1)               
    mismatch[self,which(geno_sire==1&geno_dam==2&(geno_self==0))]<-as.integer(1)            
    mismatch[self,which(geno_sire==2&geno_dam==0&geno_self!=1&geno_self!=9)]<-as.integer(1)              
    mismatch[self,which(geno_sire==2&geno_dam==1&(geno_self==0))]<-as.integer(1)              
    mismatch[self,which(geno_sire==2&geno_dam==2&geno_self!=2&geno_self!=9)]<-as.integer(1)              
    mismatch[self,which(geno_sire==9&geno_dam==0&(geno_self==2))]<-as.integer(1)               
    mismatch[self,which(geno_sire==9&geno_dam==2&(geno_self==0))]<-as.integer(1)               
    mismatch[self,which(geno_sire==0&geno_dam==9&(geno_self==2))]<-as.integer(1)               
    mismatch[self,which(geno_sire==2&geno_dam==9&(geno_self==0))]<-as.integer(1)
  }
  if((length(sire_pos)==1)&(length(dam_pos)==0)) {                                  #only sire genotyped
    n_evaluated<-n_evaluated+1
    geno_sire<-genom[sire_pos,]
    mismatch[self,which(geno_sire==0&geno_self==2)]<-as.integer(1)              
    mismatch[self,which(geno_sire==2&geno_self==0)]<-as.integer(1)              
  }
  if((length(sire_pos)==0)&(length(dam_pos)==1)) {                                  #only dam genotyped
    n_evaluated<-n_evaluated+1
    geno_dam<-genom[dam_pos,]
    mismatch[self,which(geno_dam==0&geno_self==2)]<-as.integer(1)              
    mismatch[self,which(geno_dam==2&geno_self==0)]<-as.integer(1)              
  } 
  if((length(sire_pos)==0)&(length(dam_pos)==0)) {                                  #neither parent genotyped
    n_notevaluated<-n_notevaluated+1
  } 
}
rm(geno_self,geno_dam,geno_sire,sire_pos,dam_pos)

colnames(mismatch)<-colnames(genom)
mismatch[,1]<-genom[,1]
rownames(mismatch)<-mismatch[,1]

n_evaluated                                                                          
n_notevaluated                														#animals without parent genotyped

n_mismatch_locus<-apply(mismatch[,2:ncol(mismatch)],2,sum)                         	#total number of mismatches for each locus
n_mismatch_animal<-apply(mismatch[,2:ncol(mismatch)],1,sum) 

#################################2a##############################
n_mismatch_animal<-as.data.frame(n_mismatch_animal)
n_mismatch_animal$self<-genom[,1]
n_mismatch_animal$father<-ped[,2]
many_mismatches<-subset(n_mismatch_animal, n_mismatch_animal > 0.001*(ncol(genom)-1))
many_mismatches_markers<-subset(n_mismatch_locus,n_mismatch_locus >0.005*nrow(genom))
bad_dads<- names(which(table(many_mismatches$father)>=3))

genom2<-genom[!genom[,1] %in% bad_dads, ]
ped2<-ped[!ped[,1]%in%bad_dads,]

genom<-genom2
ped<-ped2

###############################2b##############################
bad_markers<-names(which(n_mismatch_locus>=0.005*nrow(genom)))
bad_IDs<-names(which(n_mismatch_animal>=0.001*(ncol(genom)-1)))

genom3<-genom[!genom[,1] %in% bad_IDs,]
genom3<-genom3[,!colnames(genom)%in%bad_markers]

mismatch2<-mismatch[!mismatch[,1] %in% bad_IDs, ]
mismatch2<-mismatch2[,!colnames(mismatch2) %in% bad_markers]

genom3[which(mismatch2==1)]<-as.integer(9)							#final file
map2<-map[map$V2%in%colnames(genom3),]
ped3<-ped[!ped[,1]%in% bad_IDs,]

#############################3###############################
#remove offspring with only 1 parent
id_2<-numeric(nrow(ped3))
for (self in 1:nrow(ped3)){
sire_pos<-which(ped3$self==ped3$father[self])
dam_pos<-which(ped3$self==ped3$mother[self])
  if((length(sire_pos)==1)&(length(dam_pos)==1)) {  
 id_2[self]<-ped3$self[self]}}
id_2<-as.integer(id_2)

Cpar<-subset(genom3, genom3[,1]<8000000,drop=FALSE)
Coff<-subset(genom3, genom3[,1]>8000000,drop=FALSE)
Coff2<-subset(Coff, Coff[,1] %in% id_2)
genom4<-rbind(Cpar,Coff2)

ped4<-ped3[ped3[,1] %in% genom4[,1],]
#############################4################################
write.table(map2, "Output_files/W1WD/Data_preparation/map2.txt", sep="\t",quote=FALSE)
write.table(ped4, "Output_files/W1WD/Data_preparation/ped4.txt", sep="\t",quote=FALSE)
write.table(genom4,"Output_files/W1WD/Data_preparation/genom4.txt", sep="\t",quote=FALSE)

for(i in seq_len(28)) {
  filename = paste("Output_files/W1WD/Data_preparation/Chromosomes/Chr", i,"_geno.txt", sep="")
  idx <- match(subset(map2,map2[,1]==i)[,2], colnames(genom4))
Subgenom <- genom4[,c(1,idx)]
  write.table(Subgenom, filename, col.names=FALSE, row.names=FALSE, quote=FALSE)
}

for(i in seq_len(28)) {
filename = paste("Output_files/W1WD/Data_preparation/Map_files/Chr",i,"_map.txt",sep="")
write.table(subset(map2,map2[,1]==i),filename,row.names=FALSE,quote=FALSE)
}
