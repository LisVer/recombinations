#load data
C<- read.table("Input_files/C1_W1WC_70_65_NO.txt", sep=" ",header=FALSE) ##finalphase.txt output from AlphaPhase1.2
ped<-read.table("Input_files/ped4.txt")  ##output from Mendelian_errors.R
map2<-read.table("Input_files/map2.txt") ##output from Mendelian_errors.R
Chr1_map<-read.table("Input_files/Chr1_map.txt", header=TRUE) ##output from Mendelian errors.R
genom4<-read.table("Input_files/genom4.txt")  ##output from Mendelian errors.R 

#change data
C<-C[,-c(2:14)]
C<-as.matrix(C)
colnames(Chr1_map)<-c('chr','marker','pos')
colnames(C)<-c('self', as.character(Chr1_map$marker))

########################1#########################################
Coff<-subset(C, C[,1]>8000000,drop=FALSE)
even_indexes<-seq(2,nrow(Coff),2)
odd_indexes<-seq(1,nrow(Coff)-1,2)
Cfat<-Coff[odd_indexes,]				#different files to separate the father and mother line
Cmot<-Coff[even_indexes,]

#see where allele in offspring comes from
n_IDs<-nrow(ped)
PhasMisFat<-(matrix(as.integer(0),nrow=dim(ped)[1],ncol=dim(Cfat)[2]))   
for (self in 1:n_IDs) {                                                     #paar minutjes
  sire_pos1<-which(ped$self==ped$father[self])*2-1
  sire_pos2<-which(ped$self==ped$father[self])*2
  if (length(sire_pos1)>0 & length(sire_pos2)>0){  	#if the parental haplotypes exist
  geno_self<-C[seq(1,nrow(C),by=2),][self,]
  geno_sire1<-C[sire_pos1,]
  geno_sire2<-C[sire_pos2,]
	PhasMisFat[self,which(geno_self!=geno_sire1&geno_self!=geno_sire2)]<-as.integer(5)   #allele cannot come from one of the parent haplotypes
	PhasMisFat[self,which(geno_self==9|geno_sire1==9|geno_sire2==9)]<-as.integer(9)   	#missing value, origin of allele unknown
	} else{
	PhasMisFat[self,]<-as.integer(NA)	# no information on parents
	}}
	
#remove parent animals
PhasMisFat2<-PhasMisFat
PhasMisFat2[,1]<-as.integer(as.character(ped[,1]))
colnames(PhasMisFat2)<-colnames(C)
PhasMisFat2<-subset(PhasMisFat2, PhasMisFat2[,1]>8000000)
PhasMisFat2<-PhasMisFat2[complete.cases(PhasMisFat2),]

###mothers
PhasMisMot<-(matrix(as.integer(0),nrow=dim(ped)[1],ncol=dim(Cmot)[2]))   

for (self in 1:n_IDs) {                                                    
  dam_pos1<-which(ped$self==ped$mother[self])*2-1
  dam_pos2<-which(ped$self==ped$mother[self])*2
  if (length(dam_pos1)>0 & length(dam_pos2)>0){  	#if the parental haplotypes exist
  geno_self<-C[seq(2,nrow(C),by=2),][self,]
  geno_dam1<-C[dam_pos1,]
  geno_dam2<-C[dam_pos2,]
	PhasMisMot[self,which(geno_self!=geno_dam1&geno_self!=geno_dam2)]<-as.integer(5)   #allele cannot come from one of the parent haplotypes
	PhasMisMot[self,which(geno_self==9|geno_dam1==9|geno_dam2==9)]<-as.integer(9)   	#missing value, origin of allele unknown
	} else{
	PhasMisMot[self,]<-NA	# no information on parents
	}}

#remove parent animals
PhasMisMot2<-PhasMisMot
PhasMisMot2[,1]<-as.integer(as.character(ped$self))
colnames(PhasMisMot2)<-colnames(C)
PhasMisMot2<-subset(PhasMisMot2, PhasMisMot2[,1]>8000000)
PhasMisMot2<-PhasMisMot2[complete.cases(PhasMisMot2),]	 

##################################2a######################################
#count frequencies of phasing mistakes per ID fathers
CountsPerIDFat<-matrix(as.integer(0),nrow=dim(PhasMisFat2)[1],ncol=1)
for (i in 1:nrow(PhasMisFat2)){
CountsPerIDFat[i,]<-length(which(PhasMisFat2[i,]==5))/(ncol(PhasMisFat2)-1)}
rownames(CountsPerIDFat)<-PhasMisFat2[,1]
#count frequencies of phasing mistakes per ID mothers
CountsPerIDMot<-matrix(as.integer(0),nrow=dim(PhasMisMot2)[1],ncol=1)
for (i in 1:nrow(PhasMisMot2)){
CountsPerIDMot[i,]<-length(which(PhasMisMot2[i,]==5))/(ncol(PhasMisMot2)-1)}
rownames(CountsPerIDMot)<-PhasMisMot2[,1]
#individuals with many mistakes
ind_mis_fat<-rownames(subset(CountsPerIDFat, CountsPerIDFat[,1]>=0.05))
ind_mis_mot<-rownames(subset(CountsPerIDMot, CountsPerIDMot[,1]>=0.05))
ind_mis<-c(ind_mis_fat, ind_mis_mot)

genom5<-genom4[-which(genom4[,1] %in% ind_mis),]
ped5<-ped[ped[,1] %in% genom5[,1],]
PhasMisFat2<-PhasMisFat2[-which(PhasMisFat2[,1] %in% ind_mis),]
PhasMisMot2<-PhasMisMot2[-which(PhasMisMot2[,1] %in% ind_mis),]

###############################2b##########################################
#count phasing mistakes per SNP fathers
CountsPerSNPFat<-matrix(as.integer(0),nrow=dim(PhasMisFat2)[2],ncol=1)
for (i in 1:ncol(PhasMisFat2)){
CountsPerSNPFat[i,]<-length(which(PhasMisFat2[,i]==5))/(nrow(PhasMisFat2))}
rownames(CountsPerSNPFat)<-colnames(PhasMisFat2)
#count phasing mistakes per SNP mothers
CountsPerSNPMot<-matrix(as.integer(0),nrow=dim(PhasMisMot2)[2],ncol=1)
for (i in 1:ncol(PhasMisMot2)){
CountsPerSNPMot[i,]<-length(which(PhasMisMot2[,i]==5))/(nrow(PhasMisMot2))}
rownames(CountsPerSNPMot)<-colnames(PhasMisMot2)
#snps with many mistakes
CountsPerSNP<-(CountsPerSNPFat+CountsPerSNPMot)/2
snps_mis<-rownames(subset(CountsPerSNP, CountsPerSNP[,1]>=0.05))

genom6<-genom5[,-which(colnames(genom5) %in% snps_mis)]
map3<-map2[-which(as.character(map2[,2]) %in% snps_mis),]
idx <- match(subset(map3,map3[,1]==1)[,2], colnames(genom6))
Chr1_geno2 <- genom6[,c(1,idx)]														#deze nog even checken
Chr1_map2<-Chr1_map[Chr1_map$marker %in%colnames(genom6),]

##############################3#############################################
write.table(Chr1_geno2, 'Output_files/W1WC/Data_preparation/ChromosomesC1/Chr1_geno2.txt',row.names=FALSE,col.names=FALSE)
write.table(Chr1_map2, 'Output_files/W1WC/Data_preparation/Map_filesC1/Chr1_map2.txt')
write.table(ped5, 'Output_files/W1WC/Data_preparation/ped5.txt')
write.table(snps_mis, 'Output_files/W1WC/Data_preparation/snps_mis.txt')
write.table(ind_mis, 'Output_files/W1WC/Data_preparation/ind_mis.txt')
write.table(genom6, 'Output_files/W1WC/Data_preparation/genom6.txt')