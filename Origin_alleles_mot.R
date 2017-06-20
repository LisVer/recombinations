#####################################1##################################
#load data
C<- read.table("Input_files/C1_W1WC_70_65_NO2.txt", sep=" ",header=FALSE) ##finalphase.txt output from AlphaPhase1.2 #2
ped<-read.table("Input_files/ped5.txt")  ##output from Phasing_mistakes.R
map<-read.table("Input_files/Chr1_map2.txt", header=TRUE) ##output from Phasing_mistakes.R

#change data
C<-C[,-c(2:14)]
C<-as.matrix(C)
map<-map[,-3]
colnames(map)<-c('chr','marker','pos')
colnames(C)<-c('self', as.character(map$marker))

Coff<-subset(C, C[,1]>8000000,drop=FALSE)
even_indexes<-seq(2,nrow(Coff),2)
odd_indexes<-seq(1,nrow(Coff)-1,2)
Cfat<-Coff[odd_indexes,]				#different files to separate the father and mother line
Cmot<-Coff[even_indexes,]

#see where allele in offspring comes from
n_IDs<-nrow(ped)
PhasMisMot<-(matrix(as.integer(0),nrow=dim(ped)[1],ncol=dim(Cmot)[2]))   

for (self in 1:n_IDs) {                                                     #paar minutjes
  dam_pos1<-which(ped$self==ped$mother[self])*2-1
  dam_pos2<-which(ped$self==ped$mother[self])*2
  if (length(dam_pos1)>0 & length(dam_pos2)>0){  	#if the parental haplotypes exist
  geno_self<-C[seq(2,nrow(C),by=2),][self,]
  geno_dam1<-C[dam_pos1,]
  geno_dam2<-C[dam_pos2,]
	PhasMisMot[self,which(geno_self==geno_dam1&geno_self!=geno_dam2)]<-as.integer(1)  #allele from paternal haplotype 
	PhasMisMot[self,which(geno_self==geno_dam2&geno_self!=geno_dam1)]<-as.integer(2)   #allele from maternal haplotype
	PhasMisMot[self,which(geno_self==geno_dam2&geno_self==geno_dam1)]<-as.integer(0)   #allele can come from both haplotypes (homozygous)
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
write.table(PhasMisMot2, 'Output_files/W1WC/PhasMisMot2.txt')

###########################2########################################	 
#count frequencies of allele origins per SNP
CountsPerSNP<-apply(PhasMisMot2, 2, function(x){
  (as.data.frame(table(factor(x, lev=c(0,1,2,5,9)))))/nrow(PhasMisMot2)})
CountsPerSNP.df<-as.data.frame(CountsPerSNP)
even_indexes<-seq(2,ncol(CountsPerSNP.df),2)
CountsPerSNP.df<-CountsPerSNP.df[,c(1,even_indexes)]  
CountsPerSNP.df[,1]<-c(0,1,2,5,9)
#count sum of frequencies
CountsPerSNP.df[6,]<-apply(CountsPerSNP.df,2, sum)
CountsPerSNP.df[6,1]<-"sum"
colnames(CountsPerSNP.df)<-c("number",colnames(PhasMisMot2))
CountsPerSNP.df$self<-NULL

#frequencies
sums<-apply(CountsPerSNP.df[,2:ncol(CountsPerSNP.df)], 1, function(x){
sum(x)/(ncol(PhasMisMot2)-1)})
sums<-t(t(sums))
colnames(sums)<-"proportions"
sums<-as.data.frame(sums)
sums[7,]<-sums[4,]/(sums[2,]+sums[3,])
rownames(sums)<-c(0,1,2,5,9,"total","5/(1+2)")

#count frequencies of allele origins per ID
CountsPerID<-apply(PhasMisMot2, 1, function(x){
  (as.data.frame(table(factor(x, lev=c(0,1,2,5,9)))))/(ncol(PhasMisMot2)-1)})
CountsPerID.df<-as.data.frame(CountsPerID)
even_indexes<-seq(2,ncol(CountsPerID.df),2)
CountsPerID.df<-CountsPerID.df[,c(1,even_indexes)]  
CountsPerID.df[,1]<-c(0,1,2,5,9)
#count sum of frequencies
CountsPerID.df[6,]<-apply(CountsPerID.df,2, sum)
CountsPerID.df[7,]<-CountsPerID.df[4,]/(CountsPerID.df[2,]+CountsPerID.df[3,])
CountsPerID.df[6,1]<-"sum"
CountsPerID.df[7,1]<-"5/(1+2)"
colnames(CountsPerID.df)<-c("number",PhasMisMot2[,1])

#look for individuals with most mistakes (code 5)
Mistakes_ID<-as.data.frame(t(CountsPerID.df[4,2:ncol(CountsPerID.df)]))
colnames(Mistakes_ID)<-'mistakes'
Mistakes_ID<-Mistakes_ID[with(Mistakes_ID, order(-Mistakes_ID$mistakes)), ,drop=FALSE]
Mistakes_ID_scaled<-as.data.frame(t(CountsPerID.df[7,2:ncol(CountsPerID.df)]))
colnames(Mistakes_ID_scaled)<-'mistakes'
Mistakes_ID_scaled<-Mistakes_ID_scaled[with(Mistakes_ID_scaled, order(-Mistakes_ID_scaled$mistakes)), ,drop=FALSE]

#look for markers with most mistakes (code 5)
Mistakes_snp_mot<-as.data.frame(t(CountsPerSNP.df[4,2:ncol(CountsPerSNP.df)]))
colnames(Mistakes_snp_mot)<-'mistakes'
Mistakes_snp2<-Mistakes_snp_mot[with(Mistakes_snp_mot, order(-Mistakes_snp_mot$mistakes)), ,drop=FALSE]

#phasing mistake where other haplotype has a phasing mistake
PhasMisMot3<-PhasMisMot2
for (i in 1:nrow(PhasMisMot3)){
PhasMisMot3[i,which(PhasMisFat2[i,]==9)]<-as.integer(9)}

########make plots############
#make histograms
pdf('Output_files/W1WC/WC/hist_mistakes_mot.pdf')
par(mfrow=c(2,1))
hist(as.vector(CountsPerSNP.df[4,2:ncol(CountsPerSNP.df)],mode="numeric"),breaks=100,xlab="wrong phase", main='hist mistakes per marker',cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
hist(as.vector(CountsPerID.df[4,2:ncol(CountsPerID.df)],mode="numeric"),breaks=100,xlab="wrong phase",  main='hist mistakes per individual',cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
dev.off()

#make boxplots
pdf('Output_files/W1WC/WC/boxplots_mistakes_mot.pdf')
par(mfrow=c(1,2))
boxplot(t(CountsPerSNP.df[4,2:ncol(CountsPerSNP.df)]), main='mistakes marker',cex.lab=2, cex.axis=2, cex.main=2)
boxplot(t(CountsPerID.df[4,2:ncol(CountsPerID.df)]), main='mistakes ind.',cex.lab=2, cex.axis=2, cex.main=2)
dev.off()

#make barplot
pdf('Output_files/W1WC/WC/barplot_mistakes_mot.pdf')
barplot(as.vector(CountsPerSNP.df[4,2:ncol(CountsPerSNP.df)], mode='numeric'),
        names.arg=map$pos,main="Chr1 phasing mistakes",xlab="position in bp",ylab="frequency",cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
dev.off()

#scaled plots
pdf('Output_files/W1WC/WC/scaled_mistakes.pdf')
par(mfrow=c(2,1))
hist(as.vector(CountsPerID.df[7,2:ncol(CountsPerID.df)],mode="numeric"),breaks=100,xlab="wrong phase", main='scaled hist mistakes per individual',cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
boxplot(t(CountsPerID.df[7,2:ncol(CountsPerID.df)]), main='scaled boxplot mistakes per individual',cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
dev.off()

#write tables
write.table(Mistakes_ID, "Output_files/W1WC/WC/Mistakes_ID.txt")
write.table(Mistakes_snp2, "Output_files/W1WC/WC/Mistakes_snp.txt")
write.table(Mistakes_ID_scaled, "Output_files/W1WC/WC/Mistakes_ID_scaled.txt")
write.table(PhasMisMot3, "Output_files/W1WC/WC/PhasMisMot3.txt")
write.table(sums, 'Output_files/W1WC/WC/frequencies.txt')

