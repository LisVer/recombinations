#load data
PhasMisMot3<-read.table('Input_files/PhasMisMot3.txt') ##output file Origin_alleles_fat.R
RecMot2<-as.matrix(PhasMisMot3)
map<-read.table("Input_files/Chr1_map2.txt", header=TRUE) ##output file Phasing_mistakes.R
ped<-read.table("Input_files/ped5.txt") ##output file Phasing_mistakes.R

#load packages
library(data.table)
library(plyr)
library(dplyr)

#prepare data
colnames(map)<-c('chr','marker','pos')
RecMot2<-RecMot2[,-1]

#count subsequent snps with the same origin
seventy<-c(seq(0,ncol(RecMot2),by=70),ncol(RecMot2))

tor<-list()
tor2<-list()
for (j in 1:nrow(RecMot2)){
vec<-RecMot2[j,]
for (k in 1:(ncol(RecMot2)/70+1)){
x1<-0
x2<-0
tor[[k]]<-0
for(i in (1+seventy[k]):(seventy[k+1])){
if (vec[i]==1){
if (x2>0){
tor[[k]]<-c(tor[[k]],x2)}
x2<-0
x1<-x1+1
} else if (vec[i]==2){
if (x1>0){
tor[[k]]<-c(tor[[k]],x1)}
x1<-0
x2<-x2+1
}
}
if (x1>0){
tor[[k]]<-c(tor[[k]],x1)}
else if (x2>0){
tor[[k]]<-c(tor[[k]],x2)}
}
tor2[[j]]<-tor}

#count phase shifts
recomb<-lapply(tor2, function(x){
lapply(x, function(y){
if (length(y)>1){
length(y)-2
}else {y<-0}})})

recomb1b<-lapply(recomb, unlist)
recomb1c<-lapply(recomb1b,sum)
recomb1c<-as.matrix(recomb1c)
rownames(recomb1c)<-PhasMisMot3[,1]

#only count phase shifts when there are at least 6 consecutive snps from the same parent haplotype
##set shifts <5 to 0
tor3<-tor2
for(i in 1:length(tor3)){
for(j in 1: length(tor3[[i]])){
for(k in 1: length(tor3[[i]][[j]])){
if (tor3[[i]][[j]][[k]]<6){
tor3[[i]][[j]][[k]]<-0}}}}
##add up shifts from the same parental haplotype with zeroes in between
tor3b<-lapply(tor3, function(x){
lapply(x, function(y){
y<-c(y,0,0)})})
for (h in 1:length(tor3b)){
for (i in 1:length(tor3b[[h]])){
for (j in 2:length(tor3[[h]][[i]])){
if (length(tor3[[h]][[i]])>2){
if (tor3b[[h]][[i]][[j]]>0 &tor3b[[h]][[i]][[j+1]]==0){
tor3b[[h]][[i]][[j+2]]<-tor3b[[h]][[i]][[j]]+tor3b[[h]][[i]][[j+2]]
tor3b[[h]][[i]][[j]]<-0}}
}}}
##remove zeroes
tor3b<-lapply(tor3b, function(x){
lapply(x, function(y){
y[y!=0]})})
##count phase shifts
recomb3<-lapply(tor3b, function(x){
lapply(x, function(y){
if (length(y)>1){
length(y)-1
}else {y<-0}})
})
recomb3b<-lapply(recomb3, unlist)
recomb3c<-lapply(recomb3b,sum)
recomb3c<-as.matrix(recomb3c)
rownames(recomb3c)<-PhasMisMot3[,1]

#location recombinations
recomb3d<-lapply(recomb3b, function(x){
as.data.frame(t(as.data.frame(x)))})
rec_per_core<-colSums(rbindlist(recomb3d, fill=T), na.rm=T)
 names(rec_per_core)<-c(1:length(rec_per_core))
rec_per_core<-as.data.frame(rec_per_core)
rec_per_core<-rec_per_core/nrow(recomb_Mot)*100

#location corrected for uninformative markers
inf_list<-list()
for (k in 1:(ncol(RecMot2)/cl+1)){
inf_list[[k]]<-table(factor(RecMot2[,(1+seventy[k]):(seventy[k+1])],lev=c(0,1,2,5,9)))}
##estimate fraction of informative markers per core for each individual
inf_list2<-lapply(inf_list, function(x){
(x[[2]]+x[[3]])/sum(x)})
inf_vec<-unlist(inf_list2)
## corrected location recombinations
rec_per_core2<-rec_per_core/inf_vec
rec_per_core2<-as.data.frame(rec_per_core2)
rec_per_core2<-rec_per_core2/nrow(recomb_Mot)*100

#location recombinations corrected for core length (Mb)
##core length in bp
core_length<-numeric((ncol(RecMot2)/cl+1))
for (j in 1:(ncol(RecMot2)/cl+1)){
for (i in seventy[j+1]){
core_length[j]<-map[i,3]-map[seventy[j]+1,3]}}
##location recombinations corrected for core length (Mb)
rec_per_Mb<-rec_per_core/core_length*1000000
rec_per_Mb_mot<-as.data.frame(rec_per_Mb)
rec_per_Mb2<-rec_per_core2/core_length*1000000
rec_per_Mb_mot2<-as.data.frame(rec_per_Mb2)

#middle of cores
core_middle<-numeric((ncol(RecMot2)/cl+1))
for (j in 1:(ncol(RecFat2)/cl+1)){
for(i in (seventy[j]+seventy[j+1])/2){
core_middle[j]<-map[i,3]}}

#recombinations per Mother
recomb_Mot<-as.data.frame(recomb3c)
recomb_Mot$self<-as.numeric(rownames(recomb_Mot))
recomb_Mot<-merge(recomb_Mot, ped[,1:2], by='self')
recomb_Mot$V1<-as.numeric(recomb_Mot$V1)
recomb_Mot2<-recomb_Mot[,-1]
colnames(recomb_Mot)<-c('individual','recombinations','Mother')
meansum <- recomb_Mot %>% group_by(mother) %>% summarise(meansum = mean(recombinations)) 
meansum<- as.data.frame(meansum)

###############figures#################
pdf('Output_files/W1WC/WC/hist_rec_per_individual.pdf')
par(mfrow=c(1,1))
hist(as.vector(recomb3c[,1],mode="numeric"),xlab="recombinations", main='W1 recombinations per individual',right=FALSE,cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
dev.off()

pdf('Output_files/W1WC/WC/hist_rec_per_mother.pdf')
par(mfrow=c(1,1))
hist(as.vector(meansum[,2],mode='numeric'),xlab='recombinations', main='W1 recombinations per father',cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
dev.off()

pdf('Output_files/W1WC/rec_per_Mb_both_parents.pdf')
par(mfrow=c(3,1))
plot(core_middle/1000000, rec_per_Mb_fat[,1],main= 'recombination rate per Mb', xlab='Mb',ylab='rec. rate (cM/Mb)',type='n',cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
lines(core_middle/1000000, rec_per_Mb_fat[,1],type='l',col=1)
lines(core_middle/1000000, rec_per_Mb_mot[,1],type='l',col=3)
dev.off()

pdf('Output_files/W1WC/corrected_rec_per_Mb.pdf')
par(mfrow=c(3,1))
plot(core_middle/1000000, rec_per_Mb_fat2[,1],main= 'corrected recombination rate per Mb', xlab='Mb',ylab='W1 corrected recombination rate (cM/Mb)',type='n',cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
lines(core_middle/1000000, rec_per_Mb_fat2[,1],type='l',col=1)
lines(core_middle/1000000, rec_per_Mb_mot2[,1],type='l',col=4)
dev.off()

pdf('Output_files/W1WC/W1/inf_markers_vs_rec.pdf')
plot(inf_vec, rec_per_Mb[,1], xlab='fraction of informative markers', ylab='recombination rate per Mb',cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
dev.off()

###########tables###########
##recombinations
capture.output(tor2, file='Output_files/W1WC/WC/tor2.txt')
capture.output(recomb1b, file='Output_files/W1WC/WC/recomb1b.txt')
write.table(recomb1c, 'Output_files/W1WC/WC/recomb1c.txt')
capture.output(tor3b, file='Output_files/W1WC/WC/tor3b.txt')
capture.output(recomb3b, file='Output_files/W1WC/WC/recomb3b.txt')
write.table(recomb3c, file='Output_files/W1WC/WC/recomb3c.txt')
###per mother
write.table(recomb_mot, 'Output_files/W1WC/WC/recomb_mot.txt')
write.table(meansum, 'Output_files/W1WC/WC/meansum.txt')
##locations
write.table(rec_per_core, 'Output_files/W1WC/WC/rec_per_core.txt')
write.table(rec_per_core2, 'Output_files/W1WC/WC/rec_per_core2.txt')
write.table(rec_per_Mb, 'Output_files/W1WC/WC/rec_per_Mb.txt')
write.table(rec_per_Mb2, 'Output_files/W1WC/WC/rec_per_Mb2.txt')
capture.output(inf_vec, file='Output_files/W1WC/WC/inf_vec.txt')