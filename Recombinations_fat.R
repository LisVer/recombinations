#load data
PhasMisFat3<-read.table('Input_files/PhasMisFat3.txt') ##output file Origin_alleles_fat.R
RecFat2<-as.matrix(PhasMisFat3)
map<-read.table("Input_files/Chr1_map2.txt", header=TRUE) ##output file Phasing_mistakes.R
ped<-read.table("Input_files/ped5.txt") ##output file Phasing_mistakes.R

#load packages
library(data.table)
library(plyr)
library(dplyr)

#prepare data
map<-map[,-3]
colnames(map)<-c('chr','marker','pos')
RecFat2<-RecFat2[,-1]
cl<-70 ##core_length

#count subsequent snps with the same origin
seventy<-c(seq(0,ncol(RecFat2),by=cl),ncol(RecFat2))

tor<-list()
tor2<-list()
for (j in 1:nrow(RecFat2)){
vec<-RecFat2[j,]
for (k in 1:(ncol(RecFat2)/cl+1)){
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
rownames(recomb1c)<-PhasMisFat3[,1]

#only count phase shifts when there are at least 5 consecutive snps from the same parent haplotype
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
rownames(recomb3c)<-PhasMisFat3[,1]

#location recombinations
recomb3d<-lapply(recomb3b, function(x){
as.data.frame(t(as.data.frame(x)))})
rec_per_core<-colSums(rbindlist(recomb3d, fill=T), na.rm=T)
 names(rec_per_core)<-c(1:length(rec_per_core))
rec_per_core<-as.data.frame(rec_per_core)
rec_per_core<-rec_per_core/nrow(recomb_fat)*100

#location corrected for uninformative markers
inf_list<-list()
for (k in 1:(ncol(RecFat2)/cl+1)){
inf_list[[k]]<-table(factor(RecFat2[,(1+seventy[k]):(seventy[k+1])],lev=c(0,1,2,5,9)))}
##estimate fraction of informative markers per core for each individual
inf_list2<-lapply(inf_list, function(x){
(x[[2]]+x[[3]])/sum(x)})
inf_vec<-unlist(inf_list2)
## corrected location recombinations
rec_per_core2<-rec_per_core/inf_vec
rec_per_core2<-as.data.frame(rec_per_core2)
rec_per_core2<-rec_per_core2/nrow(recomb_fat)*100

#location recombinations corrected for core length (Mb)
##core length in bp
core_length<-numeric((ncol(RecFat2)/cl+1))
for (j in 1:(ncol(RecFat2)/cl+1)){
for (i in seventy[j+1]){
core_length[j]<-map[i,3]-map[seventy[j]+1,3]}}
##location recombinations corrected for core length (Mb)
rec_per_Mb<-rec_per_core/core_length*1000000
rec_per_Mb_fat<-as.data.frame(rec_per_Mb)
rec_per_Mb2<-rec_per_core2/core_length*1000000
rec_per_Mb_fat2<-as.data.frame(rec_per_Mb2)
#middle of cores
core_middle<-numeric((ncol(RecFat2)/cl+1))
for (j in 1:(ncol(RecFat2)/cl+1)){
for(i in (seventy[j]+seventy[j+1])/2){
core_middle[j]<-map[i,3]}}

#recombinations per father
recomb_fat<-as.data.frame(recomb3c)
recomb_fat$self<-as.numeric(rownames(recomb_fat))
recomb_fat<-merge(recomb_fat, ped[,1:2], by='self')
recomb_fat[,4]<-numeric(nrow(recomb_fat))
for(i in 1:nrow(recomb_fat)){
recomb_fat[i,4]<-length(which(recomb_fat[,1]==recomb_fat[i,1]))} #number of repeated observations per father
colnames(recomb_fat)<-c('father', 'individual','recombinations','repeats')
recomb_fat$recombinations<-as.numeric(recomb_fat$recombinations)
meansum <- recomb_fat %>% group_by(father) %>% summarise(meansum = mean(recombinations)) 
meansum<- as.data.frame(meansum)


###############figures#################
pdf('Output_files/W1WC/W1/hist_rec_per_individual.pdf')
par(mfrow=c(1,1))
hist(as.vector(recomb3c[,1],mode="numeric"),xlab="recombinations", main='W1 recombinations per individual',right=FALSE,cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
dev.off()

pdf('Output_files/W1WC/W1/hist_rec_per_father.pdf')
par(mfrow=c(1,1))
hist(as.vector(meansum[,2],mode='numeric'),xlab='recombinations', main='W1 recombinations per father',cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
dev.off()

pdf('Output_files/W1WC/W1/inf_markers_vs_rec.pdf')
plot(inf_vec, rec_per_Mb[,1], xlab='fraction of informative markers', ylab='recombination rate per Mb',cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
dev.off()

###########tables###########
##recombinations
capture.output(tor2, file='Output_files/W1WC/W1/tor2.txt')
capture.output(recomb1b, file='Output_files/W1WC/W1/recomb1b.txt')
write.table(recomb1c, 'Output_files/W1WC/W1/recomb1c.txt')
capture.output(tor3b, file='Output_files/W1WC/W1/tor3b.txt')
capture.output(recomb3b, file='Output_files/W1WC/W1/recomb3b.txt')
write.table(recomb3c, file='Output_files/W1WC/W1/recomb3c.txt')
###per father
write.table(recomb_fat, 'Output_files/W1WC/W1/recomb_mot.txt')
write.table(meansum, 'Output_files/W1WC/W1/meansum.txt')
##locations
write.table(rec_per_core, 'Output_files/W1WC/W1/rec_per_core.txt')
write.table(rec_per_core2, 'Output_files/W1WC/W1/rec_per_core2.txt')
write.table(rec_per_Mb, 'Output_files/W1WC/W1/rec_per_Mb.txt')
write.table(rec_per_Mb2, 'Output_files/W1WC/W1/rec_per_Mb2.txt')
capture.output(inf_vec, file='Output_files/W1WC/W1/inf_vec.txt')