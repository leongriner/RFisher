#Script for filtering relevant mutations from COSMIC list of mutations haematopoietic tissue.
#Filtering is necessary as there are a number of publications that only examine a few genes (or only 1 gene).
#This will cause problems for determining associations.

#Requirements:
#plyr
#dplyr

#import dataset
V85_38_MUTANT <- read.csv("C:/Users/lgri018/Downloads/V85_38_MUTANT.csv") #321335 R, 35 C

#order by PMID. Ordering just for inspection purposes. Some publications have only 1 gene reported.
PMID_ord<-V85_38_MUTANT[order(V85_38_MUTANT$PUBMED_PMID),]

library(plyr)

#Get number of times this gene/publication combo appears, a bit slow
num_pmid_gene<- ddply(PMID_ord,.(GENE_NAME,PUBMED_PMID),'nrow')

#Get the number of times this publication appears
num_pmid_only<- ddply(PMID_ord,.(PUBMED_PMID),'nrow')

#Stick these together.
joined_pmid_gene<-join(num_pmid_gene,num_pmid_only, by="PUBMED_PMID")

MUT_PUB_RATIO<-joined_pmid_gene[,3]/joined_pmid_gene[,4]
joined_counted<-cbind(MUT_PUB_RATIO,joined_pmid_gene)

#Join by both PMID and Gene name to original file.
library(dplyr)
megajoin<-left_join(PMID_ord,joined_counted, by=c("PUBMED_PMID","GENE_NAME"))

#Filter by MDS, not NS subtype, not neutral FATHMM, Ratio of genes per publication not >0.9.
stage1_filtered<-data.frame();
stage1_filtered<-subset(megajoin,HISTOLOGY_SUBTYPE_1=="myelodysplastic_syndrome" & HISTOLOGY_SUBTYPE_2!="NS" & FATHMM_PREDICTION!="NEUTRAL" & MUT_PUB_RATIO<0.9, select=c(1,2,5:7,12:22,26,28,30,31,36))
rm(megajoin,PMID_ord,joined_counted)

#make a list of subtypes of interest
subtype_list<-c("5q-myelodysplastic_syndrome", "refractory_anaemia", "refractory_anaemia_with_excess_blasts", "refractory_anaemia_with_ringed_sideroblasts", "refractory_cytopenia_with_multilinear_dysplasia", "refractory_cytopenia_with_multilinear_dysplasia_with_ringed_sideroblasts","refractory_cytopenia_with_unilineage_dysplasia", "refractory_thrombocytopenia")

#Get subset from list
stage1_filtered<-stage1_filtered[stage1_filtered$HISTOLOGY_SUBTYPE_2 %in% subtype_list,]

#List of indexes that are rcud
rcudlist<-data.frame(); for (i in 1:dim(stage1_filtered)[1]){if (stage1_filtered[i,8]=="refractory_cytopenia_with_unilineage_dysplasia") rcudlist<-rbind(rcudlist,i)}

#change rcud to RA/RT etc
stage1_filtered[,c(8,9)]<-sapply(stage1_filtered[,c(8,9)],as.character)
for (j in rcudlist){stage1_filtered[j,8]<-stage1_filtered[j,9]}
rm(rcudlist)

#Filter out NS for this again.
COSMIC_comut<-subset(stage1_filtered,HISTOLOGY_SUBTYPE_2!="NS",select=c(1:21))


#for (j in subtype_list){
# namehold<-sprintf("%s_mat",j)
#  assign(namehold[1],matrix(data = 0, nrow = dim(single_filtered)[1]));
#  for (i in 1:dim(stage1_filtered)[1]){
#   if (stage1_filtered$HISTOLOGY_SUBTYPE_2[i]=="j"){`paste(namehold)`[i,1]<-1}
#   }
#  }
