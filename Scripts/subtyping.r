##Generate a MDS subtype table from output from filtering.r

subtype.mat<-matrix(ncol = 7)
subtype.mat<-subtype.mat[-1,]
colnames(subtype.mat) <- c("RA","RT","RARS","RCMD","RCMD-RS","RAEB","-5q")
for(i in COSMIC_comut[,8]){
  ifelse(i == "refractory_anaemia" && nchar(i) == 18,subtype.mat<-rbind(subtype.mat,c(1,0,0,0,0,0,0)),
    ifelse (i == "refractory_thrombocytopenia",subtype.mat<-rbind(subtype.mat,c(0,1,0,0,0,0,0)),
      ifelse (i == "refractory_anaemia_with_ringed_sideroblasts",subtype.mat<-rbind(subtype.mat,c(0,0,1,0,0,0,0)),
        ifelse (i == "refractory_cytopenia_with_multilinear_dysplasia" && nchar(i) == 47,subtype.mat<-rbind(subtype.mat,c(0,0,0,1,0,0,0)),
          ifelse (i == "refractory_cytopenia_with_multilinear_dysplasia_with_ringed_sideroblasts",subtype.mat<-rbind(subtype.mat,c(0,0,0,0,1,0,0)),
            ifelse (i == "refractory_anaemia_with_excess_blasts",subtype.mat<-rbind(subtype.mat,c(0,0,0,0,0,1,0)),
              ifelse (i == "5q-myelodysplastic_syndrome",subtype.mat<-rbind(subtype.mat,c(0,0,0,0,0,0,1)), subtype.mat<-rbind(subtype.mat,c(NA,NA,NA,NA,NA,NA,NA))
                    )
                  )
                )
              )
            )
          )
        )
      }
COSMIC_expand<-cbind(COSMIC_comut,subtype.mat)

#"Gene x gene and gene x condition"

gene.mat <- matrix(0,ncol = 3); colnames(gene.mat) <- c("Comparison","pval","OR")
gene.mat <- gene.mat[-1,]
for (i in seq(from = 8,to = 57)) {
     for (j in seq(from = 8,to = 57)) {
         Tab<-table(COSMIC_expand[,i],COSMIC_expand[,j]);
         fish<-fisher.test(Tab);
         fishnames<-paste(colnames(COSMIC_expand[i]),(colnames(COSMIC_expand[j])));
         fishpval <- fish$p.value;
         fishOR <- fish[["estimate"]][["odds ratio"]];
         gene.mat <- rbind(gene.mat, c(fishnames, fishpval, fishOR))
     }
     for (k in seq(from = 58, to = 64)) {
         Tab<-table(COSMIC_expand[,i],COSMIC_expand[,k]);
         fish<-fisher.test(Tab);
         fishnames<-paste(colnames(COSMIC_expand[i]),(colnames(COSMIC_expand[k])));
         fishpval <- fish$p.value;
         fishOR <- fish[["estimate"]][["odds ratio"]];
         gene.mat <- rbind(gene.mat, c(fishnames, fishpval, fishOR))
         }
 }

write.table(gene.mat,"COSMIC_Fish.txt",sep="\t",row.names=FALSE)


#"Gene x gene and gene x condition"

gene.mat <- matrix(0,ncol = 3); colnames(gene.mat) <- c("Comparison","pval","OR")
gene.mat <- gene.mat[-1,]
for (i in seq(from = 5,to = 54)) {
     for (j in seq(from = 5,to = 54)) {
         Tab<-table(COSMIC_single[,i],COSMIC_single[,j]);
         fish<-fisher.test(Tab);
         fishnames<-paste(colnames(COSMIC_single[i]),(colnames(COSMIC_single[j])));
         fishpval <- fish$p.value;
         fishOR <- fish[["estimate"]][["odds ratio"]];
         gene.mat <- rbind(gene.mat, c(fishnames, fishpval, fishOR))
     }
     for (k in seq(from = 55, to = 61)) {
         Tab<-table(COSMIC_single[,i],COSMIC_single[,k]);
         fish<-fisher.test(Tab);
         fishnames<-paste(colnames(COSMIC_single[i]),(colnames(COSMIC_single[k])));
         fishpval <- fish$p.value;
         fishOR <- fish[["estimate"]][["odds ratio"]];
         gene.mat <- rbind(gene.mat, c(fishnames, fishpval, fishOR))
         }
 }

write.table(gene.mat,"COSMIC_single_Fish.txt",sep="\t",row.names=FALSE)

# Cytogenetics
#Generate a MDS subtype table

subtype.mat<-matrix(ncol = 7)
subtype.mat<-subtype.mat[-1,]
colnames(subtype.mat) <- c("RA","RT","RARS","RCMD","RCMD-RS","RAEB","-5q")
for(i in cosmic_cyto_counts[,2]){
  ifelse(i == "refractory_anaemia" && nchar(i) == 18,subtype.mat<-rbind(subtype.mat,c(1,0,0,0,0,0,0)),
    ifelse (i == "refractory_thrombocytopenia",subtype.mat<-rbind(subtype.mat,c(0,1,0,0,0,0,0)),
      ifelse (i == "refractory_anaemia_with_ringed_sideroblasts",subtype.mat<-rbind(subtype.mat,c(0,0,1,0,0,0,0)),
        ifelse (i == "refractory_cytopenia_with_multilinear_dysplasia" && nchar(i) == 47,subtype.mat<-rbind(subtype.mat,c(0,0,0,1,0,0,0)),
          ifelse (i == "refractory_cytopenia_with_multilinear_dysplasia_with_ringed_sideroblasts",subtype.mat<-rbind(subtype.mat,c(0,0,0,0,1,0,0)),
            ifelse (i == "refractory_anaemia_with_excess_blasts",subtype.mat<-rbind(subtype.mat,c(0,0,0,0,0,1,0)),
              ifelse (i == "5q-myelodysplastic_syndrome",subtype.mat<-rbind(subtype.mat,c(0,0,0,0,0,0,1)), subtype.mat<-rbind(subtype.mat,c(NA,NA,NA,NA,NA,NA,NA))
                    )
                  )
                )
              )
            )
          )
        )
      }
cosmic_cyto_counts<-cbind(cosmic_cyto_counts,subtype.mat)


cyto.mat <- matrix(0,ncol = 3); colnames(cyto.mat) <- c("Comparison","pval","OR")
cyto.mat <- cyto.mat[-1,]
for (i in seq(from = 5,to = 15)) {
     for (j in seq(from = 17,to = 23)) {
         Tab<-table(cosmic_cyto_counts[,i],factor(cosmic_cyto_counts[,j], levels = c(0,1))); # had to coerce 2 levels when one column had no 1s
         fish<-fisher.test(Tab);
         fishnames<-paste(colnames(cosmic_cyto_counts[i]),(colnames(cosmic_cyto_counts[j])));
         fishpval <- fish$p.value;
         fishOR <- fish[["estimate"]][["odds ratio"]];
         cyto.mat <- rbind(cyto.mat, c(fishnames, fishpval, fishOR))
     }
}

# Cytogenetics vs Variants
#Generate a MDS subtype table

cyto.mat <- matrix(0,ncol = 3); colnames(cyto.mat) <- c("Comparison","pval","OR")
cyto.mat <- cyto.mat[-1,]
for (i in seq(from = 8,to = 57)) {
     for (j in seq(from = 68,to = 78)) {
         Tab<-table(COSMIC_mut_cyto[,i],factor(COSMIC_mut_cyto[,j], levels = c(0,1))); # had to coerce 2 levels when one column had no 1s
         fish<-fisher.test(Tab);
         fishnames<-paste(colnames(COSMIC_mut_cyto[i]),(colnames(COSMIC_mut_cyto[j])));
         fishpval <- fish$p.value;
         fishOR <- fish[["estimate"]][["odds ratio"]];
         cyto.mat <- rbind(cyto.mat, c(fishnames, fishpval, fishOR))
     }
}
