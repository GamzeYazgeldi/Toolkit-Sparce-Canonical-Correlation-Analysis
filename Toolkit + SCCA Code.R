Sys.setenv(LANG="en")
# PART I: DRUG-PROTEIN ANALYSIS TOOLKIT
setwd("C:/Users/Gamze Yazgeldi/Desktop//CanonicalCorrelationAnalysis")
list.files();

# These files were retrieved from DrugBank without any changing.
drug_links <- read.csv2("drug links.csv", header=TRUE, sep=",");
target_links <-read.csv2("target links.csv",header=TRUE,sep=",");
# For seperation DrugIDs column by ;
library(tidyverse)
target_links <- target_links %>% separate_rows(Drug.IDs)


# Drug Protein Analysis Toolkit for finding host-based and pathogen-based drugs.
# Also, for comparing drug targets and non-drug targets by topological parameters analysis.
# Interested pathogen and its host interspecies PPI are retrived from PHISTOTEST and is uploaded to here.
# Host based approach drugs were detected. (Or pathogen based approach where in other code script)
# Degree and betweeness centrality values were also retrieved from PHISTOTEST througly Graph Analysis option.
func_hostbased <- function(phi_data, ga_data){
  phi_data <- read.csv2(file.choose(),header=TRUE,sep=',')     # to choose interested phi data
  ga_data <- read.csv2(file.choose(),header=TRUE,sep=',')   # to choose interested degree & btw centrality data
  
  phi_data <- phi_data[!duplicated(phi_data[,"Host.Uniprot.ID"]),]
  
  pro_drug <- merge(phi_data, target_links, by.x =c("Host.Uniprot.ID"), by.y=c("UniProt.ID"), all.x = FALSE, all.y = FALSE)
  length(unique(pro_drug$Drug.IDs))
  
  # Proteins that have not topological values are ignored. Because these values are important for network analysis
  pro_top <- merge(ga_data, pro_drug, by="Host.Uniprot.ID",all = FALSE)
  
  # Other drug IDs in different sources are retrieved.
  drug_other <- merge(pro_top,drug_links,by.x = "Drug.IDs",by.y="DrugBank.ID",all=FALSE)
  
  drug_analysis <- data.frame(host_protein_ids=drug_other$Host.Uniprot.ID, host_protein_name=drug_other$Name.x,
                              pdb_id=drug_other$PDB.ID, species=drug_other$Species, drug_ids=drug_other$Drug.IDs,
                              drug_name=drug_other$Name.y, drug_type=drug_other$Drug.Type, kegg_id=drug_other$KEGG.Drug.ID,
                              chebi_id=drug_other$ChEBI.ID, degree=drug_other$Degree, 
                              betweenness_centrality=drug_other$Betweenness.Centrality)
}
result_hostbased <- func_hostbased(phi_data,ga_data);
write.csv2(result_hostbased, file="Host Based Drug Analysis Result.csv")


# PART II: CONSTRUCTING INPUT PROFILES THAT ARE DRUG MATRIX AND PROCESS/PATHWAY MATRIX.
# Firstly, detection common have both knowledge of drug and biological pathway:
result_hostbased <- read.csv2("Host Based Drug Analysis Result.csv", header=TRUE, sep=";");
#pathway <- read.csv2("t_pathway_pc2_testlive.csv",header=TRUE, sep = ";")  #it was so low data

# Downloaded from reactome, pathway table has 881.427 rows and 6 columns
pathway <- read.csv2("UniProt2Reactome_All_Levels.txt",header=FALSE,sep="\t");  #new pathway table
# 5th column is IEA and TAS (other organisms and human?)
# Extracted only human pathways from this.
pathway <- pathway[,-5];
pathway <- pathway[which(pathway$V6=='Homo sapiens'),]; #for only human pathway table is 134.440 x 5

# Detection common proteins in order to find correlation based (common) proteins
commonproteins <- as.data.frame(intersect(result_hostbased$host_protein_ids,pathway$V1));

# Number of associated proteins
drugass <- result_hostbased[which(result_hostbased$host_protein_ids %in% commonproteins[,1]), c(1,2,6)]
pathwayass <- pathway[which(pathway$V1 %in% commonproteins[,1]), c(1,2,4)]


# Goal is to find targeting drug. If targets relevant protein, put 1, otherwise 0.
drugs <- as.data.frame(unique(drugass$drug_ids))
proteindrugmat <- matrix(data=0, nrow=nrow(commonproteins), ncol=nrow(drugs));
rownames(proteindrugmat) <- commonproteins[,1];
colnames(proteindrugmat) <- drugs[,];
drug_list <-  list();
for (i in 1:nrow(commonproteins)){   #for finding drugs (list) that targets commponproteins
  drug_list[[i]] <-  t(as.data.frame(result_hostbased[which(result_hostbased$host_protein_ids %in% commonproteins[i,]), 6]));
  for (j in 1:length(drug_list[[i]])){   #componentteki her bir ilacýn drugs'daki yerini bulup 1 yazýlmasý için loop.
    d_index <- which(drugs[,1] %in% drug_list[[i]][j]);
    proteindrugmat[i,d_index] <- 1;
  }
}
dim(proteindrugmat)
write.csv2(proteindrugmat, file="Protein-Drug Sparse Matrix.csv")


pathways <- as.data.frame(unique(pathwayass$V4))
proteinpathwaymat <- matrix(data=0, nrow=nrow(commonproteins), ncol=nrow(pathways));
rownames(proteinpathwaymat) <- commonproteins[,1];
colnames(proteinpathwaymat) <- pathways[,];
pathway_list <- list();
for (i in 1:nrow(commonproteins)){ #for finding pathways (list) that belongs to commponproteins
  pathway_list[[i]] <-  t(as.data.frame(pathway[which(pathway$V1 %in% commonproteins[i,]), 4]));
  for (j in 1:length(pathway_list[[i]])){
    p_index <- which(pathways[,1] %in% pathway_list[[i]][j]);
    proteinpathwaymat[i,p_index] <- 1;
  }                              
}
dim(proteinpathwaymat)
write.csv2(proteinpathwaymat, file="Protein-Pathway Sparse Matrix.csv")

# PART III: CREATE FUNCTIONS OF PROPOSED METHOD 
library(ROCR)
library(PMA)
#library(kernlab) is not necessary for scca method
library(CCA)

scca <- function (X, Y, c1, c2, ncomp=4) 
{
  Xraw <- X; Yraw <- Y
  colindex.nonzero.x <- (1:ncol(X))[apply(X,2,sd)!=0]
  colindex.nonzero.y <- (1:ncol(Y))[apply(Y,2,sd)!=0]
  X <- X[,colindex.nonzero.x]; Y <- Y[,colindex.nonzero.y]
  #X <- scale(X); Y <- scale(Y)
  nx <- nrow(X);  px <- ncol(X)
  ny <- nrow(Y);  py <- ncol(Y)
  n  <- nx
  # canonical correlation
  dvec <- rep(0, ncomp)                       # zero vector in ncomp size
  # weight matrix features of non-zeros
  #umat <- matrix(0, px, ncomp)
  #vmat <- matrix(0, py, ncomp)
  # weight matrix for all input features
  umatall <- matrix(0, ncol(Xraw), ncomp)
  vmatall <- matrix(0, ncol(Yraw), ncomp)
  result <- CCA(x=X, z=Y, typex = "standard", typez = "standard", penaltyx = c1, penaltyz = c2, K = ncomp, niter = 15, v = NULL, trace = FALSE, standardize = FALSE, xnames = NULL, znames = NULL, chromx = NULL, chromz = NULL, upos = FALSE, uneg = FALSE, vpos = FALSE, vneg = FALSE, outcome = NULL, y = NULL, cens = NULL) 
  dvec <- result$d
  # adjustment of positive and negatieve
  for (k in 1:ncomp){
    wmaxorderx <- order(abs(result$u[,k]), decreasing=T)[1]
    if (result$u[wmaxorderx,k] < 0){
      result$u[,k] <- - result$u[,k]
      result$v[,k] <- - result$v[,k]
    }
  }
  umatall[colindex.nonzero.x, ] <- result$u
  vmatall[colindex.nonzero.y, ] <- result$v
  rownames(umatall) <- colnames(Xraw)
  rownames(vmatall) <- colnames(Yraw)
  scorex <- X %*% result$u
  scorey <- Y %*% result$v
  rownames(scorex) <- rownames(Xraw)
  rownames(scorey) <- rownames(Yraw)
  corvec <- rep(0,ncomp)
  for (j in 1:ncomp){
    corvec[j] <- cor(scorex[,j],scorey[,j])
  } 
  #list(u=umatall, v=vmatall, d=dvec)
  list(u=umatall, v=vmatall, d=dvec, scorex=scorex, scorey=scorey, rho=corvec)
}


# PART IV: FIND CANONICAL COMPONENTS AND CREATE NETWORK REPRESENTATION
# OCCA
X=scale(proteindrugmat) #data standardization. #gold standard data set #its not necessary
X[is.na(X)] <- 0 #IS IT TRUE? NaN values is replaced by 0.
Y=scale(proteinpathwaymat) #also data transformation technique, used in regression and clustering analysis, as here.
Y[is.na(Y)] <- 0
#protein.occa <- scca(X=proteindrugmat, Y=proteinpathwaymat, c1=1, c2=1, ncomp=10) 

# SCCA
protein.scca <- scca(X, Y, c1=0.1, c2=0.1, ncomp=20) 
ncomp <- 20
u <- rownames(protein.scca[["u"]]) #each weights for each drugs;
u <- cbind(u,protein.scca[["u"]])
v <- rownames(protein.scca[["v"]]) #each weights for each pathways
v <- cbind(v,protein.scca[["v"]])

# Number of associated proteins
drugass <- result_hostbased[which(result_hostbased$host_protein_ids %in% commonproteins[,1]), c(1,2,6)]
pathwayass <- pathway[which(pathway$V1 %in% commonproteins[,1]), c(1,2,4)]

# CONSTRUCTING A NETWORK with Procedure I
# According to protein.scca results, u and v linear combinations/canonical components were detected.
# Highly scoring objects were extracted for network representation.
# FOR CYTOSCAPE INPUT:
library(dplyr)
liste <- list()
for(i in 1:ncomp){ #for listing all list
  for (i in 1:ncomp){  #for making list for each component
    # Component n: positive weights on targeting drugs
    # componentlerin içerisindeki weightleri sýfýr olan druglarý bulunuyor
    nthcompd <- protein.scca[["u"]][which(abs(protein.scca[["u"]][,i]) !=0),i];
    a <- rownames(as.data.frame(nthcompd)); #each weights for each drugs
    nthcompd <- cbind(a,nthcompd); #row names are present in one column anymore
    # bulunan druglarý drugass(ilgili ortak proteinlerin ilaçlarýyla olan ikili tablo)de arýyorum
    commonnd <- intersect(nthcompd[,1],drugass[,3]);
    # bu druglarýn hedeflediði proteinlerini çekiyorum
    commonnd<- drugass[which(drugass[,3] %in% commonnd),c(2,3)];
    
    # positive weights on belonging pathways
    nthcompp <- protein.scca[["v"]][which(abs(protein.scca[["v"]][,i]) !=0),i];
    b <- rownames(as.data.frame(nthcompp)); #each weights for each pathways
    nthcompp <- cbind(b,nthcompp);
    commonnp <- intersect(nthcompp[,1],pathwayass[,3]);
    commonnp<- pathwayass[which(pathwayass[,3] %in% commonnp),c(1,3)];
    commonnp <- unique(commonnp)   #for making unique
    
    bipartiten <- merge.data.frame(commonnd,commonnp,by.x = "host_protein_ids", by.y ="V1");
    q <- pathway[which(pathway[,4] %in% intersect(bipartiten[,3],pathway[,4])),c(2,4)]; #add the pathway ids
    q <- unique(q)
    #bipartite10 <- cbind(q[,2], bipartite10)
    liste[[i]] <- merge.data.frame(q,bipartiten, by.x = "V4", by.y = "V4")
    liste[[i]] <- cbind(i, liste[[i]])
    colnames(liste[[i]]) <- c("component number","pathway name","pathway id","host protein ids","drugids")
    #rownames(liste[[i]]) <- NULL
  }
  return(liste[[i]])  
} 
all_bipartite_cytoscape <- do.call(rbind, liste)
write.csv2(all_bipartite_cytoscape,file="Human-AllPathogen_20Components-CYTOSCAPE_BINARY.csv")

# Network representation for each components:
# Saved each component in different .csv files.
for (i in 1:ncomp){
  a <- all_bipartite_cytoscape[which(all_bipartite_cytoscape$`component number`==i),]
  file=paste(i, "Component-Cytoscape.csv")
  write.csv2(a, file=file)
}

# Construction components table for both drug and pathway contents:
# Component n: positive weights on targeting drugs
library(plyr)
liste_d <- list()
for (i in 1:ncomp){
  for (i in 1:ncomp){
    nthcompd <- protein.scca[["u"]][which(abs(protein.scca[["u"]][,i]) !=0),i] #each component drug contents
    a <- rownames(as.data.frame(nthcompd)); 
    nthcompd <- cbind(a,nthcompd); #row names are present in one column anymore
    commonnd <- intersect(nthcompd[,1],drugass[,3]);
    commonnd<- drugass[which(drugass[,3] %in% commonnd),c(2,3)]; 
    count_d <- count(commonnd[,2])
    
    weight_d <- as.data.frame(nthcompd[,2])
    a <- rownames(nthcompd); 
    weight_d <- cbind(a,weight_d);
    
    b <- count_d[order(match(weight_d$a,count_d)),2];
    c <- weight_d[order(match(weight_d$a,count_d)),c(2,1)]
    liste_d[[i]] <- data.frame(b,c)
    liste_d[[i]] <- cbind(i, liste_d[[i]])
    colnames(liste_d[[i]]) <- c("component number","number of associated protein",
                                "weights","drugids")
  }  
  return(liste_d[[i]])
} 
all_bipartite_drug <- do.call(rbind, liste_d)
write.csv2(all_bipartite_drug,file="Human-AllPathogen_20Components-DRUG.csv")


# Component n: positive weights on targeting pathways
liste_p <- list()
for (i in 1:ncomp){
  nthcompp <- protein.scca[["v"]][which(abs(protein.scca[["v"]][,i]) !=0),i]; #each component pathway contents
  b <- rownames(as.data.frame(nthcompp)); #each weights for each pathways
  nthcompp <- cbind(b,nthcompp); #row names are present in one column anymore
  commonnp <- intersect(nthcompp[,1],pathwayass[,3]);
  commonnp<- pathwayass[which(pathwayass[,3] %in% commonnp),c(1,3)];
  commonnp <- unique(commonnp)   #for making unique
  count_p <- count(commonnp[,2])
  
  weight_p <- as.data.frame(nthcompp[,2])
  a <- rownames(as.data.frame(weight_p)); 
  weight_p <- cbind(a,weight_p);
  
  liste_p[[i]] <- merge(count_p,weight_p,by.x="x",by.y = "a")
  
  t <- match(liste_p[[i]][,1],pathwayass[,3])
  t3 <- as.data.frame(unique(pathwayass[t,2]))
  
  liste_p[[i]] <- cbind(t3,liste_p[[i]])
  
  liste_p[[i]] <- cbind(i, liste_p[[i]])
  liste_p[[i]] <- liste_p[[i]][,c(1,4,5,2,3)]
  colnames(liste_p[[i]]) <- c("component number","number of associated proteins","weight",
                              "pathwayids","pathway names")
}
all_bipartite_pathway <- do.call(rbind, liste_p)
write.csv2(all_bipartite_pathway,file="Human-AllPathogen_20Components-PATHWAY.csv")


# Component n: high scoring proteins for both targeting drugs and belongs pathways
library(tidyverse)
liste_protein <- list()
for(i in 1:ncomp){   #loop for creating liste_protein
  for (i in 1:ncomp){  #loop for each proteins in the components
    e <- unique(liste[[i]][,4])  # unique proteins in the components
    
    common_pro_path <- intersect(e[1:length(e)],pathwayass[,1]);   #pathwayassde o proteinleri buluyor
    # bu proteinlere ait pathwayler bulunacak
    common_pro_path <- pathwayass[which(pathwayass[,1] %in% common_pro_path),c(1,3)];
    count_associated_pathways <- count(common_pro_path[,1])  #bu protein pathwayass tablosunda kaç pathway ile tanýmlanmýþ
    
    common_pro_drug <- intersect(e[1:length(e)],drugass[,2]); # bu proteini drugassde bul
    common_pro_drug <- drugass[which(drugass[,2] %in% common_pro_drug),c(2,3)]; #bulunca protein idsi ve ilgili ilacý al
    #kaç ilaçla etkileþmiþ onu count et
    count_associated_drugs <- count(common_pro_drug[,1]) #bu protein drugass tablosunda kaç drug ile tanýmlanmýþ
    
    liste_protein[[i]] <- data.frame(i,"",count_associated_pathways[,2],count_associated_drugs[,2],e[1:length(e)])
    colnames(liste_protein[[i]]) <- c("component number","component_score","number of associated pathways",
                                      "number of associated drugs","protein ids")
  }
  return(liste_protein[[i]])  
} 
all_bipartite_protein <- do.call(rbind, liste_protein)
write.csv2(all_bipartite_protein,file="Human-AllPathogen_20Components-PROTEIN.csv")



# PART V: INTERPRETATION OF RESULTS
# rho is correlation coefficient of each drug's varible set and pathway's variable set component 
# 0 < rho < 1 
# Drug ve pathway variable setleri arasýndaki iliþkiyi maksimum gösteren component:
max(protein.scca[["rho"]])
max_corr_component <- which(protein.scca[["rho"]]==max(protein.scca[["rho"]])) 
#8.ci component aralarýndaki iliþkiyi en iyi (maximum) gösteren component.

length(which(all_bipartite_cytoscape$`component number`=="1"))  #29 tane pathway var; bu componentte
length(which(all_bipartite_drug$`component number`=="2"))   #38 tane ilaç var; bu componentte

#2.ci componentteki protein:P07900
length(which(pathway$V1=="P07900"))  #bu proteine ait 101 tane pathway bilgisi var
length(which(drugass$host_protein_ids=="P07900"))  #bu proteine ait 49 tane ilaç var.

#yani componentte çýkan proteine göre yapýlmamýþ o ilaç ve pathway. bu okey.

length(unique(all_bipartite_cytoscape$`pathway name`))
length(unique(all_bipartite_cytoscape$`host protein ids`)) #956 tane çýktý.


#her componentteki ilaç ve pathwayler 1den fazla proteinle iliþkili.
b <- unique(all_bipartite_cytoscape[which(all_bipartite_cytoscape$`component number`=="1"),4])
View(b)

# TO FIND DRUG NAME, TARGETING PROTEIN NAME OF MAX CORRELATED COMPONENT.
# For Herpesviridae:
x <- all_bipartite_cytoscape[all_bipartite_cytoscape$`component number`==14,]
herpes_table <- result_hostbased[which(result_hostbased$drug_ids %in% x$drugids), c(2,3,6,7)]
View(unique(herpes_table$drug_name))
write.csv2(unique(herpes_table$drug_name),file="herpes-drug.csv")

y <-  all_bipartite_cytoscape[all_bipartite_cytoscape$`component number`==14,]
write.csv2(unique(y$`pathway name`),file="herpes-pathwayname.csv")

# For Orthomyxoviridae:
x <- all_bipartite_cytoscape[all_bipartite_cytoscape$`component number`==4,]
ortho_table <- result_hostbased[which(result_hostbased$drug_ids %in% x$drugids), c(2,3,6,7)]
View(unique(ortho_table$drug_name))
write.csv2(unique(ortho_table$drug_name),file="ortho-drug.csv")

y <-  all_bipartite_cytoscape[all_bipartite_cytoscape$`component number`==4,]
write.csv2(unique(y$`pathway name`),file="ortho-pathwayname.csv")
