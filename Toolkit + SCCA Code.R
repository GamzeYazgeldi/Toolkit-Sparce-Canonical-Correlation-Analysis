# TODO: To obtain drug information with using Drug-Protein Interaction Analysis.
# TODO: To use a machine learning method for obtaining correlated components.


# INPUT FILES: SARS-CoV-2-Human, t_analysis, UniProt2Reactome_All_Levels.

# About input files:
# Interested SARS-CoV-2-Human PHI data were retrieved from IntAct. 
# Pathway information of human proteins were retrieved from Reactome.


# HOST-BASED DRUGS:
func_hostbased <- function(phi_data, ga_data){
  phi_data <- read.csv2(file.choose(),header=TRUE,sep=';')     # to choose interested phi data
  ga_data <- read.csv2(file.choose(),header=FALSE,sep=',')   # to choose interested degree & btw centrality data
  
  phi_data <- phi_data[!duplicated(phi_data[,"Host.Uniprot.ID"]),]
  
  pro_drug <- merge(phi_data, target_links, by.x =c("Host.Uniprot.ID"), by.y=c("UniProt.ID"), all.x = FALSE, all.y = FALSE)
  length(unique(pro_drug$Drug.IDs))
  
  pro_top <- merge(ga_data, pro_drug, by.x="V1", by.y="Host.Uniprot.ID",all = FALSE)
  
  # Other drug IDs in different sources are retrieving.
  drug_other <- merge(pro_top,drug_links,by.x = "Drug.IDs",by.y="DrugBank.ID",all=FALSE)
  
  drug_analysis <- data.frame(host_protein_ids=drug_other$V1, host_protein_name=drug_other$Name.x,
                              pdb_id=drug_other$PDB.ID, species=drug_other$Species, drug_ids=drug_other$Drug.IDs,
                              drug_name=drug_other$Name.y, drug_type=drug_other$Drug.Type, kegg_id=drug_other$KEGG.Drug.ID,
                              chebi_id=drug_other$ChEBI.ID, degree=drug_other$V2, 
                              betweenness_centrality=drug_other$V3)
}
result_hostbased <- func_hostbased(phi_data,ga_data);
write.csv2(result_hostbased, file="Host Based Drug Analysis Result.csv")


# STEP I: Constructing Input Profiles that are Drug Matrix and Pathway Matrix for Sparse Canonical Correlation Analysis 
# Finding common have both knowledge of drug and biological pathway:
result_hostbased <- read.csv2("Host Based Drug Analysis Result.csv", header=TRUE, sep=";");
# Downloaded pathway information from Reactome, pathway table has 881.427 rows and 6 columns
pathway <- read.csv2("UniProt2Reactome_All_Levels.txt",header=FALSE,sep="\t");  
# Extracting only human pathways from this.
pathway <- pathway[,-5];
pathway <- pathway[which(pathway$V6=='Homo sapiens'),]; #for only human pathway table is 134.440 x 5
# Finding common proteins in order to find correlation based (common) proteins
commonproteins <- as.data.frame(intersect(result_hostbased$host_protein_ids,pathway$V1));
# Finding drug and pathway information of common proteins to create input matrix.
drugass <- result_hostbased[which(result_hostbased$host_protein_ids %in% commonproteins[,1]), c(1,2,6)]
pathwayass <- pathway[which(pathway$V1 %in% commonproteins[,1]), c(1,2,4)]


# If drug targets relevant protein, putting 1, otherwise 0.
drugs <- as.data.frame(unique(drugass$drug_ids))
proteindrugmat <- matrix(data=0, nrow=nrow(commonproteins), ncol=nrow(drugs));
rownames(proteindrugmat) <- commonproteins[,1];
colnames(proteindrugmat) <- drugs[,];
drug_list <-  list();
for (i in 1:nrow(commonproteins)){   #for finding drugs (list) that targets commponproteins
  drug_list[[i]] <-  t(as.data.frame(result_hostbased[which(result_hostbased$host_protein_ids %in% commonproteins[i,]), 6]));
  for (j in 1:length(drug_list[[i]])){   #componentteki her bir ilacın drugs'daki yerini bulup 1 yazılması için loop.
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


# STEP II: Using Functions of Sparse Canonical Correlation Method
# TODO: To find a relationship between pathway and drug based on proteins.
# NOTE: This method code was taken from Yamanishi Group. Parameters were set according to this study datasets.

library(ROCR) 
library(PMA) #penalized multivariate analysis library
library(CCA) #canonical correlation analysis library

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
  # weight matrix for all input features
  umatall <- matrix(0, ncol(Xraw), ncomp)
  vmatall <- matrix(0, ncol(Yraw), ncomp)
  result <- CCA(x=X, z=Y, typex = "standard", typez = "standard", penaltyx = c1, penaltyz = c2, K = ncomp, niter = 15, v = NULL, trace = FALSE, standardize = FALSE, xnames = NULL, znames = NULL, chromx = NULL, chromz = NULL, upos = FALSE, uneg = FALSE, vpos = FALSE, vneg = FALSE, outcome = NULL, y = NULL, cens = NULL) 
  dvec <- result$d
  # adjustment of positive and negative
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


# STEP III: Finding Canonical Components 
# TODO: To find correlated drug-pathway components. To find u and v linear combinations/canonical components.
# Making data standardization.
X=scale(proteindrugmat) 
X[is.na(X)] <- 0 
Y=scale(proteinpathwaymat) 
Y[is.na(Y)] <- 0

# Applying SCCA method:
protein.scca <- scca(X, Y, c1=0.1, c2=0.1, ncomp=20) 
ncomp <- 20
u <- rownames(protein.scca[["u"]]) #each weights for each drugs;
u <- cbind(u,protein.scca[["u"]])
v <- rownames(protein.scca[["v"]]) #each weights for each pathways
v <- cbind(v,protein.scca[["v"]])


# STEP IV: Constructing a Network for Network Representation in Cytoscape
# TODO: To find correlated drug-pathway components. 

# Creating a drug-pathway interaction network file:
library(dplyr)
liste <- list()
for (i in 1:ncomp){ # listing all list
  for (i in 1:ncomp){ #making list for each component
    # Finding non-zero elements in the each component.
    nthcompd <- protein.scca[["u"]][which(abs(protein.scca[["u"]][,i]) !=0),i]; # Add the row names as column.
    a <- rownames(as.data.frame(nthcompd));
    nthcompd <- cbind(a,nthcompd);
    # Finding common drugs from drugass table and drugs which have non-zero weights in the each components.
    commonnd <- intersect(nthcompd[,1],drugass[,3]);
    # Retrieving targeted proteins.
    commonnd<- drugass[which(drugass[,3] %in% commonnd),c(2,3)];
    # Finding non-zero elements in the each component.
    nthcompp <- protein.scca[["v"]][which(abs(protein.scca[["v"]][,i]) !=0),i]; # Adding the row names as column.
    b <- rownames(as.data.frame(nthcompp)); #each weights for each pathways
    nthcompp <- cbind(b,nthcompp); 
    # Finding common pathways from pathwayass table and pathways which have non-zero weights in the each components.
    commonnp <- intersect(nthcompp[,1],pathwayass[,3]);
    commonnp<- pathwayass[which(pathwayass[,3] %in% commonnp),c(1,3)];
    commonnp <- unique(commonnp) # making unique 
    # Merging drugs and pathways which have non-zero elements in the each components.
    bipartiten <- merge.data.frame(commonnd,commonnp,by.x = "host_protein_ids", by.y ="V1");
    q <- pathway[which(pathway[,4] %in% intersect(bipartiten[,3],pathway[,4])),c(2,4)]; #adding the pathway ids
    q <- unique(q)
    liste[[i]] <- merge.data.frame(q,bipartiten, by.x = "V4", by.y = "V4")
    liste[[i]] <- cbind(i, liste[[i]])
    colnames(liste[[i]]) <- c("component number","pathway name","pathway id","host protein ids","drugids")
    #rownames(liste[[i]]) <- NULL
    }
    return(liste[[i]])
  }
all_bipartite_cytoscape <- do.call(rbind, liste)
write.csv2(all_bipartite_cytoscape,file="20Components-CytoScape.csv")

# Constructing components table for both drug and pathway contents:
# Component n: positive weights on targeting drugs.
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
    
    b <- count_d[order(match(weight_d$a,count_d)),2]
    c <- weight_d[order(match(weight_d$a,count_d)),c(2,1)]
    liste_d[[i]] <- data.frame(b,c)
    liste_d[[i]] <- cbind(i, liste_d[[i]])
    colnames(liste_d[[i]]) <- c("component number","number of associated protein",
                                "weights","drugids")
  }  
  return(liste_d[[i]])
} 
all_bipartite_drug <- do.call(rbind, liste_d)
write.csv2(all_bipartite_drug,file="20Components-DRUG.csv")


# Making result of each components in the different excel sheets.
for (i in 1:ncomp){
  a <- all_bipartite_cytoscape[which(all_bipartite_cytoscape$`component number`==i),]
  file=paste(i, "Component-Cytoscape.csv")
  write.csv2(a, file=file)}


# Putting pathway information table for each component.
liste_p <- list()
for (i in 1:ncomp){
  nthcompp <- protein.scca[["v"]][which(abs(protein.scca[["v"]][,i]) !=0),i];
  b <- rownames(as.data.frame(nthcompp));
  nthcompp <- cbind(b,nthcompp);
  commonnp <- intersect(nthcompp[,1],pathwayass[,3]);
  commonnp<- pathwayass[which(pathwayass[,3] %in% commonnp),c(1,3)];
  commonnp <- unique(commonnp) #for making unique
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
  colnames(liste_p[[i]]) <- c("component number","number of associated proteins","weight", "pathwayids","pathway names")
  }
all_bipartite_pathway <- do.call(rbind, liste_p)
write.csv2(all_bipartite_pathway,file=" 20Components-PATHWAY.csv")

# Creating protein information table for each component.
library(tidyverse)
liste_protein <- list()
for(i in 1:ncomp){ #loop for creating liste_protein
  for (i in 1:ncomp){ #loop for each proteins in the components
    e <- unique(liste[[i]][,4]) # unique proteins in the components # Find pathways from pathwayass table to calculate frequency.
    common_pro_path <- intersect(e[1:length(e)],pathwayass[,1]);
    common_pro_path <- pathwayass[which(pathwayass[,1] %in% common_pro_path),c(1,3)];
    count_associated_pathways <- count(common_pro_path[,1])
    # Ordered by alphabetically, so this caused a problem.
    x<- match(count_associated_pathways[,1],e)
    count_associated_pathways <- count_associated_pathways[x,]
    count_associated_pathways <- cbind(e,count_associated_pathways[x,]) # Find drugs from drugass table to calculate frequency.
    common_pro_drug <- intersect(e[1:length(e)],drugass[,2]);
    common_pro_drug <- drugass[which(drugass[,2] %in% common_pro_drug),c(2,3)];
    count_associated_drugs <- count(common_pro_drug[,1])
    liste_protein[[i]] <- data.frame(i,"",count_associated_pathways[,3],count_associated_drugs[,2],e[1:length(e)])
    colnames(liste_protein[[i]]) <- c("component number", "component_score", "number of associated pathways", "number of associated drugs","protein ids")
    }
  return(liste_protein[[i]])
  }
all_bipartite_protein <- do.call(rbind, liste_protein)
write.csv2(all_bipartite_protein,file=" 20Components-PROTEIN.csv")
