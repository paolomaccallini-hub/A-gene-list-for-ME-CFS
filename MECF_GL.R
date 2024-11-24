# file name: MECFS_GL
#
# Rome 25th August 2024
#
path<-getwd() # current directory
#
#-------------------------------------------------------------------------------
# Packages and files
#-------------------------------------------------------------------------------
#
library(dplyr) # for data frame joining
library(rentrez) # for NCBI database
library(httr)
library(jsonlite)
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("ReactomePA",force=T)
library(ReactomePA)
library(biomaRt)
#
#-------------------------------------------------------------------------------
# STRING threshold
#-------------------------------------------------------------------------------
#
STRING.co<-0.7 # cut-off for gene interaction in STRING API 
#
#-------------------------------------------------------------------------------
# Convert gene symbol to NCBI ID (when possible)
#-------------------------------------------------------------------------------
#
Symbol2NCBI<-function(gene.symbol) {
  retries<-100
  for (i in 1:retries) {
    # Use tryCatch to catch errors and warnings
    result<-tryCatch(
      {
        # Perform the search with entrez_search
        NCBI<-lapply(gene.symbol,function(gene) {
          search<-entrez_search(db="gene",term=paste(gene,"[Gene Name] AND human[Organism]"))
        })
      },
      error = function(e) {
        # Check if the error is HTTP 500
        if (grepl("HTTP failure: 500", e$message)) {
          message(paste("Attempt", i, "failed with HTTP 500. Retrying..."))
          return(NULL)
        } else {
          stop(e)  # Stop the loop for other errors
        }
      }
    )
    # If result is not NULL, the search was successful
    if (!is.null(result)) {
      if (is.list(result$ids)) {
        return(NA)
      } else {
        NCBI<-NCBI[[1]]$ids[1]
        return(NCBI)
      }
    }
    # If we've reached this point, the search failed with HTTP 500.
    # Wait before retrying
    Sys.sleep(0.1)
  }
}
#
#-------------------------------------------------------------------------------
# Convert NCBI ID to ENSG ID (when possible)
#-------------------------------------------------------------------------------
#
NCBI2ENSG<-function(NCBI.id) {
  Sys.sleep(0.1)
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  result <- getBM(attributes = c("entrezgene_id", "ensembl_gene_id"),
                  filters = "entrezgene_id",
                  values = NCBI.id,
                  mart = mart)
  return(result$ensembl_gene_id[1])
}
#
#-------------------------------------------------------------------------------
# This function find genes that interact with the input gene, 
# according to STRING API
#-------------------------------------------------------------------------------
#
STRING<-function(gene,threshold) {
  #
  Sys.sleep(0.1)
  species_id<-9606 # Homo sapiens
  #
  # find STRING's identifier for the inputted NCBI ID
  #
  base_url<-"https://string-db.org/api/json/get_string_ids"
  identifiers<-paste(gene,collapse="%0d")
  response<-GET(base_url,query=list(identifiers=identifiers,species=species_id))
  #
  if (status_code(response)==200) {
    response_content<-rawToChar(response$content)
    Encoding(response_content)<-"UTF-8"
    data<-fromJSON(response_content) # parse
  } else {
    stop("Failed to retrieve data. Please check the NCBI IDs and your internet connection.")
  }
  #
  # find interacting genes
  #
  gene_name<-data$preferredName
  base_url<-"https://string-db.org/api/json/network"
  response<-GET(base_url,query=list(identifiers=gene_name,species=species_id))
  #
  if (status_code(response)==200) {
    # 
    response_content<-rawToChar(response$content)
    Encoding(response_content)<-"UTF-8"
    data<-fromJSON(response_content) # parse
    #
    if (nrow(data)>0) {
      interacting_genes<-unique(data$preferredName_B)
      return(interacting_genes)
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}
#
#-------------------------------------------------------------------------------
# This function add the haploinsufficiency and triplosensitivity probability 
# of a gene (if available) accepting as input its NCBI id 
#-------------------------------------------------------------------------------
#
Dosage<-function(NCBI.id) {
  #
  file.name<-paste0(path,"/Data/Dosage_NCBI.txt")
  Dosage.score<-read.csv(file=file.name,sep="\t")
  #
  df<-subset.data.frame(Dosage.score,NCBI==NCBI.id)
  if (nrow(df)==1) {
    return(df)
  } else {
    return(NA)  
  }
  #
}   
#
#-------------------------------------------------------------------------------
# This function add the probability of being a dominant gene
#-------------------------------------------------------------------------------
#
DOMINO<-function(NCBI.id) {
  #
  file.name<-paste0(path,"/Data/DOMINO_NCBI_feb_2019.txt")
  Domino.score<-read.csv(file=file.name,sep="\t")
  #
  if (!is.na(NCBI.id)) {
    df<-subset.data.frame(Domino.score,NCBI==NCBI.id)
    if (nrow(df)==1) {
      PAD<-df$Domino
      return(PAD)
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
  #
}
#
#-------------------------------------------------------------------------------
# Build gene lists
#-------------------------------------------------------------------------------
#
Gene.lists<-function() {
  #
  #-------------------------------------------------------------------------------
  # Custom gene lists
  #-------------------------------------------------------------------------------
  #
  # PrecisionLife and My gene prioritization form UK Biobank data for ME/CFS, from:
  #
  # https://pubmed.ncbi.nlm.nih.gov/36517845/
  #
  UKBB.PL.MC<-c("S100PBP","ATP9A","KCNB1","CLOCK","SLC15A4","TMEM232","GPC5",
                "PHACTR2","AKAP1","USP6NL","CDON","INSR","SLC6A11","SULF2",
                "EBF3","TPTE2P5","ATP8B1","SLC25A15") 
  UKBB.PL.MC<-c(UKBB.PL.MC,"TPTE2") # manual expansion
  #
  # Expand this list using interacting genes
  #
  vector<-UKBB.PL.MC
  for (i in 1:length(vector)) {
    vector<-c(vector,STRING(vector[i],STRING.co))
    print(vector[i])
  }
  UKBB.PL.MC<-unique(vector)
  #
  # Add NCBI ID 
  #
  id<-c()
  for (i in 1:length(UKBB.PL.MC)) {
    id[i]<-Symbol2NCBI(UKBB.PL.MC[i])
  }
  UKBB.PL.MC<-data.frame(id=id,name=UKBB.PL.MC)
  #
  # Severely Ill Patients, from: 
  #
  # https://doi.org/10.1101/2024.09.26.24314417
  #
  SIP<-c("ACADL","BRCA1","CFTR","COX10","HABP2","MFRP","PCLO","PRKN","ZFPM2")
  #
  # Expand this list using interacting genese
  #
  vector<-SIP
  for (i in 1:length(vector)) {
    vector<-c(vector,STRING(vector[i],STRING.co))
    print(vector[i])
  }
  SIP<-unique(vector)
  #
  # Add NCBI ID 
  #
  id<-c()
  for (i in 1:length(SIP)) {
    id[i]<-Symbol2NCBI(SIP[i])
  }
  SIP<-data.frame(id=id,name=SIP)
  #
  # Sporadic MECFS. From: 
  #
  # https://doi.org/10.1161/CIRCGEN.124.004741, 
  # https://doi.org/10.1186/s12967-023-04711-5,
  # https://pubmed.ncbi.nlm.nih.gov/35943990/ 
  # https://doi.org/10.1155/2024/6475425
  #
  MMECFS<-c("NOS3","AKR1C1","AKR1C2","SMPD1","ATP6","COX1")
  MMECFS<-c(MMECFS,"GABRA1","GABRA2","GABRA3","GABRA4","GABRA5","GABRB1","GABRB2",
            "GABRB3","GABRG2","GABRD") # manual expansion
  #
  # Expand this list using interacting genes
  #
  vector<-MMECFS
  for (i in 1:length(vector)) {
    vector<-c(vector,STRING(vector[i],STRING.co))
    print(vector[i])
  }
  MMECFS<-unique(vector)
  #
  # Add NCBI ID 
  #
  id<-c()
  for (i in 1:length(MMECFS)) {
    id[i]<-Symbol2NCBI(MMECFS[i])
  }
  MMECFS<-data.frame(id=id,name=MMECFS)
  #
  # Klimas's study with 23andME data. From: 
  #
  # https://doi:10.21203/rs.3.rs-3171709/v1
  #
  Klimas23<-c("MUC16","MUC19","MUC22")
  #
  # Expand this list using interacting genes
  #
  vector<-Klimas23
  for (i in 1:length(vector)) {
    vector<-c(vector,STRING(vector[i],STRING.co))
    print(vector[i])
  }
  Klimas23<-unique(vector)
  #
  # Add NCBI ID 
  #
  id<-c()
  for (i in 1:length(Klimas23)) {
    id[i]<-Symbol2NCBI(Klimas23[i])
  }
  Klimas23<-data.frame(id=id,name=Klimas23)
  #
  #-------------------------------------------------------------------------------
  # Build a list
  #-------------------------------------------------------------------------------
  #
  gene_lists<-list()
  i<-1
  gene_lists[[i]]<-UKBB.PL.MC
  names(gene_lists)[i]<-"UKBB.PL.MC"
  #
  i<-i+1
  gene_lists[[i]]<-SIP
  names(gene_lists)[i]<-"SIP"
  #
  i<-i+1
  gene_lists[[i]]<-MMECFS
  names(gene_lists)[i]<-"MMECFS"
  #
  i<-i+1
  gene_lists[[i]]<-Klimas23
  names(gene_lists)[i]<-"Klimas23"
  #
  #-------------------------------------------------------------------------------
  # Build a data frame with all the genes and the number of lists they are in 
  #-------------------------------------------------------------------------------
  #
  all.genes<-gene_lists[[1]]
  for (i in 2:length(gene_lists)) {
    all.genes<-rbind(all.genes,gene_lists[[i]])
  }
  all.genes<-unique(all.genes)
  #
  # For each gene, indicate how many gene lists include it and the names of the lists
  #
  count<-rep(0,nrow(all.genes))
  list.name<-rep(NA,nrow(all.genes))
  for (i in 1:nrow(all.genes)) { # select the gene
    for (j in 1:length(gene_lists)) { # select the list
      df<-subset.data.frame(gene_lists[[j]],name==all.genes$name[i])
      if (nrow(df)>0) {
        count[i]<-count[i]+1
        if (count[i]==1) {
          list.name[i]<-names(gene_lists)[j] 
        } else {
          list.name[i]<-paste0(list.name[i],",",names(gene_lists)[j])
        }
      } 
    }
  }
  all.genes$list.count<-count
  all.genes$list.name<-list.name
  colnames(all.genes)<-c("NCBI.id","name","list.count","list.name")
  #
  # Add ENSG ID 
  #
  ENSG.id<-c()
  for (i in 1:nrow(all.genes)) {
    ENSG.id[i]<-NCBI2ENSG(all.genes$NCBI.id[i])
  }
  all.genes$ENSG.id<-ENSG.id
  #
  # Add DOMINO probability 
  #
  pDI<-c()
  for (i in 1:nrow(all.genes)) {
    pDI[i]<-DOMINO(all.genes$NCBI.id[i])
  }
  all.genes$pDI<-pDI
  #
  # Add dosage sensitivity 
  #
  all.genes$pHaplo<-rep(NA,nrow(all.genes))
  all.genes$pTriplo<-rep(NA,nrow(all.genes))
  for (i in 1:nrow(all.genes)) {
    df<-Dosage(all.genes$NCBI.id[i])
    if (is.data.frame(df)) {
      all.genes$pHaplo[i]<-df$pHaplo
      all.genes$pTriplo[i]<-df$pTriplo
    } 
  }
  #
  # Edit 
  #
  all.genes<-subset.data.frame(all.genes,select=c("NCBI.id","ENSG.id","name",
                                                  "pDI","pHaplo","pTriplo",
                                                  "list.count","list.name"))
  all.genes<-subset.data.frame(all.genes,list.count!=0)
  #
  # Save a copy 
  #
  file.name<-file.path(path,"All_genes.tsv")
  write.table(all.genes,file=file.name,quote=F,row.names=F,col.names=T,sep="\t")
  #
  # Return the output to the main unit
  #
  return(all.genes)
}
#
#-------------------------------------------------------------------------------
# Build the gene list
#-------------------------------------------------------------------------------
#
all.genes<-Gene.lists()
#
#-------------------------------------------------------------------------------
# Cluster analysis
#-------------------------------------------------------------------------------
#
# Load gene list
#
file.name<-file.path(path,"All_genes.tsv")
all.genes<-read.csv(file=file.name,sep="\t")
ncbi_ids<-all.genes$NCBI.id
#
# Perform pathway enrichment analysis using Reactome
#
pathway_results<-enrichPathway(gene=ncbi_ids,organism="human")
#
head(pathway_results) # View the top results
#
# Plot the top 10 enriched pathways
#
barplot(pathway_results,showCategory=10,title="Top Enriched Pathways")
#
# Optional: Dot plot for visualization
#
dotplot(pathway_results,showCategory=10,title="Dotplot of Enriched Pathways")


