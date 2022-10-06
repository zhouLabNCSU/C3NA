#' Initiate the C3NA object with Correlation Calculation and Topology Overlap Matrix
#' 
#' @description Generate the initiate C3NA object, including the count matrix, the taxonomic table, and the multi-taxa stacked count matrix for the C3NA analysis. This step also includes the sparcc correlation calculation. Note: sparcc correlation can be extremely computationally expensive. The results will be summarized along with the signed network calcualted from the topology overlap matrix. Finally, different clusters of taxa will be calculated based on a minimal module size range from 3 to 40 (default). 
#' 
#' @param phyloseqObj (Required) Phyloseq \linkS4class{phyloseq} object. This should first undergo validatePhyloseq to ensure the diagnosis column are present.  
#' @param phenotype (Required) The desired phenotype that present under the diagnosis column in the metadata from phyloseqObj
#' @param prevTrh (Required) Prevalence threshold of the samples, which is a number between 0 and 1. E.g., the default 0.1 represents 10% of the samples need to have given taxa. 
#' @param nBootstrap (Required) Number of bootstrap for the sparcc command. Warning: this is a very computational step. 
#' @param nMinTotalCount (Required) The Minimal number of reads per sample. Default: 1,000.
#' @param nCPUs (Optional) Parallel computation for the sparccboot function. Default: 1. 
#' @param minModuleSize (Optional) The lowest number for the minimal size for each cluster. Default: 3. 
#' @param maxModuleSize (Optional) The highest number for the minimal size for each cluster. Default: 40. 
#' @param seed (Optional)
#' 
#' @importFrom magrittr %>%
#' @importFrom dynamicTreeCut cutreeDynamic 
#' @importFrom phyloseq prune_samples tax_table otu_table sample_data phyloseq subset_samples sample_sums tax_glom taxa_are_rows taxa_sums taxa_are_rows<-
#' @importFrom WGCNA TOMsimilarity labels2colors 
#' @importFrom dplyr mutate group_by summarise_if distinct filter rename left_join summarise n arrange ungroup pull row_number select across
#' @importFrom tibble column_to_rownames 
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom reshape2 melt
#' @importFrom metagMisc phyloseq_filter_prevalence
#' @importFrom SpiecEasi sparccboot pval.sparccboot
#' @importFrom stats na.omit p.adjust cutree as.hclust hclust as.dist 
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom cluster silhouette agnes 
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @return A list object of C3NA Objects
#' @export
#' @examples
#' data(CRC_Phyloseq)
#' curPhyloseq = validatePhloseq(phyloseqObj = CRC_Phyloseq)
#' 
#' # These steps are commented out due to time consuming step. The post sparcc correlation data will be 
#' # to avoid the step
#' phyloseq_Normal = phyloseq::subset_samples(physeq = CRC_Phyloseq, diagnosis == "Normal")
#' # C3NAObj_Normal = initiateC3NA(phyloseqObj = phyloseq_Normal, prevTrh = 0.1, 
#' #                               nCPUs = 10, nBootstrap = 1000, 
#' #                               nMinTotalCount = 1000, phenotype = "Normal", 
#' #                               seed = 100)
#' 
#' phyloseq_Cancer = phyloseq::subset_samples(physeq = CRC_Phyloseq, diagnosis == "Cancer")
#' # C3NAObj_Cancer = initiateC3NA(phyloseqObj = phyloseq_Cancer, prevTrh = 0.1, 
#' #                               nCPUs = 10, nBootstrap = 1000, 
#' #                               nMinTotalCount = 1000, phenotype = "Cancer", 
#' #                               seed = 100)
#' 
initiateC3NA <- function(phyloseqObj = phyloseqObj,
                         prevTrh = 0.1, 
                         nCPUs = 1, nBootstrap = 1000,
                         nMinTotalCount = 1000, phenotype = NA,
                         minModuleSize = 3, maxModuleSize = 40,
                         seed = 100){
  # Text Bar
  pb = txtProgressBar(min = 0, max = 4, initial = 0, style = 3)
  
  # Misc Variables
  taxaLvls <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  taxaHeaders <- c("p", "c", "o", "f", "g", "s")
  
  # Error Checking
  if(is.na(phenotype)){stop("Please provide a phenotype name.")} 
  if(!class(phyloseqObj) == "phyloseq"){
    stop(paste0("The object is not a phyloseq object, please correctly import the ",
                "OTUs/ASVs table, taxa table and sample data as phyloseq obejcts."))}
  
  # Filter out low OTUs count samples, default > 1000
  phyloseqObj = prune_samples(sample_sums(phyloseqObj)>=nMinTotalCount, phyloseqObj)

  setTxtProgressBar(pb, 1)
  # Extract Info from Phyloseq Object
  tempTaxa = as.data.frame(as.matrix(tax_table(phyloseqObj)))
  tempTaxa$OTUs = rownames(otu_table(phyloseqObj))
  curTaxaTable <- as.data.frame(as.matrix(tax_table(phyloseqObj)))
  for(t in seq_along(taxaLvls)){
    curLevel = taxaLvls[t]
    psTempOri = tax_glom(phyloseqObj, curLevel, NArm = FALSE)
    psTemp = phyloseq_filter_prevalence(psTempOri, 
                                        prev.trh = prevTrh, abund.trh = NULL,
                                        threshold_condition = "AND", abund.type = "total")
    curOTUs = as.data.frame(as.matrix(otu_table(psTemp)))
    curOTUs$TEMP = as.character(as.data.frame(as.matrix(tax_table(psTemp)))[, curLevel])
    curOTUs <- curOTUs %>%
      mutate(TEMP = ifelse(is.na(TEMP), "NA", TEMP)) %>%
      group_by(TEMP) %>%
      summarise_if(is.numeric, sum) %>%
      column_to_rownames("TEMP")
    
    # Replace NAs and Remove Useless Taxa
    curOTUs[is.na(curOTUs)] <- 0
    curOTUs <- curOTUs[rowSums(curOTUs)>0, ]
    
    if(t == 1){
      # OTUs Table
      rownames(curOTUs) = paste0(taxaHeaders[t], "_", rownames(curOTUs))
      rawCountMatrix = curOTUs
      # Accessory Information
      accInfoTable = data.frame(nTotalDisease = ncol(curOTUs),
                                nNonZeroDisease = rowSums(curOTUs != 0),
                                nCurRowSumsDisease = rowSums(curOTUs),
                                oriNTaxa = dim(otu_table(psTempOri))[1],
                                filteredNTaxa = dim(curOTUs)[1],
                                nFullRowSums = rowSums(curOTUs))
    } else{
      rownames(curOTUs) = paste0(taxaHeaders[t], "_", rownames(curOTUs))
      rawCountMatrix = rbind(rawCountMatrix, curOTUs)
      accInfoTable_temp = data.frame(nTotalDisease = ncol(curOTUs),
                                     nNonZeroDisease = rowSums(curOTUs != 0),
                                     nCurRowSumsDisease = rowSums(curOTUs),
                                     oriNTaxa = dim(otu_table(psTempOri))[1],
                                     filteredNTaxa = dim(curOTUs)[1],
                                     nFullRowSums = rowSums(curOTUs))
      accInfoTable = rbind(accInfoTable, accInfoTable_temp)
    }
  }

  # SparCC 
  setTxtProgressBar(pb, 2)
  cat("\nSparcc Caluclation: Time-Consuming Step...")
  set.seed(seed)
  allTaxaMatrix_BSCor <- SpiecEasi::sparccboot(data = t(rawCountMatrix), 
                                               R = nBootstrap, ncpus = nCPUs)
  allTaxaMatrix_BSCorP <- SpiecEasi::pval.sparccboot(x = allTaxaMatrix_BSCor)
  
  # Summarize SparCC Results
  cors <- allTaxaMatrix_BSCorP$cors
  pvals <- allTaxaMatrix_BSCorP$pvals
  sparCCpcors <- diag(0.5, nrow = dim(rawCountMatrix)[1], ncol = dim(rawCountMatrix)[1])
  sparCCpcors[upper.tri(sparCCpcors, diag=FALSE)] <- cors
  sparCCpcors <- sparCCpcors + t(sparCCpcors)
  
  sparCCpval <- diag(0.5, nrow = dim(rawCountMatrix)[1], ncol = dim(rawCountMatrix)[1])
  sparCCpval[upper.tri(sparCCpval, diag=FALSE)] <- pvals
  sparCCpval <- sparCCpval + t(sparCCpval)

  rownames(sparCCpcors) <- rownames(rawCountMatrix)
  colnames(sparCCpcors) <- rownames(rawCountMatrix)
  rownames(sparCCpval) <- rownames(rawCountMatrix)
  colnames(sparCCpval) <- rownames(rawCountMatrix)
  
  sparccCor_processed <- sparCCpcors  %>% .getUpperTri() %>% reshape2::melt() %>% na.omit() %>% rename(cor = value)
  sparccP_processed <- sparCCpval  %>% .getUpperTri() %>% reshape2::melt() %>% na.omit() %>% rename(p = value)
  
  sparccTable <- left_join(sparccCor_processed, sparccP_processed, by = c("Var1", "Var2")) %>%
    filter(Var1 != Var2) %>%
    mutate(fdr = p.adjust(p, method = "BH"))
  
  # Signed TOM Similairty 
  setTxtProgressBar(pb, 3)
  signedNetwork = WGCNA::TOMsimilarity(adjMat = sparCCpcors, TOMType = "signed", verbose = FALSE)
  colnames(signedNetwork) = rownames(signedNetwork) = colnames(sparCCpcors)
  
  # Calculate the Modular Data Based on Minimal Cluster Size
  module_df <- data.frame(
    taxaID = NULL, colors = NULL,
    minSize = NULL, dataset = NULL
  )
  for(minSize in (minModuleSize:maxModuleSize)){
    consTree = hclust(as.dist(1-signedNetwork), method = "complete");
    # Module identification using dynamic tree cut:
    unmergedLabels = dynamicTreeCut::cutreeDynamic(dendro = consTree, 
                                   distM = 1-signedNetwork,
                                   deepSplit = 2, 
                                   cutHeight = 0.995,
                                   minClusterSize = minSize,
                                   pamRespectsDendro = FALSE, verbose = FALSE);
    unmergedColors = WGCNA::labels2colors(unmergedLabels)
    module_df_temp <- data.frame(
      taxaID = consTree$labels,
      colors = unmergedColors,
      minSize = minSize,
      phenotype = phenotype
    )
    module_df = rbind(module_df, module_df_temp)   
  }
  
  # Return
  setTxtProgressBar(pb, 4)
  returnlist = list(
    "oriTaxTable" = as.data.frame(as.matrix(tax_table(phyloseqObj))),
    "filteredCountMatrix" = rawCountMatrix,
    "sparCCTable" = sparccTable,
    "moduleData" = module_df,
    "misc" = list(
      "phenotype" = phenotype,
      "accInfoTable" = accInfoTable,
      "prevTrh" = prevTrh, 
      "nBootstrap" = nBootstrap,
      "nMinTotalCount" = nMinTotalCount,
      "minModuleSize" = minModuleSize,
      "maxModuleSize" = maxModuleSize
    )
  )
  return(returnlist)
}


#' Incoporate the C3NA object with Different Correlation Methods
#' 
#' @description Generate the initiate C3NA object, including the count matrix, the taxonomic table, and the multi-taxa stacked count matrix for the C3NA analysis. This step also includes the sparcc correlation calculation. Note: sparcc correlation can be extremely computationally expensive. The results will be summarized along with the signed network calcualted from the topology overlap matrix. Finally, different clusters of taxa will be calculated based on a minimal module size range from 3 to 40 (default). 
#' 
#' @param corMatrix (Required) Symmetric correlation matrix, required column and rownames to be in the format of taxonomic name and level, e.g. "p_bacterialName".  
#' @param phyloseqObj (Required) Phyloseq \linkS4class{phyloseq} object. This should first undergo validatePhyloseq to ensure the diagnosis column are present.  
#' @param phenotype (Required) The desired phenotype that present under the diagnosis column in the metadata from phyloseqObj
#' @param prevTrh (Required) Prevalence threshold of the samples, which is a number between 0 and 1. E.g., the default 0.1 represents 10% of the samples need to have given taxa. 
#' @param nBootstrap (Required) Number of bootstrap for the sparcc command. Warning: this is a very computational step. 
#' @param nMinTotalCount (Required) The Minimal number of reads per sample. Default: 1,000.
#' @param nCPUs (Optional) Parallel computation for the sparccboot function. Default: 1. 
#' @param minModuleSize (Optional) The lowest number for the minimal size for each cluster. Default: 3. 
#' @param maxModuleSize (Optional) The highest number for the minimal size for each cluster. Default: 40. 
#' @param seed (Optional)
#' 
#' @importFrom magrittr %>%
#' @importFrom dynamicTreeCut cutreeDynamic 
#' @importFrom phyloseq prune_samples tax_table otu_table sample_data phyloseq subset_samples sample_sums tax_glom taxa_are_rows taxa_sums taxa_are_rows<-
#' @importFrom WGCNA TOMsimilarity labels2colors 
#' @importFrom dplyr mutate group_by summarise_if distinct filter rename left_join summarise n arrange ungroup pull row_number select across
#' @importFrom tibble column_to_rownames 
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom reshape2 melt
#' @importFrom metagMisc phyloseq_filter_prevalence
#' @importFrom SpiecEasi sparccboot pval.sparccboot
#' @importFrom stats na.omit p.adjust cutree as.hclust hclust as.dist 
#' @importFrom randomcoloR distinctColorPalette
#' @importFrom cluster silhouette agnes 
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @return A list object of C3NA Objects
#' @export
#' @examples
#' data(CRC_Phyloseq)
#' curPhyloseq = validatePhloseq(phyloseqObj = CRC_Phyloseq)
#' 
#' # These steps are commented out due to time consuming step. The post sparcc correlation data will be 
#' # to avoid the step
#' phyloseq_Cancer = phyloseq::subset_samples(physeq = CRC_Phyloseq, diagnosis == "Cancer")
#' stackedTaxaMatrix = getStackedTaxaMatrix(phyloseqObj = phyloseq_Cancer, phenotype = "Cancer")
#' 
#' # Correlation method
#' testCorMatrix = cor(t(stackedTaxaMatrix))
#' # C3NAObj_Normal = initiateC3NA_DiffCorr(corMatrix = testCorMatrix,
#' #                                        phyloseqObj = phyloseq_Cancer, 
#' #                                        nMinTotalCount = 1000, phenotype = "Cancer" 
#' #                                        )
#' 
initiateC3NA_DiffCorr <- function(corMatrix = corMatrix,  
                                  prevTrh = 0.1, 
                                  phyloseqObj = phyloseqObj,
                                  nCPUs = "None", nBootstrap = "None",
                                  nMinTotalCount = 1000, phenotype = NA,
                                  minModuleSize = 3, maxModuleSize = 40,
                                  seed = "None"){
  # Text Bar
  pb = txtProgressBar(min = 0, max = 4, initial = 0, style = 3)
  
  # Misc Variables
  taxaLvls <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  taxaHeaders <- c("p", "c", "o", "f", "g", "s")
  
  # Error Checking
  if(!isSymmetric(corMatrix)){stop("Please provide a symmetric correlation matrix")}
  if(is.na(phenotype)){stop("Please provide a phenotype name.")} 
  if(!class(phyloseqObj) == "phyloseq"){
    stop(paste0("The object is not a phyloseq object, please correctly import the ",
                "OTUs/ASVs table, taxa table and sample data as phyloseq obejcts."))}
  
  # Filter out low OTUs count samples, default > 1000
  phyloseqObj = phyloseq::prune_samples(phyloseq::sample_sums(phyloseqObj)>=nMinTotalCount, phyloseqObj)
  
  setTxtProgressBar(pb, 1)
  # Extract Info from Phyloseq Object
  tempTaxa = as.data.frame(as.matrix(tax_table(phyloseqObj)))
  tempTaxa$OTUs = rownames(otu_table(phyloseqObj))
  curTaxaTable <- as.data.frame(as.matrix(tax_table(phyloseqObj)))
  for(t in seq_along(taxaLvls)){
    curLevel = taxaLvls[t]
    psTempOri = tax_glom(phyloseqObj, curLevel, NArm = FALSE)
    psTemp = phyloseq_filter_prevalence(psTempOri, 
                                        prev.trh = prevTrh, abund.trh = NULL,
                                        threshold_condition = "AND", abund.type = "total")
    curOTUs = as.data.frame(as.matrix(otu_table(psTemp)))
    curOTUs$TEMP = as.character(as.data.frame(as.matrix(tax_table(psTemp)))[, curLevel])
    curOTUs <- curOTUs %>%
      mutate(TEMP = ifelse(is.na(TEMP), "NA", TEMP)) %>%
      group_by(TEMP) %>%
      summarise_if(is.numeric, sum) %>%
      column_to_rownames("TEMP")
    
    # Replace NAs and Remove Useless Taxa
    curOTUs[is.na(curOTUs)] <- 0
    curOTUs <- curOTUs[rowSums(curOTUs)>0, ]
    
    if(t == 1){
      # OTUs Table
      rownames(curOTUs) = paste0(taxaHeaders[t], "_", rownames(curOTUs))
      rawCountMatrix = curOTUs
      # Accessory Information
      accInfoTable = data.frame(nTotalDisease = ncol(curOTUs),
                                nNonZeroDisease = rowSums(curOTUs != 0),
                                nCurRowSumsDisease = rowSums(curOTUs),
                                oriNTaxa = dim(otu_table(psTempOri))[1],
                                filteredNTaxa = dim(curOTUs)[1],
                                nFullRowSums = rowSums(curOTUs))
    } else{
      rownames(curOTUs) = paste0(taxaHeaders[t], "_", rownames(curOTUs))
      rawCountMatrix = rbind(rawCountMatrix, curOTUs)
      accInfoTable_temp = data.frame(nTotalDisease = ncol(curOTUs),
                                     nNonZeroDisease = rowSums(curOTUs != 0),
                                     nCurRowSumsDisease = rowSums(curOTUs),
                                     oriNTaxa = dim(otu_table(psTempOri))[1],
                                     filteredNTaxa = dim(curOTUs)[1],
                                     nFullRowSums = rowSums(curOTUs))
      accInfoTable = rbind(accInfoTable, accInfoTable_temp)
    }
  }
  
  setTxtProgressBar(pb, 2)
  corTable = corMatrix %>% .getUpperTri() %>%
    reshape2::melt() %>% na.omit() %>% mutate(Var1 = as.character(Var1), 
                                              Var2 = as.character(Var2)) %>%
    
    rename(cor = value, from = Var1, to = Var2) %>% rowwise() %>%
    mutate(source = sort(c(from, to))[1], target = sort(c(from, to))[2]) %>%
    mutate(taxaComboName = paste0(source, "__", target))
  
  # Signed TOM Similairty 
  setTxtProgressBar(pb, 3)
  corMatrix2 = corMatrix
  corMatrix2[which(corMatrix2< 0)] = 0
  corMatrix2[which(corMatrix2>1)] = 1
  signedNetwork = WGCNA::TOMsimilarity(adjMat = (corMatrix2), TOMType = "signed", verbose = FALSE)
  colnames(signedNetwork) = rownames(signedNetwork) = colnames(corMatrix2)
  
  # Calculate the Modular Data Based on Minimal Cluster Size
  module_df <- data.frame(
    taxaID = NULL, colors = NULL,
    minSize = NULL, dataset = NULL
  )
  for(minSize in (minModuleSize:maxModuleSize)){
    consTree = hclust(as.dist(1-signedNetwork), method = "complete");
    # Module identification using dynamic tree cut:
    unmergedLabels = dynamicTreeCut::cutreeDynamic(dendro = consTree, 
                                                   distM = 1-signedNetwork,
                                                   deepSplit = 2, 
                                                   cutHeight = 0.995,
                                                   minClusterSize = minSize,
                                                   pamRespectsDendro = FALSE, verbose = FALSE);
    unmergedColors = WGCNA::labels2colors(unmergedLabels)
    module_df_temp <- data.frame(
      taxaID = consTree$labels,
      colors = unmergedColors,
      minSize = minSize,
      phenotype = phenotype
    )
    module_df = rbind(module_df, module_df_temp)   
  }
  
  # Output Data 
  setTxtProgressBar(pb, 4)
  returnlist = list(
    "oriTaxTable" = as.data.frame(as.matrix(tax_table(phyloseqObj))),
    "filteredCountMatrix" = rawCountMatrix,
    "sparCCTable" = corTable,
    "moduleData" = module_df,
    "misc" = list(
      "phenotype" = phenotype,
      "accInfoTable" = accInfoTable,
      "prevTrh" = prevTrh, 
      "nBootstrap" = nBootstrap,
      "nMinTotalCount" = nMinTotalCount,
      "minModuleSize" = minModuleSize,
      "maxModuleSize" = maxModuleSize
    )
  )
}

#' Obtained the stacked-taxa table
#' 
#' @description Generate the stacked-taxa tabel for other correlation methods calculation outside C3NA. 
#' 
#' @param phyloseqObj (Required) Phyloseq \linkS4class{phyloseq} object. This should first undergo validatePhyloseq to ensure the diagnosis column are present.  
#' @param phenotype (Required) The desired phenotype that present under the diagnosis column in the metadata from phyloseqObj
#' @param prevTrh (Required) Prevalence threshold of the samples, which is a number between 0 and 1. E.g., the default 0.1 represents 10% of the samples need to have given taxa. 
#' @param nMinTotalCount (Required) The Minimal number of reads per sample. Default: 1,000.
#' 
#' @return countMatrix
#' @export
#' @examples
#' data(CRC_Phyloseq)
#' curPhyloseq = validatePhloseq(phyloseqObj = CRC_Phyloseq)
#' 
#' # These steps are commented out due to time consuming step. The post sparcc correlation data will be 
#' # to avoid the step
#' phyloseq_Cancer = phyloseq::subset_samples(physeq = CRC_Phyloseq, diagnosis == "Cancer")
#' stackedTaxaMatrix = getStackedTaxaMatrix(phyloseqObj = phyloseq_Cancer, phenotype = "Cancer")
#' 
#' # Correlation method
#' testCorMatrix = cor(t(stackedTaxaMatrix))
getStackedTaxaMatrix <- function(prevTrh = 0.1, 
                                 phyloseqObj = phyloseqObj,
                                 nMinTotalCount = 1000, phenotype = NA){
  # Misc Variables
  taxaLvls <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  taxaHeaders <- c("p", "c", "o", "f", "g", "s")
  
  # Error Checking
  if(is.na(phenotype)){stop("Please provide a phenotype name.")} 
  if(!class(phyloseqObj) == "phyloseq"){
    stop(paste0("The object is not a phyloseq object, please correctly import the ",
                "OTUs/ASVs table, taxa table and sample data as phyloseq obejcts."))}
  
  # Filter out low OTUs count samples, default > 1000
  phyloseqObj = phyloseq::prune_samples(phyloseq::sample_sums(phyloseqObj)>=nMinTotalCount, phyloseqObj)
  
  # Extract Info from Phyloseq Object
  tempTaxa = as.data.frame(as.matrix(tax_table(phyloseqObj)))
  tempTaxa$OTUs = rownames(otu_table(phyloseqObj))
  curTaxaTable <- as.data.frame(as.matrix(tax_table(phyloseqObj)))
  for(t in seq_along(taxaLvls)){
    curLevel = taxaLvls[t]
    psTempOri = tax_glom(phyloseqObj, curLevel, NArm = FALSE)
    psTemp = phyloseq_filter_prevalence(psTempOri, 
                                        prev.trh = prevTrh, abund.trh = NULL,
                                        threshold_condition = "AND", abund.type = "total")
    curOTUs = as.data.frame(as.matrix(otu_table(psTemp)))
    curOTUs$TEMP = as.character(as.data.frame(as.matrix(tax_table(psTemp)))[, curLevel])
    curOTUs <- curOTUs %>%
      mutate(TEMP = ifelse(is.na(TEMP), "NA", TEMP)) %>%
      group_by(TEMP) %>%
      summarise_if(is.numeric, sum) %>%
      column_to_rownames("TEMP")
    
    # Replace NAs and Remove Useless Taxa
    curOTUs[is.na(curOTUs)] <- 0
    curOTUs <- curOTUs[rowSums(curOTUs)>0, ]
    
    if(t == 1){
      # OTUs Table
      rownames(curOTUs) = paste0(taxaHeaders[t], "_", rownames(curOTUs))
      rawCountMatrix = curOTUs
    } else{
      rownames(curOTUs) = paste0(taxaHeaders[t], "_", rownames(curOTUs))
      rawCountMatrix = rbind(rawCountMatrix, curOTUs)
    }
  }
  
  return(rawCountMatrix)
}

# Remove all lower triangle of a corMatrix
.getUpperTri <- function(corMatrix){
  corMatrix[lower.tri(corMatrix)]<- NA
  return(corMatrix)
}

#' Obtain the Module Information from Consensus-Based Investigation
#' 
#' @description Extract the module information based on the consensus results. 
#' 
#' @param C3NAObj (Required) The C3NAObj with a single phenotype and after the initiateC3NA function.
#' @param selectedPatterns (Required) The selected unique patterns identified in the shiny application. 
#' @param nModules (Required) The optimal number of clusters identified by the shiny application. 
#' 
#' @importFrom dplyr rowwise relocate 
#' @importFrom WGCNA standardColors checkAdjMat goodSamplesGenes mtd.apply overlapTable prepComma spaste multiData
#' @importFrom igraph graph.data.frame as_adjacency_matrix as_data_frame
#' @importFrom dynamicTreeCut indentSpaces printFlush 
#' @export
#' 
getOptMods <- function(C3NAObj = C3NAObj,
                       selectedPatterns = NA,
                       nModules = NA){
  if(!is.numeric(selectedPatterns) | !is.numeric(nModules)){
    stop("Please copy and paste the codes directly from the moduleEvals shiny application.\n")
  }
  ## Text Bar
  pb = txtProgressBar(min = 0, max = 4, initial = 0, style = 3)

  ## Calculate the moduleStat
  module_df = C3NAObj$moduleData
  # Check for duplicated patterns
  setTxtProgressBar(pb, 1)
  ## Googing over each of the minSizes and retrieve cluster module connection information
  moduleData = module_df %>%
    filter(minSize %in% selectedPatterns)
  for(minN in selectedPatterns){
    moduleData_temp = moduleData %>% filter(minSize == minN)
    consensusMatrix_temp = expand.grid(from = moduleData_temp$taxaID, to = moduleData_temp$taxaID)
    consensusMatrix_temp$Category = NA
    consensusMatrix_temp$BinarySimilarity = NA
    for(color in unique(moduleData_temp$colors)){
      tempD = moduleData_temp %>% filter(colors == color)
      consensusMatrix_temp = consensusMatrix_temp %>%
        mutate(BinarySimilarity = ifelse(from %in% tempD$taxaID & to %in% tempD$taxaID, 1, BinarySimilarity))
    }
    consensusMatrix_temp$Category = paste0("MinSize_", minN)
    consensusMatrix_temp$BinarySimilarity[is.na(consensusMatrix_temp$BinarySimilarity)] = 0
    if(minN == 3){
      consensusMatrix = consensusMatrix_temp
    } else{
      consensusMatrix = rbind(consensusMatrix, consensusMatrix_temp)
    }
  }

  # Expands the consensMatrix and caclaulte the average consensus
  consensusMatrix_wide = consensusMatrix %>%
    pivot_wider(names_from = Category, values_from = BinarySimilarity)
  consensusMatrix_wide$rowMeans = rowMeans(consensusMatrix_wide[, 3:ncol(consensusMatrix_wide)])
  consensusMatrix_Final = consensusMatrix_wide %>% dplyr::select(from, to, rowMeans)
  consensusMatrix_Final_wide = consensusMatrix_Final %>%
    pivot_wider(names_from = to, values_from = rowMeans)
  consensusMatrix_Final_wide <- consensusMatrix_Final_wide %>% column_to_rownames("from")
  
  ## Calculate the Consensus Score Matrix
  setTxtProgressBar(pb, 2)
  consensusMatrix_Final = consensusMatrix_Final %>%
    mutate(from = as.character(from), 
           to = as.character(to)) %>%
    rowwise() %>% 
    mutate(source = sort(c(from, to))[1],
           target = sort(c(from, to))[2]) %>%
    mutate(taxaComboName = paste0(source, "__", target)) %>% 
    dplyr::select(-from, -to) %>%
    distinct() %>%
    rename(ConsensusProp = rowMeans)
  
  ## Calculation of the optimal cluseter number based on the consensus matrix
  setTxtProgressBar(pb, 3)
  curHclust = cluster::agnes(x = (1-consensusMatrix_Final_wide), 
                             method = "complete", metric = "euclidean")
  ## Obtain the clustering information for the taxa
  curTree = cutree(as.hclust(curHclust), k = nModules) 
  groups <- as.data.frame(curTree)
  order = rownames(consensusMatrix_Final_wide[curHclust[["order"]],])
  groupFreqTable = as.data.frame(table(groups[order,]))
  curClusterTable = data.frame(taxaName = order, clusterID = groups[order,])
  
  ## Extract the Intramodular Correlations
  
  # Update then Return the C3NA Object
  setTxtProgressBar(pb, 4)
  C3NAObj$misc$modPrep = TRUE
  C3NAObj$consensus$modules = curClusterTable
  C3NAObj$consensus$selectedPatterns = selectedPatterns
  C3NAObj$consensus$consensusProp = consensusMatrix_Final
  C3NAObj$consensus$nModules = nModules
  return(C3NAObj)
}

#' Obtain the Module Information from Consensus-Based Investigation
#' 
#' @description Extract the module information based on the consensus results. 
#' 
#' @param C3NAObj_Comparison (Required) The C3NAObj with a single phenotype as Disease.
#' @param C3NAObj_Reference (Required) The C3NAObj with a single phenotype as Control.
#' @param corCutoff (Required) A positive correlation cut-off point, default is 0.2. 
#' @param fdr (Required) BH-Adjusted p-value, default is 0.05.
#' @param unusefulTaxa (Optional) A list of taxa names that could not be interpreated biologically, e.g. human_gut, ambiguous_taxa
#' @param nBootstrap (Optional) Number of bootstrap for the module preservation. 
#' @param verbose (Optinal) Verbose for the module preservation analysis, the higher the more detailed output. 
#' 
#' @importFrom dplyr rowwise relocate 
#' @importFrom WGCNA standardColors checkAdjMat goodSamplesGenes mtd.apply overlapTable prepComma spaste multiData
#' @importFrom igraph graph.data.frame as_adjacency_matrix as_data_frame
#' @importFrom dynamicTreeCut indentSpaces printFlush 
#' @export
#' 
comparePhenotypes <- function(C3NAObj_Comparison = C3NAObj1,
                              C3NAObj_Reference = C3NAObj2,
                              corCutoff = 0.2, fdr = 0.05,
                              unusefulTaxa = NA, verbose = 0, 
                              nBootstrap = 300
                              ){
  if(!is.numeric(corCutoff) | (corCutoff > 1 & corCutoff < 0)){
    stop("Please provide a valid correlation cut-off value, a numeric number between 0 and 1. \n")
  }
  if(!is.numeric(corCutoff) | (corCutoff > 1 & corCutoff < 0)){
    stop("Please provide a valid correlation cut-off value, a numeric number between 0 and 1. \n")
  }
  pb = txtProgressBar(min = 0, max = 7, initial = 0, style = 3)
  
  ## Misc Generate the Meaningless Taxa Table
  taxaLvls_Abbre = c("p", "c", "o", "f", "g", "s")
  if(length(unusefulTaxa) == 0 | is.na(unusefulTaxa)){
    unusefulTaxa = c("NA", "ambiguous", "uncultured", "Incertae_Sedis", 
                     "gut_metagenome", "metagenome", "unidentified", "unidentified_rumen", "unidentified_marine",
                     "uncultured_bacterium", "uncultured_organism", "uncultured_rumen", "human_gut")
  }
  removedTaxaTable = expand.grid(taxaLvls_Abbre = taxaLvls_Abbre, unusefulTaxa = unusefulTaxa)
  removedTaxaTable$removedTaxa = paste0(removedTaxaTable$taxaLvls_Abbre, "_", 
                                        removedTaxaTable$unusefulTaxa)
  
  ## Step 1. Loading the associated processed data
  setTxtProgressBar(pb, 1)
  curDiseasePhenotype = C3NAObj_Comparison$misc$phenotype
  curDisease_SparccP = C3NAObj_Comparison$sparCCTable
  curDisease_ModuleData = C3NAObj_Comparison$consensus$modules
  diseaseTaxaTable = C3NAObj_Comparison$oriTaxTable
  diseaseTaxaTable$Phylum = paste0("p_", diseaseTaxaTable$Phylum)
  diseaseTaxaTable$Class = paste0("c_", diseaseTaxaTable$Class)
  diseaseTaxaTable$Order = paste0("o_", diseaseTaxaTable$Order)
  diseaseTaxaTable$Family = paste0("f_", diseaseTaxaTable$Family)
  diseaseTaxaTable$Genus = paste0("g_", diseaseTaxaTable$Genus)
  diseaseTaxaTable$Species = paste0("s_", diseaseTaxaTable$Species)
  
  curControlPhenotype = C3NAObj_Reference$misc$phenotype
  curControl_SparccP = C3NAObj_Reference$sparCCTable
  curControl_ModuleData = C3NAObj_Reference$consensus$modules
  controlTaxaTable = C3NAObj_Reference$oriTaxTable
  controlTaxaTable$Phylum = paste0("p_", controlTaxaTable$Phylum)
  controlTaxaTable$Class = paste0("c_", controlTaxaTable$Class)
  controlTaxaTable$Order = paste0("o_", controlTaxaTable$Order)
  controlTaxaTable$Family = paste0("f_", controlTaxaTable$Family)
  controlTaxaTable$Genus = paste0("g_", controlTaxaTable$Genus)
  controlTaxaTable$Species = paste0("s_", controlTaxaTable$Species)
  
  ## Step 2. Replace WGCNA color names with Cluster ID with Standard Colors
  setTxtProgressBar(pb, 2)
  colorTable = data.frame(
    clusterID = unique(c(curDisease_ModuleData$clusterID, 
                         curControl_ModuleData$clusterID)),
    color = WGCNA::standardColors(n = length(unique(c(curDisease_ModuleData$clusterID, 
                                                      curControl_ModuleData$clusterID))))
  )
  
  curDisease_ModuleData = curDisease_ModuleData %>%
    mutate(ClusterID = paste0("Cluster ", sprintf("%02d", clusterID))) %>%
    left_join(colorTable, by = "clusterID") %>%
    dplyr::select(taxaName, ClusterID, color)
  
  curControl_ModuleData = curControl_ModuleData %>%
    mutate(ClusterID = paste0("Cluster ", sprintf("%02d", clusterID))) %>%
    left_join(colorTable, by = "clusterID") %>%
    dplyr::select(taxaName, ClusterID, color)
  
  ## Step 3. Parsing the Diseae and Control Correlations, kep only statistically significant cor > 0
  setTxtProgressBar(pb, 3)
  curDisease_SparccP_Filtered = curDisease_SparccP %>%
    dplyr::rename(from = Var1, to = Var2) %>%
    filter(fdr <= fdr & cor > corCutoff) %>%
    mutate(from = as.character(from), to = as.character(to)) %>%
    left_join(curDisease_ModuleData[, c("ClusterID", "taxaName")], by = c("from" = "taxaName")) %>%
    rename(ClusterID_From = ClusterID) %>%
    left_join(curDisease_ModuleData[, c("ClusterID", "taxaName")], by = c("to" = "taxaName")) %>%
    rename(ClusterID_To = ClusterID) %>%
    mutate(Module = ifelse(ClusterID_From == ClusterID_To, "Intra-Modular", "Inter-Modular")) %>%
    mutate(ClusterID = ifelse(ClusterID_From == ClusterID_To, ClusterID_From, NA)) %>%
    rowwise() %>%
    mutate(source = sort(c(from, to))[1],
           target = sort(c(from, to))[2]) %>%
    mutate(taxaComboName = paste0(source, "__", target)) %>%
    dplyr::select(taxaComboName, source, target, cor, ClusterID, Module) %>%
    dplyr::rename(cor_Disease = cor,
                  clusterID_Disease = ClusterID,
                  Module_Disease = Module)
  curDisease_SparccP_FilteredV2 = curDisease_SparccP_Filtered
  colnames(curDisease_SparccP_FilteredV2) = c("taxaComboName", "source", "target", "cor", "clusterID", "Module")
  curDisease_SparccP_FilteredV2$Diagnosis = curDiseasePhenotype

  
  curControl_SparccP_Filtered = curControl_SparccP %>%
    dplyr::rename(from = Var1, to = Var2) %>%
    filter(fdr <= fdr & cor > corCutoff) %>%
    mutate(from = as.character(from), to = as.character(to)) %>%
    left_join(curControl_ModuleData[, c("ClusterID", "taxaName")], by = c("from" = "taxaName")) %>%
    rename(ClusterID_From = ClusterID) %>%
    left_join(curControl_ModuleData[, c("ClusterID", "taxaName")], by = c("to" = "taxaName")) %>%
    rename(ClusterID_To = ClusterID) %>%
    mutate(Module = ifelse(ClusterID_From == ClusterID_To, "Intra-Modular", "Inter-Modular")) %>%
    mutate(ClusterID = ifelse(ClusterID_From == ClusterID_To, ClusterID_From, NA)) %>%
    rowwise() %>%
    mutate(source = sort(c(from, to))[1],
           target = sort(c(from, to))[2]) %>%
    mutate(taxaComboName = paste0(source, "__", target)) %>%
    dplyr::select(taxaComboName, source, target, cor, ClusterID, Module) %>%
    dplyr::rename(cor_Control = cor,
                  clusterID_Control = ClusterID,
                  Module_Control = Module)
  curControl_SparccP_FilteredV2 = curControl_SparccP_Filtered
  colnames(curControl_SparccP_FilteredV2) = c("taxaComboName", "source", "target", "cor", "clusterID", "Module")
  curControl_SparccP_FilteredV2$Diagnosis = curControlPhenotype
  
  sparccP_Filtered_rbind = rbind(curDisease_SparccP_FilteredV2, curControl_SparccP_FilteredV2)
  sparccP_Filtered_Combined = merge(curDisease_SparccP_Filtered, curControl_SparccP_Filtered, 
                                    by = c("taxaComboName", "source", "target"), all = TRUE)

  ## Step 4. Saving the shared taxa
  setTxtProgressBar(pb, 4)
  diseaseTaxa = data.frame(
    uniqueTaxaName = unique(unlist(curDisease_SparccP[, 1:2]))
  )
  diseaseTaxa = diseaseTaxa %>%
    filter(!(uniqueTaxaName %in% removedTaxaTable$removedTaxa)) %>%
    rowwise() %>%
    mutate(taxaLvls_Abbre = substring(uniqueTaxaName, 1, 1))

  controlTaxa = data.frame(
    uniqueTaxaName = unique(unlist(curControl_SparccP[, 1:2]))
  )
  controlTaxa = controlTaxa %>%
    filter(!(uniqueTaxaName %in% removedTaxaTable$removedTaxa)) %>%
    rowwise() %>%
    mutate(taxaLvls_Abbre = substring(uniqueTaxaName, 1, 1))
  commonTaxa = data.frame(
    uniqueTaxaName = intersect(diseaseTaxa$uniqueTaxaName, controlTaxa$uniqueTaxaName)
  )

  diseaseTaxa_only = diseaseTaxa %>% filter(!(uniqueTaxaName %in% commonTaxa$uniqueTaxaName))
  controlTaxa_only = controlTaxa %>% filter(!(uniqueTaxaName %in% commonTaxa$uniqueTaxaName))

  taxaConverstionTable = data.frame(
    taxaLvls_Abbre = c("p", "c", "o", "f", "g", "s"),
    taxaLvls_Full = c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  )

  statsTable = commonTaxa %>%
    rowwise() %>%
    mutate(taxaLvls_Abbre = substring(uniqueTaxaName, 1, 1)) %>%
    group_by(taxaLvls_Abbre) %>%
    mutate(n = 1) %>%
    left_join(taxaConverstionTable, by = "taxaLvls_Abbre") %>% ungroup() %>%
    dplyr::select(-taxaLvls_Abbre) %>% mutate(Group = "Shared Taxa") %>%
    dplyr::relocate(Group, taxaLvls_Full, n)

  diseaseTaxa_only = diseaseTaxa_only %>%
    rowwise() %>%
    mutate(taxaLvls_Abbre = substring(uniqueTaxaName, 1, 1)) %>%
    group_by(taxaLvls_Abbre) %>%
    mutate(n = 1) %>%
    left_join(taxaConverstionTable, by = "taxaLvls_Abbre") %>% ungroup() %>%
    dplyr::select(-taxaLvls_Abbre) %>% mutate(Group = "Disease-Only") %>%
    dplyr::relocate(Group, taxaLvls_Full, n)

  controlTaxa_only = controlTaxa_only %>%
    rowwise() %>%
    mutate(taxaLvls_Abbre = substring(uniqueTaxaName, 1, 1)) %>%
    group_by(taxaLvls_Abbre) %>%
    mutate(n = 1) %>%
    left_join(taxaConverstionTable, by = "taxaLvls_Abbre") %>% ungroup() %>%
    dplyr::select(-taxaLvls_Abbre) %>% mutate(Group = "Control-Only") %>%
    dplyr::relocate(Group, taxaLvls_Full, n)

  statsTable = rbind(statsTable, diseaseTaxa_only)
  statsTable = rbind(statsTable, controlTaxa_only)

  statsTable = statsTable %>%
    mutate(taxaLvls_Full = factor(taxaLvls_Full, levels = taxaConverstionTable$taxaLvls_Full),
           Group = factor(Group, levels = c("Shared Taxa", "Disease-Only", "Control-Only"))) %>%
    arrange(Group, taxaLvls_Full) %>%
    ungroup() %>%
    dplyr::select(Group, taxaLvls_Full, uniqueTaxaName, n) %>%
    dplyr::rename(`Taxonomic Levels` = taxaLvls_Full,
                  `Number of Unique Taxa` = n,
                  `Full Taxa Name` = uniqueTaxaName)

  ## Step 5. Generate the Taxa Inference Table to Determine taxa in the same phylogenetic branch
  setTxtProgressBar(pb, 5)
  taxaLvls_FullName = c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  allTaxaTable = rbind(diseaseTaxaTable, controlTaxaTable) %>% distinct()
  refTaxaTable = data.frame(
    taxaName = NULL, description = NULL
  )
  ## Phylum
  tempTable = allTaxaTable %>%
    dplyr::select(Phylum) %>%
    distinct()
  refTaxaTable_temp = data.frame(
    taxaName = tempTable$Phylum,
    description = paste0("<b>Phylum: </b>", tempTable$Phylum, "<br>"),
    Phylum = tempTable$Phylum,
    Class = NA, Order = NA, Family = NA, Genus = NA, Species = NA
  )
  refTaxaTable = rbind(refTaxaTable, refTaxaTable_temp)
  ## Class
  tempTable = allTaxaTable %>%
    dplyr::select(Phylum, Class) %>%
    distinct()
  refTaxaTable_temp = data.frame(
    taxaName = tempTable$Class,
    description = paste0("<b>Phylum: </b>", tempTable$Phylum, "<br>",
                         "<b>Class: </b>", tempTable$Class, "<br>"),
    Phylum = tempTable$Phylum, Class = tempTable$Class,
    Order = NA, Family = NA, Genus = NA, Species = NA
  )
  refTaxaTable = rbind(refTaxaTable, refTaxaTable_temp)
  
  ## Order
  tempTable = allTaxaTable %>%
    dplyr::select(Phylum, Class, Order) %>%
    distinct()
  refTaxaTable_temp = data.frame(
    taxaName = tempTable$Order,
    description = paste0("<b>Phylum: </b>", tempTable$Phylum, "<br>",
                         "<b>Class: </b>", tempTable$Class, "<br>",
                         "<b>Order: </b>", tempTable$Order, "<br>"),
    Phylum = tempTable$Phylum, Class = tempTable$Class,
    Order = tempTable$Order,
    Family = NA, Genus = NA, Species = NA
  )
  refTaxaTable = rbind(refTaxaTable, refTaxaTable_temp)
  
  
  ## Family
  tempTable = allTaxaTable %>%
    dplyr::select(Phylum, Class, Order, Family) %>%
    distinct()
  refTaxaTable_temp = data.frame(
    taxaName = tempTable$Family,
    description = paste0("<b>Phylum: </b>", tempTable$Phylum, "<br>",
                         "<b>Class: </b>", tempTable$Class, "<br>",
                         "<b>Order: </b>", tempTable$Order, "<br>",
                         "<b>Family: </b>", tempTable$Family, "<br>"),
    Phylum = tempTable$Phylum, Class = tempTable$Class,
    Order = tempTable$Order, Family = tempTable$Family,
    Genus = NA, Species = NA
  )
  refTaxaTable = rbind(refTaxaTable, refTaxaTable_temp)
  
  ## Genus
  tempTable = allTaxaTable %>%
    dplyr::select(Phylum, Class, Order, Family, Genus) %>%
    distinct()
  refTaxaTable_temp = data.frame(
    taxaName = tempTable$Genus,
    description = paste0("<b>Phylum: </b>", tempTable$Phylum, "<br>",
                         "<b>Class: </b>", tempTable$Class, "<br>",
                         "<b>Order: </b>", tempTable$Order, "<br>",
                         "<b>Family: </b>", tempTable$Family, "<br>",
                         "<b>Genus: </b>", tempTable$Genus, "<br>"),
    Phylum = tempTable$Phylum, Class = tempTable$Class,
    Order = tempTable$Order, Family = tempTable$Family,
    Genus = tempTable$Genus,
    Species = NA
  )
  refTaxaTable = rbind(refTaxaTable, refTaxaTable_temp)
  
  ## Species
  tempTable = allTaxaTable %>%
    dplyr::select(Phylum, Class, Order, Family, Genus, Species) %>%
    distinct()
  refTaxaTable_temp = data.frame(
    taxaName = tempTable$Species,
    description = paste0("<b>Phylum: </b>", tempTable$Phylum, "<br>",
                         "<b>Class: </b>", tempTable$Class, "<br>",
                         "<b>Order: </b>", tempTable$Order, "<br>",
                         "<b>Family: </b>", tempTable$Family, "<br>",
                         "<b>Genus: </b>", tempTable$Genus, "<br>",
                         "<b>Species: </b>", tempTable$Species, "<br>"),
    Phylum = tempTable$Phylum, Class = tempTable$Class,
    Order = tempTable$Order, Family = tempTable$Family,
    Genus = tempTable$Genus, Species = tempTable$Species
  )
  refTaxaTable = rbind(refTaxaTable, refTaxaTable_temp)
  refTaxaTable = refTaxaTable %>%
    filter(!(taxaName %in% removedTaxaTable$removedTaxa))

  if(verbose){
    dupData = as.data.frame(table(refTaxaTable$taxaName)) %>%
      dplyr::rename(`Taxa Name` = Var1) %>%
      filter(Freq > 1)
    cat(paste0("Uncertain/Duplicated taxa names detected: ", 
               nrow(dupData),
               ". \nThis will affect the interactive annotation display as only ", 
               "the first occurance will be used. \n"))
    print(dupData)
  }
  refTaxaTable = refTaxaTable %>%
    group_by(taxaName) %>%
    filter(row_number() == 1)
  if(verbose){
    cat(paste0("After removing the duplicated taxa names, the updated number of ",
               "post-processed taxa is: ", nrow(refTaxaTable), "\n"))
  }

  samePhyloTable = data.frame(
    from = NULL, to = NULL
  )
  for(i in seq(6)){
    for(j in seq(6)){
      if(i < j){
        tempPhyloTable = allTaxaTable[, c(taxaLvls_FullName[i], taxaLvls_FullName[j])] %>% distinct()
        colnames(tempPhyloTable) = c("from", "to")
        samePhyloTable = rbind(samePhyloTable, tempPhyloTable)
      }
    }
  }
  ## Generate the taxaComboName ordered alphebetically
  samePhyloTable = samePhyloTable %>%
    rowwise() %>%
    mutate(source = sort(c(from, to))[1],
           target = sort(c(from, to))[2]) %>%
    mutate(taxaComboName = paste0(source, "__", target)) %>%
    dplyr::select(taxaComboName) %>%
    mutate(samePhylo = TRUE)

  ## Step 6. Preservation Analysis
  setTxtProgressBar(pb, 6)
  ## Prep data
  diseaseOnlyModules = sparccP_Filtered_Combined %>%
    filter(cor_Disease >= corCutoff & Module_Disease == "Intra-Modular") %>%
    filter(!(source %in% removedTaxaTable$removedTaxa) & !(target %in% removedTaxaTable$removedTaxa)) %>%
    distinct() %>%
    left_join(samePhyloTable, by = "taxaComboName", all.x = TRUE) %>%
    mutate(samePhylo = ifelse(is.na(samePhylo), FALSE, TRUE))
  
  controlOnlyModules = sparccP_Filtered_Combined %>%
    filter(cor_Control >= corCutoff & Module_Control == "Intra-Modular") %>%
    filter(!(source %in% removedTaxaTable$removedTaxa) & !(target %in% removedTaxaTable$removedTaxa)) %>%
    distinct() %>%
    left_join(samePhyloTable, by = "taxaComboName", all.x = TRUE) %>%
    mutate(samePhylo = ifelse(is.na(samePhylo), FALSE, TRUE))
  
  disease_Adj = igraph::graph.data.frame(d = diseaseOnlyModules[, c("source", "target")],
                                         directed = FALSE)
  igraph::E(disease_Adj)$weight = diseaseOnlyModules$cor_Disease
  disease_Adj = igraph::as_adjacency_matrix(disease_Adj, type = "both", attr = "weight")
  disease_Adj = as.data.frame(as.matrix(disease_Adj))
  diag(disease_Adj) = 1
  
  control_Adj = igraph::graph.data.frame(d = controlOnlyModules[, c("source", "target", "cor_Control")],
                                         directed = FALSE)
  igraph::E(control_Adj)$weight = controlOnlyModules$cor_Control
  control_Adj = igraph::as_adjacency_matrix(control_Adj, type = "both", attr = "weight")
  control_Adj = as.data.frame(as.matrix(control_Adj))
  diag(control_Adj) = 1
  
  ## Data List
  multiAdj = list();
  multiAdj[[1]] = list(data = as.matrix(disease_Adj))
  multiAdj[[2]] = list(data = as.matrix(control_Adj))
  setLabels = c("Disease", "Control");
  names(multiAdj) = setLabels

  nodesDisease = curDisease_ModuleData[, c("taxaName", "ClusterID")] %>%
    filter(taxaName %in% colnames(disease_Adj))
  nodesControl = curControl_ModuleData[, c("taxaName", "ClusterID")] %>%
    filter(taxaName %in% colnames(control_Adj))
  
  nodesDisease_Abbrev = nodesDisease %>% select(ClusterID) 
  nodesControl_Abbrev = nodesControl %>% select(ClusterID) 
  nodesAbbre = rbind(nodesDisease_Abbrev, nodesControl_Abbrev) %>%
    distinct()
  nodesAbbre$colors = standardColors(nrow(nodesAbbre))
  if("grey" %in% nodesAbbre$colors){print("Found Grey"); greyFound = TRUE} else{greyFound = FALSE}
  nodesDisease = nodesDisease %>%
    left_join(nodesAbbre, by = "ClusterID") %>%
    column_to_rownames("taxaName")
  nodesDisease = nodesDisease[colnames(disease_Adj), ]
  nodesControl = nodesControl %>%
    left_join(nodesAbbre, by = "ClusterID") %>%
    column_to_rownames("taxaName")
  nodesControl = nodesControl[colnames(control_Adj), ]
  
  colorList = list(nodesDisease$colors, nodesControl$colors);
  names(colorList) = setLabels;
  
  nodesDiseaseV2 = nodesDisease 
  colnames(nodesDiseaseV2) = c("ClusterID_Disease", "colors_Disease")
  nodesDiseaseV2$TaxaName = rownames(nodesDisease)
  
  nodesControlV2 = nodesControl 
  colnames(nodesControlV2) = c("ClusterID_Control", "colors_Control")
  nodesControlV2$TaxaName = rownames(nodesControl)

  nodesAll = merge(nodesDiseaseV2, nodesControlV2, by = "TaxaName", all = TRUE)
  modulePreservation = suppressWarnings(.updated_modulePreservation(
    multiData = multiAdj, multiColor = colorList, 
    multiWeights = NULL,dataIsExpr = FALSE, 
    referenceNetworks = 2, testNetworks = 1, 
    networkType = "unsigned", 
    quickCor = 0,
    nPermutations = nBootstrap, # Default: 300
    verbose = verbose
  ))
  
  mp=modulePreservation
  stats= mp$preservation$observed[[1]][[1]]
  labelsX = rownames(stats)
  labelsX[labelsX=="gold"] = "orangeX"
  modColors = labelsX;
  plotMods = !(modColors %in% c("orangeX"));
  moduleSizes = stats[plotMods, 1];
  textLabels = modColors[plotMods]
  colorLabels = labelsX[plotMods];
  
  nModules = sum(plotMods);
  
  modulePlotData_Long = data.frame(
    colorType = NULL,
    varCategory = NULL,
    value = NULL,
    moduleSize = NULL,
    controlType = NULL
  )
  ## Z Summary
  modulePlotData_temp = data.frame(
    color = rownames(mp$preservation$Z[[1]][[1]])[plotMods],
    varCategory = "ZSummary",
    value = mp$preservation$Z[[1]][[1]]$Zsummary.pres[plotMods],
    moduleSize = moduleSizes,
    controlType = curControlPhenotype
  )
  modulePlotData_Long = rbind(modulePlotData_Long, modulePlotData_temp)
  ## medianRank
  modulePlotData_temp = data.frame(
    color = rownames(mp$preservation$Z[[1]][[1]])[plotMods],
    varCategory = "medianRank",
    value = mp$preservation$observed[[1]][[1]]$medianRank.pres[plotMods],
    moduleSize = moduleSizes,
    controlType = curControlPhenotype
  )
  modulePlotData_Long = rbind(modulePlotData_Long, modulePlotData_temp)
  
  modulePlotData_Wide = data.frame(
    color = rownames(mp$preservation$Z[[1]][[1]])[plotMods],
    ZSummary = mp$preservation$Z[[1]][[1]]$Zsummary.pres[plotMods],
    medianRank = mp$preservation$observed[[1]][[1]]$medianRank.pres[plotMods],
    moduleSize = moduleSizes,
    controlType = curControlPhenotype
  )
  modulePlotData_Wide$nonPreserved = ifelse((modulePlotData_Wide$ZSummary < 5 & 
                                               modulePlotData_Wide$medianRank >=8), 
                                            TRUE, FALSE)
  setTxtProgressBar(pb, 7)
  C3NAObj = list()
  C3NAObj[["curDisease_TaxaTable"]] = diseaseTaxaTable
  C3NAObj[["curControl_TaxaTable"]] = controlTaxaTable
  C3NAObj[["curDisease_OriCountTable"]] = C3NAObj_Comparison$filteredCountMatrix
  C3NAObj[["curControl_OriCountTable"]] = C3NAObj_Reference$filteredCountMatrix
  
  C3NAObj[["sparccP_Filtered_Combined"]] = sparccP_Filtered_Combined
  C3NAObj[["sparccP_Filtered_rbind"]] = sparccP_Filtered_rbind
  C3NAObj[["statsTable"]] = statsTable
  
  C3NAObj$modulePreservation$modulePlotData_Wide = modulePlotData_Wide
  C3NAObj$nodes$nodesAll = nodesAll
  C3NAObj$nodes$nodesAbbre = nodesAbbre
  C3NAObj$nodes$samePhyloTable = samePhyloTable
  C3NAObj$nodes$refTaxaTable = refTaxaTable
  
  C3NAObj$miscDisease = C3NAObj_Comparison$misc
  C3NAObj$miscControl = C3NAObj_Reference$misc
  
  C3NAObj$misc$corCut = corCutoff
  C3NAObj$misc$fdr = fdr
  
  C3NAObj$consensusDisease = C3NAObj_Comparison$consensus
  C3NAObj$consensusControl = C3NAObj_Reference$consensus
  
  C3NAObj$removedTaxaTable = removedTaxaTable
  
  return(C3NAObj)
}

