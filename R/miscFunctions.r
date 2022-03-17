#' Verify Phyloseq Objects before Running the C3NA Pipeline
#' 
#' @description This validation serve to ensure the data are prepared is ready for the C3NA. 
#' 1. Ensure the OTUs/ASVs table, Taxa table and sample data with a diagnosis column present. 
#' 2. The user can provide the correct diagnosis column required for the C3NA phenotype data processing. 
#' 
#' @param phyloseqObj (Required) A phyloseq object.
#' @param diagColName (Optional) Identify the dianogis/phenotype column, and create a 'diagnosis' column 
#' the column was named differently.   
#' 
#' @importFrom phyloseq sample_data phyloseq
#' @importFrom metagMisc phyloseq_filter_prevalence
#' @importFrom WGCNA moduleEigengenes signedKME qvalue 
#' @importFrom impute impute.knn 
#' @importFrom foreach %dopar% foreach 
#' @importFrom stats pnorm median 
#' @importFrom grDevices pdf dev.off colors colorRamp 
#' @importFrom graphics abline lines par text title
#' @importFrom stats cor lm predict sd wilcox.test
#' @return A phyloseq object
#' @export
#' @examples
#' data(CRC_Phyloseq)
#' curPhyloseq = validatePhloseq(phyloseqObj = CRC_Phyloseq)
#' 
validatePhloseq <- function(phyloseqObj = phyloseqObj, 
                            diagColName = NA) {
  if(!class(phyloseqObj) == "phyloseq"){
    stop(paste0("The object is not a phyloseq object, please correctly import the ",
                "OTUs/ASVs table, taxa table and sample data as phyloseq obejcts."))
  }
  if(!("diagnosis" %in% colnames(sample_data(phyloseqObj))) & is.na(diagColName)){
    stop(paste0("The diagnosis column is missing in from the sample_data() of the ",
                "phyloseq object. You can provide the correct phenotype column ",
                "name using the arg 'diagColName = phenotype.'"))
  } else if (!("diagnosis" %in% colnames(sample_data(phyloseqObj))) & is.na(diagColName)){
    if(diagnoColName %in% colnames(sample_data(phyloseqObj))){
      sample_data(phyloseqObj)$diagnosis = sample_data(phyloseqObj)[[diagnoColName]]
    } else {
      stop(paste0("The provided diagnosis column is not in the phyloseq object."))
    }
  }
  return(phyloseqObj)
} 

 
#' C3NA Required Package Check
#'
#' @description Check if required packages are installed and its current version
#' 
#' @importFrom utils installed.packages packageVersion
#' @export
#' @examples
#' C3NAPackageCheck()
C3NAPackageCheck <- function() {
  ## CRAN Repo
  .cran_packages <- c("stringr", "tibble", "dplyr", "readr", "stats", "visNetwork", 
                      "igraph", "scales", "magrittr", "tidyr", "randomcoloR", 
                      "ggpubr", "ggplot2",  "plotly", "colorspace", 
                      "reshape2", "shiny", "shinyWidgets", "shinydashboard", 
                      "shinyjs", "reactable", "DT", "pheatmap", "dynamicTreeCut", 
                      "cluster", "WGCNA", "RColorBrewer", "base64enc", "htmlwidgets",
                       "foreach")
  packageResults = suppressWarnings(suppressMessages(as.data.frame(sapply(.cran_packages, require, character.only = TRUE))))
  colnames(packageResults) = "Found"
  packageResults$Library = "CRAN"
  packageResults$Version = NA
  packageResults$PackageName = rownames(packageResults)
  for(i in seq(nrow(packageResults))){
    if(packageResults$Found[i] == TRUE){
      packageResults$Version[i] = as.character(packageVersion(packageResults$PackageName[i]))
    }
  }
  
  ## Bioconductor Repo
  .bioc_packages <- c("phyloseq", "qvalue", "impute")
  packageResults_temp = suppressWarnings(suppressMessages(as.data.frame(sapply(.bioc_packages, require, character.only = TRUE))))
  colnames(packageResults_temp) = "Found"
  packageResults_temp$Library = "Bioconductor"
  packageResults_temp$Version = NA
  packageResults_temp$PackageName = rownames(packageResults_temp)
  for(i in seq(nrow(packageResults_temp))){
    if(packageResults_temp$Found[i] == TRUE){
      packageResults_temp$Version[i] = as.character(packageVersion(packageResults_temp$PackageName[i]))
    }
  }
  packageResults = rbind(packageResults, packageResults_temp)
  
  ## GitHub Repo
  .git_packages <- c("SpiecEasi", "metagMisc")
  .inst <- .git_packages %in% installed.packages()
  packageResults_temp = suppressWarnings(suppressMessages(as.data.frame(sapply(.git_packages, require, character.only = TRUE))))
  colnames(packageResults_temp) = "Found"
  packageResults_temp$Library = "GitHub"
  packageResults_temp$Version = NA
  packageResults_temp$PackageName = rownames(packageResults_temp)
  for(i in seq(nrow(packageResults_temp))){
    if(packageResults_temp$Found[i] == TRUE){
      packageResults_temp$Version[i] = as.character(packageVersion(packageResults_temp$PackageName[i]))
    }
  }
  packageResults = rbind(packageResults, packageResults_temp)
  rownames(packageResults) = seq(nrow(packageResults))
  return(packageResults)
}

#' Extract C3NA Results
#'
#' @description The C3NAObj after comparePhenotypes() command is quiet large due to storing accessory data for Shiny plotting. 
#' This command is developed to extract the following information: 1. Taxa-Taxa Correlations that passed pre-set 
#' BH-adjusted p-value cut-off, here we did NOT filter by any correlation threshold. 
#' 2. All the node data for both phenotypes, it will contain the differential abundance and C3NA
#' results if avaliable. This can then be used for other network plotting software or analysis. 
#' 
#' @param twoPhenoC3NAObj C3NA objects after running comparePhenotypes()
#' @export
#' 
extractResults <- function(twoPhenoC3NAObj = twoPhenoC3NAObj){
  if(is.null(twoPhenoC3NAObj$misc$corCut)){
    stop("Please run the comparePhenotypes() before running the interactive Shiny application. ")
  }
  returnData = list()
  ## Extract the taxa-taxa correlations
  wideCorTable = twoPhenoC3NAObj$sparccP_Filtered_Combined
  
  ## Extract Nodes and add DA/C3NA if avaliable
  curNodes = twoPhenoC3NAObj$nodes$nodesAll
  if(!is.null(twoPhenoC3NAObj$C3NA_Wilcoxon)){
    curNodes = curNodes %>%
      left_join(twoPhenoC3NAObj$C3NA_Wilcoxon[, c("TaxaName", "C3NA")], 
                by = "TaxaName", all.x = TRUE) %>%
      mutate(C3NA = ifelse(is.na(C3NA), FALSE, C3NA))
  }
  if(!is.null(twoPhenoC3NAObj$DA)){
    curNodes = curNodes %>%
      left_join(twoPhenoC3NAObj$DA, 
                by = "TaxaName", all.x = TRUE) %>%
      mutate(C3NA = ifelse(is.na(C3NA), FALSE, C3NA))
  }
  
  ### Replace names with the correct phenotypes
  colnames(wideCorTable) = gsub("Disease", twoPhenoC3NAObj$miscDisease$phenotype, colnames(wideCorTable))
  colnames(curNodes) = gsub("Disease", twoPhenoC3NAObj$miscDisease$phenotype, colnames(curNodes))
  
  colnames(wideCorTable) = gsub("Control", twoPhenoC3NAObj$miscControl$phenotype, colnames(wideCorTable))
  colnames(curNodes) = gsub("Control", twoPhenoC3NAObj$miscControl$phenotype, colnames(curNodes))

  returnData[["taxaTaxaCorrelations"]] = wideCorTable
  returnData[["taxaTable"]] = curNodes
  
  return(returnData)
}

















