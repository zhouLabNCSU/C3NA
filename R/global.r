## Set global variables to prevent no visible binding for global variable warning
utils::globalVariables(c(
  # validatePhloseq
  "diagnoColName",
  # initiateC3NA
  "diagnosis",
  # Variables
  "BinarySimilarity", "C3NA", 
  "C3NAObj1", "C3NAObj2", "Category", "Class", 
  "ClusterID", "ClusterID_Control", "ClusterID_Disease", "ClusterID_From", 
  "ClusterID_To", "Colors", "ConsensusProp", "ControlConsensusProp", "DA", 
  "Diagnosis", "DiseaseConsensusProp", "Family", "Freq", "Full", "Taxa", 
  "Name", "Genus", "Group", "Membership", "Module", "Module_Control", 
  "Module_Disease", "Order", "Phylum", "Silouette", "Species", "TEMP", 
  "TaxaName", "Taxonomic", "Levels", "Var1", "Var2", "Var3", "ZSummary", 
  "abline", "betweenessWilcoxonP","clusterID", 
  "clusterID_Control", "clusterID_Disease", "color", "color.background", "color.border", 
  "colorRamp", "colors", "colors_Control", "colors_Disease", "cor", "cor_Control", "cor_Disease", 
  "degreeWilcoxonP", "dev.off", "errbar", "fdr", "foreach", "from", "head", "index", "inputTaxaName",
  "lines", "lm", "medianRank", "minSize", "nTaxa", "original", "par", "pdf", "phenotype", "predict", 
  "presence", "samePhylo", "samePhyloV2", "sample_data<-", "sd", "shape", "target", "taxaComboName",
  "taxaLevelAbbrev", "taxaLvls_Abbre", "taxaLvls_Full", "taxaName", "text", "title", "to", 
  "totalUniqueColors", "transitivityWilcoxonP", "uniqueTaxaName", "value", "verboseScatterplot", "wilcox.test",
  "Full Taxa Name", "Taxonomic Levels", "."
))

#' Correlation and Consensus-Based Cross-Taxonmony Network Analysis
#'
#' Correlation and Consensus-based Cross-taxonomy Network Analysis (C3NA) is a 
#' user-friendly tool for investigating compositional microbial sequencing data 
#' to identify and compare co-occurrence patterns across different taxonomic 
#' levels between two phenotypes. 
#'
#' @import methods
#' @name C3NA-package
#' @author Kuncheng Song 
#' @docType package
#' @keywords package
NA

