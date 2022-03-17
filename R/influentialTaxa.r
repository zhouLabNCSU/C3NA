#' Adding Influential Taxa Evaluation to the C3NA Object
#' 
#' @description Using Paired Wilcoxon Test to evaluate each taxa based on the degree, 
#' betweenness and transitivity changes. 
#' 
#' @param C3NAObj (Required) The C3NAObj after the comparePhenotypes() function
#' @param C3NAFDR (Optional) This allows to change the BH-Adjuted p-value for the paired 
#' Wilcoxon test for degree, betweenness and transitivity
#' 
#' @importFrom igraph degree betweenness transitivity delete_vertices components
#' @export
#' 
getInfluentialTaxa <- function(C3NAObj = C3NAObj, C3NAFDR = 0.05){
  ## Extract the correlation threshold
  corCut = C3NAObj$misc$corCut

  ## Obtain the Edge Data 
  edgeData = C3NAObj$sparccP_Filtered_Combined %>%
    filter((cor_Disease >= corCut) | (cor_Control >= corCut)) %>%
    mutate(
      edgeColor = case_when(
        (Module_Control=="Inter-Modular" | is.na(Module_Control)) & 
          Module_Disease=="Intra-Modular" ~ "#77BBEB70",
        (Module_Disease=="Inter-Modular" | is.na(Module_Disease)) & 
          Module_Control=="Intra-Modular" ~ "#EB787470",
        (Module_Disease=="Intra-Modular" & Module_Control=="Intra-Modular") ~ "#CCCFCB70"
      )
    ) %>%
    filter((source %in% C3NAObj$nodes$nodesAll$TaxaName)) %>%
    filter((target %in% C3NAObj$nodes$nodesAll$TaxaName))
  
  ## Obtain the Node Data 
  nodesData = C3NAObj$nodes$nodesAll %>%
    filter(TaxaName %in% unique(unlist(edgeData[, c("source", "target")])))
  
  ## Paired Wilcoxon  
  wilcoxonTest = data.frame(
    TaxaName = NULL, degreeWilcoxonP = NULL,
    betweenessWilcoxonP = NULL, transitivityWilcoxonP = NULL
  )
  pb = txtProgressBar(min = 0, max = length(unique(nodesData$TaxaName)), initial = 0, style = 3)
  counter = 0
  for(t in unique(nodesData$TaxaName)){
    setTxtProgressBar(pb, counter); counter = counter + 1;
    curTaxa = t
    diseaseClutser = subset(nodesData, TaxaName == curTaxa)$ClusterID_Disease
    controlClutser = subset(nodesData, TaxaName == curTaxa)$ClusterID_Control
    
    curNet = edgeData %>%
      filter(clusterID_Disease %in% diseaseClutser | clusterID_Control %in% controlClutser) %>%
      filter((cor_Disease >= 0.2) | (cor_Control >= 0.2)) 
    curG = graph.data.frame(curNet[, c("source", "target")],
                            directed = FALSE)
    ## Obtain the subgraph with the taxa in it
    comp_g <- components(curG)
    componentTable = as.data.frame(comp_g$membership)
    colnames(componentTable) = "Membership"
    componentTable$TaxaName = rownames(componentTable)
    
    ## Selected Membership
    curMembership = subset(componentTable, TaxaName == curTaxa)$Membership
    if(length(curMembership) >= 1){
      removedSubgraph = subset(componentTable, Membership != curMembership)
      curG = delete_vertices(curG, removedSubgraph$TaxaName)
      curG_Removed = delete_vertices(curG, curTaxa)
      curG_Data = data.frame( 
        TaxaName = names((igraph::degree(curG, normalized = TRUE))),
        degree = as.numeric(igraph::degree(curG, normalized = TRUE)),
        betweeness = as.numeric(betweenness(curG, directed = FALSE, normalized = TRUE)),
        transitivity = as.numeric(transitivity(curG, type = "local", isolates = "zero")),
        type = "Original"
      )
      curG_Data_V2 = curG_Data %>%
        filter(TaxaName != curTaxa)
      curG_Removed_Data = data.frame( 
        TaxaName = names((igraph::degree(curG_Removed, normalized = TRUE))),
        degree = as.numeric(igraph::degree(curG_Removed, normalized = TRUE)),
        betweeness = as.numeric(betweenness(curG_Removed, directed = FALSE, normalized = TRUE)),
        transitivity = as.numeric(transitivity(curG_Removed, type = "local", isolates = "zero")),
        type = "Removed"
      )
      
      fullData = rbind(curG_Data_V2, curG_Removed_Data)
      fullData[is.na(fullData)] = 0
      degreeResult = suppressWarnings(wilcox.test(degree ~ type, data = fullData, paired = TRUE))
      betweenessResult = suppressWarnings(wilcox.test(betweeness ~ type, data = fullData, paired = TRUE))
      transitivityResult = suppressWarnings(wilcox.test(transitivity ~ type, data = fullData, paired = TRUE))
      wilcoxonTest_temp = data.frame(
        TaxaName = curTaxa, 
        degreeWilcoxonP = degreeResult$p.value,
        betweenessWilcoxonP = betweenessResult$p.value,
        transitivityWilcoxonP = transitivityResult$p.value
      )
      wilcoxonTest = rbind(wilcoxonTest, wilcoxonTest_temp)  
    }
  } 
  
  ## BH-Adjusted P and C3NA 
  wilcoxonTest = wilcoxonTest %>%
    replace(is.na(.), 1) %>%
    ## Remove non-contributing nodes/taxa
    filter((degreeWilcoxonP + betweenessWilcoxonP + transitivityWilcoxonP) != 3) 
  wilcoxonTest$degreeWilcoxonPAdj = p.adjust(wilcoxonTest$degreeWilcoxonP, method = "BH")
  wilcoxonTest$betweenessWilcoxonPAdj = p.adjust(wilcoxonTest$betweenessWilcoxonP, method = "BH")
  wilcoxonTest$transitivityWilcoxonPAdj = p.adjust(wilcoxonTest$transitivityWilcoxonP, method = "BH")
  wilcoxonTest$C3NA = ifelse((wilcoxonTest$degreeWilcoxonPAdj <= C3NAFDR & 
                                wilcoxonTest$betweenessWilcoxonPAdj <= C3NAFDR & 
                                wilcoxonTest$transitivityWilcoxonPAdj <= C3NAFDR) 
                             , TRUE, FALSE)
  ## return
  C3NAObj[["C3NA_Wilcoxon"]] = wilcoxonTest
  return(C3NAObj)
}

