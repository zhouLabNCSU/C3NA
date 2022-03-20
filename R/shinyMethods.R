#' Shiny Application for Optimal Module Evaluations
#' 
#' @description Shiny application that includes two panels. On the first panel, the user can conduct a step-by-step investigation of the consensus-based module selection, and the results will be summarized into an R command in which the user copies and pastes to the R console. The second panel includes a download window for extracting the module pattern, silhouette plot, consensus, and correlation plots. 
#' 
#' @param C3NAObj (Required) The C3NAObj with a single phenotype and after the initiateC3NA function. 
#' @param colorSeed (Optional) The seed set for reproducible random colors for the bar plot fill. Default: 9.
#' 
#' @importFrom shinydashboard box tabBox
#' @importFrom shinyjs useShinyjs
#' @importFrom tibble as_tibble
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr case_when 
#' @rawNamespace import(shiny, except = c(renderDataTable, dataTableOutput))
#' @import ggplot2
#' @import pheatmap
#' @importFrom shinyWidgets setShadow actionBttn useShinydashboard awesomeCheckboxGroup updateAwesomeCheckboxGroup radioGroupButtons downloadBttn
#' @export
#' 

moduleEvals = function(C3NAObj = C3NAObj, colorSeed = 9){
  ## Extract the needed data and basic processing
  module_df = C3NAObj$moduleData
  for(minS in unique(module_df$minSize)){
    module_df_temp = module_df %>% dplyr::filter(minSize == minS)
    if(minS == 3){
      moduleStat =
        module_df_temp %>%
        group_by(phenotype, minSize, colors) %>%
        dplyr::summarise(nTaxa = n(), .groups = "drop") %>%
        mutate(totalNTaxa = length(unique(module_df_temp$taxaID)),
               totalUniqueColors = length(unique(colors)))
    }else {
      moduleStat_temp =
        module_df_temp %>%
        group_by(phenotype, minSize, colors) %>%
        dplyr::summarise(nTaxa = n(), .groups = "drop") %>%
        mutate(totalNTaxa = length(unique(module_df_temp$taxaID)),
               totalUniqueColors = length(unique(colors)))
      moduleStat = rbind(moduleStat, moduleStat_temp)
    }
  }
  
  # Check for duplicated patterns
  moduleCheck = moduleStat %>% ungroup() %>%
    dplyr::select(minSize, colors, nTaxa) %>% distinct() %>%
    pivot_wider(names_from = colors, values_from = nTaxa) %>% mutate(minSize = as.numeric(minSize)) %>%
    arrange(minSize)
  moduleCheck2 = moduleCheck %>%
    group_by(across(-c(minSize)), .drop = F) %>% mutate(n = n())
  moduleCheckRemoveDup = moduleCheck %>%
    group_by(across(-c(minSize)), .drop = F) %>% filter(row_number() == 1) %>%
    pull(minSize)

  set.seed(colorSeed)
  colorTable = data.frame(
    colors = unique(moduleStat$colors),
    randomColors = randomcoloR::distinctColorPalette(k = length(unique(moduleStat$colors)))
  )
  oriMaxClusters = max(moduleStat$totalUniqueColors)
  moduleStatV2 = moduleStat %>%
    left_join(colorTable, by = "colors")
  grayScale = moduleStatV2 %>%
    filter(!(minSize %in% moduleCheckRemoveDup)) %>%
    mutate(randomColorsGrayScale = "#000000",
           colors = "Removed")
  maxSizeCluster = moduleStatV2 %>% dplyr::select(minSize, totalUniqueColors) %>% distinct() %>% 
    filter(totalUniqueColors <= 10) %>% pull(minSize)
  
  shinyApp(
    ui = navbarPage(
      title = "C3NA", 
      # theme = bs_theme
      tags$head(
        tags$style(
          HTML(
            ".checkbox-inline { 
                    margin-left: 0px;
                    margin-right: 10px;}
             .checkbox-inline+.checkbox-inline {
                        margin-left: 0px;
                        margin-right: 10px;}
            .shiny-notification {
                  height: 100px;
                  width: 800px;
                  position:fixed;
                  top: calc(50% - 50px);
                  left: calc(50% - 400px);}
            .btn-custom-status-downloadPlots:hover, .btn-custom-status-downloadPlots:active, 
            .btn-custom-status-downloadPlots:focus, .btn-custom-status-downloadPlots:visited, 
            .btn-custom-status-downloadPlots.active, .open .dropdown-toggle.btn-custom-status-downloadPlots { 
              color: #ffffff;
              background-color: #318dc3;
              border-color: #07486f;
              box-shadow: inset 4px 4px 4px 0 #104667,
                          inset -4px -4px 4px 0 #5ebdf7;
            }
            .btn-custom-status-downloadPlots { 
              color: black; 
              background-color: #88bcdc ; 
              border-color: #07486f;
              box-shadow: 4px 4px 4px 0 rgba(0, 0, 0, 0.2),
                          -4px -4px 4px 0 rgba(255, 255, 255, 0.5),
                          inset 3px 3px 3px 0 #bdd5e4,
                          inset -3px -3px 3px 0 #0c3a56;
            } 
            "
          )
        ) 
      ),
      #####  Panel 1 #####
      tabPanel(
        ## First Panel
        title = HTML("<b>Panel 1: Consensus-Based Evaluations</b>"),
        value = "panel1",
        useShinyjs(),
        useShinydashboard(),
        setShadow(class = "box"),
        
        ##### Section 2 Taxa Statistics #####
        ##### Section 2.1 Pattern Plot #####
        fluidRow(
          box(
            width = 12,
            title = HTML("Pattern Visualization", "<font size='3'>", 
                    as.character(actionLink(inputId = "MS_PatternVis", 
                                            label = "", 
                                            icon = icon("question-circle"))), "</font>"), 
            status = "warning",
            solidHeader = TRUE, 
            shiny::plotOutput("uniquePatterns")
          )
        ),        
        ##### Section 2.2 Pattern Selector #####
        fluidRow(
          box(
            width = 6,
            title = HTML("Select the Unique Patterns from Differnt Minimal Number of Taxa per Module", "<font size='3'>", 
                         as.character(actionLink(inputId = "MS_UniquePatternSelector", 
                                                 label = "", 
                                                 icon = icon("question-circle"))), "</font>"), 
            status = "danger",
            solidHeader = TRUE, 
            uiOutput("patternSelector"),
            actionBttn(
              inputId = "uniquePatternsConfirm",
              label = "Confirm Selected Patterns",
              style = "unite", 
              color = "danger"
            )
          ),
          box(
            width = 3,
            title = HTML("Select the Optimal Number of Modules", "<font size='3'>", 
                         as.character(actionLink(inputId = "MS_OptimalModuleSelector", 
                                                 label = "", 
                                                 icon = icon("question-circle"))), "</font>"), 
            status = "primary",
            solidHeader = TRUE, 
            conditionalPanel(
              condition = "input.uniquePatternsConfirm==0",
              h6("Please Calculate the Silhouette Plot on the Left.")
            ),
            conditionalPanel(
              condition = "input.uniquePatternsConfirm>0",
              h6("Manually adjust the optimal number of modules based on the Silhouette, Consensus and Correlation Plot. "),
              uiOutput("optNModules"),
              actionBttn(
                inputId = "updateClusters",
                label = "Confirm Selected Patterns",
                style = "unite", 
                color = "primary"
              )
            )
          ),
          box(
            width = 3,
            title = HTML("Confirm the Optimal Clusters", "<font size='3'>", 
                         as.character(actionLink(inputId = "MS_OptimalModuleCodeGenerator", 
                                                 label = "", 
                                                 icon = icon("question-circle"))), "</font>"), 
            status = "success",
            solidHeader = TRUE, 
            conditionalPanel(
              condition = "input.updateClusters==0",
              h6("Please Confirm the Optimanl Number of Modules on the Left.")
            ),
            conditionalPanel(
              condition = "input.updateClusters>0",
              h6("Please run the following code after confirming the optimal number of modules."),
              h6("Make sure you replace the variables. "),
              verbatimTextOutput("optClusterCode")
            ),
          )
        ),
        ##### Section 2.3 Tabset Panel #####
        conditionalPanel(
          condition = "input.uniquePatternsConfirm>0",
          fluidRow(
            box(
              width = 12,
              title = HTML("Consensus Evaluation", "<font size='3'>", 
                           as.character(actionLink(inputId = "MS_ConsensusEval", 
                                                   label = "", 
                                                   icon = icon("question-circle"))), "</font>"), 
              status = "info",
              solidHeader = TRUE, 
              tabBox(
                title = "", 
                width = 12,
                tabPanel(title = "Silouette", 
                         shiny::plotOutput("silouettePlot", height = "600px")
                ),
                tabPanel(title = "Consensus Matrix", 
                         shiny::plotOutput("consensusPlot", height = "1600px")
                ),
                tabPanel(title = "Correlation Matrix", 
                         shiny::plotOutput("correlationPlot", height = "1600px")
                )
              )
            )
          )
        )

      ),
      tabPanel(
        ## First Panel
        title = HTML("<b>Panel 2: Download Panel </b>"),
        value = "panel2",
        fluidRow(
          box( ### 5.1 Download Plot List and Button ### 
            width = 6, 
            title = HTML("Download Plots", "<font size='3'>", 
                         as.character(actionLink(inputId = "MS_DownloadPanel", 
                                                 label = "", 
                                                 icon = icon("question-circle"))), "</font>"), 
            status = "primary", solidHeader = TRUE, 
            uiOutput("plotSelector"), 
            hr(),
            ## Conditional Panels 
            conditionalPanel(
              condition = "input.plotRadioButtons=='PDF'|input.plotRadioButtons=='PNG'",
              fluidRow(
                column(
                  width = 6,
                  textInput(inputId = "heightInput",
                            label = "Height in Inch",
                            placeholder =  "Height in Inch",
                            value = "8"),
                ),
                column(
                  width = 6,
                  textInput(inputId = "widthInput",
                            label = "Width in Inch",
                            placeholder =  "Width in Inch",
                            value = "10"),
                )
              ),
              h6(HTML("Please select the heigth and width for the plot, unit: Inch. <br>
                      For both Consensus and Correlation Plot, the dimension should started around 
                      Height: 20, Width: 20, to avoid overlapping of the taxa names. 
                             ")), # The default dpi is set to 300.
              hr()
            ),
            radioGroupButtons(
              inputId = "plotRadioButtons",
              label = "Select a format for the plot",
              choices = c("PDF", "PNG"),
              justified = TRUE,
              checkIcon = list(
                yes = icon("ok", 
                           lib = "glyphicon")),
              status = "custom-status-downloadPlots"
            ), hr(),
            fluidRow(
              column(
                12,
                downloadBttn(
                  outputId = "downloadPlots",
                  label = "Download the Selected Plot", 
                  style = "gradient",
                  color = "primary"
                )
              )
            )
          )
        )
      )
    ),
    
    ##### Server #####
    server = function(input, output, session) {
      # Reactive values
      values <- reactiveValues(
        selectedMinClust = NULL, 
        selectedPatterns = NULL
      )
      plots = reactiveValues()
      
      
      ##### 2.1 Pattern UI ######
      output$uniquePatterns <- renderPlot({
        moduleCheckRemoveDup = moduleCheckRemoveDup[moduleCheckRemoveDup <= maxSizeCluster[1]]
        output$patternSelector <- renderUI({
          awesomeCheckboxGroup(
            inputId = "patternSelector",
            label = "Default selection includes all unique patterns with ≥ ten unique modules", 
            choices = formatC(C3NAObj$misc$minModuleSize:C3NAObj$misc$maxModuleSize,flag=0,width=2), 
            selected = formatC(moduleCheckRemoveDup,flag=0,width=2),
            inline = TRUE, 
            status = "danger"
          )
        })
        
        modulePatternPlot <- 
          ggplot() +
            geom_bar(data = moduleStatV2, aes(fill = colors, x = minSize, y = nTaxa),
                     position = "fill", stat = "identity") + 
            geom_bar(data = grayScale, aes(fill = colors, x = minSize, y = nTaxa),
                     alpha = 0.5,
                     position = "fill", stat = "identity", show.legend = FALSE) +
            geom_text(data = moduleStatV2 %>% dplyr::select(minSize, totalUniqueColors), 
                      aes(x = minSize, y = 1.05, label = totalUniqueColors)) +
            geom_text(data = moduleStatV2, aes(x = 1, y = 1.05, label = paste0("#Modules:")),  
                      color = "black") +
            scale_fill_manual(breaks = c(colorTable$colors, "Removed"),
                              values = c(colorTable$randomColors, "#000000")) + theme_bw() +
            theme(legend.position = "none") + 
            labs(x = "Minimal # Taxa per Module", y ="Relative Abundance of Number of Taxa per Module",
                 caption = "Darkend bars represents duplicated patterns") + 
            scale_x_continuous(breaks= seq(C3NAObj$misc$minModuleSize,C3NAObj$misc$maxModuleSize,2), 
                               labels = seq(C3NAObj$misc$minModuleSize,C3NAObj$misc$maxModuleSize,2))
        plots[["modulePatternPlot"]] = modulePatternPlot
        modulePatternPlot
      })
      
      
      ##### 2.1.1.1 Pattern Selector ######
      output$patternSelector <- renderUI({
        awesomeCheckboxGroup(
          inputId = "patternSelector",
          label = "", 
          choices = formatC(C3NAObj$misc$minModuleSize:C3NAObj$misc$maxModuleSize,flag=0,width=2), 
          selected = formatC(moduleCheckRemoveDup,flag=0,width=2),
          inline = TRUE, 
          status = "danger"
        )
      })
      
      ##### 2.1.1.2 Optimal number of Cluster ######
      output$optNModules <- renderUI({
        sliderInput(
          inputId = "optNModules", label = "", 
          min = C3NAObj$misc$minModuleSize, max = C3NAObj$misc$maxModuleSize, 
          step = 1, value = 10
        )
      })
      
      ##### 2.1.1.3 Verbatim Text output ######
      output$optClusterCode <- renderText({ 
        if(!is.null(values[["selectedMinClust"]])){
          print(paste0("newC3NAObj = getOptMods(C3NAObj = oldC3NAObj,\n", 
                       "selectedPatterns = c(", paste0(values[["curModuleSizes"]], collapse = ","),"),\n", 
                       "nModules = ", values[["selectedMinClust"]], ")")) 
        }
      })  

      ##### 2.1.1 Update the Selected Patterns ######
      observeEvent(input$uniquePatternsConfirm, {
        values[["selectedPatterns"]] = as.numeric(input$patternSelector)
      })
      
      ##### 2.2 Silhouette Plot ######
      output$silouettePlot <- renderPlot({
        if(!is.null(values[["selectedPatterns"]])){
          withProgress(message = 'Evaluating Modular Information. Current procedure: ', value = 0, {
            ## Total Parts
            n = 3 + length(values[["selectedPatterns"]])
            counter = 1
            moduleCheckRemoveDup = as.numeric(values[["selectedPatterns"]])
            
            curminSilWitdth = 0.5
            incProgress(1/n, detail = paste("Part", counter, "Calculate the Consensus Matrix based on Selected Patterns"))
            ## Going over each of the minSizes and retrieve cluster module connection information
            curModuleSizes = moduleCheckRemoveDup[moduleCheckRemoveDup <= maxSizeCluster[1]]
            moduleData = module_df %>%
              filter(minSize %in% as.numeric(curModuleSizes)) %>%
              arrange(minSize)
            uniqueSizes = sort(unique(moduleData$minSize))
            values[["curModuleSizes"]]= curModuleSizes
            for(minN in curModuleSizes){
              incProgress(1/n, detail = paste("Part", counter, "Calculate the Consensus Matrix based on MinSize: ", minN))
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
              counter = counter + 1
            }
            ## Expands the consensMatrix and calculate the average consensus
            consensusMatrix_wide = consensusMatrix %>%
              pivot_wider(names_from = Category, values_from = BinarySimilarity)
            consensusMatrix_wide$rowMeans = rowMeans(consensusMatrix_wide[, 3:ncol(consensusMatrix_wide)])
            consensusMatrix_Final = consensusMatrix_wide %>% dplyr::select(from, to, rowMeans)
            consensusMatrix_Final_wide = consensusMatrix_Final %>%
              pivot_wider(names_from = to, values_from = rowMeans)
            consensusMatrix_Final_wide <- consensusMatrix_Final_wide %>% column_to_rownames("from")
            values[["consensusMatrix_Final_wide"]] = consensusMatrix_Final_wide
            incProgress(1/n, detail = paste("Part", counter+1, "Generating the Silhouette Plot Data"))
            ## Calculation of the optimal cluster number based on the consensus matrix
            curHclust = cluster::agnes(x = (1-consensusMatrix_Final_wide),
                                       method = "complete", metric = "euclidean")
            values[["curHclust"]] = curHclust
            # Calculate the evals with silhouette scores
            ks = C3NAObj$misc$minModuleSize:C3NAObj$misc$maxModuleSize
            withinClusterPropZeros = rep(NA, length(ks))
            names(withinClusterPropZeros) = ks
            prop_zeros = rep(NA, length(ks))
            names(prop_zeros) = ks
            silWidths = c()
            for(k in seq_along(ks)){
              sil <- cluster::silhouette(cutree(curHclust, ks[k]), 1-consensusMatrix_Final_wide)
              silSum = summary(sil)
              avgSilWidth = silSum$avg.width
              silWidths = c(silWidths, avgSilWidth)
              
              curTree = cutree(as.hclust(curHclust), k = ks[k])
              curTreeData =  as.data.frame(curTree)
              curTreeData$Taxa = rownames(curTreeData)
              entries = c()
              nonZeros = 0
              totalWithinClusterCounts = 0
              for(x in seq(ks[k])){
                curTaxa = subset(curTreeData, curTree == x)$Taxa
                tempM = consensusMatrix_Final_wide[curTaxa, curTaxa]
                entries = c(entries, tempM[lower.tri(tempM)])
                nonZeros = nonZeros + sum(tempM!=0)
                totalWithinClusterCounts = totalWithinClusterCounts + length(as.matrix(tempM))
              }
              # calculate proportion of zeroes in M
              nTotal = length(as.matrix(consensusMatrix_Final_wide))
              withinClusterPropZeros[k] = (totalWithinClusterCounts - nonZeros)/totalWithinClusterCounts
              prop_zeros[k] = (nTotal - nonZeros)/nTotal
            }
            propZerosCut = min(as.integer(names((withinClusterPropZeros[which(withinClusterPropZeros<0.1)]))))
            
            incProgress(1/n, detail = paste("Part", counter+2, "Plotting.. "))
            # Optimal Cluster
            for(i in seq(length(silWidths)-1)){
              j = i+1
              if(silWidths[i] >= curminSilWitdth & silWidths[j] < silWidths[i] & i >= propZerosCut){
                localMax = silWidths[i]
                minClust = i+2 # i value started at 1 and k started at 3
                break
              }
            }
            values[["selectedMinClust"]] = minClust
            updateSliderInput(
              session = session, 
              inputId = "optNModules", 
              value = minClust
            )
            
            clusterInfo = data.frame(
              index = C3NAObj$misc$minModuleSize:C3NAObj$misc$maxModuleSize,
              Silouette = silWidths,
              withinClusterPropZeros = withinClusterPropZeros
            ) 
            silPlot = 
              suppressWarnings(suppressMessages(
                ggplot(clusterInfo, aes(x = index)) +
                  geom_hline(yintercept = curminSilWitdth, linetype = "dashed", color = "#add8e6", alpha = 0.7,  size = 1.5) +
                  geom_hline(yintercept = 0.1, linetype = "dotted", color = "#d94c51", alpha = 0.7,  size = 1.5) +
                  ## Silhouette - Blue
                  geom_line(aes(y = Silouette), color = "#2E40DB") +
                  geom_point(aes(y = Silouette), color = "#2736BA") +
                  ## Prop Zeros - Red
                  geom_line(aes(y = withinClusterPropZeros), color = "#DB383E") +
                  geom_point(aes(y = withinClusterPropZeros), color = "#C23237") +
                  theme_bw() +
                  labs(x = "Number of Clusters") +
                  geom_vline(xintercept = minClust, linetype = "solid", alpha = 0.3,
                             color = "#2E994C", size=1.5) +
                  geom_text(aes(x = 35, y = curminSilWitdth-0.05,
                                label = paste0("Minimal Silhouette Cutoff: ", curminSilWitdth)),
                            color = "#add8e6") +
                  geom_text(aes(x = 35, y = 0.07,
                                label = paste0("10% Proportion of Zero Threshold")),
                            color = "#d94c51") +
                  
                  labs(title = paste0("Initial Estimate for the Optimal Number of Modules: ", minClust),
                       caption = paste0("Consensus Matrix Generated based on Selected Module Patterns")) +
                  ylim(0, 1) +
                  scale_y_continuous(
                    name = "Average Silhouette Width",
                    sec.axis = sec_axis(~., name = "Within Consensus k's Module Proprotion of Zeros",
                                        labels = function(b) { paste0(round(b * 100, 0), "%")})
                  ) +
                  theme(axis.title.y.left = element_text(colour = "#2736BA",
                                                         face = "bold"),
                        axis.title.y.right = element_text(colour = "#C23237",
                                                          face = "bold"),
                        axis.text.y.left = element_text(color = "#2736BA"),
                        axis.text.y.right = element_text(color = "#C23237"),
                        plot.title = element_text(color = "#2E994C", face = "italic")) +
                  scale_x_continuous(breaks= seq(5,40,5),
                                     labels = seq(5,40,5))
              ))
            plots[["silPlot"]] = silPlot
            silPlot
          })  
        }
      })
      

      ##### 2.1.1.4 MultiListener ######
      consensusListeners <- reactive({
        list(input$uniquePatternsConfirm,
             input$updateClusters)
      })
      
      observeEvent(input$optNModules, {
        values[["selectedMinClust"]] = input$optNModules
      })
      
      observeEvent(consensusListeners(), {
        ##### 2.1.2 Consensus Matrix ######
        output$consensusPlot <- renderPlot({
          withProgress(message = 'Generating Consensus Plot, please wait...', value = 0, {
            n=3
            incProgress(1/n, detail = paste("Part", 1, "calculate Consensus Data"))
            if(!is.null(values[["selectedMinClust"]])){
              minClust = values[["selectedMinClust"]]
              curHclust = values[["curHclust"]]
              consensusMatrix_Final_wide = values[["consensusMatrix_Final_wide"]]

              # Obtain the clustering information for the taxa
              curTree = cutree(as.hclust(curHclust), k = minClust)
              ## Generate Pheatmap
              my_taxa_colData_NbClust <- as.data.frame(curTree)
              colnames(my_taxa_colData_NbClust)[1] = "Cluster"
              my_taxa_colData_NbClust$Cluster = as.character(paste0("Cluster ", 
                                                                    formatC(my_taxa_colData_NbClust$Cluster,flag=0,width=2)))
              my_taxa_colData_NbClust$TaxaName = rownames(my_taxa_colData_NbClust)
              my_taxa_colData_NbClust = my_taxa_colData_NbClust %>% as_tibble() %>%
                mutate(taxaLevelAbbrev = sapply(strsplit(TaxaName, "_"), "[", 1)) %>%
                mutate(taxaLevel = case_when(
                  taxaLevelAbbrev =="p" ~ "Phylum",
                  taxaLevelAbbrev =="c" ~ "Class",
                  taxaLevelAbbrev =="o" ~ "Order",
                  taxaLevelAbbrev =="f" ~ "Family",
                  taxaLevelAbbrev =="g" ~ "Genus", 
                  taxaLevelAbbrev =="s" ~ "Species",
                )) %>% dplyr::select(-taxaLevelAbbrev) %>%
                column_to_rownames("TaxaName")
              
              my_taxa_colData_NbClust_Left = my_taxa_colData_NbClust[, c("Cluster", "taxaLevel")]
              my_taxa_colData_NbClust_Left$Phylum = as.character(ifelse(my_taxa_colData_NbClust_Left$taxaLevel=="Phylum", TRUE, FALSE))
              my_taxa_colData_NbClust_Left$Class = as.character(ifelse(my_taxa_colData_NbClust_Left$taxaLevel=="Class", TRUE, FALSE))
              my_taxa_colData_NbClust_Left$Order = as.character(ifelse(my_taxa_colData_NbClust_Left$taxaLevel=="Order", TRUE, FALSE))
              my_taxa_colData_NbClust_Left$Family = as.character(ifelse(my_taxa_colData_NbClust_Left$taxaLevel=="Family", TRUE, FALSE))
              my_taxa_colData_NbClust_Left$Genus = as.character(ifelse(my_taxa_colData_NbClust_Left$taxaLevel=="Genus", TRUE, FALSE))
              my_taxa_colData_NbClust_Left$Species = as.character(ifelse(my_taxa_colData_NbClust_Left$taxaLevel=="Species", TRUE, FALSE))
              my_taxa_colData_NbClust_Right = my_taxa_colData_NbClust_Left
              
              colors50 = c(
                "#C1F242", "#68C5E5", "#A3E86B", "#E2DE48", "#CFE5AE",
                "#EBF1D8", "#55B99F", "#C8EEEC", "#B7A7DF", "#E1CBA3",
                "#A9CBB4", "#E597DE", "#929776", "#C3E67B", "#B08BEA",
                "#82F23A", "#E9B6D6", "#A38B97", "#D12FEE", "#54ECB1",
                "#DCD263", "#DB7341", "#E37EA0", "#B053E2", "#EAF425",
                "#9CEEC8", "#9453A6", "#DC4B9C", "#E0506A", "#ACA550",
                "#69A5DE", "#63EB73", "#7CA3A9", "#E7CCCA", "#E33ED1",
                "#5C74EC", "#C6CEE6", "#A4E9A1", "#E0948A", "#7E3C62",
                "#E5AD35", "#66F0E3", "#6CB061", "#8DBB28", "#DD72E1",
                "#613CD6", "#6382D5", "#7ADCE0", "#E5EA96", "#E2AE74"
              ) 
              ## Max 15 Cluster Colors
              curClusterColors = colors50[1:minClust]
              taxaLevel_colors = c("#FFD919", "#0B701E", "#FF7300", "#120C94", "#FF748C", "#009BB3")
              names(taxaLevel_colors) = c("Phylum", "Class", "Order", "Family", "Genus", "Species")
              
              taxaLevel_p = c("#FFEC9D","#000000") # Yellow
              names(taxaLevel_p) = c("TRUE", "FALSE")  
              
              taxaLevel_c = c("#93C797", "#000000") # Green
              names(taxaLevel_c) = c("TRUE", "FALSE")
              
              taxaLevel_o = c("#FFC19F", "#000000") # Orange 
              names(taxaLevel_o) = c("TRUE", "FALSE")
              
              taxaLevel_f = c("#8194BA", "#000000") # Navy
              names(taxaLevel_f) = c("TRUE", "FALSE")
              
              taxaLevel_g = c("#FFC1DA", "#000000") # Pink
              names(taxaLevel_g) = c("TRUE", "FALSE")
              
              taxaLevel_s = c("#89F3F5", "#000000") # Cerulean
              names(taxaLevel_s) = c("TRUE", "FALSE")
              
              clusterColors = curClusterColors
              names(clusterColors) = sort(unique(my_taxa_colData_NbClust_Left$Cluster))
              colors = list(
                # Altered = altered_colors, 
                taxaLevel = taxaLevel_colors, 
                Phylum = taxaLevel_p, 
                Class = taxaLevel_c, 
                Order = taxaLevel_o, 
                Family = taxaLevel_f, 
                Genus = taxaLevel_g,
                Species = taxaLevel_s,
                Cluster = clusterColors
              )
              
              groups <- as.data.frame(curTree)
              order = rownames(consensusMatrix_Final_wide[curHclust[["order"]],])
              groupFreqTable = as.data.frame(table(groups[order,]))
              curClusterTable = data.frame(taxaName = order, clusterID = groups[order,])
              clusterFreq = groupFreqTable[unique(groups[order,]), ]$Freq
              
              ## Reorder 
              my_taxa_colData_NbClust_Right = my_taxa_colData_NbClust_Left = my_taxa_colData_NbClust_Left[order,]
              
              rowIndices = c()
              for(i in seq_along(clusterFreq)){
                rowIndices = c(rowIndices, sum(clusterFreq[1:i]))
              }
              rowIndices = rowIndices[1:(length(rowIndices)-1)]
              colorHex = colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(51)
              colorScale = seq(-1, 1, length.out = 51)
              reorderdConsensus = as.matrix(consensusMatrix_Final_wide)[order, order]
              diag(reorderdConsensus) = NA
              incProgress(1/n, detail = paste("Part", 2, "Generate Consensus Plot"))
              consensusPlot = pheatmap::pheatmap(mat = reorderdConsensus, 
                                                 main = paste0("# Optimal Clusters: ", minClust, " | ",
                                                               "Phenotype: ", C3NAObj$misc$phenotype, "; "),
                                                 annotation_row = my_taxa_colData_NbClust_Left, 
                                                 annotation_colors = colors, 
                                                 annotation_col = my_taxa_colData_NbClust_Right, 
                                                 fontsize_col = 4, fontsize_row = 4,
                                                 gaps_row = rowIndices,gaps_col = rowIndices,
                                                 cluster_rows = FALSE, cluster_cols = FALSE)
              plots[["consensusPlot"]] = consensusPlot
              consensusPlot
              incProgress(1/n, detail = paste("Part", 3, "Saving..."))
            }
          })
        })
        
        ##### 2.1.3 Correlation Matrix ######
        output$correlationPlot <- renderPlot({
          withProgress(message = 'Generating Correlation Plot, please wait...', value = 0, {
            n=3
            incProgress(1/n, detail = paste("Part", 1, "calculate Correlation Data"))
            if(!is.null(values[["selectedMinClust"]])){
              minClust = values[["selectedMinClust"]]
              curHclust = values[["curHclust"]]
              consensusMatrix_Final_wide = values[["consensusMatrix_Final_wide"]]
              ## Plot the same plot with correlation
              SparccP2 = C3NAObj$sparCCTable %>% 
                mutate(Var3 = Var1) %>% 
                select(-Var1) %>%
                rename(Var1 = Var2) %>%
                rename(Var2 = Var3) %>%
                filter(Var1 != Var2)
              
              reordered_sparccCor = C3NAObj$sparCCTable %>%
                rbind(SparccP2) %>%
                select(-p, -fdr) %>%
                pivot_wider(names_from = Var2, values_from = cor) %>%
                column_to_rownames("Var1")
              # Obtain the clustering information for the taxa
              curTree = cutree(as.hclust(curHclust), k = minClust)
              
              ## Generate Pheatmap
              my_taxa_colData_NbClust <- as.data.frame(curTree)
              colnames(my_taxa_colData_NbClust)[1] = "Cluster"
              my_taxa_colData_NbClust$Cluster = as.character(paste0("Cluster ", 
                                                                    formatC(my_taxa_colData_NbClust$Cluster,flag=0,width=2)))
              my_taxa_colData_NbClust$TaxaName = rownames(my_taxa_colData_NbClust)
              my_taxa_colData_NbClust = my_taxa_colData_NbClust %>% as_tibble() %>%
                mutate(taxaLevelAbbrev = sapply(strsplit(TaxaName, "_"), "[", 1)) %>%
                mutate(taxaLevel = case_when(
                  taxaLevelAbbrev =="p" ~ "Phylum",
                  taxaLevelAbbrev =="c" ~ "Class",
                  taxaLevelAbbrev =="o" ~ "Order",
                  taxaLevelAbbrev =="f" ~ "Family",
                  taxaLevelAbbrev =="g" ~ "Genus", 
                  taxaLevelAbbrev =="s" ~ "Species",
                )) %>% dplyr::select(-taxaLevelAbbrev) %>%
                column_to_rownames("TaxaName")
              
              my_taxa_colData_NbClust_Left = my_taxa_colData_NbClust[, c("Cluster", "taxaLevel")]
              my_taxa_colData_NbClust_Left$Phylum = as.character(ifelse(my_taxa_colData_NbClust_Left$taxaLevel=="Phylum", TRUE, FALSE))
              my_taxa_colData_NbClust_Left$Class = as.character(ifelse(my_taxa_colData_NbClust_Left$taxaLevel=="Class", TRUE, FALSE))
              my_taxa_colData_NbClust_Left$Order = as.character(ifelse(my_taxa_colData_NbClust_Left$taxaLevel=="Order", TRUE, FALSE))
              my_taxa_colData_NbClust_Left$Family = as.character(ifelse(my_taxa_colData_NbClust_Left$taxaLevel=="Family", TRUE, FALSE))
              my_taxa_colData_NbClust_Left$Genus = as.character(ifelse(my_taxa_colData_NbClust_Left$taxaLevel=="Genus", TRUE, FALSE))
              my_taxa_colData_NbClust_Left$Species = as.character(ifelse(my_taxa_colData_NbClust_Left$taxaLevel=="Species", TRUE, FALSE))
              my_taxa_colData_NbClust_Right = my_taxa_colData_NbClust_Left
              
              colors50 = c(
                "#C1F242", "#68C5E5", "#A3E86B", "#E2DE48", "#CFE5AE",
                "#EBF1D8", "#55B99F", "#C8EEEC", "#B7A7DF", "#E1CBA3",
                "#A9CBB4", "#E597DE", "#929776", "#C3E67B", "#B08BEA",
                "#82F23A", "#E9B6D6", "#A38B97", "#D12FEE", "#54ECB1",
                "#DCD263", "#DB7341", "#E37EA0", "#B053E2", "#EAF425",
                "#9CEEC8", "#9453A6", "#DC4B9C", "#E0506A", "#ACA550",
                "#69A5DE", "#63EB73", "#7CA3A9", "#E7CCCA", "#E33ED1",
                "#5C74EC", "#C6CEE6", "#A4E9A1", "#E0948A", "#7E3C62",
                "#E5AD35", "#66F0E3", "#6CB061", "#8DBB28", "#DD72E1",
                "#613CD6", "#6382D5", "#7ADCE0", "#E5EA96", "#E2AE74"
              ) 
              ## Max 15 Cluster Colors
              curClusterColors = colors50[1:minClust]
              taxaLevel_colors = c("#FFD919", "#0B701E", "#FF7300", "#120C94", "#FF748C", "#009BB3")
              names(taxaLevel_colors) = c("Phylum", "Class", "Order", "Family", "Genus", "Species")
              
              taxaLevel_p = c("#FFEC9D","#000000") # Yellow
              names(taxaLevel_p) = c("TRUE", "FALSE")  
              
              taxaLevel_c = c("#93C797", "#000000") # Green
              names(taxaLevel_c) = c("TRUE", "FALSE")
              
              taxaLevel_o = c("#FFC19F", "#000000") # Orange 
              names(taxaLevel_o) = c("TRUE", "FALSE")
              
              taxaLevel_f = c("#8194BA", "#000000") # Navy
              names(taxaLevel_f) = c("TRUE", "FALSE")
              
              taxaLevel_g = c("#FFC1DA", "#000000") # Pink
              names(taxaLevel_g) = c("TRUE", "FALSE")
              
              taxaLevel_s = c("#89F3F5", "#000000") # Cerulean
              names(taxaLevel_s) = c("TRUE", "FALSE")
              
              clusterColors = curClusterColors
              names(clusterColors) = sort(unique(my_taxa_colData_NbClust_Left$Cluster))
              colors = list(
                # Altered = altered_colors, 
                taxaLevel = taxaLevel_colors, 
                Phylum = taxaLevel_p, 
                Class = taxaLevel_c, 
                Order = taxaLevel_o, 
                Family = taxaLevel_f, 
                Genus = taxaLevel_g,
                Species = taxaLevel_s,
                Cluster = clusterColors
              )
              
              groups <- as.data.frame(curTree)
              order = rownames(consensusMatrix_Final_wide[curHclust[["order"]],])
              groupFreqTable = as.data.frame(table(groups[order,]))
              curClusterTable = data.frame(taxaName = order, clusterID = groups[order,])
              clusterFreq = groupFreqTable[unique(groups[order,]), ]$Freq
              
              ## Reorder 
              my_taxa_colData_NbClust_Right = my_taxa_colData_NbClust_Left = my_taxa_colData_NbClust_Left[order,]
              
              rowIndices = c()
              for(i in seq_along(clusterFreq)){
                rowIndices = c(rowIndices, sum(clusterFreq[1:i]))
              }
              rowIndices = rowIndices[1:(length(rowIndices)-1)]
              colorHex = colorRampPalette(rev(brewer.pal(n=11, name="RdBu")))(51)
              colorScale = seq(-1, 1, length.out = 51)
              reorderdConsensus = as.matrix(reordered_sparccCor)[order, order]
              diag(reorderdConsensus) = NA
              incProgress(1/n, detail = paste("Part", 2, "Generate Correlation Plot"))
              correlationPlot = pheatmap::pheatmap(mat = reorderdConsensus, 
                                                   main = paste0("# Optimal Clusters: ", minClust, " | ",
                                                                 "Phenotype: ", C3NAObj$misc$phenotype, "; "),
                                                   annotation_row = my_taxa_colData_NbClust_Left, 
                                                   annotation_colors = colors, 
                                                   color = colorHex, breaks = colorScale,
                                                   annotation_col = my_taxa_colData_NbClust_Right, 
                                                   fontsize_col = 4, fontsize_row = 4,
                                                   gaps_row = rowIndices,gaps_col = rowIndices,
                                                   cluster_rows = FALSE, cluster_cols = FALSE)
              plots[["correlationPlot"]] = correlationPlot
              correlationPlot
              incProgress(1/n, detail = paste("Part", 3, "Saving..."))
            }
          })
        })
      })
      
      
      ##### Download Panel ###### 
      output$plotSelector <- renderUI({
        if(input$plotRadioButtons %in% c("PNG", "PDF")){
          ## HTML Plots
          avaliablePlots <- sort(names(plots))
          listOfImportantData <- data.frame(
            original = c("modulePatternPlot", "silPlot", 
                         "consensusPlot", "correlationPlot"),
            full = c("Module Pattern Plot", "Silhouette Evaluation Plot",
                     "Consensus Plot", "Correlation Plot"))
          checker <- avaliablePlots %in% listOfImportantData$original
          tempData <- avaliablePlots[checker]
          listData <- listOfImportantData %>%
            filter(original %in% tempData)
          
          ## Obtain the list of plots available
          if(nrow(listData) < 1){
            h4("No plot was generated in Panel 1. Please check your data. ")
          } else {
            radioGroupButtons(
              inputId = "plotSelector",
              choiceNames = listData$full, 
              choiceValues = listData$original,
              individual = TRUE,
              checkIcon = list(
                yes = tags$i(class = "fa fa-check-square", 
                             style = "color: white"),
                no = tags$i(class = "fa fa-square-o", 
                            style = "color: black")),
              status = "custom-status-downloadPlots"
            )
          }
        }
      })
      
      ## Observe the download button and handle it
      selectedPlotInput <- function(){
        req(input$plotSelector)
        plots[[input$plotSelector]]
      }
      
      # create filename
      selectedPlotName <- reactive({
        name <- input$plotSelector
        if(input$plotRadioButtons == "PDF")   filename <- paste0(name, "-", Sys.Date(), ".pdf",  sep="")
        if(input$plotRadioButtons == "PNG")   filename <- paste0(name, "-", Sys.Date(), ".png",  sep="")
        return(filename)
      })
      
      output$downloadPlots <- downloadHandler(
        filename = selectedPlotName,
        content = function(file) {
          name <- input$plotSelector
          
          if(input$plotRadioButtons %in% c("PDF")) {
            ggsave(file, plot = selectedPlotInput(), 
                   height = as.numeric(input$heightInput), 
                   width = as.numeric(input$widthInput), units = "in")
          } else if(input$plotRadioButtons == "PNG") {
            ggsave(file, plot = selectedPlotInput(), dpi = 300, 
                   height = as.numeric(input$heightInput), 
                   width = as.numeric(input$widthInput), units = "in")
          }
        }
      )
      
      ##### 3.0 All Showmodels ######
      observeEvent(input$MS_PatternVis, {
        showModal(
          modalDialog(
            title = HTML('<b>Visualization of Module Patterns across Different Minimal Module Sizes</b>'),
            HTML("
The 'Pattern Visualization' section hosts the patterns in turns of the unique module due to the different minimal numbers of taxa per module. <hr>
The <b>x-axis</b> represents the different range of a minimal number of taxa per module; by default, it investigated between 3 and 40, and this range should be sufficient for observing all the necessary patterns. 
<br>
The <b>y-axis</b> represents The relative abundance of taxa per module.The colors are assigned using 'randomcolorR' package with a seed, which can be adjusted by setting the 'colorSeed' arguments when calling this Shiny application. 
<br>
The number of unique modules is displayed on top of the bar plot. 
<hr>
The plot can be saved from 'Panel 2' on the top bar. 
              "),
easyClose = TRUE
          )
        )
      })
      
      observeEvent(input$MS_UniquePatternSelector, {
        showModal(
          modalDialog(
            title = HTML('<b>Selection of Unique Patterns that will be included in Consensus Clustering</b>'),
            HTML("
The unique pattern selector enables you to manually adjust which patterns will be included to generate the consensus matrix. 
The optimal last pattern is generally around when the number of unique patterns drops to 10, with more <p style='background-color:gray;display:inline;color:white;'>prevalent duplicated patterns (darkened bar)</p> compared to patterns with more unique modules. 
Please check the supplement results from our publication for a more detailed investigation of the impact of different selected patterns.               "),
easyClose = TRUE
          )
        )
      })
      
      observeEvent(input$MS_OptimalModuleSelector, {
        showModal(
          modalDialog(
            title = HTML('<b>Manually Confirm the Number of Modules based on Consensus Clustering Results</b>'),
            HTML("
The Optimal Module Selector should be determined by the Silhouette, Consensus, and Correlation Plots below, with the Silhouette plot being the main factor. <hr>
The <b><font color='#2736BA'>Silhouette curve</font></b>  (based on the consensus matrix, in blue) should follow the pattern of quick incrementation with the increased number of clusters followed by a plateau region where the incrementation rate is much slower. 
The optimal number of clusters should be around the turning point between these two clusters with a lower than 10% proportion of intra-modular zeros, which is shown in the <b><font color='#C23237'>intra-modular proportion of zero curve</font></b>. 
<hr>
Please check the supplement results from our publication for a more detailed investigation of this topic. 

              "),
easyClose = TRUE
          )
        )
      })
      
      observeEvent(input$MS_OptimalModuleCodeGenerator, {
        showModal(
          modalDialog(
            title = HTML('<b>Manually Confirm the Number of Modules based on Consensus Clustering Results</b>'),
            HTML(
              "The code generator displays the code for the user to run after confirming the optimal pattern selection and the optimal number of clusters of taxa based on the consensus matrix. The user should select the display text and copy&paste it to R to run after closing the Shiny application. Please replace the 'oldC3NAObj' and 'newC3NAObj' with the correct and preferred name. <hr><b>Note: </b> The code is also displayed on the console. 
"
            ),
easyClose = TRUE
          )
        )
      })
      
      observeEvent(input$MS_ConsensusEval, {
        showModal(
          modalDialog(
            title = HTML('<b>Optimal Number of Module Selection Based on the Consensus Results</b>'),
            HTML(
              "The consensus eval section hosts the three module evaluation tools, the Silhouette Plot, Consensus Plot, and Correlation Plot. These should be used to determine the optimal number of clusters based on the consensus matrix generated using the selected patterns. <hr>
The <b><font color='#2736BA'>Silhouette curve</font></b>  (based on the consensus matrix, in blue) should follow the pattern of quick incrementation with the increased number of clusters followed by a plateau region where the incrementation rate is much slower. 
The optimal number of clusters should be around the turning point between these two clusters with a lower than 10% proportion of intra-modular zeros, which is shown in the <b><font color='#C23237'>intra-modular proportion of zero curve</font></b>. <br>
The <b>Consensus Plot</b> includes the consensus heatmap with the taxa name, taxonomic levels, and clusters shown. <br>
The <b>Correlation Plot</b> includes the correlation heatmap with the taxa name, taxonomic levels, and clusters shown. <br>
<hr>
To update the number of clusters for the consensus and correlation plot, please use the 'Confirm Selected Patterns' above. The legend for the clusters might not be displayed fully within the shiny display due to limited dimensions; please visit the download panel to select a wider width for the correct display of the cluster legend. The plot can be saved from 'Panel 2' on the top bar. "

            ),
easyClose = TRUE
          )
        )
      })
      
      observeEvent(input$MS_DownloadPanel, {
        showModal(
          modalDialog(
            title = HTML('<b>Download Displayed Figures in Panel 1</b>'),
            HTML(
"The download panel enables you to download the currently displayed plots from Panel 1 in PDF or PNG format. For the PNG format, the DPI is set to 300. The file should be saved to your default download folder with the plot type name followed by the date and time to avoid overwriting any plot. "            
            ),
easyClose = TRUE
          )
        )
      })
      
      # stop App on closing the browser
      session$onSessionEnded(function() {
        stopApp()
      })
      
    },
    # launch App in a browser
    options = list(launch.browser = TRUE)
  )
}


#' Shiny Application for C3NAobj Comparison 
#' 
#' @description This function calls the interactive two phenotypes after running the 
#' comparePhenotypes() step. 
#' 
#' @param C3NAObj (Required) The C3NAObj with Compared Phenotypes
#' 
#' @rawNamespace import(plotly, except = last_plot)
#' @import shinyWidgets
#' @import shinydashboard
#' @import visNetwork 
#' @import reactable
#' @importFrom DT renderDataTable dataTableOutput datatable
#' @importFrom base64enc dataURI
#' @importFrom htmlwidgets saveWidget
#' @export
#' 
compareTwoPhenoShiny = function(C3NAObj){
  ## Quick Check
  if(is.null(C3NAObj$misc$corCut)){
    stop("Please run the comparePhenotypes() before running the interactive Shiny application. ")
  }
  ## Extract the needed data for Shiny
  sparccP_Filtered_rbind = C3NAObj$sparccP_Filtered_rbind
  sparccP_Filtered_Combined = C3NAObj$sparccP_Filtered_Combined
  modulePlotData_Wide = C3NAObj$modulePreservation$modulePlotData_Wide
  nodesAll = C3NAObj$nodes$nodesAll
  ## Filter nodes All 
  uniqueTaxa = unique(unlist(sparccP_Filtered_Combined[, c("target", "source")]))
  nodesAll = subset(nodesAll, TaxaName %in% uniqueTaxa)
  nodesAll_Pheno = merge(nodesAll, C3NAObj$nodes$refTaxaTable, 
                   by.x = "TaxaName", by.y = "taxaName", all.x = TRUE)
  nodesAll_Pheno = nodesAll_Pheno %>%
    mutate(taxaLevelAbbrev = sapply(strsplit(TaxaName, "_"), "[", 1)) %>%
    mutate(`Taxonomic Levels` = case_when(
      taxaLevelAbbrev =="p" ~ "Phylum",
      taxaLevelAbbrev =="c" ~ "Class",
      taxaLevelAbbrev =="o" ~ "Order",
      taxaLevelAbbrev =="f" ~ "Family",
      taxaLevelAbbrev =="g" ~ "Genus", 
      taxaLevelAbbrev =="s" ~ "Species"
    ),
    `Number of Unique Taxa` = 1)
  nodesAll_PhenoOri = nodesAll_Pheno
  statsTable = C3NAObj$statsTable %>%
    filter(`Full Taxa Name` %in% nodesAll$TaxaName)
  nodesAbbre = C3NAObj$nodes$nodesAbbre
  samePhyloTable = C3NAObj$nodes$samePhyloTable
  C3NAResults = C3NAObj$C3NA_Wilcoxon
  data1Name = C3NAObj$miscDisease$phenotype
  data2Name = C3NAObj$miscControl$phenotype
  curDisease_OriCountTable = C3NAObj$curDisease_OriCountTable
  curControl_OriCountTable = C3NAObj$curControl_OriCountTable
  corCut = C3NAObj$misc$corCut
  fdr = C3NAObj$misc$fdr
  removedTaxaTable = C3NAObj$removedTaxaTable
  consensusProp_Disease = C3NAObj$consensusDisease$consensusProp
  consensusProp_Disease = consensusProp_Disease %>%
    select(ConsensusProp, taxaComboName) %>%
    rename(DiseaseConsensusProp = ConsensusProp) %>%
    mutate(DiseaseConsensusProp = round(DiseaseConsensusProp, digits = 4))
  consensusProp_Control = C3NAObj$consensusControl$consensusProp 
  consensusProp_Control = consensusProp_Control %>%
    select(ConsensusProp, taxaComboName) %>%
    rename(ControlConsensusProp = ConsensusProp) %>%
    mutate(ControlConsensusProp = round(ControlConsensusProp, digits = 4))
  phenotypes = c(data1Name, data2Name)
  ## DA Methods if present
  daData = C3NAObj$DA
  daMethods = names(sapply(daData, class))[as.vector(as.vector(sapply(daData, class))) == "logical"]
  daMethods = c(daMethods, "None")
  
  shinyApp(
    ui = navbarPage(
      title = "C3NA",
      tags$head(
        tags$style(
          HTML(
            ".checkbox-inline { 
                    margin-left: 0px;
                    margin-right: 10px;
              }
             .checkbox-inline+.checkbox-inline {
                        margin-left: 0px;
                        margin-right: 10px;
             }
            .multi-wrapper {height: 415px !important;}
            .multi-wrapper .non-selected-wrapper, .multi-wrapper .selected-wrapper {height: 370px !important;}
            img {float: left;}
            .myHeader {
              font-weight: 900; font-size: 20px; margin: 0 -600rem; padding: 0.15rem 600rem;
              background: rgba(250, 205, 15, 0.2); /*Original: 666*/
              border-bottom: 5rem; -webkit-transition: all .5s ease-in-out;
              -moz-transition: all .5s ease-in-out; -o-transition: all .5s ease-in-out;
              transition: all .5s ease-in-out; text-decoration: underline solid transparent;
              box-shadow: 0 0 11px;
            }
            .btn-custom-status-downloadPlots:hover, .btn-custom-status-downloadPlots:active, 
            .btn-custom-status-downloadPlots:focus, .btn-custom-status-downloadPlots:visited, 
            .btn-custom-status-downloadPlots.active, .open .dropdown-toggle.btn-custom-status-downloadPlots { 
              color: #ffffff;
              background-color: #318dc3;
              border-color: #07486f;
              box-shadow: inset 4px 4px 4px 0 #104667,
                          inset -4px -4px 4px 0 #5ebdf7;
            }
            .btn-custom-status-downloadPlots { 
              color: black; 
              background-color: #88bcdc ; 
              border-color: #07486f;
              box-shadow: 4px 4px 4px 0 rgba(0, 0, 0, 0.2),
                          -4px -4px 4px 0 rgba(255, 255, 255, 0.5),
                          inset 3px 3px 3px 0 #bdd5e4,
                          inset -3px -3px 3px 0 #0c3a56;
            } 
            "
          )
        ) 
      ),
      tabPanel(
        ## First Panel
        title = HTML("<b>Panel 1: Compare Two Conditions</b>"),
        value = "panel1",
        useShinyjs(),
        useShinydashboard(),
        setShadow(class = "box"),
        
        tags$div(class = "Banner", 
                 tags$div(class = "myHeader", id = "infoPanel1Sec1",
                          (HTML("Section 1: Dataset Information")
                          )
                 )
        ),
        br(),
        ##### Section 1 General Taxa #####
        fluidRow(
          box(
            width = 3,
            title = HTML("Comparison Phenotype", "<font size='3'>", 
                         as.character(actionLink(inputId = "MS_Pheno1", 
                                                 label = "", 
                                                 icon = icon("question-circle"))), "</font>"), 
            
            status = "warning",
            solidHeader = TRUE, # height = "33em", 
            plotlyOutput("pheno1Sunburst")
          ),
          box(
            width = 3,
            title = HTML("Reference Phenotype", "<font size='3'>", 
                         as.character(actionLink(inputId = "MS_Pheno2", 
                                                 label = "", 
                                                 icon = icon("question-circle"))), "</font>"), 
            
            status = "warning",
            solidHeader = TRUE, # height = "33em", 
            plotlyOutput("pheno2Sunburst")
          ),
          box(
            width = 6,
            title = HTML("Phenotype-based Taxa Comparison", "<font size='3'>", 
                         as.character(actionLink(inputId = "MS_PhenoComparison", 
                                                 label = "", 
                                                 icon = icon("question-circle"))), "</font>"), 
            
            status = "warning",
            solidHeader = TRUE, # height = "33em", 
            reactable::reactableOutput("taxaRelatedStatTable")
          )
        ),
        
        tags$div(class = "Banner", 
                 tags$div(class = "myHeader", id = "infoPanel1Sec1",
                          (HTML("Section 2: Preservation Statistics")
                          )
                 )
        ),
        br(),
        ##### Section 2 Taxa Statistics #####
        fluidRow(
          box(
            width = 4,
            title = HTML("Taxa Comparison among Modules from Two Phenotypes", "<font size='3'>", 
                         as.character(actionLink(inputId = "MS_ModModComparison", 
                                                 label = "", 
                                                 icon = icon("question-circle"))), "</font>"), 

            status = "info",
            solidHeader = TRUE, height = "33em", 
            plotlyOutput("moduleModuleComparison")
          ),
          box(
            width = 4,
            title = HTML("Module Preservation Statistics", "<font size='3'>", 
                         as.character(actionLink(inputId = "MS_ModPreservation", 
                                                 label = "", 
                                                 icon = icon("question-circle"))), "</font>"), 
            
            status = "info",
            solidHeader = TRUE,  height = "33em", 
            plotlyOutput("modulePreservation")
          ),
          box(
            width = 4,
            title = HTML("Module Membership", "<font size='3'>", 
                         as.character(actionLink(inputId = "MS_ModMembership", 
                                                 label = "", 
                                                 icon = icon("question-circle"))), "</font>"), 
            status = "info",
            solidHeader = TRUE, height = "33em", 
            DT::dataTableOutput("clickedModuleInfo")
          )
        ),
        #####  Section 3 Taxa Specific investigation #####
        tags$div(class = "Banner", 
                 tags$div(class = "myHeader", id = "infoPanel1Sec1",
                          (HTML("Section 3: Network Plots")
                          )
                 )
        ),
        br(),
        fluidRow(
          box(
            width = 4,
            title = HTML("Taxa Selector", "<font size='3'>", 
                         as.character(actionLink(inputId = "MS_TaxaSelector", 
                                                 label = "", 
                                                 icon = icon("question-circle"))), "</font>"), 
            status = "danger",
            solidHeader = TRUE,  height = "36em", 
            uiOutput("taxaSelector")
          ), 
          box(
            width = 5,
            title = HTML("Phenotype-Specific Taxa Membership", "<font size='3'>", 
                         as.character(actionLink(inputId = "MS_ModBarPlot", 
                                                 label = "", 
                                                 icon = icon("question-circle"))), "</font>"), 
            status = "danger",
            solidHeader = TRUE, height = "36em", 
            ##### . 4.1 VisNetwork Drop Down #####
            dropdown(
              radioGroupButtons(
                inputId = "barPlotPhenoSelector",
                label = "Phenotype Selector:", 
                choices = phenotypes,
                selected = phenotypes[1],
                status = "custom-status"
              ),
              uiOutput("barPlotTaxaSelector"),
              circle = TRUE, status = "danger", 
              icon = icon("filter"), width = "600px",
              animate = animateOptions(
                enter = shinyWidgets::animations$bouncing_entrances$bounceInLeft,
                exit = shinyWidgets::animations$bouncing_exits$bounceOutRight
              ),
              tooltip = tooltipOptions(title = "Click to Switch Phenotypes and Selected Taxa")
            ),
            plotlyOutput("barPlot")
          ),
          box(
            width = 3,
            title = HTML("Selected Taxa Modular Information", "<font size='3'>", 
                         as.character(actionLink(inputId = "MS_TaxaModTable", 
                                                 label = "", 
                                                 icon = icon("question-circle"))), "</font>"), 
            status = "danger",
            solidHeader = TRUE, height = "36em", 
            
            dropdown(
              uiOutput("taxaModularTaxaSelector"),
              circle = TRUE, status = "danger", 
              icon = icon("filter"), width = "600px",
              animate = animateOptions(
                enter = shinyWidgets::animations$bouncing_entrances$bounceInRight,
                exit = shinyWidgets::animations$bouncing_exits$bounceOutRight
              ),
              tooltip = tooltipOptions(title = "Click to Switch Taxon")
            ),
            hr(),
            DT::dataTableOutput("intramodularCorTable")
          )
        ),
        ##### Section 3 Network Plot #####
        conditionalPanel(
          condition = "input.refreshPanel1Section3==0",
          fluidRow(
            box(
              width = 4,
              title = HTML("Proceed to Network", "<font size='3'>", 
                           as.character(actionLink(inputId = "MS_ClickToProceeed", 
                                                   label = "", 
                                                   icon = icon("question-circle"))), "</font>"), 
              status = "success",
              solidHeader = TRUE,
              actionBttn(
                inputId = "refreshPanel1Section3",
                label = "Plot Network based on Selection",
                style = "material-flat",
                color = "primary"
              )
            )
          )
        ),
        conditionalPanel(
          condition = "input.refreshPanel1Section3>0",
          fluidRow(
            box(
              width = 8,
              title = HTML("Network Plot", "<font size='3'>", 
                           as.character(actionLink(inputId = "MS_Network", 
                                                   label = "", 
                                                   icon = icon("question-circle"))), "</font>"), 
              status = "success",
              height = "63em",
              solidHeader = TRUE,
              ##### . 4.1 VisNetwork Drop Down #####
              dropdown(
                prettyRadioButtons(
                  inputId = "nodesMode",
                  label = "Node-Node Connections Mode:",
                  selected = "selectedTaxa",
                  choiceNames =  c("Only the selected taxa",
                                   "All taxa directly related to selected taxa",
                                   "All taxa in related modules"),
                  choiceValues = c("selectedTaxa",
                                   "relatedTaxa",
                                   "allRelatedTaxa"),
                  status = "primary",
                  fill = TRUE
                ),
                radioGroupButtons(
                  inputId = "samePhylo",
                  label = "Remove Connection from Same Phylogenenetic Branch?",
                  choices  = c("Remove", "Keep"),
                  selected = c("Keep"),
                  checkIcon = list(
                    yes = tags$i(class = "fa fa-check-square",
                                 style = "color: steelblue"),
                    no = tags$i(class = "fa fa-square-o",
                                style = "color: steelblue"))
                ),
                checkboxGroupButtons(
                  inputId = "corMode",
                  label = "Edges Mode:",
                  choiceNames = c("Disease-Only",
                                  "Control-Only",
                                  "Shared"),
                  choiceValues = c("diseaseOnly",
                                   "controlOnly",
                                   "both"),
                  selected = c("diseaseOnly",
                               "controlOnly",
                               "both"),
                  checkIcon = list(
                    yes = tags$i(class = "fa fa-check-square",
                                 style = "color: steelblue"),
                    no = tags$i(class = "fa fa-square-o",
                                style = "color: steelblue"))
                ),
                checkboxGroupButtons(
                  inputId = "taxaLvls",
                  label = "Taxonomic Levels:",
                  choiceNames = c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
                  choiceValues = c("p", "c", "o", "f", "g", "s"),
                  selected = c("p", "c", "o", "f", "g", "s"),
                  checkIcon = list(
                    yes = tags$i(class = "fas fa-virus",
                                 style = "color: steelblue"),
                    no = tags$i(class = "fas fa-virus-slash",
                                style = "color: steelblue"))
                ),
                hr(),
                prettyRadioButtons(
                  inputId = "DASelector",
                  label = "Differential Abundance: ", inline = TRUE,
                  selected = "None", choices = daMethods,
                  status = "primary", fill = TRUE
                ),
                sliderInput(
                  inputId = "corSelector",
                  label = "Select a Correlation Cut-off Values",
                  min = C3NAObj$misc$corCut, max = 1, step = 0.05,
                  value = 0.2
                ),
                sliderTextInput(
                  inputId = "nodeBorderWidth",
                  label = "Node Border Color Width: ",  choices = 1:15,
                  selected = 6, grid = TRUE
                ),
                sliderTextInput(
                  inputId = "nodeLabelSize",
                  label = "Node Label Size: ", choices = 1:15,
                  selected = 6, grid = TRUE
                ),
                sliderTextInput(
                  inputId = "edgeLabelSize",
                  label = "Edge Label Size: ", choices = 1:15,
                  selected = 6, grid = TRUE
                ),
                circle = TRUE, status = "success",
                icon = icon("filter"), width = "600px",
                animate = animateOptions(
                  enter = shinyWidgets::animations$bouncing_entrances$bounceInLeft,
                  exit = shinyWidgets::animations$bouncing_exits$bounceOutLeft
                ),
                tooltip = tooltipOptions(title = "Click to change Network Plot Settings")
              ),
              visNetworkOutput("SelectedTaxaPlots",
                               height = "800px")
            ),
            column(
              width = 4,
              box(
                width = 12,
                title = HTML("Selected Taxon Count Distribution", "<font size='3'>", 
                             as.character(actionLink(inputId = "MS_SingelTaxaViolin", 
                                                     label = "", 
                                                     icon = icon("question-circle"))), "</font>"), 
                status = "success",
                height = "31em",
                solidHeader = TRUE,
                plotly::plotlyOutput("countDistributionPlot", height = "380px")

              ),
              box(
                width = 12,
                title = HTML("Sequentially Selected Taxa - Comparison", "<font size='3'>", 
                             as.character(actionLink(inputId = "MS_TwoTaxaViolin", 
                                                     label = "", 
                                                     icon = icon("question-circle"))), "</font>"),
                status = "success",
                height = "31em",
                solidHeader = TRUE,
                plotly::plotlyOutput("countDistributionPlot_Compare", height = "380px")
              )
            )
          ),
          ###### Section 5.1 Empty Space at the Bottom ######
          fluidRow(
            hr()
          )
        )
      ),
      tabPanel(
        ## First Panel
        title = HTML("<b>Panel 2: Download Panel </b>"),
        value = "panel2",
        fluidRow(
          box( ### 5.1 Download Plot List and Button ### 
            width = 6, 
            title = HTML("Download Plots", "<font size='3'>", 
                         as.character(actionLink(inputId = "MS_DownloadPanel", 
                                                 label = "", 
                                                 icon = icon("question-circle"))), "</font>"), 
            status = "primary", solidHeader = TRUE, 
            uiOutput("plotSelector"), 
            hr(),
            ## Conditional Panels 
            conditionalPanel(
              condition = "input.plotRadioButtons=='HTML'",
              fluidRow(
                column(
                  width = 6,
                  textInput(inputId = "heightInput",
                            label = "Height in Pixel",
                            placeholder =  "Height in Pixel",
                            value = "1920"),
                ),
                column(
                  width = 6,
                  textInput(inputId = "widthInput",
                            label = "Width in Pixel",
                            placeholder =  "Width in Pixel",
                            value = "1200"),
                )
              ),
              h6(HTML("Please select the heigth and width for the plot, unit: Pixel. ")), # The default dpi is set to 300.
              hr()
            ),
            radioGroupButtons(
              inputId = "plotRadioButtons",
              label = "Select a format for the plot",
              choices = c("HTML"),
              justified = TRUE,
              checkIcon = list(
                yes = icon("ok", 
                           lib = "glyphicon")),
              status = "custom-status-downloadPlots"
            ), hr(),
            fluidRow(
              column(
                12,
                downloadBttn(
                  outputId = "downloadPlots",
                  label = "Download the Selected Plot", 
                  style = "gradient",
                  color = "primary"
                )
              )
            )
          )
        )
      )
    ),
    
    server = function(input, output, session) {
      # Reactive values
      values <- reactiveValues(
        panel1SelectedIDs = NULL,
        click1 = NULL, click1_TaxaLvl = NULL,
        click2 = NULL, click2_TaxaLvl = NULL,
        fullegend = FALSE, selectedFuncs = NULL
      )
      plots_Plotly <- reactiveValues()
      
      # First Panel
      ###### Section 1 General Stat ###### 
      ### Taxa Related StatTable
      output$taxaRelatedStatTable <- reactable::renderReactable({
          reactable(
              statsTable,
              groupBy = c("Group", "Taxonomic Levels"),
              columns = list(
                  `Taxonomic Levels` = colDef(aggregate = "unique"),
                  `Number of Unique Taxa` = colDef(aggregate = "sum", format = colFormat(separators = TRUE))
              )
          )
      })
      
      ##### Sunburst Disease #####
      output$pheno1Sunburst <- renderPlotly({
        nodeData_Disease = subset(nodesAll, !is.na(ClusterID_Disease)) %>%
          mutate(taxaLevelAbbrev = sapply(strsplit(TaxaName, "_"), "[", 1)) %>%
          mutate(`Taxonomic Levels` = case_when(
            taxaLevelAbbrev =="p" ~ "Phylum",
            taxaLevelAbbrev =="c" ~ "Class",
            taxaLevelAbbrev =="o" ~ "Order",
            taxaLevelAbbrev =="f" ~ "Family",
            taxaLevelAbbrev =="g" ~ "Genus", 
            taxaLevelAbbrev =="s" ~ "Species"
          )) %>% 
          group_by(`Taxonomic Levels`) %>%
          dplyr::summarise(n = n())
        nodeData_Disease$`Taxonomic Levels` = factor(nodeData_Disease$`Taxonomic Levels`, 
                                                     levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))
        nodeData_Disease = nodeData_Disease %>% arrange(`Taxonomic Levels`)
        diseaseData = data.frame(
          labels = c(paste0(data1Name, " | nTaxa: ", sum(nodeData_Disease$n)), 
                     as.character(nodeData_Disease$`Taxonomic Levels`)),
          parents = c("", rep(paste0(data1Name, " | nTaxa: ", sum(nodeData_Disease$n)), nrow(nodeData_Disease))),
          values = c(0, nodeData_Disease$n),
          colors = c("#ccaf1d", "#0B701E", "#FF7300", "#120C94", "#db2545", "#009BB3", "#EBEBEB")
        )
        sunburst_Disease = 
          plot_ly(
            diseaseData, 
            labels = ~labels, parents = ~parents, 
            values = ~values,# ids = ~ids, 
            type = "sunburst"
          ) %>%
            layout(colorway  = ~colors,
                   margin = list(l = 0, r = 0, b = 0, t = 0))
        sunburst_Disease$sizingPolicy$padding <- "0"
        plots_Plotly[["sunburst_Disease"]] <- sunburst_Disease
        
        sunburst_Disease
      })
      
      ##### Sunburst Reference #####
      output$pheno2Sunburst <- renderPlotly({
        nodeData_Control = subset(nodesAll, !is.na(ClusterID_Control)) %>%
          mutate(taxaLevelAbbrev = sapply(strsplit(TaxaName, "_"), "[", 1)) %>%
          mutate(`Taxonomic Levels` = case_when(
            taxaLevelAbbrev =="p" ~ "Phylum",
            taxaLevelAbbrev =="c" ~ "Class",
            taxaLevelAbbrev =="o" ~ "Order",
            taxaLevelAbbrev =="f" ~ "Family",
            taxaLevelAbbrev =="g" ~ "Genus", 
            taxaLevelAbbrev =="s" ~ "Species"
          )) %>% 
          group_by(`Taxonomic Levels`) %>%
          dplyr::summarise(n = n())
        nodeData_Control$`Taxonomic Levels` = factor(nodeData_Control$`Taxonomic Levels`, 
                                                     levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species"))
        nodeData_Control = nodeData_Control %>% arrange(`Taxonomic Levels`)
        controlData = data.frame(
          labels = c(paste0(data2Name, " | nTaxa: ", sum(nodeData_Control$n)), 
                     as.character(nodeData_Control$`Taxonomic Levels`)),
          parents = c("", rep(paste0(data2Name, " | nTaxa: ", sum(nodeData_Control$n)), nrow(nodeData_Control))),
          values = c(0, nodeData_Control$n),
          colors = c("#ccaf1d", "#0B701E", "#FF7300", "#120C94", "#db2545", "#009BB3", "#EBEBEB")
        )
        sunburst_Control = 
          plot_ly(
            controlData, 
            labels = ~labels, parents = ~parents, 
            values = ~values,# ids = ~ids, 
            type = "sunburst"
          ) %>%
          layout(colorway  = ~colors,
                 margin = list(l = 0, r = 0, b = 0, t = 0))
        sunburst_Control$sizingPolicy$padding <- "0"
        plots_Plotly[["sunburst_Control"]] <- sunburst_Control
        sunburst_Control
      })
      
      ###### Section 2 Compartive Stat ######
      ###### . 2.1 Module-Module Preservation######
      output$moduleModuleComparison <- renderPlotly({
        df = nodesAll %>%
          mutate(presence = 1) %>%
          pivot_wider(names_from = colors_Control, values_from = presence, values_fn = list(presence = sum), 
                      id_cols = colors_Disease) %>%
          mutate(colors_Disease = ifelse(is.na(colors_Disease), "Not Found", colors_Disease)) %>%
          column_to_rownames("colors_Disease")
        if("NA" %in% colnames(df)){
          colnames(df)[match("NA", colnames(df))] = "Not Found"
        }
        curTitle = paste0(data1Name, " Vs ", data2Name)
        df[is.na(df)] = 0
        df_melt = df %>%
          mutate(colors_Disease = rownames(df)) %>%
          pivot_longer(names_to = "colors_Control", values_to = "count", cols = colnames(df))
        df_melt$text = paste0("<b>", data1Name, " Module Color: </b>", df_melt$colors_Disease, "<br>",
                              "<b>", data2Name, " Module Color: </b>", df_melt$colors_Control, "<br>",
                              "<b>Number of Shared Taxa: </b>", df_melt$count, "<br>")
        df_Wide = df_melt %>%
          pivot_wider(names_from = colors_Control, values_from = text, 
                      id_cols = colors_Disease) %>%
          mutate(colors_Disease = ifelse(is.na(colors_Disease), "Not Found", colors_Disease)) %>%
          column_to_rownames("colors_Disease")
        df_Wide = df_Wide[rownames(df), colnames(df)]
        modModPlot = 
          plot_ly(x = colnames(df), y = rownames(df),
                z = as.matrix(df), 
                hoverinfo = 'text',
                colors=colorRamp(c("#083D7F","#845E40","#E27516")),
                colorbar=list(title="Count",len=1),
                text = as.matrix(df_Wide),
                type = "heatmap") %>%
          layout(
            title = paste0(curTitle, " Module Comparison"),
            xaxis = list(
              title = data2Name,
              showspikes = TRUE,
              spikethickness = 1
            ),
            yaxis = list(
              title = data1Name,
              showspikes = TRUE,
              spikethickness = 1,
              spikedash = "10px"
            ),
            hovermode  = 'x'
          ) %>%
          add_annotations(x = df_melt$colors_Control,
                          y = df_melt$colors_Disease,
                          text = df_melt$count, 
                          showarrow = FALSE, 
                          font = list(
                            color = "white"
                          ),
                          ax = 20,
                          ay = -20)
        plots_Plotly[["modModPlot"]] <- modModPlot
        modModPlot
      })
      
      ###### . 2.2 Module Preservation Table ###### 
      output$modulePreservation <- renderPlotly({
        modulePlotData_Wide$text = paste0("<b>Module Color: </b>", modulePlotData_Wide$color, "<br>",
                                          "<b>Module Size: </b>", modulePlotData_Wide$moduleSize)
        modPreservationPlot =
          plot_ly(data = modulePlotData_Wide,
                x = ~ZSummary, y = ~medianRank,
                
                marker = list(line = list(color = "#0E1012"),
                              size = ~moduleSize,
                              color = ~color, 
                              colors = ~color),
                type = "scatter", hovertemplate = ~text) %>%
          layout(
            shapes = list(
              ## Preserved Region ZSummary < 10 and median rank > 8 
              list(type = "rect",
                   fillcolor = "#F3C6C9", 
                   line = list(color = "#B59395"), 
                   opacity = 0.3,
                   x0 = 5, 
                   x1 = max(C3NAObj$modulePreservation$modulePlotData_Wide$ZSummary
                            [which(is.finite(C3NAObj$modulePreservation$modulePlotData_Wide$ZSummary))])+2, 
                   xref = "x",
                   y0 = 0, y1 = 8, yref = "y"),
              ## Non Preserved Region ZSummary < 10 and median rank > 8 
              list(type = "rect",
                   fillcolor = "#F0EBB9", 
                   line = list(color = "#9AB2B8"),
                   opacity = 0.3,
                   x0 = 0, x1 = 10, xref = "x",
                   y0 = 8, 
                   y1 = max(C3NAObj$modulePreservation$modulePlotData_Wide$medianRank
                            [which(is.finite(C3NAObj$modulePreservation$modulePlotData_Wide$medianRank))])+2, 
                   yref = "y"),
              ## Intermediately Preserved Region 
              list(type = "rect",
                   fillcolor = "#CEEDF5", 
                   line = list(color = "#9AADB8"),
                   opacity = 0.3,
                   x0 = 10, 
                   x1 = max(C3NAObj$modulePreservation$modulePlotData_Wide$ZSummary
                            [which(is.finite(C3NAObj$modulePreservation$modulePlotData_Wide$ZSummary))])+2, 
                   xref = "x",
                   y0 = 8, 
                   y1 = max(C3NAObj$modulePreservation$modulePlotData_Wide$medianRank
                            [which(is.finite(C3NAObj$modulePreservation$modulePlotData_Wide$medianRank))])+2, 
                   yref = "y")
            ),
            annotations = list(
              list(
                x = 0.02,
                y = 0.975,
                text = "Non-Preserved",
                xref = "paper", yref = "paper",
                showarrow = FALSE, 
                font = list(color = "#A3A26D", size = 13)
              ),
              list(
                x = 0.975,
                y = 0.025,
                text = "Strongly Preserved",
                xref = "paper", yref = "paper",
                showarrow = FALSE, 
                font = list(color = "#A65358", size = 13)
              ),
              list(
                x = 0.975,
                y = 0.975,
                text = "Intermediately Preserved",
                xref = "paper", yref = "paper",
                showarrow = FALSE, 
                font = list(color = "#676D99", size = 13)
              )
            ), 
            xaxis = list(title = "<i>Z<sub>Summary</sub></i>"),
            yaxis = list(title = "<i>medianRank</i>")
          )
        
        plots_Plotly[["modPreservationPlot"]] <- modPreservationPlot
        modPreservationPlot
      })
      
      ###### . 2.3 Click Taxa ######
      output$clickedModuleInfo <- DT::renderDataTable({
        tmp <- event_data("plotly_click", priority = "event")
        closestZSum = which.min(abs(modulePlotData_Wide$ZSummary-tmp$x))
        matchedZSum = modulePlotData_Wide$ZSummary[closestZSum]
        closestMedRank = which.min(abs(modulePlotData_Wide$medianRank-tmp$y))
        matchedMedRank = modulePlotData_Wide$medianRank[closestMedRank]
        tempMod = subset(modulePlotData_Wide, ZSummary == matchedZSum & medianRank == matchedMedRank)
        if(nrow(tempMod)>0){
          selectedColor = nodesAll %>%
            filter(colors_Control == tempMod$color[1]) %>%
            dplyr::select(TaxaName, colors_Disease, colors_Control) %>%
            dplyr::rename("Comparison Color" = colors_Disease,
                          "Reference Color" = colors_Control)
          DT::datatable(selectedColor,
                        options = list(
                          rownames = FALSE,
                          autoWidth = TRUE,
                          searching = TRUE,
                          pageLength = 7
                        ))
        } else {
          placeholderData = data.frame(
            "Instruction" = paste0("Please click on the Module in the \n", 
                                   "'Module Preservation Statistics' plot\n to ", 
                                   "display the members of the module, which is \n", 
                                   "based on the modules from Phenotype: '", data2Name, 
                                   "'. (Left Figure)")
          )
          rownames(placeholderData) = ""
          datatable(
            placeholderData,              
            rownames = FALSE,
            options = list(
              dom = 't',
              autoWidth = TRUE,
              searching = FALSE,
              pageLength = 1
            )
          )
        }
      })
      
      ###### Section 3 Extract Specific Network ######
      ###### . 3.1 Taxa Selector ###### 
      output$taxaSelector <- renderUI({
        multiInput(
          inputId = "taxaSelector",
          label = "Selected Taxa:", 
          choices = NULL, 
          choiceNames = nodesAll$TaxaName, 
          choiceValues = nodesAll$TaxaName
        )
      })
      
      ###### . 3.2 Intra Modular Corr Table ######
      observeEvent(input$taxaSelector, {
        matchedID = input$taxaSelector
        values[["panel1SelectedIDs"]] = matchedID
        
        updatePickerInput(
          session, 
          inputId = "barPlotTaxaSelector",
          choices = c("All", matchedID)
        )
        
        updatePickerInput(
          session, 
          inputId = "taxaModularTaxaSelector",
          choices = c(matchedID)
        )
      })

      output$intramodularCorTable <- DT::renderDataTable({
        if(length(input$taxaSelector) > 0){
          matchedID = input$taxaModularTaxaSelector
          if(length(matchedID) > 0){
            outputDataTable = subset(sparccP_Filtered_rbind, source %in% matchedID | target %in% matchedID) %>%
              filter(cor >= corCut & fdr <= fdr) %>%
              group_by(Diagnosis, clusterID) %>% filter(clusterID != "Inter-Modular") %>%
              left_join(nodesAbbre, by = c("clusterID" = "ClusterID")) %>% ungroup() %>%
              rename(Colors = colors) %>%
              group_by(Diagnosis, Colors) %>%
              dplyr::summarise(n = n()) %>%
              dplyr::rename("Number of Important Correlations" = "n")
            if(nrow(outputDataTable) > 0){
              datatable(
                outputDataTable,
                rownames = FALSE, 
                options = list(
                  autoWidth = TRUE,
                  searching = FALSE,
                  pageLength = 7
                )
              )
            } else {
              placeholderData = data.frame(
                "Information" = paste0("No Significant Taxa-Taxa Correlations Found. ", 
                                                     "Please select another taxon.")
              )
              rownames(placeholderData) = ""
              datatable(
                placeholderData,              
                rownames = FALSE,
                options = list(
                  dom = 't',
                  autoWidth = TRUE,
                  searching = FALSE,
                  pageLength = 1
                )
              )
            }
          }
        }else {
          placeholderData = data.frame(
            "Information" = paste0("No taxon selected from the 'Taxa Selector' on the Left. ")
          )
          rownames(placeholderData) = ""
          datatable(
            placeholderData,              
            rownames = FALSE,
            options = list(
              dom = 't',
              autoWidth = TRUE,
              searching = FALSE,
              pageLength = 1
            )
          )
        }
      })
          
      output$barPlotTaxaSelector <- renderUI({
        pickerInput(
          inputId = "barPlotTaxaSelector",
          label = "Taxa Selector:", 
          choices = c("All"), 
          selected = "All"
        )
      })
      
      output$taxaModularTaxaSelector <- renderUI({
        pickerInput(
          inputId = "taxaModularTaxaSelector",
          label = "Taxa Selector:", 
          choices = NULL
        )
      })
      
      ###### . 3.3 Bar Plot ######
      observeEvent(input$barPlotTaxaSelector, {
        req(input$barPlotPhenoSelector)
        output$barPlot <- renderPlotly({
          if(input$barPlotTaxaSelector == "All"){
            nodesAll_Pheno = nodesAll_PhenoOri
            nodesAll_Pheno$`Taxonomic Levels` = 
              factor(nodesAll_Pheno$`Taxonomic Levels`,
                     levels = c("Phylum", "Class", "Order", "Family", "Genus", 
                                "Species", "Unselected"))
            if(input$barPlotPhenoSelector == phenotypes[1]){
              curPheno = "ClusterID_Disease"
            } else{
              curPheno = "ClusterID_Control"
            }
            nodesAll_Pheno$`Cluster ID` = nodesAll_Pheno[, curPheno]
            barPlotSelected =
              plot_ly(nodesAll_Pheno, 
                    x = ~`Cluster ID`,
                    y = ~`Number of Unique Taxa`,
                    type = "bar", 
                    hovertext = ~description,
                    color = ~`Taxonomic Levels`,
                    colors = c("#ccaf1d", "#0B701E", "#FF7300", "#120C94", 
                               "#db2545", "#009BB3", "#EBEBEB"),
                    marker = list(line = list(color = '#565656',
                                              width = 1))
              ) %>%
              layout(
                barmode = "relative",
                xaxis = list(title = paste0("Module ID from Phenotype: ", 
                                            input$barPlotPhenoSelector)),
                yaxis = list(title = "Number of Unique Taxa per Module")
              )
            plots_Plotly[["barPlotAll"]] <- barPlotSelected
            barPlotSelected
          } else {
            nodesAll_Pheno = nodesAll_PhenoOri
            matchedID = input$barPlotTaxaSelector
            print(matchedID)
            matchedTaxonomicLevel = subset(nodesAll_Pheno, TaxaName == matchedID)$`Taxonomic Levels`
            print(matchedTaxonomicLevel)
            
            relatedTaxaName = subset(nodesAll_Pheno, TaxaName == matchedID)
            nodesAll_Pheno$`Taxonomic Levels` = ifelse(nodesAll_Pheno[, matchedTaxonomicLevel] == matchedID,
                                                       as.character(nodesAll_Pheno$`Taxonomic Levels`), "Unselected")
            nodesAll_Pheno$`Taxonomic Levels` = 
              factor(nodesAll_Pheno$`Taxonomic Levels`,
                     levels = c("Phylum", "Class", "Order", "Family", "Genus", 
                                "Species", "Unselected"))
            if(input$barPlotPhenoSelector == phenotypes[1]){
              curPheno = "ClusterID_Disease"
            } else{
              curPheno = "ClusterID_Control"
            }
            nodesAll_Pheno$`Cluster ID` = nodesAll_Pheno[, curPheno]
            print(head(nodesAll_Pheno))
            
            barPlotAll = 
              plot_ly(nodesAll_Pheno, 
                    x = ~`Cluster ID`,
                    y = ~`Number of Unique Taxa`,
                    type = "bar", 
                    hovertext = ~description,
                    color = ~`Taxonomic Levels`,
                    colors = c("#ccaf1d", "#0B701E", "#FF7300", "#120C94", 
                               "#db2545", "#009BB3", "#EBEBEB"),
                    marker = list(line = list(color = '#565656',
                                              width = 1))
              ) %>%
                layout(
                  barmode = "relative",
                  xaxis = list(title = paste0("Module ID from Phenotype: ", 
                                              input$barPlotPhenoSelector)),
                  yaxis = list(title = "Number of Unique Taxa per Module")
                )
            plots_Plotly[["barPlotAll"]] <- barPlotAll
            barPlotAll
          }
        })
      })
      
      ## Bar plot placeholder
      output$barPlot <- renderPlotly({
        plotly_empty(type = "scatter", mode = "markers") %>%
          config(
            displayModeBar = FALSE
          ) %>%
          layout(
            title = list(
              text = paste0("Please click on the filter button on the top left."),
              yref = "paper",
              y = 0.5
            )
          )
      })
      
      ##### Section 4 Network Plots and Individual Violin Plot #####
      networkMultiListener <- reactive({
        list(input$refreshPanel1Section3, input$functionSelector)
      })
      
      ##### . 4.1 Network Plot #####
      observeEvent(networkMultiListener(), {
        req(input$nodesMode)
        req(input$corMode)
        req(input$nodeBorderWidth)
        req(input$nodeLabelSize)
        req(input$edgeLabelSize)
        req(input$samePhylo)
        req(input$DASelector)
        req(input$corSelector)
        
        output$SelectedTaxaPlots = visNetwork::renderVisNetwork({
          ## Obtain all the matched ID within 
          curCorMode = input$corMode
          curNodeBorderWidth = input$nodeBorderWidth
          curNodeLabelSize = input$nodeLabelSize
          curEdgeLabelSize = input$edgeLabelSize
          curTaxaLvls = input$taxaLvls
          cursamePhylo = input$samePhylo
          curCorFilter = input$corSelector
          
          if(!is.null(values[["panel1SelectedIDs"]])){
            if(input$nodesMode == "selectedTaxa"){
              print("Single Taxa")
              relatedTaxa = values[["panel1SelectedIDs"]]
            } else if (input$nodesMode == "relatedTaxa"){
              print("Related Taxa")
              corCut = curCorFilter
              matchedID = values[["panel1SelectedIDs"]]
              filteredData = subset(sparccP_Filtered_rbind, source %in% matchedID | target %in% matchedID) %>%
                filter(cor >= corCut) %>% filter(Module != "Inter-Modular") 
              clusterIDsDisease = subset(filteredData, Diagnosis == data1Name) %>% pull(clusterID) %>% unique()
              clusterIDsControl = subset(filteredData, Diagnosis == data2Name) %>% pull(clusterID) %>% unique()
              filteredData_Disease = subset(filteredData, Diagnosis == data1Name & clusterID %in% clusterIDsDisease)
              filteredData_Control = subset(filteredData, Diagnosis == data2Name & clusterID %in% clusterIDsControl)
              
              relatedTaxa = unique(c(filteredData_Disease$source, 
                                     filteredData_Disease$target, 
                                     filteredData_Control$source, 
                                     filteredData_Control$target))
            } else {
              print("All Taxa in Related Modules")
              matchedID = values[["panel1SelectedIDs"]]
              filteredData = subset(sparccP_Filtered_rbind, source %in% matchedID | target %in% matchedID) %>%
                filter(cor >= corCut) %>% filter(Module != "Inter-Modular") 
              clusterIDsDisease = subset(filteredData, Diagnosis == data1Name) %>% pull(clusterID) %>% unique()
              clusterIDsControl = subset(filteredData, Diagnosis == data2Name) %>% pull(clusterID) %>% unique()
              filteredData_Disease = subset(sparccP_Filtered_rbind, Diagnosis == data1Name & clusterID %in% clusterIDsDisease) %>%
                filter(cor >= corCut)
              filteredData_Control = subset(sparccP_Filtered_rbind, Diagnosis == data2Name & clusterID %in% clusterIDsControl) %>%
                filter(cor >= corCut)
              relatedTaxa = unique(c(filteredData_Disease$source, 
                                     filteredData_Disease$target, 
                                     filteredData_Control$source, 
                                     filteredData_Control$target))
            }
            print("1")
            if(input$nodesMode == "selectedTaxa"){
              sparccP_Filtered_Combined_Filtered = subset(sparccP_Filtered_Combined, source %in% relatedTaxa | target %in% relatedTaxa) %>% 
                filter(cor_Disease >= corCut | cor_Control >= corCut) %>%
                filter(!is.na(clusterID_Disease) | !is.na(clusterID_Control)) %>%
                filter(!(source %in% removedTaxaTable$removedTaxa)) %>%
                filter(!(target %in% removedTaxaTable$removedTaxa)) 
            } else{
              sparccP_Filtered_Combined_Filtered = subset(sparccP_Filtered_Combined, source %in% relatedTaxa & target %in% relatedTaxa) %>% 
                filter(cor_Disease >= corCut | cor_Control >= corCut) %>%
                filter(!is.na(clusterID_Disease) | !is.na(clusterID_Control)) %>%
                filter(!(source %in% removedTaxaTable$removedTaxa)) %>%
                filter(!(target %in% removedTaxaTable$removedTaxa)) 
            }

            print("2")
            ## Calculate Edge and Nodes Table
            nodesTable = nodesAll %>%
              filter(TaxaName %in% relatedTaxa) %>%
              mutate(color.background = colors_Disease, 
                     color.border = colors_Control,
                     borderWidth = as.integer(curNodeBorderWidth),
                     id = TaxaName,
                     label = TaxaName,
                     font.size = as.integer(curNodeLabelSize),
                     shape = "dot"
              ) %>%
              mutate(shape = ifelse(TaxaName %in% values[["panel1SelectedIDs"]], "star", shape)) %>%
              mutate(title = paste0("TaxaName: ", TaxaName, "<br><hr>",
                                    "Disease Type: ", data1Name, "<br>",
                                    "Disease Color Module: ", colors_Disease, "<br><hr>",
                                    "Control Type: ", data2Name, "<br>",
                                    "Control Color Module: ", colors_Control, "<br><hr>"
                                    ),
                     size = 10
                     ) %>%
              mutate(color.background = ifelse(color.background == "grey60", "#808080", color.background),
                     color.border = ifelse(color.border == "grey60", "#808080", color.border))
            if(!is.null(C3NAObj$C3NA_Wilcoxon)){
              nodesTable = nodesTable %>%
                left_join(C3NAResults, 
                        by = "TaxaName", all.x = TRUE)
            }
            print("3")
            ## Filter Vertices based on taxonomic levels
            nodesTable = nodesTable %>%
              rowwise() %>%
              mutate(taxaLvls_Abbre = substring(TaxaName, 1, 1)) %>%
              filter(taxaLvls_Abbre %in% curTaxaLvls) %>%
              group_by(TaxaName) %>% filter(row_number() == 1)
            print("4")
            ## Adding Differential Abundance Analysis Results if Avaliable
            ###### !!!!!!DA  ######
            print(input$DASelector)
            print(length(daMethods))
            if(length(daMethods) > 1 & input$DASelector != "None"){
              daData_temp = daData[, c("TaxaName", input$DASelector)]
              colnames(daData_temp)[2] = "DA"
              nodesTable = nodesTable %>%
                left_join(daData_temp, 
                          by = "TaxaName", all.x = TRUE) %>%
                as.data.frame() %>%
                mutate(DA = as.character(DA)) %>%
                mutate(DA = ifelse(is.na(DA), FALSE, DA)) %>%
                mutate(title = paste0(title, 
                                      "Differential Abundance: ", DA, "<br><hr>"), 
                       shape = ifelse(DA == TRUE, "square", shape)) %>%
                replace(is.na(.), values = FALSE)
              lnodes <- data.frame(
                label = c("DA", "Selected", "Others"),
                shape = c("square", "star", "dot"),
                shadow = c(FALSE, FALSE, FALSE)
              )
            } else {
              print("Only Star")
              lnodes <- data.frame(
                label = c("Selected", "Others"),
                shape = c("star", "dot"),
                shadow = c(FALSE, FALSE)
              )
            }
            print("5")
            ## C3NA Results
            print("C3NA Results")
            if(!is.null(C3NAObj$C3NA_Wilcoxon)){
              print("C3NA_Wilcoxon")
              nodesTable = nodesTable %>%
                mutate(C3NA = ifelse(is.na(C3NA), FALSE, C3NA)) %>%
                mutate(title = paste0(title, 
                                      "C3NA Wilcoxon: ", C3NA, ""), 
                       shape = ifelse(C3NA == TRUE, "triangle", shape)) 
              
              ## Update Legend
              selectedFuncs = values[["selectedFuncs"]]
              avaliableShape = c("diamond", "triangle", "ellipse")
              lnodesV2 <- data.frame(
                label = "C3NA",
                shape = "triangle",
                shadow = FALSE
              )
              lnodes = rbind(lnodesV2, lnodes)
            } else {
              print("No C3NA")
              lnodes = lnodes %>%
                filter(!(shape %in% c("triangle")))
            }
            print("6")
            ## Remove Duplicates Caused by inputTaxaName
            if("inputTaxaName" %in% colnames(nodesTable)){
              nodesTable = nodesTable %>%
                select(-inputTaxaName) %>%
                distinct()
            }
            
            print("7")
            # Check if both C3NA & DA 
            if(!is.null(C3NAObj$C3NA_Wilcoxon) & (length(daMethods) > 1 & input$DASelector != "None")){
              nodesTable = nodesTable %>%
                mutate(shape = ifelse(C3NA == TRUE & DA == TRUE, "diamond", shape))
              
              lnodesV2 <- data.frame(
                label = "Both",
                shape = "diamond",
                shadow = FALSE
              )
              lnodes = rbind(lnodesV2, lnodes)
            }
            
            ledges <- data.frame(
              label = c("Disease-Only", "Control-Only", "Shared"),
              color = c( "#77BBEB", "#EB7874", "#CCCFCB"),
              dashes = c(FALSE, FALSE, FALSE),
              arrows.from.type = c("dot", "dot", "dot"),
              arrows.to.type = c(NA, NA, NA),
              font.align = "top"
            )
            
            print("8")
            # Edge
            edgesTable = sparccP_Filtered_Combined_Filtered %>%
              filter(cor_Disease >= corCut | cor_Control >= corCut) %>%
              mutate(cor_Disease_Filtered = ifelse(cor_Disease < corCut | is.na(cor_Disease) | Module_Disease == "Inter-Modular" | is.na(Module_Disease), 0, cor_Disease),
                     cor_Control_Filtered = ifelse(cor_Control < corCut | is.na(cor_Control) | Module_Control == "Inter-Modular" | is.na(Module_Control), 0, cor_Control)) %>%
              mutate(from = source, 
                     to = target,
                     color = case_when(
                       (Module_Control=="Inter-Modular" | is.na(Module_Control)) & Module_Disease=="Intra-Modular" ~ "#77BBEB70",
                       (Module_Disease=="Inter-Modular" | is.na(Module_Disease)) & Module_Control=="Intra-Modular" ~ "#EB787470",
                       (Module_Disease=="Intra-Modular" & Module_Control=="Intra-Modular") ~ "#CCCFCB70"
                     )                          ,
                     label = case_when(
                       (Module_Disease=="Intra-Modular" & Module_Control=="Intra-Modular") ~ paste0(round(cor_Disease_Filtered, 2), " | ", round(cor_Control_Filtered, 2)),
                       (Module_Control=="Inter-Modular" | is.na(Module_Control)) & Module_Disease=="Intra-Modular" ~ as.character(round(cor_Disease_Filtered, 2)),
                       (Module_Disease=="Inter-Modular" | is.na(Module_Disease)) & Module_Control=="Intra-Modular" ~ as.character(round(cor_Control_Filtered, 2))
                     ),
                     font.color = case_when(
                       (Module_Disease=="Intra-Modular" & Module_Control=="Intra-Modular") ~ "#CCCFCB85",
                       (Module_Control=="Inter-Modular" | is.na(Module_Control)) & Module_Disease=="Intra-Modular" ~ "#77BBEB85",
                       (Module_Disease=="Inter-Modular" | is.na(Module_Disease)) & Module_Control=="Intra-Modular" ~ "#EB787485"
                     )
              ) %>%
              filter(!(cor_Disease <= corCut & Module_Control == "Inter-Modular") & 
                       !(cor_Control <= corCut & Module_Disease == "Inter-Modular")) %>%
              mutate(font.size = as.integer(curEdgeLabelSize)) %>%
              left_join(consensusProp_Disease, by = "taxaComboName", all.x = TRUE) %>%
              left_join(consensusProp_Control , by = "taxaComboName", all.x = TRUE) %>%
              mutate(title = paste0("Connection between: <br>", from, " & ", to, "<hr>",
                                    "Disease Consensus:    ", DiseaseConsensusProp, "<br>", 
                                    "Control Consensus:    ", ControlConsensusProp, "<br><hr>",
                                    "Disease Cor:    ", round(cor_Disease,3), "<br>", 
                                    "Control Cor:    ", round(cor_Control,3), "<br><hr>",
                                    "Disease Module: ", Module_Disease, "<br>",
                                    "Control Module: ", Module_Control)) %>%
              left_join(samePhyloTable, by = "taxaComboName", all.x = TRUE) %>%
              mutate(samePhylo = ifelse(is.na(samePhylo), FALSE, samePhylo))
            print(dim(edgesTable))
            save(edgesTable, nodesTable,
                 file = "../tempEdgeTable.rdata")
            
            ## Re-format based on Current mode
            if(length(curCorMode) < 3){
              if(!("diseaseOnly" %in% curCorMode)){
                ## Update Edges
                edgesTable = edgesTable %>%
                  filter(color != "#77BBEB70")
                ## Update Nodes
                nodesTable = nodesTable %>%
                  filter(TaxaName %in% unique(unlist(edgesTable[, c("source", "target")])))
              }
              if(!("controlOnly" %in% curCorMode)){
                ## Update Edges
                edgesTable = edgesTable %>%
                  filter(color != "#EB787470")
                ## Update Nodes
                nodesTable = nodesTable %>%
                  filter(TaxaName %in% unique(unlist(edgesTable[, c("source", "target")])))
              }
              if(!("diseaseOnly" %in% curCorMode)){
                ## Update Edges
                edgesTable = edgesTable %>%
                  filter(color != "#EB787470")
                ## Update Nodes
                nodesTable = nodesTable %>%
                  filter(TaxaName %in% unique(unlist(edgesTable[, c("source", "target")])))
              }
            }
            print(dim(edgesTable))

            print("9")
            if(cursamePhylo == "Remove"){
              edgesTable$samePhyloV2 = ifelse(is.na(edgesTable$samePhylo), FALSE, TRUE)
              edgesTable = edgesTable %>% dplyr::filter(samePhyloV2 != TRUE)
            }
            values[["nodesTable"]] = nodesTable
            
            print("10")
            
            networkPlot = 
              visNetwork(nodes = nodesTable, edges = edgesTable
              ) %>%
                visNodes(
                  shadow = list(enabled = FALSE)
                ) %>%
                visOptions(highlightNearest = TRUE) %>%
                visInteraction(navigationButtons = TRUE) %>%
                visEdges(smooth = FALSE) %>%
                visLegend(width=0.1, addEdges = ledges,
                          useGroups = FALSE, addNodes = lnodes,
                          zoom = FALSE, main = "Legend") %>% 
                visPhysics(solver = "forceAtlas2Based", 
                           forceAtlas2Based = list(avoidOverlap = 1, 
                                                   centralGravity = 0.001),
                           minVelocity = 20) %>%
                visEvents(click = "function(nodes){
                                Shiny.onInputChange('click', nodes.nodes[0]);
                                Shiny.onInputChange('node_selected', nodes.nodes.length);
                            ;}"
                ) 
            plots_Plotly[["networkPlot"]] <- networkPlot
            networkPlot
          }
        })
      })
  
      ##### . 4.2 Violin Click Plot #####
      observeEvent(input$click, {
        selectedID = input$click
        taxaLvls <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
        taxaHeaders <- c("p", "c", "o", "f", "g", "s")
        splitName = strsplit(selectedID, "_") 
        selectedTaxaName = selectedID
        selectedTaxaLvl = taxaLvls[match(splitName[[1]][1], taxaHeaders)]
        selectedTaxaNameFull = gsub(paste0("^", splitName[[1]][1], "_"), "", selectedTaxaName)
        curDiseaseCount = data.frame(
          count = as.numeric(curDisease_OriCountTable[selectedTaxaName, ]),
          sampleID = colnames(curDisease_OriCountTable[selectedTaxaName, ]),
          condition = data1Name
        )
        curControlCount = data.frame(
          count = as.numeric(curControl_OriCountTable[selectedTaxaName, ]),
          sampleID = colnames(curControl_OriCountTable[selectedTaxaName, ]),
          condition = data2Name
        )

        combinedCountTable = rbind(curControlCount, curDiseaseCount)
        combinedCountTable$hoverText = paste0("<b>Sample ID: </b>", combinedCountTable$sampleID, "<br>",
                                              "<b>Count: </b>", combinedCountTable$count, "<br>")

        output$countDistributionPlot <- renderPlotly({
          singleTaxonPlot = 
            plot_ly(data = combinedCountTable, 
                    x = ~condition, 
                    y = ~log(count), 
                    text  = ~hoverText,
                    points = "all", 
                    split = ~condition,
                    line = list(color = 'rgb(7,40,89)'),
                    marker = list(line = list(color = 'rgb(7,40,89)',
                                              width = 2)),
                    
                    type = 'violin',
                    box = list(
                      visible = T
                    ),
                    meanline = list(
                      visible = T
                    )) %>%
              layout(
                xaxis = list(title="<b>Condition</b>"), 
                yaxis = list(title="<b>Log<sub>2</sub>-Transformed Count</b>"),
                title = paste0(selectedTaxaLvl, " | ", selectedTaxaNameFull)
              )
          
          plots_Plotly[["singleTaxonPlot"]] <- singleTaxonPlot
          singleTaxonPlot
        })

        
        if(is.null(values[["click1"]]) & is.null(values[["click2"]])){
          values[["click1"]] = selectedTaxaName
        } else if (!is.null(values[["click1"]]) & is.null(values[["click2"]])){
          if(selectedTaxaName != values[["click1"]]){
            print("Update second Taxon")
            values[["click2"]] = selectedTaxaName
          }
        } else if (!is.null(values[["click1"]]) & !is.null(values[["click2"]])){
          if(values[["click1"]] != values[["click2"]] &
             values[["click1"]] != selectedTaxaName &
             values[["click2"]] != selectedTaxaName){
            values[["click1"]] = values[["click2"]]
            values[["click2"]] = selectedTaxaName
          }
        }

        if(!is.null(values[["click1"]]) & !is.null(values[["click2"]])){
          taxaName1 = values[["click1"]]
          taxaName2 = values[["click2"]]
          curControlCount1 = data.frame(
            count = as.numeric(curControl_OriCountTable[values[["click1"]], ]),
            sampleID = colnames(curControl_OriCountTable[values[["click1"]], ]),
            condition = data2Name,
            taxaName = taxaName1
          )
          curDiseaseCount1 = data.frame(
            count = as.numeric(curDisease_OriCountTable[values[["click1"]], ]),
            sampleID = colnames(curDisease_OriCountTable[values[["click1"]], ]),
            condition = data1Name,
            taxaName = taxaName1
          )
          combinedCountTable2 = rbind(curControlCount1, curDiseaseCount1)
          
          curControlCount2 = data.frame(
            count = as.numeric(curControl_OriCountTable[values[["click2"]], ]),
            sampleID = colnames(curControl_OriCountTable[values[["click2"]], ]),
            condition = data2Name,
            taxaName = taxaName2
          )
          curDiseaseCount2 = data.frame(
            count = as.numeric(curDisease_OriCountTable[values[["click2"]], ]),
            sampleID = colnames(curDisease_OriCountTable[values[["click2"]], ]),
            condition = data1Name,
            taxaName = taxaName2
          )
          combinedCountTable2 = rbind(combinedCountTable2, curControlCount2)
          combinedCountTable2 = rbind(combinedCountTable2, curDiseaseCount2)
          
          combinedCountTable2$hoverText = paste0("<b>Sample ID: </b>", combinedCountTable2$sampleID, "<br>",
                                                 "<b>Count: </b>", combinedCountTable2$count, "<br>")
          output$countDistributionPlot_Compare <- renderPlotly({
            fig <- plot_ly(type = 'violin')
            conditions = unique(combinedCountTable2$condition)
            taxaNames = c(values[["click2"]], values[["click1"]])
            df = combinedCountTable2
            showLegend <- c(T,F)
            
            for(i in seq(2)){
              fig <- add_trace(
                fig,
                x = df$taxaName[df$taxaName == taxaNames[i] & df$condition == conditions[1]],
                y = log(df$count[df$taxaName == taxaNames[i] & df$condition == conditions[1]]),
                text = df$hoverText[df$taxaName == taxaNames[i] & df$condition == conditions[1]],
                legendgroup = conditions[1],
                scalegroup = conditions[1],
                name = conditions[1],
                side = 'negative',
                box = list(
                  visible = T
                ),
                points = 'all',
                pointpos = -0.4,
                scalemode = 'count',
                meanline = list(
                  visible = T
                ),
                color = I("#B33A3B"),
                marker = list(
                  line = list(
                    width = 2,
                    color = "#B33A3B"
                  )
                ),
                showlegend = showLegend[i]
              )
              fig <- fig %>%
                add_trace(
                  x = df$taxaName[df$taxaName == taxaNames[i] & df$condition == conditions[2]],
                  y = log(df$count[df$taxaName == taxaNames[i] & df$condition == conditions[2]]),
                  text = df$hoverText[df$taxaName == taxaNames[i] & df$condition == conditions[2]],
                  legendgroup = conditions[2],
                  scalegroup = conditions[2],
                  name = conditions[2],
                  side = 'positive',
                  box = list(
                    visible = T
                  ),
                  points = 'all',
                  pointpos = 0.4,
                  scalemode = 'count',
                  meanline = list(
                    visible = T
                  ),
                  color = I("#2688B3"),
                  marker = list(
                    line = list(
                      width = 2,
                      color = "#2688B3"
                    )
                  ),
                  showlegend = showLegend[i]
                )
            }
            fig <- layout(
              fig,
              xaxis = list(title="<b>Condition</b>"),
              yaxis = list(title="<b>Log<sub>2</sub>-Transformed Count</b>"),
              violingap = 0,
              violingroupgap = 0,
              violinmode = 'overlay',
              legend = list(
                tracegroupgap = 0
              )
            )
            
            plots_Plotly[["twoTaxaPlot"]] <- fig
            fig
          })
        }
      })
      
      ##### 5.0 Download Panel ###### 
      output$plotSelector <- renderUI({
        if(input$plotRadioButtons %in% c("HTML")){
          ## Obtain the list of plots available
          avaliablePlots <- sort(names(plots_Plotly))
          listOfImportantData <- data.frame(
            original = c("sunburst_Disease", "sunburst_Control",   
                         "modModPlot", "modPreservationPlot", 
                         "barPlotAll", "networkPlot",        
                         "singleTaxonPlot", "twoTaxaPlot"),
            full = c("Comparison Sunburst Plot", "Reference Sunburst Plot", 
                     "Module Membership Comparison Plot", "Module Preservation Statistics Plot",
                     "Phenotype-Based Bar Plot", "Network Plot", 
                     "Single Taxon Violin Plot", "Two Taxa Vilin Plot"
                     ))
          checker <- avaliablePlots %in% listOfImportantData$original
          tempData <- avaliablePlots[checker]
          listData <- listOfImportantData %>%
            filter(original %in% tempData)
          
          ## Obtain the list of plots available
          if(nrow(listData) < 1){
            h4("No Plot was generated from Panel 1. Please check your data.")
          } else {
            radioGroupButtons(
              inputId = "plotSelector",
              choiceNames = listData$full, 
              choiceValues = listData$original,
              individual = TRUE,
              checkIcon = list(
                yes = tags$i(class = "fa fa-check-square", 
                             style = "color: white"),
                no = tags$i(class = "fa fa-square-o", 
                            style = "color: black")),
              status = "custom-status-downloadPlots"
            )
          }
        }
      })
      
      ## Observe the download button and handle it
      selectedPlotInput <- function(){
        req(input$plotSelector)
        plots_Plotly[[input$plotSelector]]
      }
      
      # create filename
      selectedPlotName <- reactive({
        name <- input$plotSelector
        if(input$plotRadioButtons == "HTML")   filename <- paste0(name, "-", Sys.Date(), ".html",  sep="")
        return(filename)
      })
      
      output$downloadPlots <- downloadHandler(
        filename = selectedPlotName,
        content = function(file) {
          name <- input$plotSelector
          
          if(input$plotRadioButtons == "HTML") {
            savedData <- selectedPlotInput()
            savedData$height = as.numeric(input$heightInput)
            savedData$width = as.numeric(input$widthInput)
            saveWidget(widget = savedData, file = file)
          }}
      )
      
      ##### 6.0 All Showmodels ######
      observeEvent(input$MS_Pheno1, {
        showModal(
          modalDialog(
            title = HTML('<b>Peek of the Avaliable Taxa Inforamation from Comparison Phentype</b>'),
            HTML("
This is the panel for the relative proportion of taxa found among the <b><i>comparison</i></b> phenotype among the Phylum, Class, Order, Family, Genus, and Species Levels. 
<hr>
<b>Notes: </b> This only includes taxa with correlations, and BH-adjusted p-values passed the thresholds defined in comparePhenotypes() command.               
"),
easyClose = TRUE
          )
        )
      })
      
      observeEvent(input$MS_Pheno2, {
        showModal(
          modalDialog(
            title = HTML('<b>Peek of the Avaliable Taxa Inforamation from Reference Phentype</b>'),
            HTML("
This is the panel for the relative proportion of taxa found among the <b><i>reference/control</i></b> phenotype among the Phylum, Class, Order, Family, Genus, and Species Levels. 
<hr>
<b>Notes: </b> This only includes taxa with correlations, and BH-adjusted p-values passed the thresholds defined in comparePhenotypes() command. 
"),
easyClose = TRUE
          )
        )
      })
      
      observeEvent(input$MS_PhenoComparison, {
        showModal(
          modalDialog(
            title = HTML('<b>Detailed Comparison between the Two Phenotypes</b>'),
            HTML("
This is an interactive table that includes all the shared and unique taxa between the comparison and reference phenotypes. <hr>
You can click the black triangle to expand the interactive table. <br>
The 'Full Taxa Name' will only display when a specific 'Taxonomic Level' is selected. 
"),
easyClose = TRUE
          )
        )
      })
      
observeEvent(input$MS_ModModComparison, {
        showModal(
          modalDialog(
            title = HTML('<b>Module-Module Taxa Comparison between Phenotypes</b>'),
            HTML("
The interactive heatmap highlights the shared taxa between the two phenotypes, the number of colors highlights the shared taxa. 
The taxa included here only includes taxa that have a correlation greater or equal to the correlation satisfied the correlation and BH-adjusted p-value cut-off values specified in comparePhenotypes() command. 
<br>
If a taxon is not found in either of the phenotypes, it will be labeled as 'Not Found.' The plot can be saved in HTML from 'Panel 2: Download Panel', or you can use the 'Download plot as a png' button on the top right corner of the plot (a default function from plotly). 
              "),
easyClose = TRUE
          )
        )
      })
      
      observeEvent(input$MS_ModPreservation, {
        showModal(
          modalDialog(
            title = HTML('<b>Module Preservation Statistics</b>'),
            HTML(
"The module preservation statistics plot shows the <i>medianRank</i> Vs. <i>Z<sub>Summary</sub></i> with the reference specified in comparePhenotypes() command.  
<br> The size of the dot represents the number of taxa in the module of the reference. "

            ),
easyClose = TRUE
          )
        )
      })
      
      observeEvent(input$MS_ModMembership, {
        showModal(
          modalDialog(
            title = HTML('<b>Module Taxa Memebrship between the Two Phenotypes</b>'),
            HTML(
"
This is an interactive table showing the taxa in each of the modules from both the disease and reference modules. In order to activate the table, you will need to click the module from the 'module preservation statistics' plot on the left.  <br> 
If the module overlapped each other, you could zoom in on the plot to click on the module.
"

            ),
easyClose = TRUE
          )
        )
      })

      observeEvent(input$MS_TaxaSelector, {
        showModal(
          modalDialog(
            title = HTML('<b>Select the Taxa You Want to Visualize on the Network Plot</b>'),
            HTML(
              "
This is an interactive module where you can click the taxa on the left column to include the selected taxon for investigation on the network plot below. <hr>
Multiple selections of taxa are allowed. <br>
If you want to remove a selected taxon, click its name on the right column. <br>
You can also search for the taxon you want by inputting the taxon in the search box. 
"

            ),
easyClose = TRUE
          )
        )
      })      
      
      observeEvent(input$MS_ModBarPlot, {
        showModal(
          modalDialog(
            title = HTML('<b>Module-based Taxa Inforamtion Bar Plot</b>'),
            HTML(
              "
<b>Note: </b> Please click on the filter button to activate the plot. 
<hr>
The box plot highlights the taxa membership per module from both phenotypes as well as its full taxonomic information, including all its parent taxonomic levels. <br>
By default, it displays all taxa, and you can also select a single taxon (selected from the Taxa Selector panel on the left). In this mode, all its children's taxonomic taxa will be shown, and all other taxa will be color in gray.  
<br>
You switch between the phenotypes and between the 'All' and individual taxa. <br>
If no selected taxa was in the dropdown menu, please random click another taxa on the left panel.
"

            ),
easyClose = TRUE
          )
        )
      })      
      
      observeEvent(input$MS_TaxaModTable, {
        showModal(
          modalDialog(
            title = HTML('<b>Number of Valid Taxa-Taxa Interactions for Selected Taxa</b>'),
            HTML(
              "
This table includes all the taxa-taxa interactions and their respecitve related module from the selected taxa. 
If no selected taxa was in the dropdown menu, please random click another taxa on the left panel.
"

            ),
easyClose = TRUE
          )
        )
      })      
      
      observeEvent(input$MS_Network, {
      showModal(modalDialog(
          title = "Interactive Network Plot Settings ", 
            tags$img(
              src = base64enc::dataURI(file = "www/NetworkLegend.png", mime = "image/png"),
              alt = "Italian Trulli",
              width = "200px"
            ),
          HTML("
The Network Plot section allows an interactive display of the selected taxa from the 'Taxa Selector' with a range of filtering settings accessible from the top left filter button. <br><b>Note:</b> <br>
<ul>
  <li>The node shape has priority shape from the top to the bottom of the legend order. E.g., if a taxon is significant from differential abundance analysis, it will be displayed as 'square' instead of 'star' even if it is a selected taxon. </li>
  <li>The border color represent the <b>Reference</b> module color, the internal color represent the <b>Comparison</b> module color as shown in the figure on the left. </li>
  <li> In the description below, we define a <i>significant connection</i> as a taxa-taxa connection that passed the correlation threshold (from both the slider input and the comparePhenotypes() function) and with a BH-adjusted p-value from comparePhenotypes(). </li>
</ul>
<hr>
<b>Node-Node Connections Mode:</b> <br>
<b>Only the selected taxa:</b> Only display the selected taxa. In this mode, only the selected taxa from the 'Taxa Selector' will render, and they will remain as nodes without connections if there are no significant connections. <br>
<b>All taxa directly related to selected taxa:</b> This is the <b><i>Suggested</i></b> mode for display. It will render all the taxa that shared significant connections with the selected taxa and intra-modular taxa that shared a significant connection with the selected taxa. <br>
<b>All taxa in related modules: </b> This is the mode to display all taxa that shared the same module membership from both the <b>Comparison</b> and <b>Reference</b> phenotypes. This mode will be likely to create multiple sub-graphs that might be harder to view. <hr>
<b>Connections from Share Common Phylogenetic Ancestor?</b> This is an optional method to remove/keep all significant connections that link between the taxa that share a common ancestor. E.g. The connection between the Family Bacteroidaceae and its children Genus <i>Bacteroides</i> will be  removed if the selection is 'Remove'. <hr>
<b>Edge Mode:</b>  The shard mode can filter the correlation based on their definition from the significant connections. <hr>
<b>Taxonomic Levels:</b> Filtering the taxa based on its taxonomic levels.<hr>
<b>Node Border Color Width, Node Label Size, Edge Label Size</b> Adjust the display setting for Nodes and Connections display settings. 


"),
          size = "l",
          
          easyClose = TRUE
        ))
      })
      
      observeEvent(input$MS_SingelTaxaViolin, {
        showModal(
          modalDialog(
            title = HTML('<b>Single Taxon Violin Plot</b>'),
            HTML(
              "This plot represents a single taxon and its respective log<sub>2</sub>-transformed count from both phenotypes. If one of them is missing, it means such taxon is missed from the respective phenotype. "            
            ),
            easyClose = TRUE
          )
        )
      })
      
      observeEvent(input$MS_TwoTaxaViolin, {
        showModal(
          modalDialog(
            title = HTML('<b>Two Taxa Violin Plot</b>'),
            HTML(
              "
Similar to the single taxon display above, this plot displays the log<sub>2</sub>-transformed count from both phenotypes for the sequentially clicked taxa. "            
            ),
            easyClose = TRUE
          )
        )
      })
      
            
      observeEvent(input$MS_DownloadPanel, {
        showModal(
          modalDialog(
            title = HTML('<b>Download Displayed Figures in Panel 1</b>'),
            HTML(
              "The download panel enables you to download the currently displayed plots from Panel 1 in PDF or PNG format. For the PNG format, the DPI is set to 300. The file should be saved to your default download folder with the plot type name followed by the date and time to avoid overwriting any plot. "            
            ),
            easyClose = TRUE
          )
        )
      })
      
      # stop App on closing the browser
      session$onSessionEnded(function() {
        stopApp()
      })
    },
    options = list(launch.browser = TRUE)
  )
}








