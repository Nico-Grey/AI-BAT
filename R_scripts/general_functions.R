

############################################### GSEA_calc_gene ##############################################

library(biomaRt)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(clusterProfiler)
library(DOSE)

gsea_pathway_dir = gsea_pathway_dir

GSEA_calc_gene <- function(gene_list,
                           DEG_list,
                           comparison = NULL,
                           species = "M", # H for human
                           genes_down = NULL,
                           genes_all = NULL,
                           GeneSets =c("GO","KEGG","DO","Hallmark","cannonicalPathways","Motifs","ImmunoSignatures"),
                           GOntology = "BP", #Alternative "MF" or "CC"
                           pCorrection = "bonferroni", # choose the p-value adjustment method
                           pvalueCutoff = 0.1, # set the unadj. or adj. p-value cutoff (depending on correction method)
                           qvalueCutoff = 0.2 # set the q-value cutoff (FDR corrected)
) {
  
  # define universe
  universe <- as.character(gene_list)
  # change symbols to ENTREZ IDs (necessary for ClusterProfiler)
  if(species == "M") {
    universe = getLDS(attributes = c("mgi_symbol"), 
                      filters = "mgi_symbol", 
                      values = universe, 
                      mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/"), 
                      attributesL = c("hgnc_symbol"), 
                      martL = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/"), 
                      uniqueRows=T)[,2]
  }else{}
  
  universe_Entrez <- bitr(universe, 
                          fromType="SYMBOL", 
                          toType="ENTREZID", 
                          OrgDb="org.Hs.eg.db")$ENTREZID
  
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  # human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",host = "https://www.ensembl.org")
  # mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",host = "https://www.ensembl.org")
  
  
  genes_up = DEG_list 
  # genes_down = DEG_down_list 
  # genes_all = DEG_all_list
  
  if(species == "M") {
    genes_up <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes_up, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)[,2]
    # genes_down <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes_down, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)[,2]
    # genes_all <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes_all, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)[,2]
  }
  
  cannonicalPathway_genes <- clusterProfiler::read.gmt(paste(gsea_pathway_dir,"c2.all.v7.5.entrez.gmt",sep=""))
  immuno_genes <- clusterProfiler::read.gmt(paste(gsea_pathway_dir,"c7.all.v7.5.entrez.gmt",sep=""))
  hallmark_genes <- clusterProfiler::read.gmt(paste(gsea_pathway_dir,"h.all.v7.5.entrez.gmt",sep=""))
  motifs <- clusterProfiler::read.gmt(paste(gsea_pathway_dir,"c3.all.v7.5.entrez.gmt",sep=""))
  
  # if(!is.null(group_by)){
  # SetIdents(object = tmp) <- group_by
  # } 
  
  # if(!is.null(genes_down)) {
  #   top_down <- genes_down
  #   entrez_down <- bitr(top_down, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  # }
  
  # if(!is.null(genes_all)) {
  #   top_all <- genes_all
  #   entrez_all <- bitr(top_all, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  # }
  
  top_up <- genes_up
  entrez_up <- bitr(top_up, fromType = "SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)$ENTREZID
  
  OrgDb = org.Hs.eg.db
  
  results <- list()
  
  
  # GO enrichment
  if("GO" %in% GeneSets){
    print("Performing GO enrichment")
    results$GOup <- as.data.frame(enrichGO(gene = entrez_up,
                                           universe = universe_Entrez,
                                           OrgDb = OrgDb,
                                           ont = GOntology,
                                           pAdjustMethod = pCorrection,
                                           pvalueCutoff  = pvalueCutoff,
                                           qvalueCutoff  = qvalueCutoff,
                                           readable      = T))
    
    if(nrow(results$GOup)>0){results$GOup$Enrichment <- paste("GO enrichment for genes upregulated in comparison ",comparison,sep="")}
    
    # if(!is.null(genes_down)) {
    #   results$GOdown <- as.data.frame(enrichGO(gene = entrez_down,
    #                                            universe = universe_Entrez,
    #                                            OrgDb = OrgDb,
    #                                            ont = GOntology,
    #                                            pAdjustMethod = pCorrection,
    #                                            pvalueCutoff  = pvalueCutoff,
    #                                            qvalueCutoff  = qvalueCutoff,
    #                                            readable      = T))
    #   if(nrow(results$GOdown)>0){results$GOdown$Enrichment <- paste("GO enrichment for genes downregulated in comparison ",comparison,sep="")}
    # }
    
  #   if(!is.null(genes_all)) {
  #     results$GOall <- as.data.frame(enrichGO(gene = entrez_all,
  #                                             universe = universe_Entrez,
  #                                             OrgDb = OrgDb,
  #                                             ont = GOntology,
  #                                             pAdjustMethod = pCorrection,
  #                                             pvalueCutoff  = pvalueCutoff,
  #                                             qvalueCutoff  = qvalueCutoff,
  #                                             readable      = T))
  #     if(nrow(results$GOall)>0){results$GOall$Enrichment <- paste("GO enrichment for genes allregulated in comparison ",comparison,sep="")}
  #   }
  }
  
  # KEGG enrichment
  if("KEGG" %in% GeneSets){
    print("Performing KEGG enrichment")
    
    org = "hsa"
    
    results$KEGGup <- as.data.frame(enrichKEGG(gene = entrez_up, 
                                               organism = org,
                                               universe = universe_Entrez, 
                                               pAdjustMethod = pCorrection,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff = qvalueCutoff))
    if(nrow(results$KEGGup)>0){results$KEGGup$Enrichment <- paste("KEGG enrichment for genes upregulated in comparison ",comparison,sep="")}
    
    # if(!is.null(genes_down)) {
    #   results$KEGGdown <- as.data.frame(enrichKEGG(gene = entrez_down, 
    #                                                organism = org,
    #                                                universe = universe_Entrez, 
    #                                                pAdjustMethod = pCorrection,
    #                                                pvalueCutoff  = pvalueCutoff,
    #                                                qvalueCutoff = qvalueCutoff))
    #   if(nrow(results$KEGGdown)>0){results$KEGGdown$Enrichment <- paste("KEGG enrichment for genes downregulated in comparison ",comparison,sep="")}
    # }
  }
  
  # DO enrichment
  if("DO" %in% GeneSets){
    print("Performing Disease Ontology enrichment")
    
    results$DOup <- as.data.frame(enrichDO(gene = entrez_up, 
                                           universe = universe_Entrez, 
                                           pAdjustMethod = pCorrection,
                                           pvalueCutoff  = pvalueCutoff,
                                           qvalueCutoff = qvalueCutoff,
                                           minGSSize     = 5,
                                           maxGSSize     = 500,
                                           readable=TRUE))
    if(nrow(results$DOup)>0){results$DOup$Enrichment <- paste("DO enrichment for genes upregulated in comparison ",comparison,sep="")}
    
    # if(!is.null(genes_down)) {
    #   results$DOdown <- as.data.frame(enrichDO(gene = entrez_down, 
    #                                            universe = universe_Entrez, 
    #                                            pAdjustMethod = pCorrection,
    #                                            pvalueCutoff  = pvalueCutoff,
    #                                            qvalueCutoff = qvalueCutoff,
    #                                            minGSSize     = 5,
    #                                            maxGSSize     = 500,
    #                                            readable=TRUE))
    #   if(nrow(results$DOdown)>0){results$DOdown$Enrichment <- paste("DO enrichment for genes downregulated in comparison ",comparison,sep="")}
    # }
  }
  
  # Hallmark enrichment
  if("Hallmark" %in% GeneSets){
    print("Performing Hallmark enrichment")
    
    results$HALLMARKup <- as.data.frame(enricher(entrez_up,
                                                 TERM2GENE=hallmark_genes,
                                                 universe = universe_Entrez,  
                                                 pAdjustMethod = pCorrection,
                                                 pvalueCutoff  = pvalueCutoff,
                                                 qvalueCutoff = qvalueCutoff))
    if(nrow(results$HALLMARKup)>0){results$HALLMARKup$Enrichment <- paste("HALLMARK enrichment for genes upregulated in comparison ",comparison,sep="")}
    
    # if(!is.null(genes_down)) {
    #   results$HALLMARKdown <- as.data.frame(enricher(entrez_down,
    #                                                  TERM2GENE=hallmark_genes,
    #                                                  universe = universe_Entrez,  
    #                                                  pAdjustMethod = pCorrection,
    #                                                  pvalueCutoff  = pvalueCutoff,
    #                                                  qvalueCutoff = qvalueCutoff))
    #   if(nrow(results$HALLMARKdown)>0){results$HALLMARKdown$Enrichment <- paste("HALLMARK enrichment for genes downregulated in comparison ",comparison,sep="")}
    # }
  }
  
  # Cannonical Pathway enrichment
  if("cannonicalPathways" %in% GeneSets){
    print("Performing Cannonical Pathway (C2) enrichment")
    
    results$cannonicalPathwaysup <- as.data.frame(enricher(entrez_up,
                                                           TERM2GENE=cannonicalPathway_genes,
                                                           universe = universe_Entrez,  
                                                           pAdjustMethod = pCorrection,
                                                           pvalueCutoff  = pvalueCutoff,
                                                           qvalueCutoff = qvalueCutoff))
    if(nrow(results$cannonicalPathwaysup)>0){results$cannonicalPathwaysup$Enrichment <- paste("Cannonical pathway enrichment for genes upregulated in comparison ",comparison,sep="")}
    
    # if(!is.null(genes_down)) {
    #   results$cannonicalPathwaysdown <- as.data.frame(enricher(entrez_down,
    #                                                            TERM2GENE=cannonicalPathway_genes,
    #                                                            universe = universe_Entrez,  
    #                                                            pAdjustMethod = pCorrection,
    #                                                            pvalueCutoff  = pvalueCutoff,
    #                                                            qvalueCutoff = qvalueCutoff))
    #   if(nrow(results$cannonicalPathwaysdown)>0){results$cannonicalPathwaysdown$Enrichment <- paste("Cannonical pathway enrichment for genes downregulated in comparison ",comparison,sep="")}
    # }
  }
  
  # Motif enrichment
  if("Motifs" %in% GeneSets){
    print("Performing Motif enrichment")
    
    results$Motifup <- as.data.frame(enricher(entrez_up,
                                              TERM2GENE=motifs,
                                              universe = universe_Entrez,  
                                              pAdjustMethod = pCorrection,
                                              pvalueCutoff  = pvalueCutoff,
                                              qvalueCutoff = qvalueCutoff))
    if(nrow(results$Motifup)>0){results$Motifup$Enrichment <- paste("TF binding motif enrichment for genes upregulated in comparison ",comparison,sep="")}
    
    # if(!is.null(genes_down)) {
    #   results$Motifdown <- as.data.frame(enricher(entrez_down,
    #                                               TERM2GENE=motifs,
    #                                               universe = universe_Entrez,  
    #                                               pAdjustMethod = pCorrection,
    #                                               pvalueCutoff  = pvalueCutoff,
    #                                               qvalueCutoff = qvalueCutoff))
    #   if(nrow(results$Motifdown)>0){results$Motifdown$Enrichment <- paste("TF binding motif enrichment for genes downregulated in comparison",comparison,sep="")}
    # }
  }
  
  # Immunosignatures enrichment
  if("ImmunoSignatures" %in% GeneSets){
    print("Performing immunesignature enrichment")
    
    results$ImmSigup <- as.data.frame(enricher(entrez_up,
                                               TERM2GENE=immuno_genes,
                                               universe = universe_Entrez,  
                                               pAdjustMethod = pCorrection,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff = qvalueCutoff))
    if(nrow(results$ImmSigup)>0){results$ImmSigup$Enrichment <- paste("Immunosignature enrichment for genes upregulated in comparison ",comparison,sep="")}
    
    # if(!is.null(genes_down)) {
    #   results$ImmSigdown <- as.data.frame(enricher(entrez_down,
    #                                                TERM2GENE=immuno_genes,
    #                                                universe = universe_Entrez,  
    #                                                pAdjustMethod = pCorrection,
    #                                                pvalueCutoff  = pvalueCutoff,
    #                                                qvalueCutoff = qvalueCutoff))
    #   if(nrow(results$ImmSigdown)>0){results$ImmSigdown$Enrichment <- paste("Immunosignature enrichment for genes downregulated in comparison ",comparison,sep="")}
    # }
  }
  return(results)
}


############################################### GSEA dotplot ##############################################


dotplotGSEA <- function(x,
                        show=25,
                        font.size=10,
                        title.size=10,
                        title.width=100,
                        order="count"){
  if(nrow(x)<1){
    print("No enrichment found.")
  }else{
    x <- if(nrow(x)>show){x[c(1:show),]}else{x}
    if(order=="padj"){
      x <- x[order(x$Count,decreasing=FALSE),] 
      x <- x[order(x$p.adjust,decreasing=TRUE),]
      x$Description <- ifelse(nchar(x$Description)<500,
                              paste(substr(x$Description, 1, 500),"[...]",sep=""),
                              x$Description)
      x$Description <- factor(x$Description, levels = unique(x$Description))
    }
    if(order=="count"){
      x <- x[order(x$Count,decreasing=FALSE),]
      x$Description <- ifelse(nchar(x$Description)>500,
                              paste(substr(x$Description, 1, 500),"[...]",sep=""),
                              x$Description)
      x$Description <- factor(x$Description, levels = unique(x$Description))
      x$GeneRatio <- factor(x$GeneRatio, levels = unique(x$GeneRatio))
    }
    
    ggplot(x, aes(x = GeneRatio, y = Description, color = p.adjust)) +
      geom_point(aes(size = Count)) +
      scale_colour_gradientn(colours=c('red', 
                                       'orange', 
                                       'darkblue',
                                       'darkblue'),
                             limits=c(0,1),
                             values   = c(0,0.05,0.2,0.5,1),
                             breaks   = c(0.05,0.2,1),
                             labels = format(c(0.05,0.2,1))) +
      ylab(NULL) +
      ggtitle(paste(strwrap(unique(x$Enrichment), width=title.width), collapse = "\n"))+
      theme_classic() +
      scale_y_discrete(position = "right") +
      theme(text = element_text(size=font.size),
            plot.title = element_text(size=title.size)) 
  }
}


