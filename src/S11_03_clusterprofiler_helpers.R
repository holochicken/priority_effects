#!/usr/bin/env Rscript

clusterprofiler_helpers <- new.env()

clusterprofiler_helpers$get_up_genes <- function(tidybulk_de_annotated) {
  tidybulk_de_annotated %>%
    tidybulk::filter(significant & direction > 0) %>%
    dplyr::pull(gene)
}

clusterprofiler_helpers$get_down_genes <- function(tidybulk_de_annotated) {
  tidybulk_de_annotated %>%
    tidybulk::filter(significant & direction < 0) %>%
    dplyr::pull(gene)
}


clusterprofiler_helpers$compute_go_terms <-
  function(tidybulk_de_annotated, up_or_down = "up", organism = "org.Gg.eg.db") {
    de_genes <-
      dplyr::if_else(
        up_or_down == "up",
        tidybulk_de_annotated %>%
          get_up_genes() %>%
          list(),
        tidybulk_de_annotated %>%
          get_down_genes() %>%
          list()
      ) %>%
      unlist()
    
    all_genes <- tidybulk_de_annotated %>% dplyr::pull(gene)
    
    clusterProfiler::enrichGO(
      minGSSize = 10, maxGSSize = 500,
      gene = de_genes,
      OrgDb = organism,
      universe = all_genes,
      keyType = "SYMBOL",
      ont = "ALL",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = T,
      pool = T
    )
  }



clusterprofiler_helpers$get_pathway_to_symbol <-
  function(universe,
           organism = "org.Gg.eg.db",
           kegg_organism = "gga",
           keggType = "KEGG") {
    pathway_to_entrezid <-
      clusterProfiler::download_KEGG(
        species = kegg_organism,
        keggType = keggType
      ) %>%
      .[["KEGGPATHID2EXTID"]] %>%
      tibble::as_tibble() %>%
      dplyr::rename(pathway = from, ENTREZID = to)
    
    symbol_to_entrezid <-
      universe %>%
      clusterProfiler::bitr(from = "SYMBOL", to = "ENTREZID", OrgDb = organism)
    
    pathway_to_entrezid %>%
      dplyr::left_join(symbol_to_entrezid) %>%
      dplyr::select(pathway, gene = SYMBOL) %>%
      dplyr::filter(!is.na(gene))
  }


clusterprofiler_helpers$get_pathway_to_name <- function(kegg_organism = "gga", keggType = "KEGG") {
  clusterProfiler::download_KEGG(kegg_organism, keggType) %>%
    .[["KEGGPATHID2NAME"]] %>%
    tibble::as_tibble() %>%
    dplyr::rename(pathway = from, name = to) %>%
    dplyr::distinct()
}


clusterprofiler_helpers$compute_kegg_terms <- function(tidybulk_de_annotated,
                                                       up_or_down = "up",
                                                       organism = "org.Gg.eg.db",
                                                       kegg_organism = "gga",
                                                       keggType = "KEGG") {
  de_genes <-
    if_else(
      up_or_down == "up",
      tidybulk_de_annotated %>% get_up_genes() %>% list(),
      tidybulk_de_annotated %>% get_down_genes() %>% list()
    ) %>%
    unlist()
  
  universe <- tidybulk_de_annotated %>% pull(gene)
  
  term2gene <- get_pathway_to_symbol(universe, organism, kegg_organism, keggType)
  term2name <- get_pathway_to_name(kegg_organism, keggType)
  
  clusterProfiler::enricher(
    gene = de_genes,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    universe = universe,
    minGSSize = 10,
    maxGSSize = 500,
    qvalueCutoff = 0.2,
    TERM2GENE = term2gene,
    TERM2NAME = term2name
  )
}


attach(clusterprofiler_helpers, name = "clusterprofiler_helpers")
rm(clusterprofiler_helpers)