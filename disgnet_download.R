# Load necessary libraries
library(disgenet2r)
library(dplyr)
library(tidyr)
library(ggplot2)

# Set API key
api_key <- "cdbfc97d-eff7-4f56-83fe-b558a7e34c25"
Sys.setenv(DISGENET_API_KEY = api_key)

# Create a named vector of psychiatric disorders with their UMLS CUIs
psychiatric_disorders <- c(
  "UMLS_C0036341",  # Schizophrenia
  "UMLS_C0005586",  # Bipolar Disorder
  "UMLS_C0011570",  # Major Depressive Disorder
  "UMLS_C0003469",  # Anxiety Disorders
  "UMLS_C0004352",  # Autism Spectrum Disorder
  "UMLS_C1263846",  # Attention Deficit Disorder
  "UMLS_C0038436",  # PTSD
  "UMLS_C0028768"   # Obsessive-Compulsive Disorder
)

# Query GDAs for psychiatric disorders
results <- disease2gene(
  disease = psychiatric_disorders,
  database = "CURATED",
  score = c(0.3, 1),
  verbose = TRUE
)

if (!is.null(results) && inherits(results, "DataGeNET.DGN")) {
  # Convert list columns to character
  process_list_column <- function(x) {
    if(is.list(x)) {
      sapply(x, function(y) paste(unlist(y), collapse = "; "))
    } else {
      x
    }
  }
  
  # Convert results to a data frame with comprehensive columns
  gda_data <- as_tibble(results@qresult) %>%
    mutate(across(everything(), process_list_column)) %>%
    select(
      gene_symbol,
      geneid,
      ensemblid,
      geneNcbiType,
      uniprotids,
      protein_class_name,
      disease_name,
      diseaseType,
      diseaseUMLSCUI,
      score,
      yearInitial,
      yearFinal,
      numPMIDs,
      evidence_index
    ) %>%
    distinct()
  
  # Create comprehensive gene summary statistics
  gene_summary <- gda_data %>%
    group_by(gene_symbol) %>%
    summarise(
      n_diseases = n_distinct(disease_name),
      mean_score = mean(score),
      max_score = max(score),
      ensembl_id = first(ensemblid),
      entrez_id = first(geneid),
      uniprot_ids = first(uniprotids),
      protein_classes = first(protein_class_name),
      diseases = paste(unique(disease_name), collapse = "; "),
      total_pmids = sum(numPMIDs),
      mean_evidence = mean(evidence_index, na.rm = TRUE),
      years_studied = paste(min(yearInitial), max(yearFinal), sep = "-"),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_score))
  
  # Create comprehensive disease summary statistics
  disease_summary <- gda_data %>%
    group_by(disease_name, diseaseType, diseaseUMLSCUI) %>%
    summarise(
      n_genes = n_distinct(gene_symbol),
      mean_score = mean(score),
      total_pmids = sum(numPMIDs),
      mean_evidence = mean(evidence_index, na.rm = TRUE),
      top_genes = paste(gene_symbol[order(-score)][1:min(5, n())], collapse = "; "),
      years_studied = paste(min(yearInitial), max(yearFinal), sep = "-"),
      .groups = "drop"
    ) %>%
    arrange(desc(n_genes))
  
  # Create STRING-compatible format
  string_compatible <- gda_data %>%
    select(gene_symbol, geneid, ensemblid, uniprotids, score) %>%
    filter(!is.na(uniprotids) & uniprotids != "") %>%
    distinct()
  
  # Save all processed datasets
  write.csv(gda_data, "psychiatric_disorders_full_GDA.csv", row.names = FALSE)
  write.csv(gene_summary, "psychiatric_disorders_gene_summary.csv", row.names = FALSE)
  write.csv(disease_summary, "psychiatric_disorders_disease_summary.csv", row.names = FALSE)
  write.csv(string_compatible, "string_compatible_genes.csv", row.names = FALSE)
  
  # Create visualization of top genes
  top_genes_plot <- ggplot(
    head(gene_summary, 20), 
    aes(x = reorder(gene_symbol, mean_score), y = mean_score)
  ) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      title = element_text(size = 12, face = "bold")
    ) +
    labs(
      title = "Top 20 Genes Associated with Psychiatric Disorders",
      x = "Gene Symbol",
      y = "Mean Association Score",
      caption = paste("Based on", nrow(gda_data), "gene-disease associations")
    )
  
  # Save plot
  ggsave("top_genes_plot.png", top_genes_plot, width = 12, height = 8)
  
  # Create high-confidence genes list with comprehensive info
  high_confidence_genes <- gene_summary %>%
    filter(mean_score > 0.5) %>%
    select(gene_symbol, mean_score, ensembl_id, entrez_id, uniprot_ids, 
           n_diseases, total_pmids, mean_evidence)
  
  # Save high-confidence genes
  write.csv(high_confidence_genes, "high_confidence_genes.csv", row.names = FALSE)
  
  # Print detailed summary statistics
  cat("\nDetailed Summary Statistics:\n")
  cat("--------------------------------\n")
  cat("Total number of genes:", nrow(gene_summary), "\n")
  cat("High-confidence genes (score > 0.5):", sum(gene_summary$mean_score > 0.5), "\n")
  cat("Number of diseases analyzed:", n_distinct(gda_data$disease_name), "\n")
  cat("Total PMIDs supporting associations:", sum(gda_data$numPMIDs), "\n")
  cat("Year range of studies:", min(gda_data$yearInitial), "-", max(gda_data$yearFinal), "\n")
  
  # Print top genes per disease
  cat("\nTop Genes by Disease:\n")
  cat("--------------------------------\n")
  for(disease in unique(disease_summary$disease_name)) {
    cat("\n", disease, ":\n")
    top_genes <- gda_data %>%
      filter(disease_name == disease) %>%
      arrange(desc(score)) %>%
      slice_head(n = 5) %>%
      select(gene_symbol, score)
    print(top_genes)
  }
  
  cat("\nFiles created:\n")
  cat("1. psychiatric_disorders_full_GDA.csv - Complete gene-disease associations\n")
  cat("2. psychiatric_disorders_gene_summary.csv - Gene-level statistics\n")
  cat("3. psychiatric_disorders_disease_summary.csv - Disease-level statistics\n")
  cat("4. string_compatible_genes.csv - Format ready for STRING database\n")
  cat("5. high_confidence_genes.csv - High confidence genes (score > 0.5)\n")
  cat("6. top_genes_plot.png - Visualization of top 20 genes\n")
  
} else {
  stop("No results were found for the query. Please check the identifiers or parameters.")
}

cat("\nAnalysis complete! Check the created files for detailed results.\n")