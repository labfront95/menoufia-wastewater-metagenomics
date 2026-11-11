# Menoufia Wastewater Metagenomics, Antimicrobial Resistance, and Ecological Network Analysis

# This script consolidates the analyses supplied for the completed study.
# Input files are expected in the relative directory: data/
# Output tables are written to the current working directory

# Load Required Libraries
# R packages required to run this code. 
# Refer SessionInfo file to check package version information
library(tidyverse)
library(dplyr)
library(tibble)
library(janitor)
library(reshape2)
library(RColorBrewer)
library(paletteer)
library(ggpubr)
library(ggalt)
library(vegan)
library(ecodist)
library(igraph)
library(ggraph)
library(ggforce)
# Custom Visualization and Analysis Functions
# Generate top-feature plots for each city
bar_plot_top10 <- function(dataframe) {

  top_all <- dataframe %>%
    group_by(city_state, id) %>%
    summarize(sum_normalized_count = sum(normalized_count), .groups = "drop") %>%
    group_by(id) %>%
    mutate(sum_sum_normalized_count = sum(sum_normalized_count)) %>%
    ungroup() %>%
    arrange(desc(sum_sum_normalized_count)) %>%
    mutate(id = factor(id, levels = unique(id))) %>%
    group_by(city_state) %>%
    top_n(10, wt = sum_normalized_count) %>%
    ungroup() %>%
    droplevels()

  top_all_20 <- dataframe %>%
    group_by(city_state, id) %>%
    summarize(sum_normalized_count = sum(normalized_count), .groups = "drop") %>%
    group_by(city_state) %>%
    top_n(20, wt = sum_normalized_count) %>%
    ungroup() %>%
    dplyr::select(city_state, id) %>%
    mutate(present = "yes") %>%
    filter(id %in% top_all$id)

  top_all_color <- setNames(
    paletteer_d("ggsci::default_igv")[1:length(levels(top_all$id))],
    levels(top_all$id)
  )

  table_plot <- top_all %>%
    ggplot(aes(y = id, x = city_state)) +
    geom_point(size = 5, shape = 21, aes(fill = id), stroke = 0.75) +
    geom_point(data = top_all_20, size = 5, shape = 21,
               aes(color = present), stroke = 0.75) +
    theme_bw(25) +
    scale_fill_manual(values = top_all_color, name = "") +
    scale_color_manual(values = c("black", "white"), name = "") +
    scale_x_discrete(
      position = "top",
      labels = c("Delhi" = "D", "Kolkata" = "K",
                 "Chennai" = "C", "Mumbai" = "M")
    ) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    ) +
    scale_y_discrete(
      limits = rev,
      position = "right",
      label = function(x) str_trunc(x, width = 25)
    ) +
    xlab("") +
    ylab("")

  city_top <- list()

  for (city in levels(as.factor(dataframe$city_state))) {

    top <- top_all %>%
      filter(city_state == city) %>%
      arrange(desc(sum_normalized_count)) %>%
      mutate(id = factor(id, levels = id))

    city_top[[city]] <- dataframe %>%
      filter(city_state == city) %>%
      group_by(month_label) %>%
      mutate(fraction_normalized_count = normalized_count / sum(normalized_count)) %>%
      ungroup() %>%
      filter(id %in% top$id) %>%
      mutate(id = factor(id, levels = levels(top$id))) %>%
      ggplot(aes(x = month_label, y = fraction_normalized_count,
                 fill = fct_inorder(id), group = id)) +
      geom_bar(position = "stack", stat = "identity") +
      scale_fill_manual(values = top_all_color, name = "") +
      scale_x_discrete(drop = FALSE) +
      theme_bw(25) +
      labs(x = "", y = "", title = city) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.3, size = 20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")
      ) +
      guides(fill = guide_legend(ncol = 1))
  }

  return(list(table_plot, city_top))
}

#  Alpha-diversity curves across sampling months
alpha_diversity_plot <- function(dataframe) {

  dataframe %>%
    group_by(month_label, city_state) %>%
    summarise(
      richness = n_distinct(id[normalized_count > 0]),
      shannon = diversity(normalized_count, index = "shannon"),
      simpson = diversity(normalized_count, index = "simpson"),
      .groups = "drop"
    ) %>%
    melt(id.vars = c("month_label", "city_state")) %>%
    ggplot(aes(x = as.factor(month_label), y = value,
               fill = city_state, group = city_state)) +
    geom_point(size = 5, shape = 21, show.legend = TRUE, alpha = 1) +
    geom_xspline(size = 1, aes(color = city_state),
                 show.legend = TRUE, alpha = 0.7) +
    scale_color_brewer(palette = "Dark2", name = "City") +
    scale_fill_brewer(palette = "Dark2", name = "City") +
    theme_minimal(25) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    labs(x = "Collection month", y = "", title = "") +
    facet_wrap(~variable, scales = "free_y", nrow = 3)
}

# PCoA plot based on Bray-Curtis dissimilarity
bray_curtis_plot <- function(dataframe) {

  dataframe_wider <- dataframe %>%
    pivot_wider(
      id_cols = id,
      names_from = fileName,
      values_from = normalized_count,
      values_fill = 0
    ) %>%
    column_to_rownames(var = "id")

  bray_curtis_dist <- vegdist(t(dataframe_wider), method = "bray")
  bray_curtis_pcoa <- pco(bray_curtis_dist)
  eigenvalues <- bray_curtis_pcoa$values
  percentVar <- (eigenvalues / sum(eigenvalues)) * 100

  bray_curtis_tibble <- tibble(
    pcoa1 = bray_curtis_pcoa$vectors[, 1],
    pcoa2 = bray_curtis_pcoa$vectors[, 2]
  ) %>%
    mutate(fileName = rownames(bray_curtis_pcoa$vectors)) %>%
    left_join(dataframe, by = "fileName")

  bray_curtis_tibble_centroid <- bray_curtis_tibble %>%
    group_by(city_state) %>%
    summarize(pcoa1 = mean(pcoa1), pcoa2 = mean(pcoa2), .groups = "drop")

  bray_curtis_tibble %>%
    left_join(
      dataframe %>%
        group_by(month_label, city_state) %>%
        summarise(
          shannon = diversity(normalized_count, index = "shannon"),
          .groups = "drop"
        ),
      by = c("month_label", "city_state")
    ) %>%
    ggplot(aes(x = pcoa1, y = pcoa2)) +
    geom_point(shape = 21, aes(fill = as.factor(city_state), size = shannon)) +
    geom_point(
      data = bray_curtis_tibble_centroid,
      aes(x = pcoa1, y = pcoa2, fill = city_state),
      size = 15, color = "black", shape = 22, stroke = 2, alpha = 0.75
    ) +
    scale_size_continuous(range = c(2, 10), transform = "reverse") +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    xlab(paste0("PCo1: (", round(percentVar[1], 2), "%)")) +
    ylab(paste0("PCo2: (", round(percentVar[2], 2), "%)")) +
    theme_bw(30) +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 7))) +
    theme(
      axis.text.x = element_text(angle = 0, vjust = 0.3),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(fill = "City", size = "Shannon") +
    stat_ellipse(aes(color = as.factor(city_state)),
                 alpha = 1, show.legend = FALSE)
}


# Load Datasets and Prepare Metadata

meta <- read_delim(file = "data/metadata.tsv", delim = "\t", col_names = TRUE) %>%
  clean_names() %>%
  arrange(year, month) %>%
  mutate(month_label = factor(month_label, levels = unique(month_label))) %>%
  mutate(city_state = factor(
    city_state,
    levels = c("Kolkata", "Delhi", "Chennai", "Mumbai",
               "Ashmoon", "El-Bagoor", "Menoof", "Shibin El-Kom")
  ))

k2_nt <- read_delim(
  file = "data/k2_nt_20230502.inspect.species.krona.phylotree.txt",
  delim = "\t"
)
bracken_phylum <- read_delim("data/bracken_phylum.tsv", delim = "\t", col_names = TRUE)
bracken_genus <- read_delim("data/bracken_genus.tsv", delim = "\t", col_names = TRUE)
bracken_species <- read_delim("data/bracken_species.tsv", delim = "\t", col_names = TRUE)
bacterial_read_df <- read_delim("data/bacterial_read_df.tsv", delim = "\t", col_names = TRUE)
rgi_bowtie <- read_delim("data/rgi_bowtie.tsv", delim = "\t", col_names = TRUE)
checkm2 <- read_delim("data/checkm2_quality_report.csv", delim = "\t")
coverage <- read_delim(
  file = "data/sample_to_dereplicated_bin_coverage.tsv",
  delim = "\t", col_names = TRUE
)
arg_matrix <- read_delim("data/ARG_counts_427_aligned.tsv", delim = "\t", col_names = TRUE)


#  ARG Audit and Sequencing-Depth Correlation

arg_features <- setdiff(colnames(arg_matrix), "sample_code")

arg_audit <- arg_matrix %>%
  column_to_rownames("sample_code") %>%
  mutate(
    ARG_richness = rowSums(. > 0),
    ARG_read_sum = rowSums(.)
  ) %>%
  rownames_to_column("sample_code") %>%
  left_join(meta, by = "sample_code")

# Spearman correlation between observed ARG richness and ARG read counts
cor_depth <- cor.test(
  arg_audit$ARG_richness,
  arg_audit$ARG_read_sum,
  method = "spearman"
)

cat(
  "Spearman Correlation (Richness vs ARG Reads): rho =",
  cor_depth$estimate,
  ", p =", cor_depth$p.value, "\n"
)


#  ARG Rarefaction Sensitivity Analysis
#    1,000 repeated multinomial draws at a target depth of 8,000 reads

TARGET_DEPTH <- 8000
N_REPS <- 1000

eligible_samples <- arg_audit %>%
  filter(ARG_read_sum >= TARGET_DEPTH)

cat(
  "Samples Retained for Rarefaction:",
  nrow(eligible_samples), "/", nrow(arg_audit), "\n"
)

rarefaction_res <- list()
set.seed(20261013)

for (i in seq_len(nrow(eligible_samples))) {

  s_code <- eligible_samples$sample_code[i]
  counts <- as.numeric(
    arg_matrix %>%
      filter(sample_code == s_code) %>%
      dplyr::select(-sample_code)
  )

  # Repeated sampling using the multinomial procedure in the supplied analysis
  sampled_richness <- replicate(N_REPS, {
    draw <- rmultinom(1, size = TARGET_DEPTH, prob = counts)
    sum(draw > 0)
  })

  rarefaction_res[[i]] <- tibble(
    sample_code = s_code,
    city_state = eligible_samples$city_state[i],
    observed_richness = sum(counts > 0),
    rarefied_richness_mean = mean(sampled_richness)
  )
}

rarefaction_df <- bind_rows(rarefaction_res)

# Kruskal-Wallis test across cities for rarefied ARG richness
kw_test <- kruskal.test(
  rarefied_richness_mean ~ city_state,
  data = rarefaction_df
)

cat(
  "Kruskal-Wallis Test across Cities (Rarefied): H =",
  kw_test$statistic,
  ", p =", kw_test$p.value, "\n"
)


#  Read-Level Metagenomics and Taxonomic Profiling

bacterial_read_df_select <- bacterial_read_df %>%
  left_join(meta, by = c("fileName" = "sample_code")) %>%
  group_by(city_state, month_label) %>%
  summarise(
    total_reads = sum(total_reads),
    bacterial_reads = sum(bacterial_reads),
    .groups = "drop"
  ) %>%
  mutate(fileName = paste0(city_state, "_", month_label))

bracken_genus_select_bacterial <- bracken_genus %>%
  left_join(meta, by = c("fileName" = "sample_code")) %>%
  group_by(city_state, month_label, name) %>%
  summarise(new_est_reads = sum(new_est_reads), .groups = "drop") %>%
  mutate(fileName = paste0(city_state, "_", month_label)) %>%
  inner_join(
    k2_nt %>%
      filter(kingdom == "Bacteria") %>%
      dplyr::select(genus) %>%
      distinct(),
    by = c("name" = "genus")
  ) %>%
  left_join(
    bacterial_read_df_select,
    by = c("city_state", "month_label", "fileName")
  ) %>%
  mutate(normalized_count = (new_est_reads / bacterial_reads) * 1e06) %>%
  dplyr::rename(id = name)

# Genus-level alpha diversity and Bray-Curtis PCoA
alpha_diversity_plot(bracken_genus_select_bacterial)
bray_curtis_plot(bracken_genus_select_bacterial)


#  Co-occurrence Network Analysis and Community Detection

coverage_bin_wide <- coverage %>%
  left_join(meta, by = c("sampleName" = "sample_code")) %>%
  group_by(bin, city_state, month_label) %>%
  summarise(
    binLen = sum(contigLen),
    nucleotide_mapped = sum(nucleotide_mapped),
    .groups = "drop"
  ) %>%
  mutate(avg_cov = nucleotide_mapped / binLen) %>%
  mutate(sampleName = paste0(city_state, "_", month_label)) %>%
  pivot_wider(
    names_from = sampleName,
    id_cols = bin,
    values_from = avg_cov,
    values_fill = 0
  ) %>%
  column_to_rownames("bin")

# Retain MAGs with mean coverage greater than 0.5
coverage_filtered <- coverage_bin_wide[rowMeans(coverage_bin_wide) > 0.5, ]

target_cities <- c("Ashmoon", "El-Bagoor", "Menoof", "Shibin El-Kom")
final_network_objects <- list()

for (city in target_cities) {

  city_cols <- grep(paste0("^", city), colnames(coverage_filtered), value = TRUE)
  sub_mat <- as.matrix(coverage_filtered[, city_cols])

  # Calculate Pearson correlations among MAG coverage profiles
  cor_mat <- cor(t(sub_mat), method = "pearson")

  # Calculate correlation p-values and apply Benjamini-Hochberg adjustment
  n_obs <- ncol(sub_mat)
  t_stat <- cor_mat * sqrt((n_obs - 2) / (1 - cor_mat^2))
  p_mat <- 2 * (1 - pt(abs(t_stat), df = n_obs - 2))
  diag(p_mat) <- 1

  adj_p_mat <- matrix(
    p.adjust(p_mat, method = "BH"),
    nrow = nrow(p_mat)
  )
  rownames(adj_p_mat) <- colnames(adj_p_mat) <- rownames(cor_mat)

  # Build adjacency matrix: Pearson r >= 0.70 and FDR < 0.05
  adj_matrix <- (cor_mat >= 0.70) & (adj_p_mat < 0.05)
  diag(adj_matrix) <- FALSE

  g <- graph_from_adjacency_matrix(
    adj_matrix * 1,
    mode = "undirected",
    diag = FALSE
  )
  g <- delete_vertices(g, which(degree(g) == 0))

  # Community detection using the Louvain algorithm
  set.seed(20261102)
  comm <- cluster_louvain(g)
  obs_Q <- modularity(comm)

  # Null-model permutations for modularity validation (1,000 repetitions)
  null_Q <- replicate(N_REPS, {
    g_null <- rewiring.edges(g, prob = 0.8)
    modularity(cluster_louvain(g_null))
  })

  z_score <- (obs_Q - mean(null_Q)) / sd(null_Q)
  p_val <- (1 + sum(null_Q >= obs_Q)) / (N_REPS + 1)

  metrics <- tibble(
    City = city,
    Nodes = vcount(g),
    Edges = ecount(g),
    Edge_Density = edge_density(g),
    Mean_Degree = mean(degree(g)),
    Observed_Q = obs_Q,
    Null_Mean_Q = mean(null_Q),
    Null_SD_Q = sd(null_Q),
    Modularity_Z = z_score,
    Permutation_P = p_val,
    Number_of_Communities = length(unique(membership(comm)))
  )

  final_network_objects[[city]] <- list(
    graph = g,
    membership = membership(comm),
    metrics = metrics,
    observed_Q = obs_Q,
    null_modularities = null_Q
  )
}

# Combine and export network topology summary (Supplementary Table S9)
table_S9 <- bind_rows(lapply(final_network_objects, function(x) x$metrics))
print(table_S9)
write_csv(table_S9, "Supplementary_Table_S9_Network_Topology.csv")


#  Correlation-Threshold Sensitivity Analysis (r = 0.70, 0.80, 0.90)

r_thresholds <- c(0.70, 0.80, 0.90)
sensitivity_list <- list()

for (city in target_cities) {

  city_cols <- grep(paste0("^", city), colnames(coverage_filtered), value = TRUE)
  sub_mat <- as.matrix(coverage_filtered[, city_cols])
  cor_mat <- cor(t(sub_mat), method = "pearson")

  n_obs <- ncol(sub_mat)
  t_stat <- cor_mat * sqrt((n_obs - 2) / (1 - cor_mat^2))
  p_mat <- 2 * (1 - pt(abs(t_stat), df = n_obs - 2))
  diag(p_mat) <- 1

  adj_p_mat <- matrix(
    p.adjust(p_mat, method = "BH"),
    nrow = nrow(p_mat)
  )

  for (rt in r_thresholds) {

    adj_matrix <- (cor_mat >= rt) & (adj_p_mat < 0.05)
    diag(adj_matrix) <- FALSE

    g <- graph_from_adjacency_matrix(
      adj_matrix * 1,
      mode = "undirected",
      diag = FALSE
    )
    g <- delete_vertices(g, which(degree(g) == 0))

    set.seed(20260902)
    comm <- cluster_louvain(g)

    sensitivity_list[[length(sensitivity_list) + 1]] <- tibble(
      City = city,
      R_threshold = rt,
      Active_Nodes = vcount(g),
      Edges = ecount(g),
      Density = edge_density(g),
      Modularity = modularity(comm)
    )
  }
}

sensitivity_df <- bind_rows(sensitivity_list)
print(sensitivity_df)
write_csv(sensitivity_df, "Network_Correlation_Threshold_Sensitivity.csv")



