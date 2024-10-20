library(dplyr)
library(ggplot2)
library(reshape2)
library(scales)
library(vegan)
library(ggrepel)
library(zCompositions)
library(compositions)
library(ANCOMBC)
library(phyloseq)
library(TreeSummarizedExperiment)
library(SummarizedExperiment)
library(pairwiseAdonis)
library(ape)
library(stringr)

# User input for file names
input_file <- "mealworm_genus.csv" 
output_prefix <- "mealworm"  
metadata_file <- "metadata_meal.csv"
#Checking for unwanted charachters
check_unwanted_characters <- function(df, file_name) {
  if (any(str_detect(names(df), "^\\s+|\\s+$"))) {
    warning(paste("Warning: Leading or trailing spaces found in column names of", file_name))
  }
  problematic_columns <- sapply(df, function(col) {
    any(str_detect(as.character(col), "^\\s+|\\s+$")) || any(str_detect(as.character(col), "[^[:alnum:][:space:].-_]"))
  })
  if (any(problematic_columns)) {
    cat(paste("Columns with unwanted characters or leading/trailing spaces found in", file_name, ":\n"))
    cat(names(df)[problematic_columns], "\n")
  } else {
    cat(paste("No issues found in", file_name, "\n"))
  }
}
abundance_data <- read.csv(input_file, row.names = 1, check.names = FALSE)
metadata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
cat("Checking input file for issues...\n")
check_unwanted_characters(abundance_data, input_file)
cat("Checking metadata file for issues...\n")
check_unwanted_characters(metadata, metadata_file)

# User input for columns of interest; select any one before proceeding
# If you want to analyze only specific columns, define them as shown below
# Otherwise, set columns_of_interest to NULL for analyzing the entire dataset
columns_of_interest <- NULL
columns_of_interest <-  c("Standard_1a", "Standard_2a", "Standard_2b", "Standard_3a") 
if (!is.null(columns_of_interest)) {
  abundance_data <- abundance_data[, columns_of_interest, drop = FALSE]
}

total_abundance <- rowSums(abundance_data)
threshold_percentage <- 0.005
threshold_value <- sum(total_abundance) * threshold_percentage
abundant_mags <- total_abundance[total_abundance >= threshold_value]
abundant_mag_names <- names(abundant_mags)
filtered_abundance <- abundance_data[abundant_mag_names, , drop = FALSE]
filtered_abundance["Others", ] <- colSums(abundance_data[!rownames(abundance_data) %in% abundant_mag_names, , drop = FALSE])
write.csv(filtered_abundance, paste0(output_prefix, "_Filtered_Abundance_genus.csv"), row.names = TRUE)

#Stacked bars
filtered_abundance <- read.csv(paste0(output_prefix, "_Filtered_Abundance_genus.csv"), row.names = 1)
abundance_long <- melt(filtered_abundance, variable.name = "Sample", value.name = "Abundance")
abundance_long$MAG <- rep(rownames(filtered_abundance), times = ncol(filtered_abundance))
unique_mags <- unique(abundance_long$MAG)
color_palette <- hue_pal()(length(unique_mags) - 1)  
color_palette <- c(color_palette, "black")  
abundance_long$MAG <- factor(abundance_long$MAG, levels = unique_mags)
stacked_bar_plot <- ggplot(abundance_long, aes(x = Sample, y = Abundance, fill = MAG)) +
  geom_bar(stat = "identity", position = "stack", color = "black", size = 0.2) +  
  scale_fill_manual(values = setNames(color_palette, unique_mags)) +
  labs(x = "Sample", y = "Abundance", fill = "Genera") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, size = 12, hjust = 1),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 10)
  )
stacked_bar_file <- paste0(output_prefix, "_stacked_bar_plot.png")
ggsave(stacked_bar_file, plot = stacked_bar_plot, width = 10, height = 8)

# Bubble plot
bubble_plot <- ggplot(abundance_long, aes(x = Sample, y = MAG, size = Abundance, fill = MAG)) +
  geom_point(shape = 21, stroke = 1, color = "black", alpha = 0.7) +  
  scale_size_continuous(range = c(3, 15), breaks = c(0.1, 1, 3, 10, 30)) +
  scale_fill_manual(values = color_palette) +
  labs(x = "Genera", y = "Sample", size = "Abundance") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, size = 14, hjust = 1),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 8),
    legend.text = element_text(size = 10)
  ) +
  guides(fill = FALSE)
bubble_plot_file <- paste0(output_prefix, "_bubble_plot.png")
ggsave(bubble_plot_file, plot = bubble_plot, width = 10, height = 8)

#alpha Diversity metrics
abundance_table <- read.csv(input_file, row.names = 1)
#user defined samples of interest
columns_of_interest <- NULL
columns_of_interest <-  c("Standard_1a", "Standard_2a", "Standard_2b", "Standard_3a") 
if (!is.null(columns_of_interest)) {
  abundance_table <- abundance_table[, columns_of_interest, drop = FALSE]
}

sample_types <- gsub("_[0-9]+[a-z]$", "", colnames(abundance_table))
results <- data.frame(
  SampleType = sample_types,
  Shannon = numeric(ncol(abundance_table)),
  Simpson = numeric(ncol(abundance_table)),
  InvSimpson = numeric(ncol(abundance_table)),
  Richness = numeric(ncol(abundance_table)),
  Evenness = numeric(ncol(abundance_table))
)
for (i in 1:ncol(abundance_table)) {
  sample_data <- abundance_table[, i]
  shannon <- diversity(sample_data, index = "shannon")
  simpson <- diversity(sample_data, index = "simpson")
  invsimpson <- diversity(sample_data, index = "invsimpson")
  richness <- specnumber(sample_data)
  evenness <- ifelse(richness > 0, shannon / log(richness), NA)
  
  results[i, "Shannon"] <- shannon
  results[i, "Simpson"] <- simpson
  results[i, "InvSimpson"] <- invsimpson
  results[i, "Richness"] <- richness
  results[i, "Evenness"] <- evenness
}
alpha_diversity_output_file <- paste0(output_prefix, "_alpha_diversity_metrics.csv")
write.csv(results, alpha_diversity_output_file, row.names = FALSE)
results_melt <- melt(results, id.vars = "SampleType")
stats <- results_melt %>%
  group_by(SampleType, variable) %>%
  summarise(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )
alpha_div_plot <- ggplot(results_melt, aes(x = SampleType, y = value, fill = SampleType)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 2, color = "black", position = position_dodge(width = 0.8)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, position = position_dodge(width = 0.8)) +
  labs(title = "Alpha Diversity Indices with Standard Deviation",
       x = "Sample Type", y = "Value") +
  theme_minimal() +
  facet_wrap(~ variable, scales = "free_y") +
  theme(legend.title = element_blank())
alpha_diversity_plot_file <- paste0(output_prefix, "_alpha_diversity_plot.png")
ggsave(alpha_diversity_plot_file, plot = alpha_div_plot, width = 10, height = 8)

#PCoA
mag_abundances <- read.csv(input_file, row.names = 1)
#user defined samples of interest
columns_of_interest <- NULL  
columns_of_interest <-  c("Standard_1a", "Standard_2a", "Standard_2b", "Standard_3a")
if (!is.null(columns_of_interest)) {
  mag_abundances <- mag_abundances[, columns_of_interest, drop = FALSE]
}
mag_abundances_t <- t(mag_abundances)
sample_types <- gsub("_[0-9]+[a-z]$", "", rownames(mag_abundances_t))

bray_curtis_dist <- vegdist(mag_abundances_t, method = "bray")
pcoa_result <- pcoa(bray_curtis_dist)
pcoa_df <- data.frame(Sample = rownames(pcoa_result$vectors),
                      Axis1 = pcoa_result$vectors[, 1],
                      Axis2 = pcoa_result$vectors[, 2],
                      Site = sample_types)

plot_untransformed <- ggplot(pcoa_df, aes(x = Axis1, y = Axis2, shape = Site, color = Site, fill = Site)) +
  geom_point(size = 3, stroke = 1.5) +
  xlab(paste0("PCoA Axis 1 (", round(pcoa_result$values$Relative_eig[1] * 100, 2), "%)")) +
  ylab(paste0("PCoA Axis 2 (", round(pcoa_result$values$Relative_eig[2] * 100, 2), "%)")) +
  ggtitle("PCoA (Bray-Curtis, Untransformed)") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  coord_fixed() +
  scale_shape_manual(values = 21:25) +
  scale_color_manual(values = rainbow(length(unique(sample_types)))) +
  scale_fill_manual(values = rainbow(length(unique(sample_types)))) +
  guides(shape = guide_legend(title = "Site"), color = guide_legend(title = "Site"), fill = guide_legend(title = "Site"))
print(plot_untransformed)
pcoa_untransformed_file <- paste0(output_prefix, "_pcoa_untransformed.png")
ggsave(pcoa_untransformed_file, plot = plot_untransformed, width = 10, height = 8)
mag_abundances_nonzero <- zCompositions::cmultRepl(mag_abundances_t, method = "CZM")
mag_abundances_clr <- clr(mag_abundances_nonzero)
summary(mag_abundances_clr)
mag_abundances_clr[mag_abundances_clr < 0] <- 0

bray_curtis_dist_clr <- vegdist(mag_abundances_clr, method = "bray")
pcoa_result_clr <- pcoa(bray_curtis_dist_clr)
pcoa_df_clr <- data.frame(Sample = rownames(pcoa_result_clr$vectors),
                          Axis1 = pcoa_result_clr$vectors[, 1],
                          Axis2 = pcoa_result_clr$vectors[, 2],
                          Site = sample_types)

plot_clr <- ggplot(pcoa_df_clr, aes(x = Axis1, y = Axis2, shape = Site, color = Site, fill = Site)) +
  geom_point(size = 3, stroke = 1.5) +
  xlab(paste0("PCoA Axis 1 (", round(pcoa_result_clr$values$Relative_eig[1] * 100, 2), "%)")) +
  ylab(paste0("PCoA Axis 2 (", round(pcoa_result_clr$values$Relative_eig[2] * 100, 2), "%)")) +
  ggtitle("PCoA (Bray-Curtis, CLR Transformed)") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  coord_fixed() +
  scale_shape_manual(values = 21:25) +
  scale_color_manual(values = rainbow(length(unique(sample_types)))) +
  scale_fill_manual(values = rainbow(length(unique(sample_types)))) +
  guides(shape = guide_legend(title = "Site"), color = guide_legend(title = "Site"), fill = guide_legend(title = "Site"))
print(plot_clr)
pcoa_clr_file <- paste0(output_prefix, "_pcoa_clr.png")
ggsave(pcoa_clr_file, plot = plot_clr, width = 10, height = 8)

#PERMANOVA_pairwise 
#This will test the null hypothesis that the spread (centroids) are the same for all groups.
#Betadisper and Tukey's Test is performed to see significant differences in dispersion/mean distance.
metadata <- read.csv(metadata_file)
abundance_data <- read.csv(input_file, row.names = 1)
abundance_cleaned <- abundance_data[rowSums(abundance_data) > 0, ]
abundance_cleaned <- abundance_cleaned[, metadata$SampleID]
metadata$Condition <- as.factor(metadata$Condition)
output_file <- paste0(output_prefix, "_combined_results_genus.txt")
sink(output_file)
cat("=== PERMANOVA and Pairwise PERMANOVA Results ===\n\n")
pairwise_results <- pairwise.adonis2(t(abundance_cleaned) ~ Condition, data = metadata, method = "bray")
print(pairwise_results)
cat("\n\n=== Betadisper and ANOVA Results ===\n\n")
betadisper_result <- betadisper(vegdist(t(abundance_cleaned), method = "bray"), metadata$Condition)
anova_result <- anova(betadisper_result)
print(anova_result)
tukey_result <- TukeyHSD(betadisper_result)
print(tukey_result)
sink()
extract_r2_p_values <- function(file_path) {
  file_content <- readLines(file_path)
  r2_values <- list()
  p_values <- list()
  for (i in seq_along(file_content)) {
    if (grepl("^\\$", file_content[i])) {
      comparison_name <- gsub("^\\$", "", file_content[i])
      for (j in (i + 1):(i + 5)) {
        if (grepl("Model", file_content[j])) {
          r2_value <- as.numeric(strsplit(file_content[j], "\\s+")[[1]][5])
          r2_values[[comparison_name]] <- r2_value
          p_value <- as.numeric(strsplit(file_content[j], "\\s+")[[1]][6])
          p_values[[comparison_name]] <- p_value
        }
      }
    }
  }
  r2_p_df <- data.frame(
    Comparison = names(r2_values),
    R_squared = unlist(r2_values),
    P_value = unlist(p_values)
  )
  return(r2_p_df)
}
r2_p_df <- extract_r2_p_values(output_file)
print(r2_p_df)
r2_p_df$Significant <- ifelse(r2_p_df$P_value < 0.05, "Significant", "Not Significant")
#plot PERMANOVA
permanova_plot <- ggplot(r2_p_df, aes(x = Comparison, y = R_squared, fill = Significant)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = c("Significant" = "yellow", "Not Significant" = "lightgrey")) +
  labs(title = "PERMANOVA R² Values", y = "R²", x = "Pairwise Comparison") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
permanova_plot_file <- paste0(output_prefix, "_permanova_r2_plot.png")
ggsave(permanova_plot_file, plot = permanova_plot, width = 10, height = 8)
#plot Tukey HSD
extract_tukey_hsd <- function(file_path) {
  file_content <- readLines(file_path)
  comparison <- c()
  diff <- c()
  lwr <- c()
  upr <- c()
  p_adj <- c()
  tukey_section <- FALSE
  for (line in file_content) {
    if (grepl("\\$group", line)) {
      tukey_section <- TRUE
    }
    if (tukey_section) {
      if (grepl("^\\s*[A-Za-z]", line)) {
        split_line <- strsplit(line, "\\s+")[[1]]
        comparison <- c(comparison, split_line[1])
        diff <- c(diff, as.numeric(split_line[2]))
        lwr <- c(lwr, as.numeric(split_line[3]))
        upr <- c(upr, as.numeric(split_line[4]))
        p_adj <- c(p_adj, as.numeric(split_line[5]))
      }
    }
  }
  tukey_data <- data.frame(
    Comparison = comparison,
    diff = diff,
    lwr = lwr,
    upr = upr,
    p_adj = p_adj
  )
  
  return(tukey_data)
}
file_path <- paste0(output_prefix, "_combined_results_genus.txt")
tukey_data <- extract_tukey_hsd(file_path)
tukey_data_cleaned <- tukey_data[complete.cases(tukey_data), ]
print(tukey_data_cleaned)
tukey_plot <- ggplot(tukey_data_cleaned, aes(x = Comparison, y = diff)) +
  geom_point(size = 3, color = "blue") +
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Tukey HSD - Group Comparisons", y = "Difference", x = "Comparison") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), axis.title.x = element_text(size = 14))
ggsave(paste0(output_prefix, "_tukey_hsd_plot.png"), plot = tukey_plot, width = 10, height = 8)


#ANCOMBC2
#LogFC: represents log-transformed fold change between X and Y. A positive value indicates genus is more abundant in X
#Standard Error: A smaller error indicates more precise estimation.
#P-value: <0.05 indicates a strong evidence that the abundance of genus is significantly different. 
#Q-value: Adjusted P value to control false discovery rate
#Diff:TRUE states the taxa is differentially abundant
#Passed_SS: TRUE states that the taxa is present in one and absent in the other 

abundance_data <- read.csv(input_file, row.names = 1, check.names = FALSE)
metadata <- data.frame(
  group = rep(c("Control", "PB", "Standard", "Styrofoam"),
              times = c(6, 6, 4, 6)),  
  row.names = colnames(abundance_data),
  stringsAsFactors = FALSE
)
metadata <- DataFrame(metadata)
tax_tab <- data.frame(Genus = rownames(abundance_data))
rownames(tax_tab) <- rownames(abundance_data)
tax_tab <- DataFrame(tax_tab)
tse <- TreeSummarizedExperiment(assays = SimpleList(counts = as.matrix(abundance_data)),
                                colData = metadata,
                                rowData = tax_tab)
colData(tse)$group <- factor(colData(tse)$group, levels = c("Control", "PB", "Standard", "Styrofoam"))
out <- ancombc2(
  data = tse, 
  assay_name = "counts", 
  tax_level = "Genus",  
  fix_formula = "group",
  rand_formula = NULL,
  p_adj_method = "holm", 
  pseudo_sens = TRUE,
  prv_cut = 0,          
  lib_cut = 0,        
  s0_perc = 0.05,
  group = "group", 
  struc_zero = FALSE,  
  neg_lb = TRUE,
  alpha = 0.01, 
  n_cl = 1, 
  verbose = TRUE,
  global = TRUE, 
  pairwise = TRUE,      
  dunnet = FALSE, 
  trend = FALSE,
  iter_control = list(tol = 1e-5, max_iter = 100, verbose = TRUE),
  em_control = list(tol = 1e-5, max_iter = 100),
  lme_control = lme4::lmerControl(),
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 1000)
)
res_pair <- out$res_pair
output_file_res_pair <- paste0(output_prefix, "_res_pair.csv")
write.csv(res_pair, output_file_res_pair)
print(res_pair)
colnames(res_pair) <- c(
  "Genus",
  "Log_FC_PB_vs_Control",                
  "Log_FC_Standard_vs_Control",          
  "Log_FC_Styrofoam_vs_Control",        
  "Log_FC_Standard_vs_PB",               
  "Log_FC_Styrofoam_vs_PB",             
  "Log_FC_Styrofoam_vs_Standard",        
  "SE_PB_vs_Control",                    
  "SE_Standard_vs_Control",              
  "SE_Styrofoam_vs_Control",            
  "SE_Standard_vs_PB",                   
  "SE_Styrofoam_vs_PB",                 
  "SE_Styrofoam_vs_Standard",           
  "W_PB_vs_Control",                     
  "W_Standard_vs_Control",              
  "W_Styrofoam_vs_Control",             
  "W_Standard_vs_PB",                    
  "W_Styrofoam_vs_PB",                   
  "W_Styrofoam_vs_Standard",             
  "P_value_PB_vs_Control",               
  "P_value_Standard_vs_Control",         
  "P_value_Styrofoam_vs_Control",        
  "P_value_Standard_vs_PB",              
  "P_value_Styrofoam_vs_PB",             
  "P_value_Styrofoam_vs_Standard",       
  "Q_value_PB_vs_Control",               
  "Q_value_Standard_vs_Control",         
  "Q_value_Styrofoam_vs_Control",        
  "Q_value_Standard_vs_PB",              
  "Q_value_Styrofoam_vs_PB",             
  "Q_value_Styrofoam_vs_Standard",       
  "Diff_PB_vs_Control",                  
  "Diff_Standard_vs_Control",            
  "Diff_Styrofoam_vs_Control",           
  "Diff_Standard_vs_PB",                 
  "Diff_Styrofoam_vs_PB",                
  "Diff_Styrofoam_vs_Standard",          
  "Passed_SS_PB_vs_Control",             
  "Passed_SS_Standard_vs_Control",       
  "Passed_SS_Styrofoam_vs_Control",      
  "Passed_SS_Standard_vs_PB",           
  "Passed_SS_Styrofoam_vs_PB",          
  "Passed_SS_Styrofoam_vs_Standard"      
)
output_file_descriptive <- paste0(output_prefix, "_res_pair_with_descriptive_names.csv")
write.csv(res_pair, output_file_descriptive)
df <- read.csv(output_file_descriptive)
diff_columns <- c(
  "Diff_PB_vs_Control", "Diff_Standard_vs_Control", "Diff_Styrofoam_vs_Control",
  "Diff_Standard_vs_PB", "Diff_Styrofoam_vs_PB", "Diff_Styrofoam_vs_Standard"
)
q_value_columns <- c(
  "Q_value_PB_vs_Control", "Q_value_Standard_vs_Control", "Q_value_Styrofoam_vs_Control",
  "Q_value_Standard_vs_PB", "Q_value_Styrofoam_vs_PB", "Q_value_Styrofoam_vs_Standard"
)
filter_condition <- rep(FALSE, nrow(df))
for (i in seq_along(diff_columns)) {
  diff_col <- diff_columns[i]
  q_val_col <- q_value_columns[i]
  filter_condition <- filter_condition | (df[[diff_col]] == 1 & df[[q_val_col]] < 0.05)
}
differentially_abundant_genera <- df[filter_condition, ]
output_file_differential <- paste0(output_prefix, "_differentially_abundant_genera.csv")
write.csv(differentially_abundant_genera, output_file_differential)
print(differentially_abundant_genera)
