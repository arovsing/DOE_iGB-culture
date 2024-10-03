## Rovsing et al. 2024
## scripts for generating heatmap in figure 6

library(data.table)
library(tidyverse)
library(pheatmap)

# Load the data
df = fread("FlowSOM_heatmap.csv")

# Preprocess the data: Remove the last two rows and V27
df_prep = df %>%
    select(-V27) %>%
    slice(1:(n() - 2))

# Normalize the data (subtract median and divide by standard deviation)
df_normalized <- df_prep %>%
    mutate(across(-V1, ~ (. - median(.)) / sd(.)))

# Prepare the data for clustering
df_wide <- df_normalized %>%
    select(-V1) %>% # Remove the V1 column for clustering
    as.data.frame()

# Generate new row names in the format Pop_0 to Pop_19
new_row_names <- paste0("Pop_", seq_len(nrow(df_wide)) - 1)
rownames(df_wide) <- new_row_names

# Clean up column names by removing the specific substring
colnames(df_wide) <- gsub(" \\| Freq. of Parent \\(%\\)", "", colnames(df_wide))

# Transpose the data for flipping the heatmap visually
df_wide_t <- t(df_wide)

# Compute the distance matrix and apply hierarchical clustering on columns (which are original rows)
dist_matrix_col <- dist(df_wide_t)
hc_col <- hclust(dist_matrix_col)

# Define the number of clusters for columns (original rows)
num_col_clusters <- 3

# Assign clusters to columns (original rows)
col_clusters <- cutree(hc_col, k = num_col_clusters)

# Create annotation data frame for columns (original rows)
col_annotation <- data.frame(Cluster = factor(col_clusters[rownames(df_wide_t)]))
col_annotation <- col_annotation[match(rownames(df_wide_t), rownames(col_annotation)), , drop = FALSE]

# Create color palette for column annotations
cluster_colors <- c(
    "1" = "#a49250",
    "2" = "#a24131",
    "3" = "#5f85ae"
)

# Map cluster numbers to colors
col_annotation_colors <- setNames(cluster_colors, levels(col_annotation$Cluster))

# Extract values for x-axis labels and include 'Pop_'
x_axis_labels <- colnames(df_wide_t)  # These are originally row names in df_wide

# Create the heatmap with pheatmap and annotations
pheatmap(
    df_wide_t,
    cluster_rows = hc_col, # Cluster rows (which are original columns)
    cluster_cols = FALSE, # Do not cluster columns (which are original rows)
    color = colorRampPalette(c("#f9a042", "#f6e8a7", "#a19c93", "#1a3f71"))(100),
    show_rownames = TRUE,
    show_colnames = TRUE,
    cutree_rows = num_col_clusters,
    annotation_row = col_annotation,
    annotation_colors = list(
        Cluster = col_annotation_colors
    ),
    labels_col = x_axis_labels # Update column labels to include 'Pop_' with numbers
)
