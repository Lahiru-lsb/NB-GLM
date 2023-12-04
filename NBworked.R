

#Installing packages


library(tidyverse)
library(spatialLIBD)
library(dbscan)
library(scran)

library(igraph)
library(CellBench)
library(scater)
library(MASS)

library(SpatialExperiment)
library(sctransform)

library(magrittr)
library(ggplot2)
library(gridExtra)
library(cluster)


#Attaching the data

data = readRDS(file = "C:\\Users\\wboyagoda\\OneDrive - The University of Melbourne\\Documents\\Spatial Data Analysis\\Data sets\\ST_for_Lahiru.rds")
seo=data[[1]]
seo <- seo[, !is.na(colData(seo)$AnnotatedCluster)]
cleanData=seo




filtering=function(spatial_obj, lib_size_thresh, gene_min_thresh, gene_max_thresh)
{
  qc=perCellQCMetrics(spatial_obj)
  qc_df=as.data.frame(qc)
  
  size_threshold=lib_size_thresh
  min_threshold=gene_min_thresh
  max_threshold=gene_max_thresh
  
  filtered_cells=qc$sum>lib_size_thresh &
    qc$detected>=gene_min_thresh &
    qc$detected<= gene_max_thresh
  
  #filter by column
  filtered_seo=spatial_obj[,filtered_cells]
  return(filtered_seo)
  
}




ls_normalization <- function(SEO) {
  if (!inherits(SEO, 'SpatialExperiment')) {
    stop("Input is not a SpatialExperiment object.")
  }
  
  # Calculate library sizes (sum of counts per spot)
  ls <- colSums(counts(SEO))
  
  # Calculate mean library size
  mean_ls <- mean(ls)
  
  # Calculate size factors
  size_fac <- ls / mean_ls
  
  # Apply normalization
  normalized_counts <- log1p(sweep(counts(SEO), 2, size_fac, FUN = "/"))
  
  # Update the counts in the SpatialExperiment object's 'logcounts' slot
  assay(SEO, "logcounts") <- normalized_counts
  
  return(SEO)
}




select_features <- function(SEO, n = 1000, assay_name = "logcounts") {
  if (!is(SEO, 'SpatialExperiment')) {
    stop("Not a Spatial Experiment object.")
  }
  
  if (!assay_name %in% assayNames(SEO)) {
    stop("Specified assay not found in the SpatialExperiment object.")
  }
  
  # Compute model var. of gene's expression across spatial locations using the specified assay
  gene_var <- scran::modelGeneVar(SEO, assay.type = assay_name)
  
  # Select top n most variable genes
  hvg <- getTopHVGs(gene_var, n = n)
  
  # Subset the SpatialExperiment object
  spe_hvg <- SEO[hvg, ]
  
  return(spe_hvg)
}




data <- cleanData[rowSums(counts(cleanData)) > 0, ]


data=filtering(data, 500, 500, 2500)
data_lsnorm=ls_normalization(data)
#feature selections

top_500=select_features(data_lsnorm)



counts_matrix=assay(top_500,'counts')

# Calculate library sizes (sum of counts per spot)
library_sizes <- colSums(counts(top_500))
logLS <- log(library_sizes )



AnnotatedCluster <- colData(top_500)$AnnotatedCluster



# Initialize a list to store models
nb_models <- vector("list", nrow(counts_matrix))

# Iterate over each gene
for (i in 1:nrow(counts_matrix)) {
  
  
  # Extract gene expression counts for the current gene
  gene_counts <- counts_matrix[i, ]
  
  # Create a dataframe for the model
  df <- data.frame(gene_expression = gene_counts,
                   AnnotatedCluster = factor(AnnotatedCluster),
                   logLS = logLS)
  
  # Check for non-finite values in each column of the data frame
  if (any(!is.finite(df$gene_expression)) || any(!is.finite(df$logLS))) {
    # Skip this iteration if non-finite values are found
    next
  }
  
  # Fit the Negative Binomial GLM with error handling
  nb_models[[i]] <- tryCatch({
    glm.nb(gene_expression ~  AnnotatedCluster+logLS, data = df,control = glm.control(maxit = 1000))
  }, error = function(e) {
    NULL  # Return NULL in case of an error
  })
}


# Count how many models were successfully fitted
successful_models <- sum(sapply(nb_models, function(x) !is.null(x)))

cat('successful_models (without spatial data)',successful_models)







top_500$Lib_size=log(colSums(counts(top_500)))
library(splines)

# Assuming 'x' and 'y' are your spatial coordinates
x <- spatialCoords(top_500)[,1]  
y <- spatialCoords(top_500)[,2]

# Generate B-Splines for x and y coordinates
bs.x <- ns(x, df=7) 
bs.y <- ns(y, df=7) 


bs.xy<- matrix(0,nrow=nrow(bs.x),ncol=ncol(bs.x)*ncol(bs.y))
for(i in 1:ncol(bs.x)) {
  for(j in 1:ncol(bs.y)) {
    bs.xy[,(i-1)*ncol(bs.x)+j] <- bs.x[,i]*bs.y[,j]
  }
}

nb_spatial_models <- vector("list", nrow(counts_matrix))

for (i in 1:nrow(counts_matrix)) {
  
  gene_counts <- counts_matrix[i, ]
  
  
  nb_spatial_models[[i]] <- tryCatch({
    
    glm.nb(assays(top_500)$count[i,] ~ AnnotatedCluster+Lib_size*bs.xy,
           data = colData(top_500),control = glm.control(maxit = 1000))
  }, error = function(e) {
    # Print an error message and return NULL
    message("Error in model fitting for gene index ", i, ": ", e$message)
    NULL
  })
}


successful_models <- sum(sapply(nb_spatial_models, function(x) !is.null(x)))
successful_models




#handeling the different number of model fittings

gene_names <- rownames(counts_matrix)  

successful_indices_model1 <- which(!sapply(nb_models, is.null))
successful_indices_model2 <- which(!sapply(nb_spatial_models, is.null))

common_indices <- intersect(successful_indices_model1, successful_indices_model2)



nb_models_filtered <- nb_models[common_indices]
nb_spatial_models_filtered <- nb_spatial_models[common_indices]
filtered_gene_names <- gene_names[common_indices]



# Set a significance threshold
significance_threshold <- 0.05

# Load necessary library
library(lmtest)

# Initialize a vector to store p-values
p_values <- numeric(length(common_indices))

# Iterate over each gene to perform the likelihood-ratio test
for (i in 1:length(common_indices)) {
  model1 <- nb_models_filtered[[i]]
  model2 <- nb_spatial_models_filtered[[i]]
  
  # Perform the likelihood-ratio test if both models are fitted
  if (!is.null(model1) && !is.null(model2)) {
    test_result <- lrtest(model1, model2)
    p_values[i] <- test_result$'Pr(>Chisq)'[2]  # Extract p-value
  } else {
    p_values[i] <- NA  # Assign NA if either model is NULL
  }
}

# Load the qvalue library
library(qvalue)

# Filter out non-finite values from p-values
valid_p_values <- p_values[is.finite(p_values)]

# Calculate q-values and analyze results if there are valid p-values
if (length(valid_p_values) > 0) {
  # Calculate q-values
  q_values <- qvalue(valid_p_values)$qvalues
  
  # Set a threshold for significance
  qvalue_threshold <- 0.05
  
  # Count significant q-values and calculate proportion
  significant_qvalues <- sum(q_values < qvalue_threshold, na.rm = TRUE)
  proportion_significant_qvalues <- significant_qvalues / length(q_values)
  
  # Print the proportion of significant q-values
  cat('proportion_significant_qvalues', proportion_significant_qvalues, '\n')
  
  # Count how many p-values are below the threshold and calculate proportion
  significant_genes <- sum(valid_p_values < significance_threshold, na.rm = TRUE)
  proportion_significant <- significant_genes / length(valid_p_values)
  
  # Print the proportion of significant genes
  cat("proportion_significant ", proportion_significant, '\n')
} else {
  # Handle the case where there are no valid p-values
  warning("No valid p-values available for q-value analysis.")
}

