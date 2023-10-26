### Q3.R ###########################################################################################
# This script corresponds to Question 3 of the Boutros Lab R Training: Statistical Analysis 

# Read the files
input1_path <- "/Users/amaanjsattar/Downloads/input1.txt"

input2_path <- "/Users/amaanjsattar/Downloads/input2.txt"

tumour.data.p1 <- read.table(input1_path, 
                             header = TRUE)
tumour.data.p2 <- read.table(input2_path, 
                             header = TRUE)

# Merge the tumour data
tumour.data <- merge(
    tumour.data.p1, 
    tumour.data.p2, 
    by = "GeneID")

# Verify that the first column in your tumour data is "GeneID"
colnames(tumour.data[1]) == "GeneID"

# Get the indices corresponding to subtype 1
subtype1.inds <- seq(
    from = 2, 
    to = ncol(tumour.data.p1))

# Get the indices corresponding to subtype 2
subtype2.inds <- seq(
    from = (ncol(tumour.data.p1) + 1),
    to = ncol(tumour.data))

### FUNCTION #######################################################################################
# Input variables:
# gene.row:  a row vector with gene expression data
# inds1:    a vector of indices corresponding to the columns in the first comparison group
# inds2:    a vector of indices corresponding to the columns in the second comparison group
# Output variables:
# p.vals:   a vector of p-values generated from the Wilcoxon Signed-Rank Test
# Description: 
# Returns p-values generated from a Wilcoxon Signed-Rank Test applied to
# each row of a table, with comparison groups denoted by user-provided indices

# Create wrapper function
tumour.data.wilcox <- function(gene.row, inds1, inds2) {
    # Run wilcox.test with comparison groups specified by inds1, inds2
    wilcox.compare <- wilcox.test(
        as.numeric(gene.row[inds1]), 
        as.numeric(gene.row)[inds2],
        paired = FALSE)
    
    # Return a vector of p-values
    p.vals <- wilcox.compare$p.value
    return(p.vals)
}

# Apply the function to tumour.data and store p-values in a vector
wilcox.p.values <- apply(
    tumour.data, 
    MARGIN = 1, 
    FUN = tumour.data.wilcox, 
    subtype1.inds,
    subtype2.inds)

### FUNCTION #######################################################################################
# Input Variables:
# gene.row: a row vector with gene expression data
# inds1:    a vector of indices corresponding to the columns in the first comparison group
# inds2:    a vector of indices corresponding to the columns in the second comparison group
# Output variables:
# p.vals:   a vector of p-values generated from the T-Test
# Description: 
# Returns p-values generated from a T-Test applied to
# each row of a table, with comparison groups denoted by user-provided indices


# Create wrapper function
tumour.data.t.test <- function(gene.row, inds1, inds2) {
    # Run t.test with comparison groups specified by inds1, inds2
    t.test.compare <- t.test(
        as.numeric(gene.row[inds1]), 
        as.numeric(gene.row)[inds2],
        paired = FALSE)
    
    # Return a vector of p-values
    p.vals <- t.test.compare$p.value
    return(p.vals)
}

# Apply the function to tumour.data and store p-values in a vector
t.test.p.values <- apply(
    tumour.data, 
    MARGIN = 1, 
    FUN = tumour.data.t.test,
    subtype1.inds, 
    subtype2.inds)

### FUNCTION #######################################################################################
tumour.data.fold.changes <- function(gene.row, inds1, inds2) {
    # Calculate mean expression levels for each gene, separated by subtype
    subtype1.mean.exp <- mean(as.numeric(gene.row[inds1]))
    subtype2.mean.exp <- mean(as.numeric(gene.row[inds2]))
    
    # Calculate the fold changes 
    subtype.fold.changes <- (subtype1.mean.exp / subtype2.mean.exp)
    
    # Return a vector of fold changes
    return(subtype.fold.changes)
}
    
tumour.fold.changes <- apply(
    tumour.data, 
    MARGIN = 1, 
    FUN = tumour.data.fold.changes,
    subtype1.inds, 
    subtype2.inds)



# Show the non-transformed histograms
hist(tumour.p.values)

hist(wilcox.p.values)

hist(tumour.fold.changes)

# Apply log(base 2) transformations to each histogram
hist(log(t.test.p.values))

hist(log2(wilcox.p.values))

hist(log2(tumour.fold.changes))
