### Q4.R ###########################################################################################
# This script corresponds to Question 4 of the Boutros Lab R Training: Statistical Analysis

# Use the FDR and Bonferroni corrections to adjust p-values for multiple comparisons
library(BoutrosLab.utilities);
library(BoutrosLab.statistics.general);
library(BoutrosLab.plotting.general);

input1_path <- '/Users/amaanjsattar/Downloads/input1.txt';

input2_path <- '/Users/amaanjsattar/Downloads/input2.txt';

tumour.data.p1 <- read.table(input1_path,
                             header = TRUE);
tumour.data.p2 <- read.table(input2_path,
                             header = TRUE);
# Merge the tumour data
tumour.data <- merge(
    tumour.data.p1,
    tumour.data.p2,
    by = "GeneID"
    );

# Get the indices corresponding to subtype 1
subtype1.inds <- seq(
    from = 2,
    to = ncol(tumour.data.p1)
    );

# Get the indices corresponding to subtype 2
subtype2.inds <- seq(
    from = (ncol(tumour.data.p1) + 1),
    to = ncol(tumour.data)
    );

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

wilcox.p.values <- apply(
    tumour.data,
    MARGIN = 1,
    FUN = tumour.data.wilcox,
    subtype1.inds,
    subtype2.inds
    );

# Implement FDR Correction
wilcox.p.values.fdr <- p.adjust(
    p = wilcox.p.values,
    method = 'fdr',
    n = length(wilcox.p.values)
    );

# Create FDR-corrected histogram
create.histogram(
    x = wilcox.p.values.fdr,
    main = 'FDR-Adjusted Wilcoxon Signed-Rank Test P-Values',
    main.cex = 1.3,
    xlab.label = 'P-Value',
    xlab.cex = 1,
    xaxis.cex = 1,
    ylab.label = 'Frequency (%)',
    ylab.cex = 1,
    yaxis.cex = 1,
    ylimits = c(0, 90),
    col = 'red3',
    );

# Implement Bonferroni Correction
wilcox.p.values.bonf <- p.adjust(
    p = wilcox.p.values,
    method = 'bonferroni',
    n = length(wilcox.p.values)
    );

# Create Bonferroni-corrected histogram
create.histogram(
    x = wilcox.p.values.bonf,
    main = 'Bonferroni-Adjusted Wilcoxon Signed-Rank Test P-Values',
    main.cex = 1.25,
    xlab.label = 'P-Value',
    xlab.cex = 1,
    xaxis.cex = 1,
    ylab.label = 'Frequency (%)',
    ylab.cex = 1,
    yaxis.cex = 1,
    col = 'red3'
    );

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
    subtype2.inds);

# Implement FDR Correction
t.test.p.values.fdr <- p.adjust(
    p = t.test.p.values,
    method = 'fdr',
    n = length(t.test.p.values));

# Create FDR-corrected histogram
create.histogram(
    x = t.test.p.values.fdr,
    main = 'FDR-Adjusted T-Test P-Values',
    main.cex = 2,
    xlab.label = 'P-Value',
    xlab.cex = 1.5,
    xaxis.cex = 1,
    ylab.label = 'Frequency (%)',
    ylab.cex = 1.5,
    yaxis.cex = 1,
    ylimits = c(0, 90),
    col = 'orange2'
    );

# Implement Bonferroni Correction
t.test.p.values.bonf <- p.adjust(
    p = t.test.p.values,
    method = 'bonferroni',
    n = length(t.test.p.values));

# Create Bonferroni-corrected histogram
create.histogram(
    x = t.test.p.values.bonf,
    main = 'Bonferroni-Adjusted T-Test P-Values',
    main.cex = 2,
    xlab.label = 'P-Value',
    xlab.cex = 1.5,
    xaxis.cex = 1,
    ylab.label = 'Frequency (%)',
    ylab.cex = 1.5,
    yaxis.cex = 1,
    col = 'orange2'
    );

# Written Portion

save.session.profile(generate.filename('RTraining', 'Q4-statistics', 'txt'), stdout = TRUE);
