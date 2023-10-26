### Q2.R ###########################################################################################
# This script corresponds to Question 2 of the Boutros Lab R Training: Statistical Analysis 

# 1) Read the files
input1.path <- "/Users/amaanjsattar/Downloads/input1.txt";
input2.path <- "/Users/amaanjsattar/Downloads/input2.txt";
tumour.data.p1 <- read.table(input1.path, header = TRUE);
tumour.data.p2 <- read.table(input2.path, header = TRUE);

# Create vectors for the column names corresponding to subtypes 1 and 2, respectively
subtype1.cols <- colnames(tumour.data.p1[, -1])
subtype2.cols <- colnames(tumour.data.p2[, -1])
    
# 2: Combine the files into one file with data on all 12 tumors 
# [order: tumour.data.p1, tumour.data.p2]

# Sort the tumor data based on GeneID 
sorted.tumour.data.p1 <- tumour.data.p1[order(tumour.data.p1[, 1]), ];
sorted.tumour.data.p2 <- tumour.data.p2[order(tumour.data.p2[, 1]), ];



# 2a) Combination method 1 [using cbind]
tumour.data.combined <- cbind(sorted.tumour.data.p1, sorted.tumour.data.p2[, -1]);

# 2b) Combination method 2 [using merge]
tumour.data.merged <- merge(sorted.tumour.data.p1, sorted.tumour.data.p2, by = "GeneID");

# Verify that the two tables are identical 
identical(tumour.data.combined, tumour.data.merged);

# 3) Conduct a t-test comparing the first 3 tumors to the last 9 for *each* gene using a for-loop
# Create an empty table to store our final p values for each gene
gene.p.values <- data.frame("GeneID" = tumour.data.merged$GeneID, 
                            "P.Value" = numeric(nrow(tumour.data.merged)));

for (gene in 1:nrow(tumour.data.merged)) {
  
    # Get data for each individual gene 
    gene.data <- tumour.data.merged[gene, ];
  
    # subset the table by subtype

    subtype.1 <- gene.data[subtype1.cols];
    subtype.2 <- gene.data[subtype2.cols];
  
    # conduct t test comparing subtypes
    gene.wise.t.test <- t.test(subtype.1, subtype.2);
  
    # extract p-values and store in table
    gene.p.values[gene, ] <- gene.wise.t.test$p.value;
}

# 4. plot a histogram of the p-values
hist(gene.p.values$P.Value);

# 5. fix axes labels if they are rotated 90 degrees
hist(gene.p.values$P.Value, las = 1);

# 6. Try plotting the histogram in log-space
hist(log(gene.p.values$P.Value),
     las = 1);

# 7. interpret the distribution;
