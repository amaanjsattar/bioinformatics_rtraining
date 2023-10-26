### Q1.R ###########################################################################################
# This script corresponds to Question 1 of the Boutros Lab R Training: Statistical Analysis 

# Read the file
AHR.file.path <- "/Users/amaanjsattar/Downloads/AHR-test-file.txt";
mouse.data <- read.table(AHR.file.path, header = TRUE);

# Perform a Paired Sample T-Test on the control and treated samples
t.test(mouse.data$Control, mouse.data$Treated, paired = TRUE);

# Perform a Wilcoxon Signed-Rank test on the control and treated samples
wilcox.test(mouse.data$Control, mouse.data$Treated, paired = TRUE);

# Calculate a fold change between control and treated
# The line below creates a vector of log FC for each mouse (using base 2)
log2(mouse.data$Treated / mouse.data$Control);




