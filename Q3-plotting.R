### Q3-Plotting.R ###########################################################################################

# This script corresponds to Question 3 of the Boutros Lab R Training: Plotting

# Load packages
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Step 1:
# Read in the data
seqcontrol.data <- read.csv(
    file = '/Users/amaanjsattar/Downloads/Q3_SeqControl_data.tsv',
    sep = '\t',
    );

# Reorder data based on Yes Votes in decreasing order
seqcontrol.ordered <- seqcontrol.data[order(seqcontrol.data$yes.votes,
                                            decreasing = TRUE), ];

# Step 2
# Create a heatmap displaying cpcgene sample names
sample.colour.scheme <- c('rosybrown1', 'rosybrown4', 'red',
                          'darkred', 'darkorange', 'gold',
                          'darkolivegreen3', 'darkgreen', 'aquamarine',
                          'cyan4', 'dodgerblue', 'darkblue'
                          );

sample.names <- c('CPCG0003P', 'CPCG0005P', 'CPCG0007P',
                  'CPCG0040P', 'CPCG0047P', 'CPCG0603P',
                  'CPCG0098P', 'CPCG0102P', 'CPCG0103P',
                  'CPCG0123P', 'CPCG0183P', 'CPCG0184P'
                  );


cpcg <- seqcontrol.ordered$CPCG;

# Factorize CPCG indices
cpcg.factors <- as.numeric(
    as.factor(cpcg)
    );

# Create a matrix (heatmap input) of cpcg factors
cpcg.matrix <- matrix(
    data = cpcg.factors,
    nrow = 1
    );

# Create heatmap of cpcg
# Add 1 to total.colours to account for a "NULL" colour
cpcg.heatmap <- create.heatmap(
    x = cpcg.matrix,
    clustering.method = 'none',
    colour.scheme = sample.colour.scheme,
    grid.col = TRUE,
    force.grid.col = TRUE,
    col.colour = 'black',
    total.colours = length(unique(cpcg.factors)) + 1,
    yaxis.tck = 0,
    xaxis.tck = 0,
    print.colour.key = FALSE
    );
    
# Step 3, Part 1
average.reads.start <- seqcontrol.ordered$Average.reads.start;

# Create a matrix (heatmap input) of average.reads.start
average.reads.start.matrix <- matrix(
    data = average.reads.start,
    nrow = 1
    );

# Create heatmap of average.reads.start
average.reads.start.heatmap <- create.heatmap(
    x = average.reads.start.matrix,
    clustering.method = 'none',
    colour.scheme = c('white', 'deeppink'),
    print.colour.key = FALSE,
    grid.col = TRUE,
    force.grid.col = TRUE,
    xaxis.tck = 0,
    yaxis.tck = 0,
    );

# Step 3, Part 2
unique.start.points <- seqcontrol.ordered$Unique.start.points;

# Create a matrix (heatmap input) of unique.start.points
unique.start.points.matrix <- matrix(
    data = unique.start.points,
    nrow = 1
    );

# Create heatmap of unique.start.points
unique.start.points.heatmap <- create.heatmap(
    x = unique.start.points.matrix,
    clustering.method = 'none',
    colour.scheme = c('white', 'darkblue'),
    print.colour.key = FALSE,
    grid.col = TRUE,
    force.grid.col = TRUE,
    xaxis.tck = 0,
    yaxis.tck = 0
    );

# Step 3, Part 3
x.bases.0.quality <- seqcontrol.ordered$X..Bases...0.quality;

# Create a matrix (heatmap input) of X..Bases...O.quality
x.bases.0.quality.matrix <- matrix(
    data = x.bases.0.quality,
    nrow = 1
    );

# Create a heatmap of X..Bases...0.quality
x.bases.0.quality.heatmap <- create.heatmap(
    x = x.bases.0.quality.matrix,
    clustering.method = 'none',
    colour.scheme = c('white', 'darkorange'),
    print.colour.key = FALSE,
    grid.col = TRUE,
    force.grid.col = TRUE,
    xaxis.tck = 0,
    yaxis.tck = 0
    );

# Step 4
# CPCG0102P and CPCG0103P are known FFPE samples
FFPE.samples <- c('CPCG0102P', 'CPCG0103P');

# Apply boolean filtering to differentiate between FFPE and non-FFPE samples
is.FFPE <- ifelse(
    seqcontrol.ordered$CPCG %in% FFPE.samples,
    1, 0
    );

# Create a matrix (heatmap input) of is.FFPE
is.FFPE.matrix <- matrix(
    data = is.FFPE,
    nrow = 1
    );

# Create a heatmap of is.FFPE
is.FFPE.heatmap <- create.heatmap(
    x = is.FFPE.matrix,
    clustering.method = 'none',
    colour.scheme = c('white', 'darkslategrey'),
    print.colour.key = FALSE,
    grid.col = TRUE,
    force.grid.col = TRUE,
    xaxis.tck = 0,
    yaxis.tck = 0
    );

# Step 5
prop.yes.votes <- seqcontrol.ordered$yes.votes;

# Apply boolean filtering to create bar colour scheme
bar.colour.scheme <- ifelse(
    1 == seqcontrol.ordered$outcome,
    'black', 'grey'
    );
    
# Create barplot of fraction of yes votes
yes.vote.barplot <- create.barplot(
    formula = prop.yes.votes ~ 1:length(prop.yes.votes),
    data = seqcontrol.ordered,
    ylab.label = 'Fraction of yes votes',
    ylab.cex = 1.5,
    yaxis.cex = 1,
    xlab.label = NULL,
    xaxis.lab = FALSE,
    xaxis.tck = 0,
    xaxis.cex = 0,
    yaxis.tck = 1,
    yat = seq(0, 1, 0.25),
    ylimits = c(0, 1.05),
    box.ratio = 2,
    col = bar.colour.scheme,
    abline.h = 0.5,
    abline.lty = 2,
    abline.col = 'grey',
    add.grid = FALSE,
    border.col = NULL,
    );

# Step 6
# Create legend.grob input as a list of individual legends
legends <- legend.grob(
    list(
        legend = list(
            title = expression(underline('Sample')),
            colours = sample.colour.scheme,
            labels = sample.names,
            size = 1.5,
            border = 'black'
        ),
        legend = list(
            title = expression(underline('Sample Preparation')),
            colours = c('white','darkslategrey'),
            labels = c('Frozen', 'FFPE'),
            size = 1.5,
            border = 'black'
        ),
        legend = list(
            colours = c('white', 'darkorange'),
            labels = c('83.0','97.0'),
            border = 'black',
            title = expression(underline('% Bases > 0 quality')),
            continuous = TRUE
        ),
        legend = list(
            colours = c('white', 'darkblue'),
            labels = scientific.notation(c(100000000, 360000000)),
            border = 'black',
            title = expression(underline('Unique start points')),
            continuous = TRUE
        ),
        legend = list(
            colours = c('white','deeppink'),
            labels = c('1.070', '1.190'),
            border = 'black',
            title = expression(underline('Average reads/start')),
            continuous = TRUE
        )
    ),
    title.just = 'left',
    label.cex = 0.7,
    title.cex = 0.7
    );

# Create barplot legend separately
barplot.legend <- legend.grob(
    list(
        legend = list(
            colours = c('grey', 'black'),
            labels = c('< 50x', '\u2265 50x'),
            border = 'black',
            title = expression(underline('Observed'))
        )
    ),
    title.just = 'left',
    label.cex = 1,
    title.cex = 1,
    x = 0.65,
    y = 0.8
    );

# Step 7
# Combine all individual plots and legend using create.multipanelplot

create.multipanelplot(
    filename = generate.filename(
        'RTraining',
        'SeqControl',
        'tiff'),
    height = 8,
    width = 8,
    x.spacing = 0,
    y.spacing = -0.5,
    ylab.axis.padding = 2,
    plot.objects = list(
        yes.vote.barplot,
        cpcg.heatmap,
        is.FFPE.heatmap,
        x.bases.0.quality.heatmap,
        unique.start.points.heatmap,
        average.reads.start.heatmap
        ),
    plot.objects.heights = c(1, 0.15, 0.15, 0.15, 0.15, 0.15),
    legend = list(
        right = list(fun = legends),
        inside = list(fun = barplot.legend)
        )
    );

save.session.profile(generate.filename('RTraining', 'Q3-Plotting-Session', 'txt'));
