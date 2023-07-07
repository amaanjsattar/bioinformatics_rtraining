### Q2-Plotting.R ###########################################################################################

# This script corresponds to Question 2 of the Boutros Lab R Training: Plotting

# Load plotting package
library(BoutrosLab.statistics.general);
library(BoutrosLab.plotting.general);
library(BoutrosLab.utilities);

# Load cars as data.frame
cars.data <- as.data.frame(cars);

# Create a simple scatterplot
create.scatterplot(
    filename = generate.filename(
        'RTraining',
        'CarsScatter',
        'tiff'
        ),
    formula = dist ~ speed,
    data = cars.data,
    main = 'Car Stopping Distance by Speed',
    main.cex = 1.5,
    xlab.label = 'Speed (mph)',
    ylab.label = 'Distance (ft)',
    xlab.cex = 2,
    ylab.cex = 2,
    xlimits = 0:30,
    ylimits = 0:130,
    xaxis.cex = 1.5,
    yaxis.cex = 1.5,
    xaxis.fontface = 1,
    yaxis.fontface = 1,
    xaxis.tck = c(1,0),
    yaxis.tck = c(1,0)
    );

# Load Loblolly as data.frame
loblolly.data <- as.data.frame(Loblolly);

# Reshape dataframe: long -> wide
loblolly.wide <- reshape(
    loblolly.data,
    idvar = 'Seed',
    v.names = 'height',
    timevar = 'age',
    direction = 'wide'
    );

# Remove the categorical "Seed" column
loblolly.age.height <- loblolly.wide[, 2:ncol(loblolly.wide)];

# Regex to remove "height" from column names and correct for special character (.)
new.names <- sub('^height\\.', '', colnames(loblolly.age.height));

colnames(loblolly.age.height) <- new.names;

# Utilize the transpose of the dataframe, such that axes are swapped: y = Seed Number, x = Age (yrs)
loblolly.age.height.t <- t(loblolly.age.height);

# Create heatmap
create.heatmap(
    filename = generate.filename(
        'RTraining',
        'LoblollyHeatmap',
        'tiff'
        ),
    x = loblolly.age.height.t,
    clustering.method = 'none',
    main = 'Loblolly Pine Tree Height Variation (ft.) by Age (yrs.) Across Seed Numbers',
    main.cex = 1,
    xlab.label = 'Age (Yrs)',
    xlab.cex = 0.75,
    xaxis.cex = 0.75,
    ylab.label = 'Seed Number',
    ylab.cex = 0.75,
    yaxis.cex = 0.75,
    yaxis.lab = unique(Loblolly$Seed),
    xaxis.fontface = 0.75,
    yaxis.fontface = 0.75,
    xaxis.lab = colnames(loblolly.age.height),
    xaxis.rot = 0,
    grid.row = TRUE,
    grid.col = TRUE,
    row.colour = 'black',
    col.colour = 'black',
    row.lwd = 0.75,
    col.lwd = 0.75,
    colourkey.cex = 0.75,
    same.as.matrix = FALSE,
    colour.scheme = c('white', 'midnightblue')
    );

# Load ChickWeight as data.frame
chicks <- as.data.frame(ChickWeight);

chicks$Time <- as.factor(chicks$Time);

# Check for NA values
any(is.na(ChickWeight));

# Subset by diet, creating four distinct data.frames
diet1.chicks <- subset(chicks, 1 == Diet);

diet2.chicks <- subset(chicks, 2 == Diet);

diet3.chicks <- subset(chicks, 3 == Diet);

diet4.chicks <- subset(chicks, 4 == Diet);

boxplot.diet1 <- create.boxplot(
    formula = weight ~ Time,
    data = diet1.chicks,
    main = 'Diet 1',
    main.cex = 1,
    outliers = FALSE,
    main.y = -1,
    xaxis.cex = 0.5,
    xlab.cex = 0.75,
    xlab.label = NULL,
    ylab.cex = 0.75,
    ylab.label = NULL,
    ylimits = c(0, 350),
    yaxis.cex = 0.5,
    col = 'violet'
    );

boxplot.diet2 <- create.boxplot(
    formula = weight ~ Time,
    data = diet2.chicks,
    main = 'Diet 2',
    main.cex = 1,
    outliers = FALSE,
    main.y = -1,
    main.x = 0.6,
    xlab.cex = 0.75,
    xlab.label = NULL,
    xaxis.cex = 0.5,
    ylab.cex = 0.75,
    ylab.label = NULL,
    ylimits = c(0, 350),
    yaxis.cex = 0.5,
    col = 'forestgreen'
    );

boxplot.diet3 <- create.boxplot(
    formula = weight ~ Time,
    data = diet3.chicks,
    main = 'Diet 3',
    main.cex = 1,
    outliers = FALSE,
    main.y = -1,
    xlab.cex = 0.75,
    xlab.label = NULL,
    xaxis.cex = 0.5,
    ylab.cex = 0.75,
    ylab.label = NULL,
    yaxis.cex = 0.5,
    ylimits = c(0, 350),
    col = 'goldenrod'
    );

boxplot.diet4 <- create.boxplot(
    formula = weight ~ Time,
    data = diet2.chicks,
    main = 'Diet 4',
    main.cex = 1,
    main.y = -1,
    main.x = 0.6,
    outliers = FALSE,
    xlab.cex = 0.75,
    xlab.label = NULL,
    xaxis.cex = 0.5,
    ylab.cex = 0.75,
    ylab.label = NULL,
    ylimits = c(0, 350),
    yaxis.cex = 0.5,
    col = 'purple'
    );


# create panel of 4 boxplots
create.multipanelplot(
    filename = generate.filename(
        'RTraining',
        'ChicksBoxPlots',
        'tiff'
    ),
    plot.objects = list(
        boxplot.diet1,
        boxplot.diet2,
        boxplot.diet3,
        boxplot.diet4),
    layout.width = 2,
    layout.height = 2,
    main = 'Chicken Weights Across 4 Experimental Diets',
    main.cex = 1.5,
    top.padding = 0.5,
    plot.objects.heights = c(1000, 1000),
    x.spacing = 1,
    y.spacing = 1,
    xlab.label = 'Time (days since birth)', xlab.cex = 1.5,
    ylab.label = 'Weight (gm)', ylab.cex = 1.5
    );

save.session.profile(generate.filename('RTraining', 'Q2-Plotting-Session', 'txt'));
