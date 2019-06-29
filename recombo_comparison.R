#!/applications/R/R-3.5.0/bin/Rscript

### Plot and compare recombination frequencies (RF) for two genotypes in various intervals

## Install the ggplot2 and ggbeeswarm packages if not already installed
# (first check if already installed with library(ggplot2) and library(ggbeeswarm);
#  if there are no error messages, they're already installed.
#  If not installed, use the following commented-out commands)
#install.pacakges("ggplot2")
#install.packages("ggbeeswarm")
# Load packages to be used for data processing and plotting
library(ggplot2)
library(ggbeeswarm)

## Set working directory to directory containing RF CSV files
setwd("/path/to/directory/containing/CSVs/")
setwd("/home/ajt200/recombo_comparison/")
print(getwd())

## Specify genotype names to be used for plotting
geno1Name <- "Col-0"
geno2Name <- "Tanz-1"
## Read in CSV files as data.frames
# "./" represents the current working directory, as set above
geno1 <- read.csv("./SeedScoringResultsTanz-1_FTLs_Col_Col.csv")
geno2 <- read.csv("./SeedScoringResultsTanz-1_FTLs_Col_Tanz-1.csv")

## Have a look at the first and last 6 lines of each data.frame
print(head(geno1))
print(tail(geno1))
print(head(geno2))
print(tail(geno2))

## Retain relevant columns
df <- data.frame(genotype = c(rep(geno1Name, dim(geno1)[1]),
                              rep(geno2Name, dim(geno2)[1])),
                 # Use a regular expression pattern to extract
                 # the "interval" part of the Sample ID:
                 # ".*_" matches any character up to and including the first underscore
                 # "(.*)" represents the interval name,
                 # with brackets enabling extraction using "\\1"
                 # "_\\d*" matches the number after the second underscore
                 interval = c(sub(pattern = ".*_(.*)_\\d*",
                                  replacement = "\\1",
                                  x = geno1$Sample),
                              sub(pattern = ".*_(.*)_\\d*",
                                  replacement = "\\1",
                                  x = geno2$Sample)),
                 # Use a regular expression pattern to extract
                 # the "individual" part of the Sample ID:
                 # ".*_" matches any character up to and including the first underscore
                 # The subsequent ".*_" matches any character up to and including the second underscore
                 # "(\\d*)" matches the individual numbers,
                 # with brackets enabling extraction using "\\1"
                 individual = c(sub(pattern = ".*_.*_(\\d*)",
                                    replacement = "\\1",
                                    x = geno1$Sample),
                                sub(pattern = ".*_.*_(\\d*)",
                                    replacement = "\\1",
                                    x = geno2$Sample)),
                 RF = c(geno1$RF...,
                        geno2$RF...))

## Have a look at df
print(df)

## Order the levels of df$genotype (geno1Name followed by geno2Name)
## so that geno1Name is plotted before geno2Name
df$genotype <- factor(df$genotype,
                      levels = c(geno1Name, geno2Name))

## Compare RF values using Mann-Whitney U tests
# Here lapply() iterates over each unique interval name in df
# and applies a U test comparing geno1 and geno2 RFs 
# This creates a list object (Utests) in which each list element
# (accessible with double brackets; e.g., Utests[[1]])
# contains the results of a U test for a given interval
Utests <- lapply(seq_along(levels(df$interval)), function(x) {
  wilcox.test(x = df[df$interval == levels(df$interval)[x] &
                     df$genotype == geno1Name,]$RF,
              y =  df[df$interval == levels(df$interval)[x] &
                      df$genotype == geno2Name,]$RF,
              alternative = "two.sided")
})
# Note warning message: "cannot compute exact p-value with ties"
# This occurred because there are identical RF values within a given interval
# for a given genotype. These values therefore receive the same ranking
# in the U test, which prevents calculation of exact P-values

## Name list elements of "Utests" with interval names to make them
## accessible by interval name (e.g., Utests$`420-CEN3`)
names(Utests) <- levels(df$interval)
print(Utests)

## Extract U test P-values
# sapply() works similarly to lapply(), but the resulting object is
# a vector rather than a list of elements
# (specific elements of a vector are accessible with single brackets;
# e.g., UtestPvals[1])
UtestPvals <- sapply(seq_along(Utests), function(x) {
  Utests[[x]]$p.val
})
names(UtestPvals) <- levels(df$interval)
print(UtestPvals)

## Disable scientific notation
options(scipen=999)

## Simplify P-values
UtestPvalsChar <- sapply(seq_along(UtestPvals), function (x) {
  if(UtestPvals[x] < 0.0001) {
    "< 0.0001"
  } else {
    paste0("= ", as.character(round(UtestPvals[x], digits = 4)))
  }
})
names(UtestPvalsChar) <- levels(df$interval)
print(UtestPvalsChar)

## RF plotting function
RFplotFun <- function(dataFrame, intervalName, Pval, genoColours) {
  RFplot <- ggplot(data = dataFrame[dataFrame$interval == intervalName,],
                   mapping = aes(x = genotype,
                                 y = RF,
                                 colour = genotype)) +
            scale_colour_manual(values = genoColours) +
            #geom_boxplot(mapping = aes(colour = genotype),
            #             varwidth = TRUE) +
            geom_violin(scale = "count",
                        trim = FALSE,
                        draw_quantiles = c(0.25, 0.50, 0.75)) +
            geom_beeswarm(cex = 6,
                          size = 4) +
            labs(x = "",
                 y = "Recombination frequency (cM)") +
            theme_bw() +
            theme(axis.line.y = element_line(size = 1.0, colour = "black"),
                  axis.ticks.y = element_line(size = 1.0, colour = "black"),
                  axis.ticks.x = element_blank(),
                  axis.text.y = element_text(size = 18, colour = "black"),
                  axis.text.x = element_text(size = 18, colour = genoColours),
                  axis.title = element_text(size = 20, colour = "black"),
                  legend.position = "none",
                  panel.grid = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  plot.margin = unit(c(0.3,0.9,0.0,0.3), "cm"),
                  plot.title = element_text(hjust = 0.5, size = 30)) +
            ggtitle(bquote(italic(.(intervalName)) ~ ~ ~ ~
                           "MWW" ~ italic("P") ~ .(Pval)))
  ggsave(RFplot,
         file = paste0("./", intervalName, "_",
                       geno1Name, "_vs_", geno2Name,
                       "_RF_Utest_and_plot.pdf"),
         width = 18, height = 18, units = "cm")
}

## Iterate over each unique interval name in df
# and apply plotting function 
sapply(seq_along(levels(df$interval)), function(x) {
  RFplotFun(dataFrame = df,
            intervalName = levels(df$interval)[x],
            Pval = UtestPvalsChar[x],
            genoColours = c("darkorange", "dodgerblue2"))
})
