###Plotting Twisst Results###

#Based on the plotting script provided with Twisst from Simon H. Martin (https://github.com/simonhmartin/twisst/blob/master/example_plot.R)
rm(list = ls())

library(ape)
library(data.table)
library(tools)
library(ggplot2)
library(cowplot)
library(magick)
library(dplyr)
library(R.utils)

#Load the plotting functions, this script also listed in this repo with some modifications from the original
source("plot_twisst.R")

############################## input files ######################################

# It is possible to import one or more weights files at a time.
# Here we just import 1
#weights file with a column for each topology
weights_file <- "outputislands1.weights100.csv.gz"
#Also read in the weights as a dataframe
weightsdf<-read.csv("outputislands1.weights100.csv.gz", header = F)
# It is not necessary to import window data files, but if you do there should be one for
# each weights file

#coordinates file for each window
window_data_file <- "output.islands1phyml_bioml.w100.data.tsv"

################################# import data ##################################

# The function import.twisst reads the weights, window data  files into a list object
# If there are multiple weights files, or a single file with different chromosomes/scaffolds/contigs
# in the window data file, these will be separated when importing.

twisst_data <- import.twisst(weights_files=weights_file,
                             window_data_files=window_data_file)

#Change the island names in the topologies to abbreviations
for (i in 1:length(twisst_data$topos)){
  twisst_data$topos[[i]]$tip.label<-as.character(ifelse(twisst_data$topos[[i]]$tip.label == "Outgroup", 'O',
                                                        ifelse(twisst_data$topos[[i]]$tip.label == "Darwin", 'D',
                                                               ifelse(twisst_data$topos[[i]]$tip.label == "Gardner", 'Ga',
                                                                      ifelse(twisst_data$topos[[i]]$tip.label == "Isabela", 'I',
                                                                             ifelse(twisst_data$topos[[i]]$tip.label == "Pinta", 'P',
                                                                                    ifelse(twisst_data$topos[[i]]$tip.label == "SanCristobal", 'SC',
                                                                                           ifelse(twisst_data$topos[[i]]$tip.label == "Genovesa", 'Ge',
                                                                                                  ifelse(twisst_data$topos[[i]]$tip.label == "Wolf", 'W', 'Other')))))))))
}
  
#################### subset to only the most abundant topologies #################
sum(twisst_data$weights_overall_mean)
#get list of the most abundant topologies (top 10 in this case)
top10_topos <- order(twisst_data$weights_overall_mean, decreasing=T)[1:10]

#subset twisst object for these
twisst_data_top10topos <- subset.twisst.by.topos(twisst_data, top10_topos)
#this can then be used in all the plotting functions.

#Get the weighting values for the top 10 topos
top10weights <- 100*twisst_data_top10topos$weights_overall_mean
top10weights

############################## combined plots ##################################
# there are a functions available to plot both the weightings and the topologies
#a summary plot shows all the topologies and a bar plot of their relative weightings
plot.twisst.summary(twisst_data, lwd=3, cex=0.7)
plot.twisst.summary.boxplot(twisst_data, lwd=3, cex=0.7)

#Same with top 10 topos
tiff("~/galmocks/Top10Topos.tiff", width = 2200, height = 1500, res = 300)
plot.twisst.summary(twisst_data_top10topos, lwd=3, cex=0.7)
dev.off()

#Make a histogram of the weightings
#First, convert the weights to a df
weightings <- as.data.frame(twisst_data$weights_overall_mean)
#Then rename the column to "Weights"
colnames(weightings)[1] = "Weights"
#Plot with ggplot
WeightsHist<-ggplot(weightings, aes(x = Weights))+
  geom_histogram(binwidth = .000001, color = "black")+
  theme_bw()+
  labs(x = "Topology Weight", y = "Count")
#Get the number of occurrences of each topo

window<-read.delim(window_data_file, sep = "\t")
weightings<-mutate(weightings, occurrences = round(weightings$Weights*length(window$scaffold)))
ggplot(weightings, aes(x =occurrences))+
  geom_histogram(binwidth = 1, color = "black", fill = "grey")+
  theme_bw()+
  labs(x = "Topology Weight", y = "Count")

#Boxplot isn't really relevant here
plot.twisst.summary.boxplot(twisst_data_top10topos, lwd=3, cex=0.7)

#Put the plots together in cowplot
Topos <- ggdraw()+draw_image("~/galmocks/Top10Topos.tiff")
ToposAndWeights <- plot_grid(Topos, WeightsHist, labels = "AUTO", ncol = 1)
ToposAndWeights

tiff("~/galmocks/ToposWeights.tiff", res = 300, height = 2500, width = 2000)
plot_grid(Topos, WeightsHist, labels = "AUTO", ncol = 1)
dev.off()
