"""
Created on 12 Jun 2019

@author: 2138645C
"""
import os
import pandas as pd
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage


class FigureGen(object):

    r_function = '''
#### General Information ####
# This code is used to generate figures from codon usage statistics (RSCU, ENC, FOP, GC3). 
# Code developed by Antonia Chalka
# General Figure Settings (for reference):
# TODO ADD 
#Load/Install required packages
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
if (!require("reshape2")) {
  install.packages("reshape2")
  library(reshape2)
}
if (!require("ggplot2")) {
  install.packages("ggplot2")
  library(ggplot2)
}
if (!require("pheatmap")) {
  install.packages("pheatmap")
  library(pheatmap)
}
if (!require("grid")) {
  install.packages("grid")
  library(grid)
}

#Create custom color palletes for heatmaps. For each one, a value of 0 is set to white (#F8F8F8) 
rscu.palette <- c("#F8F8F8", colorRampPalette(c("white","yellow", "orange", "red"))(n=35))
enc.palette <- c("#F8F8F8", colorRampPalette(c("white","yellow", "orange", "red"))(n=44))
fop.palette <- c("#F8F8F8", colorRampPalette(c("white","yellow", "orange", "red"))(n=20))

Gene_RSCU_Clustermap <- function(data, outpath, ann_data)
{ 
  #Convert input data to matrix
  dat <- as.matrix(data)
  
  #Remove rows with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store values that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Call appropriate figure method 
  if(missing(ann_data)) { #Make non-annotated heatmap if there is no annotation data
    make_heatmap(dat, rscu.palette, 8, 3, c(seq(0,4,length.out=36)), paste(outpath, "figure.png", sep="_"))
    print('Making non-annotated')
  } else { #Make annotated heatmap if annotation data has been provided 
    make_ann_heatmap(dat, ann_data, rscu.palette,8, 3, c(seq(0,4,length.out=36)), paste(outpath, "figure.png", sep="_"))
    print('Making annotated')
  }
}

Strain_RSCU_Clustermap <- function(data, outpath, ann_data)
{
  #Convert input data to matrix
  dat <- as.matrix(data)
  
  #Remove rows with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store values that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data  
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Call appropriate figure method 
  if(missing(ann_data)) { #Make non-annotated heatmap if there is no annotation data
    make_heatmap(dat, rscu.palette, 8, 3, c(seq(0,4,length.out=36)), paste(outpath, "figure.png", sep="_"))
  } else { #Make annotated heatmap if annotation data has been provided 
    make_ann_heatmap(dat, ann_data, rscu.palette,8, 3, c(seq(0,4,length.out=36)), paste(outpath, "figure.png", sep="_"))
  }
}

ENC_Heatmap <- function(data, outpath, gene_ann_data, strain_ann_data)
{
  # Transpose ENC data
  dat <- acast(data, data[,1]~data[,2], value.var='ENC')
  
  #Remove rows and columns with over half of their values as nan.
  colremoved <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store column values that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  rowremoved <- dat[,colSums(is.na(dat)) >= nrow(dat)-nrow(dat)/2] #Store row that will be removed
  dat <- dat[,colSums(is.na(dat)) < nrow(dat)-nrow(dat)/2] #Remove values
  
  # removed <- rbind(colremoved, rowremoved) # Combined removed values
  
  #Export data    
  # write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Call appropriate figure method 
  if(missing(ann_data) && missing(strain_ann_data)) { #Make non-annotated heatmap if there is no annotation data
    make_heatmap(dat, enc.palette, 3, 3, c(seq(21,64,length.out=45)), paste(outpath, "figure.png", sep="_"))
  } else { #Make annotated heatmap if annotation data has been provided 
    png(outname, width = 300, height = 300, units='mm', res = 300)
    pheatmap(dat, 
             #Color Options (general palette, nan color, do not display border)
             color=figure_palette, na_col ='grey', border_color=NA,
             
             #Clustering Options
             clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
             clustering_method = "ward.D2",
             
             #Annotation (type, colors)
             annotation_col=gene_ann_data,annotation_row=strain_ann_data, gp = gpar(fill = "grey"),
             
             #Other (size, font size, min/max value(breaks))
             treeheight_row = 100, treeheight_col = 100, fontsize_col = colsize, fontsize_row = rowsize, 
             breaks = limits
    )
    dev.off()
  }
}

Gene_RSCU_Clustermap_Di <- function(data, outpath, ann_data)
{
  #Convert input data to matrix
  dat <- as.matrix(data)
  
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Call appropriate figure method 
  if(missing(ann_data)) { #Make non-annotated heatmap if there is no annotation data
    make_heatmap(dat, rscu.palette, 8, 3, c(seq(0,2,length.out=36)), paste(outpath, "figure.png", sep="_"))
  } else { #Make annotated heatmap if annotation data has been provided 
    make_ann_heatmap(dat, ann_data, rscu.palette,8, 3, c(seq(0,2,length.out=36)), paste(outpath, "figure.png", sep="_"))
  }
}

Gene_RSCU_Clustermap_Tri <- function(data, outpath, ann_data)
{
  #Convert input data to matrix
  dat <- as.matrix(data)
  
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Call appropriate figure method 
  if(missing(ann_data)) { #Make non-annotated heatmap if there is no annotation data
    make_heatmap(dat, rscu.palette, 8, 3, c(seq(0,3,length.out=36)), paste(outpath, "figure.png", sep="_"))
  } else { #Make annotated heatmap if annotation data has been provided 
    make_ann_heatmap(dat, ann_data, rscu.palette,8, 3, c(seq(0,3,length.out=36)), paste(outpath, "figure.png", sep="_"))
  }
}

Gene_RSCU_Clustermap_Tetra <- function(data, outpath, ann_data)
{
  #Convert input data to matrix
  dat <- as.matrix(data)
  
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Call appropriate figure method 
  if(missing(ann_data)) { #Make non-annotated heatmap if there is no annotation data
    make_heatmap(dat, rscu.palette, 8, 3, c(seq(0,4,length.out=36)), paste(outpath, "figure.png", sep="_"))
  } else { #Make annotated heatmap if annotation data has been provided 
    make_ann_heatmap(dat, ann_data, rscu.palette,8, 3, c(seq(0,4,length.out=36)), paste(outpath, "figure.png", sep="_"))
  }
}

Gene_RSCU_Clustermap_Hexa <- function(data, outpath, ann_data)
{
  #Convert input data to matrix
  dat <- as.matrix(data)
  
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Call appropriate figure method 
  if(missing(ann_data)) { #Make non-annotated heatmap if there is no annotation data
    make_heatmap(dat, rscu.palette, 8, 3, c(seq(0,6,length.out=36)), paste(outpath, "figure.png", sep="_"))
  } else { #Make annotated heatmap if annotation data has been provided 
    make_ann_heatmap(dat, ann_data, rscu.palette,8, 3, c(seq(0,6,length.out=36)), paste(outpath, "figure.png", sep="_"))
  }
}

Strain_RSCU_Clustermap_Di <- function(data, outpath, ann_data)
{
  #Convert input data to matrix
  dat <- as.matrix(data)
  
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Call appropriate figure method 
  if(missing(ann_data)) { #Make non-annotated heatmap if there is no annotation data
    make_heatmap(dat, rscu.palette, 8, 3, c(seq(0,2,length.out=36)), paste(outpath, "figure.png", sep="_"))
  } else { #Make annotated heatmap if annotation data has been provided 
    make_ann_heatmap(dat, ann_data, rscu.palette,8, 3, c(seq(0,2,length.out=36)), paste(outpath, "figure.png", sep="_"))
  }
}

Strain_RSCU_Clustermap_Tri <- function(data, outpath, ann_data)
{
  #Convert input data to matrix
  dat <- as.matrix(data)
  
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Call appropriate figure method 
  if(missing(ann_data)) { #Make non-annotated heatmap if there is no annotation data
    make_heatmap(dat, rscu.palette, 8, 3, c(seq(0,3,length.out=36)), paste(outpath, "figure.png", sep="_"))
  } else { #Make annotated heatmap if annotation data has been provided 
    make_ann_heatmap(dat, ann_data, rscu.palette,8, 3, c(seq(0,3,length.out=36)), paste(outpath, "figure.png", sep="_"))
  }
}

Strain_RSCU_Clustermap_Tetra <- function(data, outpath, ann_data)
{
  #Convert input data to matrix
  dat <- as.matrix(data)
  
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Call appropriate figure method 
  if(missing(ann_data)) { #Make non-annotated heatmap if there is no annotation data
    make_heatmap(dat, rscu.palette, 8, 3, c(seq(0,4,length.out=36)), paste(outpath, "figure.png", sep="_"))
  } else { #Make annotated heatmap if annotation data has been provided 
    make_ann_heatmap(dat, ann_data, rscu.palette,8, 3, c(seq(0,4,length.out=36)), paste(outpath, "figure.png", sep="_"))
  }
}

Strain_RSCU_Clustermap_Hexa <- function(data, outpath, ann_data)
{
  #Convert input data to matrix
  dat <- as.matrix(data)
  
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Call appropriate figure method 
  if(missing(ann_data)) { #Make non-annotated heatmap if there is no annotation data
    make_heatmap(dat, rscu.palette, 8, 3, c(seq(0,6,length.out=36)), paste(outpath, "figure.png", sep="_"))
  } else { #Make annotated heatmap if annotation data has been provided 
    make_ann_heatmap(dat, ann_data, rscu.palette,8, 3, c(seq(0,6,length.out=36)), paste(outpath, "figure.png", sep="_"))
  }
}

FOP_Gene_Clustermap <- function(data, outpath, ann_data)
{
  #Convert input data to matrix
  dat <- as.matrix(data)
  
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Call appropriate figure method 
  if(missing(ann_data)) { #Make non-annotated heatmap if there is no annotation data
    make_heatmap(dat, fop.palette, 8, 3, c(seq(0,1,length.out=21)), paste(outpath, "figure.png", sep="_"))
  } else { #Make annotated heatmap if annotation data has been provided 
    make_ann_heatmap(dat, ann_data, fop.palette,8, 3, c(seq(0,1,length.out=21)), paste(outpath, "figure.png", sep="_") )
  }
}

FOP_Strain_Clustermap <- function(data, outpath, ann_data)
{
  #Convert input data to matrix
  dat <- as.matrix(data)
  
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Call appropriate figure method 
  if(missing(ann_data)) { #Make non-annotated heatmap if there is no annotation data
    make_heatmap(dat, fop.palette, 8, 3, c(seq(0,1,length.out=21)), paste(outpath, "figure.png", sep="_"))
  } else { #Make annotated heatmap if annotation data has been provided 
    make_ann_heatmap(dat, ann_data, fop.palette,8, 3, c(seq(0,1,length.out=21)), paste(outpath, "figure.png", sep="_") )
  }
  
}

FOP_Ref_Clustermap  <- function(data, outpath, ann_data)
{
  #Convert input data to matrix
  dat <- as.matrix(data)
  
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Create figure and save as to pre-specified dir (outpath)
  png(paste(outpath, "figure.png", sep="_"), width = 300, height = 300, units='mm', res = 300)
  pheatmap(dat, 
           #Color Options (general palette, nan color, do not display border)
           na_col ='grey', border_color=NA,
           #Clustering Options
           clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
           clustering_method = "ward.D2",
           #Other (size, )
           treeheight_row = 100, treeheight_col = 100, fontsize_col = 8, fontsize_row = 3
  )
  dev.off()
}

FOP_SamplesGenes_Clustermap  <- function(data, outpath, ann_data)
{
  #Convert input data to matrix
  dat <- as.matrix(data)
  
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Call appropriate figure method 
  if(missing(ann_data)) { #Make non-annotated heatmap if there is no annotation data
    make_heatmap(dat, fop.palette, 1, 3, c(seq(0,1,length.out=21)), paste(outpath, "figure.png", sep="_"))
  } else { #Make annotated heatmap if annotation data has been provided 
    make_ann_heatmap(dat, ann_data, fop.palette,1, 3, c(seq(0,1,length.out=21)), paste(outpath, "figure.png", sep="_") )
  }
  
}

FOP_SamplesStrains_Clustermap  <- function(data, outpath,ann_data)
{
  #Convert input data to matrix
  dat <- as.matrix(data)
  
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  #Call appropriate figure method 
  if(missing(ann_data)) { #Make non-annotated heatmap if there is no annotation data
    make_heatmap(dat, fop.palette, 1, 3, c(seq(0,1,length.out=21)), paste(outpath, "figure.png", sep="_"))
  } else { #Make annotated heatmap if annotation data has been provided 
    make_ann_heatmap(dat, ann_data, fop.palette,1, 3, c(seq(0,1,length.out=21)), paste(outpath, "figure.png", sep="_") )
  }
}

ENC_GC3_Genes  <- function(dat, outpath)
{
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  # Set up expected values of ENC vs GC3
  encmodel <- data.frame(GC3=seq(0.0, 1, by=0.01)) # Set up gc3 column
  encmodel$ENC <- 2 + encmodel$GC3 + 29/(encmodel$GC3^2 + (1-encmodel$GC3)^2) # calculate enc based on gc3
  
  #Create figure and save as to pre-specified dir (outpath)
  png(paste(outpath, "figure.png", sep="_"), width = 200, height = 200, units='mm', res = 300)
  p <-ggplot(dat, aes(GC3, ENC)) + 
    geom_point(aes(color="Observed")) +
    geom_line(data=encmodel,aes(color="Expected")) +
    labs(color="Legend text") +
    geom_text(aes(label=rownames(dat)),hjust=0,vjust=0,size=1)+
    theme_bw()
  print(p)
  dev.off()
}

ENC_GC3_Strains  <- function(dat, outpath)
{
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  # Set up expected values of ENC vs GC3
  encmodel <- data.frame(GC3=seq(0.0, 1, by=0.01)) # Set up gc3 column
  encmodel$ENC <- 2 + encmodel$GC3 + 29/(encmodel$GC3^2 + (1 - encmodel$GC3)^2) # calculate enc based on gc3
  
  #Create figure and save as to pre-specified dir (outpath)
  png(paste(outpath, "figure.png", sep="_"), width = 200, height = 200, units='mm', res = 300)
  p<-ggplot(dat, aes(GC3, ENC)) + 
    geom_point(aes(color="Observed")) +
    geom_line(data=encmodel,aes(color="Expected")) +
    labs(color="Legend text") +
    geom_text(aes(label=rownames(dat)),hjust=0,vjust=0,size=1) +
    theme_bw()
  print(p)
  dev.off()
}

ENC_GC3_ALL  <- function(dat, outpath)
{
  #Remove rows and columns with over half of their values as nan.
  removed <- dat[rowSums(is.na(dat)) >= ncol(dat)-ncol(dat)/2,] #Store row that will be removed
  dat <- dat[rowSums(is.na(dat)) < ncol(dat)-ncol(dat)/2,] #Remove values
  
  #Export data    
  write.table(removed, file = paste(outpath, "nan","removed.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  write.table(dat, file = paste(outpath, "data.csv", sep="_"), quote = FALSE, sep = ",", col.names =NA)
  
  # Set up expected values of ENC vs GC3
  encmodel <- data.frame(GC3=seq(0.0, 1, by=0.01)) # Set up gc3 column
  encmodel$ENC <- 2 + encmodel$GC3 + 29/(encmodel$GC3^2 + (1 - encmodel$GC3)^2) # calculate enc based on gc3
  
  #Create figure and save as to pre-specified dir (outpath)
  png(paste(outpath, "figure.png", sep="_"), width = 200, height = 200, units='mm', res = 300)
  p<-ggplot(dat, aes(GC3, ENC)) + 
    geom_point(aes(color="Observed")) +
    geom_line(data=encmodel,aes(color="Expected")) +
    labs(color="Legend text") +
    geom_text(aes(label=paste(Strain_ID, Gene, sep='_')),hjust=0,vjust=0,size=1) +
    theme_bw()
  print(p)
  dev.off()
}

make_heatmap <- function(dat, figure_palette, colsize, rowsize, limits, outname) {
  png(outname, width = 300, height = 300, units='mm', res = 300)
  pheatmap(dat, 
           #Color Options (general palette, nan color, do not display border)
           color=figure_palette, na_col ='grey', border_color=NA,
           
           #Clustering Options
           clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
           clustering_method = "ward.D2",
           
           #Other (size, font size, min/max value(breaks))
           treeheight_row = 100, treeheight_col = 100, fontsize_col = colsize, fontsize_row = rowsize, 
           breaks = limits
  )
  dev.off()
}

make_ann_heatmap <- function(dat, ann_data, figure_palette,  colsize, rowsize, limits, outname) {
  print(ann_data)
  png(outname, width = 300, height = 300, units='mm', res = 300)
  pheatmap(dat, 
           #Color Options (general palette, nan color, do not display border)
           color=figure_palette, na_col ='grey', border_color=NA,
           
           #Clustering Options
           clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
           clustering_method = "ward.D2",
           
           #Annotation (type, colors)
           annotation_row=ann_data, gp = gpar(fill = "grey"),
           
           #Other (size, font size, min/max value(breaks))
           treeheight_row = 100, treeheight_col = 100, fontsize_col = colsize, fontsize_row = rowsize, 
           breaks = limits
  )
  dev.off()
}
    '''

    # Compile above R code into a pseudo-package, allowing its function to be called for generating figures
    RFunction = SignatureTranslatedAnonymousPackage(r_function, "RFunction")

    # Name of RSCU columns to be used for filtering. Does not include start, stop, Tryp and digenerate
    # (C,D,E,F,H,K,N,Q,Y) aa with T/A at 3rd position
    # Codons removed: C (UGU), D (GAU), E (GAA), F (UUU), H (CAU), K (AAA), START (AUG), N (AAU),
    # Q (CAA), W (UGG), Y (UAU), STOP (UAA UAG UGA)
    rscu_columns = ['S_UCU_RSCU',
                    'F_UUC_RSCU', 'S_UCC_RSCU', 'Y_UAC_RSCU', 'C_UGC_RSCU',
                    'L_UUA_RSCU', 'S_UCA_RSCU',
                    'L_UUG_RSCU', 'S_UCG_RSCU',
                    'L_CUU_RSCU', 'P_CCU_RSCU', 'R_CGU_RSCU',
                    'L_CUC_RSCU', 'P_CCC_RSCU', 'H_CAC_RSCU', 'R_CGC_RSCU',
                    'L_CUA_RSCU', 'P_CCA_RSCU', 'R_CGA_RSCU',
                    'L_CUG_RSCU', 'P_CCG_RSCU', 'Q_CAG_RSCU', 'R_CGG_RSCU',
                    'I_AUU_RSCU', 'T_ACU_RSCU', 'S_AGU_RSCU',
                    'I_AUC_RSCU', 'T_ACC_RSCU', 'N_AAC_RSCU', 'S_AGC_RSCU',
                    'I_AUA_RSCU', 'T_ACA_RSCU', 'R_AGA_RSCU',
                    'T_ACG_RSCU', 'K_AAG_RSCU', 'R_AGG_RSCU',
                    'V_GUU_RSCU', 'A_GCU_RSCU', 'G_GGU_RSCU',
                    'V_GUC_RSCU', 'A_GCC_RSCU', 'D_GAC_RSCU', 'G_GGC_RSCU',
                    'V_GUA_RSCU', 'A_GCA_RSCU', 'G_GGA_RSCU',
                    'V_GUG_RSCU', 'A_GCG_RSCU', 'E_GAG_RSCU', 'G_GGG_RSCU']
    di_aa = ['C_UGC_RSCU',
             'D_GAC_RSCU',
             'E_GAG_RSCU',
             'F_UUC_RSCU',
             'H_CAC_RSCU',
             'K_AAG_RSCU',
             'N_AAC_RSCU',
             'Q_CAG_RSCU',
             'Y_UAC_RSCU']
    tri_aa = ['I_AUU_RSCU', 'I_AUC_RSCU', 'I_AUA_RSCU']
    tetra_aa = ['A_GCU_RSCU', 'A_GCC_RSCU', 'A_GCA_RSCU', 'A_GCG_RSCU',
                'G_GGU_RSCU', 'G_GGC_RSCU', 'G_GGA_RSCU', 'G_GGG_RSCU',
                'P_CCU_RSCU', 'P_CCC_RSCU', 'P_CCA_RSCU', 'P_CCG_RSCU',
                'T_ACU_RSCU', 'T_ACC_RSCU', 'T_ACA_RSCU', 'T_ACG_RSCU',
                'V_GUU_RSCU', 'V_GUC_RSCU', 'V_GUA_RSCU', 'V_GUG_RSCU']
    hexa_aa = ['L_UUA_RSCU', 'L_UUG_RSCU', 'L_CUU_RSCU', 'L_CUC_RSCU', 'L_CUA_RSCU', 'L_CUG_RSCU',
               'R_CGU_RSCU', 'R_CGC_RSCU', 'R_CGA_RSCU', 'R_CGG_RSCU', 'R_AGA_RSCU', 'R_AGG_RSCU',
               'S_UCU_RSCU', 'S_UCC_RSCU', 'S_UCA_RSCU', 'S_UCG_RSCU', 'S_AGU_RSCU', 'S_AGC_RSCU']

    def __init__(self, data, gene_annotation_path=None, strain_annotation_path=None, regex="Possible|note|RL5A_+|RL6_+"):
        self.masterdf = data

        if gene_annotation_path is not None:
            gene_ann = pd.read_csv(gene_annotation_path)
            gene_ann.set_index(gene_ann.columns[0], inplace=True)
            self.gene_ann = gene_ann
            print(self.gene_ann)
        else:
            self.gene_ann = None

        if strain_annotation_path is not None:
            strain_ann = pd.read_csv(strain_annotation_path)
            strain_ann.set_index(strain_ann.column[0], inplace=True)
            self.strain_ann = strain_ann
        else:
            self.strain_ann = None

        self.filter_regex = regex
        pandas2ri.activate()

    def make_graph(self, graph_type, out_path=os.path.dirname(os.path.realpath(__file__)) + '//Output//'):
        # Regex filter
        print(graph_type)
        if graph_type != 'FOP_Ref_Clustermap':
            filter_rows = self.masterdf['Gene'].str.contains(self.filter_regex)
            removed = self.masterdf[filter_rows]
            removed.to_csv(out_path + graph_type + "_regex_removed.csv")
            self.masterdf = self.masterdf[~filter_rows]

        outvar = out_path + graph_type
        if graph_type == 'Gene_RSCU_Clustermap':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Gene')[self.rscu_columns].mean()
            graph_values.to_csv(outvar, index=True)

            # If annotation file has been loaded, add it when calling R function
            if self.gene_ann is None:
                self.RFunction.Gene_RSCU_Clustermap(graph_values, outvar)
            elif self.gene_ann is not None:
                self.RFunction.Gene_RSCU_Clustermap(graph_values, outvar, self.gene_ann)

        elif graph_type == 'Strain_RSCU_Clustermap':

            # group by Strain, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Strain_ID')[self.rscu_columns].mean()
            graph_values.to_csv(outvar, index=True)

            # Check number of strains in data, because clustering algorithm requires 2 rows/columns
            if len(graph_values.index) > 1:
                # If annotation file has been loaded, add it when calling R function
                if self.strain_ann is None:
                    self.RFunction.Strain_RSCU_Clustermap(graph_values, outvar)
                elif self.strain_ann is not None:
                    self.RFunction.Strain_RSCU_Clustermap(graph_values, outvar, self.strain_ann)

        elif graph_type == 'ENC_Heatmap':

            graph_values = self.masterdf[['Strain_ID', 'Gene', 'ENC']]  # Obtain relevant columns
            graph_values.to_csv(outvar, index=True)

            # Check number of strains + genes (rows and columns), because heatmap.2 requires 2 rows/columns
            if graph_values.Strain_ID.nunique() > 1 and graph_values.Gene.nunique() > 1:
                # 'Melt' ENC values (transpose into a matrix), with Strains as rows and Genes as columns
                # graph_values = graph_values.pivot_table(index='Strain_ID', columns='Gene', values='ENC')

                # If annotation files have been loaded, add it when calling R function
                if (self.strain_ann is None) and (self.gene_ann is None):
                    self.RFunction.ENC_Heatmap(graph_values, outvar)
                elif (self.strain_ann is not None) and (self.gene_ann is not None):
                    self.RFunction.ENC_Heatmap(graph_values, outvar, self.gene_ann, self.strain_ann)

        elif graph_type == 'Gene_RSCU_Clustermap_Di':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Gene')[self.di_aa].mean()
            graph_values.to_csv(outvar, index=True)

            # If annotation file has been loaded, add it when calling R function
            if self.gene_ann is None:
                self.RFunction.Gene_RSCU_Clustermap_Di(graph_values, outvar)
            elif self.gene_ann is not None:
                self.RFunction.Gene_RSCU_Clustermap_Di(graph_values, outvar, self.gene_ann)

        elif graph_type == 'Gene_RSCU_Clustermap_Tri':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Gene')[self.tri_aa].mean()
            graph_values.to_csv(outvar, index=True)

            # If annotation file has been loaded, add it when calling R function
            if self.gene_ann is None:
                self.RFunction.Gene_RSCU_Clustermap_Tri(graph_values, outvar)
            elif self.gene_ann is not None:
                self.RFunction.Gene_RSCU_Clustermap_Tri(graph_values, outvar, self.gene_ann)

        elif graph_type == 'Gene_RSCU_Clustermap_Tetra':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Gene')[self.tetra_aa].mean()
            graph_values.to_csv(outvar, index=True)

            # If annotation file has been loaded, add it when calling R function
            if self.gene_ann is None:
                self.RFunction.Gene_RSCU_Clustermap_Tetra(graph_values, outvar)
            elif self.gene_ann is not None:
                self.RFunction.Gene_RSCU_Clustermap_Tetra(graph_values, outvar, self.gene_ann)

        elif graph_type == 'Gene_RSCU_Clustermap_Hexa':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Gene')[self.hexa_aa].mean()
            graph_values.to_csv(outvar, index=True)

            # If annotation file has been loaded, add it when calling R function
            if self.gene_ann is None:
                self.RFunction.Gene_RSCU_Clustermap_Hexa(graph_values, outvar)
            elif self.gene_ann is not None:
                self.RFunction.Gene_RSCU_Clustermap_Hexa(graph_values, outvar, self.gene_ann)

        elif graph_type == 'Strain_RSCU_Clustermap_Di':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Strain_ID')[self.di_aa].mean()
            graph_values.to_csv(outvar, index=True)

            if len(graph_values.index) > 1:
                # If annotation file has been loaded, add it when calling R function
                if self.strain_ann is None:
                    self.RFunction.Strain_RSCU_Clustermap_Di(graph_values, outvar)
                elif self.strain_ann is not None:
                    self.RFunction.Strain_RSCU_Clustermap_Di(graph_values, outvar, self.strain_ann)

        elif graph_type == 'Strain_RSCU_Clustermap_Tri':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Strain_ID')[self.tri_aa].mean()
            graph_values.to_csv(outvar, index=True)

            if len(graph_values.index) > 1:
                # If annotation file has been loaded, add it when calling R function
                if self.strain_ann is None:
                    self.RFunction.Strain_RSCU_Clustermap_Tri(graph_values, outvar)
                elif self.strain_ann is not None:
                    self.RFunction.Strain_RSCU_Clustermap_Tri(graph_values, outvar, self.strain_ann)

        elif graph_type == 'Strain_RSCU_Clustermap_Tetra':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Strain_ID')[self.tetra_aa].mean()
            graph_values.to_csv(outvar, index=True)

            if len(graph_values.index) > 1:
                # If annotation file has been loaded, add it when calling R function
                if self.strain_ann is None:
                    self.RFunction.Strain_RSCU_Clustermap_Tetra(graph_values, outvar)
                elif self.strain_ann is not None:
                    self.RFunction.Strain_RSCU_Clustermap_Tetra(graph_values, outvar, self.strain_ann)

        elif graph_type == 'Strain_RSCU_Clustermap_Hexa':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Strain_ID')[self.hexa_aa].mean()
            graph_values.to_csv(outvar, index=True)

            if len(graph_values.index) != 1:
                # If annotation file has been loaded, add it when calling R function
                if self.strain_ann is None:
                    self.RFunction.Strain_RSCU_Clustermap_Hexa(graph_values, outvar)
                elif self.strain_ann is not None:
                    self.RFunction.Strain_RSCU_Clustermap_Hexa(graph_values, outvar, self.strain_ann)

        elif graph_type == 'FOP_Gene_Clustermap':
            fop_columns = self.masterdf.filter(regex="^FOP_*").columns.tolist()
            graph_values = self.masterdf.groupby('Gene')[fop_columns].mean()
            graph_values.to_csv(outvar, index=True)

            # If annotation file has been loaded, add it when calling R function
            if self.gene_ann is None:
                self.RFunction.FOP_Gene_Clustermap(graph_values, outvar)
            elif self.gene_ann is not None:
                self.RFunction.FOP_Gene_Clustermap(graph_values, outvar, self.gene_ann)

        elif graph_type == 'FOP_Strain_Clustermap':
            fop_columns = self.masterdf.filter(regex="^FOP_*").columns.tolist()
            graph_values = self.masterdf.groupby('Strain_ID')[fop_columns].mean()
            graph_values.to_csv(outvar, index=True)

            if len(graph_values.index) != 1:
                # If annotation file has been loaded, add it when calling R function
                if self.strain_ann is None:
                    self.RFunction.FOP_Strain_Clustermap(graph_values, outvar)
                elif self.strain_ann is not None:
                    self.RFunction.FOP_Strain_Clustermap(graph_values, outvar, self.strain_ann)

        elif graph_type == 'FOP_Ref_Clustermap':
            graph_values = self.masterdf
            graph_values = graph_values.select_dtypes(include='number')

            self.RFunction.FOP_Ref_Clustermap(graph_values, outvar)

        elif graph_type == 'FOP_SamplesGenes_Clustermap':
            fop_columns = self.masterdf.filter(regex="^FOP_*").columns.tolist()
            graph_values = self.masterdf.groupby('Gene')[fop_columns].mean()
            graph_values.to_csv(outvar, index=True)

            # If annotation file has been loaded, add it when calling R function
            if self.gene_ann is None:
                self.RFunction.FOP_SamplesGenes_Clustermap(graph_values, outvar)
            elif self.gene_ann is not None:
                self.RFunction.FOP_SamplesGenes_Clustermap(graph_values, outvar, self.gene_ann)

        elif graph_type == 'FOP_SamplesStrains_Clustermap':
            fop_columns = self.masterdf.filter(regex="^FOP_*").columns.tolist()
            graph_values = self.masterdf.groupby('Strain_ID')[fop_columns].mean()
            graph_values.to_csv(outvar, index=True)

            if len(graph_values.index) != 1:
                # If annotation file has been loaded, add it when calling R function
                if self.strain_ann is None:
                    self.RFunction.FOP_SamplesStrains_Clustermap(graph_values, outvar)
                elif self.strain_ann is not None:
                    self.RFunction.FOP_SamplesStrains_Clustermap(graph_values, outvar, self.strain_ann)

        elif graph_type == 'ENC_GC3_Genes':
            graph_values = self.masterdf.groupby('Gene')[['ENC', 'GC3']].mean()
            graph_values.to_csv(outvar, index=True)

            self.RFunction.ENC_GC3_Genes(graph_values, outvar)

        elif graph_type == 'ENC_GC3_Strains':
            graph_values = self.masterdf.groupby('Strain_ID')[['ENC', 'GC3']].mean()
            graph_values.to_csv(outvar, index=True)

            if len(graph_values.index) != 1:
                self.RFunction.ENC_GC3_Strains(graph_values, outvar)

        elif graph_type == 'ENC_GC3_ALL':
            graph_values = self.masterdf[['Gene', 'Strain_ID', 'ENC', 'GC3']]
            graph_values.to_csv(outvar, index=True)

            self.RFunction.ENC_GC3_ALL(graph_values, outvar)


'''
# Make fop ref clustermap
fop_ref = pd.read_csv("fop_ref.csv")
fop_ref.set_index('Tissue', inplace=True)
fop_ref.drop(['SMTS','SMTSD'], axis=1, inplace=True)
print(fop_ref)
new = FigureGen(fop_ref)
new.make_graph('FOP_Ref_Clustermap')

new = FigureGen(pd.read_csv("C:\\Users\\anni1\\PycharmProjects\\MScProject\\CodonUsage_Output\\fopmode_masterfile.csv"),
                pd.read_csv("C:\\Users\\anni1\\PycharmProjects\\MScProject\\CodonUsage\\ref_files\\time_class_name.csv",
                 index_col='Gene')
                )
new.make_graph('Gene_RSCU_Clustermap')
new.make_graph('Strain_RSCU_Clustermap')

new.make_graph('ENC_Heatmap')

new.make_graph('Gene_RSCU_Clustermap_Di')
new.make_graph('Gene_RSCU_Clustermap_Tri')
new.make_graph('Gene_RSCU_Clustermap_Tetra')
new.make_graph('Gene_RSCU_Clustermap_Hexa')
new.make_graph('Strain_RSCU_Clustermap_Di')
new.make_graph('Strain_RSCU_Clustermap_Tri')
new.make_graph('Strain_RSCU_Clustermap_Tetra')
new.make_graph('Strain_RSCU_Clustermap_Hexa')

new.make_graph('FOP_Gene_Clustermap')
new.make_graph('FOP_Strain_Clustermap')

new.make_graph('FOP_SamplesStrains_Clustermap')
new.make_graph('FOP_SamplesGenes_Clustermap')

new.make_graph('ENC_GC3_ALL')
new.make_graph('ENC_GC3_Genes')
new.make_graph('ENC_GC3_Strains')
'''
