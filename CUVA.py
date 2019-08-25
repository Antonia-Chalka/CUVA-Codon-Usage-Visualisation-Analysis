"""
Created on 20 May 2019
@author: Antonia Chalka
"""
import pandas as pd
import os
import numpy as np
import subprocess as sb
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
import PySimpleGUI as sg


class CodonData(object):
    """
    Class representing codon data across strains and genes. Relies on codonW to perform RSCU and ENC analysis.
    All data is stored and manipulated in masterdf, a pandas data frame.
    """
    # List of codon names, used for filtering RSCU input file
    CodonNames = ['UUU', 'UCU', 'UAU', 'UGU',
                  'UUC', 'UCC', 'UAC', 'UGC',
                  'UUA', 'UCA', 'UAA', 'UGA',
                  'UUG', 'UCG', 'UAG', 'UGG',
                  'CUU', 'CCU', 'CAU', 'CGU',
                  'CUC', 'CCC', 'CAC', 'CGC',
                  'CUA', 'CCA', 'CAA', 'CGA',
                  'CUG', 'CCG', 'CAG', 'CGG',
                  'AUU', 'ACU', 'AAU', 'AGU',
                  'AUC', 'ACC', 'AAC', 'AGC',
                  'AUA', 'ACA', 'AAA', 'AGA',
                  'AUG', 'ACG', 'AAG', 'AGG',
                  'GUU', 'GCU', 'GAU', 'GGU',
                  'GUC', 'GCC', 'GAC', 'GGC',
                  'GUA', 'GCA', 'GAA', 'GGA',
                  'GUG', 'GCG', 'GAG', 'GGG']
    # List of aa names, used for filtering RSCU input file
    AANames = ['Ala', 'Arg', 'Asn', 'Asp',
               'Cys', 'Glu', 'Gln', 'Gly',
               'His', 'Ile', 'Leu', 'Lys',
               'Met', 'Phe', 'Pro', 'Ser',
               'Thr', 'Trp', 'Tyr', 'Val', 'TER']
    # Setting up columns of masterdf
    col_names = ['Strain_ID', 'Gene', 'Total_Num_Codons',
                 'F_UUU_Num', 'F_UUU_RSCU',
                 'S_UCU_Num', 'S_UCU_RSCU',
                 'Y_UAU_Num', 'Y_UAU_RSCU',
                 'C_UGU_Num', 'C_UGU_RSCU',
                 'F_UUC_Num', 'F_UUC_RSCU',
                 'S_UCC_Num', 'S_UCC_RSCU',
                 'Y_UAC_Num', 'Y_UAC_RSCU',
                 'C_UGC_Num', 'C_UGC_RSCU',
                 'L_UUA_Num', 'L_UUA_RSCU',
                 'S_UCA_Num', 'S_UCA_RSCU',
                 '*_UAA_Num', '*_UAA_RSCU',
                 '*_UGA_Num', '*_UGA_RSCU',
                 'L_UUG_Num', 'L_UUG_RSCU',
                 'S_UCG_Num', 'S_UCG_RSCU',
                 '*_UAG_Num', '*_UAG_RSCU',
                 'W_UGG_Num', 'W_UGG_RSCU',
                 'L_CUU_Num', 'L_CUU_RSCU',
                 'P_CCU_Num', 'P_CCU_RSCU',
                 'H_CAU_Num', 'H_CAU_RSCU',
                 'R_CGU_Num', 'R_CGU_RSCU',
                 'L_CUC_Num', 'L_CUC_RSCU',
                 'P_CCC_Num', 'P_CCC_RSCU',
                 'H_CAC_Num', 'H_CAC_RSCU',
                 'R_CGC_Num', 'R_CGC_RSCU',
                 'L_CUA_Num', 'L_CUA_RSCU',
                 'P_CCA_Num', 'P_CCA_RSCU',
                 'Q_CAA_Num', 'Q_CAA_RSCU',
                 'R_CGA_Num', 'R_CGA_RSCU',
                 'L_CUG_Num', 'L_CUG_RSCU',
                 'P_CCG_Num', 'P_CCG_RSCU',
                 'Q_CAG_Num', 'Q_CAG_RSCU',
                 'R_CGG_Num', 'R_CGG_RSCU',
                 'I_AUU_Num', 'I_AUU_RSCU',
                 'T_ACU_Num', 'T_ACU_RSCU',
                 'N_AAU_Num', 'N_AAU_RSCU',
                 'S_AGU_Num', 'S_AGU_RSCU',
                 'I_AUC_Num', 'I_AUC_RSCU',
                 'T_ACC_Num', 'T_ACC_RSCU',
                 'N_AAC_Num', 'N_AAC_RSCU',
                 'S_AGC_Num', 'S_AGC_RSCU',
                 'I_AUA_Num', 'I_AUA_RSCU',
                 'T_ACA_Num', 'T_ACA_RSCU',
                 'K_AAA_Num', 'K_AAA_RSCU',
                 'R_AGA_Num', 'R_AGA_RSCU',
                 'M_AUG_Num', 'M_AUG_RSCU',
                 'T_ACG_Num', 'T_ACG_RSCU',
                 'K_AAG_Num', 'K_AAG_RSCU',
                 'R_AGG_Num', 'R_AGG_RSCU',
                 'V_GUU_Num', 'V_GUU_RSCU',
                 'A_GCU_Num', 'A_GCU_RSCU',
                 'D_GAU_Num', 'D_GAU_RSCU',
                 'G_GGU_Num', 'G_GGU_RSCU',
                 'V_GUC_Num', 'V_GUC_RSCU',
                 'A_GCC_Num', 'A_GCC_RSCU',
                 'D_GAC_Num', 'D_GAC_RSCU',
                 'G_GGC_Num', 'G_GGC_RSCU',
                 'V_GUA_Num', 'V_GUA_RSCU',
                 'A_GCA_Num', 'A_GCA_RSCU',
                 'E_GAA_Num', 'E_GAA_RSCU',
                 'G_GGA_Num', 'G_GGA_RSCU',
                 'V_GUG_Num', 'V_GUG_RSCU',
                 'A_GCG_Num', 'A_GCG_RSCU',
                 'E_GAG_Num', 'E_GAG_RSCU',
                 'G_GGG_Num', 'G_GGG_RSCU',
                 'ENC', 'GC3']
    # Used for rscu quality controls (qc_rscu)
    SynonymousCodons = {
        'Ala': ['A_GCA', 'A_GCC', 'A_GCG', 'A_GCU'],
        'Cys': ['C_UGU', 'C_UGC'],
        'Asp': ['D_GAU', 'D_GAC'],
        'Glu': ['E_GAG', 'E_GAA'],
        'Phe': ['F_UUU', 'F_UUC'],
        'Gly': ['G_GGU', 'G_GGG', 'G_GGA', 'G_GGC'],
        'His': ['H_CAU', 'H_CAC'],
        'Ile': ['I_AUC', 'I_AUA', 'I_AUU'],
        'Lys': ['K_AAG', 'K_AAA'],
        'Leu': ['L_UUA', 'L_UUG', 'L_CUC', 'L_CUU', 'L_CUG', 'L_CUA'],
        'Met': ['M_AUG'],
        'Asn': ['N_AAC', 'N_AAU'],
        'Pro': ['P_CCU', 'P_CCG', 'P_CCA', 'P_CCC'],
        'Gln': ['Q_CAA', 'Q_CAG'],
        'Arg': ['R_CGA', 'R_CGC', 'R_CGG', 'R_CGU', 'R_AGG', 'R_AGA'],
        'Ser': ['S_UCU', 'S_UCG', 'S_UCA', 'S_UCC', 'S_AGC', 'S_AGU'],
        'Thr': ['T_ACC', 'T_ACA', 'T_ACG', 'T_ACU'],
        'Val': ['V_GUA', 'V_GUC', 'V_GUG', 'V_GUU'],
        'Trp': ['W_UGG'],
        'Tyr': ['Y_UAU', 'Y_UAC']
    }

    def __init__(self, fop_type=None, fop_ref_df=None, fop_dic=None,
                 codonw_dir=os.path.dirname(os.path.realpath(__file__)) + "\\ref_files\\codonW",
                 codonw_out_dir=os.path.dirname(os.path.realpath(__file__)) + "\\codonw_out\\"):
        print("Initialising model data frame...")
        self.masterdf = pd.DataFrame(columns=self.col_names)  # Set up masterdf
        self.fop_type = fop_type
        self.fop_ref_df = fop_ref_df
        self.fop_dic = fop_dic
        self.codonw_dir = codonw_dir
        self.codonw_out_dir = codonw_out_dir

    def make_rscu_file(self, input_file, input_dir):
        print("Calculating RSCU values for file:" + input_file)

        # command summary = codonw inputfile -silent -nomenu -rscuoutname1 outname2
        cmd = self.codonw_dir + " " + input_dir + input_file + ".fasta" + \
            ' -silent -nomenu -rscu ' + \
            self.codonw_out_dir + input_file + '.none ' + \
            self.codonw_out_dir + input_file + '.blk'
        print("Running CodonW command: " + cmd)
        out = sb.run(cmd, stdout=sb.PIPE).stdout.decode('utf-8')

        # print("Parsing RSCU file...")
        rscu_file = open(self.codonw_out_dir + input_file + '.blk', "r")
        sections = rscu_file.read().split("\n\n")  # split when 2 newlines (1 gene has 5 of such sections)
        sections = [x.split() for x in sections]  # remove whitespace
        # remove Amino Acid and Codon Names
        sections = [[i for i in x if i not in self.AANames and i not in self.CodonNames] for x in sections]
        rscu_file.close()

        return sections, out

    def make_enc_file(self, input_file, input_dir):
        print("Calculating ENC values for file:" + input_file)
        # command summary = codonw inputfile (options) rscuoutname1 outname2
        cmd = self.codonw_dir + " " + input_dir + input_file + ".fasta" + \
            ' -silent -nomenu -enc -machine ' + \
            self.codonw_out_dir + input_file + '.enc ' + \
            self.codonw_out_dir + input_file + '.blk'
        print("Running CodonW command: " + cmd)
        out = sb.run(cmd, stdout=sb.PIPE).stdout.decode('utf-8')

        print("Parsing ENC file...")
        enc_file = open(self.codonw_out_dir + input_file + ".enc", "r")
        enc_rows = enc_file.read().split("\n")  # Extract rows
        enc_file.close()

        return enc_rows, out

    def make_gc3_file(self, input_file, input_dir):
        print("Calculating GC3 values for file: " + input_file)
        # command summary = codonw inputfile (options) rscuoutname1 outname2
        cmd = self.codonw_dir + " " + input_dir + input_file + ".fasta" + \
            ' -silent -nomenu -gc3s -machine ' + \
            self.codonw_out_dir + input_file + '.gc3 ' + \
            self.codonw_out_dir + input_file + '.blk'
        print("Running CodonW command: " + cmd)
        out = sb.run(cmd, stdout=sb.PIPE).stdout.decode('utf-8')

        print("Parsing GC3 file...")
        gc3_file = open(self.codonw_out_dir + input_file + ".gc3", "r")
        gc3_rows = gc3_file.read().split("\n")  # Extract rows
        gc3_file.close()

        return gc3_rows, out

    def parse_masterfile_values(self, strain_name, rscu_sections, enc_rows, gc3_rows):
        print("Gathering values to add to masterdf for Strain: " + strain_name + "...")

        # Populate master file
        print("Gathering RSCU values...")
        rscu_values = []
        rscu_gene_list = []

        # Each loop represents a gene and is used to create a row to be added to data frame
        for i in range(0, len(rscu_sections) - 1, 5):
            # Collate RSCU Sections
            rcsu_row = [int(rscu_sections[i + 4][0])]  # Number of codons
            rcsu_row.extend((float(j) for j in rscu_sections[i]))
            rcsu_row.extend((float(j) for j in rscu_sections[i + 1]))
            rcsu_row.extend((float(j) for j in rscu_sections[i + 2]))
            rcsu_row.extend((float(j) for j in rscu_sections[i + 3]))

            rscu_values.append(rcsu_row)
            rscu_gene_list.append(rscu_sections[i + 4][3])

        print("Parsing ENC values...")
        enc_values = []
        enc_gene_list = []
        for i in range(1, len(enc_rows) - 1):
            enc_fields = enc_rows[i].split()
            enc_gene_list.append(enc_fields[0])
            if enc_fields[1] == '*****':
                enc_fields[1] = np.nan
            enc_values.append(float(enc_fields[1]))

        print("Parsing GC3 values...")
        gc3_values = []
        gc3_gene_list = []
        for i in range(1, len(gc3_rows) - 1):
            gc3_fields = gc3_rows[i].split()
            gc3_gene_list.append(gc3_fields[0])
            if len(gc3_fields) == 1:
                gc3_fields.append(np.nan)
            gc3_values.append(float(gc3_fields[1]))

        print("Collating values in masterdf...")
        for i in range(0, len(rscu_values)):
            masterfile_row = [strain_name, enc_gene_list[i]]
            masterfile_row.extend(rscu_values[i])
            masterfile_row.append(enc_values[i])
            masterfile_row.append(gc3_values[i])

            self.masterdf.loc[len(self.masterdf)] = masterfile_row

    def iterate_populate_master_file(self, file_path, infile):
        input_file = os.path.splitext(infile)[0]

        rscu_sections, rscu_out = self.make_rscu_file(input_file, file_path)
        enc_rows, enc_out = self.make_enc_file(input_file, file_path)
        gc3_rows, gc3_out = self.make_gc3_file(input_file, file_path)

        self.parse_masterfile_values(input_file, rscu_sections, enc_rows, gc3_rows)
        return rscu_out + "\n" + enc_out + "\n" + gc3_out + "\n"

    def qc_rscu(self, rscuqc_thres=15):  # Quality control RSCU  values
        # rscuqc_thres: Threshold for how many times a codon column must have a non-NaN value.
        # Above threshold counts are retained and vice versa.

        print("Beginning quality control of RSCU values across amino acids...")
        for i in self.masterdf.index:
            # print("Checking row " + str(i))
            for aa in self.SynonymousCodons:
                degeneracy = len(self.SynonymousCodons[aa])
                aa_count = 0

                for codon in self.SynonymousCodons[aa]:
                    aa_count += self.masterdf.at[i, codon + "_Num"]
                # print("AA: " + aa + ". Codon (column): " + codon +
                #      ". Degeneracy: " + str(degeneracy) + ". Count: " + str(aa_count))
                if degeneracy > aa_count:
                    # print("Num of codons less than degeneracy. Filling with Nan values...")
                    for codon in self.SynonymousCodons[aa]:
                        self.masterdf.at[i, codon + "_RSCU"] = np.nan

        print("Beginning quality control of RSCU values across strains...")
        grouped = self.masterdf.groupby("Strain_ID")
        for strain, group in grouped:
            # print("Checking codons for Strain:" + strain)
            for aa in self.SynonymousCodons:
                for codon in self.SynonymousCodons[aa]:
                    rscu_column = codon + "_RSCU"
                    if group[rscu_column].notna().count() < rscuqc_thres:
                        # print("Codon with > " + str(rscuqc_thres) + " reps ( " + strain + " " + codon +
                        # " ). Setting to Nan...")
                        self.masterdf.loc[self.masterdf.Strain_ID == strain, rscu_column] = np.nan

    def set_fop_ref(self, fop_ref_path):
        print("reading fop file...")
        self.fop_ref_df = pd.read_csv(fop_ref_path)

    def set_fop_dict(self, fop_type, type_column_name='SMTS', tissue_column_name='Tissue'):
        self.fop_type = fop_type
        print("setting reference..." + fop_type)
        columns = ['A_GCA', 'A_GCC', 'A_GCG', 'A_GCU', 'C_UGC', 'C_UGU', 'E_GAA', 'E_GAG', 'D_GAC', 'D_GAU',
                   'G_GGA', 'G_GGC', 'G_GGG', 'G_GGU', 'F_UUC', 'F_UUU', 'I_AUA', 'I_AUC', 'I_AUU', 'H_CAC',
                   'H_CAU', 'K_AAA', 'K_AAG', 'L_CUA', 'L_CUC', 'L_CUG', 'L_CUU', 'L_UUA', 'L_UUG', 'N_AAC',
                   'N_AAU', 'Q_CAA', 'Q_CAG', 'P_CCA', 'P_CCC', 'P_CCG', 'P_CCU', 'S_AGC', 'S_AGU', 'S_UCA',
                   'S_UCC', 'S_UCG', 'S_UCU', 'R_AGA', 'R_AGG', 'R_CGA', 'R_CGC', 'R_CGG', 'R_CGU', 'T_ACA',
                   'T_ACC', 'T_ACG', 'T_ACU', 'V_GUA', 'V_GUC', 'V_GUG', 'V_GUU', 'Y_UAC', 'Y_UAU']
        if fop_type == 'mode':
            self.fop_dic = self.fop_ref_df.groupby(type_column_name)[columns].agg(pd.Series.mode)
        elif fop_type == 'random':
            self.fop_dic = self.fop_ref_df.groupby(type_column_name)[columns].first()
        elif fop_type == 'all':
            # Rename Tissue column to include it's type (e.g. Whole Blood, EBV infected etc)
            self.fop_dic = self.fop_ref_df.assign(Tissue=self.fop_ref_df[tissue_column_name]
                                                  + "_" + self.fop_ref_df[type_column_name])
            columns.insert(0, tissue_column_name)
            self.fop_dic = self.fop_dic[columns]  # Filter columns

            self.fop_dic.set_index('Tissue', drop=True, inplace=True)  # Set Tissue as index

    def calculate_fop(self, tissues):
        codons = ['A_GCA', 'A_GCC', 'A_GCG', 'A_GCU', 'C_UGC', 'C_UGU', 'E_GAA', 'E_GAG', 'D_GAC', 'D_GAU',
                  'G_GGA', 'G_GGC', 'G_GGG', 'G_GGU', 'F_UUC', 'F_UUU', 'I_AUA', 'I_AUC', 'I_AUU', 'H_CAC',
                  'H_CAU', 'K_AAA', 'K_AAG', 'L_CUA', 'L_CUC', 'L_CUG', 'L_CUU', 'L_UUA', 'L_UUG', 'N_AAC',
                  'N_AAU', 'Q_CAA', 'Q_CAG', 'P_CCA', 'P_CCC', 'P_CCG', 'P_CCU', 'S_AGC', 'S_AGU', 'S_UCA',
                  'S_UCC', 'S_UCG', 'S_UCU', 'R_AGA', 'R_AGG', 'R_CGA', 'R_CGC', 'R_CGG', 'R_CGU', 'T_ACA',
                  'T_ACC', 'T_ACG', 'T_ACU', 'V_GUA', 'V_GUC', 'V_GUG', 'V_GUU', 'Y_UAC', 'Y_UAU']
        prefix = 'FOP'
        if tissues == 'all':
            tissues = self.fop_dic.index.tolist()

        fop_dict = self.fop_dic.to_dict(orient='index')  # Transform dataframe to dictionary for quicker access
        fop_all_values = {}
        for row in self.masterdf.itertuples():
            total_index = self.masterdf.columns.get_loc('Total_Num_Codons') + 1  # Get index of column with codon count
            fop_row_values = {}
            for tissue in tissues:
                fop = 0
                for codon in codons:
                    if fop_dict[tissue][codon] == 1:
                        col_index = self.masterdf.columns.get_loc(codon + '_Num') + 1
                        fop += row[col_index]
                if row[total_index] == 0:
                    fop = np.nan
                else:
                    fop = fop / row[total_index]
                fop_row_values[prefix + "_" + tissue] = fop
            fop_all_values[row[0]] = fop_row_values

        print('making dataframe fop')
        fop_df = pd.DataFrame.from_dict(fop_all_values, orient='index')
        print('joining')
        self.masterdf = self.masterdf.join(fop_df)

    def export_csv(self, out_name, out_path=os.path.dirname(os.path.realpath(__file__)) + '//Output//'):
        self.masterdf.to_csv(out_path + out_name, index=False)

    def import_csv(self, csv_path):
        self.masterdf = pd.read_csv(csv_path)


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
  if(missing(gene_ann_data) && missing(strain_ann_data)) { #Make non-annotated heatmap if there is no annotation data
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

    def __init__(self, data, gene_annotation_path=None, strain_annotation_path=None,
                 regex="Possible|note|RL5A_+|RL6_+"):
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

            # If annotation file has been loaded, add it when calling R function
            if self.gene_ann is None:
                self.RFunction.Gene_RSCU_Clustermap(graph_values, outvar)
            elif self.gene_ann is not None:
                self.RFunction.Gene_RSCU_Clustermap(graph_values, outvar, self.gene_ann)

        elif graph_type == 'Strain_RSCU_Clustermap':

            # group by Strain, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Strain_ID')[self.rscu_columns].mean()

            # Check number of strains in data, because clustering algorithm requires 2 rows/columns
            if len(graph_values.index) > 1:
                # If annotation file has been loaded, add it when calling R function
                if self.strain_ann is None:
                    self.RFunction.Strain_RSCU_Clustermap(graph_values, outvar)
                elif self.strain_ann is not None:
                    self.RFunction.Strain_RSCU_Clustermap(graph_values, outvar, self.strain_ann)

        elif graph_type == 'ENC_Heatmap':

            graph_values = self.masterdf[['Strain_ID', 'Gene', 'ENC']]  # Obtain relevant columns

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

            # If annotation file has been loaded, add it when calling R function
            if self.gene_ann is None:
                self.RFunction.Gene_RSCU_Clustermap_Di(graph_values, outvar)
            elif self.gene_ann is not None:
                self.RFunction.Gene_RSCU_Clustermap_Di(graph_values, outvar, self.gene_ann)

        elif graph_type == 'Gene_RSCU_Clustermap_Tri':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Gene')[self.tri_aa].mean()

            # If annotation file has been loaded, add it when calling R function
            if self.gene_ann is None:
                self.RFunction.Gene_RSCU_Clustermap_Tri(graph_values, outvar)
            elif self.gene_ann is not None:
                self.RFunction.Gene_RSCU_Clustermap_Tri(graph_values, outvar, self.gene_ann)

        elif graph_type == 'Gene_RSCU_Clustermap_Tetra':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Gene')[self.tetra_aa].mean()

            # If annotation file has been loaded, add it when calling R function
            if self.gene_ann is None:
                self.RFunction.Gene_RSCU_Clustermap_Tetra(graph_values, outvar)
            elif self.gene_ann is not None:
                self.RFunction.Gene_RSCU_Clustermap_Tetra(graph_values, outvar, self.gene_ann)

        elif graph_type == 'Gene_RSCU_Clustermap_Hexa':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Gene')[self.hexa_aa].mean()

            # If annotation file has been loaded, add it when calling R function
            if self.gene_ann is None:
                self.RFunction.Gene_RSCU_Clustermap_Hexa(graph_values, outvar)
            elif self.gene_ann is not None:
                self.RFunction.Gene_RSCU_Clustermap_Hexa(graph_values, outvar, self.gene_ann)

        elif graph_type == 'Strain_RSCU_Clustermap_Di':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Strain_ID')[self.di_aa].mean()

            if len(graph_values.index) > 1:
                # If annotation file has been loaded, add it when calling R function
                if self.strain_ann is None:
                    self.RFunction.Strain_RSCU_Clustermap_Di(graph_values, outvar)
                elif self.strain_ann is not None:
                    self.RFunction.Strain_RSCU_Clustermap_Di(graph_values, outvar, self.strain_ann)

        elif graph_type == 'Strain_RSCU_Clustermap_Tri':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Strain_ID')[self.tri_aa].mean()

            if len(graph_values.index) > 1:
                # If annotation file has been loaded, add it when calling R function
                if self.strain_ann is None:
                    self.RFunction.Strain_RSCU_Clustermap_Tri(graph_values, outvar)
                elif self.strain_ann is not None:
                    self.RFunction.Strain_RSCU_Clustermap_Tri(graph_values, outvar, self.strain_ann)

        elif graph_type == 'Strain_RSCU_Clustermap_Tetra':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Strain_ID')[self.tetra_aa].mean()

            if len(graph_values.index) > 1:
                # If annotation file has been loaded, add it when calling R function
                if self.strain_ann is None:
                    self.RFunction.Strain_RSCU_Clustermap_Tetra(graph_values, outvar)
                elif self.strain_ann is not None:
                    self.RFunction.Strain_RSCU_Clustermap_Tetra(graph_values, outvar, self.strain_ann)

        elif graph_type == 'Strain_RSCU_Clustermap_Hexa':
            # group by Gene, get only RSCU columns and calculate each one's mean
            graph_values = self.masterdf.groupby('Strain_ID')[self.hexa_aa].mean()

            if len(graph_values.index) != 1:
                # If annotation file has been loaded, add it when calling R function
                if self.strain_ann is None:
                    self.RFunction.Strain_RSCU_Clustermap_Hexa(graph_values, outvar)
                elif self.strain_ann is not None:
                    self.RFunction.Strain_RSCU_Clustermap_Hexa(graph_values, outvar, self.strain_ann)

        elif graph_type == 'FOP_Gene_Clustermap':
            fop_columns = self.masterdf.filter(regex="^FOP_*").columns.tolist()
            graph_values = self.masterdf.groupby('Gene')[fop_columns].mean()

            # If annotation file has been loaded, add it when calling R function
            if self.gene_ann is None:
                self.RFunction.FOP_Gene_Clustermap(graph_values, outvar)
            elif self.gene_ann is not None:
                self.RFunction.FOP_Gene_Clustermap(graph_values, outvar, self.gene_ann)

        elif graph_type == 'FOP_Strain_Clustermap':
            fop_columns = self.masterdf.filter(regex="^FOP_*").columns.tolist()
            graph_values = self.masterdf.groupby('Strain_ID')[fop_columns].mean()

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

            # If annotation file has been loaded, add it when calling R function
            if self.gene_ann is None:
                self.RFunction.FOP_SamplesGenes_Clustermap(graph_values, outvar)
            elif self.gene_ann is not None:
                self.RFunction.FOP_SamplesGenes_Clustermap(graph_values, outvar, self.gene_ann)

        elif graph_type == 'FOP_SamplesStrains_Clustermap':
            fop_columns = self.masterdf.filter(regex="^FOP_*").columns.tolist()
            graph_values = self.masterdf.groupby('Strain_ID')[fop_columns].mean()

            if len(graph_values.index) != 1:
                # If annotation file has been loaded, add it when calling R function
                if self.strain_ann is None:
                    self.RFunction.FOP_SamplesStrains_Clustermap(graph_values, outvar)
                elif self.strain_ann is not None:
                    self.RFunction.FOP_SamplesStrains_Clustermap(graph_values, outvar, self.strain_ann)

        elif graph_type == 'ENC_GC3_Genes':
            graph_values = self.masterdf.groupby('Gene')[['ENC', 'GC3']].mean()

            self.RFunction.ENC_GC3_Genes(graph_values, outvar)

        elif graph_type == 'ENC_GC3_Strains':
            graph_values = self.masterdf.groupby('Strain_ID')[['ENC', 'GC3']].mean()

            if len(graph_values.index) != 1:
                self.RFunction.ENC_GC3_Strains(graph_values, outvar)

        elif graph_type == 'ENC_GC3_ALL':
            graph_values = self.masterdf[['Gene', 'Strain_ID', 'ENC', 'GC3']]

            self.RFunction.ENC_GC3_ALL(graph_values, outvar)


Input_Layout = [
    [sg.Text('Import coding sequences:', font=('Any', 15, 'underline bold'))],
    [sg.Radio('Directory', "FileRadio", enable_events=True, key='_DirRadio_'),
     sg.In(key='_DirPath_',size=(50,None)),
     sg.FolderBrowse(button_text='Choose File Directory', target='_DirPath_',
                     tooltip='Choose a directory with fasta files that are to be analysed.')],
    [sg.Radio('Load Fasta Files', "FileRadio", enable_events=True, key='_FastaRadio_'),
     sg.In(key='_FastaPath_', size=(50,None)),
     sg.FilesBrowse(button_text='Choose fasta files', file_types=(('FASTA', '*.fasta'), ('FASTA', '*.fa'),
                                                                  ('Text', '*.txt')),
                    target='_FastaPath_', tooltip='Choose 1 or more fasta files that are to be analysed')],
    [sg.Text('If importing fasta files, make sure the name of the file is the genome/strain identifier.')],
    # [sg.Text('Alternatively, import csv file previously generated by CUVA:', font=('Any', 15, 'underline bold'))],
    # [sg.Radio('Load pre-made CSV File', "FileRadio", enable_events=True, key='_CSVRadio_'),
    #  sg.In(key='_CSVDir_'),
    #  sg.FileBrowse(button_text='Load CSV file', file_types=(('CSV', '*.csv'),), target='_CSVDir_',
    #                tooltip='Choose 1 csv file already generated by this program.')],
    # [sg.Text('If importing a csv file, make sure to check/uncheck FOP calculation depending on your needs.')],
    # [sg.Checkbox('Imported file already contains FOP values', key='_masterfileFOP_')]
]
RSCU_Layout = [
    [sg.Checkbox('Calculate RSCU', default=True, disabled=True, key='_RSCUCheck_')],
    [sg.Checkbox('Enable RSCU Quality Control', default=True, key='_RSCUQC_')],
    [sg.Frame(layout=[
        [sg.Checkbox('Discard RSCU if degeneracy of amino acid is larger than amino acid count in a gene.',
                     default=True, disabled=True, key='_RSCUQCGene_')],
        [sg.Checkbox('Discard RSCU if the count of non-empty (NaN) values of a codon is less than...',
                     default=True, disabled=True),
         sg.In('15', key='_RSCUQCThr_')]],
        title='Analysis Settings', relief='raised')]
]
ENC_Layout = [
    [sg.Checkbox('Calculate ENC', default=True, disabled=True, key='_ENCCheck_')],
    [sg.Checkbox('Calculate GC3', default=True, disabled=True, key='_GC3Check_')]
]
FOP_Layout = [
    [sg.Checkbox('Calculate FOP', default=False, key='_FOPCheck_'),
     sg.Text('WARNING: Enabling FOP will significantly increase analysis time', background_color='red')],
    [sg.Frame(layout=[
        [sg.Text('Optimal Codon Table (including type) Path:'),
         sg.In('None selected', key='_FOPRefPath_',size=(40,None)),
         sg.FileBrowse(button_text='Choose file', file_types=(('CSV', '*.csv'),), target='_FOPRefPath_')],
        [sg.Text('Sample ID Column Name: '), sg.In('Tissue', key='_FOPTissueColumn_')],
        [sg.Text('Tissue Type Column Name: '), sg.In('SMTS', key='_FOPSampleColumn_')]
    ], title='Reference File Settings', relief='raised')],
    [sg.Frame(layout=[
        # [sg.Radio('Group by type and use modal result as optimal codon value', 'FOPcalcRadio', key='_FOPModal_')],
        # [sg.Radio('Use optimal values of 1 sample from each type', 'FOPcalcRadio', key='_FOPRandom_')],
        # [sg.Radio('Calculate FOP for each sample', 'FOPcalcRadio',
        #          key='_FOPAll_')],
        [sg.Text('Which tissues to calculate FOP for (default all). Separate tissues by comma. '),
         sg.In('all', key='_FOPTissues_')]
    ], title='FOP Calculation Settings', relief='raised')]
]
Output_Layout = [
    [sg.Text('Output Directory'),
     sg.InputText(default_text=os.path.dirname(os.path.realpath(__file__)) + '//Output//', key='_OutPath_',
                  size=(50,None)),
     sg.FolderBrowse(button_text='Save Directory', target='_OutPath_',
                     tooltip='Choose where the output will be saved. Defaults to script folder.')],
    [sg.Text('Export masterfile (will be saved at output directory) as...'),
     sg.InputText(default_text='masterfile.csv', key='_masterdfName_')],
    # [sg.Checkbox('If calculating FOP, export extra masterfile without FOP values.', default=True,
    #              key='_Extramasterdf_'),
    # sg.InputText(default_text='nofop_masterfile.csv', key='_ExtramasterdfName_')],
    [sg.Checkbox('Generate & save all possible figures.', default=True,
                 key='_GenFigures_')]
]
GenFigures_Layout = [
    [sg.Text('Figure regex filter'), sg.InputText('Possible|note|RL5A_+|RL6_+', key='_FigRegex_')],
    [sg.Checkbox('Annotate Genes in heatmaps based on file:', key='_GeneAnnot_'),
     sg.In('None selected', key='_GeneAnnotFile_', size=(50,None)),
     sg.FileBrowse(button_text='Choose file', file_types=(('CSV', '*.csv'),), target='_GeneAnnotFile_')],
    [sg.Checkbox('Annotate Genomes in heatmaps based on file:', key='_GenomeAnnot_'),
     sg.In('None selected', key='_GenomeAnnotFile_', size=(50,None)),
     sg.FileBrowse(button_text='Choose file', file_types=(('CSV', '*.csv'),), target='_GenomeAnnotFile_')],
]
CodonW_Layout = [
    [sg.Text('CodonW Location'),
     sg.InputText(os.path.dirname(os.path.realpath(__file__))+"\\ref_files\\codonW", key='_CodonWFile_',
                  size=(50,None)),
     sg.FileBrowse(file_types=(('exe', '*.exe'),), target='_CodonWFile_')],
    [sg.Text('CodonW Output Directory'),
     sg.InputText(os.path.dirname(os.path.realpath(__file__))+"\\codonw_out\\", key='_CodonWOutDir_', size=(50,None)),
     sg.FolderBrowse(target='_CodonWOutDir_')]
]

Window_Layout = [
    [sg.TabGroup(layout=[[
        sg.Tab('Input Settings', layout=Input_Layout),
        sg.Tab('RSCU Settings', layout=RSCU_Layout),
        sg.Tab('ENC Settings', layout=ENC_Layout),
        sg.Tab('FOP Settings', layout=FOP_Layout),
        sg.Tab('Figure Settings', layout=GenFigures_Layout),
        sg.Tab('Output Settings', layout=Output_Layout),
        sg.Tab('CodonW Settings', layout=CodonW_Layout)]])],
    [sg.Multiline('Awaiting user input...\n', size=(50, 5), disabled=True, autoscroll=True, key='_Console_')],
    [sg.Submit(button_text='RUN', key='_RUN_', size=(30, 2)), sg.Button('Reset to Default', key='_Reset_')]
]
window = sg.Window('CUVA (Codon Usage Visualisation & Analysis', Window_Layout, default_element_size=(12, 1))

figure_types = [
    'Gene_RSCU_Clustermap', 'Strain_RSCU_Clustermap',
    'ENC_Heatmap',
    'Gene_RSCU_Clustermap_Di', 'Gene_RSCU_Clustermap_Tri', 'Gene_RSCU_Clustermap_Tetra', 'Gene_RSCU_Clustermap_Hexa',
    'Strain_RSCU_Clustermap_Di', 'Strain_RSCU_Clustermap_Tri', 'Strain_RSCU_Clustermap_Tetra',
    'Strain_RSCU_Clustermap_Hexa',
    'FOP_Gene_Clustermap', 'FOP_Strain_Clustermap',
    'FOP_SamplesGenes_Clustermap', 'FOP_SamplesStrains_Clustermap',
    'ENC_GC3_Genes', 'ENC_GC3_Strains', 'ENC_GC3_ALL',
    'FOP_Ref_Clustermap'
]
# print('Beginning event handling...')
while True:
    event, values = window.Read()
    print(event, values)
    if (event is None) or (event == 'Exit'):
        break
    if event == '_Reset_':
        # Clear Input
        window.Element('_DirPath_').Update('')
        window.Element('_FastaPath_').Update('')
        # window.Element('_CSVDir_').Update('')

        # Reset RSCU
        window.Element('_RSCUQC_').Update(True)
        window.Element('_RSCUQCThr_').Update('15')

        # Reset FOP
        window.Element('_FOPCheck_').Update(False)
        window.Element('_FOPRefPath_').Update('None selected')
        window.Element('_FOPTissueColumn_').Update('SMTS')
        window.Element('_FOPSampleColumn_').Update('Tissue')
        window.Element('_FOPTissues_').Update('all')

        # Reset Output
        window.Element('_OutPath_').Update(os.path.dirname(os.path.realpath(__file__)) + '//Output//')
        window.Element('_GenFigures_').Update(True)
        window.Element('_masterdfName_').Update('masterfile.csv')
        # window.Element('_Extramasterdf_').Update(True)
        # window.Element('_ExtramasterdfName_').Update('nofop_masterfile.csv')

        # Reset Figure Settings
        window.Element('_FigRegex_').Update('Possible|note|RL5A_+|RL6_+')
        window.Element('_GeneAnnot_').Update(False)
        window.Element('_GenomeAnnot_').Update(False)

        # Reset CodonW
        window.Element('_CodonWFile_').Update(os.path.dirname(os.path.realpath(__file__))+"\\ref_files\\codonW")
        window.Element('_CodonWOutDir_').Update(os.path.dirname(os.path.realpath(__file__))+"\\codonw_out\\")

    if event == '_RUN_':
        # Obtaining Input files
        files = []
        if values['_DirRadio_'] is True:  # Directory Input
            for file in os.listdir(values['_DirPath_']):
                if file.endswith(".fasta"):
                    files.append(os.path.abspath(os.path.join(values['_DirPath_'], file)))
        elif values['_FastaRadio_'] is True:  # Fasta File(s) Input
            files = values['_FastaPath_'].split(';')
        elif values['_CSVRadio_'] is True:  # CSV Input
            files = values['_CSVDir_']
        # Create model class with CodonW dir and output
        main_model = CodonData(codonw_dir=values['_CodonWFile_'], codonw_out_dir=values['_CodonWOutDir_'])

        # Calculate RSCU, ENC, GC3 if Input is *not* csv
        if (values['_DirRadio_'] is True) or (values['_FastaRadio_'] is True):
            for file in files:
                out = main_model.iterate_populate_master_file(os.path.dirname(file) + "\\", os.path.basename(file))
        elif values['_CSVRadio_'] is True:
            main_model.import_csv(values['_CSVDir_'])

        # Perform quality controls for RSCU  (optional)
        if values['_RSCUCheck_'] is True:
            main_model.qc_rscu(int(values['_RSCUQCThr_']))

        # Calculate FOP  (optional)
        if values['_FOPCheck_'] is True:
            # Export non-fop masterfile (optional)
            if values['_Extramasterdf_'] is True:
                main_model.export_csv(values['_ExtramasterdfName_'], values['_OutPath_'] + "\\")

            main_model.set_fop_ref(values['_FOPRefPath_'])

            fop_type = None
            if values['_FOPModal_'] is True:
                fop_type = 'mode'
            elif values['_FOPRandom_'] is True:
                fop_type = 'random'
            elif values['_FOPAll_'] is True:
                fop_type = 'all'
            main_model.set_fop_dict(fop_type, type_column_name=values['_FOPTissueColumn_'],
                                    tissue_column_name=values['_FOPSampleColumn_'])

            if values['_FOPTissues_'] != 'all':
                tissues = values['_FOPTissues_'].split(",")
            else:
                tissues = 'all'
            main_model.calculate_fop(tissues=tissues)

        # Export file
        main_model.export_csv(values['_masterdfName_'], values['_OutPath_'] + "\\")

        # Export Figures
        if values['_GenFigures_'] is True:
            if values['_GeneAnnot_'] is True:
                gene_annot = values['_GeneAnnotFile_']
            else:
                gene_annot = None
            if values['_GenomeAnnot_'] is True:
                genome_annot = values['_GenomeAnnotFile_']
            else:
                genome_annot = None
            figures = FigureGen(data=main_model.masterdf, regex=values['_FigRegex_'],
                                gene_annotation_path=gene_annot, strain_annotation_path=genome_annot)
            if values['_FOPCheck_'] is True:
                pass
            elif values['_masterfileFOP_'] is True:
                figure_types = [figure for figure in figure_types if figure not in {'FOP_Ref_Clustermap'}]
            else:
                figure_types = [figure for figure in figure_types if figure not in {'FOP_Gene_Clustermap',
                                                                                    'FOP_Strain_Clustermap',
                                                                                    'FOP_SamplesGenes_Clustermap',
                                                                                    'FOP_SamplesStrains_Clustermap',
                                                                                    'FOP_Ref_Clustermap'}]
            for figure in figure_types:
                if figure == 'FOP_Ref_Clustermap':
                    fop_ref = main_model.fop_ref_df
                    fop_ref.set_index(values['_FOPSampleColumn_'])

                    fop_ref_figgen = FigureGen(fop_ref)
                    fop_ref_figgen.make_graph(figure, out_path=values['_OutPath_'] + "\\")
                else:
                    try:
                        figures.make_graph(figure, out_path=values['_OutPath_'] + "\\")
                    except:
                        print("Could not make figure: " + figure)

window.Close()

