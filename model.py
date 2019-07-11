"""
Created on 20 May 2019

@author: Antonia Chalka

"""
import pandas as pd
import os
import numpy as np
from subprocess import call

'''
Class representing codon data across strains and genes. Relies on codonW to perform RSCU and ENC analysis.
All data is stored and manipulated in masterdf, a pandas data frame.
Its methods include:
 make_rscu_file + make_enc_file: methods that call CodonW and retrieve codon information data.
 make_cai_file: tbd
 makeFOPvalues: tbd
 parse_masterfile_values: gather codon data from the file-making methods, and use them to create rows to poulate master
 qc_rscu: quality control RSCU. eplaces some RSCU cells with nan depending on quality criteria
 
 export_csv: export masterdf to CSV file. May be moved
 make_graph: make Graph depending on keyword, May be moved
 iteratePopulatemastefile: take a directory and iterate through fasta files to add them to masterdf using above methods. 
 May be moved.

'''


class CodonData(object):
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
                 codonw_dir=os.path.dirname(os.path.realpath(__file__))+"\\ref_files\\codonW",
                 codonw_out_dir=os.path.dirname(os.path.realpath(__file__))+"\\codonw_out\\"):
        print("Initialising model data frame...")
        self.masterdf = pd.DataFrame(columns=self.col_names)  # Set up masterdf
        self.fop_type = fop_type
        self.fop_ref_df = fop_ref_df
        self.fop_dic = fop_dic
        self.codonw_dir = codonw_dir
        self.codonw_out_dir = codonw_out_dir + "\\"

    def make_rscu_file(self, input_file, input_dir):
        print("Calculating RSCU values for file:" + input_file)

        # command summary = codonw inputfile -silent -nomenu -rscuoutname1 outname2
        cmd = self.codonw_dir + " " + input_dir + input_file + ".fasta" +\
            ' -silent -nomenu -rscu ' +\
            self.codonw_out_dir + input_file + '.none ' +\
            self.codonw_out_dir + input_file + '.blk'
        print("Running CodonW command: " + cmd)
        call(cmd, shell=True)
        
        # print("Parsing RSCU file...")
        rscu_file = open(self.codonw_out_dir + input_file + '.blk', "r")
        sections = rscu_file.read().split("\n\n")  # split when 2 newlines (1 gene has 5 of such sections)
        sections = [x.split() for x in sections]  # remove whitespace
        # remove Amino Acid and Codon Names
        sections = [[i for i in x if i not in self.AANames and i not in self.CodonNames] for x in sections]
        rscu_file.close()

        return sections
    
    def make_enc_file(self, input_file, input_dir):
        print("Calculating ENC values for file:" + input_file)
        # command summary = codonw inputfile (options) rscuoutname1 outname2
        cmd = self.codonw_dir + " " + input_dir + input_file + ".fasta" +\
            ' -silent -nomenu -enc -machine ' +\
            self.codonw_out_dir + input_file + '.enc ' +\
            self.codonw_out_dir + input_file + '.blk'
        print("Running CodonW command: " + cmd)
        call(cmd, shell=True)

        print("Parsing ENC file...")
        enc_file = open(self.codonw_out_dir + input_file + ".enc", "r")
        enc_rows = enc_file.read().split("\n")  # Extract rows
        enc_file.close()
        
        return enc_rows

    def make_gc3_file(self, input_file, input_dir):
        print("Calculating GC3 values for file: " + input_file)
        # command summary = codonw inputfile (options) rscuoutname1 outname2
        cmd = self.codonw_dir + " " + input_dir + input_file + ".fasta" + \
            ' -silent -nomenu -gc3s -machine ' + \
            self.codonw_out_dir + input_file + '.gc3 ' + \
            self.codonw_out_dir + input_file + '.blk'
        print("Running CodonW command: " + cmd)
        call(cmd, shell=True)

        print("Parsing GC3 file...")
        gc3_file = open(self.codonw_out_dir + input_file + ".gc3", "r")
        gc3_rows = gc3_file.read().split("\n")  # Extract rows
        gc3_file.close()

        return gc3_rows

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

        rscu_sections = self.make_rscu_file(input_file, file_path + "\\")
        enc_rows = self.make_enc_file(input_file, file_path + "\\")
        gc3_rows = self.make_gc3_file(input_file, file_path + "\\")

        self.parse_masterfile_values(input_file, rscu_sections, enc_rows, gc3_rows)

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
                print(fop)
                fop_row_values[prefix + "_" + tissue] = fop
            fop_all_values[row[0]] = fop_row_values

        print('making dataframe fop')
        fop_df = pd.DataFrame.from_dict(fop_all_values, orient='index')
        print('joining')
        self.masterdf = self.masterdf.join(fop_df)

        
    def export_csv(self, out_name, out_path=os.path.dirname(os.path.realpath(__file__)) + '//Output//'):
        self.masterdf.to_csv(out_path+out_name, index=False)
        
    def import_csv(self, csv_path):
        self.masterdf = pd.read_csv(csv_path)


model = CodonData()

# filepath= os.path.dirname(os.path.realpath(__file__))
filepath = os.path.abspath("C:\\Users\\anni1\\PycharmProjects\\MScProject\\CodonUsage\\fasta_cds")
for file in os.listdir(filepath):
    if file.endswith(".fasta"):
        model.iterate_populate_master_file(filepath, file)

'''
model.calculate_fop(tissues)
print(model.masterdf)
model.export_csv('fopmode_masterfile.csv')

model.import_csv('no_masterfile.csv')
model.set_fop_dict('random')
model.calculate_fop(tissues)
print(model.masterdf)
model.export_csv('foprandom_masterfile.csv')


# model.import_csv('foptissue_masterfile.csv')

# model.set_fop_dict('mode')
# model.calculate_fop('all')
# model.export_csv('fopall_masterfile.csv')

'''
