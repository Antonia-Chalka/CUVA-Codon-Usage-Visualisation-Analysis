"""
Created on 20 May 2019

@author: Antonia Chalka

"""
import pandas as pd
import os
import numpy as np

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
                 'ENC', 'CAI', 'GC3']
    
    # Used for qc_rscu
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
    # Threshold for how many times a codon column must have a non-NaN value.
    # Above threshold counts are retained and vice versa.
    RSCUQC_Thres = 15

    def __init__(self):
        print("Initialising model data frame...")
        self.masterdf = pd.DataFrame(columns=self.col_names)  # Set up masterdf

    def make_rscu_file(self, input_file):
        print("Calculating RSCU values for file:" + input_file)
        
        cmd = 'CodonW ' + input_file + '.fasta -nomenu -silent -rscu'
        print("Running CodonW command: " + cmd)
        os.system(cmd)
        
        print("Parsing RSCU file...")
        rscu_file = open(input_file + ".blk", "r")
        sections = rscu_file.read().split("\n\n")  # split when 2 newlines (1 gene has 5 of such sections)
        sections = [x.split() for x in sections]  # remove whitespace
        # remove Amino Acid and Codon Names
        sections = [[i for i in x if i not in self.AANames and i not in self.CodonNames] for x in sections]
        rscu_file.close()
        
        return sections
    
    def make_enc_file(self, input_file):
        print("Calculating ENC values for file:" + input_file)
        cmd = "codonW " + input_file + ".fasta -nomenu -silent -enc -machine "
        print("Running CodonW command: " + cmd)
        os.system(cmd)
        print("Parsing ENC file...")
        enc_file = open(input_file + ".out", "r")
        enc_rows = enc_file.read().split("\n")  # Extract rows
        enc_file.close()
        
        return enc_rows

    def make_gc3_file(self, input_file):
        print("Calculating GC3 values for file: " + input_file)
        cmd = "codonW " + input_file + ".fasta -nomenu -gc3s -machine -silent"
        print("Running CodonW command: " + cmd)
        os.system(cmd)
        print("Parsing GC3 file...")
        gc3_file = open(input_file + ".out", "r")
        gc3_rows = gc3_file.read().split("\n")  # Extract rows
        gc3_file.close()

        return gc3_rows
    # TODO ADD CALCULATE GC3 files

    def make_cai_file(self, input_file):
        # TBD LATER
        cai_rows = []
        return cai_rows

    def make_fop_ref(self, codon_file, info_file):
        fop_codon_df = pd.read_csv(codon_file, sep='\t', header=0)
        sample_ref_df = pd.read_csv(info_file, sep='\t', header=0, usecols=['SAMPID', 'SMTS', 'SMTSD'])
        self.fop_ref_df = fop_codon_df.set_index('Tissue').join(sample_ref_df.set_index('SAMPID'))
        self.fop_ref_df = self.fop_ref_df.dropna(axis=1)
        self.fop_ref_df.to_csv('fop_ref.csv')

    def set_fop_ref(self, fopref_path):
        print("reading fop file...")
        self.fop_ref_df = pd.read_csv(fopref_path)

    def set_fop_dict(self, foptype):
        self.foptype = foptype  # TODO MOVE
        print("setting reference..." + foptype)
        columns = ['A_GCA', 'A_GCC', 'A_GCG', 'A_GCU', 'C_UGC', 'C_UGU', 'E_GAA', 'E_GAG', 'D_GAC', 'D_GAU',
                   'G_GGA', 'G_GGC', 'G_GGG', 'G_GGU', 'F_UUC', 'F_UUU', 'I_AUA', 'I_AUC', 'I_AUU', 'H_CAC',
                   'H_CAU', 'K_AAA', 'K_AAG', 'L_CUA', 'L_CUC', 'L_CUG', 'L_CUU', 'L_UUA', 'L_UUG', 'N_AAC',
                   'N_AAU', 'Q_CAA', 'Q_CAG', 'P_CCA', 'P_CCC', 'P_CCG', 'P_CCU', 'S_AGC', 'S_AGU', 'S_UCA',
                   'S_UCC', 'S_UCG', 'S_UCU', 'R_AGA', 'R_AGG', 'R_CGA', 'R_CGC', 'R_CGG', 'R_CGU', 'T_ACA',
                   'T_ACC', 'T_ACG', 'T_ACU', 'V_GUA', 'V_GUC', 'V_GUG', 'V_GUU', 'Y_UAC', 'Y_UAU']
        if foptype == 'mode':
            self.fopdic = self.fop_ref_df.groupby('SMTS')[columns].agg(pd.Series.mode)
        elif foptype == 'random':
            self.fopdic = self.fop_ref_df.groupby('SMTS')[columns].first()
        elif foptype == 'detail_all':
            columns.append('Tissue')
            ebv = self.fop_ref_df.loc[self.fop_ref_df['SMTSD'] == 'Cells - EBV-transformed lymphocytes'][columns]
            ebv = ebv.assign(Tissue=ebv.Tissue + '_EBV')

            fibro = self.fop_ref_df.loc[self.fop_ref_df['SMTSD'] == 'Cells - Transformed fibroblasts'][columns]
            fibro = fibro.assign(Tissue=fibro.Tissue + '_FIBRO')

            blood = self.fop_ref_df.loc[self.fop_ref_df['SMTSD'] == 'Whole Blood'][columns].head(100)
            blood = blood.assign(Tissue=blood.Tissue + "_WB")

            skin = self.fop_ref_df[self.fop_ref_df['SMTSD'].str.contains("Skin")][columns].head(100)
            skin = skin.assign(Tissue=skin.Tissue + "_Skin")

            self.fopdic = pd.concat([ebv, fibro, blood, skin], ignore_index=True)
            self.fopdic.set_index('Tissue', drop=True, inplace=True)
        elif foptype == 'temp':

            lung = self.fop_ref_df.loc[self.fop_ref_df['SMTS'] == 'Lung'][columns].agg(pd.Series.mode)
            saliva = self.fop_ref_df.loc[self.fop_ref_df['SMTS'] == 'Salivary Gland'][columns].agg(pd.Series.mode)
            sintestine = self.fop_ref_df.loc[self.fop_ref_df['SMTS'] == 'Small Intestine'][columns].agg(pd.Series.mode)
            kidney = self.fop_ref_df.loc[self.fop_ref_df['SMTS'] == 'Kidney'][columns].agg(pd.Series.mode)
            liver = self.fop_ref_df.loc[self.fop_ref_df['SMTS'] == 'Liver'][columns].agg(pd.Series.mode)
            pancreas = self.fop_ref_df.loc[self.fop_ref_df['SMTS'] == 'Pancreas'][columns].agg(pd.Series.mode)


            # 3 tissues WB,
            pass


    def calculate_fop(self, tissues):
        codons = ['A_GCA', 'A_GCC', 'A_GCG', 'A_GCU', 'C_UGC', 'C_UGU', 'E_GAA', 'E_GAG', 'D_GAC', 'D_GAU',
                  'G_GGA', 'G_GGC', 'G_GGG', 'G_GGU', 'F_UUC', 'F_UUU', 'I_AUA', 'I_AUC', 'I_AUU', 'H_CAC',
                  'H_CAU', 'K_AAA', 'K_AAG', 'L_CUA', 'L_CUC', 'L_CUG', 'L_CUU', 'L_UUA', 'L_UUG', 'N_AAC',
                  'N_AAU', 'Q_CAA', 'Q_CAG', 'P_CCA', 'P_CCC', 'P_CCG', 'P_CCU', 'S_AGC', 'S_AGU', 'S_UCA',
                  'S_UCC', 'S_UCG', 'S_UCU', 'R_AGA', 'R_AGG', 'R_CGA', 'R_CGC', 'R_CGG', 'R_CGU', 'T_ACA',
                  'T_ACC', 'T_ACG', 'T_ACU', 'V_GUA', 'V_GUC', 'V_GUG', 'V_GUU', 'Y_UAC', 'Y_UAU']
        if tissues == 'all':
            tissues = self.fopdic.index.tolist()
        if self.foptype == 'detail_all':
            prefix = 'SFOP'
        elif self.foptype == 'mode' or self.foptype == 'random':
            prefix = 'FOP'

        fop_dict = self.fopdic.to_dict(orient='index')
        fop_all_values = {}
        i = 0
        for row in self.masterdf.itertuples():
            totaln_index = self.masterdf.columns.get_loc('Total_Num_Codons') + 1
            print(i)
            i += 1
            fop_row_values = {}
            for tissue in tissues:
                fop = 0
                for codon in codons:
                    if fop_dict[tissue][codon] == 1:
                        col_index = self.masterdf.columns.get_loc(codon + '_Num') + 1
                        fop += row[col_index]
                if row[totaln_index] == 0:
                    fop = np.nan
                else:
                    fop = fop / row[totaln_index]
                print(fop)
                fop_row_values[prefix + "_" + tissue] = fop
            fop_all_values[row[0]] = fop_row_values

        print('making dataframe fop')
        fop_df = pd.DataFrame.from_dict(fop_all_values, orient='index')
        print('joining')
        self.masterdf = self.masterdf.join(fop_df)

    def parse_masterfile_values(self, input_file, rscu_sections, enc_rows, cai_rows, gc3_rows):
        
        strain_name = os.path.splitext(input_file)[0]  # Get Strain Name
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
        print(enc_gene_list)
        print(enc_values)

        #GC3 section
        if len(enc_rows) != len(gc3_rows):
            print("Different lengths: " + strain_name + " ENC: " + str(len(enc_rows)) + " GC3: " + len(gc3_rows))

        print("Parsing GC3 values...")
        gc3_values = []
        gc3_gene_list = []
        for i in range(1, len(gc3_rows) - 1):
            gc3_fields = gc3_rows[i].split()
            gc3_gene_list.append(gc3_fields[0])

            if len(gc3_fields) == 1:
                print(gc3_fields)
                gc3_fields.append(np.nan)
                print(gc3_fields)
            gc3_values.append(float(gc3_fields[1]))
        print(str(gc3_values))

        # CAI section - TODO
        cai_values = [0]*len(enc_rows)
        
        print("Collating values in masterdf...")
        # Collate into masterfile rows
        for i in range(0, len(rscu_values)):
            masterfile_row = [strain_name, enc_gene_list[i]]
            masterfile_row.extend(rscu_values[i])
            masterfile_row.append(enc_values[i])
            masterfile_row.append(cai_values[i])
            masterfile_row.append(gc3_values[i])

            self.masterdf.loc[len(self.masterdf)] = masterfile_row

    def qc_rscu(self):  # Quality control RSCU  values
        # For aa
        print("Beginning quality control of RSCU values across amino acids...")
        
        for i in self.masterdf.index:
            print("Checking row " + str(i))
            
            for aa in self.SynonymousCodons:
                degeneracy = len(self.SynonymousCodons[aa]) 
                aa_count = 0
                
                for codon in self.SynonymousCodons[aa]:
                    numcolumn = codon+"_Num"
                    aa_count += self.masterdf.at[i, numcolumn]
                
                print("AA: " + aa + ". Codon (column): " + codon +
                      ". Degeneracy: " + str(degeneracy) + ". Count: " + str(aa_count))

                if degeneracy > aa_count:
                    print("Num of codons less than degeneracy. Filling with Nan values...")
                    # replace with missing values
                    for codon in self.SynonymousCodons[aa]:
                        rscu_column = codon + "_RSCU"
                        self.masterdf.at[i, rscu_column] = np.nan
        # Across strains
        print("Beginning quality control of RSCU values across strains...")
        
        grouped = self.masterdf.groupby("Strain_ID")  
        for strain, group in grouped:
            print("Checking codons for Strain:" + strain)
            
            for aa in self.SynonymousCodons:
                for codon in self.SynonymousCodons[aa]:
                    rscu_column = codon + "_RSCU"
                    if group[rscu_column].notna().count() < self.RSCUQC_Thres:
                        print("Codon with >15 reps ( " + strain + " " + codon + " ). Setting to Nan...")
                        self.masterdf.loc[self.masterdf.Strain_ID == strain, rscu_column] = np.nan
    
    def replace_gene_names(self, info_filepath):
        # Get values from external file
        print("Changing gene names...")
        for line in open(info_filepath, "r"):
            values = line.split(";")
            for value in values:
                if "gene=" in value:
                    gene = value[5:]
                elif "product=" in value:
                    extra_fields = value[8:]
            # Iterate through gene from external file (each line should have a gene)
            old_gene_name = gene + "_+"
            new_gene_name = gene + "_" + extra_fields
            print("Changing gene name from " + old_gene_name + " to " + new_gene_name)
            self.masterdf.loc[self.masterdf.Gene == old_gene_name, 'Gene'] = new_gene_name
    
    def iterate_populate_master_file(self, file_path):
        print("Changing directory to " + file_path)
        os.chdir(file_path)
        
        for file in os.listdir(file_path):
            if file.endswith(".fasta"):
                input_file = os.path.splitext(file)[0]
                rscu_sections = self.make_rscu_file(input_file)
                enc_rows = self.make_enc_file(input_file)
                cai_rows = self.make_cai_file(input_file)
                gc3_rows = self.make_gc3_file(input_file)

                self.parse_masterfile_values(input_file, rscu_sections, enc_rows, cai_rows, gc3_rows)
        self.qc_rscu()
        self.replace_gene_names("Extra_geneInfoField.txt")
        
    def export_csv(self, outfilepath):
        self.masterdf.to_csv(outfilepath, index=False)
        
    def import_csv(self, csv_path):
        self.masterdf = pd.read_csv(csv_path)
        

#  RUNNING IT

model = CodonData()

# filepath= os.path.dirname(os.path.realpath(__file__))
filepath = os.path.abspath("C:\\Users\\anni1\\PycharmProjects\\MScProject\\CodonUsage\\fasta_cds")
# model.iterate_populate_master_file(filepath)
# model.export_csv('masterfile.csv')

tissues = ['Lung', 'Salivary Gland', 'Small Intestine', 'Kidney', 'Liver', 'Pancreas']

model.import_csv('no_masterfile.csv')

model.set_fop_ref('fop_ref.csv')

model.set_fop_dict('mode')
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
Code used to make fop reference file:
fop_codon_df = pd.read_csv(codon_file, sep='\t', header=0)
sample_ref_df = pd.read_csv(info_file, sep='\t', header=0, usecols=['SAMPID', 'SMTS', 'SMTSD'])
self.fop_ref_df = fop_codon_df.set_index('Tissue').join(sample_ref_df.set_index('SAMPID'))
self.fop_ref_df = self.fop_ref_df.dropna(axis=1)
self.fop_ref_df.to_csv('fop_ref.csv')
'''
