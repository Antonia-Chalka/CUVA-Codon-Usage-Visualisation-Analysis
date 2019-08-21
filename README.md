# CUVA (Codon Usage Visualisation & Analysis)

## Requirements:
 
* CodonW installation 
* R installation in path

## Manual

The CUVA GUI is split into multiple tabs, separated for input/output, figure generation, codonW settings and codon usage bias metric settings (RSCU, ENC-GC3, FOP). The information here can also be found in the main report.

### Input 

CUVA takes one or more fasta file of coding sequences as input, with an example given below:

>Gene Name 1

ATGG......

>Gene Name 2

ATGTGACCC...

 The file name of each fasta file is used an identifier (can be strain/genome/cultivar name, or any other variable of choice).
 
### Codon Usage Bias Settings
The RSCU, ENC & FOP tabs contain settings for each calculation. RSCU, ENC  GC3 calculations are performed by default, FOP calculations require setup by the user.

RSCU settings allow the modification of quality controls (enabled by default).
The ENC settings cannot be modified in the current version.
In the FOP settings, the user can enable FOP calculations, but must provide an optimal codon reference file beforehand, with optimality being denoted as 1/0 (optimal/not optimal). The format of the optimal codon reference file must be as such (columns may be in any order):

| Sample ID | Tissue (Optional Column) | AUG | CGG | ... (Other Codon Columns) |
| --- | --- | --- | --- |--- |
| Sample 1 | Blood |0| 0| |
| Sample 2 | Kidney|0 | 1| |
| Sample 3 | Blood |1 | 1| |
| Sample 4 | Lung | 1| 0| |
| Sample 5 | Lung | 0| 1| |
...

The Tissue name column is optional, and can be used to select samples for which to calculate FOP for (FOP is calculated based on each optimal codons of each sample in the group).

### Figure Generation
The Figure settings allow the user to select whether CUVA will output.
CUVA will group metrics by gene or by genome/strain/cultivar. The user may filter out genes using the regex filter. The following graphs are generated:
* Heatmap of RSCU scores across all codons, averaged by gene or genome
* Heatmaps of RSCU scores of di/tri/tetra/hexa-degenerated codons, averaged by gene or genome
* ENC heatmap with genomes as rows and genes as columns
* ENC-plot, with a curve of expected ENC scores and with observed ENC scores averaged by genome or gene
* FOP heatmap, with FOP scores averaged by gene or genome
* Heatmap of optimal codons of the FOP reference file
The user can provide annotation files for the genes or genomes to create annotation colorbars for the heatmaps.
In the file, the first column must contain the gene or genomes, as shown below:

|Gene/Genome|	Annotation |
| --- | --- |
|Gene 1|	Group 1|
|Gene 2|	Group 2|
|Gene 3|	Group 3|
...

If the annotation fails, CUVA will output non-annotated figures.

### Output
The user can specify the directory where the file containing all the codon usage bias metrics as well as the figures (if the user selects to generate them).

### CodonW Settings
The user may specify the CodonW installation directory, as well as the directory of the outputted CodonW files.


