# PaGeSearch (Pathway Genome Search)
### A Tool for Identifying Genes within Pathways in Unannotated Genomes

PaGeSearch can identify a list of genes or genes in a pathway efficiently from the genome without any annotation information.
PaGeSearch combines sequence similarity search, gene prediction, and neural network models to identify the most possible orthologs of the genes in the list of pathway.

## Installation
You can download the source code of PaGeSearch and the filtering model according to your species from the github repositories.
Do not change the name of the model file for further use.
### Dependencies
#### Python libraries
* python	3.7
* pandas	1.3.5
* Bio	1.79
* keras	2.11.0
#### Softwares
* seqkit	2.4
* mmseqs2	14.7
* samtools	1.17
* bedtools	2.31
* augustus	3.5
* exonerate	2.4

## Running PaGeSearch
PaGeSearch is designed to be applicable for five taxonomic classes including an archetype species.

|   Species   |    Scientific Name    | Applicable Class |
|-------------|-----------------------|------------------|
|    human    | Homo sapiens          | Mammals          |
|  zebrafish  | Danio rerio           | Fish             |
|   chicken   | Gallus gallus         | Birds            |
| arabidopsis | Arabidopsis thaliana  | Eudicotyledons   |
|    wheat    | Triticum aestivum     | Liliopsida       |

You need to provide the paths of the genome and the query gene sequences, the archetype species specified, and the path of the model downloaded to run PaGeSearch. 
You can specify the output directory and prefix, and the number of threads.

### Run PaGeSearch
Assume that you have downloaded the query genes in a folder namned "pathway", and have the fasta file of the genome at genome.fa 
Also, assume that you have downloaded the Homo sapiens model at a folder "names path_to_model/Homo_sapiens_NN_model.rds". 
```
python pagesearch.py -g genome.fa -p pathway -od pagesearch_results -op test -s human -m path_to_model -t 4
```
Your results will be saved at the pagesearch_results folder, as files named test.txt and testb.bed.

### Download gene sequences 
You can download the seqeucnes of genes in a pathway or a list, that can directly be used as the query of PaGeSearch.
The pathway information is based on the Reactome database, and the pathway name must be entered exactly as it appears in the Reactome database.
Refer to https://reactome.org/ for animal species and https://plantreactome.gramene.org/index.php?lang=en for plant species.

#### Download gene sequences from pathway name
Download sequences of genes in the pathway 'Metabolism of nucleotides' of Homo sapiens and save to a folder named metabolism_of_nucleotides_gene_sequences.
The sequences of each gene is saved to separate fasta files.
```
python download_gene_sequences.py -p "Metabolism of nucleotides" -s human -t 9606 -o metabolism_of_nucleotides_gene_sequences
```
You can run PaGeSearch using the downloaded gene sequences as:
```
python pagesearch.py -g genome.fa -p metabolism_of_nucleotides_gene_sequences -od pagesearch_results -op metabolism_of_nucleotides -s human -t 4
```
Download all sequences that are orthologs of the genes in the humen 'Metabolism of nucleotides' pathway within mammali species.
The taxonomy ID is based on NCBI taxonomy (https://www.ncbi.nlm.nih.gov/taxonomy).
```
python download_gene_sequences.py -p "Metabolism of nucleotides" -s human -t 40674 -o metabolism_of_nucleotides_mammals_orthologs
```
#### Download gene sequences from a user defined list
The gene list is a text file of Ensembl gene IDs separated by line breaks.
```
python download_gene_sequences.py -l genelist.txt -s human -t 9606 -o gene_sequences
```

### Options
#### pagesearch.py
|           Option           |                    Description                    |  Default   |
|----------------------------|---------------------------------------------------|------------|
| -g, --genome               | File path of the genome sequence                  | None       |
| -p. --pathway-geneseq-dir  | Folder path where the query gene sequences are in | None       |
| -od, --outdir              | Folder path where outputs will be saved           | ./         |
| -op, --outprefix           | Prefix of output files                            | pagesearch |
| -s, --species              | Archetype species                                 | human      |
| -t, --threads              | Number of threads used                            | 4          |

#### download_gene_sequences.py
|     Option     |               Description                |      Default       |
|----------------|------------------------------------------|--------------------|
| -p, --pathway  | Reactome pathway name                    | None               |
| -l, --genelist | Path to user-defined gene list txt file  | None               |
| -s, --species  | Archetype species                        | human              |
| -t, --taxon    | NCBI taxon ID to download orthologs from | 9606               |
| -o, --outdir   | Folder path where outputs will be saved  | ./pagesearch_query |


## Contact
wsy415@snu.ac.kr
