# PaGeSearch (Pathway Genome Search)
### A Tool for Identifying Genes within Pathways in Unannotated Genomes

PaGeSearch can identify a list of genes or genes in a pathway efficiently from the genome without any annotation information.
PaGeSearch combines sequence similarity search, gene prediction, and neural network models to identify the most possible orthologs of the genes in the list of pathway.

## Installation
You can download the source code of PaGeSearch and the filtering model according to your species from the github repositories.
Do not change the name of the model file for further use.
### Codes:
### Models:

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
### Usage
#### Run PaGeSearch
Assume that you have downloaded the query genes in a folder namned "pathway", and have the fasta file of the genome at genome.fa 
Also, assume that you have downloaded the Homo sapiens model at a folder "names path_to_model/Homo_sapiens_NN_model.rds". 
```
python pathway_gene_search.py -g genome.fa -p pathway -od pagesearch_results -op test -s human -m path_to_model -t 4
```
Your results will be saved at the pagesearch_results folder, as files named test.txt and testb.bed.

You can download the seqeucnes of genes in a pathway or a list, that can directly be used as the query of PaGeSearch.
The pathway information is based on the Reactome database, and the pathway name must be entered exactly as it appears in the Reactome database.
Refer to https://reactome.org/ for animal species and https://plantreactome.gramene.org/index.php?lang=en for plant species.
#### Download gene sequences from pathway name
Download sequences of glycolysis related genes of Homo sapiens and save to a folder named Glycolysis_gene_sequences.
The sequences of each gene is saved to separate fasta files.
```
python download_gene_sequences.py -p Glycolysis -s human -t 9606 -o Glycolysis_gene_sequences
```
You can run PaGeSearch using the downloaded gene sequences as:
```
python pathway_gene_search.py -g genome.fa -p Glycolysis_gene_sequences -od pagesearch_results -op glycolysis -s human -m path_to_model -t 4
```
Download all sequences that are orthologs of human glycolysis genes in mammalia.
The taxonomy ID is based on NCBI taxonomy (https://www.ncbi.nlm.nih.gov/taxonomy).
```
python download_gene_sequences.py -p Glycolysis -s human -t 40674 -o Glycolysis_mammals_orthologs
```
#### Download gene sequences from a user defined list
The gene list is a text file of Ensembl gene IDs separated by line breaks.
```
python download_gene_sequences.py -l genelist.txt -s human -t 9606 -o gene_sequences
```

### Options

