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
Assume that you have downloaded the query genes in a folder namned "pathway", and have the fasta file of the genome at genome.fa 
Also, assume that you have downloaded the Homo sapiens model at a folder "names path_to_model/Homo_sapiens_NN_model.rds". 
```
python pathway_gene_search.py -g genome.fa -p pathway -od pagesearch_results -op test -s human -m path_to_model -t 4
```
Your results will be saved at the pagesearch_results folder, as files named test.txt and testb.bed.
### Options


## pathway_gene_sequence_download.py


## pathway_gene_search.py
