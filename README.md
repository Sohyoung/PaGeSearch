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
* tensorflow 2.11.0
#### Softwares
* seqkit	2.4
* mmseqs2	14.7
* samtools	1.17
* bedtools	2.31
* augustus	3.5
* exonerate	2.4

## Running PaGeSearch
#### Choosing archetype species model
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

#### Preparing Your Input Files
First, you need two types of input files:
Genome Sequence File: This is a fasta file (.fa) containing the genome sequence you are interested in. For example, genome.fa.
Query Gene Sequences: These are the nucleotide sequences of the genes you want to find. Each gene's sequence should be in a separate fasta file.
Ensure all the fasta files for your query genes are placed in a directory, such as /directory/of/query/gene/sequences/folder.

#### Setting up PaGeSearch
Make sure you have PaGeSearch installed and know the path to its folder, for instance, /path/to/PaGeSearch.
Change your working directory to the PaGeSearch folder using the command line:
```
cd /path/to/PaGeSearch
```

#### Running PaGeSearch
Now run PaGeSearch with the following command:
```
python Codes/pagesearch.py -g /path/to/genome.fa -p /directory/of/query/gene/sequences/folder -od ../pagesearch_results -op test -s human -t 4
```
In this command:
-g specifies the path to your genome sequence file.
-p specifies the path to the folder containing your query gene sequences. Each gene should be in a separate fasta file within this folder.
-od determines where the results will be saved. In this case, results will be in the pagesearch_results folder located one level up from the current directory.
-op sets the prefix for the output files. Here, the output files will be named with the prefix 'test', resulting in test.bed and test.gff.
-s indicates the archetype species (in this case, 'human').
-t sets the number of threads to use (here, it's 4).

#### Understanding the Output
After running PaGeSearch, the results will be saved in the specified output directory. You'll find two files:
##### test.bed
This file contains information on the locations of the found genes in BED format.
Each row contains information of the chromosome number, start codon, stop codon, and gene name. An example is as follows.
```
1	3934136	3943753	FBA6
1	3934136	3938303	FBA8
1	3934136	3938303	FBA8
1	20740849	20743255	GAPC1
1	20740849	20743255	GAPC2
```
##### test.gff
This file provides detailed annotations in GFF format. For detailed information, please refer to the following link: https://en.wikipedia.org/wiki/General_feature_format. 
The first column indicates the chromosome number, the third column the feature type, the third and fourth columns represent the start and stop codons, and the last column contains information related to gene annotation. 
In the last column, the content following 'gene=' indicates the name of the gene being searched for. 
Although PaGeSearch aims to present the most probable gene model for each gene, in cases where prediction probabilities are similar (within 5% difference), it considers it a tie and may report multiple models. In such cases, the gene_id is indicated with _1, _2, etc. 
The probability represents the prediction probability calculated by PaGeSearch’s algorithm. An example is as follows.
```
1	PaGeSearch	gene	3934136	3938303	.	+	.	gene_id=FBA8_1;gene=FBA8;probability=0.4765756
1	PaGeSearch	transcript	3934136	3938303	.	+	.	gene_id=FBA8_1
1	PaGeSearch	start_codon	3934136	3934138	.	+	.	gene_id=FBA8_1
1	PaGeSearch	intron	3934298	3935799	.	+	.	gene_id=FBA8_1
1	PaGeSearch	intron	3936070	3936187	.	+	.	gene_id=FBA8_1
1	PaGeSearch	intron	3936298	3936499	.	+	.	gene_id=FBA8_1
1	PaGeSearch	intron	3936591	3937221	.	+	.	gene_id=FBA8_1
1	PaGeSearch	intron	3937492	3938024	.	+	.	gene_id=FBA8_1
1	PaGeSearch	CDS	3934136	3934297	.	+	.	gene_id=FBA8_1
1	PaGeSearch	CDS	3935800	3936069	.	+	.	gene_id=FBA8_1
1	PaGeSearch	CDS	3936188	3936297	.	+	.	gene_id=FBA8_1
1	PaGeSearch	CDS	3936500	3936590	.	+	.	gene_id=FBA8_1
1	PaGeSearch	CDS	3937222	3937491	.	+	.	gene_id=FBA8_1
1	PaGeSearch	CDS	3938025	3938303	.	+	.	gene_id=FBA8_1
1	PaGeSearch	stop_codon	3938301	3938303	.	+	.	gene_id=FBA8_1
```


### Download gene sequences 
You can download the seqeucnes of genes in a pathway or a list, that can directly be used as the query of PaGeSearch.
The pathway information is based on the Reactome database, and the pathway name must be entered exactly as it appears in the Reactome database.
Refer to https://reactome.org/ for animal species and https://plantreactome.gramene.org/index.php?lang=en for plant species.

#### Download gene sequences from pathway name
Download sequences of genes in the pathway 'Metabolism of nucleotides' of Homo sapiens and save to a folder named metabolism_of_nucleotides_gene_sequences.
The sequences of each gene is saved to separate fasta files.
```
python Codes/download_gene_sequences.py -p "Metabolism of nucleotides" -s human -t 9606 -o ../metabolism_of_nucleotides_gene_sequences
```
You can run PaGeSearch using the downloaded gene sequences as:
```
python Codes/pagesearch.py -g ../path/to/genome.fa -p ../metabolism_of_nucleotides_gene_sequences -od ../pagesearch_results -op metabolism_of_nucleotides -s human -t 4
```
Download all sequences that are orthologs of the genes in the humen 'Metabolism of nucleotides' pathway within mammali species.
The taxonomy ID is based on NCBI taxonomy (https://www.ncbi.nlm.nih.gov/taxonomy).
```
python Codes/download_gene_sequences.py -p "Metabolism of nucleotides" -s human -t 40674 -o metabolism_of_nucleotides_mammals_orthologs
```
#### Download gene sequences from a user defined list
The gene list is a text file of Ensembl gene IDs separated by line breaks.
```
python Codes/download_gene_sequences.py -l genelist.txt -s human -t 9606 -o gene_sequences
```

### Example run
Test data and code is available at the Github 'example' folder. 


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
