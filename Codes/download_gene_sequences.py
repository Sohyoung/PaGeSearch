import http.client
import json
import os
import sys
import time
import argparse



def download_gene_sequences(gene_id, taxon_id, outdir):
    # Get the gene symbol from the Ensembl API
    gene_url = f"/lookup/id/{gene_id}?format=full"
    conn = http.client.HTTPSConnection("rest.ensembl.org")
    conn.request("GET", gene_url, headers={"Content-Type": "application/json"})
    gene_response = conn.getresponse()
    gene_data = json.loads(gene_response.read().decode("utf-8"))
    gene_name = gene_data.get("display_name", gene_id)

    # Construct the Ensembl API request URL for the gene's sequence
    gene_sequence_url = f"/sequence/id/{gene_id}?type=protein;multiple_sequences=1"

    # Send the API request to download the gene sequence
    conn.request("GET", gene_sequence_url, headers={"Content-Type": "text/plain"})
    gene_sequence_response = conn.getresponse()
    gene_sequence = gene_sequence_response.read().decode("utf-8") if gene_sequence_response.status in range(200, 300) else None

    # Write the gene's sequences to a file
    if gene_sequence is not None:
        gene_sequence_list = gene_sequence.split('\n')
        if ' ' in gene_name or '/' in gene_name or ';' in gene_name:
            gene_file_path = os.path.join(outdir, f"{gene_name.replace(' ', '_').replace('/', '_').replace(';', '_')}.fa")
        else:
            gene_file_path = os.path.join(outdir, f"{gene_name}.fa")
        print(f"Write sequence of {gene_id}")
        if len(gene_sequence_list) > 1:
            for i in range(len(gene_sequence_list) - 1):
                with open(gene_file_path, 'a') as outfile:
                    outfile.write(f">{gene_id}_{str(i)}\n{gene_sequence_list[i]}\n")
        else:
            with open(gene_file_path, 'a') as outfile:
                outfile.write(f">{gene_id}\n{gene_sequence_list[0]}\n")

    # Construct the Ensembl API request URL for the gene's orthologs
    ortholog_url = f"/homology/id/{gene_id}?type=orthologues;target_taxon={taxon_id}"
    # Send the API request to download the orthologs
    try:
        conn = http.client.HTTPSConnection("rest.ensembl.org")
        conn.request("GET", ortholog_url, headers={"Content-Type": "application/json"})
        resp = conn.getresponse()
        if resp.status == 200:
            # Parse the API response as JSON data
            data = json.loads(resp.read().decode("utf-8"))
            if data.get("data") and len(data["data"]) > 0 and len(data["data"][0].get("homologies", [])) > 0:
                with open(outdir + '/download.log', 'a') as f:
                    f.write(gene_id + '\t' + gene_name + '\n')
                # Loop through each ortholog and download its sequence
                for ortholog in data["data"][0]["homologies"]:
                    # Construct the Ensembl API request URL for the ortholog's sequence
                    sequence_url = f"/sequence/id/{ortholog['target']['id']}?type=protein;multiple_sequences=1"
                    # Send the API request to download the ortholog's sequence
                    conn.request("GET", sequence_url, headers={"Content-Type": "text/plain"})
                    sequence_resp = conn.getresponse()
                    if sequence_resp.status == 200:
                        # Write the ortholog gene name and sequence to the file
                        sequence_list = sequence_resp.read().decode("utf-8").split('\n')
                        if ' ' in gene_name or '/' in gene_name or ';' in gene_name:
                            gene_file_path = os.path.join(outdir, f"{gene_name.replace(' ', '_').replace('/', '_').replace(';', '_')}.fa")
                        else:
                            gene_file_path = os.path.join(outdir, f"{gene_name}.fa")
                        if (len(sequence_list)) > 1:
                            for i in range(len(sequence_list) - 1):
                                with open(gene_file_path, 'a') as outfile:
                                    outfile.write(f">{ortholog['target']['id']}_{str(i)}\n{sequence_list[i]}\n")
                        else:
                            with open(gene_file_path, 'a') as outfile:
                                outfile.write(f">{ortholog['target']['id']}\n{sequence_list[0]}\n")
                # Print the success message to the console
                print(f"All sequences downloaded and written to {os.path.join(outdir, f'{gene_name}.fa')}")
    except Exception as e:
        print(f"Failed to retrieve orthologs for {gene_id}: {e}")
    finally:
        conn.close()
    return gene_name

def download_pathway_genes(gene_pathway_file_path, pathway_name, species, taxon_id, outdir):
    if os.path.exists(outdir) == False:
        os.mkdir(outdir)
    with open(f'{outdir}/download.log', 'w') as f:
        f.write(f'GeneName\tGeneID\n')
    # Define lists to store the Gene IDs
    gene_list = []
    # Open the file and read the contents
    with open(gene_pathway_file_path, 'r') as file:
        # Skip the header row
        next(file)
        # Iterate through each line in the file
        for line in file:
            # Split the line into columns
            columns = line.strip().split('\t')
            # Check if the Reactome ID of the row matches the Reactome ID of interest
            if columns[3] == pathway_name:
                gene_id = columns[0]
                # If gene_id is not in list, download gene sequences
                if not gene_id in gene_list:
                    gene_list.append(gene_id)
                    print("Download orthologs of " + gene_id)
                    gene_name = download_gene_sequences(gene_id, taxon_id, outdir)
                    with open(f'{outdir}/download.log', 'a') as f:
                        f.write(f'{gene_name}\t{gene_id}\n')


def download_genelist_genes(genelist_fname, species, taxon_id, outdir):
    # Define lists to store the Gene IDs
    gene_list = []
    # Open the file and read the contents
    with open(genelist_fname, 'r') as file:
        # Skip the header row
        next(file)
        # Iterate through each line in the file
        for line in file:
            gene_list.append(line.split('\t')[0].strip())
    for gene_id in gene_list:
        print("Download orthologs of " + gene_id)
        download_gene_sequences(gene_id, taxon_id, outdir)


def main(argv):
    stime = time.time()
    parser = argparse.ArgumentParser(prog='PROG', usage='%(prog)s [options]')
    parser.add_argument('-p', '--pathway', action="store", default='')
    parser.add_argument('-l', '--genelist', action="store", default='')
    parser.add_argument('-s', '--species', dest='species', type=str, action="store", default='human')
    parser.add_argument('-t', '--taxon', dest='taxonid', type=str, action="store", default='9606')
    parser.add_argument('-o', '--outdir', dest='outdir', type=str, action="store", default='./pagesearch_query')

    args = parser.parse_args()

    if os.path.exists(args.outdir) == False:
        os.mkdir(args.outdir)

    if args.pathway != '':
        gene_pathway_file_path = f"Reactome/{args.species}.txt"
        download_pathway_genes(gene_pathway_file_path, args.pathway, args.species, args.taxonid, args.outdir)
    if args.genelist != '':
        download_genelist_genes(args.genelist, args.species, args.taxonid, args.outdir)

    print('\nGene sequence download finished.')

main(sys.argv)
