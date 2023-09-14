import http.client
import json
import os


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
                    else:
                        print(
                            f"Failed to download sequence for {ortholog['target']['id']}: {sequence_resp.status} {sequence_resp.reason}")
                # Print the success message to the console
                print(f"All sequences downloaded and written to {os.path.join(outdir, f'{gene_name}.fa')}")
            else:
                print(f"Failed to retrieve orthologs for {gene_id}")
        else:
            print(f"Failed to retrieve orthologs for {gene_id}: {resp.status} {resp.reason}")
    except Exception as e:
        print(f"Failed to retrieve orthologs for {gene_id}: {e}")
    finally:
        conn.close()


def download_pathway_genes(gene_pathway_file_path, pathway_name, species, taxon_id, outdir):
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
            if columns[5] == species:
                if columns[3] == pathway_name:
                    gene_id = columns[0]
                    # If gene_id is not in list, download gene sequences
                    if not gene_id in gene_list:
                        print("Download orthologs of " + gene_id)
                        download_gene_sequences(gene_id, taxon_id, outdir)
                        gene_list.append(gene_id)


def download_genelist_genes(genelist_fname, species, taxon_id, outdir):
    # Define lists to store the Gene IDs
    gene_list = []
    # Open the file and read the contents
    with open(genelist_fname, 'r') as file:
        # Skip the header row
        next(file)
        # Iterate through each line in the file
        for line in file:
            gene_list.append(line.split('\t')[0])
    for gene_id in gene_list:
        print("Download orthologs of " + gene_id)
        download_gene_sequences(gene_id, taxon_id, outdir)



def main():
    # Open the file for reading
    with open('/disk1/1.Sohyoung_Pipeline/PathwayTest/Testset_taxons.txt', 'r') as f:
        # Skip the header row
        next(f)
        # Initialize an empty dictionary to store the parsed information
        taxon_dic = {}
        # Loop through the lines in the file
        for line in f:
            # Split the line into columns based on the tab delimiter
            columns = line.strip().split('\t')
            # Get the scientific name of species 1 and the NCBI Taxon ID of the largest taxon
            scientific_name = columns[1]
            largest_taxon_id = columns[8]
            # Add the information to the dictionary using species 1 as the key
            taxon_dic[scientific_name] = largest_taxon_id

    '''
    for pathway_tmp in pathway_dic.keys():
        for species, taxon_id in taxon_dic.items():
            if species == 'Arabidopsis thaliana' or species == 'Triticum aestivum':
                pathway = pathway_dic[pathway_tmp].split(' (')[0].replace(' ', '_')
                gene_pathway_file_path = '/disk1/1.Sohyoung_Pipeline/Data/Ensembl/' + species.replace(' ',  '_') + '_reactome_genes_only.txt'
                if not os.path.exists('/disk1/1.Sohyoung_Pipeline/PathwayTest/' + species.replace(' ', '_')):
                    os.makedirs('/disk1/1.Sohyoung_Pipeline/PathwayTest/' + species.replace(' ', '_'))
                if not os.path.exists('/disk1/1.Sohyoung_Pipeline/PathwayTest/' + species.replace(' ', '_') + '/Data/'):
                    os.makedirs('/disk1/1.Sohyoung_Pipeline/PathwayTest/' + species.replace(' ', '_') + '/Data/')
                outdir = '/disk1/1.Sohyoung_Pipeline/PathwayTest/' + species.replace(' ', '_') + '/Data/' + pathway_tmp.split(' (')[0].replace(' ', '_')
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                else:
                    os.system('rm ' + outdir + '/*')
                os.chdir(outdir)
                print(species + '\t' + pathway)
                download_pathway_genes(gene_pathway_file_path, pathway_dic[pathway_tmp], species, taxon_id, outdir)

            else:
                pathway = pathway_tmp.split(' (')[0].replace(' ', '_')
                gene_pathway_file_path = '/disk1/1.Sohyoung_Pipeline/Data/Ensembl/' + species.replace(' ',  '_') + '_reactome_genes_only.txt'
                if not os.path.exists('/disk1/1.Sohyoung_Pipeline/PathwayTest/' + species.replace(' ', '_')):
                    os.makedirs('/disk1/1.Sohyoung_Pipeline/PathwayTest/' + species.replace(' ', '_'))
                if not os.path.exists('/disk1/1.Sohyoung_Pipeline/PathwayTest/' + species.replace(' ', '_') + '/Data/'):
                    os.makedirs('/disk1/1.Sohyoung_Pipeline/PathwayTest/' + species.replace(' ', '_') + '/Data/')
                outdir = '/disk1/1.Sohyoung_Pipeline/PathwayTest/' + species.replace(' ', '_') + '/Data/' +  pathway_tmp.split(' (')[0].replace(' ', '_')
                if not os.path.exists(outdir):
                    os.makedirs(outdir)
                else:
                    os.system('rm ' + outdir + '/*')
                os.chdir(outdir)
                print(species + '\t' + pathway)
                download_pathway_genes(gene_pathway_file_path, pathway_tmp, species, taxon_id, outdir)
    '''
    species = 'Arabidopsis thaliana'
    taxon_id = taxon_dic[species]
    pathway = "TCA cycle (plant)"
    gene_pathway_file_path = '/disk1/1.Sohyoung_Pipeline/Data/Ensembl/' + species.replace(' ', '_') + '_reactome_genes_only.txt'
    outdir = '/disk1/1.Sohyoung_Pipeline/PathwayTest/' + species.replace(' ', '_') + '/Data/' + "Citric acid cycle".replace(' ', '_')
    download_pathway_genes(gene_pathway_file_path, pathway, species, taxon_id, outdir)

if __name__ == '__main__':
    main()
    