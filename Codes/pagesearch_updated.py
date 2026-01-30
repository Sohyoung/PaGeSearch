import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import subprocess
import pandas as pd
import time
import glob
import sys
import getopt
from subprocess import DEVNULL, STDOUT, check_call
import argparse
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
import math
from tensorflow import keras


# -------------------------
# Utilities
# -------------------------
def sanitize_fasta_ambiguous_to_X(in_fa, out_fa=None, ambiguous_set=set("BJOUZOU")):
    """
    Replace ambiguous/rare amino acid letters in FASTA sequences with 'X'.
    Only modifies sequence lines, leaves headers unchanged.

    Returns:
        (n_seqs, n_seqs_modified, n_chars_modified)
    """
    if out_fa is None:
        out_fa = in_fa

    n_seqs = 0
    n_seqs_modified = 0
    n_chars_modified = 0

    with open(in_fa, 'r') as fin, open(out_fa + ".san.tmp", 'w') as fout:
        header = None
        seq_chunks = []
        seq_modified_flag = False

        def flush_record(h, chunks):
            nonlocal n_seqs, n_seqs_modified, n_chars_modified
            if h is None:
                return
            n_seqs += 1
            seq = "".join(chunks)
            if not seq:
                # allow empty; downstream will catch if file becomes empty
                fout.write(h)
                fout.write("\n")
                return
            # Replace ambiguous letters
            original = seq
            # Count changes efficiently
            # Do a simple loop to avoid regex overhead
            seq_list = list(seq)
            changed = 0
            for i, ch in enumerate(seq_list):
                if ch in ambiguous_set:
                    seq_list[i] = 'X'
                    changed += 1
            if changed > 0:
                n_seqs_modified += 1
                n_chars_modified += changed
            sanitized = "".join(seq_list)

            fout.write(h + "\n")
            # Write sequence in same style (single line here; safe for downstream)
            fout.write(sanitized + "\n")

        for line in fin:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # flush previous
                flush_record(header, seq_chunks)
                header = line
                seq_chunks = []
            else:
                if line.strip() != "":
                    seq_chunks.append(line.strip())

        # flush last
        flush_record(header, seq_chunks)

    os.replace(out_fa + ".san.tmp", out_fa)
    return n_seqs, n_seqs_modified, n_chars_modified


# 1. Seed sequence search by MMseqs2
def mmseqs(genome, db_dir, pathway_gene_dir, mmseqs_results, nthreads, sanitize_ambiguous=True):
    """
    Build query fasta from pathway_gene_dir/*.fa, optionally sanitize ambiguous letters -> X,
    deduplicate, then run mmseqs search against genome db.

    sanitize_ambiguous=True:
        Replace BJOUZOU with X in query sequences (recommended).
    """
    genome_db = db_dir + '/genome'
    print('\nCandidate region identification through sequence similarity search')

    # Log file for mmseqs (stdout+stderr) (B)
    mmseqs_log = mmseqs_results.replace('.txt', '.mmseqs.log')

    def run_logged(cmd_list):
        with open(mmseqs_log, 'a') as log:
            log.write("\n# CMD: " + " ".join(cmd_list) + "\n")
            subprocess.check_call(cmd_list, stdout=log, stderr=log)

    # Make MMseqs2 genome database
    if os.path.exists(genome + '.dbtype'):
        os.system('rm ' + genome + '.dbtype*')
    run_logged(['mmseqs', 'createdb', genome, genome_db])

    # Prepare paths
    tmp_fa = mmseqs_results.replace('.txt', '.tmp.fa')
    query_fa = mmseqs_results.replace('.txt', '.fa')
    uniq_fa = mmseqs_results.replace('.txt', '_unique.fa')
    query_db_prefix = mmseqs_results.replace('.txt', '')

    # Clean old intermediates
    for fp in [tmp_fa, query_fa, uniq_fa]:
        if os.path.exists(fp):
            os.remove(fp)

    # Ensure there are .fa files
    fa_list = glob.glob(os.path.join(pathway_gene_dir, '*.fa'))
    if len(fa_list) == 0:
        raise RuntimeError(f"No .fa files found in pathway_gene_dir: {pathway_gene_dir}")

    # Concatenate
    os.system('cat ' + pathway_gene_dir + '/*.fa >> ' + tmp_fa)

    # Convert to 2-line FASTA without deleting sequences (no grep deletion)
    # (Keeps headers+sequence joined; downstream tools handle 2-line fine.)
    line = (
        "perl -pe '$. > 1 and /^>/ ? print " + '"\\n"' + " : chomp' " + tmp_fa +
        ' | paste - - | tr "\\t" "\\n" > ' + query_fa
    )
    os.system(line)

    if os.path.exists(tmp_fa):
        os.remove(tmp_fa)

    # C) Sanitize ambiguous letters -> X (sequence lines only)
    if sanitize_ambiguous:
        n_seqs, n_mod_seqs, n_mod_chars = sanitize_fasta_ambiguous_to_X(query_fa)
        with open(mmseqs_log, 'a') as log:
            log.write(
                f"\n# SANITIZE: replaced BJOUZOU -> X in {query_fa}\n"
                f"# SANITIZE: sequences={n_seqs}, modified_sequences={n_mod_seqs}, modified_chars={n_mod_chars}\n"
            )

    # Deduplicate sequences by sequence content
    line = 'seqkit rmdup -s ' + query_fa + ' 2> tmp.log | cat > ' + uniq_fa
    os.system(line)

    # A) Guard: unique fasta must exist and be non-empty
    if (not os.path.exists(uniq_fa)) or (os.path.getsize(uniq_fa) == 0):
        log_tail = ""
        if os.path.exists('tmp.log'):
            try:
                with open('tmp.log', 'r') as f:
                    log_tail = "".join(f.readlines()[-80:])
            except:
                pass
        raise RuntimeError(
            f"Empty or missing query FASTA after dedup: {uniq_fa}\n"
            f"- Check query/*.fa exists and is valid FASTA.\n"
            f"- Check seqkit rmdup stderr (tmp.log tail below).\n\n"
            f"{log_tail}"
        )

    # Create MMseqs2 query database (B logged)
    run_logged(['mmseqs', 'createdb', uniq_fa, query_db_prefix])

    # Remove query_fa (keep uniq_fa for debugging; comment out if you want to delete)
    if os.path.exists(query_fa):
        os.remove(query_fa)

    # Run MMseqs2 search
    tmp_dir = db_dir + '/tmp'
    os.makedirs(tmp_dir, exist_ok=True)

    run_logged([
        'mmseqs', 'search',
        query_db_prefix, genome_db, f'{mmseqs_results}.m8', tmp_dir,
        '-a', '--start-sens', '1', '--sens-steps', '3', '-s', '10',
        '--threads', str(nthreads), '--alignment-mode', '1'
    ])

    run_logged([
        'mmseqs', 'convertalis',
        query_db_prefix, genome_db, f'{mmseqs_results}.m8', mmseqs_results,
        '--threads', str(nthreads),
        '--format-output', 'query,target,qstart,qend,tstart,tend,alnlen,raw,nident,qlen'
    ])

    os.system('rm -r ' + tmp_dir)


# 2. Extended from mmseqs seed regions and extract the sequences
# Merge nearby hsps in the same strand
def merge_results(grouped_data):
    _, group = grouped_data  # Unpack the tuple
    merged_regions = []
    current_region = group.iloc[0].copy()
    for idx, row in group.iloc[1:].iterrows():
        if abs(current_region['end'] - row['start']) <= 100000:
            current_region['end'] = max(current_region['end'], row['end'])
        else:
            merged_regions.append(current_region)
            current_region = row.copy()
    merged_regions.append(current_region)
    return pd.DataFrame(merged_regions)

# Merge overlapping hsp regions
def find_largest_overlapping_region(grouped_data):
    gene, group = grouped_data  # Unpack the tuple
    group.sort_values(['chr', 'start'], inplace=True)
    merged_regions = []
    current_region = group.iloc[0].copy()
    for _, row in group.iloc[1:].iterrows():
        if row['start'] <= current_region['end'] and row['geneid'] != current_region['geneid']:
            current_region['end'] = max(current_region['end'], row['end'])
        else:
            merged_regions.append(current_region.copy())
            current_region = row.copy()
    merged_regions.append(current_region)
    return pd.DataFrame(merged_regions)

# Write bed and sequence files
def process_group(task_params):
    group, gene, genome, seq_extract_dir, pathway_gene_dir = task_params
    group['chr_start_end'] = group['chr'].astype(str) + ':' + group['start'].astype(str) + '-' + group['end'].astype(str)
    bed_dir = os.path.join(seq_extract_dir, 'bed')
    os.makedirs(bed_dir, exist_ok=True)
    bed_filename = os.path.join(bed_dir, f"{gene}.bed")
    group[['chr', 'start', 'end', 'chr_start_end']].to_csv(bed_filename, sep='\t', index=False, header=False)

    seq_dir = os.path.join(seq_extract_dir, 'sequence_fasta')
    os.makedirs(seq_dir, exist_ok=True)
    os.system(f'bedtools getfasta -fi {genome} -bed {bed_filename} -fo {os.path.join(seq_dir, f"{gene}.fa")} > /dev/null')

    id_dir = os.path.join(seq_extract_dir, 'gene_fasta')
    os.makedirs(id_dir, exist_ok=True)
    geneid_filename = os.path.join(id_dir, f"{gene}_id.txt")
    group['geneid'].drop_duplicates().to_csv(geneid_filename, index=False, header=False)

    gene_fasta = os.path.join(pathway_gene_dir, f"{gene}.fa")
    with open(geneid_filename) as f:
        gene_ids = [line.strip() for line in f]

    os.system(f'samtools faidx {gene_fasta} {" ".join(gene_ids)} > {geneid_filename.replace("_id.txt", ".fa")}')

# Parse mmseqs2 dataframe and use multiprocessing to extend seeds and extract sequences
def seq_extract(genome, mmseqs_results, seq_extract_dir, nextension, chr_len_fname, pathway_gene_dir, gene_id_fname, nthreads):
    nthreads = int(nthreads)
    columns = ['geneid', 'chr', 'qstart', 'qend', 'start', 'end', 'alnlen', 'score', 'nident', 'qlen']
    data = pd.read_csv(mmseqs_results, sep='\t', header=None, names=columns, dtype={'chr': 'str'})
    geneid_dic = pd.read_csv(gene_id_fname, sep='\t', header=None, names=['gene', 'geneid'])
    with open(chr_len_fname, 'r') as file:
        chrlen_dic = {line.split()[0]: int(line.split()[1].strip()) for line in file}

    data = data.merge(geneid_dic, on='geneid', how='left')
    data['strand'] = data.apply(lambda row: '+' if row['start'] < row['end'] else '-', axis=1)
    data.loc[data['strand'] == '-', ['start', 'end']] = data.loc[data['strand'] == '-', ['end', 'start']].values
    data = data.sort_values(by=['gene', 'geneid', 'chr', 'strand', 'start'])
    grouped = list(data.groupby(['gene', 'geneid', 'chr', 'strand']))

    with ProcessPoolExecutor(max_workers=nthreads) as executor:
        result = executor.map(merge_results, grouped)
    merged_df = pd.concat(result)
    merged_df = merged_df.sort_values(by=['gene', 'chr', 'start'])

    grouped_df = merged_df.groupby('gene')
    tasks = [(gene, group) for gene, group in grouped_df]

    with ProcessPoolExecutor(max_workers=nthreads) as executor:
        result = executor.map(find_largest_overlapping_region, tasks)
    overlapping_regions_df = pd.concat(result)

    overlapping_regions_df['start'] = merged_df['start'].apply(lambda x: max(0, x - int(nextension)))

    def adjust_end(row, chrlen_dic):
        return min(chrlen_dic[str(row['chr'])], row['end'] + int(nextension))

    overlapping_regions_df['end'] = overlapping_regions_df.apply(lambda row: adjust_end(row, chrlen_dic), axis=1)

    grouped_df = overlapping_regions_df.groupby('gene')
    tasks = [(group, gene, genome, seq_extract_dir, pathway_gene_dir) for gene, group in grouped_df]

    with ProcessPoolExecutor(max_workers=nthreads) as executor:
        executor.map(process_group, tasks)


# 3. Gene prediction
def exonerate_hints(gene_fasta, fasta_extended, exonerate_hints_gff):
    line = 'exonerate ' + gene_fasta + ' ' + fasta_extended + ' --model protein2genome --showtargetgff T > ' + exonerate_hints_gff.replace('.gff', '_tmp.txt')
    os.system(line)
    line = 'exonerate2hints.pl --in=' + exonerate_hints_gff.replace('.gff', '_tmp.txt') + ' --source=M --out=' + exonerate_hints_gff
    os.system(line)
    os.remove(exonerate_hints_gff.replace('.gff', '_tmp.txt'))

def augustus(fasta_extended, exonerate_hints_gff, gene_prediction_gff, species_model):
    if os.path.exists(exonerate_hints_gff) and os.path.getsize(exonerate_hints_gff):
        line = 'augustus --genemodel=complete --noInFrameStop=true --UTR=off --species=' + species_model + ' --hintsfile=' + exonerate_hints_gff + ' ' + fasta_extended + ' > ' + gene_prediction_gff
        os.system(line)
    else:
        line = 'augustus --genemodel=complete --noInFrameStop=true --UTR=off --species=' + species_model + ' ' + fasta_extended + ' > ' + gene_prediction_gff
        os.system(line)
    line = 'getAnnoFasta.pl ' + gene_prediction_gff
    os.system(line)

def exonerate_hints_task(params):
    gene_fasta_dir, fasta_extended, exonerate_hints_dir = params
    genename = os.path.basename(fasta_extended).replace('.fa', '')
    gene_fasta = os.path.join(gene_fasta_dir, genename + '.fa')
    exonerate_hints_gff = os.path.join(exonerate_hints_dir, genename + '_hints.gff')
    exonerate_hints(gene_fasta, fasta_extended, exonerate_hints_gff)

def augustus_task(params):
    fasta_extended, exonerate_hints_dir, gene_prediction_dir, species_model = params
    genename = os.path.basename(fasta_extended).replace('.fa', '')
    exonerate_hints_gff = os.path.join(exonerate_hints_dir, genename + '_hints.gff')
    gene_prediction_gff = os.path.join(gene_prediction_dir, genename + '.gff')
    augustus(fasta_extended, exonerate_hints_gff, gene_prediction_gff, species_model)

def gene_prediction(gene_fasta_dir, fasta_extended_dir, exonerate_hints_dir, gene_prediction_dir, species_model, nthreads):
    nthreads = int(nthreads)
    fasta_extended_list = glob.glob(fasta_extended_dir + '/*.fa')
    print('Making protein hints')
    with ProcessPoolExecutor(max_workers=nthreads) as executor:
        params = [(gene_fasta_dir, fasta_extended, exonerate_hints_dir) for fasta_extended in fasta_extended_list]
        executor.map(exonerate_hints_task, params)
    print('Ab initio gene prediction')
    with ProcessPoolExecutor(max_workers=nthreads) as executor:
        params = [(fasta_extended, exonerate_hints_dir, gene_prediction_dir, species_model) for fasta_extended in fasta_extended_list]
        executor.map(augustus_task, params)


# 4. Align reference gene sequences to the gene predicted sequences using exonerate
def exonerate(gene_fasta, fasta_genepredicted, protein_alignment_output):
    line = 'exonerate --target ' + gene_fasta + ' --query ' + fasta_genepredicted + ' --querytype protein --targettype protein --model ungapped --bestn 1 --showalignment no --showvulgar no --ryo "%qi\t%ti\t%s\t%ql\t%tl\t%et\t%ei\t%em\t%pi\t%ps\t%s\n" > ' + protein_alignment_output
    os.system(line)

def exonerate_task(params):
    fasta_genepredicted, pathway_gene_dir, exonerate_dir = params
    genename = os.path.basename(fasta_genepredicted).replace('.aa', '')
    gene_fasta = os.path.join(pathway_gene_dir, genename + '.fa')
    line = "perl -pe '$. > 1 and /^>/ ? print " + '"\\n"' + " : chomp' " + gene_fasta + ' | paste - - | grep -P -v ".+\\t.*[BJOUZ]+.*" | tr "\\t" "\\n" > ' + gene_fasta.replace('.fa', '.tmp.fa')
    os.system(line)
    protein_alignment_output = os.path.join(exonerate_dir, genename + '.txt')
    exonerate(gene_fasta.replace('.fa', '.tmp.fa'), fasta_genepredicted, protein_alignment_output)
    os.remove(gene_fasta.replace('.fa', '.tmp.fa'))

def run_exonerate_multithreading(fasta_genepredicted_dir, pathway_gene_dir, exonerate_dir, nthreads):
    nthreads = int(nthreads)
    fasta_genepredicted_list = glob.glob(fasta_genepredicted_dir + '/*.aa')
    with ProcessPoolExecutor(max_workers=nthreads) as executor:
        params = [(fasta_genepredicted, pathway_gene_dir, exonerate_dir) for fasta_genepredicted in fasta_genepredicted_list]
        executor.map(exonerate_task, params)


# 5. Parse and filter results
def parse_and_merge_onegene(genename, gene_prediction_gff, exonerate_results, summary_fname):
    data = []
    with open(gene_prediction_gff, 'r') as gff_file:
        for line in gff_file:
            if line.startswith('# start gene'):
                fields = next(gff_file).split('\t')
                chromosome = fields[0].split(':')[0]
                chr_start_end = fields[0]
                gene_start = int(fields[0].split(':')[1].split('-')[0]) + int(fields[3])
                gene_end = int(fields[0].split(':')[1].split('-')[0]) + int(fields[4])
                transcript_id = next(gff_file).split('\t')[8].strip()
            elif line.startswith('# % of transcript supported by hints'):
                hints_percent = float(line.split(':')[1].strip())
                data.append([chromosome, gene_start, gene_end, transcript_id, hints_percent])

    augustus_df = pd.DataFrame(data, columns=['Chromosome', 'Gene_Start', 'Gene_End', 'Transcript_ID', 'Percent_of_Transcript_Supported_by_Hints'])

    exonerate_df = pd.read_csv(
        exonerate_results, sep='\t', skiprows=3, header=None, skipfooter=1, engine='python', index_col=None,
        names=['Transcript_ID', 'Gene_ID_exonerate', 'Score', 'Sequence_Length_exonerate', 'Gene_Length_exonerate',
               'Alignment_Length_exonerate', 'Identity_Length', 'Match_Length', 'Percent_Identity', 'Percent_Similarity', 'Raw_Score'],
        dtype={'Sequence_Length': int, 'Gene_Length': int, 'Alignment_Length': int, 'Identity_Length': int, 'Match_Length': int, 'Raw_Score': int}
    )

    combined_df = pd.merge(exonerate_df, augustus_df, on='Transcript_ID', how='inner')
    combined_df['Gene_Cover_exonerate'] = combined_df['Alignment_Length_exonerate'] / combined_df['Gene_Length_exonerate']
    combined_df['Sequence_Cover_exonerate'] = combined_df['Alignment_Length_exonerate'] / combined_df['Sequence_Length_exonerate']
    combined_df['Normalized_Score_exonerate'] = combined_df['Raw_Score'] / combined_df['Alignment_Length_exonerate']
    combined_df['Percent_Identity_exonerate'] = combined_df['Percent_Identity'] / 100.0
    combined_df['Percent_Similarity_exonerate'] = combined_df['Percent_Similarity'] / 100.0
    combined_df['Percent_Transcript_Supported_by_Hints'] = combined_df['Percent_of_Transcript_Supported_by_Hints'] / 100.0
    combined_df['Gene'] = genename
    combined_df = combined_df.drop('Percent_Identity', axis=1)

    df_out = combined_df[['Gene', 'Transcript_ID', 'Gene_ID_exonerate', 'Percent_of_Transcript_Supported_by_Hints',
                          'Percent_Identity_exonerate', 'Percent_Similarity_exonerate', 'Gene_Cover_exonerate',
                          'Sequence_Cover_exonerate', 'Normalized_Score_exonerate', 'Chromosome', 'Gene_Start', 'Gene_End']]
    df_out.to_csv(summary_fname.replace('.txt', '.tmp.txt'), sep='\t', index=False, mode='a', header=False)

def parse_and_merge_onegene_wrapper(args):
    genename, gene_prediction_gff, exonerate_results, summary_fname = args
    try:
        parse_and_merge_onegene(genename, gene_prediction_gff, exonerate_results, summary_fname)
    except:
        pass

def filter_and_calculate(params):
    filtered_regions = []
    gene, group, df2 = params
    for _, row in group.iterrows():
        chromosome = row['Chromosome']
        gene_start = row['Gene_Start']
        gene_end = row['Gene_End']
        df2_tmp = df2[(df2['chr'] == chromosome)]
        overlapping_regions = df2_tmp[(df2_tmp['end'] >= gene_start) & (df2_tmp['start'] <= gene_end)]
        grouped_overlap_regions = list(overlapping_regions.groupby(['geneid']))
        if len(grouped_overlap_regions) > 0:
            for geneid, overlap_bygeneid in grouped_overlap_regions:
                alnset = set()
                pident_sum = 0
                alnlen_sum = 0
                for _, o_row in overlap_bygeneid.iterrows():
                    if o_row['strand'] == '+':
                        astart = o_row['qstart']
                        aend = o_row['qend']
                        start = o_row['start']
                        end = o_row['end']
                        if gene_start > start:
                            astart = o_row['qstart'] + math.ceil((gene_start - start) / 3)
                        if gene_end < end:
                            aend = o_row['qend'] - math.ceil((end - gene_end) / 3)
                        o_row['alnset'] = set(range(astart, aend + 1))
                        alnset = alnset.union(o_row['alnset'])
                        pident_sum += o_row['Percent_Identity_mmseqs'] * o_row['alnlen']
                        alnlen_sum += o_row['alnlen']
                    if o_row['strand'] == '-':
                        astart = o_row['qstart']
                        aend = o_row['qend']
                        start = o_row['start']
                        end = o_row['end']
                        if gene_start > start:
                            aend = o_row['qend'] - math.ceil((gene_start - start) / 3)
                        if gene_end < end:
                            astart = o_row['qstart'] + math.ceil((end - gene_end) / 3)
                        o_row['alnset'] = set(range(astart, aend + 1))
                        alnset = alnset.union(o_row['alnset'])
                        pident_sum += o_row['Percent_Identity_mmseqs'] * o_row['alnlen']
                        alnlen_sum += o_row['alnlen']
                pident = pident_sum / alnlen_sum
                filtered_regions.append({
                    'Gene': gene,
                    'Transcript_ID': row['Transcript_ID'],
                    'Gene_ID': geneid,
                    'Percent_of_Transcript_Supported_by_Hints': row['Percent_of_Transcript_Supported_by_Hints'],
                    'Percent_Identity_exonerate': row['Percent_Identity_exonerate'],
                    'Percent_Similarity_exonerate': row['Percent_Similarity_exonerate'],
                    'Gene_Cover_exonerate': row['Gene_Cover_exonerate'],
                    'Sequence_Cover_exonerate': row['Sequence_Cover_exonerate'],
                    'Normalized_Score_exonerate': row['Normalized_Score_exonerate'],
                    'Chromosome': chromosome,
                    'Gene_Start': gene_start,
                    'Gene_End': gene_end,
                    'alnset': alnset,
                    'Percent_Identity_mmseqs': pident
                })

    filtered_df = pd.DataFrame(filtered_regions)
    if filtered_df.shape[0] > 0:
        df2 = df2.drop_duplicates(subset=['geneid'])
        filtered_df['qlen'] = filtered_df['Gene_ID'].map(df2.set_index('geneid')['qlen'])
        filtered_df['Gene_Cover_mmseqs'] = filtered_df['alnset'].apply(len) / filtered_df['qlen']
        filtered_df['Gene_Transcript_ID'] = filtered_df.apply(lambda row: f"{row['Gene']}_{row['Transcript_ID']}", axis=1)
        filtered_df.sort_values(['Gene_Cover_mmseqs', 'Percent_Identity_mmseqs'], ascending=[False, False], inplace=True)
        filtered_df.drop_duplicates(subset=['Gene_Transcript_ID'], keep='first', inplace=True)
        filtered_df.set_index('Gene_Transcript_ID', inplace=True)
        outdf = filtered_df[['Gene', 'Transcript_ID', 'Percent_of_Transcript_Supported_by_Hints', 'Percent_Identity_exonerate',
                             'Percent_Similarity_exonerate', 'Gene_Cover_exonerate', 'Sequence_Cover_exonerate',
                             'Normalized_Score_exonerate', 'Gene_Cover_mmseqs', 'Percent_Identity_mmseqs',
                             'Chromosome', 'Gene_Start', 'Gene_End']]
    else:
        outdf = pd.DataFrame(columns=['Gene', 'Transcript_ID', 'Percent_of_Transcript_Supported_by_Hints',
                                      'Percent_Identity_exonerate', 'Percent_Similarity_exonerate', 'Gene_Cover_exonerate',
                                      'Sequence_Cover_exonerate', 'Normalized_Score_exonerate', 'Gene_Cover_mmseqs',
                                      'Percent_Identity_mmseqs', 'Chromosome', 'Gene_Start', 'Gene_End'])
    return outdf

def parseResults(gene_prediction_dir, exonerate_dir, mmseqs_results, gene_id_fname, summary_fname, nthreads):
    nthreads = int(nthreads)
    gene_prediction_gff_list = glob.glob(gene_prediction_dir + '/*.gff')
    tasks = []
    for gene_prediction_gff in gene_prediction_gff_list:
        genename = os.path.basename(gene_prediction_gff).replace('.gff', '')
        exonerate_results = os.path.join(exonerate_dir, genename + '.txt')
        tasks.append((genename, gene_prediction_gff, exonerate_results, summary_fname))

    with ProcessPoolExecutor(max_workers=nthreads) as executor:
        executor.map(parse_and_merge_onegene_wrapper, tasks)

    colnames = ['Gene', 'Transcript_ID', 'Gene_ID_exonerate', 'Percent_of_Transcript_Supported_by_Hints',
                'Percent_Identity_exonerate', 'Percent_Similarity_exonerate', 'Gene_Cover_exonerate',
                'Sequence_Cover_exonerate', 'Normalized_Score_exonerate', 'Chromosome', 'Gene_Start', 'Gene_End']
    merged_df = pd.read_csv(summary_fname.replace('.txt', '.tmp.txt'), sep='\t', header=None, names=colnames, dtype={'Chromosome': str})

    columns = ['geneid', 'chr', 'qstart', 'qend', 'start', 'end', 'alnlen', 'score', 'nident', 'qlen']
    mmseqs_df = pd.read_csv(mmseqs_results, sep='\t', header=None, names=columns, dtype={'chr': str})
    geneid_dic = pd.read_csv(gene_id_fname, sep='\t', header=None, names=['gene', 'geneid'])
    mmseqs_df = mmseqs_df.merge(geneid_dic, on='geneid', how='left')

    mmseqs_df['strand'] = mmseqs_df.apply(lambda row: '+' if row['start'] < row['end'] else '-', axis=1)
    mmseqs_df['Percent_Identity_mmseqs'] = mmseqs_df['nident'] / mmseqs_df['alnlen']
    mmseqs_df['alnset'] = mmseqs_df.apply(lambda row: set(range(row['qstart'], row['qend'] + 1)), axis=1)
    mmseqs_df.loc[mmseqs_df['strand'] == '-', ['start', 'end']] = mmseqs_df.loc[mmseqs_df['strand'] == '-', ['end', 'start']].values
    mmseqs_df = mmseqs_df.sort_values(by=['gene', 'geneid', 'chr', 'strand', 'start'])

    grouped_df = list(merged_df.groupby(['Gene']))
    tasks = [(gene, group, mmseqs_df.loc[(mmseqs_df['gene'] == gene)]) for gene, group in grouped_df]

    with ProcessPoolExecutor(max_workers=nthreads) as executor:
        results2 = executor.map(filter_and_calculate, tasks)

    outdf = pd.concat(results2)
    outdf.to_csv(summary_fname.replace('.txt', '_tmp.txt'), sep='\t', index=False, mode='w', header=True)

def keep_highest_score_non_overlapping(df_sorted):
    non_overlapping_regions = df_sorted.iloc[[0]]
    for i in range(1, len(df_sorted)):
        current_region = df_sorted.iloc[i]
        previous_region = non_overlapping_regions.iloc[-1]
        if current_region['Chromosome'] != previous_region['Chromosome'] or current_region['Gene_Start'] >= previous_region['Gene_End']:
            non_overlapping_regions = non_overlapping_regions.append(current_region)
    return non_overlapping_regions

def filterResults(results_path, model_path):
    model = keras.models.load_model(model_path)
    species = os.path.basename(model_path).replace('_NN_model.h5', '')

    threshold_df = pd.read_csv(model_path.replace(os.path.basename(model_path), 'Threshold95.txt'), header=None, sep=' ')
    threshold = threshold_df[threshold_df[0] == species].iloc[0, 1]

    df_tmp = pd.read_csv(results_path.replace('.txt', '_tmp.txt'), sep='\t').dropna(subset=['Chromosome', 'Gene_Start', 'Gene_End', 'Transcript_ID'])
    df_tmp['Percent_of_Transcript_Supported_by_Hints'] /= 100
    X = df_tmp[['Percent_of_Transcript_Supported_by_Hints', 'Percent_Identity_exonerate', 'Percent_Similarity_exonerate',
                'Gene_Cover_exonerate', 'Gene_Cover_mmseqs', 'Percent_Identity_mmseqs']]
    df_tmp['Pred_prob'] = model.predict(X)

    df_uniq = df_tmp.groupby(['Chromosome', 'Gene_Start', 'Gene_End'])
    df_uniq = df_uniq.apply(lambda x: x.sort_values('Pred_prob', ascending=False))
    df_uniq = df_uniq.reset_index(drop=True)
    df_uniq = df_uniq[df_uniq['Pred_prob'] >= df_uniq.groupby(['Chromosome', 'Gene_Start', 'Gene_End'])['Pred_prob'].transform('max') * 0.99]

    df_uniq = df_uniq.groupby(['Gene']).apply(lambda x: x.sort_values('Pred_prob', ascending=False)).reset_index(drop=True)
    df_uniq = df_uniq[df_uniq['Pred_prob'] >= df_uniq.groupby('Gene')['Pred_prob'].transform('max') * 0.95]
    df_uniq = df_uniq[df_uniq['Pred_prob'] >= threshold]
    df_uniq = df_uniq.sort_values(by=['Gene', 'Pred_prob'], ascending=[True, False])

    df_uniq = df_uniq.rename(columns={
        'Percent_Identity_exonerate': 'Percent_Identity_Protein_Alignment',
        'Percent_Similarity_exonerate': 'Percent_Similarity_Protein_Alignment',
        'Gene_Cover_exonerate': 'Gene_Coverage_Protein_Alignment',
        'Sequence_Cover_exonerate': 'Sequence_Coverage_Protein_Alignment',
        'Normalized_Score_exonerate': 'Normalized_Score_Protein_Alignement',
        'Gene_Cover_mmseqs': 'Gene_Coverage_Similarity_Search',
        'Percent_Identity_mmseqs': 'Percent_Identity_Similarity_Search',
        'Pred_prob': 'Prediction_Probability'
    })
    df_uniq.to_csv(results_path, sep='\t', index=False)
    df_uniq[['Chromosome', 'Gene_Start', 'Gene_End', 'Gene']].to_csv(results_path.replace('.txt', '.bed'), sep='\t', index=False, header=False)
    os.system('rm ' + results_path.replace('.txt', '.tmp.txt'))
    os.system('rm ' + results_path.replace('.txt', '_tmp.txt'))

def make_gff_output(results_file, gene_prediction_dir, out_gff):
    result = pd.read_csv(results_file, sep='\t', header=0)
    for i in range(result.shape[0]):
        gene = result.loc[i, 'Gene']
        transcript = result.loc[i, 'Transcript_ID'].replace('.t1', '')
        pred = result.loc[i, 'Prediction_Probability']
        count = result['Gene'].value_counts().get(gene, 0)
        if count == 1:
            j = ''
        else:
            j = 0
        with open(f"{gene_prediction_dir}/{gene}.gff") as gff_file:
            chr_start_end = 'stringnotfound'
            for line in gff_file:
                if line.startswith(f'# start gene {transcript}'):
                    if j != '':
                        j += 1
                    fields = next(gff_file).split('\t')
                    chromosome = fields[0].split(':')[0]
                    chr_start_end = fields[0]
                    gene_start = int(fields[0].split(':')[1].split('-')[0]) + int(fields[3])
                    gene_end = int(fields[0].split(':')[1].split('-')[0]) + int(fields[4])
                    with open(out_gff, 'a') as outfile:
                        if j == '':
                            outfile.write(f"{chromosome}\tPaGeSearch\tgene\t{gene_start}\t{gene_end}\t.\t{fields[6]}\t.\tgene_id={gene};gene={gene};probability={pred}\n")
                        else:
                            outfile.write(f"{chromosome}\tPaGeSearch\tgene\t{gene_start}\t{gene_end}\t.\t{fields[6]}\t.\tgene_id={gene}_{j};gene={gene};probability={pred}\n")
                if line.startswith(chr_start_end) and transcript in line:
                    fields = line.split('\t')
                    start = int(fields[0].split(':')[1].split('-')[0]) + int(fields[3])
                    end = int(fields[0].split(':')[1].split('-')[0]) + int(fields[4])
                    with open(out_gff, 'a') as outfile:
                        if j == '':
                            outfile.write(f"{chromosome}\tPaGeSearch\t{fields[2]}\t{start}\t{end}\t.\t{fields[6]}\t.\tgene_id={gene}\n")
                        else:
                            outfile.write(f"{chromosome}\tPaGeSearch\t{fields[2]}\t{start}\t{end}\t.\t{fields[6]}\t.\tgene_id={gene}_{j}\n")


def runPipeline(genome_tmp, pathway_gene_dir_tmp, outdir_tmp, outprefix, species_model, nthreads):
    cwd_full_path = os.getcwd()
    stime_tot = time.time()
    outdir = os.path.abspath(outdir_tmp)
    pathway_gene_dir = os.path.abspath(pathway_gene_dir_tmp)
    genome = os.path.abspath(genome_tmp)

    print('Start PaGeSearch')
    print('\nGenome assembly at: ' + genome)
    print('Query gene sequences at: ' + pathway_gene_dir)
    print('Archetype species: ' + species_model)
    print(f'Results will be saves at: {outdir}/{outprefix}.txt')
    print(f'Using {nthreads} threads')

    nextension = 20000

    if os.path.exists(outdir) == False:
        os.makedirs(outdir)
    os.chdir(outdir)
    if os.path.exists(outdir + '/tmp') == False:
        os.makedirs(outdir + '/tmp')
    os.environ['TMPDIR'] = outdir + '/tmp'

    if os.path.exists(outdir + '/MMseqs') == False:
        os.makedirs(outdir + '/MMseqs')
    if os.path.exists(outdir + '/MMseqs/' + outprefix) == False:
        os.makedirs(outdir + '/MMseqs/' + outprefix)
    else:
        os.system('rm -r ' + outdir + '/MMseqs/' + outprefix)
        os.makedirs(outdir + '/MMseqs/' + outprefix)

    mmseqs_dir = outdir + '/MMseqs/' + outprefix + '/'
    mmseqs_results = mmseqs_dir + outprefix + '.txt'

    # ✅ sanitize_ambiguous=True: BJOUZOU -> X 치환 적용
    mmseqs(genome, mmseqs_dir, pathway_gene_dir, mmseqs_results, nthreads, sanitize_ambiguous=True)

    if os.path.exists(outdir + '/ExtractSequences') == False:
        os.makedirs(outdir + '/ExtractSequences')
    if os.path.exists(outdir + '/ExtractSequences/' + outprefix) == False:
        os.makedirs(outdir + '/ExtractSequences/' + outprefix)
    if os.path.exists(outdir + '/ExtractSequences/' + outprefix + '/gene_fasta') == False:
        os.makedirs(outdir + '/ExtractSequences/' + outprefix + '/gene_fasta')
    else:
        os.system('rm -r ' + outdir + '/ExtractSequences/' + outprefix + '/gene_fasta')
        os.makedirs(outdir + '/ExtractSequences/' + outprefix + '/gene_fasta')

    if os.path.exists(outdir + '/ExtractSequences/' + outprefix + '/sequence_fasta') == False:
        os.makedirs(outdir + '/ExtractSequences/' + outprefix + '/sequence_fasta')
    else:
        os.system('rm -r ' + outdir + '/ExtractSequences/' + outprefix + '/sequence_fasta')
        os.makedirs(outdir + '/ExtractSequences/' + outprefix + '/sequence_fasta')

    if os.path.exists(outdir + '/ExtractSequences/' + outprefix + '/bed') == False:
        os.makedirs(outdir + '/ExtractSequences/' + outprefix + '/bed')
    else:
        os.system('rm -r ' + outdir + '/ExtractSequences/' + outprefix + '/bed')
        os.makedirs(outdir + '/ExtractSequences/' + outprefix + '/bed')

    seq_extract_dir = outdir + '/ExtractSequences/' + outprefix + '/'
    chr_len_fname = outdir + '/ExtractSequences/chrlen_dictionary.txt'
    print('Candidate region grouping and extension')

    seqlen_dic = {}
    with open(genome, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            seqlen_dic[record.id] = len(record.seq)

    with open(chr_len_fname, 'w') as f:
        for key, value in seqlen_dic.items():
            f.write('%s\t%s\n' % (key, value))

    gene_id_fname = seq_extract_dir + 'geneid_dic.txt'
    gene_fasta_list = glob.glob(pathway_gene_dir + '/*.fa')
    with open(gene_id_fname, 'a') as f:
        for gene_fasta in gene_fasta_list:
            genename = gene_fasta.split('/')[-1].replace('.fa', '')
            with open(gene_fasta, 'r') as g:
                for line in g:
                    if line.startswith('>'):
                        gene_id = line[1:].strip()
                        f.write(genename + '\t' + gene_id + '\n')

    seq_extract(genome, mmseqs_results, seq_extract_dir, nextension, chr_len_fname, pathway_gene_dir, gene_id_fname, nthreads)
    os.system('rm -r ' + seq_extract_dir + '/bed')
    for file_name in os.listdir(mmseqs_dir):
        if os.path.join(mmseqs_dir, file_name) != mmseqs_results:
            file_path = os.path.join(mmseqs_dir, file_name)
            os.remove(file_path)

    if os.path.exists(outdir + '/GenePrediction/') == False:
        os.makedirs(outdir + '/GenePrediction/')
    if os.path.exists(outdir + '/GenePrediction/' + outprefix) == False:
        os.makedirs(outdir + '/GenePrediction/' + outprefix)
    if os.path.exists(outdir + '/GenePrediction/' + outprefix + '/hints') == False:
        os.makedirs(outdir + '/GenePrediction/' + outprefix + '/hints')
    else:
        os.system('rm -r ' + outdir + '/GenePrediction/' + outprefix + '/hints')
        os.makedirs(outdir + '/GenePrediction/' + outprefix + '/hints')
    if os.path.exists(outdir + '/GenePrediction/' + outprefix + '/prediction') == False:
        os.makedirs(outdir + '/GenePrediction/' + outprefix + '/prediction')
    else:
        os.system('rm -r ' + outdir + '/GenePrediction/' + outprefix + '/prediction')
        os.makedirs(outdir + '/GenePrediction/' + outprefix + '/prediction')

    gene_fasta_dir = seq_extract_dir + '/gene_fasta/'
    fasta_extended_dir = seq_extract_dir + '/sequence_fasta/'
    exonerate_hints_dir = outdir + '/GenePrediction/' + outprefix + '/hints'
    gene_prediction_dir = outdir + '/GenePrediction/' + outprefix + '/prediction'
    fasta_extended_list = glob.glob(seq_extract_dir + '.fa')
    gene_prediction(gene_fasta_dir, fasta_extended_dir, exonerate_hints_dir, gene_prediction_dir, species_model, nthreads)
    os.system('rm -r ' + seq_extract_dir + '/*/')
    os.system('rm -r ' + exonerate_hints_dir)

    gene_prediction_dir = outdir + '/GenePrediction/' + outprefix + '/prediction'
    if os.path.exists(outdir + '/Exonerate/') == False:
        os.makedirs(outdir + '/Exonerate/')
    if os.path.exists(outdir + '/Exonerate/' + outprefix) == False:
        os.makedirs(outdir + '/Exonerate/' + outprefix)
    else:
        os.system('rm -r ' + outdir + '/Exonerate/' + outprefix)
        os.makedirs(outdir + '/Exonerate/' + outprefix)

    exonerate_dir = outdir + '/Exonerate/' + outprefix
    print('Protein sequence alignement')
    run_exonerate_multithreading(gene_prediction_dir, pathway_gene_dir, exonerate_dir, nthreads)

    summary_fname = outdir + '/' + outprefix + '.txt'
    print('Results filtering')
    parseResults(gene_prediction_dir, exonerate_dir, mmseqs_results, gene_id_fname, summary_fname, nthreads)

    species_dic = {'human': 'Homo_sapiens', 'zebrafish': 'Danio_rerio', 'chicken': 'Gallus_gallus',
                   'arabidopsis': 'Arabidopsis_thaliana', 'wheat': 'Triticum_aestivum', 'fly': 'Drosophila_melanogaster'}
    model_path = f"{cwd_full_path}/Models/{species_dic[species_model]}_NN_model.h5"
    filterResults(summary_fname, model_path)
    make_gff_output(summary_fname, gene_prediction_dir, summary_fname.replace('.txt', '.gff'))

    os.system('rm -r ' + outdir + '/GenePrediction')
    os.system('rm -r ' + outdir + '/Exonerate')
    os.system('rm -r ' + outdir + '/MMseqs')
    os.system('rm -r ' + outdir + '/ExtractSequences')
    os.system('rm -r ' + outdir + '/tmp')

    genes_found = pd.read_csv(summary_fname.replace('.txt', '.bed'), header=None, sep='\t')[3]
    genes_found = set(genes_found)
    query_genes = glob.glob(f"{pathway_gene_dir}/*.fa")
    query_genes = [i.split('/')[-1].replace('.fa', '') for i in query_genes]
    query_genes = set(query_genes)
    not_found = query_genes - genes_found
    with open(summary_fname.replace('.txt', '_notfound.txt'), 'w') as outfile:
        for gene in not_found:
            outfile.write(gene + '\n')

    print('Time consumed:  ' + str(time.time() - stime_tot) + '\n')


def file_exists(filepath):
    full_path = os.path.abspath(filepath)
    if not os.path.isfile(full_path):
        raise argparse.ArgumentTypeError(f"The file {full_path} does not exist!")
    return full_path

def directory_exists_and_not_empty(directory):
    full_directory = os.path.abspath(directory)
    if not os.path.isdir(directory):
        raise argparse.ArgumentTypeError(f"The directory {directory} does not exist!")
    if not os.listdir(directory):
        raise argparse.ArgumentTypeError(f"The directory {directory} is empty!")
    return full_directory

def validate_species(species):
    valid_species = ["arabidopsis", "wheat", "human", "chicken", "zebrafish"]
    if species.lower() not in valid_species:
        raise argparse.ArgumentTypeError(f"Invalid species: {species}. Must be one of {', '.join(valid_species)}.")
    return species

def main(argv):
    parser = argparse.ArgumentParser(
        prog='python Codes/pagesearch.py',
        usage='%(prog)s -g /path/to/genome.fa -p /directory/of/query/gene/sequences/folder/ -od ./-op pagesearch -s human -t 4',
        description='This script searches for pathway genes in a given genome.'
    )
    parser.add_argument('-g', '--genome', dest='genome', type=file_exists, required=True, help='Path to your genome sequence file. Required.')
    parser.add_argument('-p', '--pathway-geneseq-dir', dest='pathway_gene_dir', type=directory_exists_and_not_empty, required=True, help='Path to the folder containing your query gene sequences. Required.')
    parser.add_argument('-od', '--outdir', dest='outdir', type=str, action="store", default='./', help='Where the results will be saved (default "./").')
    parser.add_argument('-op', '--outprefix', dest='outprefix', type=str, action="store", default='pagesearch', help='Prefix for the output files (default "pagesearch").')
    parser.add_argument('-s', '--species', dest='modelspecies', type=validate_species, action="store", default='human', help='Archetype species (default "human").')
    parser.add_argument('-t', '--threads', dest='nthreads', type=int, action="store", default=4, help='Number of threads to use (default 4).')

    if len(argv) <= 1:
        parser.print_help()
        return

    args = parser.parse_args()

    try:
        stime = time.time()
        runPipeline(args.genome, args.pathway_gene_dir, args.outdir, args.outprefix, args.modelspecies, args.nthreads)
        print('\nPathway gene search in genome finished.')
        etime = time.time()
        print(f'Elapsed time: {etime - stime} seconds')
    except argparse.ArgumentTypeError as e:
        print(e)
        return

if __name__ == "__main__":
    main(sys.argv)
