import hashlib
import os
import requests
import subprocess


from math import log2
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from datetime import datetime
from pathlib import Path


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns



GC_TARGET = 58

# Sample barcodes, can be unique or shared with transcriptome barcodes 
BARCODE_IDS = [
    "BC001", "BC002", "BC003", "BC004", "BC005", "BC006", "BC007", "BC008",
    "BC009", "BC010", "BC011", "BC012", "BC013", "BC014", "BC015", "BC016",
]
BARCODE_SEQUENCES = [
    "ACTTTAGG", "AACGGGAA", "AGTAGGCT", "ATGTTGAC", "ACAGACCT", "ATCCCAAC",
    "AAGTAGAG", "AGCTGTGA", "ACAGTCTG", "AGTGAGTG", "AGAGGCAA", "ACTACTCA",
    "ATACGTCA", "ATCATGTG", "AACGCCGA", "ATTCGGTT",
]


def max_homopolymer_length(seq):
    """Length of the longest single-base run. Returns 0 for an empty sequence."""
    if not len(seq):
        return 0

    max_length = 1
    current_length = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i - 1]:
            current_length += 1
        else:
            current_length = 1

        max_length = max(max_length, current_length)

    return max_length


def calculate_gc_content(sequence):
    """Calculate the GC content of a DNA sequence."""
    gc_count = sequence.count("G") + sequence.count("C")
    total_bases = len(sequence)
    gc_content = (gc_count / total_bases) * 100
    return gc_content


def calculate_entropy(sequence):
    """Calculate the entropy of a DNA sequence."""
    probabilities = [float(sequence.count(base)) / len(sequence) for base in "ACGT"]
    entropy = -sum(p * log2(p) if p > 0 else 0 for p in probabilities)
    entropy = entropy / log2(4)
    return entropy


def score_sequence_repetition(sequence, min_n=2, max_n=5):
    """This function scores the degree of repetitive sequences in an input DNA seq. Default nmer range is 2-5."""
    nmer_counts = {n: {} for n in range(min_n, max_n + 1)}

    for n in range(min_n, max_n + 1):
        for i in range(len(sequence) - n + 1):
            nmer = sequence[i : i + n]
            if nmer in nmer_counts[n]:
                nmer_counts[n][nmer] += 1
            else:
                nmer_counts[n][nmer] = 1

    repetitiveness_score = sum(
        sum(len(nmer) * count for nmer, count in counts.items() if count > 1)
        for counts in nmer_counts.values()
    )

    return repetitiveness_score


def nominate_probes(
    fasta_path,
    gc_diff_threshold=14,
    homopolymer_length_threshold=4,
    entropy_threshold=0.8,
    repetitive_sequence_score_threshold=125,
):
    """
    This function identifies and nominates probes from a given FASTA file based on specified criteria.

    Parameters:
    -----------
    fasta_path : str
        The path to the input FASTA file containing DNA sequences for probe design.

    gc_diff_threshold : int
        The threshold for the difference in GC content from GC_TARGET (58%). The default value is 14%.

    homopolymer_length_threshold : int
        The threshold for the maximum allowed homopolymer length in probe sequences. The default is 4 bases.

    entropy_threshold : int
        Threshold for the minimum sequence entropy as a surrogate for sequence diversity. See calculate_entropy() for more information. Default threshold is set to 0.8.

    repetitive_sequence_score_threshold : int
         Threshold for the maximum sequence repetitiveness calculated using score_sequence_repetition(). Numerates the number of 2-5 nmer repeats in LHS and RHS probes and aggregates to an int score.
         See score_sequence_repetition() for details. Default threshold is set to 125.

    Returns:
    --------
    Creates a directory probe_nomination in the current working directory. Creates an annotation CSV with the measured probe 5' and 3' sites as well as the maximum homopolymer length and GC content.
    Also writes the probes that pass the threshold in FASTA format. Returns the annotation table as a DataFrame.

    Prints output file paths and the number of probe sites that passed the threshold.

    Note on ordering: the thresholds are tested cheapest-first, and each 25mer's
    metrics are cached, because the RHS arm of one bin is the LHS arm of the bin
    25 bases downstream. Output is identical to testing every metric up front.
    """
    annotation_rows = []
    pass_sequences_records = []
    sequence_bins_count = {}

    base_filename = os.path.basename(fasta_path).rsplit(".", 1)[0]

    output_directory = os.path.join(os.getcwd(), "probe_nomination")
    os.makedirs(output_directory, exist_ok=True)

    with open(fasta_path) as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            reverse_complement_sequence = record.seq.reverse_complement().upper()
            sequence = str(reverse_complement_sequence)
            bin_count = 0

            arm_cache = {}

            def arm_metrics(pos):
                if pos not in arm_cache:
                    arm = sequence[pos : pos + 25]
                    arm_cache[pos] = (
                        round(calculate_gc_content(arm), 3),
                        round(calculate_entropy(arm), 3),
                        max_homopolymer_length(arm),
                        score_sequence_repetition(arm),
                    )
                return arm_cache[pos]

            for i in range(0, len(sequence) - 50 + 1):
                if sequence[i + 24] != "T":
                    continue

                gc_content_lhs, lhs_entropy, lhs_max_homopolymer, lhs_repetitiveness = arm_metrics(i)
                gc_content_rhs, rhs_entropy, rhs_max_homopolymer, rhs_repetitiveness = arm_metrics(i + 25)

                if not (
                    abs(gc_content_lhs - GC_TARGET) <= gc_diff_threshold
                    and abs(gc_content_rhs - GC_TARGET) <= gc_diff_threshold
                    and lhs_max_homopolymer <= homopolymer_length_threshold
                    and rhs_max_homopolymer <= homopolymer_length_threshold
                    and lhs_entropy >= entropy_threshold
                    and rhs_entropy >= entropy_threshold
                    and lhs_repetitiveness <= repetitive_sequence_score_threshold
                    and rhs_repetitiveness <= repetitive_sequence_score_threshold
                ):
                    continue

                start_position = i + 1
                probe_id = f"{record.id}_bin_{start_position}"
                description = f"50-bp bin starting at position {start_position}"

                annotation_rows.append(
                    {
                        "id": probe_id,
                        "target": record.id,
                        "description": description,
                        "start": start_position,
                        "end": start_position + 49,
                        "lhs_gc_content": gc_content_lhs,
                        "rhs_gc_content": gc_content_rhs,
                        "lhs_max_homopolymer": lhs_max_homopolymer,
                        "rhs_max_homopolymer": rhs_max_homopolymer,
                        "lhs_entropy": lhs_entropy,
                        "rhs_entropy": rhs_entropy,
                        "lhs_repetitiveness": lhs_repetitiveness,
                        "rhs_repetitiveness": rhs_repetitiveness,
                        "lhs_sequence": sequence[i : i + 25],
                        "rhs_sequence": sequence[i + 25 : i + 50],
                    }
                )
                pass_sequences_records.append(
                    SeqRecord(
                        reverse_complement_sequence[i : i + 50],
                        id=probe_id,
                        description=description,
                    )
                )
                bin_count += 1

            sequence_bins_count[record.id] = bin_count

    annotations = pd.DataFrame(annotation_rows)

    annotations_csv = os.path.join(output_directory, f"{base_filename}_annotations.csv")
    annotations.to_csv(annotations_csv, index=False)
    print(f"Annotations written to: {annotations_csv}")

    output_fasta = os.path.join(output_directory, f"{base_filename}_probes.fasta")
    SeqIO.write(pass_sequences_records, output_fasta, "fasta")
    print(f"Probes written to: {output_fasta}")

    for record_id, bin_count in sequence_bins_count.items():
        print(
            f"Sequence {record_id} has {bin_count} potential probe binding sites passing threshold."
        )

    return annotations


def lhs_rhs_probe_binding_site_split(input_fasta, output_fasta):
    """This function ouputs the LHS and RHS probes as a FASTA from input probe-binding site nomination FASTA"""
    with open(input_fasta, "r") as input_file, open(output_fasta, "w") as output_file:
        for record in SeqIO.parse(input_file, "fasta"):
            lhs_sequence = record.seq[:25]
            rhs_sequence = record.seq[-25:]

            new_record = (
                f">{record.id}_LHS\n{lhs_sequence}\n>{record.id}_RHS\n{rhs_sequence}\n"
            )
            output_file.write(new_record)
    print(f"LHS and RHS probes written to: {output_fasta}")


def read_blast_results(file_path):
    """
    Reads a TSV file into a DataFrame, filters data, calculates columns,
    performs pattern matching, and extracts information.

    Parameters:
    - file_path (str): The path to the TSV file.

    Returns:
    - lhs_df, rhs_df (tuple of DataFrames): Two DataFrames containing processed data.
    """
    column_names = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "evalue",
        "bitscore",
        "stitle",
    ]
    blast_out = pd.read_csv(file_path, sep="\t", header=None, names=column_names)
    blast_out["off_target_binding"] = blast_out["length"] - blast_out["mismatch"]
    probe_max_offtarg = blast_out.groupby("qseqid")["off_target_binding"].idxmax()
    blast_out = blast_out.loc[probe_max_offtarg].copy()

    pattern = r"\((.*?)\),"
    blast_out["off_target_gene"] = blast_out["stitle"].str.extract(pattern, expand=False)

    lhs_df = blast_out[blast_out["qseqid"].str.endswith("_LHS")].copy()
    rhs_df = blast_out[blast_out["qseqid"].str.endswith("_RHS")].copy()
    return lhs_df, rhs_df


def process_lhs_rhs(lhs_df, rhs_df, annotation_csv_path, output_csv_path):
    """
    Merge information from 'lhs_df' and 'rhs_df' DataFrames with the annotation table
    from nominate_probes() and save the merged DataFrame to a CSV file.

    Parameters:
    -----------
    lhs_df : pandas DataFrame
        Left-hand side BLAST hits, with a 'qseqid' column of the form '<probe id>_LHS'.
    rhs_df : pandas DataFrame
        Right-hand side BLAST hits, with a 'qseqid' column of the form '<probe id>_RHS'.
    annotation_csv_path : str
        Filepath to the annotation CSV from nominate_probes().
    output_csv_path : str
        Filepath to save the resulting merged DataFrame.

    Returns:
    --------
    The merged DataFrame (also written to output_csv_path). Probes with no BLAST
    hit keep NaN off-target columns, which calculate_priority_score() reads as
    zero off-target binding.
    """
    target_table = pd.read_csv(annotation_csv_path)

    merged = target_table
    for side, side_df in (("lhs", lhs_df), ("rhs", rhs_df)):
        side_df = side_df.copy()
        side_df["id"] = side_df["qseqid"].str.replace(f"_{side.upper()}$", "", regex=True)
        merged = merged.merge(
            side_df[["id", "length", "mismatch", "off_target_binding", "off_target_gene"]],
            on="id",
            how="left",
        ).rename(
            columns={
                "length": f"{side}_off_target_length",
                "mismatch": f"{side}_mismatch",
                "off_target_binding": f"{side}_max_off_target_binding",
                "off_target_gene": f"{side}_off_target_gene",
            }
        )

    merged.to_csv(output_csv_path, index=False)
    return merged


def calculate_priority_score(anno, gc_target=GC_TARGET):
    """
    Score each probe pair on GC content, off-target binding, and sequence entropy.

    The score is the product of three bounded sub-scores, each in (0, 1], scaled
    to 0-100:

        gc_score        = 1 / (1 + (mean |GC - gc_target| / 10) ** 2)
        off_target_score= 1 / (1 + (lhs_hit / 10) ** 2 + (rhs_hit / 10) ** 2)
        entropy_score   = min(lhs_entropy, rhs_entropy)

        priority_score  = 100 * gc_score * off_target_score * entropy_score

    A probe with perfect GC, no off-target hits and maximal entropy scores 100;
    everything else scores less. Scores are directly comparable across targets.

    This replaces the earlier formula, which divided by the GC deviation. That
    made the score unbounded (a probe sitting exactly at the GC optimum divided
    by zero and returned inf, taking its target's whole set with it) and made
    small GC differences near the optimum dominate everything else. The ranking
    logic is otherwise unchanged: same three inputs, same directions, same
    quadratic off-target penalty.

    Parameters:
    - anno (DataFrame): Annotation DataFrame which is the output of process_lhs_rhs().
    - gc_target (float): GC content (%) treated as optimal, matching nominate_probes().

    Returns:
    - anno (DataFrame): DataFrame with added 'priority_score' and 'target' columns.
    """
    anno = anno.copy()

    lhs_off_target = anno["lhs_max_off_target_binding"].fillna(0)
    rhs_off_target = anno["rhs_max_off_target_binding"].fillna(0)

    gc_deviation = (
        (anno["lhs_gc_content"] - gc_target).abs()
        + (anno["rhs_gc_content"] - gc_target).abs()
    ) / 2
    gc_score = 1 / (1 + (gc_deviation / 10) ** 2)

    off_target_score = 1 / (
        1 + (lhs_off_target / 10) ** 2 + (rhs_off_target / 10) ** 2
    )

    entropy_score = anno[["lhs_entropy", "rhs_entropy"]].min(axis=1)

    anno["priority_score"] = 100 * gc_score * off_target_score * entropy_score

    if "target" not in anno.columns:
        anno["target"] = anno["id"].str.rsplit("_bin_", n=1).str[0]

    return anno


def build_compatibility_matrix(group, overlap_buffer=0):
    """
    Boolean (n x n) matrix answering "can probes i and j be in the same set?".

    Two probes are compatible if their 50-bp intervals do not overlap (padded by
    overlap_buffer bp) and they share no off-target gene on either arm.
    Indexing is positional within `group`, not by DataFrame index; the diagonal
    is False.

    Parameters:
    - group (DataFrame): rows for a single target.
    - overlap_buffer (int): extra bp required between selected probes. 0 enforces
      physical non-overlap only.

    Returns:
    - compat (np.ndarray[bool], n x n)
    """
    n = len(group)
    starts = group["start"].to_numpy()
    ends = group["end"].to_numpy()

    overlaps = (starts[:, None] <= ends[None, :] + overlap_buffer) & (
        starts[None, :] <= ends[:, None] + overlap_buffer
    )

    lhs_genes = group["lhs_off_target_gene"].to_numpy()
    rhs_genes = group["rhs_off_target_gene"].to_numpy()
    all_genes = pd.unique(
        pd.Series(np.concatenate([lhs_genes, rhs_genes])).dropna()
    )

    if len(all_genes) == 0:
        shares_gene = np.zeros((n, n), dtype=bool)
    else:
        gene_to_col = {gene: col for col, gene in enumerate(all_genes)}
        gene_matrix = np.zeros((n, len(all_genes)), dtype=np.int16)
        for row, (lhs_gene, rhs_gene) in enumerate(zip(lhs_genes, rhs_genes)):
            if pd.notna(lhs_gene):
                gene_matrix[row, gene_to_col[lhs_gene]] = 1
            if pd.notna(rhs_gene):
                gene_matrix[row, gene_to_col[rhs_gene]] = 1
        shares_gene = (gene_matrix @ gene_matrix.T) > 0

    compat = ~overlaps & ~shares_gene
    np.fill_diagonal(compat, False)
    return compat


def find_best_probe_set(group, overlap_buffer=0):
    """
    Find the highest-scoring set of 3 compatible probes for one target, falling
    back to the best pair and then the best single probe.

    This is an exhaustive search, not a greedy one: every compatible pair is
    considered, and for each the best compatible third probe is found by
    vectorized lookup. Candidates are visited in descending score order so the
    search can stop once no remaining candidate could beat the best set found,
    which keeps it fast (milliseconds) even with ~1000 candidates per target.

    Parameters:
    - group (DataFrame): candidate rows for a single target, with 'priority_score',
      'start', 'end', 'lhs_off_target_gene' and 'rhs_off_target_gene'.
    - overlap_buffer (int): passed through to build_compatibility_matrix().

    Returns:
    - list of DataFrame index labels for the selected probes (never None unless
      the group is empty).
    """
    n = len(group)
    if n == 0:
        return None

    index_labels = group.index.to_numpy()
    if n == 1:
        return [index_labels[0]]

    compat = build_compatibility_matrix(group, overlap_buffer=overlap_buffer)

    order = np.argsort(-group["priority_score"].to_numpy())
    scores = group["priority_score"].to_numpy()[order]
    compat = compat[np.ix_(order, order)]

    best_score = -np.inf
    best_set = None

    for i in range(n):
        if scores[i] * 3 <= best_score:
            break 
        for j in range(i + 1, n):
            if not compat[i, j]:
                continue
            if scores[i] + scores[j] * 2 <= best_score:
                break 
            mask = compat[i] & compat[j]
            mask[: j + 1] = False  
            if not mask.any():
                continue
            k = int(np.argmax(np.where(mask, scores, -np.inf)))
            total = scores[i] + scores[j] + scores[k]
            if total > best_score:
                best_score = total
                best_set = (i, j, k)

    if best_set is None:
        # No compatible triple. Try the best compatible pair.
        best_score = -np.inf
        for i in range(n):
            if scores[i] * 2 <= best_score:
                break
            js = np.nonzero(compat[i, i + 1 :])[0] + (i + 1)
            if len(js) == 0:
                continue
            j = int(js[np.argmax(scores[js])])
            if scores[i] + scores[j] > best_score:
                best_score = scores[i] + scores[j]
                best_set = (i, j)

    if best_set is None:
        best_set = (0,)  # scores are sorted, so position 0 is the best probe

    return [index_labels[order[p]] for p in best_set]


def nominate_top_probe_set(input_csv, output_csv, overlap_buffer=0):
    """
    Process an annotation CSV file, calculate priority scores, select the
    highest-scoring compatible set of up to 3 probes per target, and write them
    to a new CSV file.

    Expects the input CSV to have the following columns: 'id', 'start', 'end',
    'lhs_off_target_gene', 'rhs_off_target_gene', 'lhs_gc_content',
    'rhs_gc_content', 'lhs_entropy', 'rhs_entropy', 'lhs_max_off_target_binding'
    and 'rhs_max_off_target_binding'.

    Parameters:
    - input_csv (str): Path to the input annotation CSV file.
    - output_csv (str): Path to the output CSV file to store the selected probe sets.
    - overlap_buffer (int): bp required between selected probes, on top of
      non-overlap. 0 enforces non-overlap only.

    Returns:
    - The selected probes as a DataFrame, also written to output_csv.
    """
    annotated_data = calculate_priority_score(pd.read_csv(input_csv))

    selected_indices = []

    for target, group in annotated_data.groupby("target"):
        best_set = find_best_probe_set(group, overlap_buffer=overlap_buffer)

        if not best_set:
            print(f"ERROR: No probes nominated for target {target}")
            continue
        if len(best_set) < 3:
            print(
                f"WARNING: Only {len(best_set)} probe(s) nominated for target {target}"
            )
        selected_indices.extend(best_set)

    selected = annotated_data.loc[selected_indices].sort_values(["target", "start"])
    selected.to_csv(output_csv, index=False)
    return selected


def construct_custom_probe_opool(nominated_probe_csv, pool_name, lhs_r2_adaptor, n_barcodes, directory = 'opools'):
    """
    Constructs an Excel file containing probe sequences and pool names from a nominated probe CSV file. This xlsx should be able to be uploaded to IDT for ordering.

    Parameters:
    nominated_probe_csv (str): Path to the nominated probe CSV file from nominate_top_probe_set().
    pool_name (str): Name used in the pool names and the output filename.
    lhs_r2_adaptor (str): Partial TruSeq R2 adaptor prepended to each LHS probe.
    n_barcodes (int): Number of barcode sequences to consider. Typically 4 or 16.
    directory (str): Output directory for the .xlsx order sheet.

    Returns:
    None
    """
    data = pd.read_csv(nominated_probe_csv)
    sequences = []

    if n_barcodes is None:
        n_barcodes = len(BARCODE_SEQUENCES)
    else:
        n_barcodes = min(n_barcodes, len(BARCODE_SEQUENCES))

    for index, row in data.iterrows():
        lhs_sequence = row['lhs_sequence']
        rhs_sequence = row['rhs_sequence']

        lhs_seq =  lhs_r2_adaptor + lhs_sequence
        lhs_pool = pool_name + "_lhs"
        sequences.append({'Pool name': lhs_pool, 'Sequence': lhs_seq})

        for i in range(n_barcodes):
            seq = "/5Phos/" + rhs_sequence + "ACGCGGTTAGCACGTANN" + BARCODE_SEQUENCES[i] + "CGGTCCTAGCAA"
            pool = BARCODE_IDS[i] + "_" + pool_name + "_rhs"
            sequences.append({'Pool name': pool, 'Sequence': seq})

    oligo_pool = pd.DataFrame(sequences)
    oligo_pool = oligo_pool.sort_values(by='Pool name')

    date_today = datetime.now().strftime("%y%m%d")
    file_name = f"{date_today}_{pool_name}.xlsx"

    if not os.path.exists(directory):
        os.makedirs(directory)

    file_path = os.path.join(directory, file_name)

    with pd.ExcelWriter(file_path) as writer:
        oligo_pool.to_excel(writer, index=False)



def gene_to_ensembl_wta_filter(gene_names, input_csv = 'references/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv'):
    """
    Fetches Ensembl IDs for a list of gene names and filters data from a CSV file based on matching Ensembl IDs.

    Args:
    - gene_names (list): List of gene names to fetch Ensembl IDs for.
    - input_csv (str, optional): Path to the CSV file containing gene data.
      Defaults to 'references/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv'.

    Returns:
    - filtered_data (pandas DataFrame): Filtered data from the CSV based on matched Ensembl IDs.

    Note:
    - Uses the Ensembl REST batch POST endpoint, which avoids the per-request
      rate limiting that makes one-symbol-at-a-time lookups fail partway through
      a large panel.
    - Prints any gene names that could not be resolved.
    """
    ensembl_ids = {}
    endpoint = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
    gene_names = list(dict.fromkeys(gene_names))

    with requests.Session() as session:
        for offset in range(0, len(gene_names), 200):
            batch = gene_names[offset : offset + 200]
            response = session.post(
                endpoint,
                json={"symbols": batch},
                headers={"Content-Type": "application/json", "Accept": "application/json"},
                timeout=30,
            )
            if not response.ok:
                print(f"Failed to retrieve data for {len(batch)} gene name(s). Check your input or try again later.")
                continue
            for symbol, record in response.json().items():
                if isinstance(record, dict) and "id" in record:
                    ensembl_ids[symbol] = record["id"]

    missing_genes = [gene for gene in gene_names if gene not in ensembl_ids]
    if missing_genes:
        print(f"Ensembl ID not found for: {', '.join(missing_genes)}")

    data = pd.read_csv(input_csv, names=['gene_id', 'probe_seq', 'probe_id', 'included', 'region'], skiprows=6)

    filtered_data = data[data['gene_id'].isin(ensembl_ids.values())]
    return filtered_data


def generate_flex_sgrna_opool(a_position_protospacers, b_position_protospacers, pool_name, n_barcodes, directory = 'opools'):
    """
    From a list of A and B protospacers generates flex probes.

    Args:
    - a_position_protospacers (Series): Protospacer sequences for position A.
    - b_position_protospacers (Series): Protospacer sequences for position B.
    - pool_name (str): Name of the pool.
    - n_barcodes (int): Number of barcode sequences.
    - directory (str): Output directory for the .xlsx order sheet.

    Output:
    - Creates an Excel file named YYMMDD_pool_name.xlsx containing probe orders in directory opools.
    """
    sgrna_a_rhs_constant_region = 'CTGAAAC'
    sgrna_b_rhs_constant_region = 'AAAC'
    partial_capture_seqeuence_1 = 'CGGTCCTAGCAA'

    parital_truseq_r2 = "CAGACGTGTGCTCTTCCGATCT"
    lhs_constant_sequence = "ACGCGGTTAGCACGTANN"
    sgrna_a_lhs_constant_region = "CTTGCTATGCACTCTTGTGCTTAGCT"
    sgrna_b_lhs_constant_region = "GCTATGCTGTTTCCAGCTTAGCTCTT"

    a_seq = a_position_protospacers.apply(Seq)
    b_seq = b_position_protospacers.apply(Seq)

    pool_names_rhs = [f"RHS_{pool_name}"] * (len(a_seq) * 2)
    sequences_rhs = []

    for a, b in zip(a_seq, b_seq):
        sequences_rhs.extend([
            f"/5Phos/{sgrna_a_rhs_constant_region}{a.reverse_complement()}{partial_capture_seqeuence_1}",
            f"/5Phos/{sgrna_b_rhs_constant_region}{b.reverse_complement()}{partial_capture_seqeuence_1}"
        ])

    rhs_probe_order = pd.DataFrame({
        'Pool name': pool_names_rhs,
        'Sequence': sequences_rhs
    })

    pool_names_lhs = []
    sequences_lhs = []

    for i in range(n_barcodes):
        pool_names_lhs.extend([
            f"{BARCODE_IDS[i]}_LHS_{pool_name}",
            f"{BARCODE_IDS[i]}_LHS_{pool_name}"
        ])

        sequences_lhs.extend([
            f"{parital_truseq_r2}{lhs_constant_sequence}{BARCODE_SEQUENCES[i]}{sgrna_a_lhs_constant_region}",
            f"{parital_truseq_r2}{lhs_constant_sequence}{BARCODE_SEQUENCES[i]}{sgrna_b_lhs_constant_region}"
        ])

    lhs_probe_order = pd.DataFrame({
        'Pool name': pool_names_lhs,
        'Sequence': sequences_lhs
    })

    oligo_pool = pd.concat([rhs_probe_order, lhs_probe_order])

    date_today = datetime.now().strftime("%y%m%d")
    file_name = f"{date_today}_{pool_name}.xlsx"

    if not os.path.exists(directory):
        os.makedirs(directory)

    file_path = os.path.join(directory, file_name)

    with pd.ExcelWriter(file_path) as writer:
        oligo_pool.to_excel(writer, index=False)


def generate_flex_gene_opool(gene_list, pool_name, n_barcodes, directory = 'opools'):
    """
    From stock Flex reference probes generates a barcoded oPool order.

    Args:
    - gene_list (DataFrame): Rows of the Flex probe reference, with a 'probe_seq'
      column holding the full-length probe; arms are taken from each end.
    - pool_name (str): Name of the pool.
    - n_barcodes (int): Number of barcode sequences.
    - directory (str): Output directory for the .xlsx order sheet.

    Output:
    - Creates an Excel file named YYMMDD_pool_name.xlsx containing probe orders in directory opools.
    """

    partial_capture_sequence_1 = 'CGGTCCTAGCAA'
    rhs_consatant_sequence = 'ACGCGGTTAGCACGTANN'
    parital_truseq_r2 = "CAGACGTGTGCTCTTCCGATCT"

    sequences = []

    for index, row in gene_list.iterrows():
        lhs_sequence = row['probe_seq'][:25]
        rhs_sequence = row['probe_seq'][-25:]

        lhs_seq =  parital_truseq_r2 + lhs_sequence
        lhs_pool = pool_name + "_lhs"
        sequences.append({'Pool name': lhs_pool, 'Sequence': lhs_seq})

        for i in range(n_barcodes):
            seq = "/5Phos/" + rhs_sequence + rhs_consatant_sequence + BARCODE_SEQUENCES[i] + partial_capture_sequence_1
            pool = BARCODE_IDS[i] + "_" + pool_name + "_rhs"
            sequences.append({'Pool name': pool, 'Sequence': seq})

    sequences_df = pd.DataFrame(sequences)
    sequences_df = sequences_df.sort_values(by='Pool name')

    date_today = datetime.now().strftime("%y%m%d")
    file_name = f"{date_today}_{pool_name}.xlsx"

    if not os.path.exists(directory):
        os.makedirs(directory)

    file_path = os.path.join(directory, file_name)

    with pd.ExcelWriter(file_path) as writer:
        sequences_df.to_excel(writer, index=False)


def probe_seq_hash(input_sequence):
    """
    Computes a shortened SHA-256 hash of the given input sequence.

    Parameters:
    - input_sequence (str): The input DNA sequence to be hashed.

    Returns:
    - str: The first 7 characters of the SHA-256 hash of the input sequence.
    """
    sha256_hash = hashlib.sha256(input_sequence.encode()).hexdigest()
    short_hash = sha256_hash[:7]

    return short_hash


def construct_probe_reference(
    custom_probes,
    custom_reference,
    genome_reference="references/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv",
):
    """
    Constructs a custom probe reference by combining user-defined probes with a genome reference.

    Parameters:
    - custom_probes (str): Path to a CSV file containing custom probe information, including columns 'target', 'lhs_sequence', and 'rhs_sequence'.
    - custom_reference (str): Path to the output updated custom reference.
    - genome_reference (str): Path to the genome reference CSV file (default is provided GRCh38 v1.0.1 reference).

    Returns:
    - None: The function writes the updated reference to the specified custom_reference file.
    """
    custom_probes = pd.read_csv(custom_probes)

    custom_probe_ref = pd.DataFrame(
        {
            "gene_id": custom_probes["target"],
            "probe_seq": custom_probes["lhs_sequence"] + custom_probes["rhs_sequence"],
            "probe_id": custom_probes["target"]
            + "|"
            + custom_probes["target"]
            + "|"
            + (custom_probes["lhs_sequence"] + custom_probes["rhs_sequence"]).apply(
                lambda x: probe_seq_hash(x)
            ),
            "included": True,
            "region": "unspliced",
        }
    )

    with open(genome_reference, "r") as f:
        lines = f.readlines()

    comments = [line.strip() for line in lines if line.startswith("#")]
    genome_reference = pd.read_csv(genome_reference, comment="#", header=0)
    updated_reference = pd.concat(
        [genome_reference, custom_probe_ref], ignore_index=True
    )

    with open(custom_reference, "w") as f:
        f.write("\n".join(comments) + "\n")
        updated_reference.to_csv(f, index=False)