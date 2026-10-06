import pandas as pd

def xlsx_sheet_to_fasta(xlsx_path, sheet_name, id_col, seq_col, output_fasta, header_cols=None):
    """
    Read a sheet from an Excel file and write it out as a FASTA file.

    Parameters
    ----------
    xlsx_path : str
        Path to the .xlsx file.
    sheet_name : str or int
        Name or index of the sheet to read.
    id_col : str
        Column name to use as the FASTA header ID (e.g. "Probe_ID").
    seq_col : str
        Column name containing the sequence (e.g. "Sequence").
    output_fasta : str
        Path to write the output FASTA file.
    header_cols : list of str, optional
        Additional column names to append to the FASTA header, 
        space-separated after the ID (e.g. ["Gene", "Target"]).
    """
    df = pd.read_excel(xlsx_path, sheet_name=sheet_name)

    # Drop rows with missing ID or sequence
    df = df.dropna(subset=[id_col, seq_col])

    with open(output_fasta, "w") as f:
        for _, row in df.iterrows():
            header = str(row[id_col]).strip()
            if header_cols:
                extra = " ".join(str(row[c]).strip() for c in header_cols if pd.notna(row[c]))
                if extra:
                    header += " " + extra
            seq = str(row[seq_col]).strip().replace(" ", "").replace("\n", "")
            f.write(f">{header}\n{seq}\n")

    print(f"Wrote {len(df)} sequences to {output_fasta}")


def filter_probes_and_save(nominated_file, gapfill_file, output_file, gapfill_sheet):
    nominated_df = pd.read_csv(nominated_file)
    gapfill_df = pd.read_excel(gapfill_file, sheet_name=gapfill_sheet)

    indices_to_remove = []

    for _, gap_row in gapfill_df.iterrows():
        gene = gap_row['gene']
        gap_start = gap_row['start']
        gap_end = gap_row['end']
        gene_length = gap_row['gene_length']

        exclusion_start = gene_length - gap_end
        exclusion_end = gene_length - gap_start

        for idx, probe_row in nominated_df.iterrows():
            probe_id = probe_row['id']
            probe_start = probe_row['start']
            probe_end = probe_row['end']

            if f"{gene}_bin" in probe_id:
                if exclusion_start <= probe_end and probe_start <= exclusion_end:
                    indices_to_remove.append(idx)
        
    filtered_df = nominated_df.drop(indices_to_remove)

    filtered_df.to_csv(output_file, index=False)

    print(f"Filtered probes have been saved to {output_file}")
    
    removed_df = nominated_df.iloc[indices_to_remove]
    print("Rows that will be removed:")
    print(removed_df)