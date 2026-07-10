import os
import re
import gzip
import pandas as pd

RAW_DIR = "GSE1751"
OUTPUT_DIR = "datasets"

MARKERS = [
    "TUT7", "PPP1CC", "CEBPB", "CEBPD", "FTO", "METTL16",
    "IGF2BP3", "YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2",
]

MATRIX_FILE = os.path.join(RAW_DIR, "GSE1751_series_matrix.txt")
ANNOT_FILE = os.path.join(RAW_DIR, "GPL96.annot")

GENE_AVERAGED_OUT = os.path.join(OUTPUT_DIR, "huntingtons_ML_dataset_GSE1751.csv")

ID_COL = "Sample_ID"
TARGET_COL = "Huntingtons_Disease"

# Decision point: how to handle presymptomatic HD gene-carriers (P1-P5).
# "exclude" drops them entirely (clean symptomatic-vs-control comparison).
# "hd" labels them as 1 (treats "carries the HD mutation" as the target).
PRESYMPTOMATIC_HANDLING = "exclude"  # or "hd"


def load_annotation(annot_file, marker_set):
    probe_to_gene = {}
    open_func = gzip.open if annot_file.endswith(".gz") else open
    mode = "rt" if annot_file.endswith(".gz") else "r"

    try:
        with open_func(annot_file, mode, encoding="utf-8") as f:
            in_table = False
            headers = []
            id_idx, gene_idx = -1, -1

            for line in f:
                if line.startswith("!platform_table_begin"):
                    in_table = True
                    continue

                if in_table and not headers:
                    headers_line = line.strip("\n").split("\t")
                    headers = [h.strip('"') for h in headers_line]
                    if "ID" in headers and "Gene symbol" in headers:
                        id_idx = headers.index("ID")
                        gene_idx = headers.index("Gene symbol")
                    continue

                if in_table and headers:
                    if line.startswith("!platform_table_end"):
                        break
                    parts = line.strip("\n").split("\t")
                    if len(parts) > max(id_idx, gene_idx):
                        probe_id = parts[id_idx].strip()
                        gene_str = parts[gene_idx].strip()
                        if gene_str:
                            mapped_genes = [g.strip() for g in gene_str.split("///")]
                            for gene in mapped_genes:
                                if gene in marker_set:
                                    probe_to_gene[probe_id] = gene
                                    break
    except Exception as e:
        print(f"Error reading {annot_file}: {e}")

    return probe_to_gene


def label_from_description(desc):
    """Parse !Sample_description text into a label, handling the presymptomatic case explicitly."""
    d = desc.lower()
    if "presymptomatic" in d:
        if PRESYMPTOMATIC_HANDLING == "exclude":
            return None
        elif PRESYMPTOMATIC_HANDLING == "hd":
            return 1
    if "huntington" in d or "hd" in d:
        return 1
    if "control" in d or "healthy" in d or "normal" in d:
        return 0
    return None


def process_geo_matrix(matrix_file, probe_to_gene):
    sample_ids = []
    sample_data = {}

    with open(matrix_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("!Sample_geo_accession"):
                sample_ids = re.findall(r"GSM\d+", line)
                for gsm in sample_ids:
                    sample_data[gsm] = {"Sample_ID": gsm, "Huntingtons_Disease": None}

            elif line.startswith("!Sample_description"):
                parts = line.split("\t")
                for i, val in enumerate(parts[1:]):
                    if i < len(sample_ids):
                        gsm = sample_ids[i]
                        sample_data[gsm]["Huntingtons_Disease"] = label_from_description(val)

            elif line.startswith('"ID_REF"') or line.startswith("!"):
                continue

            else:
                parts = line.rstrip("\n").split("\t")
                probe_id = parts[0].strip('"')
                if probe_id in probe_to_gene:
                    gene = probe_to_gene[probe_id]
                    feature_name = f"{gene}_{probe_id}"
                    for i, val in enumerate(parts[1:]):
                        if i < len(sample_ids):
                            gsm = sample_ids[i]
                            try:
                                sample_data[gsm][feature_name] = float(val.strip('"'))
                            except ValueError:
                                pass

    return sample_data


def get_gene_symbol(col_name):
    return col_name.split("_")[0]


def main():
    probe_to_gene = load_annotation(ANNOT_FILE, MARKERS)
    found_genes = set(probe_to_gene.values())
    missing_genes = [g for g in MARKERS if g not in found_genes]
    print(f"Mapped genes found on GPL96: {len(found_genes)}/{len(MARKERS)}")
    if missing_genes:
        print(f"WARNING - genes not found on this platform: {missing_genes}")

    data = process_geo_matrix(MATRIX_FILE, probe_to_gene)

    combined = [row for row in data.values() if row.get("Huntingtons_Disease") is not None]
    excluded_count = len(data) - len(combined)
    print(f"\nSamples with a valid label: {len(combined)} / {len(data)} total ({excluded_count} excluded)")

    df = pd.DataFrame(combined)
    cols = [ID_COL, TARGET_COL] + [c for c in df.columns if c not in (ID_COL, TARGET_COL)]
    df = df[cols]
    print(f"Class balance:\n{df[TARGET_COL].value_counts()}")

    feature_cols = [c for c in df.columns if c not in (ID_COL, TARGET_COL)]
    probes_per_gene = {}
    for col in feature_cols:
        gene = get_gene_symbol(col)
        probes_per_gene.setdefault(gene, []).append(col)

    print("\nProbes per gene:")
    for gene, probes in probes_per_gene.items():
        print(f"  {gene}: {len(probes)} probe(s) -> {probes}")

    gene_avg_df = pd.DataFrame(index=df.index)
    for gene, probes in probes_per_gene.items():
        gene_avg_df[gene] = df[probes].mean(axis=1)

    final_df = pd.concat([df[[ID_COL, TARGET_COL]], gene_avg_df], axis=1)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    final_df.to_csv(GENE_AVERAGED_OUT, index=False)

    print(f"\nFinal shape: {final_df.shape}")
    print(f"Saved: {GENE_AVERAGED_OUT}")


if __name__ == "__main__":
    main()