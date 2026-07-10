
import os
import re
import gzip
import pandas as pd

RAW_DIR = "GSE3790"
OUTPUT_DIR = "datasets"

MARKERS = [
    "TUT7", "PPP1CC", "CEBPB", "CEBPD", "FTO", "METTL16",
    "IGF2BP3", "YTHDF1", "YTHDF2", "YTHDF3", "YTHDC1", "YTHDC2",
]

PLATFORMS = [
    {
        "matrix_file": os.path.join(RAW_DIR, "GSE3790-GPL96_series_matrix.txt"),
        "annot_file": os.path.join(RAW_DIR, "GPL96.annot"),
    },
    {
        "matrix_file": os.path.join(RAW_DIR, "GSE3790-GPL97_series_matrix.txt"),
        "annot_file": os.path.join(RAW_DIR, "GPL97.annot"),
    },
]

GENE_AVERAGED_OUT = os.path.join(OUTPUT_DIR, "GSE3790_Huntington_dataset.csv")

ID_COL = "Sample_ID"
TARGET_COL = "Huntingtons_Disease"


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


def process_geo_matrix(matrix_file, probe_to_gene):
    sample_ids = []
    sample_data = {}

    with open(matrix_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("!Series_sample_id"):
                sample_ids = re.findall(r"GSM\d+", line)
                for gsm in sample_ids:
                    sample_data[gsm] = {"Sample_ID": gsm, "Huntingtons_Disease": None}

            elif line.startswith("!Sample_characteristics_ch1") or line.startswith("!Sample_title"):
                parts = line.split("\t")
                for i, val in enumerate(parts[1:]):
                    if i < len(sample_ids):
                        gsm = sample_ids[i]
                        if "hd" in val.lower() or "huntington" in val.lower():
                            sample_data[gsm]["Huntingtons_Disease"] = 1
                        elif "control" in val.lower() or "normal" in val.lower():
                            sample_data[gsm]["Huntingtons_Disease"] = 0

            elif not line.startswith("!"):
                parts = line.split("\t")
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
    #extract probe-level data 
    combined_by_gsm = {}
    for platform in PLATFORMS:
        probe_to_gene = load_annotation(platform["annot_file"], MARKERS)
        data = process_geo_matrix(platform["matrix_file"], probe_to_gene)
        for gsm, row in data.items():
            combined_by_gsm.setdefault(gsm, {}).update(row)

    combined = [row for row in combined_by_gsm.values() if row.get("Huntingtons_Disease") is not None]
    df = pd.DataFrame(combined)
    cols = [ID_COL, TARGET_COL] + [c for c in df.columns if c not in (ID_COL, TARGET_COL)]
    df = df[cols]

    #average probe columns into one column per gene
    feature_cols = [c for c in df.columns if c not in (ID_COL, TARGET_COL)]
    probes_per_gene = {}
    for col in feature_cols:
        gene = get_gene_symbol(col)
        probes_per_gene.setdefault(gene, []).append(col)

    print("Probes per gene:")
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