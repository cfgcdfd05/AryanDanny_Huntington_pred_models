
import os
import pandas as pd

INPUT_PATH = "GSE64810/GSE64810_norm_counts_CLEANED.txt"
OUTPUT_DIR = "datasets"
OUTPUT_PATH = os.path.join(OUTPUT_DIR, "GSE64810_Huntington_dataset.csv")


# Gene symbol -> Ensembl gene ID (versionless). Verify against your GTF if unsure.
GENE_TO_ENSEMBL = {
    "TUT7":    "ENSG00000083223",
    "PPP1CC":  "ENSG00000186298",  
    "CEBPB":   "ENSG00000172216",
    "CEBPD":   "ENSG00000221869",
    "FTO":     "ENSG00000140718",
    "METTL16": "ENSG00000127804",
    "IGF2BP3": "ENSG00000136231",
    "YTHDF1":  "ENSG00000149658",
    "YTHDF2":  "ENSG00000198492",
    "YTHDF3":  "ENSG00000185728",
    "YTHDC1":  "ENSG00000083896",  
    "YTHDC2":  "ENSG00000047188",
}

ENSEMBL_TO_GENE = {v: k for k, v in GENE_TO_ENSEMBL.items()}


def strip_version(ensembl_id):
    """ENSG00000000003.10 -> ENSG00000000003"""
    return ensembl_id.split(".")[0]


def label_from_sample_id(sample_id):
    """H_xxxx -> 1 (Huntington's), C_xxxx -> 0 (control)."""
    prefix = sample_id.strip().split("_")[0].upper()
    if prefix == "H":
        return 1
    elif prefix == "C":
        return 0
    return None


def main():
    # First column has no header (it's the gene ID column) -> read with index_col=0
    df = pd.read_csv(INPUT_PATH, sep="\t", index_col=0)

    # Strip version suffixes from the Ensembl IDs in the index
    df.index = [strip_version(idx) for idx in df.index]

    # Find which of our target genes are present
    found = {gene: ensg for gene, ensg in GENE_TO_ENSEMBL.items() if ensg in df.index}
    missing = [gene for gene in GENE_TO_ENSEMBL if gene not in found]

    print(f"Found {len(found)}/{len(GENE_TO_ENSEMBL)} biomarkers in the matrix:")
    for gene, ensg in found.items():
        print(f"  {gene}: {ensg}")
    if missing:
        print(f"\nWARNING - not found in matrix (check gene ID mapping): {missing}")

    # Subset to target genes, transpose so samples become rows
    subset = df.loc[list(found.values())].T
    subset.columns = [ENSEMBL_TO_GENE[c] for c in subset.columns]

    subset.insert(0, "Sample_ID", subset.index)
    subset.insert(1, "Huntingtons_Disease", [label_from_sample_id(s) for s in subset.index])
    subset = subset.reset_index(drop=True)

    unlabeled = subset[subset["Huntingtons_Disease"].isna()]
    if len(unlabeled) > 0:
        print(f"\nWARNING - {len(unlabeled)} samples had unrecognized ID prefixes and got no label:")
        print(unlabeled["Sample_ID"].tolist())

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    subset.to_csv(OUTPUT_PATH, index=False)

    print(f"\nFinal shape: {subset.shape}")
    print(f"Class balance:\n{subset['Huntingtons_Disease'].value_counts()}")
    print(f"Saved: {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
