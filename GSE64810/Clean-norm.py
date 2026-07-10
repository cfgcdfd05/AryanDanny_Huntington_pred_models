INPUT_PATH = "./GSE64810/GSE64810_mlhd_DESeq2_norm_counts_adjust.txt"
CLEANED_PATH = "./GSE64810/GSE64810_norm_counts_CLEANED.txt"

### This file removes corrupted entries from the original files

TARGET_ENSG = {
    "ENSG00000083223", "ENSG00000186298", "ENSG00000172216", "ENSG00000221869",
    "ENSG00000140718", "ENSG00000127804", "ENSG00000136231", "ENSG00000149658",
    "ENSG00000198492", "ENSG00000185728", "ENSG00000083896", "ENSG00000047188",
}

with open(INPUT_PATH, encoding="utf-8") as f:
    lines = f.readlines()

header_fields = lines[0].rstrip("\n").split("\t")
expected_count = len(header_fields)
print(f"Header field count: {expected_count}")
print(f"Total lines: {len(lines)}")

good_lines = [lines[0]]
bad_lines = []

for i, line in enumerate(lines[1:], start=2):
    fields = line.rstrip("\n").split("\t")
    if len(fields) == expected_count:
        good_lines.append(line)
    else:
        bad_lines.append((i, len(fields), fields))

print(f"\nGood lines: {len(good_lines) - 1}")
print(f"Bad lines: {len(bad_lines)}")

for line_num, count, fields in bad_lines[:5]:
    gene_id = fields[0].split(".")[0] if fields else "?"
    flag = " <-- ONE OF YOUR 12 TARGET GENES" if gene_id in TARGET_ENSG else ""
    print(f"\n--- Line {line_num}: {count} fields (expected {expected_count}){flag} ---")
    print(f"First field (gene ID): {repr(fields[0])}")
    print(f"Last field: {repr(fields[-1])}")
    print(f"Full field list: {fields}")

affected_targets = [
    fields[0] for _, _, fields in bad_lines
    if fields and fields[0].split(".")[0] in TARGET_ENSG
]
if affected_targets:
    print(f"\n*** WARNING: {len(affected_targets)} of your 12 target genes are on malformed lines: {affected_targets} ***")
    print("These would be silently dropped by the cleanup below - do not proceed without inspecting them.")
else:
    print("\nNone of your 12 target genes are on malformed lines - safe to drop bad lines and continue.")

# Write a cleaned file with only well-formed lines
with open(CLEANED_PATH, "w", encoding="utf-8") as f:
    f.writelines(good_lines)

print(f"\nCleaned file written to: {CLEANED_PATH}")
print(f"Dropped {len(bad_lines)} malformed line(s) — line numbers: {[ln for ln, _, _ in bad_lines]}")