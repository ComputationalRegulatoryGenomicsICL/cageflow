#!/usr/bin/env bash
# Usage:
#   softclip5_se_exact.sh IN.bam OUT.bam [N] [T] [MOTIF]
# Example:
#   softclip5_se_exact.sh in.bam out.se.5p3S.ATG.bam 3 8 ATG
#
# - N: exact length of the 5' soft-clip (default 1)
# - T: number of threads (default 4)
# - motif (optional): expected 5' bases (e.g., ATG). If given, only reads
#     with 5' soft-clipped sequence that matches given motif are kept.
#     Matching is strand-aware: reverse reads compare to RC(motif).
# - For reverse-strand alignments, SAM stores SEQ reverse-complemented;
#   we therefore compare the last N bases of SEQ to revcomp(MOTIF).

IN=${1:?}; OUT=${2:?}; N=${3:-1}; T=${4:-4}; MOTIF=${5:-}

if [[ -n "$MOTIF" && ${#MOTIF} -ne $N ]]; then
  echo "Error: motif length (${#MOTIF}) must equal N (got N='$N')" >&2
  exit 2
fi

# Function to reverse-complement a DNA motif (ACGTN only).
revcomp() {
  local s="$1"
  s=$(printf "%s" "$s" | tr 'acgtn' 'ACGTN' | rev | tr 'ACGTN' 'TGCAN')
  printf "%s" "$s"
}

RC_MOTIF="$(revcomp "$MOTIF")"

samtools merge -O bam - \
  <(
    # Forward strand: CIGAR starts with ^NS
    samtools view -h -F 20 "$IN" | \
      awk -v N="$N" -v M="$MOTIF" 'BEGIN{
        OFS="\t"
        hasToMatch = (M != "")
      }
      /^@/ { print; next }
      {
        if (!match($6, "^"N"S")) next
        if (!hasToMatch) { print; next }
        seq5 = substr($10, 1, N)
        if (seq5 == M) print
      }' | \
      samtools view -b -
  ) \
  <(
    # Reverse strand: CIGAR ends with NS$
    samtools view -h -f 16 "$IN" | \
      awk -v N="$N" -v MRC="$RC_MOTIF" 'BEGIN{
        OFS="\t"
        hasToMatch = (MRC != "")
      }
      /^@/ { print; next }
      {
        if (!match($6, N "S$")) next
        if (!hasToMatch) { print; next }
        seq3 = substr($10, length($10) - N + 1, N)
        if (seq3 == MRC) print
      }' | \
      samtools view -b -
  ) | \
samtools sort -@ "$T" -O bam -o "$OUT"
