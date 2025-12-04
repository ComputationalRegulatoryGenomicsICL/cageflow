#!/usr/bin/env bash
# Usage:
#   softclip5_se_exact.sh IN.bam OUT.bam [N] [T] [FIFTH]
# Examples:
# - N>0: exact length of the 5' soft-clip, optionally with motif.
#   softclip5_se_exact.sh in.bam out.se.5p3S.ATG.bam 1 8
#   softclip5_se_exact.sh in.bam out.se.5p3S.ATG.bam 3 8 ATG
# - N==0: select reads with NO soft-clip at the 5'-end of the read. Soft-clips at the 3'-end are allowed.
#     a) with motif: if length(motif) <= length(5'-end CIGAR match), then select reads with this motif
#   softclip5_se_exact.sh in.bam out.no5S.motif.bam 0 8 ATG
#     b) without motif, but with M: minimum required length of the 5'-end CIGAR match
#   softclip5_se_exact.sh in.bam out.no5S.M20.bam 0 8 20
#     c) without motif and without M: select reads with no 5' soft-clip
#   softclip5_se_exact.sh in.bam out.no5S.bam 0 8
#
# Arguments:
# - N: length of 5' soft-clip (CIGAR) (default 1):
#   - If N>0: exact length of 5' soft-clip
#   - If N==0: select reads WITHOUT 5' soft-clip. 3' soft-clips are allowed.
# - T: number of threads (default 4)
# - FIFTH (optional):
#   - if N>0: treated as motif; length must equal N
#   - if N==0:
#        - if the argument contains only digits, then M: minimum length of the 5'-end CIGAR match;
#        - otherwise, motif string; len(motif) is checked per read: length(motif)<=length(CIGAR matches).
#
# Motif checking is strand-aware:.
# - For reverse-strand alignments, SAM stores SEQ reverse-complemented;
# - We therefore compare the last N bases of SEQ to revcomp(MOTIF).

IN=${1:?}; OUT=${2:?}; N=${3:-1}; T=${4:-4}; FIFTH=${5:-}

MOTIF=""
M_MIN=""

# Parse the 5th argument depending on N
if [[ $N -eq 0 ]]; then
  if [[ -n "$FIFTH" && "$FIFTH" =~ ^[0-9]+$ ]]; then
    # Integer to M (minimum 5'-match length in CIGAR)
    M_MIN="$FIFTH"
  else
    # Otherwise to motif
    MOTIF="$FIFTH"
  fi
else
  # N>0: always treat the 5th argument as motif
  MOTIF="$FIFTH"
fi

# For N>0, if motif is given, its length must equal N
if [[ -n "$MOTIF" && $N -ne 0 && ${#MOTIF} -ne $N ]]; then
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
        hasMotif    = (M != "")
        hasMinMatch = (M_MIN != "")
        motifLen    = length(M)
      }
      /^@/ { print; next }
      {
        if (N == 0) {
          if (match($6, "^[0-9]+S")) next

          lenMatch = 0
          if (match($6, "^([0-9]+)M", a)) {
            lenMatch = a[1] + 0
          } else {
            if (hasMotif || hasMinMatch) next
          }

          if (hasMotif) {
            if (lenMatch < motifLen) next
            seq5 = substr($10, 1, motifLen)
            if (seq5 == M) print
            next
          }

          if (hasMinMatch) {
            if (lenMatch >= M_MIN) print
            next
          }

          print
          next
        }

        if (!match($6, "^"N"S")) next
        if (!hasMotif) { print; next }

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
        hasMotif    = (MRC != "")
        hasMinMatch = (M_MIN != "")
        motifLen    = length(MRC)
      }
      /^@/ { print; next }
      {
        if (N == 0) {
          if (match($6, "[0-9]+S$")) next

          lenMatch = 0
          if (match($6, "([0-9]+)M$", a)) {
            lenMatch = a[1] + 0
          } else {
            if (hasMotif || hasMinMatch) next
          }

          if (hasMotif) {
            if (lenMatch < motifLen) next
            seq3 = substr($10, length($10) - motifLen + 1, motifLen)
            if (seq3 == MRC) print
            next
          }

          if (hasMinMatch) {
            if (lenMatch >= M_MIN) print
            next
          }

          print
          next
        }

        if (!match($6, N "S$")) next
        if (!hasMotif) { print; next }

        seq3 = substr($10, length($10) - N + 1, N)
        if (seq3 == MRC) print
      }' | \
      samtools view -b -
  ) | \
samtools sort -@ "$T" -O bam -o "$OUT"
