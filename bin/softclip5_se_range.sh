#!/usr/bin/env bash
# Usage: softclip5_se_range.sh IN.bam OUT.bam [M] [N] [T] [BASE]
#
# - M: positive number (default 1) 
# - N: positive number (default 3) 
#      Keeps reads with 5' soft-clip length x where M <= x <= N.
# - T: number of threads (default 4)
# - BASE (optional): One of A/C/G/T/N.
#   If BASE (A/C/G/T/N) is provided, the soft-clipped segment (from M to N)
#   must be a homopolymer of BASE in the original read orientation
#   (in reverse reads BASE is complement).
#
# Examples:
#   softclip5_se_range.sh in.bam out.se.5p1to3S.bam 1 3 8
#   softclip5_se_range.sh in.bam out.se.5p2to5S.Gpoly.bam 2 5 8 G

IN=${1:?}; OUT=${2:?}; M=${3:-1}; N=${4:-3}; T=${5:-4}; BASE=${6:-}

# Ensure M and N are positive integers
if ! [[ "$M" =~ ^[0-9]+$ && "$N" =~ ^[0-9]+$ && "$M" -ge 1 && "$N" -ge 1 ]]; then
  echo "Error: M and N must be positive integers (got M='$M', N='$N')." >&2
  exit 2
fi

# Ensure M <= N
if [[ "$M" -gt "$N" ]]; then
  echo "Error: Require M <= N (got M=$M, N=$N)." >&2
  exit 2
fi

# Normalise and validate BASE is uppercase (if provided)
if [[ -n "$BASE" ]]; then
  BASE=$(printf "%s" "$BASE" | tr 'a-z' 'A-Z')
  if [[ ${#BASE} -ne 1 || ! "$BASE" =~ ^[ACGTN]$ ]]; then
    echo "Error: BASE must be one of A/C/G/T/N (got '$BASE')." >&2; exit 2
  fi
fi

# Complement for reverse-strand homopolymer check
comp_base() {
  case "$1" in
    A) echo T;;
    C) echo G;;
    G) echo C;;
    T) echo A;;
    N) echo N;;
    *) echo N;;
  esac
}
BASE_RC=""
if [[ -n "$BASE" ]]; then
  BASE_RC=$(comp_base "$BASE")
fi

samtools merge -O bam - \
  <(
    # Forward strand:
    #   - Keep mapped non-reverse reads (-F 20 filters out unmapped and reverse).
    #   - Require CIGAR to start with xS, and M <= x <= N.
    #   - If BASE is set, the entire 5' soft-clip (first x bases of SEQ) must be a BASE homopolymer.
    #   - x = m[1] + 0 makes this value a number
    samtools view -h -F 20 "$IN" | \
      awk -v M="$M" -v N="$N" -v BASE="$BASE" '
        BEGIN{ OFS="\t"; hasBASE=(BASE!=""); upperBASE=toupper(BASE) }
        /^@/ { print; next }
        {
          if (!match($6, "^([0-9]+)S", m)) next
          x = m[1] + 0 
          if (x < M || x > N) next
          if (!hasBASE) { print; next }

          sc_full = substr($10, 1, x)
          t = toupper(sc_full)
          if (t ~ ("^" upperBASE "+$")) print
        }' | \
      samtools view -b -
  ) \
  <(
    # Reverse strand: CIGAR ends with xS, M<=x<=N
    # SEQ is reverse-complemented, so original 5' = last x bases; compare to complement(BASE)
    samtools view -h -f 16 "$IN" | \
      awk -v M="$M" -v N="$N" -v BASE_RC="$BASE_RC" '
        BEGIN{ OFS="\t"; hasBASE=(BASE_RC!=""); upperBASE=toupper(BASE_RC) }
        /^@/ { print; next }
        {
          if (!match($6, "([0-9]+)S$", m)) next
          x = m[1] + 0
          if (x < M || x > N) next
          if (!hasBASE) { print; next }

          sc_full = substr($10, length($10) - x + 1, x)
          t = toupper(sc_full)
          if (t ~ ("^" upperBASE "+$")) print
        }' | \
      samtools view -b -
  ) | \
samtools sort -@ "$T" -O bam -o "$OUT"