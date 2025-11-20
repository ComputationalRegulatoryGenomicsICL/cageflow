#!/usr/bin/env python3
import sys
import pysam

SOFT = 4  # soft-clip in pysam CIGAR codes

def comp_base(base: str) -> str:
    """
    Complement a of a single DNA base (A/C/G/T/N), case-insensitive.
    """
    base = base.upper()
    table = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return table.get(base, "N")

def has_5prime_softclip_homopolymer_upto(aln, M, N, base=None, rc_base=None):
    """
    Return True if this alignment has:
      - a 5' soft-clip of length x (orientation-aware: 5' of read),
      - where M <= x <= N,
      - and if base is provided (A/C/G/T/N):
          - forward reads: entire 5' soft-clip (first x bases of SEQ) is homopolymer of base
          - reverse reads: entire 5' soft-clip (last x bases of SEQ) is homopolymer of rc_base

      - For forward (not is_reverse): 5' soft-clip is the first CIGAR element (cig[0] == xS).
      - For reverse (is_reverse): 5' soft-clip is the last CIGAR element (cig[-1] == xS).
      - As SEQ for reverse reads were reverse-complemented we compare last soft-clipped bases to complement(base).
    """
    if aln.is_unmapped:
        return False

    cig = aln.cigartuples

    # get the 5' operation depending on orientation
    op, length = cig[-1] if aln.is_reverse else cig[0]

    # must be soft-clip with length in [M, N]
    if op != SOFT:
        return False
    x = int(length)
    if x < M or x > N:
        return False

    # if no base given, we return True
    if not base:
        return True

    seq = aln.query_sequence
    if not seq:
        return False

    # extract the soft-clipped segment and test homopolymer
    if aln.is_reverse:
        # reverse 5' == last x bases
        sc = seq[-x:]
        target_base = rc_base.upper()
    else:
        # forward 5' == first x bases
        sc = seq[:x]
        target_base = base.upper()

    # all bases in soft-clip must match target_base (homopolymer)
    return all(b.upper() == target_base for b in sc)


def process_group(records, M, N, base, rc_base, out_bam):
    """
    records: list of AlignedSegment with same query_name (paired-end).

    - If any read1 in this group has a 5' soft-clip of length x in [M, N]
      (and passes the optional BASE homopolymer check):
        - write only those read1 alignments that match
        - write all corresponding read2 alignments (based on RNEXT/PNEXT of R1 vs RNAME/RPOS of R2)
    - Else: write nothing for this group.
    """
    if not records:
        return

    # R1s that match the soft-clip + homopolymer criteria
    r1_soft = [
        r for r in records
        if r.is_read1 and has_5prime_softclip_homopolymer_upto(r, M, N, base, rc_base)
    ]
    if not r1_soft:
        return

    # Mate R2s coordinates for matching R1s
    # RNEXT / PNEXT in SAM => next_reference_id / next_reference_start in pysam
    mate_coords = {
        (r.next_reference_id, r.next_reference_start)
        for r in r1_soft
        if r.next_reference_id >= 0 and r.next_reference_start >= 0
    }

    # Now select only those R2 alignments that map to those coordinates
    r2_all = [
        r for r in records
        if r.is_read2 and (r.reference_id, r.reference_start) in mate_coords
    ]

    for r in r1_soft + r2_all:
        out_bam.write(r)


def main():
    if len(sys.argv) < 5:
        sys.stderr.write(
            "Usage:\n"
            "  python3 softclip5_pe_upto.py IN.bam OUT.bam [M] [N] [BASE]\n\n"
            "Description:\n"
            "  Select paired-end reads where read1 has a 5' soft-clip of length x\n"
            "  such that M <= x <= N. If BASE (A/C/G/T/N) is provided, the entire\n"
            "  soft-clipped 5' segment must be a homopolymer:\n"
            "    - forward:  5' soft-clip bases == BASE\n"
            "    - reverse:  5' soft-clip bases == complement(BASE)\n\n"
            "  - Input BAM must be name-sorted.\n"
            "  - All matching read1 alignments are written, plus their matching\n"
            "    read2 mates (matched by RNEXT/PNEXT <-> RNAME/RPOS).\n\n"
            "Example:\n"
            "  python3 softclip5_pe_upto.py in.namesort.bam out.pe.5p1to3S.bam 1 3\n"
            "  python3 softclip5_pe_upto.py in.namesort.bam out.pe.5p2to5S.Gpoly.bam 2 5 G\n"
        )
        sys.exit(1)

    in_bam_path = sys.argv[1]
    out_bam_path = sys.argv[2]

    try:
        M = int(sys.argv[3])
        N = int(sys.argv[4])
    except ValueError:
        sys.stderr.write("Error: M and N must be integers.\n")
        sys.exit(2)

    if M < 1 or N < 1:
        sys.stderr.write("Error: M and N must be positive integers.\n")
        sys.exit(2)
    if M > N:
        sys.stderr.write(f"Error: require M <= N (got M={M}, N={N}).\n")
        sys.exit(2)

    base = sys.argv[5].upper() if len(sys.argv) >= 6 else None
    rc_base = None

    if base:
        if len(base) != 1 or base not in "ACGTN":
            sys.stderr.write(
                f"Error: BASE must be one of A/C/G/T/N (got '{base}').\n"
            )
            sys.exit(2)
        rc_base = comp_base(base)

    in_bam = pysam.AlignmentFile(in_bam_path, "rb")
    out_bam = pysam.AlignmentFile(out_bam_path, "wb", template=in_bam)

    current_qname = None
    group = []

    # Assumes name-sorted BAM
    for aln in in_bam.fetch(until_eof=True):
        qn = aln.query_name
        if current_qname is None:
            current_qname = qn
            group = [aln]
        elif qn == current_qname:
            group.append(aln)
        else:
            process_group(group, M, N, base, rc_base, out_bam)
            current_qname = qn
            group = [aln]

    # process last group
    if group:
        process_group(group, M, N, base, rc_base, out_bam)

    in_bam.close()
    out_bam.close()


if __name__ == "__main__":
    main()
