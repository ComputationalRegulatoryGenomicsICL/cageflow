#!/usr/bin/env python3
import sys
import pysam

SOFT = 4  # soft-clip in pysam CIGAR codes

def revcomp(seq: str) -> str:
    """
    Reverse-complement a DNA sequence (ACGTN, case-insensitive).
    """
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]

def has_5prime_softclip_with_motif(aln, N, motif=None, rc_motif=None):
    """
    Return True if this alignment has:
      - an exact N-bp soft-clip at its 5' end (orientation-aware), and
      - if motif is provided:
          - forward: first N bases of SEQ == motif
          - reverse: last N bases of SEQ == RC(motif)
    """
    if aln.is_unmapped:
        return False

    cig = aln.cigartuples

    # Check 5' soft-clip of length N
    op, length = cig[-1] if aln.is_reverse else cig[0]

    if op != SOFT or length != N:
        return False

    # If no motif provided, only 5' soft-clip condition matters
    if not motif:
        return True

    seq = aln.query_sequence
    # Check if seq is not None
    if not seq:
        return False

    if aln.is_reverse:
        # SEQ is revcomp - 5' of read corresponds to last N bases of SEQ
        seq_part = seq[-N:]
        return seq_part.upper() == rc_motif.upper()
    else:
        # Forward: 5' of read == first N bases of SEQ
        seq_part = seq[:N]
        return seq_part.upper() == motif.upper()

def process_group(records, N, motif, rc_motif, out_bam):
    """
    records: list of AlignedSegment with same query_name

    - If any read1 has 5' N-bp soft-clip (and matches motif if given):
        - write only those read1 alignments that match
        - write all corresponding read2 alignments (based on RNEXT/PNEXT of R1 vs RNAME/RPOS of R2)
    - Else: write nothing.
    """
    if not records:
        return

    r1_soft = [
        r for r in records
        if r.is_read1 and has_5prime_softclip_with_motif(r, N, motif, rc_motif)
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
    if len(sys.argv) < 4:
        sys.stderr.write(
            "Usage:\n"
            "  python3 softclip5_pe_exact.py IN.bam OUT.bam [N] [MOTIF]\n\n"
            "  - N: exact length of the 5' soft-clip.\n"
            "  - MOTIF (optional): expected 5' bases on read1.\n"
            "    - forward:  first N bases of SEQ == MOTIF\n"
            "    - reverse:  last N bases of SEQ == revcomp(MOTIF)\n"
            "  - Input BAM must be name-sorted.\n"
            "  - All matching read1 alignments are written, plus their matching\n"
            "    read2 mates (matched by RNEXT/PNEXT <-> RNAME/RPOS).\n\n"
            "Example:\n"
            "  python3 softclip5_pe_exact.py in.namesort.bam out.pe.5p3S.ATG.bam 3 ATG\n\n"
        )
        sys.exit(1)

    in_bam_path = sys.argv[1]
    out_bam_path = sys.argv[2]
    N = int(sys.argv[3])

    motif = sys.argv[4] if len(sys.argv) >= 5 else None
    rc_motif = None

    if motif:
        if len(motif) != N:
            sys.stderr.write(
                f"Error: motif length ({len(motif)}) must equal N (got N='{N}')\n"
            )
            sys.exit(2)
        rc_motif = revcomp(motif)

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
            process_group(group, N, motif, rc_motif, out_bam)
            current_qname = qn
            group = [aln]

    # Process last group
    if group:
        process_group(group, N, motif, rc_motif, out_bam)

    in_bam.close()
    out_bam.close()

if __name__ == "__main__":
    main()
