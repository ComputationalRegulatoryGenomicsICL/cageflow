#!/usr/bin/env python3
import sys
import pysam

SOFT = 4  # 'S' soft-clip in pysam CIGAR codes
MATCH = 0  # 'M' (alignment match) in pysam CIGAR codes

def revcomp(seq: str) -> str:
    """
    Reverse-complement a DNA sequence (ACGTN, case-insensitive).
    """
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]

def has_5prime_condition(aln, N, motif=None, rc_motif=None, min_match=None):
    """
    Check 5'-end conditions on a single alignment.

    For all cases:
      - Reads must be mapped.

    5'-end definition:
      - Forward reads: 5'-end corresponds to the first CIGAR element.
      - Reverse reads: 5'-end corresponds to the last CIGAR element.

    If N == 0:
      - Require NO 5'-end soft-clip (5'-end CIGAR op != SOFT).
      - 3'-end soft-clips are allowed.
      - Additionally:
          - If 'motif' is provided:
              - Let L be the length of the 5'-end CIGAR match (M op).
              - Require L >= len(motif).
              - Forward:  first len(motif) bases of SEQ == motif (case-insensitive).
              - Reverse:  last len(motif) bases of SEQ == rc_motif (case-insensitive).
          - Else if 'min_match' is provided:
              - Require L >= min_match.
          - Else:
              - Only the "no 5'-soft-clip" condition is enforced.

    If N > 0:
      - Require an exact N-bp soft-clip at the 5'-end:
          - Forward: first CIGAR element is (SOFT, N).
          - Reverse: last  CIGAR element is (SOFT, N).
      - If 'motif' is provided:
          - Forward:  first N bases of SEQ == motif (case-insensitive).
          - Reverse:  last  N bases of SEQ == rc_motif (case-insensitive).
      - If 'motif' is not provided:
          - Only the soft-clip condition matters.
    """
    if aln.is_unmapped:
        return False

    cig = aln.cigartuples

    # Determine the 5'-end CIGAR element 
    if aln.is_reverse:
        op_5p, len_5p = cig[-1]
    else:
        op_5p, len_5p = cig[0]

    # N == 0: select reads with NO 5'-end soft-clip (but allow soft-clip elsewhere)
    if N == 0:
        # Reject reads with 5'-end soft-clip
        if op_5p == SOFT:
            return False

        # If neither motif nor min_match is specified, only enforce "no 5'-soft-clip"
        if motif is None and min_match is None:
            return True

        # We need the length of the 5'-end match (M) to apply motif or min_match
        len_match_5p = 0
        if aln.is_reverse:
            # 5'-end is the last CIGAR element; require it to be a match
            op_last, len_last = cig[-1]
            if op_last == MATCH:
                len_match_5p = len_last
        else:
            # 5'-end is the first CIGAR element; require it to be a match
            op_first, len_first = cig[0]
            if op_first == MATCH:
                len_match_5p = len_first

        # If we require motif or min_match but have no initial match, fail
        if len_match_5p == 0 and (motif is not None or min_match is not None):
            return False

        seq = aln.query_sequence
        # Check if seq is not None
        if not seq:
            return False

        if motif is not None:
            motif_len = len(motif)
            # Require 5'-end match to be at least motif length
            if len_match_5p < motif_len:
                return False

            if aln.is_reverse:
                # SEQ is stored reverse-complemented, so 5'-end bases are at the 3'-end of SEQ
                seq_part = seq[-motif_len:]
                return seq_part.upper() == rc_motif.upper()
            else:
                # Forward: 5'-end bases are at the beginning of SEQ
                seq_part = seq[:motif_len]
                return seq_part.upper() == motif.upper()

        if min_match is not None:
            # Only minimal 5'-match length is enforced
            return len_match_5p >= min_match

        # Otherwise (should not be reached)
        return False

    # N > 0: require an exact N-bp 5'-end soft-clip
    if op_5p != SOFT or len_5p != N:
        return False

    # If no motif is provided, only the soft-clip is required
    if motif is None:
        return True

    seq = aln.query_sequence
    if not seq:
        return False

    if aln.is_reverse:
        # SEQ is revcomp; 5'-end of the read is at the end of SEQ
        seq_part = seq[-N:]
        return seq_part.upper() == rc_motif.upper()
    else:
        # Forward: 5'-end of the read is at the beginning of SEQ
        seq_part = seq[:N]
        return seq_part.upper() == motif.upper()


def process_group(records, N, motif, rc_motif, out_bam):
    """
    Process a group of alignments with the same query_name.

    records: list of AlignedSegment objects with the same query_name.

    Logic:
      - Identify read1 alignments (is_read1) that satisfy the 5'-end conditions:
            - N > 0:
                - exact N-bp 5'-end soft-clip
                - and if motif is specified: matching motif requirements
            - N == 0:
                - no 5'-end soft-clip
                - and if motif is specified: motif constraints on the 5'-end match
                - or if min_match is specified: minimal 5'-end match length
      - If there is at least one R1 that passes:
        - write all passing R1 alignments
        - and write all corresponding R2 alignments whose (reference_id, reference_start)
        match any (next_reference_id, next_reference_start) of the passing R1s.
      - Otherwise, nothing is written for this group.
    """
    if not records:
        return

    r1_selected = [
        r for r in records
        if r.is_read1 and has_5prime_condition(r, N, motif, rc_motif, min_match)
    ]
    if not r1_selected:
        return

    # Mate R2s coordinates for matching R1s
    # RNEXT and PNEXT in SAM are next_reference_id and next_reference_start in pysam
    mate_coords = {
        (r.next_reference_id, r.next_reference_start)
        for r in r1_selected
        if r.next_reference_id >= 0 and r.next_reference_start >= 0
    }

    # Now select only those R2 alignments that map to those coordinates
    r2_selected = [
        r for r in records
        if r.is_read2 and (r.reference_id, r.reference_start) in mate_coords
    ]

    for r in r1_selected + r2_selected:
        out_bam.write(r)

def main():
    if len(sys.argv) < 4:
        sys.stderr.write(
            "Usage:\n"
            "  python3 softclip5_pe_exact.py IN.bam OUT.bam N [MOTIF_or_M]\n\n"
            "  IN.bam        : name-sorted input BAM.\n"
            "  OUT.bam       : output BAM.\n"
            "  N             : 5'-end soft-clip length on read1.\n"
            "                  - If N > 0: require an exact N-bp 5'-end soft-clip.\n"
            "                  - If N = 0: require NO 5'-end soft-clip on read1.\n"
            "  MOTIF_or_M    : optional fourth argument with two modes:\n"
            "                  - If N > 0:\n"
            "                       - interpreted as MOTIF (string);\n"
            "                       - len(MOTIF) must equal N;\n"
            "                       - forward R1 : first N bases of SEQ == MOTIF;\n"
            "                       - reverse R1 : last  N bases of SEQ == revcomp(MOTIF).\n"
            "                  - If N = 0:\n"
            "                       - if the argument is all digits -> M (integer):\n"
            "                             minimal 5'-end CIGAR match length on R1;\n"
            "                       - otherwise -> MOTIF (string):\n"
            "                             let L be the 5'-end CIGAR match length;\n"
            "                             require L >= len(MOTIF) and matching\n"
            "                             5'-end bases (forward) / revcomp bases (reverse).\n\n"
            "Notes:\n"
            "  - Input BAM must be name-sorted.\n"
            "  - Output contains:\n"
            "       - R1 alignments that satisfy the conditions above;\n"
            "       - their corresponding R2 alignments (matched by mate coordinates).\n\n"
            "Examples:\n"
            "  # N = 0, no 5'-soft-clip on R1, no additional motif/length filter\n"
            "  python3 softclip5_pe_exact.py in.namesort.bam out.no5S.bam 0\n\n"
            "  # N = 0, no 5'-soft-clip on R1, minimal 5'-match length = 20\n"
            "  python3 softclip5_pe_exact.py in.namesort.bam out.no5S.M20.bam 0 20\n\n"
            "  # N = 0, no 5'-soft-clip on R1, motif ATG at 5'-end (length <= 5'-match)\n"
            "  python3 softclip5_pe_exact.py in.namesort.bam out.no5S.ATG.bam 0 ATG\n\n"
            "  # N = 3, require 3-bp 5'-soft-clip on R1, no motif\n"
            "  python3 softclip5_pe_exact.py in.namesort.bam out.5p3S.bam 3\n\n"
            "  # N = 3, require 3-bp 5'-soft-clip on R1 and motif ATG\n"
            "  python3 softclip5_pe_exact.py in.namesort.bam out.5p3S.ATG.bam 3 ATG\n\n"
        )
        sys.exit(1)

    in_bam_path = sys.argv[1]
    out_bam_path = sys.argv[2]
    N = int(sys.argv[3])

    fifth = sys.argv[4] if len(sys.argv) >= 5 else None

    motif = None
    rc_motif = None
    min_match = None

    # Parse the 4th argument depending on N
    if N == 0:
        if fifth is not None:
            if fifth.isdigit():
                # Numeric argument -> minimal 5'-end match length
                min_match = int(fifth)
            else:
                # Non-numeric argument -> motif string
                motif = fifth
                rc_motif = revcomp(motif)
    else:
        # N > 0: interpret the 4th argument as motif, if provided
        if fifth is not None:
            motif = fifth
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
            process_group(group, N, motif, rc_motif, min_match, out_bam)
            current_qname = qn
            group = [aln]

    # Process last group
    if group:
        process_group(group, N, motif, rc_motif, min_match, out_bam)

    in_bam.close()
    out_bam.close()

if __name__ == "__main__":
    main()
