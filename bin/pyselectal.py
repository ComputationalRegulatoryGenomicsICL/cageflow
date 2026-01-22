#!/usr/bin/env python3

MANUAL = """\
pyselectal — filter alignments by 5'-end soft-clipping or matched prefix

Complete manual: https://github.com/nikitin-p/pyselectal

This script filters reads in a BAM file based on their 5'-end:
    - presence/absence of 5'-soft-clipped bases
    - length of the 5'-soft-clip (exact or within a range)
    - optional 5'-end sequence (motif / prefix / homopolymer)

It supports both single-end (SE) and paired-end (PE) reads.
If --paired is used, input should be name-sorted OR use --sort.

"""
################################################################################
# Usage:
#
#     python pyselectal.py [options] in.namesort.bam out.bam
#
################################################################################
# Options:
#
#   -n / --min-softclip      (int, required)
#       Minimum number of 5'-soft-clipped bases.
#
#   -m / --max-softclip      (int, required)
#       Maximum number of 5'-soft-clipped bases.
#
#   Together, -n and -m define the mode:
#       1) Soft-clip of exact length            : n = m > 0
#       2) Mapped 5'-end (no S)                 : n = m = 0
#       3) Soft-clip of length within a range   : n < m
#
#   -p / --prefix          (optional, meaning depends on mode)
#
#       If n = m > 0  (exact soft-clip):
#           -p is a prefix of length n.
#           Forward: soft-clipped seq == prefix.
#           Reverse: soft-clipped seq == revcomp(prefix).
#
#       If n = m = 0  (mapped 5'-end):
#           -p is a prefix at the mapped 5' end.
#           Forward: read sequence starts with prefix.
#           Reverse: read sequence ends with revcomp(prefix).
#
#       If n < m      (range mode):
#           -p is a single base A/C/G/T/N, defining a homopolymer.
#           Forward: all 5'-soft-clipped bases == base.
#           Reverse: all 5'-soft-clipped bases == complement(base).
#
#       Validation:
#           - Exact mode: prefix length must equal n.
#           - Range mode: prefix must be exactly one base A/C/G/T/N.
#
#   -k                     (int, optional; only meaningful when n = m = 0)
#       If --prefix is NOT given:
#           - keep reads whose 5'-end CIGAR operation is MATCH (M) and M-length >= k
#           - if k==0: keep reads with 5'-end operation == M (i.e., no 5' soft-clip)
#
#       If --prefix IS given:
#           1) if k <= len(prefix) (including k==0):
#                 - ignore k; select ONLY by full prefix match (orientation-aware)
#                 - emit a warning to stderr once
#           2) if k > len(prefix):
#               - require BOTH:
#                   - 5'-end operation == M with M-length >= k
#                   - full prefix match (orientation-aware)
#
#   -s / --sort            (optional)
#       Name-sort input BAM internally (samtools sort -n via pysam).
#       Useful for --paired so you don't need to pre-sort the BAM.
#
#   -t / --threads         (optional, default: 1)
#       Number of threads for pysam I/O (BGZF compression/decompression).
#
#   --paired               (optional)
#       Turn on paired-end mode:
#           - Group reads by query_name.
#           - Evaluate only read1 (R1) by 5'-end rules.
#           - For each passing R1:
#                 output R1
#                 output any R2 whose (reference_id, reference_start)
#                 matches R1's recorded mate coordinates.
#
################################################################################
# Examples:
#
# 1) SE: exact 3-bp soft-clip with prefix ATG
#
#     python pyselectal.py \
#         -n 3 -m 3 \
#         -p ATG \
#         in.namesort.bam \
#         out.se.5p3S.ATG.bam
#
#   Keeps reads with exact 3 bp 5'-soft-clip and soft-clipped prefix ATG/CAT.
#
# 2) SE: mapped 5'-end, prefix ATG
#
#     python pyselectal.py \
#         -n 0 -m 0 \
#         -p ATG \
#         in.namesort.bam \
#         out.se.no5S.ATG.bam
#
#   Keeps reads with no 5'-soft-clip and mapped prefix ATG (or CAT on reverse).
#
# 3) SE: range 1–3 bp soft-clip, any prefix
#
#     python pyselectal.py \
#         -n 1 -m 3 \
#         in.namesort.bam \
#         out.se.5p1to3S.bam
#
#   Keeps reads with 1–3 bp 5'-soft-clip.
#
# 4) SE: range 2–5 bp soft-clip, G homopolymer
#
#     python pyselectal.py \
#         -n 2 -m 5 \
#         -p G \
#         in.namesort.bam \
#         out.se.5p2to5S.Gpoly.bam
#
#   Keeps reads with 2–5 bp soft-clip that is all G (or C on reverse).
#
# 5) PE: exact 3-bp soft-clip with prefix on R1; output R2 mates
#
#     python pyselectal.py \
#         -n 3 -m 3 \
#         -p ATG \
#         --paired \
#         in.namesort.bam \
#         out.pe.5p3S.ATG.bam
#
#   Keeps:
#       - R1 reads passing soft-clip + prefix check
#       - Corresponding R2 mates (matched by mate coordinates)
#
################################################################################
# Notes:
#
#   - Unmapped reads are always rejected.
#
################################################################################

import argparse
import os
import sys
import tempfile
import pysam

SOFT = 4   # 'S' soft-clip in pysam CIGAR codes
MATCH = 0  # 'M' (alignment match) in pysam CIGAR codes

class HelpfulArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write(f"Error: {message}\nUse -h for help.\n")
        raise SystemExit(2)

def die(message):
    sys.stderr.write(f"Error: {message}\nUse -h for help.\n")
    raise SystemExit(2)

def warn_once(key: str, message: str):
    if key in _WARNED:
        return
    _WARNED.add(key)
    sys.stderr.write(f"Warning: {message}\n")

def revcomp(seq: str) -> str:
    """
    Reverse-complement a DNA sequence (ACGTN, case-insensitive).
    """
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def comp_base(base: str) -> str:
    """
    Complement of a single DNA base (A/C/G/T/N), case-insensitive.
    """
    base = base.upper()
    table = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return table.get(base, "N")


def get_5prime_cigar(aln):
    """
    Get the 5'-end CIGAR operation and length, orientation-aware.

    Forward reads: first CIGAR element.
    Reverse reads: last  CIGAR element.
    """
    cig = aln.cigartuples
    if not cig:
        return None, None
    return (cig[-1] if aln.is_reverse else cig[0])


def has_5prime_mapped_exact(aln, prefix, rc_prefix, k, warn_k_ignored=False):
    """
    Exact match mode (n = m = 0):
        - require that the 5'-end is mapped;
        - 3'-end soft-clips are allowed;
    If --prefix is NOT given:
        - keep reads whose 5'-end CIGAR operation is MATCH (M) and M-length >= k
           - if k==0: keep reads with 5'-end operation == M (i.e., no 5' soft-clip)

    If --prefix IS given:
        1) if k <= len(prefix) (including k==0):
            - ignore k; select ONLY by full prefix match (orientation-aware)
            - emit a warning to stderr once
        2) if k > len(prefix):
            - require BOTH:
                - 5'-end operation == M with M-length >= k
                - full prefix match (orientation-aware)
    """
    if aln.is_unmapped:
        return False

    op_5p, len_5p = get_5prime_cigar(aln)

    # Need a MATCH at 5'-end to enforce a prefix
    if op_5p != MATCH:
        return False

    if len_5p is None:
        return False
    

    # If no prefix, only enforce minimum 5' MATCH length (>=k); if k==0, any 5' M passes.
    if not prefix:
        if k <= 0:
            return True
        return int(len_5p) >= k
    
    # If prefix provided:
    # Case 1: k <= len(prefix) => ignore k, select only by FULL prefix match
    if k <= len(prefix):
        if warn_k_ignored:
            warn_once("k_ignored", f"-k={k} is <= len(--prefix)={len(prefix)} in mapped mode; ignoring -k and selecting only by full prefix.")
        seq = aln.query_sequence
        if not seq:
            return False
        if aln.is_reverse:
            return seq[-len(prefix):].upper() == rc_prefix.upper()
        return seq[:len(prefix)].upper() == prefix.upper()

    # Case 2: k > len(prefix) => require BOTH: 5' MATCH length >= k AND full prefix match
    if len_5p < k:
        return False

    seq = aln.query_sequence
    if not seq:
        return False

    if aln.is_reverse:
        return seq[-len(prefix):].upper() == rc_prefix.upper()
    return seq[:len(prefix)].upper() == prefix.upper()

def has_5prime_softclip_exact(aln, n, prefix, rc_prefix):
    """
    Exact mode: n = m > 0 (exact n-bp 5' soft-clip):
      - require 5'-end is (S, n)
      - if prefix provided (length n):
          forward: first n bases == prefix
          reverse: last  n bases == revcomp(prefix)
    """
    if aln.is_unmapped:
        return False

    op_5p, len_5p = get_5prime_cigar(aln)
    if op_5p != SOFT or len_5p != n:
        return False

    if not prefix:
        return True

    seq = aln.query_sequence
    if not seq:
        return False

    if aln.is_reverse:
        return seq[-n:].upper() == rc_prefix.upper()
    return seq[:n].upper() == prefix.upper()

def has_5prime_softclip_range(aln, n, m, base=None, rc_base=None):
    """
    Range mode: n < m:

      - require a 5' soft-clip of length x where n <= x <= m;
      - if base is None: only length is used;
      - if base is provided (single letter A/C/G/T/N):
          forward: all soft-clipped bases == base
          reverse: all soft-clipped bases == complement(base)
    """
    if aln.is_unmapped:
        return False

    op_5p, len_5p = get_5prime_cigar(aln)

    if op_5p != SOFT:
        return False

    x = int(len_5p)
    if x < n or x > m:
        return False

    if not base:
        return True

    seq = aln.query_sequence
    if not seq:
        return False

    if aln.is_reverse:
        sc = seq[-x:]
        target = rc_base.upper()
    else:
        sc = seq[:x]
        target = base.upper()

    for b in sc:
        if b.upper() != target:
            return False
    return True


def alignment_matches(aln, n, m, prefix, k, warn_k_ignored=False):
    """
    Wrapper that dispatches to the correct mode based on n_min, n_max.

      if n == m > 0, then: exact soft-clip regime with optional sequence.
      if n == m == 0,    then: mapped 5'-end regime with optional prefix.
      if n < m,          then: range soft-clip regime with optional homopolymer base.
    """
    if n == m:
        if n == 0:
            rc_prefix = revcomp(prefix) if prefix else None
            return has_5prime_mapped_exact(aln, prefix, rc_prefix, k, warn_k_ignored=warn_k_ignored)
        rc_prefix = revcomp(prefix) if prefix else None
        return has_5prime_softclip_exact(aln, n, prefix, rc_prefix)

    base = prefix.upper() if prefix else None
    rc_base = comp_base(base) if base else None
    return has_5prime_softclip_range(aln, n, m, base, rc_base)


def process_group_pe(records, n, m, prefix, k, out_bam, warn_k_ignored=False):
    """
    Paired-end mode: records = list of AlignedSegment with the same query_name.

    - Select R1 alignments passing alignment_matches.
    - For each selected R1:
        - collect its (next_reference_id, next_reference_start).
    - Write:
        - all passing R1s;
        - all R2s with (reference_id, reference_start) collected above.
    """
    if not records:
        return

    r1_selected = []
    for r in records:
        if r.is_read1 and alignment_matches(r, n, m, prefix, k, warn_k_ignored=warn_k_ignored):
            r1_selected.append(r)

    if not r1_selected:
        return

    mate_coords = []
    for r in r1_selected:
        if r.next_reference_id is not None and r.next_reference_start is not None:
            if r.next_reference_id >= 0 and r.next_reference_start >= 0:
                mate_coords.append((r.next_reference_id, r.next_reference_start))

    r2_selected = []
    for r in records:
        if r.is_read2 and (r.reference_id, r.reference_start) in mate_coords:
            r2_selected.append(r)

    for r in r1_selected + r2_selected:
        out_bam.write(r)


def process_stream_se(in_bam, out_bam, n, m, prefix, k, warn_k_ignored=False):
    """
    Single-end mode: filter each alignment independently.
    """
    for aln in in_bam.fetch(until_eof=True):
        if alignment_matches(aln, n, m, prefix, k, warn_k_ignored=warn_k_ignored):
            out_bam.write(aln)


def process_stream_pe(in_bam, out_bam, n, m, prefix, k, warn_k_ignored=False):
    """
    Paired-end mode: assumes name-sorted BAM (all alignments with the same
    query_name are contiguous).
    """
    current_qname = None
    group = []

    for aln in in_bam.fetch(until_eof=True):
        qn = aln.query_name
        if current_qname is None:
            current_qname = qn
            group = [aln]
        elif qn == current_qname:
            group.append(aln)
        else:
            process_group_pe(group, n, m, prefix, k, out_bam, warn_k_ignored=warn_k_ignored)
            current_qname = qn
            group = [aln]

    if group:
        process_group_pe(group, n, m, prefix, k, out_bam, warn_k_ignored=warn_k_ignored)


def parse_args(argv):

    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("-h", "--help", action="store_true")
    pre_args, _ = pre.parse_known_args(argv)

    parser = HelpfulArgumentParser(
        add_help=False,
        description="Filter BAM alignments by 5'-end soft-clipping / mapped-prefix pattern."
    )

    # Manual/help
    parser.add_argument("-h", "--help", action="store_true", help="Show the manual and exit.")

    parser.add_argument("in_bam", help="Input BAM path or '-' for stdin.")
    parser.add_argument("out_bam", help="Output BAM path or '-' for stdout.")

    parser.add_argument("-n", "--min-softclip", type=int, required=True,
                        help="Minimum 5' soft-clip length (n).")
    parser.add_argument("-m", "--max-softclip", type=int, required=True,
                        help="Maximum 5' soft-clip length (m).")

    parser.add_argument("-p", "--prefix", type=str, default=None,
                        help="Optional prefix/base (meaning depends on mode; see -h).")

    parser.add_argument("-k", type=int, default=0,
                        help=("Mapped mode only (n=m=0): "
                              "If --prefix not given: require >=k 5' MATCH bases (CIGAR). "
                              "If --prefix given: if k<=len(prefix), k is ignored and full prefix is required; "
                              "if k>len(prefix), require both >=k 5' MATCH bases AND full prefix."))

    parser.add_argument("-s", "--sort", action="store_true",
                        help="Name-sort input BAM internally (samtools sort -n via pysam).")

    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of BGZF threads for pysam I/O.")

    parser.add_argument("--paired", action="store_true",
                        help="Paired-end mode: select read1 and emit matching read2 mates.")

    if pre_args.help or len(argv) == 0:
        sys.stdout.write(MANUAL)
        parser.print_help(sys.stdout)
        raise SystemExit(0)
    
    args = parser.parse_args(argv)

    if args.help:
        sys.stdout.write(MANUAL)
        parser.print_help(sys.stdout)
        raise SystemExit(0)

    # basic checks on n, m
    n = args.min_softclip
    m = args.max_softclip
    prefix = args.prefix
    k = args.k

    if n < 0 or m < 0:
        parser.error("Both n and m must be non-negative integers.")
    if n > m:
        parser.error("Require n <= m (min-softclip <= max-softclip).")

    # -k is only meaningful when n=m=0
    if k != 0 and not (n == 0 and m == 0):
        parser.error("-k is only valid when n = m = 0 (mapped 5' mode).")
    if k < 0:
        parser.error("-k must be >= 0.")

    # mode-specific validation of -p/--prefix
    if n == m:
        if n == 0:
            # mapped mode: prefix can be any length if provided
            pass
        else:
            # exact soft-clip mode: prefix length must equal n, if provided
            if prefix is not None and len(prefix) != n:
                parser.error(
                    f"In exact soft-clip mode (n = m = {n}), --prefix length ({len(prefix)}) must equal n."
                )
    else:
        # range mode: prefix is homopolymer base, if provided
        if prefix is not None:
            if len(prefix) != 1 or prefix.upper() not in "ACGTN":
                parser.error(
                    "In range mode (n < m), --prefix must be a single base A/C/G/T/N if provided."
                )

    if args.threads < 1:
        parser.error("--threads must be >= 1.")

    return args


def name_sort_bam(in_bam_path, threads):
    """
    Name-sort BAM using pysam.sort (-n). Returns path to sorted BAM.
    Caller is responsible for cleanup of the returned temp file.
    """
    # Put temp next to output for large files. Using system temp by default.
    fd, tmp_path = tempfile.mkstemp(prefix="softclip5.namesort.", suffix=".bam")
    os.close(fd)

    try:
        # pysam exposes samtools sort as pysam.sort()
        if threads > 1:
            # samtools sort thread option is -@; passthrough in pysam.sort
            pysam.sort("-n", "-@", str(threads), "-o", tmp_path, in_bam_path, catch_stdout=False)
        else:
            pysam.sort("-n", "-o", tmp_path, in_bam_path, catch_stdout=False)
    except Exception as e:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        die(f"Failed to name-sort BAM: {e}")

    return tmp_path

def spool_stdin_to_temp_bam():
    """
    Copy stdin (BAM stream) to a temp file and return its path.
    """
    fd, tmp_path = tempfile.mkstemp(prefix="softclip5.stdin.", suffix=".bam")
    os.close(fd)
    try:
        with open(tmp_path, "wb") as f:
            # binary copy in chunks
            while True:
                chunk = sys.stdin.buffer.read(1024 * 1024)
                if not chunk:
                    break
                f.write(chunk)
    except Exception as e:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        die(f"Failed to spool stdin to temp BAM: {e}")
    return tmp_path

def open_alignment_files(in_path, out_path, threads):
    """
    Support '-' for stdin/stdout BAM streaming.
    """
    try:
        if in_path == "-":
            if threads > 1:
                in_bam = pysam.AlignmentFile(sys.stdin.buffer, "rb", threads=threads)
            else:
                in_bam = pysam.AlignmentFile(sys.stdin.buffer, "rb")
        else:
            if threads > 1:
                in_bam = pysam.AlignmentFile(in_path, "rb", threads=threads)
            else:
                in_bam = pysam.AlignmentFile(in_path, "rb")

        if out_path == "-":
            if threads > 1:
                out_bam = pysam.AlignmentFile(sys.stdout.buffer, "wb", template=in_bam, threads=threads)
            else:
                out_bam = pysam.AlignmentFile(sys.stdout.buffer, "wb", template=in_bam)
        else:
            if threads > 1:
                out_bam = pysam.AlignmentFile(out_path, "wb", template=in_bam, threads=threads)
            else:
                out_bam = pysam.AlignmentFile(out_path, "wb", template=in_bam)
    except Exception as e:
        die(f"Failed to open BAM(s): {e}")
    return in_bam, out_bam
                                                        
def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_args(argv)

    n = args.min_softclip
    m = args.max_softclip
    prefix = args.prefix
    k = args.k
    threads = args.threads

    # warn once if k will be ignored in mapped mode with prefix
    warn_k_ignored = False
    if (n == 0 and m == 0) and prefix and (k <= len(prefix)):
        warn_k_ignored = True
    
    temp_files = []
    in_path = args.in_bam

    # If input is stdin and we need to sort, spool stdin first.
    if args.sort and in_path == "-":
        # in_path only, remove spooled
        spooled = spool_stdin_to_temp_bam()
        temp_files.append(spooled)
        in_path = spooled
    
    # If paired-end mode, name sorting is required for correct grouping.
    # We only enforce via --sort; otherwise assume the user already provided name-sorted input.
    temp_sorted = None
    if args.sort:
        temp_sorted = name_sort_bam(in_path, threads)
        temp_files.append(temp_sorted)
        in_path = temp_sorted

    in_bam, out_bam = open_alignment_files(in_path, args.out_bam, threads)

    try:
        if args.paired:
            process_stream_pe(in_bam, out_bam, n, m, prefix, k, warn_k_ignored=warn_k_ignored)
        else:
            process_stream_se(in_bam, out_bam, n, m, prefix, k, warn_k_ignored=warn_k_ignored)
    finally:
        try:
            in_bam.close()
        except Exception:
            pass
        try:
            out_bam.close()
        except Exception:
            pass
        # temp_files is needed for removing temporary files
        for p in temp_files:
            try:
                os.remove(p)
            except OSError:
                pass

if __name__ == "__main__":
    main()
