#!/usr/bin/env python

''' Script to append to sample table with id, single_end and bam / bigwig files for CAGEr '''

import argparse
from pathlib import Path

parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
parser.add_argument(
    "data_in",
    metavar="data_in",
    type=lambda x: dict(i.split(':') for i in x.split(',')),
    help="Input in Dict[str, str] format",
)
parser.add_argument(
    "file_out",
    metavar="FILE_OUT",
    type=Path,
    help="Transformed output samplesheet in CSV format.",
)

args = parser.parse_args()

with open(args.file_out, encoding="utf-8", mode="w+") as outfile:
    line = f"{args.data_in}"
    outfile.write(line)
