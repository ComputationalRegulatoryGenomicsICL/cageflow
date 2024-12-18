#!/usr/bin/env python3
# Given an input file listing all uniquely (or multi-) mapped bigwig outputs
# creates a sample list csv file in the following format
# input:
# path/to/S10.Signal.UniqueMultiple.str1.out.wig.bw
# path/to/S10.Signal.UniqueMultiple.str2.out.wig.bw
# path/to/S14.Signal.UniqueMultiple.str1.out.wig.bw
# path/to/S14.Signal.UniqueMultiple.str2.out.wig.bw
# output:
# id,single_end,path
# S10,false,[path/to/S10.Signal.UniqueMultiple.str1.out.wig.bw path/to/S10.Signal.UniqueMultiple.str2.out.wig.bw]
# S14,false,[path/to/S14.Signal.UniqueMultiple.str1.out.wig.bw path/to/S14.Signal.UniqueMultiple.str2.out.wig.bw]
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f','--filepath', type=str, help='Path to the file with bigwigs')
args = parser.parse_args()

outdict = {}
with open(args.filepath, "r", encoding="utf-8") as filein:
    for line in filein:
        samplepath = line.strip()
        sample_name = samplepath.split("/")[-1].split(".Signal.")[0]
        if sample_name not in outdict:
            outdict[sample_name] = [samplepath]
        else:
            outdict[sample_name].append(samplepath)

with open("sample_list.csv", "w+", encoding="utf-8") as outfile:
    for sample, paths in outdict.items():
        if len(paths) > 1:
            singleend = "false"
            path_str = " ".join(paths)
        else:
            singleend = "true"
            path_str = paths
        line_to_write = f"{sample},{singleend},[{path_str}]\n"
        outfile.write(line_to_write)
