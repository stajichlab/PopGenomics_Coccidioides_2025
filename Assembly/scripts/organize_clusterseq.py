#!/usr/bin/env python3

import sys
import re
import os
import argparse
import csv
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

def parse_arguments():
    parser = argparse.ArgumentParser(description="Gather region protein files by cluster-type.")
    parser.add_argument("input_dir", help="Directory containing the AntiSMASH protein processed fasta files.")
    parser.add_argument("output_dir", help="Directory to save the renamed region files.")
    parser.add_argument("-a", "--annotations",
                        default="record_annotations.tsv",
                        help="annotations file renamed region files.")
    parser.add_argument("-t","--threads", default="2", help="threads for processing datafiles in parallel")
    
    return parser.parse_args()

def process_classtype(indata):
    """
    Process all records for a given classtype and write to output file.
    """
    import sys
    import os
    import re
    (classtype, records, inputdir, outdir) = indata

    print(f"There are {len(records)} records for {classtype}", file=sys.stderr)
    outfile = os.path.join(outdir, f"{classtype}.clusters.fas")
    n = 0
    with open(outfile, "w") as outfh:
        for record in records:
            m = re.match(r'^(\S+)_R\.bin', record)
            MAG = ""
            if m:
                MAG = m.group(1)
            else:
                m = re.match(r'^(\S+)\.bin', record)
                if m:
                    MAG = m.group(1)
                else:
                    print(f"cannot determine MAG parent from {record}", file=sys.stderr)
                    continue
            fafile = os.path.join(inputdir, MAG, f"{record}.fasta")
            if not os.path.exists(fafile):
                print(f"cannot find file {fafile} in type {classtype}", file=sys.stderr)
                continue
            with open(fafile, "rt") as fain:
                for line in fain:
                    if line.startswith(">"):
                        n += 1
                    outfh.write(line)
    print(f"Wrote {n} sequences to Class {classtype} {outfile}", file=sys.stderr)
    
def gather_cluster_files(inputdir, annotations, outdir, num_jobs):
    """
    Process region files fasta and clusters.
    """
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    annot_table = dict()
    with open(annotations, "rt") as infh:
        recordparser = csv.DictReader(infh,delimiter="\t")
        for row in recordparser:
            recordname = row["Record"]
            gbk = row["GBK"]
            clusterclass = row["Class"]
            category = row["Category"]
            if category not in annot_table:
                annot_table[category] = set()
            annot_table[category].add(gbk)

    process_args = []
    for classtype, records in annot_table.items():
        process_args.append([classtype, records, inputdir, outdir])

    print(f"Found {len(process_args)} classtypes to process", file=sys.stderr)
    # Process classtypes in parallel
    if num_jobs is None:
        num_jobs = min(len(process_args), multiprocessing.cpu_count())
    
    print(f"Processing {len(process_args)} classtypes using {num_jobs} processes", file=sys.stderr)
    
    with ProcessPoolExecutor(max_workers=num_jobs) as executor:
        results = list(executor.map(process_classtype, process_args))
    
    # Print summary
    total_sequences = sum(result[1] for result in results)
    print(f"Total sequences processed: {total_sequences}", file=sys.stderr)

if __name__ == "__main__":
    args = parse_arguments()
    gather_cluster_files(args.input_dir, args.annotations, args.output_dir, int(args.threads))
