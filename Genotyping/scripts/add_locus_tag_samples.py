#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import pandas as pd
import hashlib

def add_locus_tag_samples(input_file, output_file):
    """
    Reads a CSV file with sample information and adds a 'LOCUSTAG' column.
    The 'LOCUSTAG' is derived from the 'STRAIN' column by md5sum on the STRAIN ID.
    Writes the modified DataFrame to a new CSV file.
    """
    try:
        md5sum = lambda x: hashlib.md5(x.encode()).hexdigest()
        df = pd.read_csv(input_file)
        seentags = dict()
        # might want to fix this to only 
        # update if value is 
        print(df.columns)
        if 'LocusTag' not in df.columns:
            df['LocusTag'] = None
        for row in df.itertuples():
            if pd.isna(row.LocusTag) or row.LocusTag == "":
                df.at[row.Index, 'LocusTag'] = "C" + md5sum(row.Strain).upper()[-6:]
            tag = df.at[row.Index, 'LocusTag']
            if tag not in seentags:
                seentags[tag] = 1
            else:
                seentags[tag] += 1
                print(f"Duplicate LocusTag found: {tag} for {row.Strain}, incrementing count to {seentags[tag]}")
                # If we have a duplicate, we could append a number or something to make it unique
                df.at[row.Index, 'LocusTag'] = f"{tag}{seentags[tag]}"
    
        
        df.to_csv(output_file, index=False)
        print(f"Successfully added LOCUSTAG to {output_file}")
    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python add_locus_tag_samples.py <input_file> <output_file>", file=sys.stderr)
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    add_locus_tag_samples(input_file, output_file)
        