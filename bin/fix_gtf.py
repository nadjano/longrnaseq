#!/usr/bin/env python3
"""
Script to fix GTF files by propagating gene_name and oId from transcripts to exons,
and standardizing gene_id and transcript_id fields.

Usage:
    python fix_gtf.py input.gtf output.gtf
"""

import sys
import argparse
import pyranges as pr
import pandas as pd


def fix_gtf(input_gtf, output_gtf):
    """
    Fix GTF file by filling missing gene_name and oId in exon features
    from their corresponding transcript features.
    
    Parameters:
    -----------
    input_gtf : str
        Path to input GTF file
    output_gtf : str
        Path to output GTF file
    """
    
    print(f"Reading GTF file: {input_gtf}")
    ann = pr.read_gtf(input_gtf, ignore_bad=True)
    
    # Convert to dataframe
    ann_df = ann.df
    print(f"Total features: {len(ann_df)}")
    
    # Populate missing oid and gene_name of exons with the transcript values
    print("Filling missing gene_name and oId for exon features...")
    for col in ['gene_name', 'oId']:
        if col in ann_df.columns:
            ann_df[col] = ann_df.groupby('transcript_id')[col].transform(
                lambda x: x.ffill().bfill()
            )
        else:
            print(f"Warning: Column '{col}' not found in GTF file")
    
    # Standardize gene_id and transcript_id
    print("Standardizing gene_id and transcript_id fields...")
    if 'gene_name' in ann_df.columns:
        ann_df['gene_id'] = ann_df['gene_name']
    if 'oId' in ann_df.columns:
        ann_df['transcript_id'] = ann_df['oId']
    
    # Select relevant columns
    columns_to_keep = [
        'Chromosome', 'Source', 'Feature', 'Start', 'End', 
        'Score', 'Strand', 'Frame', 'transcript_id', 'gene_id', 'gene_name'
    ]
    
    # Only keep columns that exist in the dataframe
    columns_to_keep = [col for col in columns_to_keep if col in ann_df.columns]
    ann_df_filter = ann_df[columns_to_keep]
    
    # Convert back to PyRanges and save
    print(f"Writing fixed GTF file: {output_gtf}")
    ann_fixed = pr.PyRanges(ann_df_filter)
    ann_fixed.to_gtf(output_gtf)
    
    print("Done!")


def main():
    parser = argparse.ArgumentParser(
        description='Fix GTF file by propagating gene_name and oId from transcripts to exons'
    )
    parser.add_argument(
        'input_gtf',
        help='Input GTF file path'
    )
    parser.add_argument(
        'output_gtf',
        help='Output GTF file path'
    )
    
    args = parser.parse_args()
    
    try:
        fix_gtf(args.input_gtf, args.output_gtf)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()