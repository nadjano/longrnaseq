#!/usr/bin/env python3
"""
Script to fix GTF files by propagating gene_name and oId from transcripts to exons,
standardizing gene_id and transcript_id fields, and making transcript IDs unique per gene.

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
    from their corresponding transcript features, and making transcript IDs unique.
    
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
        # Fill gene_name with gene_id where gene_name is empty/NaN
        ann_df['gene_name'] = ann_df['gene_name'].fillna(ann_df['gene_id'])
        ann_df['gene_id'] = ann_df['gene_name']
    elif 'gene_id' in ann_df.columns:
        # If gene_name doesn't exist, create it from gene_id
        ann_df['gene_name'] = ann_df['gene_id']

    if 'oId' in ann_df.columns:
        ann_df['transcript_id'] = ann_df['oId']
    
    # Make transcript IDs unique per gene
    print("Making transcript IDs unique ...")
    if 'gene_id' in ann_df.columns and 'transcript_id' in ann_df.columns:
        # Create a unique combination identifier
        ann_df['gene_tx_combo'] = ann_df['gene_id'] + '::' + ann_df['transcript_id']
        
        # Assign a unique numeric ID to each gene+transcript combination
        ann_df['unique_id'] = pd.factorize(ann_df['gene_tx_combo'])[0] + 1
        
        # Check which transcript IDs appear in multiple genes (need suffix)
        tx_gene_count = ann_df.groupby('transcript_id')['gene_id'].transform('nunique')
        needs_suffix = tx_gene_count > 1
        
        # Add suffix to make transcript_id globally unique
        ann_df.loc[needs_suffix, 'transcript_id'] = (
            ann_df.loc[needs_suffix, 'transcript_id'] + '_' + 
            ann_df.loc[needs_suffix, 'unique_id'].astype(str)
        )
        
        # Clean up temporary columns
        ann_df = ann_df.drop(columns=['gene_tx_combo', 'unique_id'])
        
        print(f"Modified {needs_suffix.sum()} features with duplicate transcript IDs")
    else:
        print("Warning: Could not make transcript IDs unique - missing gene_id or transcript_id columns")

    
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
        description='Fix GTF file by propagating gene_name and oId from transcripts to exons and making transcript IDs unique'
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