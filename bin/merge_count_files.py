#!/usr/bin/env python3

import argparse
import pandas as pd
import sys
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(
        description='Merge count files from multiple replicates'
    )
    parser.add_argument(
        '--input_files',
        nargs='+',
        required=True,
        help='Input count files to merge'
    )
    parser.add_argument(
        '--sample_ids',
        nargs='+',
        required=True,
        help='Sample IDs corresponding to each input file'
    )
    parser.add_argument(
        '--output',
        required=True,
        help='Output merged count file'
    )
    parser.add_argument(
        '--transcript_id_col',
        default=0,
        type=int,
        help='Column index for transcript IDs (0-based, default: 0)'
    )
    parser.add_argument(
        '--length_col',
        default=1,
        type=int,
        help='Column index for transcript lengths (0-based, default: 1)'
    )
    parser.add_argument(
        '--count_col',
        default=2,
        type=int,
        help='Column index for counts (0-based, default: 2)'
    )
    parser.add_argument(
        '--sep',
        default='\t',
        help='Field separator (default: tab)'
    )
    parser.add_argument(
        '--header',
        action='store_true',
        help='Input files have header row'
    )

    return parser.parse_args()

def read_count_file(file_path, sample_id, args):
    """Read a count file and return a DataFrame with transcript info and counts"""
    try:
        # Read the file
        if args.header:
            df = pd.read_csv(file_path, sep=args.sep, header=0)
            # Get column names by position
            cols = df.columns.tolist()
            transcript_col = cols[args.transcript_id_col]
            length_col = cols[args.length_col]
            count_col = cols[args.count_col]
        else:
            df = pd.read_csv(file_path, sep=args.sep, header=None)
            transcript_col = args.transcript_id_col
            length_col = args.length_col
            count_col = args.count_col

        # Create a clean DataFrame with standardized column names
        result_df = pd.DataFrame({
            'transcript_id': df.iloc[:, args.transcript_id_col],
            'length': df.iloc[:, args.length_col],
            sample_id: df.iloc[:, args.count_col]
        })

        return result_df

    except Exception as e:
        print(f"Error reading file {file_path}: {e}", file=sys.stderr)
        sys.exit(1)

def merge_count_files(input_files, sample_ids, args):
    """Merge multiple count files into a single DataFrame"""

    if len(input_files) != len(sample_ids):
        print("Error: Number of input files must match number of sample IDs", file=sys.stderr)
        sys.exit(1)

    merged_df = None

    for file_path, sample_id in zip(input_files, sample_ids):
        print(f"Processing {file_path} for sample {sample_id}")

        # Read the count file
        df = read_count_file(file_path, sample_id, args)

        if merged_df is None:
            # First file - use as base
            merged_df = df
        else:
            # Merge with existing data
            merged_df = pd.merge(
                merged_df,
                df,
                on=['transcript_id', 'length'],
                how='outer'
            )

    # Fill NaN values with 0 (for transcripts not present in all samples)
    count_columns = [col for col in merged_df.columns if col not in ['transcript_id', 'length']]
    merged_df[count_columns] = merged_df[count_columns].fillna(0)

    # Sort by transcript_id for consistent output
    merged_df = merged_df.sort_values('transcript_id')

    return merged_df

def main():
    args = parse_args()

    print(f"Merging {len(args.input_files)} count files...")

    # Merge the count files
    merged_df = merge_count_files(args.input_files, args.sample_ids, args)

    # Write the merged file
    merged_df.to_csv(args.output, sep='\t', index=False, header = True)

    # Write the header of the merged file
    merged_df.columns.to_series().to_csv('column_names.txt', index=False, header=False)

    print(f"Merged count file written to {args.output}")
    print(f"Final dimensions: {merged_df.shape[0]} transcripts x {merged_df.shape[1]} columns")
    print(f"Sample columns: {[col for col in merged_df.columns if col not in ['transcript_id', 'length']]}")

if __name__ == '__main__':
    main()
