#!/usr/bin/env python3
"""
Add read groups or haplotype tags to BAM file based on primary alignment location.
All alignments (primary, secondary, supplementary) of a read get the same tag
based on where the primary alignment maps.
"""

import argparse
import pysam
import sys
from pathlib import Path
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description='Add RG or HP tags to all alignments based on primary alignment chromosome',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Add RG tags based on primary alignment chromosome
  %(prog)s -i input.bam -o output.bam -t RG -m chromosome

  # Add HP tags based on primary alignment position bins
  %(prog)s -i input.bam -o output.bam -t HP -m bin -b 1000000

  # Add RG tags with custom prefix and header entries
  %(prog)s -i input.bam -o output.bam -t RG -m chromosome -p haplotype --add-rg-header

Notes:
  - Requires BAM file to be sorted by read name (samtools sort -n)
  - All alignments of a read get the same tag based on primary alignment
  - Unmapped reads and reads without primary alignments are skipped
        """
    )

    parser.add_argument('-i', '--input', required=True,
                        help='Input BAM file (should be name-sorted)')
    parser.add_argument('-o', '--output', required=True,
                        help='Output BAM file')
    parser.add_argument('-t', '--tag', choices=['RG', 'HP'], default='RG',
                        help='Tag type to add (RG=read group, HP=haplotype) (default: RG)')
    parser.add_argument('-m', '--mode', choices=['bin', 'chromosome'], default='chromosome',
                        help='Assignment mode based on primary alignment (default: chromosome)')
    parser.add_argument('-b', '--bin-size', type=int, default=1000000,
                        help='Bin size in bp for position-based grouping (default: 1000000)')
    parser.add_argument('-p', '--prefix', default='group',
                        help='Prefix for read group/haplotype IDs (default: group)')
    parser.add_argument('--add-rg-header', action='store_true',
                        help='Add @RG lines to BAM header (only for RG tags)')
    parser.add_argument('--coordinate-sort', action='store_true',
                        help='Output coordinate-sorted BAM instead of name-sorted')

    return parser.parse_args()


def get_tag_value(read, mode, bin_size, prefix):
    """Generate tag value based on alignment position"""
    if mode == 'bin':
        bin_id = read.reference_start // bin_size
        return f"{prefix}_{read.reference_name}_{bin_id}"
    else:  # chromosome mode
        return f"{prefix}_{read.reference_name}"


def add_rg_to_header(header, rg_ids, prefix):
    """Add @RG lines to BAM header"""
    if 'RG' not in header:
        header['RG'] = []

    for rg_id in rg_ids:
        if not any(rg['ID'] == rg_id for rg in header['RG']):
            header['RG'].append({
                'ID': rg_id,
                'SM': prefix,
                'PL': 'UNKNOWN',
                'LB': prefix
            })

    return header


def process_bam_two_pass(args):
    """
    Two-pass approach:
    Pass 1: Build lookup table of read_name -> tag based on primary alignment
    Pass 2: Apply tags to all alignments
    """

    if not Path(args.input).exists():
        print(f"Error: Input file '{args.input}' not found", file=sys.stderr)
        sys.exit(1)

    # PASS 1: Build lookup table from primary alignments
    print("Pass 1: Scanning primary alignments...", file=sys.stderr)

    read_tags = {}  # read_name -> tag_value
    primary_count = 0

    try:
        infile = pysam.AlignmentFile(args.input, "rb")
    except Exception as e:
        print(f"Error opening input BAM: {e}", file=sys.stderr)
        sys.exit(1)

    for read in infile:
        # Only process primary alignments
        if read.is_secondary or read.is_supplementary:
            continue

        if read.is_unmapped:
            continue

        # Get tag value based on primary alignment position
        tag_value = get_tag_value(read, args.mode, args.bin_size, args.prefix)
        read_tags[read.query_name] = tag_value
        primary_count += 1

        if primary_count % 100000 == 0:
            print(f"  Processed {primary_count:,} primary alignments...", file=sys.stderr)

    infile.close()

    print(f"  Found {primary_count:,} primary alignments", file=sys.stderr)
    print(f"  Unique reads: {len(read_tags):,}", file=sys.stderr)

    # PASS 2: Apply tags to all alignments
    print("\nPass 2: Applying tags to all alignments...", file=sys.stderr)

    infile = pysam.AlignmentFile(args.input, "rb")
    header = infile.header.to_dict()

    # Add RG entries to header if requested
    if args.tag == 'RG' and args.add_rg_header:
        unique_tags = set(read_tags.values())
        header = add_rg_to_header(header, unique_tags, args.prefix)
        print(f"  Added {len(unique_tags)} @RG entries to header", file=sys.stderr)

    # Open output file
    try:
        outfile = pysam.AlignmentFile(args.output, "wb", header=header)
    except Exception as e:
        print(f"Error opening output BAM: {e}", file=sys.stderr)
        infile.close()
        sys.exit(1)

    total_reads = 0
    tagged_reads = 0
    skipped_reads = 0

    for read in infile:
        total_reads += 1

        # Look up tag for this read
        if read.query_name in read_tags:
            tag_value = read_tags[read.query_name]

            # Set tag based on type
            if args.tag == 'HP':
                # HP tags should be integers
                # Extract haplotype number from tag value (e.g., "group_chr4_1" -> 1)
                # Format expected: prefix_chrX_N where N is the haplotype number
                try:
                    parts = tag_value.split('_')
                    # Get the last part and try to convert to int
                    tag_value_int = int(parts[-1])
                except (ValueError, IndexError):
                    # If can't extract integer, use 0
                    tag_value_int = 0
                read.set_tag(args.tag, tag_value_int)
            else:
                read.set_tag(args.tag, tag_value)

            tagged_reads += 1
        else:
            # No primary alignment found for this read
            skipped_reads += 1

        outfile.write(read)

        if total_reads % 100000 == 0:
            print(f"  Processed {total_reads:,} total alignments...", file=sys.stderr)

    infile.close()
    outfile.close()

    # Summary
    print(f"\nComplete!", file=sys.stderr)
    print(f"Total alignments processed: {total_reads:,}", file=sys.stderr)
    print(f"Alignments tagged: {tagged_reads:,}", file=sys.stderr)
    print(f"Alignments skipped (no primary): {skipped_reads:,}", file=sys.stderr)
    print(f"Output written to: {args.output}", file=sys.stderr)

    # Sort and index if requested
    if args.coordinate_sort:
        print(f"\nSorting by coordinate...", file=sys.stderr)
        sorted_output = args.output.replace('.bam', '.sorted.bam')
        pysam.sort("-o", sorted_output, args.output)
        print(f"Sorted output: {sorted_output}", file=sys.stderr)

        print(f"Indexing sorted BAM...", file=sys.stderr)
        pysam.index(sorted_output)
        print(f"Index created: {sorted_output}.bai", file=sys.stderr)
    else:
        print(f"\nNote: Output is in same sort order as input", file=sys.stderr)
        print(f"For JBrowse2, you may need to coordinate-sort and index:", file=sys.stderr)
        print(f"  samtools sort -o {args.output.replace('.bam', '.sorted.bam')} {args.output}", file=sys.stderr)
        print(f"  samtools index {args.output.replace('.bam', '.sorted.bam')}", file=sys.stderr)


def main():
    args = parse_args()
    process_bam_two_pass(args)


if __name__ == '__main__':
    main()
