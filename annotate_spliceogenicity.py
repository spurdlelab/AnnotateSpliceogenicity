#!/usr/bin/env python
# coding: utf-8
###############################################################################
#
#    annotate_splicegenecity.py
#
#    Adds spliceogenicity predictions using the ENIGMA thresholds
#
#    Copyright (C) 2019 QIMR Berghofer Medical Research Institute
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

import argparse
import csv
import gzip
import logging
import re
import sys


maxentscan_fieldnames = [
    'MaxEntScan_ref',
    'MaxEntScan_alt',
    'MaxEntScan_diff',
    'MES-SWA_donor_ref',
    'MES-SWA_donor_alt',
    'MES-SWA_donor_ref_comp',
    'MES-SWA_donor_diff',
    'MES-SWA_acceptor_ref',
    'MES-SWA_acceptor_alt',
    'MES-SWA_acceptor_ref_comp',
    'MES-SWA_acceptor_diff',
    'MES-NCSS_upstream_acceptor',
    'MES-NCSS_upstream_donor',
    'MES-NCSS_downstream_acceptor',
    'MES-NCSS_downstream_donor',
]

hgvs_intron_offset_fieldnames = [
    'HGVS_IntronStartOffset',
    'HGVS_IntronEndOffset',
]

additional_output_fields = [
    'Splicing_var_type',
    'Native_loss',
    'Donor_gain',
]


def get_option_parser():

    parser = argparse.ArgumentParser()

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        'input_file',
        nargs='?',
        help='Read tab-delimited values from FILE [stdin]',
        metavar='<in.tsv>',
        type=argparse.FileType('r'),
        default=sys.stdin)

    group.add_argument(
        'output_file',
        nargs='?',
        help='Write tab-delimited values to FILE [stdout]',
        metavar='<out.tsv>',
        type=argparse.FileType('w'),
        default=sys.stdout)

    return parser


def get_data(tsv):
    for row in tsv:
        for key in row.keys():
            if row[key] == '-':
                row[key] = None
            elif key in maxentscan_fieldnames:
                row[key] = float(row[key])
            elif key in hgvs_intron_offset_fieldnames:
                row[key] = int(row[key])
        yield row


def get_splicing_var_type(row, missing=None):
    """
    Annotate splice region according Burge et al
    """
    start_offset = row['HGVS_IntronStartOffset']
    end_offset = row['HGVS_IntronEndOffset']
    ncss_upstream_acceptor = row['MES-NCSS_upstream_acceptor']
    ncss_downstream_donor = row['MES-NCSS_downstream_donor']
    swa_donor_ref = row['MES-SWA_donor_ref']
    swa_acceptor_ref = row['MES-SWA_acceptor_ref']
    if start_offset is None or end_offset is None:
        return missing
    if start_offset == 0 and end_offset == 0 and re.search('splice', row['Consequence']):
        if ncss_downstream_donor is None or ncss_upstream_acceptor is None:
            return 'last_exon'
        else:
            if swa_donor_ref == ncss_downstream_donor:
                return 'Exonic_donor_splice_region'
            if swa_acceptor_ref == ncss_upstream_acceptor:
                return 'Exonic_acceptor_splice_region'
            return 'check'
    if 0 < start_offset <= 6:
        return 'Intronic_donor_splice_region'
    if -20 <= start_offset < 0:
        return 'Intronic_acceptor_splice_region'
    if start_offset == 0 and end_offset > 0:
        return 'Intronic_donor_splice_region'
    if start_offset < 0 and end_offset == 0:
        return 'Intronic_acceptor_splice_region'
    return 'Outside_native'


def get_native_loss(row):
    if row['Splicing_var_type'] in ['Intronic_donor_splice_region', 'Exonic_donor_splice_region']:
        if row['VARIANT_CLASS'] == 'SNV':
            return annotate_native_loss(row['MaxEntScan_alt'], row['MaxEntScan_diff'])
        else:
            return annotate_native_loss(row['MES-SWA_donor_alt'], row['MES-SWA_donor_diff'])
    if row['Splicing_var_type'] in ['Exonic_acceptor_splice_region', 'Intronic_acceptor_splice_region']:
        if row['VARIANT_CLASS'] == 'SNV':
            return annotate_native_loss(row['MaxEntScan_alt'], row['MaxEntScan_diff'])
        else:
            return annotate_native_loss(row['MES-SWA_acceptor_alt'], row['MES-SWA_acceptor_diff'])
    return None


def annotate_native_loss(alt, diff):
    if alt is None or diff is None:
        return None
    if diff <= 0:
        return 'IMPROVED'
    if alt < 6.2:
        if diff >= 1.15:
            return 'HIGH'
        else:
            return 'MODERATE*'
    if alt > 8.5:
        return 'LOW'
    if diff >= 1.15:
        return 'MODERATE'
    return 'LOW'


def get_donor_gain(row):

    if row['Splicing_var_type'] == 'Intronic_donor_splice_region':
        if row['VARIANT_CLASS'] != 'SNV':
            return annotate_donor_gain(row['MES-SWA_donor_alt'])
        if row['MES-SWA_donor_ref_comp'] is None or row['MES-NCSS_upstream_donor'] is None:
            return None
        if row['MES-SWA_donor_ref_comp'] == row['MES-NCSS_upstream_donor']:
            return 'LOW'
        return annotate_donor_gain(row['MES-SWA_donor_alt'])

    if row['Splicing_var_type'] == 'Exonic_donor_splice_region':
        if row['VARIANT_CLASS'] != 'SNV':
            return annotate_donor_gain(row['MES-SWA_donor_alt'])
        if row['MES-SWA_donor_ref_comp'] is None or row['MES-NCSS_downstream_donor'] is None:
            return None
        if row['MES-SWA_donor_ref_comp'] == row['MES-NCSS_downstream_donor']:
            return 'LOW'
        return annotate_donor_gain(row['MES-SWA_donor_alt'])

    if row['Splicing_var_type'] == 'Exonic_acceptor_splice_region':
        if row['VARIANT_CLASS'] != 'SNV':
            return annotate_donor_gain(row['MES-SWA_donor_alt'])
        if row['MES-SWA_donor_ref_comp'] is None or row['MES-SWA_donor_diff'] is None:
            return None
        if row['MES-NCSS_upstream_donor'] is None or row['MES-NCSS_downstream_donor'] is None:
            return None
        if row['MES-SWA_donor_diff'] >= 0:
            return 'LOW'
        donor_gain = annotate_donor_gain(row['MES-SWA_donor_alt'])
        if donor_gain == 'MODERATE':
            if row['MES-SWA_donor_alt'] > row['MES-NCSS_downstream_donor'] and row['MES-SWA_donor_alt'] > row['MES-NCSS_upstream_donor']:
                return donor_gain + '[both]'
            if row['MES-SWA_donor_alt'] > row['MES-NCSS_downstream_donor']:
                return donor_gain + '[downstream]'
            if row['MES-SWA_donor_alt'] > row['MES-NCSS_upstream_donor']:
                return donor_gain + '[upstream]'
        return donor_gain

    if row['Splicing_var_type'] == 'Outside_native':
        if row['MES-SWA_donor_ref_comp'] is None or row['MES-SWA_donor_diff'] is None:
            return None
        if row['MES-NCSS_upstream_donor'] is None or row['MES-NCSS_downstream_donor'] is None:
            return None
        if row['MES-SWA_donor_diff'] >= 0:
            return 'LOW'
        donor_gain = annotate_donor_gain(row['MES-SWA_donor_alt'])
        if donor_gain == 'MODERATE' and re.search('intron', row['Consequence']):
            if row['MES-SWA_donor_alt'] > row['MES-NCSS_upstream_donor']:
                return donor_gain + '[intronic]'
            return 'LOW'
        return donor_gain

    return None


def annotate_donor_gain(alt):
    if alt is None:
        return None
    if alt >= 8.5:
        return 'HIGH'
    if alt >= 6.2:
        return 'MODERATE'
    return 'LOW'


def main():

    logger = logging.getLogger()
    handler = logging.StreamHandler(sys.stderr)
    logger.addHandler(handler)

    parser = get_option_parser()
    args = parser.parse_args()

    if sys.stdin.isatty() and args.input_file is sys.stdin:
        logger.error("Input could not be read from STDIN.")
        return

    if not sys.stdin.isatty() and args.input_file is not sys.stdin:
        logger.error("Cannot read from both STDIN and a FILE.")
        return

    if args.input_file.name.endswith('.gz'):
        input_file = gzip.open(args.input_file.name, 'rt')
    else:
        input_file = args.input_file

    if args.output_file.name.endswith('.gz'):
        output_file = gzip.open(args.output_file.name, 'wt')
    else:
        output_file = args.output_file

    with input_file as infile, output_file as outfile:
        lines = (row.lstrip('#') for row in infile if not row.startswith('##'))
        reader = csv.DictReader(lines, delimiter='\t')
        fieldnames = reader.fieldnames + additional_output_fields
        writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t', lineterminator='\n')
        writer.writeheader()

        for row in get_data(reader):

            row['Splicing_var_type'] = get_splicing_var_type(row)
            row['Native_loss'] = get_native_loss(row)
            row['Donor_gain'] = get_donor_gain(row)

            writer.writerow(row)


if __name__ == '__main__':
    main()
