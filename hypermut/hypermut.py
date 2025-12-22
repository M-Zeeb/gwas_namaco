#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 09:55:03 2022

@author: adapted from https://gist.github.com/philiptzou/8d6c7c61d2242f730a1a2f87ba9a2a72; https://doi.org/10.1371/journal.pone.0225352; https://doi.org/10.1093/bioinformatics/16.4.400
"""

# Hypermut


import re
import time
import sys
from scipy.stats import fisher_exact

inf = float('inf')

IUPAC_CODES = {
    'R': r'[RAG]',
    'Y': r'[YCT]',
    'B': r'[BYSKCGT]',
    'D': r'[DRWKAGT]',
    'H': r'[HYWMACT]',
    'V': r'[VRSMAGC]',
    'N': r'[NRYBDHVWSKMACGT]',
    'W': r'[WAT]',
    'S': r'[SCG]',
    'K': r'[KGT]',
    'M': r'[MAC]',
}

VALID_NA_PATTERN = re.compile(r'[NRYBDHVWSKMACGT]')


def expand_iupac(pattern):
    result = []
    for char in pattern:
        if char in IUPAC_CODES:
            result.append(IUPAC_CODES[char])
        else:
            result.append(char)
    return ''.join(result)


DEFAULT_PATTERNS = {
    'hypermut_from': re.compile(expand_iupac(r'^G(?=RD)')),
    'hypermut_to': re.compile(expand_iupac(r'^A(?=RD)')),
    'control_from': re.compile(expand_iupac(r'^G(?=YN|RC)')),
    'control_to': re.compile(expand_iupac(r'^A(?=YN|RC)')),
}


def fasta_reader(filename):
    with open(filename) as fp:
        header = None
        seq = []
        for line in fp:
            if line.startswith('#'):
                continue
            elif line.startswith('>'):
                if seq:
                    yield header, ''.join(seq).upper()
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if seq:
            yield header, ''.join(seq).upper()


def find_sites(seq, pattern, site_range):
    sites = []
    for offset in site_range:
        match = pattern.search(seq[offset:])
        if match:
            sites.append(offset)
    return sites


def get_comparable_sites(refseq, naseq):
    sites = []
    for offset, (ref, na) in enumerate(zip(refseq, naseq)):
        if not VALID_NA_PATTERN.match(ref):
            continue
        if not VALID_NA_PATTERN.match(na):
            continue
        sites.append(offset)
    return sites


def hypermut(refseq, naseq, patterns=DEFAULT_PATTERNS):
    comparable_sites = get_comparable_sites(refseq, naseq)
    potential_muts = find_sites(
        refseq, patterns['hypermut_from'], comparable_sites)
    potential_ctrls = find_sites(
        refseq, patterns['control_from'], comparable_sites)
    matched_muts = find_sites(
        naseq, patterns['hypermut_to'], potential_muts)
    matched_ctrls = find_sites(
        naseq, patterns['control_to'], potential_ctrls)
    num_potential_muts = len(potential_muts)
    num_matched_muts = len(matched_muts)
    num_potential_ctrls = len(potential_ctrls)
    num_matched_ctrls = len(matched_ctrls)
    try:
        oddsratio = (
                (num_matched_muts / num_potential_muts) /
                (num_matched_ctrls / num_potential_ctrls)
        )
    except ZeroDivisionError:
        oddsratio = inf
    _, p = fisher_exact([
        [num_matched_muts, num_potential_muts - num_matched_muts],
        [num_matched_ctrls, num_potential_ctrls - num_matched_ctrls]
    ], 'greater')
    return (
        num_matched_muts,
        num_potential_muts,
        num_matched_ctrls,
        num_potential_ctrls,
        oddsratio,
        p)

def main():
    if len(sys.argv) != 3:
        print('Usage: {} <FASTA_FILE>'.format(sys.argv[0]),
              file=sys.stderr)
        exit(1)
    fasta_filename = sys.argv[1]
    # fasta_filename = "/Users/mariuszeeb/Downloads/bamtrial/6733_ex_reads.fa"
    sequences = list(fasta_reader(fasta_filename))
    ref_filename = sys.argv[2]
    # ref_filename = "/Users/mariuszeeb/Downloads/bamtrial/6733_ex_ref_reads.fa"
    refseqs = list(fasta_reader(ref_filename))

    for refseq, naseq in zip(refseqs, sequences):
        _, refseqss = refseq
        _, naseqss = naseq
        print(hypermut(refseqss, naseqss, DEFAULT_PATTERNS))


if __name__ == '__main__':
    main()
