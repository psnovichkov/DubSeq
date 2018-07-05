import os
import sys
import argparse

__ACTG = ['A', 'C', 'G', 'T']
__RC_ACTG = ['T', 'G', 'C', 'A']
__NUCL_INDEX = [0] * 1000
__NUCL_INDEX[ord('A')] = 0
__NUCL_INDEX[ord('C')] = 1
__NUCL_INDEX[ord('G')] = 2
__NUCL_INDEX[ord('T')] = 3


if sys.version_info[0] < 3:
    def print_prefix(x): print(x,)
else:
    def print_prefix(x): print(x, end='', flush=True)


def reverse_complement(sequence):
    rc_sequence = ''
    for ch in sequence:
        rc_sequence = __RC_ACTG[__NUCL_INDEX[ord(ch)]] + rc_sequence
    return rc_sequence


def to_hex_code(sequence, rev_compl=False):
    if rev_compl:
        sequence = reverse_complement(sequence)

    # we will chop it at the end
    code = 15
    for ch in sequence:
        code <<= 2
        code += __NUCL_INDEX[ord(ch)]
    return hex(code)[3:]


def process_similar_sequences(chars, processor, *args, **kwargs):
    for i, ch in enumerate(chars):
        for _ch in __ACTG:
            if ch != _ch:
                chars[i] = _ch
                processor(chars, *args, **kwargs)
        chars[i] = ch


def process_fastq_files(source, processor, *args, **kwargs):
    process_files(source, processor, ('.fastq', '.fastq.gz'), *args, **kwargs)


def process_files(source, processor, extensions, *args, **kwargs):
    ''' Iterates over all files defined by source and invokes callback
    '''

    if os.path.isfile(source):
        # the srouce is a file
        file_name = source
        processor(file_name, *args, **kwargs)

    elif os.path.isdir(source):
        # the source is a directory
        for root, _, files in os.walk(source):
            for file in files:
                if file.endswith(extensions):
                    file_name = os.path.join(root, file)
                    processor(file_name, *args, **kwargs)


class RawDescriptionArgumentDefaultsHelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                                                  argparse.RawDescriptionHelpFormatter):
    pass
