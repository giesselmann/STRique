# \MODULE\-------------------------------------------------------------------------
#
#  CONTENTS      : Mask signal segment in nanopore fast5 file
#
#  DESCRIPTION   :
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# ---------------------------------------------------------------------------------
# Copyright (c) 2018-2019,  Pay Giesselmann, Max Planck Institute for Molecular Genetics
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Written by Pay Giesselmann
# ---------------------------------------------------------------------------------
# public imports
import os, sys
import argparse
import numpy as np
from collections import namedtuple
from signal import signal, SIGPIPE, SIG_DFL
from STRique_lib import fast5Index




if __name__ == '__main__':
    signal(SIGPIPE, SIG_DFL)
    # cmd arguments
    parser = argparse.ArgumentParser(description="Mask region in raw nanopore fast5 file. Read stdin as tab-separated ID, ticks, offset")
    parser.add_argument("index", help="Path to input fast5 index")
    parser.add_argument("counts", help="Path to STRique count output file")
    parser.add_argument("output", help="Path to output .fast5 directory with masked reads")
    parser.add_argument("--format", default='bulk', choices=['single', 'bulk'], help="Output fast5 format")
    args = parser.parse_args()
    STRique_record = namedtuple('STRique_record', ['ID', 'target', 'strand', 'count', 'score_prefix', 'score_suffix', 'log_p', 'offset', 'ticks', 'mod'])
    # read STRique count output into memory
    print("[INFO] Reading repeat counts.", file=sys.stderr)
    with open(args.counts, 'r') as fp:
        row_iter = (row.strip().split('\t') for row in fp if row and not row.startswith('ID'))
        records = [STRique_record(*row[:3], int(row[3]), *[float(x) for x in row[4:7]], int(row[7]), int(row[8]), row[9]) for row in row_iter]
    # write read IDs to batch.txt in output
    print("[INFO] Extract read IDs.", file=sys.stderr)
    os.makedirs(args.output, exist_ok=True)
    f_read_IDs = os.path.join(args.output, 'reads.txt')
    f_read_Index = os.path.join(args.output, 'reads.fofn')
    with open(f_read_IDs, 'w') as fp:
        fp.write('\n'.join([record.ID for record in records]))
        fp.write('\n')
    # load input fast5 index and extract reads from count file
    print("[INFO] Extracting evaluated reads from {input} to {output}".format(
        input=os.path.dirname(args.index), output=args.output), file=sys.stderr)
    f5 = fast5Index.fast5Index(args.index)
    f5.extract(f_read_IDs, args.output, format=args.format)
    # index extracted output
    print("[INFO] Indexing output directory.", file=sys.stderr)
    with open(f_read_Index, 'w') as fp:
        for record in fast5Index.fast5Index.index(args.output):
            fp.write(record)
        fp.write('\n')
    # truncate reads in output directory
    print("[INFO] Masking output reads in-place.", file=sys.stderr)
    f5 = fast5Index.fast5Index(f_read_Index)
    for record in records:
        raw_record = None
        try:
            raw_record = f5.get_raw(record.ID)
        except:
            print("[WARNING] Could not get raw signal for read {}".format(record.ID), file=sys.stderr)
            continue
        if raw_record is not None:
            mask = np.ones(raw_record.shape, dtype=bool)
            mask[record.offset:record.offset+record.ticks] = False
            f5.set_raw(record.ID, raw_record[mask])
