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
import shutil
import numpy as np
from signal import signal, SIGPIPE, SIG_DFL
from STRique_lib import fast5Index




if __name__ == '__main__':
    signal(SIGPIPE, SIG_DFL) 
    # cmd arguments
    parser = argparse.ArgumentParser(description="Mask region in raw nanopore fast5 file. Read stdin as tab-separated ID, ticks, offset")
    parser.add_argument("index", help="Path to input fast5 index")
    parser.add_argument("output", help="Path to output .fast5 directory with masked reads")
    args = parser.parse_args()
    if os.path.exists(args.output):
        raise RuntimeError("[ERROR] Output directory already exists")
    shutil.copytree(os.path.dirname(args.index), args.output)
    f5 = fast5Index.fast5Index(args.index)
    # read stdin as ID  ticks   offset
    for line in sys.stdin:
        ID, offset, ticks = line.strip().split('\t')
        ticks = int(ticks)
        offset = int(offset)
        try:
            record = f5.get_raw(ID)
        except:
            continue
        if not record is None:
            mask = np.ones(record.shape, dtype=bool)
            mask[offset:offset+ticks] = False
            f5.set_raw(ID, record[mask])
