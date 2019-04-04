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
import os, sys, h5py
import argparse
import shutil
import numpy as np
import matplotlib.pyplot as plt
from signal import signal, SIGPIPE, SIG_DFL




# constants
LOC_RAW = "/Raw/"
LOC_CHANNEL_ID = "/UniqueGlobalKey/channel_id/"




# index fast5 files directory
class fast5Index():
    # data types
    class fast5Record():
        def __init__(self):
            self.ID = ""
            self.raw = np.array([])
            self.attributes = {}

    # init
    def __init__(self, path2fast5, recursive=True):
        self.path2fast5 = path2fast5
        self.__records = {}
        if not recursive:
            f5files = glob.glob(os.path.join(path2fast5, '*.fast5'))
        else:
            f5files = [os.path.join(dirpath, f) for dirpath, _, files in os.walk(path2fast5) for f in files if f.endswith('.fast5')]
        for f in f5files:
            try:
                self.__records[self.__getRecordID__(f)] = os.path.relpath(f, start=path2fast5).replace("\\", "/")
            except:
                sys.stderr.write("Error reading " + f + '\n')
                
    # get unique ID of fast5 read
    def __getRecordID__(self, path2fast5):
        with h5py.File(path2fast5, 'r') as f5:
            s = f5[LOC_RAW].visit(lambda name: name if 'Signal' in name else None)
            return str(f5[LOC_RAW + '/' + s.rpartition('/')[0]].attrs['read_id'], 'utf-8')

    # extract read data from fast5 file
    def get_record(self, name):
        record = None
        if name in self.__records:
            record = self.fast5Record()
            with h5py.File(os.path.join(self.path2fast5, self.__records[name]), 'r') as f5:
                # get raw current signal
                if LOC_RAW in f5:
                    s = f5[LOC_RAW].visit(lambda name: name if 'Signal' in name else None)
                    if s is not None:
                        record.raw = f5[LOC_RAW + '/' + s].value
                        record.ID = str(f5[LOC_RAW + '/' + s.rpartition('/')[0]].attrs['read_id'], 'utf-8')
                        record.attributes['rawStart'] = f5[LOC_RAW + '/' + s].parent.attrs['start_time']
                        record.attributes['readNumber'] = f5[LOC_RAW + '/' + s].parent.attrs['read_number']
                        try:
                            record.attributes['samplingRate'] = f5[LOC_CHANNEL_ID].attrs['sampling_rate']
                            record.attributes['digitisation'] = f5[LOC_CHANNEL_ID].attrs['digitisation']
                            record.attributes['offset'] = f5[LOC_CHANNEL_ID].attrs['offset']
                            record.attributes['channel'] = int(f5[LOC_CHANNEL_ID].attrs['channel_number'])
                        except:
                            pass
        return record
        
    def set_record_signal(self, name, raw_signal):
        if name in self.__records:
            with h5py.File(os.path.join(self.path2fast5, self.__records[name]), 'r+') as f5:
                # get raw current signal
                if LOC_RAW in f5:
                    s = f5[LOC_RAW].visit(lambda name: name if 'Signal' in name else None)
                    if s is not None:
                        del f5[LOC_RAW + '/' + s]
                        ds = f5.create_dataset(LOC_RAW + '/' + s, data=raw_signal)
                        ds.parent.attrs['duration'] = len(raw_signal)




if __name__ == '__main__':
    signal(SIGPIPE, SIG_DFL) 
    # cmd arguments
    parser = argparse.ArgumentParser(description="Mask region in raw nanopore fast5 file. Read stdin as tab-separated ID, ticks, offset")
    parser.add_argument("--f5_in", required=True, help="Path to input .fast5 directory")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--f5_out", default='', help="Path to output .fast5 directory with masked reads")
    group.add_argument("--inplace", action='store_true', help="Mask signal in place. This is NOT reversible")
    args = parser.parse_args()
    if not args.inplace and args.f5_out != '':
        shutil.copytree(args.f5_in, args.f5_out)
        f5 = fast5Index(args.f5_out, recursive=True)
    elif args.inplace:
        f5 = fast5Index(args.f5_in, recursive=True)
    # read stdin as ID  ticks   offset
    for line in sys.stdin:
        ID, offset, ticks = line.strip().split('\t')
        ticks = int(ticks)
        offset = int(offset)
        record = f5.get_record(ID)
        if record:
            mask = np.ones(record.raw.shape, dtype=bool)
            mask[offset:offset+ticks] = False
            f5.set_record_signal(ID, record.raw[mask])
