# \MODULE\-------------------------------------------------------------------------
#
#  CONTENTS      : Basic fast5 file reader
#
#  DESCRIPTION   : 
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# ---------------------------------------------------------------------------------
# Copyright (c) 2018,  Pay Giesselmann, Max Planck Institute for Molecular Genetics
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
import string
import h5py




# Constants
# fast5 file data locations
LOC_ANALYSES = "/Analyses/"
LOC_RAW = "/Raw/"
LOC_CHANNEL_ID = "/UniqueGlobalKey/channel_id/"
NAME_BASECALL = "Basecall_1D"
NAME_BASECALL_TEMPLATE = "BaseCalled_template"
NAME_EVENT = "EventDetection"




# data types
class fast5Record():
    def __init__(self):
        self.name = ""
        self.ID = ""
        self.sequence = ""
        self.quality = ""
        self.raw = None
        self.rawEvents = None
        self.bcEvents = None
        self.attributes = {}




# Extract currents and quality information from ONT fast5 file
class fast5Reader(object):
    ## Methods
    # constructor
    def __init__(self):
        pass

    def __grpSelect__(self, name, obj):
        if isinstance(obj, h5py.Group):
            if "Events" in obj and "Fastq" in obj:
                return name

    def getRecordID(self, path2fast5):
        with h5py.File(path2fast5, 'r') as f5:
            s = f5[LOC_RAW].visit(lambda name: name if 'Signal' in name else None)
            return str(f5[LOC_RAW + '/' + s.rpartition('/')[0]].attrs['read_id'], 'utf-8')

    def getRecordName(self, path2fast5):
        with h5py.File(path2fast5, 'r') as f5:
            grp = f5[LOC_ANALYSES]
            s = f5[LOC_ANALYSES].visititems(self.__grpSelect__)
            if s is not None:
                fastq = f5[LOC_ANALYSES + '/' + s]["Fastq"].value.decode("utf-8").split('\n')
                return fastq[0][1:]

    # extract read data from fast5 file
    def getRecord(self, path2fast5):
        record = fast5Record()
        with h5py.File(path2fast5, 'r') as f5:          
            # get read and base calling events if available
            if LOC_ANALYSES in f5:
                grp = f5[LOC_ANALYSES]
                s = grp.visititems(self.__grpSelect__)
                if s is not None:
                    sub = grp[s]
                    fastq = sub["Fastq"].value.decode("utf-8").split('\n')
                    record.name = fastq[0][1:]
                    record.sequence = fastq[1]
                    record.quality = fastq[3]
                    record.bcEvents = sub["Events"].value
                for w in grp.keys():
                    if NAME_EVENT in w:
                        s = grp[w].visit(lambda name: name if 'Events' in name else None)
                        if s is not None:
                            record.rawEvents = grp[w + '/' + s].value
                            record.attributes['eventsStart'] = grp[w + '/' + s].parent.attrs['start_time']
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
