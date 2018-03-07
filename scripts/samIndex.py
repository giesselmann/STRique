# \MODULE\---------------------------------------------------------------
#
#  CONTENTS      : Basic sam file parser
#
#  DESCRIPTION   : 
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# -----------------------------------------------------------------------
#  All rights reserved to Max Planck Institute for Molecular Genetics
#  Berlin, Germany
#  Written by Pay Giesselmann
# -----------------------------------------------------------------------
import os, json
import string, re
import collections


# constants
class samCols():
    QNAME = 0
    FLAG = 1
    RNAME = 2
    POS = 3
    MAPQ = 4
    CIGAR = 5
    RNEXT = 6
    PNEXT = 7
    TLEN = 8
    SEQ = 9
    QUAL = 10


# data types
class samRecord():
    def __init__(self):
        self.name = "*"
        self.flags = 0
        self.reference = "*"
        self.pos = 0
        self.cigar = "*"
        self.seq = "*"
        self.mapq = 0


# decode cigar into list of edits
def decodeCigar(cigar):
    ops = [(int(op[:-1]), op[-1]) for op in re.findall('(\d*\D)',cigar)]
    return ops


# return length of recognized operations in decoded cigar
def opsLength(ops, recOps='MIS=X'):
    n = [op[0] for op in ops if op[1] in recOps]
    return sum(n)


# read sequence alignments from sam file
class samIndex(object):
    # constructor
    def __init__(self):
        self.filename = ""
        self.records = {}
        self.recordNameType = collections.namedtuple('meta', ['name', 'flag', 'reference'])

    # destructor
    def __del__(self):
        pass

    def __str__(self):
        return self.filename

    # read lines, store offsets for reads
    def analyze(self, path2sam):
        self.filename = path2sam
        self.records = {}
        offset = 0
        with open(self.filename, 'r', encoding="ascii", errors="surrogateescape")as sam:
            size = os.path.getsize(self.filename)
            for line in sam:
                if line[0] != '@':      # header line
                    cols = line.split()
                    name = self.__simpleName__(cols[samCols.QNAME])
                    flags = int(cols[samCols.FLAG])
                    reference = cols[samCols.RNAME]
                    if (flags & 0x900) == 0:
                        self.records[name] = (offset, flags, reference)
                l = len(line)
                offset += l

    # load dictionary from disk
    def load(self, path2sam):
        self.filename = path2sam
        with open(path2sam + '.json', 'r') as fp:
            self.records = json.load(fp)

    # load index from disk, create if not existent
    def loadOrAnalyze(self, path2sam):
        if os.path.isfile(path2sam + '.json'):
            self.load(path2sam)
        else:
            self.analyze(path2sam)
            self.save(path2sam)

    # save dictionary to disk
    def save(self, path2sam):
        with open(path2sam + '.json', 'w') as fp:
            json.dump(self.records, fp)

    # return list of (name, flag, reference) in current dataset
    def getRecordNames(self):
        return list(self.records.keys())

    def getRecordMeta(self):
        return [self.recordNameType(x[0], x[1][1], x[1][2]) for x in self.records.items()]

    # return specific dataset if exists in current file
    def getRecord(self, dsName):
        if dsName in self.records.keys():
            with open(self.filename, 'r', encoding="ascii", errors="surrogateescape") as sam:
                offset, flags, reference = self.records[dsName]
                sam.seek(offset)
                line = sam.readline()
                cols = line.split()
                dataset = samRecord()
                dataset.name = self.__simpleName__(cols[samCols.QNAME])
                dataset.flags = int(cols[samCols.FLAG])
                dataset.reference = cols[samCols.RNAME]
                dataset.pos = int(cols[samCols.POS])
                dataset.cigar = cols[samCols.CIGAR]
                dataset.seq = cols[samCols.SEQ]
                dataset.mapq = int(cols[samCols.MAPQ])
                return dataset

    # return specific dataset with reduced memory footprint
    def getShortRecord(self, dsName):
        if dsName in self.records.keys():
            with open(self.filename, 'r', encoding="ascii", errors="surrogateescape") as sam:
                offset, flags, reference = self.records[dsName]
                sam.seek(offset)
                line = sam.readline()
                cols = line.split()
                dataset = samRecord()
                dataset.name = self.__simpleName__(cols[samCols.QNAME])
                dataset.flags = int(cols[samCols.FLAG])
                dataset.reference = cols[samCols.RNAME]
                dataset.pos = int(cols[samCols.POS])
                dataset.seq = cols[samCols.SEQ]
                return dataset

    def __simpleName__(self, name):
        return list(filter(None, re.split("[,.: ]+", name)))[0]
