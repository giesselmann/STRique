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
# public imports
import os, sys, glob
import re
import shutil
import json
import tarfile
# private imports
import fast5Reader


class fast5Index(object):
    def __init__(self):
        self.__records = {}
        self.__path2fast5 = ""
        self.__f5Reader = fast5Reader.fast5Reader()

    # read files, store names for reads
    def analyze(self, path2fast5, recursive=False):
        self.__path2fast5 = os.path.abspath(path2fast5)
        if not recursive:
            f5files = glob.glob(os.path.join(path2fast5, '*.fast5'))
            tarFiles = glob.glob(os.path.join(path2fast5, '*.tar'))
        else:
            f5files = [os.path.join(dirpath, f) for dirpath, _, files in os.walk(path2fast5) for f in files if f.endswith('.fast5')]
            tarFiles = [os.path.join(dirpath, f) for dirpath, _, files in os.walk(path2fast5) for f in files if f.endswith('.tar')]
        for f in f5files:
            try:
                self.__records[self.__f5Reader.getRecordID(f)] = os.path.relpath(f, start=path2fast5).replace("\\", "/")
            except:
                sys.stderr.write("Error reading " + f + '\n')
        for f in tarFiles:
            print("Indexing", os.path.basename(f))
            try:
                with tarfile.open(f) as tar:
                    for member in tar.getmembers():   
                        try:
                            if member.isfile():
                                name = member.name
                                member.name = os.path.basename(member.name)
                                tar.extract(member, path=path2fast5)
                                ID = self.__f5Reader.getRecordID(os.path.join(path2fast5, member.name))
                                os.remove(os.path.join(path2fast5, member.name))
                                self.__records[ID] = '/'.join([os.path.relpath(f, start=path2fast5), name])
                        except:
                            sys.stderr.write("Error reading " + name + '\n')
                        pass
            except:
                sys.stderr.write("Error reading archive " + os.path.basename(f) + '\n')
        pass

    # load dictionary from disk
    def load(self, path2fast5):
        self.__path2fast5 = path2fast5
        self.__records = {}
        with open(path2fast5 + '/index.json', 'r') as fp:
            self.__records = json.load(fp)

    # load index from disk, create if not existent
    def loadOrAnalyze(self, path2fast5, recursive=False):
        if os.path.isfile(path2fast5 + '/index.json'):
            self.load(path2fast5)
        else:
            self.analyze(path2fast5, recursive=recursive)
            self.save(path2fast5)

    # save dictionary to disk
    def save(self, path2fast5):
        with open(path2fast5 + '/index.json', 'w') as fp:
            json.dump(self.__records, fp)

    # return names of valid records in current index
    def getRecordNames(self):
        return sorted(list(self.__records.keys()))

    # retrun name, filename
    def getFilenameDict(self):
        return dict((v,k) for k, v in self.__records.items())

    def extract(self, name, destination):
        if name in self.__records.keys():
            loc = self.__records[name]
            if not ".tar" in loc:
                shutil.copy2(os.path.join(self.__path2fast5, loc), destination)
            else:
                f, _, tarPath = loc.rpartition(".tar/")
                tar = tarfile.open(os.path.join(self.__path2fast5, f + ".tar"))
                member = tar.getmember(tarPath)
                member.name = os.path.basename(member.name)
                tar.extract(member, path=destination)
            return name

    # return fast5 record
    def getRecord(self, name):
        if name in self.__records.keys():
            loc = self.__records[name]
            if not ".tar" in loc:
                return self.__f5Reader.getRecord(os.path.join(self.__path2fast5, loc))
            else:
                f, _, tarPath = loc.rpartition(".tar/")
                tar = tarfile.open(os.path.join(self.__path2fast5, f + ".tar"))
                member = tar.getmember(tarPath)
                member.name = os.path.basename(member.name)
                tar.extract(member, path=self.__path2fast5)
                record = self.__f5Reader.getRecord(os.path.join(self.__path2fast5, member.name))
                os.remove(os.path.join(self.__path2fast5, member.name))
                return record


if __name__ == '__main__':
    f5Path = "D:/data/runs/20170626_plasmid_c9orf72_RAD002/reads"
    f5Index = fast5Index()
    f5Index.loadOrAnalyze(f5Path, recursive=True)
    pass