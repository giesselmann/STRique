# \MODULE\---------------------------------------------------------------
#
#  CONTENTS      : Raw nanopore signal repeat detection pipeline
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
import os, sys, json, argparse
import Bio.Seq
import signal
import numpy as np
import numpy.ma as ma
import scipy.signal as sp
import pomegranate as pg
from skimage.morphology import opening, closing, rectangle
from multiprocessing import Pool, Process, Event, Value

# private imports
import fast5Index
import pyseqan


# basic normalization and simulation
class pore_model():
    def __init__(self, model_file):
        model_data = np.genfromtxt(model_file, names=True, dtype=None, comments="#", encoding='UTF-8')
        self.__model = {}
        for state in model_data:
            self.__model[state['kmer']] = (state['level_mean'], state['level_stdv'])
        self.kmer = len(model_data['kmer'][0])
        self.model_median = np.median(model_data['level_mean'])
        self.model_MAD = np.mean(np.absolute(np.subtract(model_data['level_mean'], self.model_median)))
        model_dict = {}
        for state in model_data:
            state_idx = state['kmer']
            model_dict[state_idx] = (state['level_mean'], state['level_stdv'])
        min_state = min(model_dict.values(), key=lambda x:x[0])
        max_state = max(model_dict.values(), key=lambda x:x[0])
        self.model_min = min_state[0] - 6 * min_state[1]
        self.model_max = max_state[0] + 6 * max_state[1]
        self.model_dict = model_dict

    def normalize2model(self, signal, clip=True):
        rawMedian = np.median(signal)
        rawMAD = np.mean(np.absolute(np.subtract(signal, rawMedian)))
        signal = np.divide(np.subtract(signal, np.median(signal)), rawMAD)
        signal = np.add(np.multiply(signal, self.model_MAD), self.model_median)
        if clip == True:
            np.clip(signal, self.model_min + .5, self.model_max - .5, out=signal)
        return signal

    def generate_signal(self, sequence, samples=10):
        signal = []
        level_means = np.array([self.__model[kmer][0] for kmer in [sequence[i:i+self.kmer] for i in range(len(sequence)-self.kmer + 1)]])
        return np.repeat(level_means, samples)


# parse bed file
class bed_parser():
    def __init__(self, bedFile):
        self.records = {}
        with open(bedFile, 'r') as fp:
            for line in fp:
                if line:
                    cols = line.strip().split()
                    if len(cols) >= 6:
                        chr, start, stop, name, val, strand = cols[:6]
                        self.records[name] = (chr, int(start), int(stop), strand)
        
    def __getitem__(self, key): 
        if key in self.records:
            return self.records[key]
        else:
            return None
        
        
# profile HMM base class
class profileHMM(pg.HiddenMarkovModel):         
    def __init__(self, sequence, 
                 model_file, transition_probs={}, state_prefix='', no_silent=False, std_scale=1.2):
        super().__init__()
        self.pore_model = pore_model(model_file)
        self.sequence = sequence
        self.state_prefix = state_prefix
        self.no_silent = no_silent
        self.std_scale = std_scale
        self.transition_probs = {'loop': .95,
                                 'move': .042,
                                 'insert': .006,
                                 'loop_insert' : .30,
                                 'delete': .002,
                                 'delete_multiple': .005}
        for key, value in transition_probs.items():
            self.transition_probs[key] = value
        self.__init_model__()

    def __init_model__(self):
        self.match_states, insertion_states, deletetion_states = self.__extract_states__(self.sequence)
        self.insertion_states = insertion_states
        self.deletion_states = deletetion_states
        self.s1 = pg.State(None, name=self.state_prefix+'s1')
        self.s2 = pg.State(None, name=self.state_prefix+'s2')
        self.e1 = pg.State(None, name=self.state_prefix+'e1')
        self.e2 = pg.State(None, name=self.state_prefix+'e2')
        self.__connect_states__()

    def __extract_states__(self, sequence):
        match_states = []
        insertion_states = []
        deletion_states = []
        digits = np.ceil(np.log10(len(sequence) - self.pore_model.kmer + 1)).astype(np.int)
        for idx, kmer in enumerate([sequence[i:i+self.pore_model.kmer] for i in range(len(sequence) - self.pore_model.kmer + 1)]):
            state_name = self.state_prefix + str(idx).rjust(digits,'0')
            state_mean, state_std = self.pore_model.model_dict[kmer]
            match_states.append(pg.State(pg.NormalDistribution(state_mean, state_std * self.std_scale), 
                                name=state_name + 'm'))
            deletion_states.append(pg.State(None, name=state_name + 'd'))
            insertion_states.append(pg.State(pg.UniformDistribution(self.pore_model.model_min, self.pore_model.model_max),
                                    name=state_name + 'i'))                                              
        return match_states, insertion_states, deletion_states

    def __connect_states__(self):
        self.add_states(self.match_states)
        self.add_states(self.insertion_states)
        if not self.no_silent:
            self.add_states(self.deletion_states)
        self.add_states([self.s1, self.s2, self.e1, self.e2])
        # matches
        for i, state in enumerate(self.match_states):
            self.add_transition(state, state, self.transition_probs['loop'], group='loop')
            if i < len(self.match_states) - 1:
                self.add_transition(state, self.match_states[i + 1], self.transition_probs['move'], group='move')
        # insertions
        for i, state in enumerate(self.insertion_states):
            self.add_transition(state, state, self.transition_probs['loop_insert'], group='loop_insert')
            self.add_transition(self.match_states[i], state, self.transition_probs['insert'], group='insert')
            if i < len(self.deletion_states) - 1 and not self.no_silent:
                self.add_transition(state, self.deletion_states[i+1], (1 - self.transition_probs['loop_insert']) * .1, group='loop_insert_n')
            if i < len(self.match_states) - 1:
                self.add_transition(state, self.match_states[i+1], (1 - self.transition_probs['loop_insert']) * .9, group='loop_insert_n')
        # deletions
        if not self.no_silent:
            for i, state in enumerate(self.deletion_states):
                self.add_transition(state, self.insertion_states[i], (1 - self.transition_probs['delete_multiple']) * .1, group='delete_multiple_n')
                if i > 0:
                    self.add_transition(self.match_states[i-1], state, self.transition_probs['delete'], group='delete')
                if i < len(self.match_states) - 1:
                    self.add_transition(state, self.match_states[i+1], (1 - self.transition_probs['delete_multiple']) * .9, group='delete_multiple_n')
                if i < len(self.deletion_states) - 1:
                    self.add_transition(state, self.deletion_states[i+1], self.transition_probs['delete_multiple'], group='delete_multiple')
            self.add_transition(self.s1, self.deletion_states[0], 1)
            self.add_transition(self.s2, self.match_states[0], 1)
            self.add_transition(self.deletion_states[-1], self.e1, self.transition_probs['delete_multiple'])
            self.add_transition(self.deletion_states[-1], self.e2, (1 - self.transition_probs['delete_multiple']))
        else:
            for i, state in enumerate(self.match_states):
                if i < len(self.match_states) - 2:
                    self.add_transition(state, self.match_states[i+2], self.transition_probs['delete'], group='delete')
            self.add_transition(self.s1, self.insertion_states[0], 1)
            self.add_transition(self.s2, self.match_states[0], 1)
        self.add_transition(self.insertion_states[-1], self.e1, (1 - self.transition_probs['loop_insert']) * .1)
        self.add_transition(self.insertion_states[-1], self.e2, (1 - self.transition_probs['loop_insert']) * .9)
        self.add_transition(self.match_states[-1], self.e2, self.transition_probs['move'])
        self.add_transition(self.match_states[-1], self.e1, self.transition_probs['delete'])

    def bake(self, *args, **kwargs):
        self.add_transition(self.start, self.s1, .5)
        self.add_transition(self.start, self.s2, .5)
        self.add_transition(self.e1, self.end, 1)
        self.add_transition(self.e2, self.end, 1)
        super().bake(*args, **kwargs)


# repeat count profile HMM
class repeatHMM(pg.HiddenMarkovModel):
    def __init__(self, repeat, model_file, transition_probs={}, state_prefix='', std_scale=1.2):
        super().__init__()
        self.repeat = repeat
        self.model_file = model_file
        self.transition_probs = {'skip': .999,   # 99
                                 'leave_repeat': .0002,
                                 'loop': .95,
                                 'move': .042,
                                 'insert': .006,
                                 'loop_insert' : .30,
                                 'delete': .002}
        for key, value in transition_probs.items():
            self.transition_probs[key] = value
        self.state_prefix = state_prefix
        self.std_scale = std_scale
        self.__build_model__()

    def __build_model__(self):
        self.pore_model = pore_model(self.model_file)
        repeat = ''.join([self.repeat] * int(np.ceil((self.pore_model.kmer + 2) / len(self.repeat))))[:-1] # TODO check this expansion on shorter repeats
        self.repeat_hmm = profileHMM(repeat, self.model_file, transition_probs=self.transition_probs, state_prefix=self.state_prefix, 
                                     no_silent=True, std_scale=self.std_scale)
        self.add_model(self.repeat_hmm)
        self.skip_distribution = pg.NormalDistribution(self.pore_model.model_median, self.pore_model.model_MAD)
        self.dummy_distribution = pg.UniformDistribution(self.pore_model.model_min, self.pore_model.model_max)
        self.d1 = pg.State(self.dummy_distribution, name=self.state_prefix+'dummy1')
        self.d2 = pg.State(self.dummy_distribution, name=self.state_prefix+'dummy2')
        self.e1 = pg.State(None, name=self.state_prefix+'e1')
        self.e2 = pg.State(None, name=self.state_prefix+'e2')
        self.add_state(self.d1)
        self.add_state(self.d2)
        self.s1 = self.repeat_hmm.s1
        self.s2 = self.repeat_hmm.s2
        self.add_transition(self.repeat_hmm.e1, self.d1, 1)
        self.add_transition(self.repeat_hmm.e2, self.d2, 1)
        self.add_transition(self.d1, self.e1, self.transition_probs['leave_repeat'])
        self.add_transition(self.d2, self.e2, self.transition_probs['leave_repeat'])
        self.add_transition(self.d1, self.s1, 1 - self.transition_probs['leave_repeat'])
        self.add_transition(self.d2, self.s2, 1 - self.transition_probs['leave_repeat'])

    def bake(self, *args, **kwargs):
        self.add_transition(self.start, self.s1, .5)
        self.add_transition(self.start, self.s2, .5)
        self.add_transition(self.e1, self.end, 1)
        self.add_transition(self.e2, self.end, 1)
        super().bake(*args, **kwargs)

    def bake2(self, *args, **kwargs):
        self.skip = pg.State(self.skip_distribution, name=self.state_prefix+'skip')
        self.add_state(self.skip)
        self.add_transition(self.start, self.s1, .5)
        self.add_transition(self.start, self.s2, .5)
        self.add_transition(self.e1, self.skip, 1)
        self.add_transition(self.e2, self.skip, 1)
        self.add_transition(self.skip, self.skip, self.transition_probs['skip'])
        self.add_transition(self.skip, self.end, 1 - self.transition_probs['skip'])
        super().bake(*args, **kwargs)

    def count_repeats(self, states):
        states = np.array([x[1] for x in states])
        n1 = np.sum(states == self.d1)
        n2 = np.sum(states == self.d2)
        return n1 + n2


# repeat detection profile HMM
class embeddedRepeatHMM(pg.HiddenMarkovModel):
    def __init__(self, repeat, 
                 prefix, suffix,
                 model_file):
        super().__init__()
        self.transition_probs = {'skip': 1-1e-4}
        self.model_file = model_file
        self.repeat = repeat
        self.prefix = prefix
        self.suffix = suffix
        self.__build_model__()

    def __build_model__(self):
        self.pore_model = pore_model(self.model_file)
        # expand primer sequences and get profile HMMs
        prefix = self.prefix + self.repeat[:-1]
        suffix = self.repeat + self.suffix  
        self.prefix_model = profileHMM(prefix, self.model_file, self.transition_probs, state_prefix='prefix', std_scale=1.0)
        self.suffix_model = profileHMM(suffix, self.model_file, self.transition_probs, state_prefix='suffix', std_scale=1.0)
        self.repeat_model = repeatHMM(self.repeat, self.model_file, self.transition_probs, state_prefix='repeat', std_scale=1.0)
        # add sub-modules, flanking and skip states
        self.add_model(self.prefix_model)
        self.add_model(self.repeat_model)
        self.add_model(self.suffix_model)
        self.add_transition(self.start, self.prefix_model.s1, .1)
        self.add_transition(self.start, self.prefix_model.s2, .9)
        # repeat model
        self.add_transition(self.prefix_model.e1, self.repeat_model.s1, 1)
        self.add_transition(self.prefix_model.e2, self.repeat_model.s2, 1)
        # suffix model
        self.add_transition(self.repeat_model.e1, self.suffix_model.s1, 1)
        self.add_transition(self.repeat_model.e2, self.suffix_model.s2, 1)
        self.add_transition(self.suffix_model.e1, self.end, 1)
        self.add_transition(self.suffix_model.e2, self.end, 1)
        # bake and store state IDs 
        self.bake(merge='All')

    def count_repeats(self, sequence, **kwargs):
        p, path = super().viterbi(sequence, **kwargs)
        if path is not None:
            # repeat number equals loops in repeat model + 1 repeat in flanking sequences
            n = self.repeat_model.count_repeats(path) + 1
            path = np.array([x[0] for x in path])
            path = path[path < self.silent_start]
            return n, p
        else:
            return 0, 0


# main repeat detection methods
class repeatDetection(object):
    def __init__(self, model_file,
                 repeat, prefix, suffix,
                 prefix2=None, suffix2=None):
        self.algn = pyseqan.align_raw()
        self.algn.dist_offset = 8.0
        self.algn.gap_open = -2.0           # -16.0
        self.algn.gap_extension = -8.0     # -4.0
        self.algn.dist_min = -16.0
        self.pm = pore_model(model_file)
        if not prefix2:
            prefix2 = prefix
        if not suffix2:
            suffix2 = suffix
        self.score_min = 200
        self.samples = 10
        self.sim_prefix = self.pm.generate_signal(prefix, samples=self.samples)
        self.sim_suffix = self.pm.generate_signal(suffix, samples=self.samples)
        self.sim_prefix2 = self.pm.generate_signal(prefix2, samples=self.samples)
        self.sim_suffix2 = self.pm.generate_signal(suffix2, samples=self.samples)
        self.flanked_model = embeddedRepeatHMM(repeat, prefix, suffix, model_file)

    def detect_range(self, signal, segment):
        score, idx_signal, idx_segment = self.algn.align_overlap(signal, segment)
        segment_begin = np.abs(np.array(idx_signal) - idx_segment[0]).argmin()
        segment_end = np.abs(np.array(idx_signal) - idx_segment[-1]).argmin()
        return score, segment_begin, segment_end

    def detect_short(self, segment):
        return self.flanked_model.count_repeats(segment)

    def detect(self, signal):
        nrm_signal = sp.medfilt(signal, kernel_size=3)
        nrm_signal = (nrm_signal - np.median(nrm_signal)) / np.std(nrm_signal)
        nrm_signal = np.clip(nrm_signal * 42 + 127, 0, 255).astype(np.dtype('uint8')).reshape((1, len(nrm_signal)))
        flt = rectangle(1, 6)
        nrm_signal = opening(nrm_signal, flt)
        nrm_signal = closing(nrm_signal, flt)
        nrm_signal = self.pm.normalize2model(nrm_signal[0].astype(np.dtype('float')), clip=True)
        score_prefix, prefix_begin, prefix_end = self.detect_range(nrm_signal, self.sim_prefix2)
        score_suffix, suffix_begin, suffix_end = self.detect_range(nrm_signal, self.sim_suffix2)
        n = 0; p = 0        
        # if prefix_end < suffix_begin and score_prefix > self.score_min and score_suffix > self.score_min:
            # n, p = self.detect_short(nrm_signal[prefix_begin:suffix_end])
        if prefix_end < suffix_begin:
            sp2, hmm_begin, _ = self.detect_range(nrm_signal[prefix_begin:suffix_end], self.sim_prefix)
            ss2, _, hmm_end = self.detect_range(nrm_signal[prefix_begin:suffix_end], self.sim_suffix)
            if prefix_begin + hmm_begin < prefix_begin + hmm_end and score_prefix > self.score_min and score_suffix > self.score_min:
                n, p = self.detect_short(nrm_signal[prefix_begin+hmm_begin:prefix_begin+hmm_end])
        return n, score_prefix, score_suffix, p, suffix_begin - prefix_end, prefix_end


# main repeat detection worker
def repeat_detection(config, record_IDs, output_file, counter):
    try:
        # load config
        fast5_dir=config['fast5_dir']
        bed_file=config['bed_file']
        output_file=os.path.join(config['out_dir'], output_file)
        model_file=config['model_file']
        repeat=config['repeat']
        prefix=config['prefix']
        suffix=config['suffix']
        prefix2=config['prefix2']
        suffix2=config['suffix2']
        # init
        assert(isinstance(repeat, dict))
        assert(isinstance(prefix, dict))
        assert(isinstance(suffix, dict))
        bed = None
        if bed_file:
            bed = bed_parser(bed_file)
        f5 = fast5Index.fast5Index()
        f5.load(fast5_dir)
        pm = pore_model(model_file)
        # init profile HMMs, simulate signals
        detection = {}
        for key, value in repeat.items():
            detection[key] = (repeatDetection(model_file=model_file,
                                               repeat=value, 
                                               prefix=prefix[key], 
                                               suffix=suffix[key], 
                                               prefix2=prefix2[key],
                                               suffix2=suffix2[key]),
                               repeatDetection(model_file=model_file,
                                               repeat=Bio.Seq.reverse_complement(value), 
                                               prefix=Bio.Seq.reverse_complement(suffix[key]),
                                               suffix=Bio.Seq.reverse_complement(prefix[key]),
                                               prefix2=Bio.Seq.reverse_complement(suffix2[key]),
                                               suffix2=Bio.Seq.reverse_complement(prefix2[key]))
                               )
        # write header line
        with open(output_file, 'w') as fp:
            print('\t'.join(['ID', 'ref', 'flag', 'count', 'score_prefix', 'score_suffix', 'log_p', 'ticks', 'offset']), file=fp)
        # analyze each read
        for record_ID in record_IDs:
            f5_record = f5.getRecord(record_ID)
            if f5_record is None:
                continue
            raw_signal = f5_record.raw    
            if bed:
                bed_record = bed[record_ID]
                if not bed_record:
                    continue
                if not bed_record[0] in detection:
                    continue
                else:
                    if bed_record[3] == '+':
                        model = detection[bed_record[0]][0]
                    elif bed_record[3] == '-':
                        model = detection[bed_record[0]][1]
                    else:
                        continue
                n, score_prefix, score_suffix, log_p, ticks, offset = model.detect(raw_signal)
                with open(output_file, 'a+') as fp:
                    print('\t'.join([f5_record.ID, bed_record[0], bed_record[3], str(n), 
                                     str(score_prefix), str(score_suffix), str(log_p), str(ticks), str(offset)]), file=fp)
                with counter.get_lock():
                    counter.value += 1
            else:
                for key, value in repeat.items():
                    model = detection[key][0]
                    n0, score_prefix0, score_suffix0, log_p0, ticks0, offset0 = model.detect(raw_signal)
                    model = detection[key][1]
                    n1, score_prefix1, score_suffix1, log_p1, ticks1, offset1 = model.detect(raw_signal)
                    score0 = score_prefix0 + score_suffix0
                    score1 = score_prefix1 + score_suffix1
                    if score0 > score1 and ticks0 > 0:
                        with open(output_file, 'a+') as fp:
                            print('\t'.join([f5_record.ID, key, '+', str(n0), 
                                         str(score_prefix0), str(score_suffix0), str(log_p0), str(ticks0), str(offset0)]), file=fp) 
                    elif score1 > score0 and ticks1 > 0:
                        with open(output_file, 'a+') as fp:
                            print('\t'.join([f5_record.ID, key, '-', str(n1), 
                                         str(score_prefix1), str(score_suffix1), str(log_p1), str(ticks1), str(offset1)]), file=fp) 
                with counter.get_lock():
                    counter.value += 1       

    except KeyboardInterrupt:
        return

            
# parse config.json            
def parse_config(f5_dir, model_file, config_file, bed_file, out_dir):
    with open(config_file) as fp:
        ld_conf = json.load(fp)
    config = {}
    try:
        config['model_file'] = model_file
        config['fast5_dir'] = f5_dir
        if os.path.isfile(bed_file):
            config['bed_file'] = bed_file
        else:
            config['bed_file'] = None
        config['out_dir'] = out_dir
        config['repeat'] = ld_conf['repeat']
        config['prefix'] = ld_conf['prefix']
        config['suffix'] = ld_conf['suffix']
        if 'prefix2' in ld_conf:
            config['prefix2'] = ld_conf['prefix2']
        else:
            config['prefix2'] = ld_conf['prefix']
        if 'suffix2' in ld_conf:
            config['suffix2'] = ld_conf['suffix2']
        else:
            config['suffix2'] = ld_conf['suffix']
    except KeyError as e:
        print('Error loading config file, missing', e.args[0])
        exit(1)
    return config
        
 
# print progess on command line
def process_monitor(counter, max_count, stop_event):
    import time
    try:
        while not stop_event.is_set():
            time.sleep(1)
            print("\rProcessed ", counter.value, 'of ', (max_count), end=" ")
    except KeyboardInterrupt:
        return
    except:
        print('Unknown Exception in Montior Process')
        
        
if __name__ == '__main__':
    # command line
    parser = argparse.ArgumentParser(description="STR Detection in raw nanopore data")
    parser.add_argument("f5", help="fast5 file directory") 
    parser.add_argument("model", help="pore model")
    parser.add_argument("config", help="job config file")
    parser.add_argument("out", help="output directory")
    parser.add_argument("--bed", help="alignment information in BED-6 format")
    parser.add_argument("--t", type=int, default=1, help="Number of processes to use in parallel")
    args = parser.parse_args()
    # load config
    config = parse_config(args.f5, args.model, args.config, args.bed, args.out)
    # index/load reads
    f5 = fast5Index.fast5Index()
    f5.loadOrAnalyze(config['fast5_dir'], recursive=True)
    records = np.array(f5.getRecordNames())
    # split in chunks for multi-processing
    record_chunks = np.array_split(records, args.t)
    counter = Value('i', 0)
    stop_event = Event()
    monitor = Process(target=process_monitor, args=(counter, len(records), stop_event))
    monitor.start()
    worker = []
    # fork into multiple worker processes
    try:
        if len(record_chunks) > 1:
            for i, chunk in enumerate(record_chunks[1:]):
                worker.append(Process(target=repeat_detection, args=(config, chunk, str(i) + '.tsv', counter)))
            for w in worker:
               w.start()
        repeat_detection(config, record_chunks[0], '0.tsv', counter)
        for w in worker:
           w.join()
    except KeyboardInterrupt:
        print("Interrupt by user")
    stop_event.set()
    try:
        monitor.join()
    except RuntimeError:
        print('Could not join monitor thread')
    
    print("")
    
    