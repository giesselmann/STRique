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
import queue
import itertools
import numpy as np
import numpy.ma as ma
import scipy.signal as sp
import pomegranate as pg
from skimage.morphology import opening, closing, dilation, erosion, rectangle
import matplotlib.pyplot as plt
from multiprocessing import Pool, Process, Event, Value, Queue


# private imports
import fast5Index
import pyseqan




def sliding_window(a, n=3, mode='same'):
    if mode == 'mean':
        a = np.append(a, (n-1) * [np.mean(a)])
    elif mode == 'median':
        a = np.append(a, (n-1) * [np.median(a)])
    elif mode == 'mirror':
        a = np.append(a, a[-1:-1-(n-1):-1])
    else:
        a = np.append(a, (n-1) * [a[-1]])
    shape = a.shape[:-1] + (a.shape[-1] - n + 1, n)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)




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
        
    def MAD(self, signal):
        return np.mean(np.absolute(np.subtract(signal, np.median(signal))))

    def normalize2model(self, signal, clip=True, mode='median', plot=False):
        if mode == 'minmax':
            model_values = np.array([x[0] for x in self.model_dict.values()])
            q5_sig, q95_sig = np.percentile(signal, [1, 99])
            q5_mod, q95_mod = np.percentile(model_values, [1, 99])
            m5_sig = np.median(signal[signal < q5_sig])
            m95_sig = np.median(signal[signal > q95_sig])
            m5_mod = np.median(model_values[model_values < q5_mod])
            m95_mod = np.median(model_values[model_values > q95_mod])
            nrm_signal = (signal - (m5_sig + (m95_sig - m5_sig) / 2)) / ((m95_sig - m5_sig) / 2)
            nrm_signal = nrm_signal * ((m95_mod - m5_mod) / 2) + (m5_mod + (m95_mod - m5_mod) / 2)
        elif mode == 'entropy':
            #sliding_std = [np.std(x) for x in sliding_window(signal, n=1000, mode='mirror')]
            sliding_std = [self.MAD(x) for x in sliding_window(signal, n=500, mode='mirror')]
            sliding_std += [sliding_std[-1]]
            diff_signal = np.abs(np.diff(sliding_std))
            ind = np.argpartition(diff_signal, -50)[-50:]
            # q = np.percentile(diff_signal, 97)  # 99.5
            # diff_mask = np.array([0 if x < q else 1 for x in diff_signal], dtype=np.dtype('uint8'))
            diff_mask = np.zeros(len(diff_signal), dtype=np.dtype('uint8'))
            diff_mask[ind] = 1
            diff_mask = dilation(diff_mask.reshape((1, len(diff_mask))), rectangle(1, 750))[0].astype(np.bool)
            rawMedian = np.median(signal[diff_mask])
            rawMAD = self.MAD(signal[diff_mask])
            nrm_signal = np.divide(np.subtract(signal, rawMedian), rawMAD)
            nrm_signal = np.add(np.multiply(nrm_signal, self.model_MAD), self.model_median)
        else:
            rawMedian = np.median(signal)
            rawMAD = self.MAD(signal)
            nrm_signal = np.divide(np.subtract(signal, rawMedian), rawMAD)
            nrm_signal = np.add(np.multiply(nrm_signal, self.model_MAD), self.model_median)
        if clip == True:
            np.clip(nrm_signal, self.model_min + .5, self.model_max - .5, out=nrm_signal)
        #if mask and plot:
            # f, ax = plt.subplots(3, sharex=True)
            # ax[0].plot(nrm_signal, 'k-')
            # ax[0].axhline(rawMedian, color='r')
            # ax[0].axhline(np.median(signal), color='g')
            # ax[1].plot(diff_mask, 'k-')
            # ax[2].plot(diff_signal, 'k-')
            # plt.show()
            # f, ax = plt.subplots(3, sharex=False)
            # bins = np.linspace(self.model_min, self.model_max, 75)
            # ax[0].plot(nrm_signal, 'b-')
            # ax[1].hist(nrm_signal, bins=bins)
            # ax[1].hold(False)
            # ax[2].hist([x[0] for x in self.model_dict.values()], bins=bins)
            # ax[0].callbacks.connect('xlim_changed', lambda axes, ax1=ax[1], nrm_signal=nrm_signal, bins=bins: ax1.hist(nrm_signal[max([int(axes.get_xlim()[0]), 0]) : int(axes.get_xlim()[1])], bins=bins))
            # plt.show()
        return nrm_signal

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




# profile HMM
class profileHMM(pg.HiddenMarkovModel):
    def __init__(self, sequence,
                 pm_base, transition_probs={}, state_prefix='',
                 no_silent=False,
                 std_scale=1.0, std_offset=0.0
                 ):
        super().__init__()
        self.pm_base = pm_base
        self.sequence = sequence
        self.state_prefix = state_prefix
        self.no_silent = no_silent
        self.std_scale = std_scale
        self.std_offset = std_offset
        self.transition_probs = {'match_loop': .75,     #
                                 'match_match': .15,    # sum to 1
                                 'match_insert': .09,   #
                                 'match_delete': .01,   #

                                 'insert_loop' : .15,   #
                                 'insert_match_0': .40, # sum to 1
                                 'insert_match_1': .40, #
                                 'insert_delete': .05,  #

                                 'delete_delete': .005, #
                                 'delete_insert': .05,  # sum to 1
                                 'delete_match': .945   #
                                 }
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
        digits = np.ceil(np.log10(len(sequence) - self.pm_base.kmer + 1)).astype(np.int)
        for idx, kmer in enumerate([sequence[i:i+self.pm_base.kmer] for i in range(len(sequence) - self.pm_base.kmer + 1)]):
            state_name = self.state_prefix + str(idx).rjust(digits,'0')
            state_mean, state_std = self.pm_base.model_dict[kmer]
            match_states.append(pg.State(pg.NormalDistribution(state_mean, state_std * self.std_scale + self.std_offset), name=state_name + 'm'))
            if not self.no_silent:
                deletion_states.append(pg.State(None, name=state_name + 'd'))
            insertion_states.append(pg.State(pg.UniformDistribution(self.pm_base.model_min, self.pm_base.model_max),
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
            self.add_transition(state, state, self.transition_probs['match_loop'], group='match_loop')
            if i < len(self.match_states) - 1:
                self.add_transition(state, self.match_states[i + 1], self.transition_probs['match_match'], group='match_match')
        # insertions
        for i, state in enumerate(self.insertion_states):
            self.add_transition(state, state, self.transition_probs['insert_loop'], group='insert_loop')
            self.add_transition(self.match_states[i], state, self.transition_probs['match_insert'], group='match_insert')
            self.add_transition(state, self.match_states[i], self.transition_probs['insert_match_1'], group='insert_match_1')
            if i < len(self.deletion_states) - 1 and not self.no_silent:
                self.add_transition(state, self.deletion_states[i+1], self.transition_probs['insert_delete'], group='insert_delete')
            if i < len(self.match_states) - 1:
                self.add_transition(state, self.match_states[i+1], self.transition_probs['insert_match_0'], group='insert_match_0')
        # deletions
        if not self.no_silent:
            for i, state in enumerate(self.deletion_states):
                self.add_transition(state, self.insertion_states[i], self.transition_probs['delete_insert'], group='delete_insert')
                if i > 0:
                    self.add_transition(self.match_states[i-1], state, self.transition_probs['match_delete'], group='match_delete')
                if i < len(self.match_states) - 1:
                    self.add_transition(state, self.match_states[i+1], self.transition_probs['delete_match'], group='delete_match')
                if i < len(self.deletion_states) - 1:
                    self.add_transition(state, self.deletion_states[i+1], self.transition_probs['delete_delete'], group='delete_delete')
            self.add_transition(self.s1, self.deletion_states[0], 1)
            self.add_transition(self.s2, self.match_states[0], 1)
            self.add_transition(self.deletion_states[-1], self.e1, self.transition_probs['delete_delete'])
            self.add_transition(self.deletion_states[-1], self.e2, self.transition_probs['delete_match'])
        else:
            for i, state in enumerate(self.match_states):
                if i < len(self.match_states) - 2:
                    self.add_transition(state, self.match_states[i+2], self.transition_probs['match_delete'], group='match_delete')
            self.add_transition(self.s1, self.insertion_states[0], 1)
            self.add_transition(self.s2, self.match_states[0], 1)
        self.add_transition(self.insertion_states[-1], self.e1, self.transition_probs['insert_delete'], group='insert_delete')
        self.add_transition(self.insertion_states[-1], self.e2, self.transition_probs['insert_match_0'], group='insert_match_0')
        self.add_transition(self.match_states[-1], self.e2, self.transition_probs['match_match'])
        self.add_transition(self.match_states[-1], self.e1, self.transition_probs['match_delete'])

    def bake(self, *args, **kwargs):
        self.add_transition(self.start, self.s1, .5)
        self.add_transition(self.start, self.s2, .5)
        self.add_transition(self.e1, self.end, 1)
        self.add_transition(self.e2, self.end, 1)
        super().bake(*args, **kwargs)




# repeat count profile HMM
class repeatHMM(pg.HiddenMarkovModel):
    def __init__(self, repeat, model_file, transition_probs={}, state_prefix='', std_scale=1.0, std_offset=0.0):
        super().__init__()
        self.repeat = repeat
        self.model_file = model_file
        # self.transition_probs = {'skip': .999,   # 99
                                 # 'leave_repeat': .0002,
                                 # 'loop': .95,
                                 # 'move': .042,
                                 # 'insert': .006,
                                 # 'loop_insert' : .30,
                                 # 'delete': .002}
        self.transition_probs = {'skip': .999,   # 99
                                 'leave_repeat': .002
                                 }
        for key, value in transition_probs.items():
            self.transition_probs[key] = value
        self.state_prefix = state_prefix
        self.std_scale = std_scale
        self.std_offset = std_offset
        self.__build_model__()

    def __build_model__(self):
        self.pore_model = pore_model(self.model_file)
        if len(self.repeat) >= self.pore_model.kmer:
            repeat = self.repeat + self.repeat[: self.pore_model.kmer - 1]
            self.repeat_offset = 0
        else:
            ext = self.pore_model.kmer - 1 + (len(self.repeat) - 1) - ((self.pore_model.kmer - 1) % len(self.repeat))
            repeat = self.repeat + ''.join([self.repeat] * self.pore_model.kmer)[:ext]
            self.repeat_offset = int(len(repeat) / len(self.repeat)) - 1
        self.repeat_hmm = profileHMM(repeat,self.pore_model, transition_probs=self.transition_probs, state_prefix=self.state_prefix, 
                                     no_silent=True, std_scale=self.std_scale, std_offset=self.std_offset)
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
        return n1 + n2 - self.repeat_offset




# repeat detection profile HMM
class embeddedRepeatHMM(pg.HiddenMarkovModel):
    def __init__(self, repeat, 
                 prefix, suffix,
                 model_file, config=None):
        super().__init__()
        self.transition_probs = {'skip': 1-1e-4,
                                 'seq_std_scale': 1.0,
                                 'rep_std_scale': 1.0,
                                 'seq_std_offset': 0.0,
                                 'rep_std_offset': 0.0,
                                 'e1_ratio': 0.1}
        if config and isinstance(config, dict):
            for key, value in config.items():
                self.transition_probs[key] = value
        self.model_file = model_file
        self.repeat = repeat
        self.prefix = prefix
        self.suffix = suffix
        self.__build_model__()

    def __build_model__(self):
        self.pore_model = pore_model(self.model_file)
        # expand primer sequences and get profile HMMs
        prefix = self.prefix + ''.join([self.repeat] * int(np.ceil(self.pore_model.kmer / len(self.repeat))))[:-1]
        suffix = ''.join([self.repeat] * int(np.ceil(self.pore_model.kmer / len(self.repeat)))) + self.suffix
        self.flanking_count = int(np.ceil(self.pore_model.kmer / len(self.repeat))) * 2 - 1
        self.prefix_model = profileHMM(prefix, self.pore_model, self.transition_probs, state_prefix='prefix', std_scale=self.transition_probs['seq_std_scale'], std_offset=self.transition_probs['seq_std_offset'])
        self.suffix_model = profileHMM(suffix, self.pore_model, self.transition_probs, state_prefix='suffix', std_scale=self.transition_probs['rep_std_scale'], std_offset=self.transition_probs['seq_std_offset'])
        self.repeat_model = repeatHMM(self.repeat, self.model_file, self.transition_probs, state_prefix='repeat', std_scale=self.transition_probs['seq_std_scale'], std_offset=self.transition_probs['rep_std_offset'])
        # add sub-modules, flanking and skip states
        self.add_model(self.prefix_model)
        self.add_model(self.repeat_model)
        self.add_model(self.suffix_model)
        self.add_transition(self.start, self.prefix_model.s1, self.transition_probs['e1_ratio'])
        self.add_transition(self.start, self.prefix_model.s2, (1-self.transition_probs['e1_ratio']))
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
            n = self.repeat_model.count_repeats(path) + self.flanking_count
            path = np.array([x[0] for x in path])
            path = path[path < self.silent_start]
            return n, p, path
        else:
            return 0, 0, np.array([])




# main repeat detection methods
class repeatDetection(object):
    def __init__(self, model_file,
                 repeat, prefix, suffix,
                 prefix2=None, suffix2=None, align_config=None, HMM_config=None):
        default_config =  {'dist_offset': 16.0,
                             'dist_min': 0.0,
                             'gap_open_h': -1.0,
                             'gap_open_v': -16.0,
                             'gap_extension_h': -1.0,
                             'gap_extension_v': -16.0,
                             'samples': 6}
        if align_config and isinstance(align_config, dict):
            for key, value in align_config.items():
                default_config[key] = value
        self.algn = pyseqan.align_raw()
        self.algn.dist_offset = default_config['dist_offset']
        self.algn.dist_min = default_config['dist_min']
        self.algn.gap_open_h = default_config['gap_open_h']
        self.algn.gap_open_v = default_config['gap_open_v']
        self.algn.gap_extension_h = default_config['gap_extension_h']
        self.algn.gap_extension_v = default_config['gap_extension_v']
        self.pm = pore_model(model_file)
        if not prefix2:
            prefix2 = prefix
        if not suffix2:
            suffix2 = suffix
        self.samples = default_config['samples']
        self.sim_prefix = self.pm.generate_signal(prefix, samples=self.samples)
        self.sim_suffix = self.pm.generate_signal(suffix, samples=self.samples)
        self.sim_prefix2 = self.pm.generate_signal(prefix2, samples=self.samples)
        self.sim_suffix2 = self.pm.generate_signal(suffix2, samples=self.samples)
        self.flanked_model = embeddedRepeatHMM(repeat, prefix, suffix, model_file, HMM_config)

    def detect_range(self, signal, segment, pre_trim=0, post_trim=0):
        score, idx_signal, idx_segment = self.algn.align_overlap(signal, segment)
        segment_begin = np.abs(np.array(idx_signal) - idx_segment[0]).argmin()
        segment_end = np.abs(np.array(idx_signal) - idx_segment[-1]).argmin()
        score = score / (segment_end - segment_begin)
        # trim
        segment_begin = np.abs(np.array(idx_signal) - idx_segment[0 + pre_trim]).argmin()
        segment_end = np.abs(np.array(idx_signal) - idx_segment[-1 - post_trim]).argmin()
        return score, segment_begin, segment_end

    def detect_short(self, segment):
        return self.flanked_model.count_repeats(segment)

    def detect(self, signal, record_ID=''):
        flt_signal = sp.medfilt(signal, kernel_size=3)
        nrm_signal = (flt_signal - np.median(flt_signal)) / self.pm.MAD(flt_signal)
        nrm_signal = np.clip(nrm_signal * 24 + 127, 0, 255).astype(np.dtype('uint8')).reshape((1, len(nrm_signal)))
        flt = rectangle(1, 8)
        nrm_signal = opening(nrm_signal, flt)
        nrm_signal = closing(nrm_signal, flt)[0].astype(np.dtype('float'))
        nrm_signal = self.pm.normalize2model(nrm_signal.astype(np.dtype('float')), mode='minmax', plot=True)
        flt_signal = self.pm.normalize2model(flt_signal.astype(np.dtype('float')), mode='minmax')
        trim_prefix = len(self.sim_prefix2) - len(self.sim_prefix)
        trim_suffix = len(self.sim_suffix2) - len(self.sim_suffix)
        score_prefix, prefix_begin, prefix_end = self.detect_range(nrm_signal, self.sim_prefix2, pre_trim=trim_prefix)
        score_suffix, suffix_begin, suffix_end = self.detect_range(nrm_signal, self.sim_suffix2, post_trim=trim_suffix)
        n = 0; p = 0; path = np.array([])
        if prefix_end < suffix_begin:
            n, p, path = self.detect_short(flt_signal[prefix_begin:suffix_end])
        # if n == 0:
            # f, axs = plt.subplots(2, sharex=True)
            # ax = axs[0]
            # ax.plot(flt_signal, 'k-')
            # ax2 = ax.twinx()
            # ax.plot(nrm_signal, 'k-', alpha=0.3)
            # ax.axvline(prefix_begin, color='r')
            # ax.axvline(prefix_end, color='r')
            # ax.axvline(suffix_begin, color='b')
            # ax.axvline(suffix_end, color='b')
            # ax.set_title('Detected ' + str(n) + ' repeats ' + record_ID)
            # ax = axs[1]
            # ax.plot(np.arange(len(path)) + prefix_begin, path, 'b-')
            # plt.show()
        return n, score_prefix, score_suffix, p, suffix_begin - prefix_end, prefix_end




# main repeat detection worker
def process_detection(config, record_IDs, output_queue, counter):
    try:
        # load config
        fast5_dir=config['fast5_dir']
        bed_file=config['bed_file']
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
                                               suffix2=suffix2[key],
                                               align_config=config['align'], HMM_config=config['HMM']),
                               repeatDetection(model_file=model_file,
                                               repeat=Bio.Seq.reverse_complement(value), 
                                               prefix=Bio.Seq.reverse_complement(suffix[key]),
                                               suffix=Bio.Seq.reverse_complement(prefix[key]),
                                               prefix2=Bio.Seq.reverse_complement(suffix2[key]),
                                               suffix2=Bio.Seq.reverse_complement(prefix2[key]),
                                               align_config=config['align'], HMM_config=config['HMM'])
                               )

        # analyze each read
        for record_ID in record_IDs:
            f5_record = f5.getRecord(record_ID)
            if f5_record is None:
                continue
            with counter.get_lock():
                    counter.value += 1
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
                n, score_prefix, score_suffix, log_p, ticks, offset = model.detect(f5_record.raw, bed_record[3], record_ID)
                # if n < 3:
                    # continue
                output_queue.put('\t'.join([f5_record.ID, bed_record[0], bed_record[3], str(n), 
                                     str(score_prefix), str(score_suffix), str(log_p), str(ticks), str(offset)]))
            else:
                for key, value in repeat.items():
                    model = detection[key][0]
                    n0, score_prefix0, score_suffix0, log_p0, ticks0, offset0 = model.detect(f5_record.raw, bed_record[3])
                    model = detection[key][1]
                    n1, score_prefix1, score_suffix1, log_p1, ticks1, offset1 = model.detect(f5_record.raw, bed_record[3])
                    score0 = score_prefix0 + score_suffix0
                    score1 = score_prefix1 + score_suffix1
                    if score0 > score1 and ticks0 > 0:
                        if n0 < 3:
                            continue
                        output_queue.put('\t'.join([f5_record.ID, key, '+', str(n0), 
                                         str(score_prefix0), str(score_suffix0), str(log_p0), str(ticks0), str(offset0)]))
                    elif score1 > score0 and ticks1 > 0:
                        if n1 < 3:
                            continue
                        output_queue.put('\t'.join([f5_record.ID, key, '-', str(n1), 
                                         str(score_prefix1), str(score_suffix1), str(log_p1), str(ticks1), str(offset1)])) 

    except KeyboardInterrupt:
        return

        

        
# parse config.json            
def parse_config(f5_dir, model_file, repeat_config_file, bed_file, config_file=None):
    with open(repeat_config_file) as fp:
        ld_conf = json.load(fp)
    config = {}
    # parse repeat config
    try:
        config['model_file'] = model_file
        config['fast5_dir'] = f5_dir
        if bed_file and os.path.isfile(bed_file):
            config['bed_file'] = bed_file
        else:
            config['bed_file'] = None
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
        print('Error loading repeat config file, missing', e.args[0])
        exit(1)
    # parse HMM and alignment config
    if config_file:
        with open(config_file) as fp:
            ld_conf = json.load(fp)
        try:
            assert(isinstance(ld_conf, dict))
            assert(isinstance(ld_conf['align'], dict))
            assert(isinstance(ld_conf['HMM'], dict))
            # Do not check values, missing ones get defaulted, additional ones ignored
            config['align'] = ld_conf['align']
            config['HMM'] = ld_conf['HMM']
        except KeyError as e:
            print('Error loading HMM config file, missing', e.args[0])
            exit(1)
        except AssertionError as e:
            print('Config file format broken', b)
            exit(1)
    return config
        
 
 
 
# print progess on command line
def process_output(counter, return_values, output_file, max_count, stop_event):
    import time
    # write header line
    with open(output_file, 'w') as fp:
        print('\t'.join(['ID', 'ref', 'flag', 'count', 'score_prefix', 'score_suffix', 'log_p', 'ticks', 'offset']), file=fp)
        # read result queue until stop event is set
        try:
            while not stop_event.is_set():
                try:
                    record = return_values.get(True, 1)
                    print(record, file=fp)
                except queue.Empty:
                    pass
                print("\rProcessed ", counter.value, 'of ', (max_count), end=" ")
        except KeyboardInterrupt:
            return
        except Exception as e:
            print('Exception in I/O Process: ', e)
        
 

 
if __name__ == '__main__':
    # command line
    parser = argparse.ArgumentParser(description="STR Detection in raw nanopore data")
    parser.add_argument("f5", help="fast5 file directory") 
    parser.add_argument("model", help="pore model")
    parser.add_argument("config", help="job config file")
    parser.add_argument("out", help="output file name")
    parser.add_argument("--bed", help="alignment information in BED-6 format")
    parser.add_argument("--t", type=int, default=1, help="Number of processes to use in parallel")
    parser.add_argument("--hmm", help="Config file for HMM transition probabilities")
    args = parser.parse_args()
    # validate input arguments
    if not os.path.isdir(args.f5):
        print("Fast5 file path is not a directory.")
        exit(1)
    if not os.path.isfile(args.model):
        print("Model file not found.")
        exit(1)
    if not os.path.isfile(args.config):
        print("Repeat config file " + str(args.config) + " not found.")
        exit(1)
    if args.bed and not os.path.isfile(args.bed):
        print("Bed file specified but not found.")
        exit(1)
    if args.hmm and not os.path.isfile(args.hmm):
        print("HMM config file specified but not found.")
        exit(1)
    # load config
    config = parse_config(args.f5, args.model, args.config, args.bed, args.hmm)
    # index/load reads
    f5 = fast5Index.fast5Index()
    f5.loadOrAnalyze(config['fast5_dir'], recursive=True)
    records = np.array(f5.getRecordNames())
    # split in chunks for multi-processing
    record_chunks = np.array_split(records, args.t)
    counter = Value('i', 0)
    return_values = Queue()
    stop_event = Event()
    monitor = Process(target=process_output, args=(counter, return_values, args.out, len(records), stop_event))
    monitor.start()
    worker = []
    # fork into multiple worker processes
    try:
        if len(record_chunks) > 1:
            for i, chunk in enumerate(record_chunks[1:]):
                worker.append(Process(target=process_detection, args=(config, chunk, return_values, counter)))
            for w in worker:
               w.start()
        process_detection(config, record_chunks[0], return_values, counter)
        for w in worker:
           w.join()
    except KeyboardInterrupt:
        pass
    stop_event.set()
    try:
        monitor.join()
    except RuntimeError:
        print('Could not join monitor thread')
    
    print("")
    
    