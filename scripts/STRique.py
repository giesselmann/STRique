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
import re
import signal
import queue
import itertools
import numpy as np
import numpy.ma as ma
import scipy.signal as sp
import pomegranate as pg
import matplotlib.pyplot as plt
from signal import signal, SIGPIPE, SIG_DFL
from collections import namedtuple, defaultdict
from skimage.morphology import opening, closing, dilation, erosion, rectangle
from multiprocessing import Pool, Process, Event, Value, Queue


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
        
    def __sliding_window__(a, n=3, mode='same'):
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
            #sliding_std = [np.std(x) for x in self.__sliding_window__(signal, n=1000, mode='mirror')]
            sliding_std = [self.MAD(x) for x in self.__sliding_window__(signal, n=500, mode='mirror')]
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
    def __init__(self, repeat, pm, transition_probs={}, state_prefix='', std_scale=1.0, std_offset=0.0):
        super().__init__()
        self.repeat = repeat
        self.pore_model = pm
        self.transition_probs = {'skip': .999,
                                 'leave_repeat': .002
                                 }
        for key, value in transition_probs.items():
            self.transition_probs[key] = value
        self.state_prefix = state_prefix
        self.std_scale = std_scale
        self.std_offset = std_offset
        self.__build_model__()

    def __build_model__(self):
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
class flankedRepeatHMM(pg.HiddenMarkovModel):
    def __init__(self, repeat, 
                 prefix, suffix,
                 pm, config=None):
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
        self.pore_model = pm
        self.repeat = repeat
        self.prefix = prefix
        self.suffix = suffix
        self.__build_model__()
        
    def free_bake_buffers(self):
        super().free_bake_buffers()

    def __build_model__(self): 
        # expand primer sequences and get profile HMMs
        prefix = self.prefix + ''.join([self.repeat] * int(np.ceil(self.pore_model.kmer / len(self.repeat))))[:-1]
        suffix = ''.join([self.repeat] * int(np.ceil(self.pore_model.kmer / len(self.repeat)))) + self.suffix
        self.flanking_count = int(np.ceil(self.pore_model.kmer / len(self.repeat))) * 2 - 1
        self.prefix_model = profileHMM(prefix, self.pore_model, self.transition_probs, state_prefix='prefix', std_scale=self.transition_probs['seq_std_scale'], std_offset=self.transition_probs['seq_std_offset'])
        self.suffix_model = profileHMM(suffix, self.pore_model, self.transition_probs, state_prefix='suffix', std_scale=self.transition_probs['rep_std_scale'], std_offset=self.transition_probs['seq_std_offset'])
        self.repeat_model = repeatHMM(self.repeat, self.pore_model, self.transition_probs, state_prefix='repeat', std_scale=self.transition_probs['seq_std_scale'], std_offset=self.transition_probs['rep_std_offset'])
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
class repeatCounter(object):
    def __init__(self, model_file, align_config=None, HMM_config=None):
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
        self.samples = default_config['samples']
        self.HMM_config = HMM_config
        self.targets = {}
        self.target_classifier = namedtuple('target_classifier', field_names=['prefix', 'suffix', 'prefix_ext', 'suffix_ext', 'repeatHMM'])
        
    def __reverse_complement__(self, sequence):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return "".join(complement.get(base, base) for base in reversed(sequence))

    def __detect_range__(self, signal, segment, pre_trim=0, post_trim=0):
        score, idx_signal, idx_segment = self.algn.align_overlap(signal, segment)
        segment_begin = np.abs(np.array(idx_signal) - idx_segment[0]).argmin()
        segment_end = np.abs(np.array(idx_signal) - idx_segment[-1]).argmin()
        score = score / (segment_end - segment_begin)
        segment_begin = np.abs(np.array(idx_signal) - idx_segment[0 + pre_trim]).argmin()
        segment_end = np.abs(np.array(idx_signal) - idx_segment[-1 - post_trim]).argmin()
        return score, segment_begin, segment_end

    def __detect_short__(self, flanked_model, segment):
        return flanked_model.count_repeats(segment)
        
    def add_target(self, target_name, repeat, prefix, suffix):
        if not target_name in self.targets:
            prefix_ext = prefix
            prefix = prefix[:50]
            suffix_ext = suffix
            suffix = suffix[:50]
            # template model
            tc_plus = self.target_classifier(
                self.pm.generate_signal(prefix, samples=self.samples),
                self.pm.generate_signal(suffix, samples=self.samples),
                self.pm.generate_signal(prefix_ext, samples=self.samples),
                self.pm.generate_signal(suffix_ext, samples=self.samples),
                flankedRepeatHMM(repeat, prefix, suffix, self.pm, self.HMM_config) )
            # complement model
            tc_minus = self.target_classifier(
                self.pm.generate_signal(self.__reverse_complement__(suffix), samples=self.samples),
                self.pm.generate_signal(self.__reverse_complement__(prefix), samples=self.samples),
                self.pm.generate_signal(self.__reverse_complement__(suffix_ext), samples=self.samples),
                self.pm.generate_signal(self.__reverse_complement__(prefix_ext), samples=self.samples),
                flankedRepeatHMM(self.__reverse_complement__(repeat), self.__reverse_complement__(suffix), self.__reverse_complement__(prefix), self.pm, self.HMM_config))
            self.targets[target_name] = (tc_plus, tc_minus)
        else:
            raise ValueError("[repeatCounter] Target with name " + str(target_name) + " already defined")
            
    def detect(self, target_name, raw_signal, strand):
        if target_name in self.targets:
            tc_plus, tc_minus = self.targets[target_name]
            if strand == '+':
                tc = tc_plus
            elif strand == '-':
                tc = tc_minus
            else:
                ValueError("[repeatCounter] Strand must be + or -")
            flt_signal = sp.medfilt(raw_signal, kernel_size=3)
            nrm_signal = (flt_signal - np.median(flt_signal)) / self.pm.MAD(flt_signal)
            nrm_signal = np.clip(nrm_signal * 24 + 127, 0, 255).astype(np.dtype('uint8')).reshape((1, len(nrm_signal)))
            flt = rectangle(1, 8)
            nrm_signal = opening(nrm_signal, flt)
            nrm_signal = closing(nrm_signal, flt)[0].astype(np.dtype('float'))
            nrm_signal = self.pm.normalize2model(nrm_signal.astype(np.dtype('float')), mode='minmax')
            flt_signal = self.pm.normalize2model(flt_signal.astype(np.dtype('float')), mode='minmax')
            trim_prefix = len(tc.prefix_ext) - len(tc.prefix)
            trim_suffix = len(tc.suffix_ext) - len(tc.suffix)
            score_prefix, prefix_begin, prefix_end = self.__detect_range__(nrm_signal, tc.prefix_ext, pre_trim=trim_prefix)
            score_suffix, suffix_begin, suffix_end = self.__detect_range__(nrm_signal, tc.suffix_ext, post_trim=trim_suffix)
            n = 0; p = 0; path = np.array([])
            if prefix_end < suffix_begin:
                n, p, path = self.__detect_short__(tc.repeatHMM, flt_signal[prefix_begin:suffix_end])
            return n, score_prefix, score_suffix, p, prefix_end, suffix_begin - prefix_end
        else:
            raise ValueError("[repeatCounter] Target with name " + str(target_name) + " not defined")




# multi locus repeat detection 
class repeatDetector(object):
    class sam_record(object):
        def __init__(self):
            self.QNAME = ''
            self.FLAG = 0
            self.RNAME = ''
            self.POS = 0
            self.TLEN = 0
            self.CLIP_BEGIN = 0
            self.CLIP_END = 0
            
    def __init__(self, repeat_config, model_file, fast5_dir, align_config=None, HMM_config=None):
        self.repeatCounter = repeatCounter(model_file, align_config, HMM_config)
        self.repeatLoci = defaultdict(lambda : [])
        self.repeat_config = repeat_config
        self.is_init = False
        self.f5 = fast5Index.fast5Index()
        self.f5.load(fast5_dir)
        
    def __init_hmm__(self):
        for target_name, (chr, begin, end, repeat, prefix, suffix) in self.repeat_config.items():
            self.repeatCounter.add_target(target_name, repeat, prefix, suffix)
            self.repeatLoci[chr].append((target_name, begin, end))
        self.is_init = True
        
    def __decode_cigar__(self, cigar):
        ops = [(int(op[:-1]), op[-1]) for op in re.findall('(\d*\D)',cigar)]
        return ops

    def __ops_length__(self, ops, recOps='MIS=X'):
        n = [op[0] for op in ops if op[1] in recOps]
        return sum(n)
        
    def __decode_sam__(self, sam_line):
        cols = sam_line.rstrip().split('\t')
        sr = self.sam_record()
        if len(cols) >= 11:
            try:
                sr.QNAME = cols[0]
                sr.FLAG = int(cols[1])
                sr.RNAME = cols[2]
                sr.POS = int(cols[3])
                cigar_ops = self.__decode_cigar__(cols[5])
                sr.TLEN = self.__ops_length__(cigar_ops, recOps='MDN=X')
                sr.CLIP_BEGIN = sum([op[0] for op in cigar_ops[:2] if op[1] in 'SH'])
                sr.CLIP_END = sum([op[0] for op in cigar_ops[-2:] if op[1] in 'SH'])
            except:
                return self.sam_record()
        return sr
        
    def __intersect_target__(self, sam_record):
        target_names = []
        if sam_record.RNAME in self.repeatLoci.keys():
            for target_name, begin, end in self.repeatLoci[sam_record.RNAME]:
                if begin > sam_record.POS - sam_record.CLIP_BEGIN and end < sam_record.POS + sam_record.TLEN + sam_record.CLIP_END:
                    target_names.append(target_name)
        return target_names

    def detect(self, sam_line=''):
        if not self.is_init:
            self.__init_hmm__()
        sam_record = self.__decode_sam__(sam_line)
        f5_record = self.f5.getRecord(sam_record.QNAME)
        target_counts = []
        if sam_record and f5_record:
            if sam_record.FLAG & 0x10 == 0:
                strand = '+'
            else:
                strand = '-'
            target_names = self.__intersect_target__(sam_record)
            for target_name in target_names:
                repeat_count = self.repeatCounter.detect(target_name, f5_record.raw, strand)
                target_counts.append((sam_record.QNAME, target_name, strand, *repeat_count))
        return {'target_counts': target_counts}




# writes repeat detection output to file or stdout
class outputWriter(object):
    def __init__(self, output_file=None):
        self.output_file = output_file
        if self.output_file:
            with open(self.output_file, 'w') as fp:
                print('\t'.join(['ID', 'target', 'strand', 'count', 'score_prefix', 'score_suffix', 'log_p', 'offset', 'ticks']), file=fp)
        else:
            print('\t'.join(['ID', 'target', 'strand', 'count', 'score_prefix', 'score_suffix', 'log_p', 'offset', 'ticks']))

    def write_line(self, target_counts=[]):
        if self.output_file:
            with open(self.output_file, 'a') as fp:
                for target_count in target_counts:
                    print('\t'.join([str(x) for x in target_count]), file=fp)
        else:
            for target_count in target_counts:
                print('\t'.join([str(x) for x in target_count]))




# multiprocess dispatcher
class mt_dispatcher():
    def __init__(self, input_queue, threads=1, worker_callables=[], collector_callables=[]):
        self.worker_callables = worker_callables
        self.collector_callables = collector_callables
        self.n_processed = 0
        self.n_worker = threads
        self.input_queue = input_queue
        self.output_queue = Queue()
        self.collector_queue = Queue(threads * 10)
        self.worker = []
        for i in range(threads):
            self.worker.append(Process(target=self.__worker__, ))
            self.worker[-1].start()
        self.collector = Process(target=self.__collector__, )
        self.collector.start()
        
    def __worker__(self):
        try:
            while True:
                input = self.input_queue.get()
                if not input:
                    break
                try:
                    for worker_callable in self.worker_callables:
                        input = worker_callable(**input)
                    self.collector_queue.put(input)
                    self.collector_queue.put('done')
                except KeyboardInterrupt:
                    break
                # except:
                    # continue
        except KeyboardInterrupt:
            self.collector_queue.put(None)
            return
        self.collector_queue.put(None)
        self.collector_queue.close()
        
    def __collector__(self):
        poison_count = self.n_worker
        try:
            while True:
                input = self.collector_queue.get()
                if input is None:
                    poison_count -= 1
                    if poison_count <= 0:
                        break
                    continue
                elif input == 'done':
                    self.n_processed += 1
                    continue
                for collector_callable in self.collector_callables:
                    input = collector_callable(**input)
                self.output_queue.put(input)
        except KeyboardInterrupt:
            pass
        # except:
            # pass
        self.output_queue.put(None)
        self.output_queue.close()
            
    def __stop__(self):
        for w in self.worker:
            self.input_queue.put(None)
        self.input_queue.close()
        for w in self.worker:
            w.join()
        self.collector_queue.close()
        self.collector.join()
        
    def __iter__(self):
        return self
        
    def __next__(self):
        while True:
            result = self.output_queue.get()
            if result is None:
                self.output_queue.close()
                raise StopIteration()
            else:
                return result
                
    def close(self):
        self.__stop__()
    
    def n_processed(self):
        return self.n_processed




 # parse config.json            
def parse_config(repeat_config_file, param_config_file=None):
    config = {}
    # parse repeat config
    repeats = {}
    with open(repeat_config_file, 'r') as fp:
        header = next(fp)
        for line in fp:
            cols = line.rstrip().split('\t')
            if len(cols) == 7:
                repeats[cols[3]] = (cols[0], int(cols[1]), int(cols[2]), cols[4], cols[5], cols[6])
            else:
                print("[Config] Repeat config column mismatch while parsing", file=sys.stderr)
                print("[Config] " + line, file=sys.stderr)
    config['repeat'] = repeats
    config['align'] = None
    config['HMM'] = None
    # parse HMM and alignment config
    if param_config_file:
        with open(param_config_file) as fp:
            ld_conf = json.load(fp)
        try:
            assert(isinstance(ld_conf, dict))
            assert(isinstance(ld_conf['align'], dict))
            assert(isinstance(ld_conf['HMM'], dict))
            # Do not check values, missing ones get defaulted, additional ones ignored
            config['align'] = ld_conf['align']
            config['HMM'] = ld_conf['HMM']
        except KeyError as e:
            print('[Config] Error loading HMM config file, missing', e.args[0], file=sys.stderr)
            exit(1)
        except AssertionError as e:
            print('[Config] Config file format broken', file=sys.stderr)
            exit(1)
    return config




# main
if __name__ == '__main__':
    signal(SIGPIPE,SIG_DFL)
    # command line
    parser = argparse.ArgumentParser(description="STR Detection in raw nanopore data")
    parser.add_argument("f5", help="fast5 file directory")
    parser.add_argument("model", help="pore model")
    parser.add_argument("repeat", help="repeat region config file")
    parser.add_argument("--out", default=None, help="output file name, if not given print to stdout")
    parser.add_argument("--algn", default=None, help="alignment in sam format, if not given read from stdin")
    parser.add_argument("--config", help="Config file with HMM transition probabilities")
    parser.add_argument("--t", type=int, default=1, help="Number of processes to use in parallel")
    args = parser.parse_args()
    # load config
    config = parse_config(args.repeat, args.config)
    # index/load reads
    f5 = fast5Index.fast5Index()
    f5.loadOrAnalyze(args.f5, recursive=True)
    # repeat detector
    rd = repeatDetector(config['repeat'], args.model, args.f5, align_config=config['align'], HMM_config=config['HMM'])
    ow = outputWriter(args.out)
    # run repeat detection
    if args.t > 1:
        sam_queue = Queue(100)
        mt = mt_dispatcher(sam_queue, threads=args.t, worker_callables=[rd.detect], collector_callables=[ow.write_line])
        if args.algn:
            with open(args.algn, 'r') as fp:
                for line in fp:
                    sam_queue.put({'sam_line': line})
        else:
            for line in sys.stdin:
                sam_queue.put({'sam_line': line})
        mt.close()
    else:
        if args.algn:
            with open(args.algn, 'r') as fp:
                for line in fp:
                    counts = rd.detect({'sam_line': line})
                    ow.write_line(**counts)
        else:
            for line in sys.stdin:
                counts = rd.detect({'sam_line': line})
                ow.write_line(**counts)                
        
    
    