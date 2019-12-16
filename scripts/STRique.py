# \MODULE\-------------------------------------------------------------------------
#
#  CONTENTS      : STRique
#
#  DESCRIPTION   : Raw nanopore signal repeat detection pipeline
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
import os, sys, traceback, json, argparse
import re, itertools
import datetime
import threading, queue
import enum
import numpy as np
import numpy.ma as ma
import scipy.signal as sp
import pomegranate as pg
from signal import signal, SIGPIPE, SIG_DFL
from collections import namedtuple, defaultdict
from skimage.morphology import opening, closing, dilation, erosion, rectangle
from multiprocessing import Pool, Process, Event, Value, Queue
# private imports
from STRique_lib import fast5Index, pyseqan




# simple parallel logging
class logger():
    logs = [sys.stderr]
    class log_type(enum.Enum):
        Error = "[ERROR]"
        Warning = "[WARNING]"
        Info = "[INFO]"
        Debug = "[DEBUG]"
    log_types = []
    lock = threading.Lock()
    log_queue = Queue()

    def __logger__():
        while True:
            print_message = logger.log_queue.get()
            if not print_message:
                break
            for log in logger.logs:
                if isinstance(log, str):
                    with open(log, 'a') as fp:
                        print(print_message, file = fp)
                else:
                    print(print_message, file = log)
                    sys.stderr.flush()

    def init(file=None, log_level='info'):
        if log_level == 'error':
            logger.log_types = [logger.log_type.Error]
        elif log_level == 'warning':
            logger.log_types = [logger.log_type.Error, logger.log_type.Warning]
        elif log_level == 'info':
            logger.log_types = [logger.log_type.Error, logger.log_type.Warning, logger.log_type.Info]
        else:
            logger.log_types = [logger.log_type.Error, logger.log_type.Warning, logger.log_type.Info, logger.log_type.Debug]
        if file:
            if os.path.isfile(file) and os.access(file, os.W_OK) or os.access(os.path.abspath(os.path.dirname(file)), os.W_OK):
                logger.logs.append(file)

        logger.log_runner = Process(target=logger.__logger__, )
        logger.log_runner.start()
        logger.log("Logger created.")
        if file and len(logger.logs) == 1:
            logger.log("Log-file {file} is not accessible".format(file=file), logger.log_type.Error)

    def close():
        logger.log_queue.put(None)
        logger.log_queue.close()
        logger.log_runner.join()

    def log(message, type=log_type.Info):
        with logger.lock:
            if type in logger.log_types:
                print_message = ' '.join([datetime.datetime.now().strftime("%d.%m.%Y %H:%M:%S"), "[PID {}]".format(os.getpid()), str(type.value), message])
                logger.log_queue.put(print_message)




# basic normalization and simulation
class pore_model():
    def __init__(self, model_file):
        def model_iter(iterable):
            for line in iterable:
                yield line.strip().split('\t')[:3]
        with open(model_file, 'r') as fp:
            model_dict = {x[0]:(float(x[1]), float(x[2])) for x in model_iter(fp)}
        self.kmer = len(next(iter(model_dict.keys())))
        self.model_median = np.median([x[0] for x in model_dict.values()])
        self.model_MAD = np.mean(np.absolute(np.subtract([x[0] for x in model_dict.values()], self.model_median)))
        min_state = min(model_dict.values(), key=lambda x:x[0])
        max_state = max(model_dict.values(), key=lambda x:x[0])
        self.model_min = min_state[0] - 6 * min_state[1]
        self.model_max = max_state[0] + 6 * max_state[1]
        self.model_dict = model_dict

    def __sliding_window__(self, a, n=3, mode='same'):
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

    def scale2stdv(self, other):
        self_median = np.median(np.array([x[1] for x in self.model_dict.values()]))
        other_median = np.median(np.array([x[1] for x in other.model_dict.values()]))
        return other_median / self_median

    def normalize2model(self, signal, clip=True, mode='median'):
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
            sliding_std = [self.MAD(x) for x in self.__sliding_window__(signal, n=500, mode='mirror')]
            sliding_std += [sliding_std[-1]]
            diff_signal = np.abs(np.diff(sliding_std))
            ind = np.argpartition(diff_signal, -50)[-50:]
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
        return nrm_signal

    def generate_signal(self, sequence, samples=10, noise=False):
        signal = []
        level_means = np.array([self.model_dict[kmer][0] for kmer in [sequence[i:i+self.kmer] for i in range(len(sequence)-self.kmer + 1)]])
        if samples and not noise:
            sig = np.repeat(level_means, samples)
        elif not noise:
            sig = np.repeat(level_means, np.random.uniform(6, 10, len(level_means)).astype(int))
        else:
            level_stdvs = np.array([self.model_dict[kmer][1] for kmer in [sequence[i:i+self.kmer] for i in range(len(sequence)-self.kmer + 1)]])
            level_samples = np.random.uniform(6, 10, len(level_means)).astype(int)
            level_means = np.repeat(level_means, level_samples)
            level_stdvs = np.repeat(level_stdvs, level_samples)
            sig = np.random.normal(level_means, level_stdvs)
        return sig




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
        self.transition_probs = {'match_loop': .75,     # .75
                                 'match_match': .15,    # .15          sum to 1
                                 'match_insert': .09,   # .09
                                 'match_delete': .01,   # .01

                                 'insert_loop' : .15,   # .15
                                 'insert_match_0': .40, # .40          sum to 1
                                 'insert_match_1': .40, # .40
                                 'insert_delete': .05,  # .05

                                 'delete_delete': .005, # .005
                                 'delete_insert': .05,  # .05          sum to 1
                                 'delete_match': .945   # .945
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

#    def free_bake_buffers(self, *args, **kwargs):
#        print("Free bake buffers in process {id}".format(id=os.getpid()))
#        super().free_bake_buffers()

    def __build_model__(self):
        # expand primer sequences and get profile HMMs
        prefix = self.prefix + ''.join([self.repeat] * int(np.ceil(self.pore_model.kmer / len(self.repeat))))[:-1]
        suffix = ''.join([self.repeat] * int(np.ceil(self.pore_model.kmer / len(self.repeat)))) + self.suffix
        self.flanking_count = int(np.ceil(self.pore_model.kmer / len(self.repeat))) * 2 - 1
        self.prefix_model = profileHMM(prefix, self.pore_model, self.transition_probs, state_prefix='prefix', std_scale=self.transition_probs['seq_std_scale'], std_offset=self.transition_probs['seq_std_offset'])
        self.suffix_model = profileHMM(suffix, self.pore_model, self.transition_probs, state_prefix='suffix', std_scale=self.transition_probs['seq_std_scale'], std_offset=self.transition_probs['seq_std_offset'])
        self.repeat_model = repeatHMM(self.repeat, self.pore_model, self.transition_probs, state_prefix='repeat', std_scale=self.transition_probs['rep_std_scale'], std_offset=self.transition_probs['rep_std_offset'])
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
            path = [x[1].name for x in path if x[0] < self.silent_start]
            return n, p, path
        else:
            return 0, 0, []




# repeat base modification detection
class repeatModHMM(pg.HiddenMarkovModel):
    def __init__(self, repeat, pm_base, pm_mod, config=None):
        super().__init__()
        self.transition_probs = {'rep_std_scale': 1.5,
                                 'rep_std_offset': 0.0,
                                 'leave_repeat': .002}
        if config and isinstance(config, dict):
            for key, value in config.items():
                self.transition_probs[key] = value
        self.pore_model_base = pm_base
        self.pore_model_mod = pm_mod
        self.repeat = repeat
        self.__build_model__()

    def __build_model__(self):
        if len(self.repeat) >= self.pore_model_base.kmer:
            repeat = self.repeat + self.repeat[: self.pore_model_base.kmer - 1]
        else:
            ext = self.pore_model_base.kmer - 1 + (len(self.repeat) - 1) - ((self.pore_model_base.kmer - 1) % len(self.repeat))
            repeat = self.repeat + ''.join([self.repeat] * self.pore_model_base.kmer)[:ext]
        self.model_min = min(self.pore_model_base.model_min, self.pore_model_mod.model_min)
        self.model_max = max(self.pore_model_base.model_max, self.pore_model_mod.model_max)
        self.s0 = pg.State(pg.UniformDistribution(self.model_min, self.model_max), name='s0')
        self.e0 = pg.State(pg.UniformDistribution(self.model_min, self.model_max), name='e0')
        self.base_model = profileHMM(repeat, self.pore_model_base, self.transition_probs, state_prefix='base', no_silent=True, std_scale=self.transition_probs['rep_std_scale'], std_offset=self.transition_probs['rep_std_offset'])
        self.mod_model = profileHMM(repeat, self.pore_model_mod, self.transition_probs, state_prefix='mod', no_silent=True, std_scale=self.transition_probs['rep_std_scale'] * self.pore_model_mod.scale2stdv(self.pore_model_base), std_offset=self.transition_probs['rep_std_offset'])
        self.add_model(self.base_model)
        self.add_model(self.mod_model)
        self.add_state(self.s0)
        self.add_state(self.e0)
        # transitions
        self.add_transition(self.start, self.s0, 1)
        self.add_transition(self.s0, self.base_model.s1, 0.25)
        self.add_transition(self.s0, self.base_model.s2, 0.25)
        self.add_transition(self.s0, self.mod_model.s1, 0.25)
        self.add_transition(self.s0, self.mod_model.s2, 0.25)
        self.add_transition(self.base_model.e1, self.e0, 1)
        self.add_transition(self.base_model.e2, self.e0, 1)
        self.add_transition(self.mod_model.e1, self.e0, 1)
        self.add_transition(self.mod_model.e2, self.e0, 1)
        self.add_transition(self.e0, self.end, self.transition_probs['leave_repeat'])
        self.add_transition(self.e0, self.s0, 1-self.transition_probs['leave_repeat'])
        # bake
        self.bake(merge='All')

    def mod_repeats(self, signal, **kwargs):
        p, path = super().viterbi(np.clip(signal, self.model_min, self.model_max), **kwargs)
        if path is not None:
            states = [x[1].name for x in path if x[0] < self.silent_start]
            states_gr = [next(g) for k,g in itertools.groupby(states, key=lambda x : False if x in ['s0', 'e0'] else True) if k]
            pattern = ''.join(['1' if 'mod' in x else '0' for x in states_gr])
            return pattern
        else:
            return '-'



# main repeat detection methods
class repeatCounter(object):
    def __init__(self, model_file, mod_model_file=None, align_config=None, HMM_config=None):
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
        if mod_model_file:
            self.pm_mod = pore_model(mod_model_file)
        else:
            self.pm_mod = self.pm
        self.samples = default_config['samples']
        self.HMM_config = HMM_config
        self.targets = {}
        self.target_classifier = namedtuple('target_classifier', field_names=['prefix', 'suffix', 'prefix_ext', 'suffix_ext', 'repeatHMM', 'modHMM'])

    def __reverse_complement__(self, sequence):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return "".join(complement.get(base, base) for base in reversed(sequence))

    def __detect_range__(self, signal, segment, pre_trim=0, post_trim=0):
        score, idx_signal, idx_segment = self.algn.align_overlap(signal, segment)
        segment_begin = np.abs(np.array(idx_signal) - idx_segment[0]).argmin()
        segment_end = np.abs(np.array(idx_signal) - idx_segment[-1]).argmin()
        if segment_end > segment_begin:
            score = score / (segment_end - segment_begin)
        else:
            score = 0.0
        segment_begin = np.abs(np.array(idx_signal) - idx_segment[0 + pre_trim]).argmin()
        segment_end = np.abs(np.array(idx_signal) - idx_segment[-1 - post_trim]).argmin()
        return score, segment_begin, segment_end

    def __detect_short__(self, flanked_model, segment):
        return flanked_model.count_repeats(segment)

    def add_target(self, target_name, repeat, prefix, suffix):
        if not target_name in self.targets:
            prefix_ext = prefix.upper()
            prefix = prefix[-50:].upper()
            suffix_ext = suffix.upper()
            suffix = suffix[:50].upper()
            repeat = repeat.upper()
            # template model
            tc_plus = self.target_classifier(
                self.pm.generate_signal(prefix, samples=self.samples),
                self.pm.generate_signal(suffix, samples=self.samples),
                self.pm.generate_signal(prefix_ext, samples=self.samples),
                self.pm.generate_signal(suffix_ext, samples=self.samples),
                flankedRepeatHMM(repeat, prefix, suffix, self.pm, self.HMM_config),
                repeatModHMM(repeat, self.pm, self.pm_mod, config=self.HMM_config))
            # complement model
            tc_minus = self.target_classifier(
                self.pm.generate_signal(self.__reverse_complement__(suffix), samples=self.samples),
                self.pm.generate_signal(self.__reverse_complement__(prefix), samples=self.samples),
                self.pm.generate_signal(self.__reverse_complement__(suffix_ext), samples=self.samples),
                self.pm.generate_signal(self.__reverse_complement__(prefix_ext), samples=self.samples),
                flankedRepeatHMM(self.__reverse_complement__(repeat), self.__reverse_complement__(suffix), self.__reverse_complement__(prefix), self.pm, self.HMM_config),
                repeatModHMM(self.__reverse_complement__(repeat), self.pm, self.pm_mod, config=self.HMM_config))
            self.targets[target_name] = (tc_plus, tc_minus)
            logger.log("RepeatCounter: Added target {}".format(target_name), logger.log_type.Info)
        else:
            raise ValueError("RepeatCounter: Target with name " + str(target_name) + " already defined.")

    def detect(self, target_name, raw_signal, strand):
        if target_name in self.targets:
            tc_plus, tc_minus = self.targets[target_name]
            if strand == '+':
                tc = tc_plus
            elif strand == '-':
                tc = tc_minus
            else:
                raise ValueError("RepeatCounter: Strand must be + or -.")
            flt_signal = sp.medfilt(raw_signal, kernel_size=3)
            morph_signal = (flt_signal - np.median(flt_signal)) / self.pm.MAD(flt_signal)
            morph_signal = np.clip(morph_signal * 24 + 127, 0, 255).astype(np.dtype('uint8')).reshape((1, len(morph_signal)))
            flt = rectangle(1, 8)
            morph_signal = opening(morph_signal, flt)
            morph_signal = closing(morph_signal, flt)[0].astype(np.dtype('float'))
            morph_signal = self.pm.normalize2model(morph_signal.astype(np.dtype('float')), mode='minmax')
            flt_signal = self.pm.normalize2model(flt_signal.astype(np.dtype('float')), mode='minmax')
            trim_prefix = len(tc.prefix_ext) - len(tc.prefix)
            trim_suffix = len(tc.suffix_ext) - len(tc.suffix)
            score_prefix, prefix_begin, prefix_end = self.__detect_range__(morph_signal, tc.prefix_ext, pre_trim=trim_prefix)
            score_suffix, suffix_begin, suffix_end = self.__detect_range__(morph_signal, tc.suffix_ext, post_trim=trim_suffix)
            n = 0; p = 0; states = []; mod_pattern = '-'
            if prefix_begin < suffix_end and score_prefix > 0.0 and score_suffix > 0.0:
                n, p, states = self.__detect_short__(tc.repeatHMM, flt_signal[prefix_begin:suffix_end])
                if self.pm != self.pm_mod:
                    #rep_signal = flt_signal[prefix_begin:suffix_end][np.array([True if 'repeat' in x else False for x in states])]
                    nrm_signal = self.pm.normalize2model(raw_signal.astype(np.dtype('float')), mode='minmax')
                    rep_signal = nrm_signal[prefix_begin:suffix_end][np.array([True if 'repeat' in x else False for x in states])]
                    mod_pattern = tc.modHMM.mod_repeats(rep_signal)
                    # import matplotlib.pyplot as plt
                    # f, ax = plt.subplots(2, sharex=True)
                    # ax[0].plot(raw_signal[prefix_begin:suffix_end][np.array([True if 'repeat' in x else False for x in states])], 'k-')
                    # ax[1].plot(rep_signal, 'k-')
                    # ax[0].set_title('Count {count}, strand {strand}'.format(count=n, strand=strand))
                    # plt.show()
            return n, score_prefix, score_suffix, p, prefix_end, max(suffix_begin - prefix_end, 0), mod_pattern
        else:
            raise ValueError("RepeatCounter: Target with name " + str(target_name) + " not defined.")




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

    def __init__(self, repeat_config, model_file, fast5_index_file, mod_model_file=None, align_config=None, HMM_config=None):
        self.repeatCounter = repeatCounter(model_file, mod_model_file=mod_model_file, align_config=align_config, HMM_config=HMM_config)
        self.repeatLoci = defaultdict(lambda : [])
        self.repeat_config = repeat_config
        self.is_init = False
        self.f5 = fast5Index.fast5Index(fast5_index_file)

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
        target_counts = []
        sam_record = self.__decode_sam__(sam_line)
        if not sam_record.QNAME:
            logger.log("Detector: Error parsing alignment \n{}".format(sam_line), logger.log_type.Error)
            return None
        if sam_record.FLAG & 0x10 == 0:
            strand = '+'
        else:
            strand = '-'
        target_names = self.__intersect_target__(sam_record)
        if not target_names:
            logger.log("Detector: No target for {}".format(sam_record.QNAME), logger.log_type.Debug)
            return None
        f5_record = self.f5.get_raw(sam_record.QNAME)
        if f5_record is None:
            logger.log("Detector: No fast5 for ID {id}".format(id=sam_record.QNAME), logger.log_type.Warning)
            return None
        logger.log("Detector: Test {id} for targets: {targets}.".format(id=sam_record.QNAME, targets=','.join(target_names)), logger.log_type.Debug)
        for target_name in target_names:
            repeat_count = self.repeatCounter.detect(target_name, f5_record, strand)
            target_counts.append((sam_record.QNAME, target_name, strand, *repeat_count))
        return {'target_counts': target_counts}




# writes repeat detection output to file or stdout
class outputWriter(object):
    def __init__(self, output_file=None):
        self.output_file = output_file
        if self.output_file:
            with open(self.output_file, 'w') as fp:
                print('\t'.join(['ID', 'target', 'strand', 'count', 'score_prefix', 'score_suffix', 'log_p', 'offset', 'ticks', 'mod']), file=fp)
        else:
            print('\t'.join(['ID', 'target', 'strand', 'count', 'score_prefix', 'score_suffix', 'log_p', 'offset', 'ticks', 'mod']))

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
                    if input:
                        self.collector_queue.put(input)
                        self.collector_queue.put('done')
                except KeyboardInterrupt:
                    logger.log("Factory: Worker terminating on user request.", logger.log_type.Info)
                    break
                except Exception as e:
                    type, value, trace = sys.exc_info()
                    logger.log('\n'.join(["Factory: Unexpected error in Worker, proceeding wiht remaining reads."] +
                                            traceback.format_exception(*sys.exc_info())), logger.log_type.Warning)
                    pass
        except KeyboardInterrupt:
            self.collector_queue.put(None)
            self.collector_queue.close()
            logger.log("Factory: Worker terminating on user request.", logger.log_type.Info)
            return
        self.collector_queue.put(None)
        self.collector_queue.close()
        logger.log("Factory: Worker terminating.", logger.log_type.Debug)

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
            logger.log("Factory: Collector terminating on user request.", logger.log_type.Info)
            pass
        except:
            logger.log("Factory: Unexpected error in collector, terminating.", logger.log_type.Error)
            pass
        self.output_queue.put(None)
        self.output_queue.close()
        logger.log("Factory: Collector terminating.", logger.log_type.Debug)

    def __stop__(self):
        for w in self.worker:
            self.input_queue.put(None)
        self.input_queue.close()
        for w in self.worker:
            w.join()
        self.collector_queue.close()
        self.collector.join()
        logger.log("Factory: Terminated.", logger.log_type.Debug)

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
            cols = line.rstrip().split()
            if len(cols) == 7:
                repeats[cols[3]] = (cols[0], int(cols[1]), int(cols[2]), cols[4], cols[5], cols[6])
            else:
                logger.log("Config: Repeat config column mismatch while parsing \n{line}".format(line=line), logger.log_type.Error)
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
            logger.log('Config: Error loading HMM config file, missing {}'.format(e.args[0]), logger.log_type.Error)
            exit(1)
        except AssertionError as e:
            logger.log('Config: file format broken', logger.log_type.Error)
            exit(1)
    return config




# main class
class main():
    def __init__(self):
        parser = argparse.ArgumentParser(
        description='STRique: a nanopore raw signal repeat detection pipeline',
        usage='''STRique.py <command> [<args>]
Available commands are:
   index      Index batch(es) of bulk-fast5 or tar archived single fast5
   count      Count single read repeat expansions
   plot       Plot repeat signal after counting
''')
        parser.add_argument('command', help='Subcommand to run')
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print('Unrecognized command', file=sys.stderr)
            parser.print_help(file=sys.stderr)
            exit(1)
        getattr(self, args.command)(sys.argv[2:])

    def index(self, argv):
        parser = argparse.ArgumentParser(description="Fast5 raw data archive indexing")
        parser.add_argument("input", help="Input batch or directory of batches")
        parser.add_argument("--recursive", action='store_true', help="Recursively scan input")
        parser.add_argument("--out_prefix", default="", help="Prefix for file paths in output")
        parser.add_argument("--tmp_prefix", default=None, help="Prefix for temporary data")
        args = parser.parse_args(argv)
        for record in fast5Index.fast5Index.index(args.input, recursive=args.recursive, output_prefix=args.out_prefix, tmp_prefix=args.tmp_prefix):
            print(record)

    def count(self, argv):
        # command line
        parser = argparse.ArgumentParser(description="STR Detection in raw nanopore data")
        parser.add_argument("f5Index", help="Fast5 index")
        parser.add_argument("model", help="Pore model")
        parser.add_argument("repeat", help="Repeat region config file")
        parser.add_argument("--out", default=None, help="Output file name, if not given print to stdout")
        parser.add_argument("--algn", default=None, help="Alignment in sam format, if not given read from stdin")
        parser.add_argument("--mod_model", default=None, help="Base modification pore model")
        parser.add_argument("--config", help="Config file with HMM transition probabilities")
        parser.add_argument("--t", type=int, default=1, help="Number of processes to use in parallel")
        parser.add_argument("--log_level", default='warning', choices=['error', 'warning', 'info', 'debug'], help="Log level")
        args = parser.parse_args(argv)
        logger.init(log_level=args.log_level)
        # load config
        config = parse_config(args.repeat, args.config)
        logger.log("Main: Parsed config.", logger.log_type.Debug)
        # index/load reads
        if not os.path.isfile(args.f5Index):
            logger.log("Main: Fast5 index file does not exist.", logger.log_type.Error)
            exit(1)
        # model files
        if not os.path.isfile(args.model):
            logger.log("Main: Pore model file does not exist.", logger.log_type.Error)
            exit(1)
        if args.mod_model and not os.path.isfile(args.mod_model):
            logger.log("Main: Modification pore model file does not exist.", logger.log_type.Error)
            exit(1)
        # repeat detector
        rd = repeatDetector(config['repeat'], args.model, args.f5Index, mod_model_file=args.mod_model, align_config=config['align'], HMM_config=config['HMM'])
        ow = outputWriter(args.out)
        # run repeat detection
        sam_queue = Queue(100)
        mt = mt_dispatcher(sam_queue, threads=args.t, worker_callables=[rd.detect], collector_callables=[ow.write_line])
        if args.algn:
            with open(args.algn, 'r') as fp:
                for line in fp:
                    if not line.startswith('@'):
                        sam_queue.put({'sam_line': line})
        else:
            for line in sys.stdin:
                if not line.startswith('@'):
                    sam_queue.put({'sam_line': line})
        mt.close()
        logger.close()

    def plot(self, argv):
        parser = argparse.ArgumentParser(description="Signal plots over STR expansions")
        parser.add_argument("f5Index", help="Fast5 index")
        parser.add_argument("--counts", default=None, help="Repeat count output from STRique, if not given read from stdin")
        parser.add_argument("--output", default=None, help="Output directory for plots, use instead of interactive GUI")
        parser.add_argument("--format", default='png', choices={"png", "pdf", "svg"}, help="Output format when writing to files")
        parser.add_argument("--width", default=16, type=int, help="Plot width")
        parser.add_argument("--height", default=9, type=int, help="Plot height")
        parser.add_argument("--dpi", default=80, type=int, help="Resolution of plot")
        parser.add_argument("--extension", type=float, default=0.1, help="Extension as fraction of repeat signal around STR region to plot")
        parser.add_argument("--zoom", type=int, default=500, help="Region around prefix and suffix to plot")
        parser.add_argument("--log_level", default='warning', choices=['error', 'warning', 'info', 'debug'], help="Log level")
        args = parser.parse_args(argv)
        logger.init(log_level=args.log_level)
        # index/load reads
        import matplotlib.pyplot as plt
        if not os.path.isfile(args.f5Index):
            logger.log("Main: Fast5 index file does not exist.", logger.log_type.Error)
            exit(1)
        f5Index = fast5Index.fast5Index(args.f5Index)
        # create output directory if needed
        if args.output:
            os.makedirs(args.output, exist_ok=True)
        def tsv_iter(input):
            if input:
                with open(input, 'r') as fp:
                    for line in fp:
                        if not line.startswith('ID'):
                            yield line.strip().split('\t')
            else:
                for line in sys.stdin:
                    if not line.startswith('ID'):
                        yield line.strip().split('\t')
        for record in tsv_iter(args.counts):
            ID, target, strand, count, score_prefix, score_suffix, _, offset, ticks = record[:9]
            offset = int(offset)
            ticks = int(ticks)
            score_prefix = float(score_prefix)
            score_suffix = float(score_suffix)
            raw_signal = f5Index.get_raw(ID)
            if raw_signal is not None:
                flt_signal = sp.medfilt(raw_signal, kernel_size=3)
                flt_signal = (flt_signal - np.median(flt_signal)) / np.std(flt_signal)
                prefix_extend = max(0, offset - int(ticks * args.extension))
                suffix_extend = min(len(flt_signal), offset + ticks + int(ticks * args.extension))
                prefix_begin = max(offset - args.zoom, 0)
                prefix_end = prefix_begin + args.zoom * 2
                suffix_begin = offset + ticks - args.zoom
                suffix_end = min(len(flt_signal), suffix_begin + args.zoom * 2)
                plt.figure(num=None, figsize=(args.width, args.height), dpi=args.dpi, facecolor='w', edgecolor='k')
                plt.subplot(2,1,1)
                plt.plot(flt_signal[prefix_extend:suffix_extend], 'k-', linewidth=0.5, label='genome')
                plt.plot(np.arange(ticks) + (offset - prefix_extend), flt_signal[offset:offset+ticks], 'b-', linewidth=1.0, label='STR')
                plt.legend()
                plt.title("Read {} with {} repeats".format(ID, count))
                plt.subplot(2,2,3)
                plt.plot(flt_signal[prefix_begin:prefix_end], 'k-', label='prefix')
                plt.plot(np.arange(args.zoom, 2 * args.zoom), flt_signal[prefix_begin + args.zoom:prefix_end], 'b-')
                plt.axvline(args.zoom, color='red', label='STR begin')
                plt.legend()
                plt.title("Prefix region with score {:.2f}".format(score_prefix))
                plt.subplot(2,2,4)
                plt.plot(flt_signal[suffix_begin:suffix_end], 'k-', label='suffix')
                plt.plot(flt_signal[suffix_begin:suffix_end - args.zoom], 'b-')
                plt.axvline(args.zoom, color='red', label='STR end')
                plt.legend()
                plt.title("Suffix region with score {:.2f}".format(score_suffix))
                plt.tight_layout()
                if args.output:
                    f_name = os.path.join(args.output, '_'.join([target, count, ID]) + '.' + args.format)
                    plt.savefig(f_name, )
                else:
                    plt.show()
            else:
                logger.log("Plot: No fast5 for ID {id}".format(id=sam_record.QNAME), logger.log_type.Warning)
        logger.close()
        exit(0)


# main
if __name__ == '__main__':
    signal(SIGPIPE,SIG_DFL)
    main()
