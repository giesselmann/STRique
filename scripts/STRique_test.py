# \TEST\-------------------------------------------------------------------------
#
#  CONTENTS      : STRique tests
#
#  DESCRIPTION   : Raw nanopore signal repeat detection pipeline tests
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
import os
import random
import unittest
# private imports
import STRique




# Test cases
class DetectionTest(unittest.TestCase):
    # Test if pipeline is able to quantify repeat count for hexamere repeat GGCCCC
    def test_Detection(self):
        model_file = os.path.join( os.path.dirname(os.path.realpath(__file__)), '..', 'models', 'r9_4_450bps.model')
        pm = STRique.pore_model(model_file)
        backbone = ''.join(random.choice('ACTG') for _ in range(2000))
        prefix = 'CGGCAGCCGAACCCCAAACAGCCACCCGCCAGGATGCCGCCTCCTCACTCACCCACTCGCCACCGCCTGCGCCTCCGCCGCCGCGGGCGCAGGCACCGCAACCGCAGCCCCGCCCCGGGCCCGCCCCCGGGCCCGCCCCGACCACGCCCC'
        suffix = 'TAGCGCGCGACTCCTGAGTTCCAGAGCTTGCTACAGGCTGCGGTTGTTTCCCTCCTTGTTTTCTTCTGGTTAATCTTTATCAGGTCTTTTCTTGTTCACCCTCAGCGAGTACTGTGAGAGCAAGTAGTGGGGAGAGAGGGTGGGAAAAAC'
        repeat = 'GGCCCC'
        dt = STRique.repeatCounter(model_file)
        dt.add_target('c9orf72', repeat, prefix, suffix)
        print('Test GGCCCC repeat')
        for i in range(100, 301, 100):
            print('Test repeat length:', i)
            seq = backbone[:1000] + prefix + repeat * i + suffix + backbone[-1000:]
            sig = pm.generate_signal(seq, samples=8)
            n, score_prefix, score_suffix, p, offset, ticks, mod_pattern = dt.detect('c9orf72', sig, '+')
            self.assertEqual(n, i)
            #print('Expected: ', i, ' detected: ', n)

    # Test special case of repeats shorter than kmer in pore model
    # this requires interpolation and limits detection resolution
    def test_Interpolation(self):
        model_file = os.path.join( os.path.dirname(os.path.realpath(__file__)), '..', 'models', 'r9_4_450bps.model')
        pm = STRique.pore_model(model_file)
        backbone = ''.join(random.choice('ACTG') for _ in range(2000))
        prefix = 'AGCGGGCCGGGGGTTCGGCCTCAGTCAGGCGCTCAGCTCCGTTTCGGTTTCACTTCCGGTGGAGGGCCGCCTCTGAGCGGGCGGCGGGCCGACGGCGAGCGCGGGCGGCGGCGGTGACGGAGGCGCCGCTGCCAGGGGGCGTGCGGCAGC'
        suffix = 'GAGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCTGGGCCTCGAGCGCCCGCAGCCCACCTCTCGGGGGCGGGCTCCCGGCGCTAGCAGGGCTGAAGAGAAGATGGAGGAGCTGGTGGTGGAAGTGCGGGGCTCCAATGGCGCTTTCTACAA'
        repeat = 'GCG'
        dt = STRique.repeatCounter(model_file)
        dt.add_target('fmr1', repeat, prefix, suffix)
        print('Test GCG repeat')
        for i in range(100, 301, 100):
            print('Test repeat length:', i)
            seq = backbone[:1000] + prefix + repeat * i + suffix + backbone[-1000:]
            sig = pm.generate_signal(seq, samples=8)
            n, score_prefix, score_suffix, p, offset, ticks, mod_pattern = dt.detect('fmr1', sig, '+')
            self.assertEqual(n, i)
            #print('Expected: ', i, ' detected: ', n)

    # test normalization for short prefix/suffix sequences
    def test_Normalization(self):
        model_file = os.path.join( os.path.dirname(os.path.realpath(__file__)), '..', 'models', 'r9_4_450bps.model')
        pm = STRique.pore_model(model_file)
        prefix = 'CGGCAGCCGAACCCCAAACAGCCACCCGCCAGGATGCCGCCTCCTCACTCACCCACTCGCCACCGCCTGCGCCTCCGCCGCCGCGGGCGCAGGCACCGCAACCGCAGCCCCGCCCCGGGCCCGCCCCCGGGCCCGCCCCGACCACGCCCC'
        suffix = 'TAGCGCGCGACTCCTGAGTTCCAGAGCTTGCTACAGGCTGCGGTTGTTTCCCTCCTTGTTTTCTTCTGGTTAATCTTTATCAGGTCTTTTCTTGTTCACCCTCAGCGAGTACTGTGAGAGCAAGTAGTGGGGAGAGAGGGTGGGAAAAAC'
        repeat = 'GGCCCC'
        dt = STRique.repeatCounter(model_file)
        dt.add_target('c9orf72', repeat, prefix, suffix)
        print('Test normalization on short prefix/suffix sequences')
        for i in range(10, 100, 10):
            print('Test repeat length:', i)
            seq = prefix + repeat * i + suffix
            sig = pm.generate_signal(seq, samples=8)
            n, score_prefix, score_suffix, p, offset, ticks, mod_pattern = dt.detect('c9orf72', sig, '+')
            self.assertEqual(n, i)
            #print('Expected: ', i, ' detected: ', n)

    # test base modification function
    def test_Modification(self):
        model_file = os.path.join( os.path.dirname(os.path.realpath(__file__)), '..', 'models', 'r9_4_450bps.model')
        mod_model_file = os.path.join( os.path.dirname(os.path.realpath(__file__)), '..', 'models', 'r9_4_450bps_mCpG.model')
        pm = STRique.pore_model(model_file)
        pm_mod = STRique.pore_model(mod_model_file)
        backbone = ''.join(random.choice('ACTG') for _ in range(2000))
        prefix = 'CGGCAGCCGAACCCCAAACAGCCACCCGCCAGGATGCCGCCTCCTCACTCACCCACTCGCCACCGCCTGCGCCTCCGCCGCCGCGGGCGCAGGCACCGCAACCGCAGCCCCGCCCCGGGCCCGCCCCCGGGCCCGCCCCGACCACGCCCC'
        suffix = 'TAGCGCGCGACTCCTGAGTTCCAGAGCTTGCTACAGGCTGCGGTTGTTTCCCTCCTTGTTTTCTTCTGGTTAATCTTTATCAGGTCTTTTCTTGTTCACCCTCAGCGAGTACTGTGAGAGCAAGTAGTGGGGAGAGAGGGTGGGAAAAAC'
        repeat = 'GGCCCC'
        dt = STRique.repeatCounter(model_file, mod_model_file=mod_model_file)
        dt.add_target('c9orf72', repeat, prefix, suffix)
        print('Test GGCCCC repeat')
        for i in range(100, 301, 100):
            print('Test repeat length:', i)
            seq = backbone[:1000] + prefix + repeat * i + suffix + backbone[-1000:]
            sig_base = pm.generate_signal(seq, samples=8, noise=True)
            n, score_prefix, score_suffix, p, offset, ticks, mod_pattern_base = dt.detect('c9orf72', sig_base, '+')
            sig_mod = pm_mod.generate_signal(seq, samples=8, noise=True)
            n, score_prefix, score_suffix, p, offset, ticks, mod_pattern_mod = dt.detect('c9orf72', sig_mod, '+')
            #print("\n".join([mod_pattern_base, mod_pattern_mod]))
            self.assertEqual(n, i)




# main function
if __name__ == '__main__':
    unittest.main()
