# \TEST\---------------------------------------------------------------
#
#  CONTENTS      : Raw nanopore signal repeat detection pipeline test
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
import os
import random
import unittest

# private imports
import STRique

# Test cases

class DetectionTest(unittest.TestCase):
    # Test if pipeline is able to quantify repeat count for hexamere repeat GGCCCC
    def test_Detection(self):
        model_file = os.path.join( os.path.dirname(os.path.realpath(__file__)), '..', 'models', 'template_median68pA6mer.model')
        pm = STRique.pore_model(model_file)
        backbone = ''.join(random.choice('ACTG') for _ in range(2000))
        prefix = 'CGGCAGCCGAACCCCAAACAGCCACCCGCCAGGATGCCGCCTCCTCACTCACCCACTCGCCACCGCCTGCGCCTCCGCCGCCGCGGGCGCAGGCACCGCAACCGCAGCCCCGCCCCGGGCCCGCCCCCGGGCCCGCCCCGACCACGCCCC'
        suffix = 'TAGCGCGCGACTCCTGAGTTCCAGAGCTTGCTACAGGCTGCGGTTGTTTCCCTCCTTGTTTTCTTCTGGTTAATCTTTATCAGGTCTTTTCTTGTTCACCCTCAGCGAGTACTGTGAGAGCAAGTAGTGGGGAGAGAGGGTGGGAAAAAC'
        repeat = 'GGCCCC'
        dt = STRique.repeatDetection(model_file, repeat, prefix[-50:], suffix[:50], prefix, suffix)
        print('Test GGCCCC repeat')
        for i in range(100, 800, 100):
            print('Test repeat length:', i)
            seq = backbone[:1000] + prefix + repeat * i + suffix + backbone[-1000:]
            sig = pm.generate_signal(seq, samples=8)
            n, score_prefix, score_suffix, p, ticks, prefix_end = dt.detect(sig)
            self.assertEqual(n, i)
            #print('Expected: ', i, ' detected: ', n)
         
    # Test special case of repeats shorter than kmer in pore model
    # this requires interpolation and limits detection resolution          
    def test_Interpolation(self):
        model_file = os.path.join( os.path.dirname(os.path.realpath(__file__)), '..', 'models', 'template_median68pA6mer.model')
        pm = STRique.pore_model(model_file)
        backbone = ''.join(random.choice('ACTG') for _ in range(2000))
        prefix = 'AGCGGGCCGGGGGTTCGGCCTCAGTCAGGCGCTCAGCTCCGTTTCGGTTTCACTTCCGGTGGAGGGCCGCCTCTGAGCGGGCGGCGGGCCGACGGCGAGCGCGGGCGGCGGCGGTGACGGAGGCGCCGCTGCCAGGGGGCGTGCGGCAGC'
        suffix = 'GAGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCTGGGCCTCGAGCGCCCGCAGCCCACCTCTCGGGGGCGGGCTCCCGGCGCTAGCAGGGCTGAAGAGAAGATGGAGGAGCTGGTGGTGGAAGTGCGGGGCTCCAATGGCGCTTTCTACAA'
        repeat = 'GCG'
        dt = STRique.repeatDetection(model_file, repeat, prefix[-50:], suffix[:50], prefix, suffix)
        print('Test GCG repeat')
        for i in range(100, 800, 100):
            print('Test repeat length:', i)
            seq = backbone[:1000] + prefix + repeat * i + suffix + backbone[-1000:]
            sig = pm.generate_signal(seq, samples=8)
            n, score_prefix, score_suffix, p, ticks, prefix_end = dt.detect(sig)
            self.assertEqual(n, i)
            #print('Expected: ', i, ' detected: ', n)

            
# main function
if __name__ == '__main__':
    unittest.main()