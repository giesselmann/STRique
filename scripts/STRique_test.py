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
    def test_STR(self):
        model_file = os.path.join( os.path.dirname(os.path.realpath(__file__)), '..', 'models', 'template_median68pA6mer.model')
        pm = STRique.pore_model(model_file)
        backbone = ''.join(random.choice('ACTG') for _ in range(1000))
        prefix = 'CGGCAGCCGAACCCCAAACAGCCACCCGCCAGGATGCCGCCTCCTCACTCACCCACTCGCCACCGCCTGCGCCTCCGCCGCCGCGGGCGCAGGCACCGCAACCGCAGCCCCGCCCCGGGCCCGCCCCCGGGCCCGCCCCGACCACGCCCC'
        suffix = 'TAGCGCGCGACTCCTGAGTTCCAGAGCTTGCTACAGGCTGCGGTTGTTTCCCTCCTTGTTTTCTTCTGGTTAATCTTTATCAGGTCTTTTCTTGTTCACCCTCAGCGAGTACTGTGAGAGCAAGTAGTGGGGAGAGAGGGTGGGAAAAAC'
        repeat = 'GGCCCC'
        dt = STRique.repeatDetection(model_file, repeat, prefix[-50:], suffix[:50], prefix, suffix)
        for i in range(100, 800, 100):
            print('Test repeat length:', i)
            seq = backbone[:500] + prefix + repeat * i + suffix + backbone[-500:]
            sig = pm.generate_signal(seq, samples=8)
            n, score_prefix, score_suffix, p, ticks, prefix_end = dt.detect(sig)
            self.assertEqual(n, i)
            
# main function
if __name__ == '__main__':
    unittest.main()