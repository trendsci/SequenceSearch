#! /home/sershn/anaconda2/bin/python

import unittest
import seqa
import sys
import os

class test_seq_parse_uniprot(unittest.TestCase):
  """Test fetching and parsing uniprot fasta sequence from the web"""
  def test_func_seq_prase_uniprot(self):
    """Test to make sure uniprot sequences are being prased to fasta correctly
       fastaid is the fasta header, and fasta is the fasta formatted sequence"""
    fastaid, fasta = seqa.seq_parse_uniprot("P15337")
    self.assertEqual(fastaid.rstrip(), 
    ">sp|P15337|CREB1_RAT Cyclic AMP-responsive element-binding protein 1 OS=Rattus norvegicus GN=Creb1 PE=1 SV=1")
    self.assertEqual(fasta,
    "MTMDSGADNQQSGDAAVTEAESQQMTVQAQPQIATLAQVSMPAAHATSSAPTVTLVQLPNGQTVQVHGVIQAAQPSVIQSPQVQTVQSSCKDLKRLFSGTQISTIAESEDSQESVDSVTDSQKRREILSRRPSYRKILNDLSSDAPGVPRIEEEKSEEETSAPAITTVTVPTPIYQTSSGQYIAITQGGAIQLANNGTDGVQGLQTLTMTNAAATQPGTTILQYAQTTDGQQILVPSNQVVVQAASGDVQTYQIRTAPTSTIAPGVVMASSPALPTQPAEEAARKREVRLMKNREAARECRRKKKEYVKCLENRVAVLENQNKTLIEELKALKDLYCHKSD".strip()   )    


class test_correct_finding_of_roi_in_sequence(unittest.TestCase):
  """Test to check that amino acids are found in the correct location"""
  class CGI_dictionary_spoof_function():
    def __init__(self, dictToUse="normal"):
          self.dictToUse = dictToUse
          testSeq = "asdfa"
          if dictToUse == "normal":
            self.dictToReturn = {"roi":self.f("A"), "roitype":self.f("normal"), 
                     "advanced_search_term":self.f("none"), 
                     "seqsource":self.f("user"), "seq":self.f(testSeq), 
                     "dpi":self.f("96"), "filetype":self.f("png"), 
                     "colorscheme":self.f("standard")}
          if dictToUse == "regex":
            self.dictToReturn = {"roi":self.f("[A]"), "roitype":self.f("regex"), 
                     "advanced_search_term":self.f("none"), 
                     "seqsource":self.f("user"), "seq":self.f(testSeq), 
                     "dpi":self.f("96"), "filetype":self.f("png"), 
                     "colorscheme":self.f("standard")}

    def __call__(self):
      return self.dictToReturn
#    def __getitem__(self,item):
#      return self.dictToReturn[item]

    class f():
      """helps spoof the CGI_parameters["parameter"].value call"""
      def __init__(self, param):
        self.value = param
    

  def find_roi_in_seq_generic(self, kind):
    seqa.parse_CGI_param = self.CGI_dictionary_spoof_function(
    dictToUse=kind)
    consolOut = sys.stdout
    seqa.sys.stdout = open(os.devnull,'w') #open("junk.txt","w")   
    found_in_sequence = seqa.main()
    sys.stdout = consolOut
    return found_in_sequence

  def test_find_roi_in_seq_normal(self):
    found_in_sequence = self.find_roi_in_seq_generic("normal")
    self.assertEqual(found_in_sequence, {'A': [1, 5]})

  def test_find_roi_in_seq_regex(self):
    found_in_sequence = self.find_roi_in_seq_generic("regex")
    self.assertEqual(found_in_sequence, [(1, 2), (5, 6)])

