#! /home/sershn/anaconda2/bin/python


import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cgi
import sys
import urllib2
import re
import datetime

def printp(text, pre=True, new_line=True, title="",tag=""):
  """Print text surroundd by <pre> ... </pre>"""
  if pre:
    text = "<pre {title} {tag}>{text}</pre>".format(title=title,tag=tag,text=text)
  if new_line:
    print text
  else:
    print text,

def commenter(text, comments_file="seqa_comments.txt",
              append=True, new_line=True, marker=True):
  """Allowes one to print out comments or debug information to a file.
     text = text to be saved in file
     comments_file = location of file 
     append = append existing text or start off with new file
  """
  if new_line:
    text += '\n'
  if marker:
    text = '-->  ' + text
  if append:
    with open(comments_file, "a+") as f:
      f.write(text)
  else:
    with open(comments_file, "w+") as f:
      f.write(text)
  return 1 #simply return 1 to indicate all is well

def parse_CGI_param():
  """parse the cgi parameters passed via GET or POST
     this function can be overloaded to emulae CGI parameters from a script
  """
  arguments_web = cgi.FieldStorage()
  return arguments_web

def seq_parse_uniprot(uniprotID):
  """Takes uniprotID and returns sequence in 
     fasta format of first sequence in uniprot entry
     returns fasta_id, sequence_fasta
  """
  uniprot_id = uniprotID.strip().replace('\n','').replace('\r','')
  sequence_web_object = urllib2.urlopen(
      'http://www.uniprot.org/uniprot/{}.fasta'.format(uniprot_id))
  fasta_id = sequence_web_object.readline()
  sequence_web = sequence_web_object.read()
  sequence_web = sequence_web.upper().strip().replace('\n','').replace('\r','')
  #sequence_web now contains the fasta 
  commenter("Uniprot sequence: {}".format(sequence_web))
  return fasta_id, sequence_web

def seq_printer(text, block_size=10, line_size=40, numbered='left',
                seperator="<br>"):
  """Print sequence in nice block format (similar to protparam output)"""
  f_out = "" # f_out stands for: formatted_output
  for i in range(0,len(text),block_size):
    ij = i + 1
    if ij == 1:
      if numbered: f_out += " %4s "%(str(ij))
    elif (ij-1)%line_size == 0:
      if numbered:
        f_out += "%s %4s "%(seperator,str(ij))
      else:  
        f_out += "%s"%(seperator)
    f_out += text[i:i+block_size] + " "

  return f_out

class html_printer():
  def __init__(self, stdout=sys.stdout, HTML_to_add = "<br>", HTML_total=""):
    self.HTML_to_add = HTML_to_add
    self.HTML_total = HTML_total
    self.stdout = stdout
  def add_HTML(self,HTML_to_add):
    self.HTML_total += HTML_to_add
  def write(self,text):
    self.add_HTML(text)
  def print_HTML(self):
    self.stdout.write(self.HTML_total)
  def get_original_stdout(self):
    return self.stdout

## END OF SEQ_UTIL functions ##
###############################



