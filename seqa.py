#! /home/sershn/anaconda2/bin/python
#import importfasta
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cgi
import sys
import urllib2
import re
import datetime
from sequtil.sequtil import * #my own package with dependancies for this program
import signal

def printException(text, target="html"):
  if target == "html":  
    print "Content-type: text/html\n\n"
    print "<br>"
    print text


#def sig_hand(signal, frame):
#  print "all good"
#  sys.exit(0)

def main(printer):
  """Description here"""
  #time this program, set start time.
  time_start = datetime.datetime.now()
  commenter("Program started at {}".format(datetime.datetime.utcnow()))
  #save error stream to file
  sys.stderr = open('stderr_seqa_py.txt', 'w')

  ##Overload print statement, instead of directly printing to screen,
  ##we will save output to a varible and print the final page at once.
  sys.stdout = printer #html_printer(stdout=sys.stdout) #assign new reference

  #enable python web presentation on dreamhost
  #html_setup = "Content-type: text/html\n\n"
  #print "Content-type: text/html\n\n"
  #read in the POST or GET parameters passed from the HTML GUI
  try:
    arguments_web = parse_CGI_param()
    commenter("arguments supplied via CGI: {}".format(arguments_web))

  # define variables
  ########################VARIABLE DEFINITIONS###################
  ###############################################################  
  # expected variables:
  # roi - residue of interest (string)
  # roitype - normal string or regex expression ('normal'/'regex')
  # etc...
    seq_source = str(arguments_web["seqsource"].value)
    dpi_web = int(arguments_web["dpi"].value)
    filetype = arguments_web["filetype"].value
    colorscheme = arguments_web["colorscheme"].value
    roitype = arguments_web["roitype"].value
    seq_input = arguments_web["seq"].value.upper() 
  #roi is "Residues Of Interest"
  #get the search term from user input (e.g. residues or RegEx expression)
    roi_raw = arguments_web["roi"].value.upper()
  #disaply amino acid type for each residue?
    try:
      show_res_label = arguments_web["label"].value
    except:
      show_res_label = 0
  # done defining variables from CGI passed parameters
  ##############################################################
  ##############################################################
  except Exception as e:
    return "<h3>Error</h3>Error parsing provided parameters from form.<br>Did you leave a field empty or filled the form incorrectly?<br>Try to go back and fix your submission<br><br>Program terminated.<br><br>Traceback:<br>({})".format(e)
  #output the arguments the script received from the web (via POST or GET)

  ## Simple HTML openning with head section ##
  ############################################
  ############################################
  javascript_1 = """
  <script type="text/javascript">
    var mytextbox = document.getElementById('mytext');
    var mydropdown = document.getElementById('dropdown');

    mydropdown.onchange = function(){
          mytextbox.value = mytextbox.value  + this.value; //to appened
         //mytextbox.innerHTML = this.value;
    }

   function showhide(id) {
    var e = document.getElementById(id);
    e.style.display = (e.style.display == 'block') ? 'none' : 'block';
 }    

  </script>
  """  
  print "<html><head><link rel=\"stylesheet\" type=\"text/css\" href=\"/css/svg1.css\">{jscript_1}</head>".format(jscript_1=javascript_1)
  print "<body>"

 
  #check where to get ASCII sequence from (e.g. uniprot or user input)
  if seq_source == "uniprot":
    uniprot_id = ''.join(seq_input.split())
    #get fasta_id and sequence in fasta format from uniprot
    fasta_id, sequence_web = seq_parse_uniprot(uniprot_id) 
    sequence_source = "Fasta from uniprot {}<br>{}".format(uniprot_id, fasta_id[:])
    printp("Fasta from uniprot {}<br>{}".format(uniprot_id,seq_printer(fasta_id.rstrip(),block_size=60,line_size=60,numbered=False,seperator="..<br>")),tag="style='margin:0;'")
   # printp("Sequence:<br>{}".format(seq_printer(sequence_web)), new_line=False)

  elif seq_source == "user":
    sequence_web = ''.join(seq_input.split())
    sequence_source = "Sequence from user input"

  seq = sequence_web
  seq_len = len(seq)
  printp("Sequence (length: {seql}):<br>{seq}".format(seql=seq_len,seq=seq_printer(sequence_web)), new_line=False)


  #Check what type of search the user is performing (e.g. regular or regex)
  #format the search term according to search type
  if roitype == "normal":
    roi = set([letter for letter in roi_raw])
  elif roitype == "regex":
    roi = roi_raw
    roi_user_print = roi_raw
    try: #check if search term present in text file
         #if it is, fetch the correct regex to use instead of term 
      with open("advanced_search","r") as f:
        for line in f:
          if roi == line.split()[0].upper(): 
            roi = line.split()[1].upper()
            break
    except:
      commenter("Error raised in open(advanced_search)")


##############################
#######generate bar plot######
  time_before_figure = datetime.datetime.now() - time_start
  fig = plt.figure(figsize=(8,2))

  plt2 = fig.add_subplot(111)
  colors = ['blue','green','red','black','pink']
  #dictionary of color-pellate types
  colors = {"standard": plt.cm.Dark2(np.linspace(0,0.85,len(roi))),
            "bytype": {"D":"#ff6666","E":"red",
                       "R":"Blue", "H":"cyan", "K":"#6666ff"},
            "black":'black',
            "dark": plt.cm.Dark2(np.linspace(0,0.85,len(roi)))
           }
  yheights = []
  xheights = []
  if roitype == "regex":
    commenter("Advanced search (regex) query: {}".format(roi))
    resultRegEx = re.finditer(roi, (seq))
    print "Found matches at position(s):<br>"
    regexList = []
    za = resultRegEx
    for m in resultRegEx:
      regexList.append((int(m.start()+1),int(m.end()+1)))
      commenter("regex result: {}".format(m.group()))
      print "Residuese> {} - {}.<br>".format(m.start()+1,m.end())
    print "<br>"
  nt = 1 #set to 1 so sequence starts at 1
  check = 1
  
  #updated faster algorithm for "normal" search. >2x faster than
  # normal_old_slow.
  if roitype == "normal":
    #make dictionary with ROIs as Keys:
    commenter("ROI (search terms): {}".format(roi))
    roi_dict = dict((term,[]) for term in roi)
    commenter("roi_dict that was generated: {}".format(roi_dict))
    #y_locations_of_roi_in_seq = [(roi_dict[char].append(pos+1),char) for pos, char in enumerate(seq) if char in roi]
    for pos, char in enumerate(seq):
      if char in roi: 
        roi_dict[char].append(pos+1)
    
    to_print = '' 
    for key, val in roi_dict.iteritems():
      to_print += "> {} @ {}. ({}) <br>".format(key,str(val)[1:-1],len(val))
    printp("Location of search terms:<br>{}".format(to_print), title="title='Search term @ locations found. (number of occurances)'",new_line=False)

    commenter("roi_dict: {}".format(roi_dict))
    color_num = 0 #used as counter to enumerate different color to each set of bars
    for (a_roi, locations) in (roi_dict.iteritems()):
      xheights = locations
      y_height = 1 #can set this variable somewhere outside of this function
      #so it can be used by others
      yheights = [y_height for a in locations]
      if locations:
        rects = plt2.bar(xheights,yheights,color=colors["standard"][color_num],
                       alpha=1,width=0.9,linewidth=0,
                       label=a_roi, gid="ssRectTest")
      color_num += 1

  elif roitype == "regex":
    commenter("generating regex figure")
    
    regexListX = [x for x,xf in regexList]
    regexListX = []
    for element in regexList:
      for i in range(element[0],element[1]):
        regexListX.append(i)
    regexListY = [1 for x in regexListX]
    rects = plt2.bar(regexListX,regexListY,color=colors["standard"][1],alpha=1,width=0.9,
            linewidth=0,align='edge')
    print "<done>"
  plt.legend(loc='center left', bbox_to_anchor=(1.01,0.5), fontsize=10)

  #plot gray divide lines
  #for item in np.arange(1,len(seq)+1,1): 
  #  plt.axvline(item,color='gray',
  #              linewidth=0.2,linestyle='-',alpha=0.5)  
  

  plt.subplots_adjust(top=0.8,left=0.05,right=0.9,bottom=0.25)

  plt.ylim(0,1)
  plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off') # labels along the bottom edge are off
  if seq_len > 60:
    plt.xticks(np.arange(0,seq_len,round(seq_len/15,-1)))
  ax1 = plt.subplot(111)
  ax = fig.gca()
  ax.set_axis_bgcolor('#f7f7f7')
 # ax = plt.gca()
  #for x,y in ((2.5,rects[0].get_y()),(4.5,-0.2)):
  #a1 = rects[0].get_y()
  #used for residue labels
  a1 = 0
  a2 = range(1,seq_len+1,1)
  a3 = []
  for a in a2:
    a3.append((a+0.5,a1))

  if show_res_label == "AA_each":
    for i, (x,y) in enumerate(a3):
      ax1.text(x,y,seq[i],fontsize=0.75*96*3/float(seq_len),
           ha='center', va='bottom')#,transform=ax.transAxes)
##  ax1.text(3,0,"W",fontsize=0.75/1.8*dpi_web/float(len(seq)),
##           ha='center', va='bottom')#,transform=ax.transAxes)
  title_roi_format = ''.join(roi)
  if roitype == "normal":
    plt.title("Location of amino acids: %s." 
             % ', '.join(roi))
  elif roitype == "regex":
    plt.title("Location of: %s" % roi_user_print)
  plt.xlabel("Residue number")
#  plt.ylabel("AA present")
#  plt2.bar(xheights,yheights,color='blue',alpha=1,#width=0.05,
#          edgecolor='blue',align='center')
  ax1.set_xlim(1,seq_len+1)
#  print "<br> almost <br>"
  plt.savefig("../test.%s"%filetype, dpi=dpi_web)
  time_after_figure = datetime.datetime.now() - time_start
  time_to_draw_figure = time_after_figure - time_before_figure
 # print "<object data=\"test.svg\" type=\"image/svg+xml\"></object>"
  print "<img src=\"test.%s\" alt=\"test png\" width=\"650\">"%filetype
  print "<br>"
  print "<a href=\"test.%s\" target=\"_blank\">View full size image</a>"%filetype
  print "<br>"
 
  print """<a href="javascript:showhide('debug_info')">
        Show/hide debug information
    </a>
   <div id="debug_info" style="display:none;">
    """
  
  printp(sequence_source,tag="style='margin:0;'")
  printp("Timing of program:<br>Time before drawing figure: {}<br>Time after drawing figure: {}<br>Time to draw figure: {}".format(time_before_figure,time_after_figure,time_to_draw_figure),tag="style='padding:0;margin:0;'")
  print "</div>"

  print "</body></html>"

  return  sys.stdout.get_HTML()
  try: 
    if roitype == "normal": return roi_dict
    if roitype == "regex" : return regexList
  except:
    pass

if __name__ == "__main__":
  html_header = "Content-type: text/html\n\n"
  p = html_printer()
  try:
    content = main(printer=p)
  except Exception as e:
    content = "Exception occured, critical. Program terminated.<br><br>Traceback:<br>{}".format(e)
  sys.stdout = p.get_original_stdout()
  print  "{}{}".format(html_header,content)
