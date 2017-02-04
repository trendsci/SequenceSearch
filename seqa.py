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
#from scipy.interpolate import spline


def calculateRoi(seq,roi,roi_label='0'):
  #seq_w_roi = []
  roi_in_seq = 0
  for char in seq:
    if char == roi:
      roi_in_seq += 1
  return roi_in_seq

#  qres = [[i,n] for i,n in enumerate(seq)]

def generate_qplot(seq,roi,sequence_file='protein',one_plot=0,domains={},
                   colors='standard',savefig=True):
#  pbar = pb.progress_indicator()
#  pbar.max_progress = len(seq)
  val = 4 #+- residues that are effected by nearby ROI. Default=4.
  height = val #hight of plotted ROI
  domain_line_width = wid = 1 #width for domain borders
  domain_line_color = col = 'black' #color for domain borders
  if 'x' in roi: 
    plt.subplot(2, 1, 1)
  for domain in domains.keys():
    if domains[domain][2] == 1:
      plt.axvline(x=domains[domain][0], linewidth=domain_line_width,
                color=domain_line_color)
      #print domains[domain][0]
  plt.ylim([0,height*4])
  plt.xlim(0,len(seq))
  qres = [[i,n,0] for i,n in enumerate(seq)]
  #print qres

  colors_l = ['#E60606','black','yellow','orange','red','green','blue',
            'green','red','orange','green','black','gray','orange','purple']#deprecated.
  colors_list = """[230,6,6] #E60606 R - Arg - Arginine
    [198,66,0]  #C64200 K - Lys - Lysine
    [255,102,0] #FF6600 Q - Gln - Glutamine
    [255,153,0] #FF9900 N - Asn - Asparagine
    [255,204,0] #FFCC00 E - Glu - Glutamic_Acid
    [255,204,153] #FFCC99 D - Asp - Aspartic_Acid
    [255,255,153] #FFFF99 H - His - Histidine
    [255,255,0] #FFFF00 P - Pro - Proline
    [204,255,204] #CCFFCC Y - Tyr - Tyrosine
    [204,153,255] #CC99FF W - Trp - Tryptophan
    [204,255,153] #CCFF99 S - Ser - Serine
    [0,255,153] #00FF99 T - Thr - Threonine
    [0,255,0] #00FF00 G - Gly - Glycine
    [204,255,255] #CCFFFF A - Ala - Alanine
    [153,204,255] #99CCFF M - Met - Methionine
    [0,255,255] #00FFFF C - Cys - Cysteine
    [0,204,255] #00CCFF F - Phe - Phenylalanine
    [51,102,255]  #3366FF L - Leu - Leucine
    [0,0,255] #0000FF V - Val - Valine
    [0,0,128] #000080 I - Ile - Isoleucine"""
  colors_list = zip(colors_list.split()[2::7],colors_list.split()[1::7])
  colors_dict = dict(colors_list)
  if colors == 'bw': colors_dict = {}
  for c,aminoacid in enumerate(roi):
    print 'Working on residue type:',aminoacid
    if aminoacid == 'x':
      plt.legend()
      plt.ylim(0,height*4)
      plt.subplot(2, 1, 2)
      continue
    #list(aminoacid)
    if one_plot == True:
      aminoacid=roi
    else:
      for i in qres: i[2]=0 #reset qres
    for i, char in enumerate(seq):
 #     pbar.progress_percent(current_progress = i)
     # print "i",i
      if char in aminoacid:
        for k in range(i-val,i+val):
          if k >= 0 and k <= len(seq)-1:
            num = height - abs(k-i)
            if num >= 0:
              qres[k][2] += height - abs(k-i)
            else:
              qres[k][2] += 0
            #print 'char',char, qres[k][2], 'qres', 'k', k,'i',i
    #print qres
    xval = [ x for x,aa,y in qres]
    yval = [ y for x,aa,y in qres]
    if len(seq) < 500:
      xvalnew = np.linspace(0,len(seq),len(seq)*10)
      yvalnew = spline(xval,yval,xvalnew)
    else:
      yvalnew = yval
      xvalnew = xval
    plot = plt.fill_between(xvalnew, yvalnew,0, color=colors_dict.get(aminoacid,'black'), linewidth=1, label=aminoacid, alpha=1)
    plot_legened = plt.fill(0, 0, color=colors_dict.get(aminoacid,'black'), linewidth=0, label=aminoacid, alpha=1)
    if one_plot == True: break
  #plt.axvline(x=85,linewidth=wid, color=col)
  #plt.axvline(x=160,linewidth=wid, color=col)
  #plt.axvline(x=280,linewidth=wid, color=col)
  plt.ylim(0,height*4)
  plt.xlim(0,len(seq))
  plt.title('Amino acid distribution in: '+sequence_file)
  z = plt.legend()
  if savefig: 
    plt.savefig('seqfig.pdf')
  plt.savefig('seqfig.pdf')
  plt.show()

def findRoi(seq,roi,roi_label):
  seq_w_roi = []
  roi_in_seq = 0
  for char in seq:
    if char == roi:
      roi_in_seq += 1
      seq_w_roi.append(roi_label)
    else:
      seq_w_roi.append(char)
  return seq_w_roi, roi_in_seq

def findRoiDomains(seq,roi,roi_label,domains,n_term_linker = 3):
 # print "Domains used:",
 # print ', '.join(sorted(domains.keys()))
 # print 'Residue of interest', roi,'\n'
  roi_in_protein = 0

  #qyres = []
  #qplot = (qres, qyres)
  for domain in sorted(domains.keys()):
    seq_w_roi = []
    roi_in_domain = 0
    #print '\n'
  #  print domain, 'Residues:',
  #  print domains[domain]
    for char in seq[domains[domain][0]+(n_term_linker-1):domains[domain][1]+n_term_linker]:
      if char in roi or char in roi_label:
        roi_in_domain += 1
        roi_in_protein += 1
        seq_w_roi.append(roi_label[char])
        #seq_w_roi.append('\033[0;34;1;40m'+char+'\033[0m')
      else:
        seq_w_roi.append(char)
  #  print ''.join(seq_w_roi)
    perc = float(roi_in_domain)/len(seq_w_roi)*100
    output = (str(roi_in_domain)+'/'+str(len(seq_w_roi))+'  '+ '%.1f'+"%%\n") % (perc)
  #  print output
  #roiperc = float(roi_in_protein)/len(seq_w_roi)
  #print 'Total', roi,':',roi_in_protein,',','%.1f%%' % (roiperc)

def main():
  #save error stream to file
  sys.stderr = open('error0909.txt', 'w')
  #enable python web presentation on dreamhost
  print "Content-type: text/html\n\n"
  
  # define variables
  sequence_file = ""
  
  arguments_web = cgi.FieldStorage()
  seq_source = str(arguments_web["seqsource"].value)
  if seq_source == "uniprot":
    uniprot_id = arguments_web["seq"].value.strip().replace('\n','').replace('\r','')
    sequence_web_temp = urllib2.urlopen('http://www.uniprot.org/uniprot/'
                                        +uniprot_id+'.fasta')
    print "<pre>"
    print "Fasta from uniprot", uniprot_id, "<br>",sequence_web_temp.readline(),
    print "</pre>",
    sequence_web = sequence_web_temp.read()
    sequence_web = sequence_web.upper().strip().replace('\n','').replace('\r','')
    print "<pre>" 
    print "Sequence:",
    for i in range(1,len(sequence_web),10):
      if i == 1 or (i-1)%40 == 0: print "<br>",'%4s' %(str(i)),
      print sequence_web[i:i+10],

    print "</pre>"
  elif seq_source == "user":
    sequence_web = arguments_web["seq"].value.upper().strip().replace('\n','').replace('\r','')
  dpi_web = int(arguments_web["dpi"].value)
  
  #roi is "Residues Of Interest"
  roi_raw = arguments_web["roi"].value.upper()
  if arguments_web["roitype"].value == "normal":
    roi = set([letter for letter in roi_raw])
  elif arguments_web["roitype"].value == "regex":
    roi = roi_raw

  #disaply amino acid type for each residue?
  try:
    show_res_label = arguments_web["label"].value
  except:
    show_res_label = 0
  
  filetype = arguments_web["filetype"].value
  colorscheme = arguments_web["colorscheme"].value

  """ Arguments required for the program (either GET or POST):
  seq = "ATGA" (any alphanumeric letter or number, length 1-inf)
  roi = "AT" (any alphanumeric, length 1-inf)
  """

  #domains can be redifined from UNIPROT or other databases in the future
  domains = {'1_Q1'  : (-2,100,0),
             '2_KID' : (85,160,1),
             '3_Q2'  : (160,280,1),
             '4_bZip': (270,341,1),
             '5_MTSL1': (302,318,0),
             '6_DNAbinding' :(286,305,0)}
  n_term_linker = 0 #not in use, for future
  c_term_linker = 0 #not in use, for future
  
  #Unicode symbols to enter in text if needed
  possible_labels = [unicode(u'\u2588').encode('utf-8'),
                     unicode(u'\u2592').encode('utf-8'),
                     unicode(u'\u2591').encode('utf-8'),
                     unicode(u'\u2663').encode('utf-8'),
                     unicode(u'\u25B2').encode('utf-8'),
                     unicode(u'\u2206').encode('utf-8')]
  
  #generate dictionary to match roi with UTF8 symbol
  roi_label = {}
  for i, letter in enumerate(roi):
    try:
      roi_label[letter] = possible_labels[i]
    except IndexError as e:
      roi_label[letter] = '*'
  # get fasta sequence
  seq = [sequence_web]#sets seq to a list of one item which is the sequence

  ##Future possible uniprot implementation
  #if seq[0] == 'uniprot' or seq[0] == 'u':
  #  uniprot_seq = ['http://www.uniprot.org/uniprot/'+seq[1]+'.fasta']
  #  print uniprot_seq
  #  seq = importfasta.importWebFasta(uniprot_seq)[0][1]
  #else:
  #  seq = importfasta.importTextFasta(sequence_file)[0][1]
  #print 'Sequence. \nLength:', len(seq),'bp.','\n',seq

  #seq_w_roi = findRoiDomains(seq,roi,roi_label,domains)

  # roin = calculateRoi(seq,roi[0])
#!@  roin = [ calculateRoi(seq,roi1) for roi1 in roi ]
#  print '<br> ok4 <br>'
#!@  roiprec = [ '%.1f%%' % (float(roin1)/len(seq)*100) for roin1 in roin ]
#  print '#======== Residue Composition ========#'
#!@  roiprec = 'Total residue composition:\n'+str(roi)+' : '+ str(roin)+' :  ' + str(roiprec)
#  print roiprec
#  print '<br><br><br> OKKKKKKK <br>'
  #generate_qplot(seq,roi,sequence_file=sequence_file,
  #                domains=domains,colors=colors,savefig=True)
  seq_temp = [a for a in seq[0]]
  seq = seq_temp

##############################
#######generate bar plot######
 
  #fig = plt.figure(figsize=(8,2))
  fig = plt.figure(figsize=(8,2))

  plt2 = fig.add_subplot(111)
  colors = ['blue','green','red','black','pink']
  colors = {"standard": plt.cm.Dark2(np.linspace(0,0.85,len(roi))),
            "bytype": {"D":"#ff6666","E":"red",
                       "R":"Blue", "H":"cyan", "K":"#6666ff"},
            "black":'black'
           }
  #colors = plt.cm.Dark2(np.linspace(0,0.85,len(roi)))
  yheights = []
  xheights = []
  if arguments_web["roitype"].value == "regex":
    print "<mark>RegEx is an experimental feature. Use with caution.</mark><br>"
    resultRegEx = re.finditer(roi, ''.join(seq))
    print "<pre>Sequence:<br>",''.join(seq),"</pre>"
    print "Found RegEx matches at residues:<br>"
    regexList = []
    for m in resultRegEx:
      regexList.append((int(m.start()+1),int(m.end()+1)))
      print m.start()+1
    print "<br>"
  nt = 1 #set to 1 so sequence starts at 1
  check = 1
  if arguments_web["roitype"].value == "normal":
    for n, residue in enumerate(roi):
      for i, aacid in enumerate(seq):
        if aacid == residue:
          yheights.append(1)
          xheights.append(i+nt)
          #print i, nt, aacid,xheights
          if check:
            check = 0
          else:
            check = 1
        else:
          if check == 1: 
            yheights.append(0)
            check = 0
          else:
            yheights.append(0)
            check = 1
          xheights.append(i+nt)
      if colorscheme == "standard":
        rects = plt2.bar(xheights,yheights,color=colors["standard"][n],alpha=1,width=0.9,
            linewidth=0,align='edge',label=residue,gid="ssRectTest")
      
      elif colorscheme == "bytype":
        rects = plt2.bar(xheights,yheights,color=colors["bytype"].setdefault(residue,"black"),alpha=1,width=0.9,
            linewidth=0,align='edge',label=residue)
    
      xheights = []
      yheights = []
  elif arguments_web["roitype"].value == "regex":
    print "generating regex figure"
    
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
  for item in np.arange(1,len(seq)+1,1): 
    plt.axvline(item,color='gray',
                linewidth=0.2,linestyle='-',alpha=0.5)  


  plt.subplots_adjust(top=0.8,left=0.1,right=0.9,bottom=0.25)

  plt.ylim(0,1)
  plt.tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left='off',      # ticks along the bottom edge are off
    right='off',         # ticks along the top edge are off
    labelleft='off') # labels along the bottom edge are off
  if len(seq) > 60:
    plt.xticks(np.arange(0,351,20))
  ax1 = plt.subplot(111)
  ax = fig.gca()
 # ax = plt.gca()
  #for x,y in ((2.5,rects[0].get_y()),(4.5,-0.2)):
  #a1 = rects[0].get_y()
  a1 = 0
  a2 = range(1,len(seq)+1,1)
  a3 = []
  for a in a2:
    a3.append((a+0.5,a1))

  if show_res_label == "AA_each":
    for i, (x,y) in enumerate(a3):
      ax1.text(x,y,seq[i],fontsize=0.75*96*3/float(len(seq)),
           ha='center', va='bottom')#,transform=ax.transAxes)
##  ax1.text(3,0,"W",fontsize=0.75/1.8*dpi_web/float(len(seq)),
##           ha='center', va='bottom')#,transform=ax.transAxes)
  title_roi_format = ''.join(roi)
  if arguments_web["roitype"].value == "normal":
    plt.title("Location of amino acids: %s." 
             % ', '.join(roi))
  elif arguments_web["roitype"].value == "regex":
    plt.title("Location of RegEx: %s" % roi)
  plt.xlabel("Residue number")
#  plt.ylabel("AA present")
#  plt2.bar(xheights,yheights,color='blue',alpha=1,#width=0.05,
#          edgecolor='blue',align='center')
  ax1.set_xlim(1,len(seq)+1)
#  print "<br> almost <br>"
  plt.savefig("test.%s"%filetype, dpi=dpi_web)
#  print "<br><br><br> DONEEE ! <br><br>"
  print "<html><head><link rel=\"stylesheet\" type=\"text/css\" href=\"/css/svg1.css\"></head>"
  print "Figure shown below.<br>"
  print "Number of residues: ",len(seq),"<br>"
  print "<a href=\"test.%s\">Click here to view full size image</a>"%filetype
  print "<br>"
 # print "<object data=\"test.svg\" type=\"image/svg+xml\"></object>"
  print "<body><img src=\"test.%s\" alt=\"test png\" width=\"800\"></body></html>"%filetype
#  print "Post image"


if __name__ == "__main__":
    main()
