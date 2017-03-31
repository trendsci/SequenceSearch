#! /home/sershn/anaconda2/bin/python

#exception textsi
ERROR_UNKNOWN = "<h3>Error</h3>Critical exception occured. Program terminated.<br><br>Traceback:<br>{}"
ERROR_CGI_PARSING = "<h3>Error</h3>Error parsing provided parameters from form.<br>Did you leave a field empty or filled the form incorrectly?<br>Try to go back and fix your submission<br><br>Program terminated.<br><br>Traceback:<br>({})"

#HTML / javascript stuff
# all the needed javascript in one variable

HTML_HEADER = "Content-type: text/html\n\n"

JAVASCRIPT_ALL = """
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
