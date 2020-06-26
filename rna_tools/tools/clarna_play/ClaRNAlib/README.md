This is ClaRNA 3.0 , prepared by @mmagnus for rna-tools that works with Python3.

ClaRNA
===========================================================

Cite & authors:

[1]	T. Waleń, G. Chojnowski, P. Gierski, and J. M. Bujnicki, “ClaRNA: a classifier of contacts in RNA 3D structures based on a comparative analysis of various classification schemes.,” Nucleic Acids Research, vol. 42, no. 19, pp. e151–e151, Oct. 2014.

`--save-graph=some-file.json`
-----------------------------------------------------------

Save-scores is option used only for debugging/development purposes.

If you want to save output of clarna in JSON format use --save-graph=some-file.json option.

To my knowledge the clarna library (files in lib directory) are the same in repo and on the web server.

ClaRNA requires:

sudo pip install simplejson==2.6.1 networkx==1.8.1 scipy

[!] somehow I could not install scipy via pip. I run sudo apt-get install python-scipy.

See a repo of ClaRNA for more!
