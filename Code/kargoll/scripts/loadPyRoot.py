from os import getenv
from sys import exit
# check that ROOT can be imported
if not getenv('ROOTSYS')+'/lib' in getenv('LD_LIBRARY_PATH'):
    print "Environment not set for usage of PyROOT. Please configure your environment:"
    print "\texport LD_LIBRARY_PATH=$ROOTSYS/lib:$PYTHONDIR/lib:$LD_LIBRARY_PATH"
    print "\texport PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH"
    exit("ERROR: Cannot execute PyROOT.")

import ROOT