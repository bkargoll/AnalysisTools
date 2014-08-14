#!/usr/bin/env python

import pickle
from optparse import OptionParser

ptval = 0

parser = OptionParser()
parser.add_option("--pt",help="pt range to use",type="int",dest="pt")
(options,args) = parser.parse_args()

ptval = options.pt

f = open('MuonEfficiencies_ISO_Run_2012ReReco_53X.pkl','r')
dict = pickle.load(f)

eta = ['ptabseta<0.9','ptabseta0.9-1.2','ptabseta1.2-2.1','ptabseta2.1-2.4']
pt = ['10_20','20_25','25_30','30_35','35_40','40_50','50_60','60_90','90_140','140_300']

for i in range(len(eta)):
	print pt[ptval]+","+eta[i]#+"\n"
	print dict['combRelIsoPF04dBeta<012_Tight'][eta[i]][pt[ptval]]['data/mc']['err_hi']
	print dict['combRelIsoPF04dBeta<012_Tight'][eta[i]][pt[ptval]]['data/mc']['err_low']
	print "\n"
