#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TColor.h"

#include <vector>

#include "tdrstyle.C"
#include "CMS_lumi.C"

#include "userConfig.h"
#include "defineSamplesAndPlots.h"
#include "plottingHelpers.h"

void plotting(){
	if(verbose) std::cout << "--> plotting()" << std::endl;
	gROOT->LoadMacro("tdrstyle.C");
	setTDRStyle();
	gROOT->LoadMacro("CMS_lumi.C");
	writeExtraText = true;
	gStyle->SetOptStat(0);

	// define which samples to use
	std::vector<sample> samples = defineSamples();

	// define which plots to draw
	std::vector<plotInfo> plots = definePlots();

	/* THIS STUFF STILL NEEDS TO BE ADAPTED
	// todo: implement syst
	std::vector<double> syst;
	// lumi + xsec + eid + muid + pileup + trigger
	const int nsyst = 6;
	double qcd[nsyst] = {0.026,0.387,0.083,0.083,0.,0.083};
	double zz4l[nsyst] = {0.026,0.15,0.002,0.006,0.008,0.007};
	double zz2l2q[nsyst] = {0.026,0.15,0.0044,0.0042,0.004,0.007};
	double zz2l2nu[nsyst] = {0.026,0.15,0.0044,0.007,0.011,0.009};
	double wz3lnu[nsyst] = {0.026,0.15,0.0018,0.006,0.042,0.035};
	double wz2l2q[nsyst] = {0.026,0.15,0.007,0.006,0.015,0.014};
	double ww2l2nu[nsyst] = {0.026,0.15,0.0016,0.0055,0.013,0.015};
	double ttbar[nsyst] = {0.026,0.047,0.0014,0.0053,0.021,0.010};
	double tw[nsyst] = {0.026,0.09,0.0011,0.0055,0.294,0.325};
	double tbarw[nsyst] = {0.026,0.09,0.0011,0.0055,0.029,0.054};
	double dyll[nsyst] = {0.026,0.033,0.0015,0.0063,0.060,0.075};
	double dytt[nsyst] = {0.026,0.033,0.0016,0.0055,0.067,0.066};
	double dyemu[nsyst+1] = {0.026,0.055,0.0016,0.0055,0.011,0.011};
	syst.push_back(QuadraticSum(nsyst,qcd));
	syst.push_back(QuadraticSum(nsyst,zz4l));
	syst.push_back(QuadraticSum(nsyst,zz2l2q));
	syst.push_back(QuadraticSum(nsyst,zz2l2nu));
	syst.push_back(QuadraticSum(nsyst,wz3lnu));
	syst.push_back(QuadraticSum(nsyst,wz2l2q));
	syst.push_back(QuadraticSum(nsyst,ww2l2nu));
	syst.push_back(QuadraticSum(nsyst,ttbar));
	syst.push_back(QuadraticSum(nsyst,tw));
	syst.push_back(QuadraticSum(nsyst,tbarw));
	syst.push_back(QuadraticSum(nsyst,dyll));
	syst.push_back(QuadraticSum(nsyst,dytt));
	syst.push_back(QuadraticSum(nsyst+1,dyemu));
	*/

	// test validity of structs
	testInputs(conf, samples, plots);
	
	// create plots
	if(testPlotting){
		plotInfo testPlot = plots.at(0);
		TH1D* datahist = getHisto(conf, testPlot.identifier+"Data", 1, 1);
		manipulateHisto(datahist, testPlot);
		drawPlot(conf, testPlot, datahist, samples);
	}else{
		std::vector<TH1D*> datahists;
		for(unsigned p = 0; p < plots.size(); p++){
			TH1D* datahist = getHisto(conf, plots.at(p).identifier+"Data", 1, 1);
			manipulateHisto(datahist, plots.at(p));
			datahists.push_back(datahist);
		}
		for(unsigned p = 0; p < plots.size(); p++){
			drawPlot(conf, plots.at(p), datahists.at(p), samples);
		}
	}

	// convert all created eps files into pdf
	std::cout << "Plots created. Convert all eps into pdf..." << std::endl;
	system("for i in `ls -1 *.eps`; do epstopdf $i; done");

	/* THIS STUFF STILL NEEDS TO BE ADAPTED
	//TH1D* onejet_signal = getHisto("onejetMC_emu_DY",mcscale.at(mcscale.size()-1),0,infile,syst.at(syst.size()-1));
	//TH1D* onejet_dy = getHisto("onejetMC_tautau_DY",mcscale.at(mcscale.size()-2),2345,infile,syst.at(syst.size()-2));
	//onejet_signal->Scale(1./onejet_signal->Integral());
	//onejet_dy->Scale(1./onejet_dy->Integral());
	//TH1D* onejet_ratio = getDataMC(onejet_signal,onejet_dy);
	//drawPlot(onejet_signal, onejet_dy, onejet_ratio, "Custom MC", "Official MC", "Single Jet Pt", "GeV");

	TH1D* zmass_signal = getHisto("zmassMC_emu_DY",1,0,infile,0);
	TH1D* zmass_dy = getHisto("zmassMC_DY",1,2345,infile,0);
	zmass_signal->Scale(1./zmass_signal->Integral());
	zmass_dy->Scale(1./zmass_dy->Integral());
	TH1D* zmass_ratio = getDataMC(zmass_signal,zmass_dy);
	//drawPlot(zmass_signal,zmass_dy,zmass_ratio,"Custom MC","Official MC","Z mass","GeV");

	TH1D* zpt_signal = getHisto("zptMC_emu_DY",1,0,infile,0);
	TH1D* zpt_dy = getHisto("zptMC_tautau_DY",1,2345,infile,0);
	zpt_signal->Scale(1./zpt_signal->Integral());
	zpt_dy->Scale(1./zpt_dy->Integral());
	TH1D* zpt_ratio = getDataMC(zpt_signal,zpt_dy);
	//for(unsigned i=1;i<=zpt_ratio->GetNbinsX();i++){
	//	std::cout << "bin " << i << ": " << zpt_ratio->GetBinContent(i) << std::endl;
	//}
	//drawPlot(zpt_signal,zpt_dy,zpt_ratio,"Custom MC","Official MC","Z pt","GeV");

	TH1D* zeta_signal = getHisto("zetaMC_emu_DY",1,0,infile,0);
	TH1D* zeta_dy = getHisto("zetaMC_tautau_DY",1,2345,infile,0);
	zeta_signal->Scale(1./zeta_signal->Integral());
	zeta_dy->Scale(1./zeta_dy->Integral());
	TH1D* zeta_ratio = getDataMC(zeta_signal,zeta_dy);
	//drawPlot(zeta_signal,zeta_dy,zeta_ratio,"Custom MC","Official MC","Z eta","GeV");

	TH1D* njets_signal = getHisto("NJetsMC_emu_DY",mcscale.at(mcscale.size()-1),0,infile,syst.at(syst.size()-1));
	TH1D* njets_dy = getHisto("NJetsMC_tautau_DY",mcscale.at(mcscale.size()-2),2345,infile,syst.at(syst.size()-2));
	njets_signal->Scale(1./njets_signal->Integral());
	njets_dy->Scale(1./njets_dy->Integral());
	TH1D* njets_ratio = getDataMC(njets_signal,njets_dy);
	//drawPlot(njets_signal,njets_dy,njets_ratio,"Custom MC","Official MC","Number of jets","");
*/
}

