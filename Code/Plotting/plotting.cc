#include "TFile.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"

#include <vector>

bool testPlotting = false;
bool signaltop = false;

bool verbose = true;

// define lumi
double lumi = 19712.;
// enter filename here
TFile* infile = new TFile("/user/nehrkorn/analysis_new.root");
bool isFileLumiScaled = true;

void plotting(){
	gROOT->LoadMacro("tdrstyle.C");
	setTDRStyle();
	gROOT->LoadMacro("CMS_lumi.C");
	writeExtraText = true;
	gStyle->SetOptStat(0);
	
	// define samples to use
	sample s_ggHiggs	= createSample("H_{gg}"					, 632);
	sample s_vbfHiggs	= createSample("H_{VBF}"				, 633);
	sample s_apHiggs	= createSample("H_{AP}"					, 634);
	sample s_qcd		= createSample("QCD"					, 606);
	sample s_ztt		= createSample("Z#rightarrow#tau#tau"	, 400);
	sample s_zll		= createSample("Z#rightarrowll"			, 800);
	sample s_ww			= createSample("WW"						, 417);
	sample s_wz2l2q		= createSample("WZ#rightarrowllqq"		, 418);
	sample s_wz3l1nu	= createSample("WZ#rightarrowlll#nu"	, 419);
	sample s_zz4l		= createSample("ZZ#rightarrowllll"		, 420);
	sample s_zz2l2nu	= createSample("ZZ#rightarrowll#nu#nu"	, 416);
	sample s_zz2l2q		= createSample("ZZ#rightarrowllqq"		, 415);
	sample s_ttbar		= createSample("t#bar{t}"				, 600);
	sample s_tw			= createSample("tW"						, 603);
	sample s_tbarW		= createSample("#bar{t}W"				, 591);
	sample s_Wlnu		= createSample("W#rightarrowl#nu"		, 876);
	sample s_Wtaunu		= createSample("W#rightarrow#tau#nu"	, 874);

	// *** examples for sample creation when scaling manually ***
	// std::map<int, int> MnEvents = readSkimSummary("mySkimSummary.log");
	// sample s_ggHiggs =	createSample(11,"H_{gg}", 1.233664, 632, MnEvents );
	// sample s_vbfHiggs =	createSample(12,"H_{VBF}", 0.0997296, 633, MnEvents );
	// sample s_apHiggs =	createSample(13,"H_{AP}", 0.0771792, 634, MnEvents );

	// combine samples for plotting


	// define which samples to plot
	std::vector<sample> samples;

	// define crosssections and lumi
	double xsignal = 0.057;
	double xdytautau = 1966.7;
	double xdyee = 1966.7;
	double xdymumu = 1966.7;
	double xww = 5.824;
	double xtt = 239.;//245.8;
	double xtw = 11.1;
	double xtbarw = 11.1;
	double xwz3lnu = 1.058;
	double xwz2l2q = 2.207;
	double xzz4l = 0.181;
	double xzz2l2q = 2.502;
	double xzz2l2nu = 0.716;
	double xdytautaum50 = 1177.3;
	double xdyllm50 = 2354.6;

	int nsignal = 100000;
	int ndytautau = 48790562;
	int ndyee = 3297045;
	int ndymumu = 3293740;
	int nww = 1933235;
	int ntt = 21675970;
	int ntw = 497658;
	int ntbarw = 493460;
	int nwz3lnu = 2017254;
	int nwz2l2q = 3212348;
	int nzz4l = 3854021;
	int nzz2l2q = 1577042;
	int nzz2l2nu = 777964;
	int ndytautaum50 = 10152445;
	int ndyllm50 = 20307058;
	
	// define colors
	int csignal = 10;
	int cqcd = 619;
	int cdytautau = 5;
	int cdyll = 53;
	int cdyee = 8;
	int cdymumu = 9;
	int cww = 7;
	int ctt = 2;
	int ctw = 3;
	int ctbarw = 4;
	int cwz3lnu = 20;
	int cwz2l2q = 21;
	int czz4l = 22;
	int czz2l2q = 23;
	int czz2l2nu = 24;
	int cwz = 20;
	int czz = 30;
	
	// define order of backgrounds
	std::vector<TString> legendnames;
	legendnames.push_back("QCD/W+jets");
	legendnames.push_back("ZZ");
	legendnames.push_back("WZ");
	legendnames.push_back("WW");
	legendnames.push_back("t#bar{t}");
	legendnames.push_back("tW");
	legendnames.push_back("#bar{t}W");
	legendnames.push_back("Z#rightarrow#mu#mu /ee");
	legendnames.push_back("Z#rightarrow#tau#tau");
	legendnames.push_back("Z#rightarrow e#mu (Signal)");
	
	std::vector<TString> leg;
	leg.push_back("QCD/W(Z)+jets");
	leg.push_back("electroweak");
	leg.push_back("t#bar{t} + singletop");
	leg.push_back("Z#rightarrow#tau#tau");
	leg.push_back("Z#rightarrow e#mu (Signal)");
	
	std::vector<double> mcscale;
	std::vector<int> colors;
	mcscale.push_back(1);
	mcscale.push_back(lumi*xzz4l/nzz4l);
	mcscale.push_back(lumi*xzz2l2q/nzz2l2q);
	mcscale.push_back(lumi*xzz2l2nu/nzz2l2nu);
	mcscale.push_back(lumi*xwz3lnu/nwz3lnu);
	mcscale.push_back(lumi*xwz2l2q/nwz2l2q);
	mcscale.push_back(lumi*xww/nww);
	mcscale.push_back(lumi*xtt/ntt);
	mcscale.push_back(lumi*xtw/ntw);
	mcscale.push_back(lumi*xtbarw/ntbarw);
	mcscale.push_back(lumi*xdyllm50/ndyllm50);
	mcscale.push_back(lumi*xdytautaum50/ndytautaum50);
	mcscale.push_back(lumi*xsignal/nsignal);

	colors.push_back(cqcd);
	colors.push_back(czz4l);
	colors.push_back(czz2l2q);
	colors.push_back(czz2l2nu);
	colors.push_back(cwz3lnu);
	colors.push_back(cwz2l2q);
	colors.push_back(cww);
	colors.push_back(ctt);
	colors.push_back(ctw);
	colors.push_back(ctbarw);
	colors.push_back(cdyll);
	colors.push_back(cdytautau);
	colors.push_back(csignal);
	
	TColor* darkblue = new TColor(1234,0.,0.328125,0.62109375,"",1.);
	TColor* lightblue = new TColor(2345,0.5546875,0.7265625,0.89453125,"",1.);
	
	std::vector<int> reducedColors;
	reducedColors.push_back(1234);
	reducedColors.push_back(38);
	reducedColors.push_back(2345);
	reducedColors.push_back(18);
	reducedColors.push_back(0);
	
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
	
	// create plots
	if(testPlotting){
		TString plot = "NPV";
		TString unit = "";
		TH1D* datahist = getHisto(plot+"Data",1,1,infile);
		drawPlot(datahist,getHistos(plot,mcscale,colors,infile,syst),reducedColors,leg,"",unit);
	}else{
		const int nplots = 19;
		TString plots[nplots] = {"PtMu","etaMu","PtE","etaE","onejet","met","mtMu","ptbal","invmass_ptbalance_m","NPV","invmass_vetos_m","invmass_jetveto_m","zmass_zoom","nm0_met","nm0_onejet","nm0_mtmu","nm0_ptbalance","Cut_10_Nminus0_ptBalance_","mtmu_phicorr"};
		TString units[nplots] = {"GeV","","GeV","","GeV","GeV","GeV","GeV","GeV","","GeV","GeV","GeV","GeV","GeV","GeV","GeV","GeV","GeV"};
		std::vector<TH1D*> datahists;
		for(unsigned i=0;i<nplots;i++){
			datahists.push_back(getHisto(plots[i]+"Data",1,1,infile));
		}
		for(unsigned i=0;i<nplots;i++){
			drawPlot(datahists.at(i),getHistos(plots[i],mcscale,colors,infile,syst),reducedColors,leg,"",units[i]);
		}
	}

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

}

