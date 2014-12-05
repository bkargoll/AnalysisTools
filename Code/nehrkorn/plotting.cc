#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TString.h"

#include <vector>

/*
 * !!! Everything except for the helper functions needs to be customized for each user !!!
 */
bool testPlotting = true;
bool dym50 = true;
bool signaltop = true;
TString prepend = "ztoemu_default_";
TString signalName = "emu_DY";

void plotting(){
	SetStyle();
	//SetExtraText("Simulation"); // for simulation comparison plots
	
	bool verbose = true;

	// enter filename here
	TString filename = "/user/nehrkorn/analysis_dym50_jeteta52.root";
	TFile* infile = new TFile(filename);
	TFile* upfile = new TFile("/user/nehrkorn/jerjecup.root");
	TFile* downfile = new TFile("/user/nehrkorn/jerjecdown.root");

	// define crosssections and lumi
	double lumi = 19712.;
	double xsignal = 0.057;
	double xdytautaum20 = 1966.7;
	double xdyeem20 = 1966.7;
	double xdymumum20 = 1966.7;
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
	int ndytautaum20 = 48790562;
	int ndyeem20 = 3297045;
	int ndymumum20 = 3293740;
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
	std::vector<TString> names;
	names.push_back("MC_QCD");
	names.push_back("MC_ZZ_4l");
	names.push_back("MC_ZZ_2l2q");
	names.push_back("MC_ZZ_2l2nu");
	names.push_back("MC_WZ_3l1nu");
	names.push_back("MC_WZ_2l2q");
	names.push_back("MC_WW_2l2nu");
	names.push_back("MC_ttbar");
	names.push_back("MC_tw");
	names.push_back("MC_tbarw");
	names.push_back("MC_DY");
	names.push_back("MC_tautau_DY");
	names.push_back("MC_emu_DY");
	
	// vectors necessary for reduced histograms
	std::vector<int> hqcd;
	std::vector<int> htop;
	std::vector<int> hewk;
	std::vector<int> hdyt;
	std::vector<int> hsig;
	hqcd.push_back(0);
	htop.push_back(7);
	htop.push_back(8);
	htop.push_back(9);
	hewk.push_back(1);
	hewk.push_back(2);
	hewk.push_back(3);
	hewk.push_back(4);
	hewk.push_back(5);
	hewk.push_back(6);
	hewk.push_back(10);
	hdyt.push_back(11);
	hsig.push_back(12);
	std::vector<std::vector<int>> histpositions;
	histpositions.push_back(hqcd);
	histpositions.push_back(hewk);
	histpositions.push_back(htop);
	histpositions.push_back(hdyt);
	histpositions.push_back(hsig);
	std::vector<TString> histnames;
	histnames.push_back("qcd");
	histnames.push_back("ewk");
	histnames.push_back("top");
	histnames.push_back("dyt");
	histnames.push_back("sig");

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
	if(!dym50)mcscale.push_back(lumi*xdyeem20/ndyeem20);
	if(!dym50)mcscale.push_back(lumi*xdymumum20/ndymumum20);
	if(!dym50)mcscale.push_back(lumi*xdytautaum20/ndytautaum20);
	if(dym50)mcscale.push_back(lumi*xdyllm50/ndyllm50);
	if(dym50)mcscale.push_back(lumi*xdytautaum50/ndytautaum50);
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
	if(!dym50)colors.push_back(cdyee);
	if(!dym50)colors.push_back(cdymumu);
	if(!dym50)colors.push_back(cdytautau);
	if(dym50)colors.push_back(cdyll);
	if(dym50)colors.push_back(cdytautau);
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
	double qcd[nsyst] = 	{0.026,0.387,0.000,0.000,0.000,0.000};
	double zz4l[nsyst] = 	{0.026,0.150,0.000,0.000,0.038,0.000};
	double zz2l2q[nsyst] = 	{0.026,0.150,0.000,0.000,0.000,0.000};
	double zz2l2nu[nsyst] = {0.026,0.150,0.000,0.000,0.333,0.000};
	double wz3lnu[nsyst] = 	{0.026,0.150,0.007,0.007,0.047,0.003};
	double wz2l2q[nsyst] = 	{0.026,0.150,0.000,0.000,0.000,0.000};
	double ww2l2nu[nsyst] = {0.026,0.150,0.008,0.007,0.006,0.005};
	double ttbar[nsyst] = 	{0.026,0.047,0.009,0.008,0.054,0.006};
	double tw[nsyst] = 		{0.026,0.090,0.008,0.007,0.264,0.005};
	double tbarw[nsyst] = 	{0.026,0.090,0.000,0.000,0.000,0.000};
	double dyll[nsyst] = 	{0.026,0.033,0.010,0.008,0.220,0.005};
	double dytt[nsyst] = 	{0.026,0.033,0.009,0.008,0.049,0.005};
	double dyemu[nsyst+1] = {0.026,0.033,0.008,0.008,0.008,0.005};
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
	if(!dym50)syst.push_back(TMath::Sqrt(TMath::Power(0.026,2)+TMath::Power(0.055,2))); // dy_ee
	if(!dym50)syst.push_back(TMath::Sqrt(TMath::Power(0.026,2)+TMath::Power(0.055,2))); // dy_mumu
	if(!dym50)syst.push_back(TMath::Sqrt(TMath::Power(0.026,2)+TMath::Power(0.055,2))); // dy_tautau
	syst.push_back(QuadraticSum(nsyst+1,dyemu));
	
	// create plots
	if(testPlotting){
		TString plot = "NPV";
		TString unit = "";
		TH1D* datahist = getHisto(plot+"Data",1,1,infile);
		drawPlot(datahist,getHistos(plot,names,mcscale,colors,infile,syst),histpositions,histnames,reducedColors,leg,"",unit);
	}else{
		const int nplots = 18;
		TString plots[nplots] = {"PtMu","etaMu","PtE","etaE","onejet","met","mtMu","ptbal","invmass_ptbalance_m","NPV","invmass_vetos_m","invmass_jetveto_m","zmass_zoom","nm0_met","nm0_onejet","nm0_mtmu","nm0_ptbalance","Cut_10_Nminus0_ptBalance_"};
		TString units[nplots] = {"GeV","","GeV","","GeV","GeV","GeV","GeV","GeV","","GeV","GeV","GeV","GeV","GeV","GeV","GeV","GeV"};
		std::vector<TH1D*> datahists;
		for(unsigned i=0;i<nplots;i++){
			datahists.push_back(getHisto(plots[i]+"Data",1,1,infile));
		}
		for(unsigned i=0;i<nplots;i++){
			drawPlot(datahists.at(i),getHistos(plots[i],names,mcscale,colors,infile,syst),histpositions,histnames,reducedColors,leg,"",units[i]);
		}
	}
	std::vector<TH1D*> blub = getHistos("invmass_zmass",names,mcscale,colors,infile);
	TH1D* tot = produceTotal(blub);
	double toterr(0), totsum(0);
	for(unsigned i=1;i<=tot->GetNbinsX();i++){
		if(tot->GetBinContent(i)>0){
			toterr += pow(tot->GetBinError(i),2);
			totsum += tot->GetBinContent(i);
		}
	}

	TH1D* onejetdata = getHisto("onejetData",1,1,infile);
	std::vector<TH1D*> onejethists = getHistos("onejet",names,mcscale,colors,infile);
	std::vector<TH1D*> onejethistsup = getHistos("onejet",names,mcscale,colors,upfile);
	std::vector<TH1D*> onejethistsdown = getHistos("onejet",names,mcscale,colors,downfile);
	TH1D* totalhist = new TH1D("totalhist","totalhist",onejethists.at(0)->GetNbinsX(),onejethists.at(0)->GetXaxis()->GetXmin(),onejethists.at(0)->GetXaxis()->GetXmax());
	totalhist->Sumw2();
	TH1D* totalhistup = new TH1D("totalhistup","totalhistup",onejethistsup.at(0)->GetNbinsX(),onejethistsup.at(0)->GetXaxis()->GetXmin(),onejethistsup.at(0)->GetXaxis()->GetXmax());
	totalhistup->Sumw2();
	TH1D* totalhistdown = new TH1D("totalhistdown","totalhistdown",onejethistsdown.at(0)->GetNbinsX(),onejethistsdown.at(0)->GetXaxis()->GetXmin(),onejethistsdown.at(0)->GetXaxis()->GetXmax());
	totalhistdown->Sumw2();
	for(unsigned i=0;i<onejethists.size()-1;i++){
		totalhist->Add(onejethists.at(i));
		totalhistup->Add(onejethistsup.at(i));
		totalhistdown->Add(onejethistsdown.at(i));
	}
	for(unsigned i=1;i<=onejethists.at(0)->GetNbinsX();i++){
		totalhist->SetBinError(i,fabs(totalhistup->GetBinContent(i)-totalhistdown->GetBinContent(i))/2);
	}
	TH1D* totalerr = (TH1D*)totalhist->Clone();
	totalerr->SetFillStyle(3001);
	totalerr->SetFillColor(kBlack);
	totalerr->SetLineColor(18);
	totalerr->SetMarkerColor(1);
	totalerr->SetMarkerSize(0.001);
	std::vector<TH1D*> onejethistos = produceReducedHistos(onejethists,histpositions,histnames,reducedColors);
	THStack* stack = produceHistStack(onejethistos);
	TCanvas* testcan = new TCanvas();
	testcan->SetWindowSize(800,800);
	testcan->cd();
	//onejetdata->Draw("");
	//stack->Draw("hist");
	//totalerr->Draw("e2same");


	TH1D* onejet_signal = getHisto("onejetMC_emu_DY",mcscale.at(mcscale.size()-1),0,infile,syst.at(syst.size()-1));
	TH1D* onejet_dy = getHisto("onejetMC_tautau_DY",mcscale.at(mcscale.size()-2),2345,infile,syst.at(syst.size()-2));
	onejet_signal->Scale(1./onejet_signal->Integral());
	onejet_dy->Scale(1./onejet_dy->Integral());
	TH1D* onejet_ratio = getDataMC(onejet_signal,onejet_dy);
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

	TH1D* zptw_signal = getHisto("zpt_weirdbinsMC_emu_DY",1,0,infile,0);
	TH1D* zptw_dy = getHisto("zpt_weirdbinsMC_tautau_DY",1,2345,infile,0);
	zptw_signal->Scale(1./zptw_signal->Integral(),"width");
	zptw_dy->Scale(1./zptw_dy->Integral(),"width");
	TH1D* zptw_ratio = getDataMC(zptw_signal,zptw_dy);
	//drawPlot(zptw_signal,zptw_dy,zptw_ratio,"Custom MC","Official MC","Z pt","GeV");

	TH1D* zeta_signal = getHisto("zetaMC_emu_DY",1,0,infile,0);
	TH1D* zeta_dy = getHisto("zetaMC_tautau_DY",1,2345,infile,0);
	zeta_signal->Scale(1./zeta_signal->Integral());
	zeta_dy->Scale(1./zeta_dy->Integral());
	TH1D* zeta_ratio = getDataMC(zeta_signal,zeta_dy);
	//drawPlot(zeta_signal,zeta_dy,zeta_ratio,"Custom MC","Official MC","Z eta","GeV");

}

/////////////////////////////////////////////////
//
// helper functions
//
/////////////////////////////////////////////////

void SetStyle(){
	gROOT->LoadMacro("tdrstyle.C");
	setTDRStyle();
	gROOT->LoadMacro("CMS_lumi.C");
	writeExtraText = true;
	gStyle->SetOptStat(0);
}

void SetExtraText(TString extraText_){
	extraText = extraText_;
}

// Returns scaled and colored histogram from file
TH1D* getHisto(TString name, double scale, int color, TFile* file){
	TH1D* hist = ((TH1D*)file->Get(prepend+name));
	hist->Scale(scale);
	hist->SetFillColor(color);
	hist->SetLineColor(1);
	hist->SetLineWidth(1);
	return hist;
}

// Returns scaled and colored histogram from file with systematic added to statistical uncertainty
// Signal histogram gets different fill style, color and line width
TH1D* getHisto(TString name, double scale, int color, TFile* file, double systematic){
	TH1D* hist = ((TH1D*)file->Get(prepend+name));
	hist->Scale(scale);
	for(unsigned i=1;i<hist->GetNbinsX();i++){
		hist->SetBinError(i,TMath::Sqrt(TMath::Power(hist->GetBinError(i),2)+TMath::Power(systematic*hist->GetBinContent(i),2)));
	}
	hist->SetFillColor(color);
	hist->SetLineColor(1);
	hist->SetLineWidth(1);
	if(name.Contains(signalName)){
		hist->SetFillStyle(0);
		hist->SetLineStyle(7);
		hist->SetLineColor(kBlack);
		hist->SetLineWidth(3);
	}
	return hist;
}

// Returns vector of histograms scaled and colored from file. You can choose a different color for each histogram
// 'append' is a vector of strings containing the different background names you want to see, e.g, 'MC_DY'
std::vector<TH1D*> getHistos(TString name, std::vector<TString> append, std::vector<double> scale, std::vector<int> color, TFile* file){
	TString histname;
	std::vector<TH1D*> histos;
	for(unsigned i=0;i<scale.size();i++){
		histname = name+append.at(i);
		histos.push_back(getHisto(histname,scale.at(i),color.at(i),file));
	}
	return histos;
}

// Same as above but with systematics added to statistical uncertainties
std::vector<TH1D*> getHistos(TString name, std::vector<TString> append, std::vector<double> scale, std::vector<int> color, TFile* file, std::vector<double> systematics){
	TString histname;
	std::vector<TH1D*> histos;
	for(unsigned i=0;i<scale.size();i++){
		histname = name+append.at(i);
		histos.push_back(getHisto(histname,scale.at(i),color.at(i),file,systematics.at(i)));
	}
	return histos;
}

// Returns ratio histogram of two histograms (typically data/MC)
TH1D* getDataMC(TH1D* datahist, TH1D* MChist){
	int nbins = datahist->GetNbinsX();
	double xlow = datahist->GetXaxis()->GetXmin();
	double xhigh = datahist->GetXaxis()->GetXmax();
	TH1D* hist = new TH1D("hist",";;Data/MC",nbins,xlow,xhigh);
	
	for(unsigned int i=1;i<=nbins;i++){
		double data = datahist->GetBinContent(i);
		double dataerror = datahist->GetBinError(i);
		double mc = MChist->GetBinContent(i);
		double mcerror = MChist->GetBinError(i);
		if(mc>0){
			hist->SetBinContent(i,data/mc);
			hist->SetBinError(i,TMath::Sqrt(pow(dataerror/mc,2)+pow(data*mcerror/pow(mc,2),2)));
		}
	}
	double min = 999;
	for(unsigned i=1;i<nbins;i++){
		if(hist->GetBinContent(i)>0){
			if(hist->GetBinContent(i)<min)min=hist->GetBinContent(i);
		}
	}
	double max = hist->GetMaximum();
	if(max>2)max=1.5;
	if(min<0.5)min=0.5;
	hist->GetYaxis()->SetRangeUser(min-0.2,max+0.2);
	hist->GetXaxis()->SetTitle(datahist->GetXaxis()->GetTitle());
	return hist;
}

// Same as above but with ratio between one histogram and the sum of all backgrounds
TH1D* getDataMC(TH1D* datahist, std::vector<TH1D*> MChists){
	int nbins = datahist->GetNbinsX();
	double xlow = datahist->GetXaxis()->GetXmin();
	double xhigh = datahist->GetXaxis()->GetXmax();
	TH1D* hist = new TH1D("hist",";;Data/MC",nbins,xlow,xhigh);
	
	for(unsigned int i=1;i<=nbins;i++){
		double data = datahist->GetBinContent(i);
		double dataerror = datahist->GetBinError(i);
		double mc = 0;
		double mcerr = 0;
		for(unsigned j=0;j<MChists.size()-1;j++){
			mc+=MChists.at(j)->GetBinContent(i);
			mcerr+=TMath::Power(MChists.at(j)->GetBinError(i),2);
		}
		double mcerror = TMath::Sqrt(mcerr);
		if(mc>0){
			hist->SetBinContent(i,data/mc);
			hist->SetBinError(i,dataerror/mc);
		}
	}
	double min = 999;
	for(unsigned i=1;i<nbins;i++){
		if(hist->GetBinContent(i)>0){
			if(hist->GetBinContent(i)<min)min=hist->GetBinContent(i);
		}
	}
	double max = hist->GetMaximum();
	if(max>2)max=1.5;
	if(min<0.5)min=0.5;
	//hist->GetYaxis()->SetRangeUser(min-0.2,max+0.2);
	hist->GetYaxis()->SetRangeUser(0.7,1.3);
	hist->GetXaxis()->SetTitle(datahist->GetXaxis()->GetTitle());
	return hist;
}

// Returns a stack of all histograms except the last one in the vector (signal)
THStack* produceHistStack(std::vector<TH1D*> histos){
	THStack* stack = new THStack("stack","stack");
	for(unsigned i=0;i<histos.size()-1;i++){ //
		stack->Add(histos.at(i));
	}
	return stack;
}

// Returns sum of all histograms given (nice for uncertainty band of background) except the last one in vector (signal)
TH1D* produceTotal(std::vector<TH1D*> histos){
	TH1D* total = new TH1D("total","total",histos.at(0)->GetNbinsX(),histos.at(0)->GetXaxis()->GetXmin(),histos.at(0)->GetXaxis()->GetXmax());
	total->Sumw2();
	for(unsigned i=0;i<histos.size()-1;i++){
		total->Add(histos.at(i));
	}
	total->SetFillStyle(3013);
	total->SetFillColor(kGray+1);
	total->SetLineColor(18);
	total->SetMarkerColor(1);
	total->SetMarkerSize(0.001);
	return total;
}

// Produces vector of merged histograms. Useful if you want to combine several histograms to one (e.g, ZZ+WW backgrounds to EWK)
// 'histpositions' is a vector of vector of the positions of the histograms to be combined (e.g., ZZ is the second histogram, WW the fifth. Then the vector should contain a vector of the number 1 and 4)
// 'names' is a vector of TStrings and should contain the names of your combined 'categories' (e.g., EWK for ZZ and WW, top for ttbar, tW and tbarW, etc)
std::vector<TH1D*> produceReducedHistos(std::vector<TH1D*> histos, std::vector<std::vector<int> > histpositions, std::vector<TString> names, std::vector<int> colors){
	std::vector<TH1D*> reducedhistos;
	if(histos.size()==0) return reducedhistos;
	for(unsigned i=0;i<histpositions.size();i++){
		TH1D* hist = new TH1D(names.at(i),names.at(i),histos.at(0)->GetNbinsX(),histos.at(0)->GetXaxis()->GetXmin(),histos.at(0)->GetXaxis()->GetXmax());
		for(unsigned j=0;j<histpositions.at(i).size();j++){
			hist->Add(histos.at(histpositions.at(i).at(j)));
		}
		hist->SetFillColor(colors.at(i));
		reducedhistos.push_back(hist);
	}
	return reducedhistos;
}

// draw single histogram with name 'name', title 'title' and the unit on the y-axis 'unit'
void drawPlot(TH1D* data, TString name, TString title, TString unit){
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPalette(1);  
	gROOT->ForceStyle(true);
	
	TCanvas* can = new TCanvas();
	TPad* Pad1 = new TPad("Pad1","Pad1",0.,0.,1.,1.);
	Pad1->SetTopMargin(0.07);
	Pad1->SetLeftMargin(0.15);
	Pad1->SetRightMargin(0.05);
	Pad1->SetBottomMargin(0.15);
	Pad1->Draw();
	Pad1->cd();
	
	data->GetYaxis()->SetTitleOffset(1.5);
	data->SetFillColor(2345);
	TString ytit = "Events / %.2f ";
	TString yTitle = ytit+unit;
	data->GetYaxis()->SetTitle(Form(yTitle.Data(),data->GetBinWidth(1)));
	data->SetTitle(title);
	data->SetTitle("CMS preliminary, #sqrt{s}=8 TeV, L=19.7 fb^{-1}");
	data->Draw("hist");
	TLegend* legend = createLegend(data,name);
	legend->Draw("same");
	can->cd();
	can->SetWindowSize(800,800);
	CMS_lumi(Pad1,2,0);
}

// Draws plot with data, all backgrounds and combined histograms for the backgrounds
// A ratio plot is automatically created and plotted below the histograms
// Also, an uncertainty band from all backgrounds is computed and drawn
void drawPlot(TH1D* data, std::vector<TH1D*> allhistos, std::vector<std::vector<int> > histpositions, std::vector<TString> histnames, std::vector<int> reducedColors, std::vector<TString> names, TString title, TString unit){
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPalette(1);
	gROOT->ForceStyle(true);

	std::vector<TH1D*> histos = produceReducedHistos(allhistos,histpositions,histnames,reducedColors);
	TH1D* ratio = getDataMC(data,allhistos);

	TCanvas* can = new TCanvas();
	THStack* stack = produceHistStack(histos);
	TH1D* total = produceTotal(histos);
	TLine* line = new TLine(data->GetXaxis()->GetXmin(),1,data->GetXaxis()->GetXmax(),1);
	TPad* Pad1 = new TPad("Pad1","Pad1",0.,0.3,1.,1.);
	Pad1->SetTopMargin(0.07);
	Pad1->SetLeftMargin(0.15);
	Pad1->SetRightMargin(0.05);
	Pad1->SetBottomMargin(0);
	Pad1->Draw();
	Pad1->cd();
	CMS_lumi(Pad1,2,0);
	Pad1->Update();

	if(data->GetMaximum()>=total->GetMaximum()){
		data->GetYaxis()->SetRangeUser(0,data->GetMaximum()*1.2);
	}else{
		data->GetYaxis()->SetRangeUser(0,total->GetMaximum()*1.2);
	}
	data->GetYaxis()->SetLabelSize(0.07);
	data->GetYaxis()->SetTitleSize(0.07);
	data->GetYaxis()->SetTitleOffset(1.15);
	data->SetMarkerColor(kBlack);
	data->SetMarkerStyle(20);
	TString ytit = "Events / %.2f ";
	TString yTitle = ytit+unit;
	data->GetYaxis()->SetTitle(Form(yTitle.Data(),data->GetBinWidth(1)));
	data->SetTitle(title);
	data->SetTitle("CMS preliminary, #sqrt{s}=8 TeV, L=19.7 fb^{-1}");
	data->Draw("E");
	int signalhist = histos.size()-1;
	TH1D* signal = histos.at(signalhist)->Clone();
	if(signaltop)stack->Add(signal);
	stack->Draw("Histsame");
	total->Draw("E2same");
	if(!signaltop){
		signal->SetLineColor(kBlack);
		signal->Scale(1);
		signal->Draw("Histsame");
	}
	data->Draw("Esame");
	data->Draw("axissame");
	data->SetMinimum(1.001);
	histos.push_back(total);
	names.push_back("Bkg uncertainty");
	TLegend* legend = createLegend(data,histos,names);
	legend->Draw("same");
	can->cd();
	TH1D* ratioband = new TH1D("ratioband","ratioband",data->GetNbinsX(),data->GetXaxis()->GetXmin(),data->GetXaxis()->GetXmax());
	for(unsigned i=1;i<=ratioband->GetNbinsX();i++){
		ratioband->SetBinContent(i,1);
		ratioband->SetBinError(i,total->GetBinError(i)/total->GetBinContent(i));
	}
	ratioband->SetFillStyle(3013);
	ratioband->SetFillColor(kGray+1);
	ratioband->SetLineColor(18);
	ratioband->SetMarkerColor(1);
	ratioband->SetMarkerSize(0.001);
	TPad* Pad2 = new TPad("Pad1","Pad1",0.,0.,1.,0.3);
	Pad2->SetTopMargin(0);
	Pad2->SetLeftMargin(0.15);
	Pad2->SetRightMargin(0.05);
	Pad2->SetBottomMargin(0.4);
	Pad2->SetTickx(kTRUE);
	Pad2->SetGridx();
	Pad2->SetGridy();
	Pad2->Draw();
	Pad2->cd();

	ratio->GetXaxis()->SetTitleSize(0.15);
	ratio->GetXaxis()->SetLabelSize(0.15);
	ratio->GetXaxis()->SetTickLength(0.075);
	ratio->GetYaxis()->SetTitleSize(0.15);
	ratio->GetYaxis()->SetLabelSize(0.15);
	ratio->GetYaxis()->SetTitleOffset(0.35);
	ratio->GetYaxis()->CenterTitle();

	ratio->SetMarkerStyle(20);
	ratio->SetMarkerSize(0.7);
	ratio->SetLineColor(kBlack);

	ratio->GetYaxis()->SetNdivisions(4,5,0,kTRUE);
	ratio->GetXaxis()->Set(data->GetXaxis()->GetNbins(),data->GetXaxis()->GetXmin(),data->GetXaxis()->GetXmax());
	ratio->Draw("E");
	ratioband->Draw("E2same");
	line->Draw("same");
	ratio->Draw("Esame");
	can->cd();
	can->SetWindowSize(800,800);
	CMS_lumi(Pad1,2,0);
}

// Draw plot with two histograms as well as a ratio plot of the two
void drawPlot(TH1D* histo1, TH1D* histo2, TH1D* ratio, TString name1, TString name2, TString title, TString unit){
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPalette(1);  
	gROOT->ForceStyle(true);

	TH1D* total = histo2->Clone();
	total->SetFillStyle(3005);
	total->SetFillColor(1);
	total->SetLineColor(18);
	total->SetMarkerColor(1);
	total->SetMarkerSize(0.001);
	
	TCanvas* can = new TCanvas();
	TLine* line = new TLine(histo1->GetXaxis()->GetXmin(),1,histo1->GetXaxis()->GetXmax(),1);
	TPad* Pad1 = new TPad("Pad1","Pad1",0.,0.3,1.,1.);
	Pad1->SetTopMargin(0.07);
	Pad1->SetLeftMargin(0.15);
	Pad1->SetRightMargin(0.05);
	Pad1->SetBottomMargin(0);
	Pad1->Draw();
	Pad1->cd();
	
	if(histo1->GetMaximum()>=histo2->GetMaximum()){
		histo1->GetYaxis()->SetRangeUser(0,histo1->GetMaximum()*1.2);
	}else{
		histo1->GetYaxis()->SetRangeUser(0,histo2->GetMaximum()*1.2);
	}
	histo1->GetYaxis()->SetLabelSize(0.07);
	histo1->GetYaxis()->SetTitleSize(0.07);
	histo1->GetYaxis()->SetTitleOffset(1.15);
	TString ytit = "Events / %.2f ";
	TString yTitle = ytit+unit;
	histo1->GetYaxis()->SetTitle(Form(yTitle.Data(),histo1->GetBinWidth(1)));
	histo1->SetTitle(title);	
	histo1->SetTitle("CMS simulation, #sqrt{s}=8 TeV");
	histo1->SetMarkerStyle(20);
	histo1->SetMarkerSize(0.7);
	histo1->Draw("E");
	histo2->Draw("Histsame");
	total->Draw("E2same");
	histo1->Draw("Esame");
	histo1->Draw("axissame");
	histo1->SetMinimum(1.001);
	TLegend* legend = createLegend(histo1,histo2,name1,name2);
	legend->Draw("same");
	can->cd();
	TPad* Pad2 = new TPad("Pad1","Pad1",0.,0.,1.,0.3);
	Pad2->SetTopMargin(0);
	Pad2->SetLeftMargin(0.15);
	Pad2->SetRightMargin(0.05);
	Pad2->SetBottomMargin(0.4);
	Pad2->SetTickx(kTRUE);
	Pad2->SetGridx();
	Pad2->SetGridy();
	Pad2->Draw();
	Pad2->cd();
	
	ratio->GetXaxis()->SetTitleSize(0.15);
	ratio->GetXaxis()->SetLabelSize(0.15);
	ratio->GetXaxis()->SetTickLength(0.075);
	ratio->GetYaxis()->SetTitleSize(0.15);
	ratio->GetYaxis()->SetLabelSize(0.15);
	ratio->GetYaxis()->SetTitleOffset(0.35);
	ratio->GetYaxis()->CenterTitle();
	
	ratio->GetYaxis()->SetNdivisions(4,5,0,kTRUE);
	ratio->GetXaxis()->Set(histo1->GetXaxis()->GetNbins(),histo1->GetXaxis()->GetXmin(),histo1->GetXaxis()->GetXmax());
	ratio->SetMarkerStyle(20);
	ratio->SetMarkerSize(0.6);
	ratio->SetLineColor(kBlack);
	ratio->Draw("E");
	line->Draw("same");
	can->cd();
	can->SetWindowSize(800,800);
	CMS_lumi(Pad1,2,0);
}

// creates legend for plot with two histograms
TLegend* createLegend(TH1D* histo1, TH1D* histo2, TString name1, TString name2){
	TLegend* legend = new TLegend(0.7,0.77,0.90,0.87);
	legend->SetFillColor(0);
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);
	legend->AddEntry(histo1,name1,"pe");
	legend->AddEntry(histo2,name2,"F");
	return legend;
}

// creates legend for one histogram
TLegend* createLegend(TH1D* data, TString name){
	TLegend* legend = new TLegend(0.7,0.77,0.90,0.87);
	legend->SetFillColor(0);
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);
	legend->AddEntry(data,name,"F");
	return legend;
}

// creates legend for data histogram and all backgrounds
TLegend* createLegend(TH1D* data, std::vector<TH1D*> histos, std::vector<TString> names){
	TLegend* legend = new TLegend(0.73,0.37,0.93,0.87);
	legend->SetFillColor(0);
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);
	legend->AddEntry(data,"Data","pe");
	for(int i=histos.size()-1;i>=0;i--){
		if(names.at(i).Contains("Signal")){
			legend->AddEntry(histos.at(i),names.at(i),"l");
		}else{
			legend->AddEntry(histos.at(i),names.at(i),"F");
		}
	}
	return legend;
}

// helper function to add stuff (usually uncertainties) in quadrature
double QuadraticSum(int nval, double values[]){
	double sum = 0.;
	for(unsigned i=0;i<nval;i++){
		sum+=TMath::Power(values[i],2);
	}
	return TMath::Sqrt(sum);
}
