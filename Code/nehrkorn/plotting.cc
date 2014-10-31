#include "TFile.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TLegend.h"

#include <vector>

bool testPlotting = false;
bool dym50 = true;
bool signaltop = false;

void plotting(){
	gROOT->LoadMacro("tdrstyle.C");
	setTDRStyle();
	gROOT->LoadMacro("CMS_lumi.C");
	writeExtraText = true;
	gStyle->SetOptStat(0);
	
	bool verbose = true;

	// enter filename here
	if(!dym50)TFile* infile = new TFile("/user/nehrkorn/analysis_new.root");
	if(dym50)TFile* infile = new TFile("/user/nehrkorn/analysis_dym50_jeteta52.root");//m60120.root");
	// define crosssections and lumi
	double lumi = 19712.;
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
	if(!dym50)mcscale.push_back(lumi*xdyee/ndyee);
	if(!dym50)mcscale.push_back(lumi*xdymumu/ndymumu);
	if(!dym50)mcscale.push_back(lumi*xdytautau/ndytautau);
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
	if(!dym50)syst.push_back(TMath::Sqrt(TMath::Power(0.026,2)+TMath::Power(0.055,2))); // dy_ee
	if(!dym50)syst.push_back(TMath::Sqrt(TMath::Power(0.026,2)+TMath::Power(0.055,2))); // dy_mumu
	if(!dym50)syst.push_back(TMath::Sqrt(TMath::Power(0.026,2)+TMath::Power(0.055,2))); // dy_tautau
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

/////////////////////////////////////////////////
//
// helper functions
//
/////////////////////////////////////////////////

TH1D* getHisto(TString name, TFile* file){
	TString pre = "ztoemu_default_";
	TH1D* hist = ((TH1D*)file->Get(pre+name));
	return hist;
}

TH1D* getHisto(TString name, double scale, int color, TFile* file){
	TString pre = "ztoemu_default_";
	//TString pre = "ztoemu_skim_default_";
	TH1D* hist = ((TH1D*)file->Get(pre+name));
	hist->Scale(scale);
	hist->SetFillColor(color);
	hist->SetLineColor(1);
	hist->SetLineWidth(1);
	return hist;
}

TH1D* getHisto(TString name, double scale, int color, TFile* file, double systematic){
	TString pre = "ztoemu_default_";
	//TString pre = "ztoemu_skim_default_";
	TH1D* hist = ((TH1D*)file->Get(pre+name));
	hist->Scale(scale);
	for(unsigned i=1;i<hist->GetNbinsX();i++){
		hist->SetBinError(i,TMath::Sqrt(TMath::Power(hist->GetBinError(i),2)+TMath::Power(systematic*hist->GetBinContent(i),2)));
	}
	hist->SetFillColor(color);
	hist->SetLineColor(1);
	hist->SetLineWidth(1);
	if(name.Contains("DY_emu")){
	//if(name.Contains("emu_DY")){
		hist->SetFillStyle(0);
		hist->SetLineStyle(7);
		hist->SetLineColor(kBlack);
		hist->SetLineWidth(3);
	}
	return hist;
}

std::vector<TH1D*> getHistos(std::vector<TString> names, std::vector<double> scale, std::vector<int> color, TFile* file){
	std::vector<TH1D*> histos;
	for(unsigned i=0;i<names.size();i++){
		histos.push_back(getHisto(names.at(i),scale.at(i),color.at(i),file));
	}
	return histos;
}

std::vector<TH1D*> getHistos(TString name, std::vector<double> scale, std::vector<int> color, TFile* file){
	if(!dym50)TString append[14] = {"MC_QCD","MC_ZZ_4l","MC_ZZ_2l2q","MC_ZZ_2l2nu","MC_WZ_3l1nu","MC_WZ_2l2q","MC_WW_2l2nu","MC_ttbar","MC_tw","MC_tbarw","MC_ee_DY","MC_mumu_DY","MC_tautau_DY","MC_emu_DY"};
	if(dym50)TString append[13] = {"MC_QCD","MC_ZZ_4l","MC_ZZ_2l2q","MC_ZZ_2l2nu","MC_WZ_3l1nu","MC_WZ_2l2q","MC_WW_2l2nu","MC_ttbar","MC_tw","MC_tbarw","MC_DY","MC_tautau_DY","MC_emu_DY"};
	TString histname;
	std::vector<TH1D*> histos;
	for(unsigned i=0;i<scale.size();i++){
		histname = name+append[i];
		histos.push_back(getHisto(histname,scale.at(i),color.at(i),file));
	}
	return histos;
}

std::vector<TH1D*> getHistos(TString name, std::vector<double> scale, std::vector<int> color, TFile* file, std::vector<double> systematics){
	if(dym50){
		TString append[13] = {"MC_QCD","MC_ZZ_4l","MC_ZZ_2l2q","MC_ZZ_2l2nu","MC_WZ_3l1nu","MC_WZ_2l2q","MC_WW_2l2nu","MC_ttbar","MC_tw","MC_tbarw","MC_DY","MC_tautau_DY","MC_emu_DY"};
	}else{
		TString append[14] = {"MC_QCD","MC_ZZ_4l","MC_ZZ_2l2q","MC_ZZ_2l2nu","MC_WZ_3l1nu","MC_WZ_2l2q","MC_WW_2l2nu","MC_ttbar","MC_tw","MC_tbarw","MC_ee_DY","MC_mumu_DY","MC_tautau_DY","MC_emu_DY"};
	}
	TString histname;
	std::vector<TH1D*> histos;
	for(unsigned i=0;i<scale.size();i++){
		histname = name+append[i];
		histos.push_back(getHisto(histname,scale.at(i),color.at(i),file,systematics.at(i)));
	}
	return histos;
}

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
			//hist->SetBinError(i,TMath::Sqrt(pow(dataerror/mc,2)+pow(data*mcerror/pow(mc,2),2)));
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

THStack* produceHistStack(std::vector<TH1D*> histos){
	THStack* stack = new THStack("stack","stack");
	for(unsigned i=0;i<histos.size()-1;i++){ //
		stack->Add(histos.at(i));
	}
	return stack;
}

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

std::vector<TH1D*> produceReducedHistos(std::vector<TH1D*> histos, std::vector<int> colors){
	if(!dym50){
		TH1D* qcd = histos.at(0);
		TH1D* dyt = histos.at(12);
		TH1D* sig = histos.at(13);
		TH1D* top = new TH1D("top","top",histos.at(0)->GetNbinsX(),histos.at(0)->GetXaxis()->GetXmin(),histos.at(0)->GetXaxis()->GetXmax());
		top->Add(histos.at(7));
		top->Add(histos.at(8));
		top->Add(histos.at(9));
		TH1D* ewk = new TH1D("ewk","ewk",histos.at(0)->GetNbinsX(),histos.at(0)->GetXaxis()->GetXmin(),histos.at(0)->GetXaxis()->GetXmax());
		ewk->Add(histos.at(1));
		ewk->Add(histos.at(2));
		ewk->Add(histos.at(3));
		ewk->Add(histos.at(4));
		ewk->Add(histos.at(5));
		ewk->Add(histos.at(6));
		ewk->Add(histos.at(10));
		ewk->Add(histos.at(11));
	}
	if(dym50){
		TH1D* qcd = histos.at(0);
		TH1D* dyt = histos.at(11);
		TH1D* sig = histos.at(12);
		TH1D* top = new TH1D("top","top",histos.at(0)->GetNbinsX(),histos.at(0)->GetXaxis()->GetXmin(),histos.at(0)->GetXaxis()->GetXmax());
		top->Add(histos.at(7));
		top->Add(histos.at(8));
		top->Add(histos.at(9));
		TH1D* ewk = new TH1D("ewk","ewk",histos.at(0)->GetNbinsX(),histos.at(0)->GetXaxis()->GetXmin(),histos.at(0)->GetXaxis()->GetXmax());
		ewk->Add(histos.at(1));
		ewk->Add(histos.at(2));
		ewk->Add(histos.at(3));
		ewk->Add(histos.at(4));
		ewk->Add(histos.at(5));
		ewk->Add(histos.at(6));
		ewk->Add(histos.at(10));
	}
	qcd->SetFillColor(colors.at(0));
	ewk->SetFillColor(colors.at(1));
	top->SetFillColor(colors.at(2));
	dyt->SetFillColor(colors.at(3));
	sig->SetFillColor(colors.at(4));
	ewk->SetLineColor(1);
	std::vector<TH1D*> reducedhistos;
	reducedhistos.push_back(qcd);
	reducedhistos.push_back(ewk);
	reducedhistos.push_back(top);
	reducedhistos.push_back(dyt);
	reducedhistos.push_back(sig);
	return reducedhistos;
}

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

void drawPlot(TH1D* data, std::vector<TH1D*> allhistos, std::vector<int> reducedColors, std::vector<TString> names, TString title, TString unit){
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPalette(1);
	gROOT->ForceStyle(true);

	std::vector<TH1D*> histos = produceReducedHistos(allhistos,reducedColors);
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
		//signal->SetFillColor(10);
		//signal->SetFillStyle(3004);
		//signal->SetLineStyle(9);
		//signal->SetFillStyle(0);
		//signal->SetLineWidth(2);
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

void drawPlot(TH1D* data, std::vector<TH1D*> histos, TH1D* ratio, std::vector<TString> names, TString title, TString unit){
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPalette(1);  
	gROOT->ForceStyle(true);
	
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
	
	if(data->GetMaximum()>=total->GetMaximum()){
		data->GetYaxis()->SetRangeUser(0,data->GetMaximum()*1.2);
	}else{
		data->GetYaxis()->SetRangeUser(0,total->GetMaximum()*1.2);
	}
	data->GetYaxis()->SetLabelSize(0.07);
	data->GetYaxis()->SetTitleSize(0.07);
	data->GetYaxis()->SetTitleOffset(1.15);
	TString ytit = "Events / %.2f ";
	TString yTitle = ytit+unit;
	data->GetYaxis()->SetTitle(Form(yTitle.Data(),data->GetBinWidth(1)));
	data->SetTitle(title);	
	data->SetTitle("CMS preliminary, #sqrt{s}=8 TeV, L=19.7 fb^{-1}");
	data->Draw("E");
	stack->Draw("Histsame");
	total->Draw("E2same");
	int bla = histos.size()-1;
	TH1D* signal = histos.at(bla)->Clone();
	//signal->SetFillColor(10);
	//signal->SetFillStyle(3004);
	//signal->SetLineStyle(9);
	//signal->SetFillStyle(0);
	//signal->SetLineWidth(2);
	//signal->SetLineColor(kBlack);
	signal->Scale(1);
	signal->Draw("Histsame");
	data->Draw("Esame");
	data->Draw("axissame");
	data->SetMinimum(1.001);
	histos.push_back(total);
	names.push_back("Bkg uncertainty");
	TLegend* legend = createLegend(data,histos,names);
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
	ratio->GetXaxis()->Set(data->GetXaxis()->GetNbins(),data->GetXaxis()->GetXmin(),data->GetXaxis()->GetXmax());
	ratio->Draw("E");
	line->Draw("same");
	can->cd();
	can->SetWindowSize(800,800);
	CMS_lumi(Pad1,2,0);
}

TLegend* createLegend(TH1D* histo1, TH1D* histo2, TString name1, TString name2){
	TLegend* legend = new TLegend(0.7,0.77,0.90,0.87);
	legend->SetFillColor(0);
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);
	legend->AddEntry(histo1,name1,"pe");
	legend->AddEntry(histo2,name2,"F");
	return legend;
}

TLegend* createLegend(TH1D* data, TString name){
	TLegend* legend = new TLegend(0.7,0.77,0.90,0.87);
	legend->SetFillColor(0);
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);
	legend->AddEntry(data,name,"F");
	return legend;
}

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

double QuadraticSum(int nval, double values[]){
	double sum = 0.;
	for(unsigned i=0;i<nval;i++){
		sum+=TMath::Power(values[i],2);
	}
	return TMath::Sqrt(sum);
}
