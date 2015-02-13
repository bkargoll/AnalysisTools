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
bool signaltop = false;
bool emb = false;
bool generatorcomparison = true;
TString prepend = "ztoemu_default_";
TString signalName = "emu_DY";

void plotting(){
	SetStyle();
	//SetExtraText("Simulation"); // for simulation comparison plots
	
	bool verbose = true;

	// enter filename here
	TString filename;
	if(emb) filename = "/user/nehrkorn/analysis_emb_scaled.root";//dym50_jeteta52.root";
	else if (generatorcomparison) filename = "/user/nehrkorn/analysis_gencomparison.root";
	else filename = "/user/nehrkorn/analysis_dym50_jeteta52.root";
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
	if(dym50)names.push_back("MC_DY");
	if(!dym50)names.push_back("MC_ee_DY");
	if(!dym50)names.push_back("MC_mumu_DY");
	if(emb)names.push_back("MC_tautau_emb");
	if(!emb)names.push_back("MC_tautau_DY");
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
	if(dym50){
		hewk.push_back(10);
		hdyt.push_back(11);
		hsig.push_back(12);
	}else{
		hewk.push_back(10);
		hewk.push_back(11);
		hdyt.push_back(12);
		hsig.push_back(13);
	}
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
	leg.push_back("top");
	leg.push_back("Z#rightarrow#tau#tau");
	leg.push_back("Signal, B(Z#rightarrow e#mu)=1.7#times10^{-6}");
	
	std::vector<double> mcscale;
	std::vector<int> colors;
	if(emb){
		mcscale.push_back(1);
		mcscale.push_back(1);
		mcscale.push_back(1);
		mcscale.push_back(1);
		mcscale.push_back(1);
		mcscale.push_back(1);
		mcscale.push_back(1);
		mcscale.push_back(1);
		mcscale.push_back(1);
		mcscale.push_back(1);
		mcscale.push_back(1);
		mcscale.push_back(1);
		mcscale.push_back(1);
	}else{
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
	}

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
		TString plots[nplots] = {"PtMu","etaMu","PtE","etaE","onejet","met","mtMu","ptbal","invmass_ptbalance_m","NPV","invmass_vetos_m","invmass_jetveto_m","zmass_zoom","NJetsTight","Cut_08_Nminus1_oneJet_","sip","sip_nm0","ptbal_zoom"};
		TString units[nplots] = {"GeV","","GeV","","GeV","GeV","GeV","GeV","GeV","","GeV","GeV","GeV","","GeV","","","GeV"};
		std::vector<TH1D*> datahists;
		for(unsigned i=0;i<nplots;i++){
			datahists.push_back(getHisto(plots[i]+"Data",1,1,infile));
		}
		for(unsigned i=0;i<nplots;i++){
			std::cout << "Plot: " << plots[i] << std::endl;
			drawPlot(datahists.at(i),getHistos(plots[i],names,mcscale,colors,infile,syst),histpositions,histnames,reducedColors,leg,"",units[i]);
		}
	}

	if(generatorcomparison){
		const int ncomp = 12;
		TString cplot[ncomp] = {"zpt","zeta","zmass","znjets","zjetpt","zmet","zmtlead","zmttrail","znjets_rec","zjetpt_rec","zleadpt","ztrailpt"};
		TString num[ncomp] = {"MC_emu_DY","MC_emu_DY","MC_emu_DY","MC_emu_DY","MC_emu_DY","MC_emu_DY","MC_emu_DY","MC_emu_DY","MC_emu_DY","MC_emu_DY","MC_emu_DY","MC_emu_DY"};
		//TString denom[ncomp] = {"MC_ee_DY","MC_ee_DY","MC_ee_DY","MC_ee_DY","MC_ee_DY","MC_ee_DY","MC_ee_DY","MC_ee_DY","MC_ee_DY","MC_ee_DY","MC_ee_DY","MC_ee_DY"};
		TString denom[ncomp] = {"MC_mumu_DY","MC_mumu_DY","MC_mumu_DY","MC_mumu_DY","MC_mumu_DY","MC_mumu_DY","MC_mumu_DY","MC_mumu_DY","MC_mumu_DY","MC_mumu_DY","MC_mumu_DY","MC_mumu_DY"};
		TString ctitle[ncomp] = {};
		TString cunits[ncomp] = {"GeV","","GeV","","GeV","GeV","GeV","GeV","","GeV","GeV","GeV"};
		for(unsigned i=0;i<ncomp;i++){
			DrawComparison(cplot[i]+num[i],cplot[i]+denom[i],1.,1.,0.,0.,infile,ctitle[i],cunits[i]);
		}
	}

	TH1D* reddata = getHisto("invmass_ptbalance_widerangeData",1,1,infile);
	std::vector<TH1D*> gethists = getHistos("invmass_ptbalance_m",names,mcscale,colors,infile,syst);
	std::vector<TH1D*> redhists = produceReducedHistos(getHistos("invmass_ptbalance_widerange",names,mcscale,colors,infile,syst),histpositions,histnames,reducedColors);
	TFile* outfile = new TFile("invariant_mass.root","RECREATE");
	reddata->Write();
	for(unsigned i=0; i<redhists.size();i++){
		redhists.at(i)->Write();
	}
	for(unsigned i=0;i<gethists.size();i++){
		gethists.at(i)->Write();
	}
	outfile->Close();

	/*TH1D* onejetdata = getHisto("onejetData",1,1,infile);
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
	testcan->SetLeftMargin(0.15);
	testcan->SetRightMargin(0.05);
	testcan->SetBottomMargin(0.12);
	testcan->cd();
	onejetdata->Draw();
	stack->Draw("hist");
	totalerr->Draw("e2same");
	TLegend* onejetlegend = createLegend(onejetdata,onejethistos,leg);
	onejetlegend->AddEntry(totalerr,"#pm 1 #sigma_{JEC+JER}","F");
	onejetlegend->Draw("same");
	CMS_lumi(testcan,2,0);

	TH1D* totalerr2 = (TH1D*)totalerr->Clone();
	for(unsigned i=1; i<=totalerr2->GetNbinsX(); i++){
		double ibincontent = totalerr2->GetBinContent(i);
		double bincontent = ibincontent;
		for(unsigned j=0; j<onejethistos.size(); j++){
			bincontent -= onejethistos.at(j)->GetBinContent(i);
		}
		if(ibincontent>0){
			totalerr2->SetBinContent(i,1);
			totalerr2->SetBinError(i,bincontent/ibincontent);
		}else{
			totalerr2->SetBinContent(i,1);
			totalerr2->SetBinError(i,0);
		}
	}
	TLine* line2 = new TLine(0,1,200,1);
	TLegend* leg2 = new TLegend(0.73,0.37,0.93,0.87);
	leg2->SetFillColor(0);
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);
	leg2->AddEntry(totalerr2,"#pm1#sigma_{JEC+JER}","F");
	TCanvas* can2 = new TCanvas();
	can2->SetWindowSize(800,800);
	can2->SetLeftMargin(0.15);
	can2->SetRightMargin(0.05);
	can2->SetBottomMargin(0.12);
	can2->cd();
	totalerr2->Draw("e2");
	line2->Draw("same");
	leg2->Draw("same");
	CMS_lumi(can2,2,0);*/

	/*TH1D* zptw_signal = getHisto("zptMC_emu_DY",1,0,infile,0);
	TH1D* zptw_dy = getHisto("zptMC_tautau_DY",1,2345,infile,0);
	const int n = 20;
	double ynlo[n] = {37.9705,198.21,101.814,62.5264,41.8615,29.5613,21.9275,16.7949,13.1418,10.4161,8.3696,6.80754,5.51898,4.57045,3.78405,3.16415,2.64978,2.22681,1.87878,1.61495};
	double ynloe[n] = {0.473081,0.232424,0.0941277,0.053511,0.0364891,0.0278376,0.0225998,0.0187661,0.0162394,0.0145898,0.0130264,0.0119675,0.0108968,0.0100309,0.00930257,0.00889576,0.00811673,0.00781533,0.00709309,0.00794441};
	//double ynlo[n] = {73.6553,197.43,103.449,63.321,42.3373,30.0333,21.8391,16.6274,12.9576,10.2494,8.28162,6.70545,5.41327,4.47742,3.70678,3.12096,2.58732,2.19813,1.84185,1.59675};
	//double ynloe[n] = {0.504075,0.238295,0.0963858,0.0544992,0.0368917,0.0286056,0.0223235,0.0187628,0.0162989,0.0142999,0.012791,0.0119629,0.011009,0.00990161,0.00927022,0.00867596,0.00798123,0.00749499,0.00710057,0.00806605};
	//double ynlo[n] = {155.67941,181.61194,83.011525,48.252917,31.402339,21.780034,15.632109,11.795680,9.1478508,7.2054868,5.7628290,4.6664257,3.7911657,3.1082682,2.5719705,2.1375079,1.7863222,1.4998253,1.2645822,1.0717867};
	//double ynloe[n] = {0.13140797,0.74189260E-01,0.36551345E-01,0.22589270E-01,0.16174990E-01,0.12449918E-01,0.10051555E-01,0.84226652E-02,0.72210897E-02,0.62163384E-02,0.54222855E-02,0.47661320E-02,0.41873075E-02,0.37040839E-02,0.32836813E-02,0.29176814E-02,0.25966418E-02,0.23202844E-02,0.20788599E-02,0.19005943E-02};
	TH1D* zptnnlo = new TH1D("zptnnlo","Z p_{T} from FEWZ 3.1 (NNLO);p_{T}^{Z} / GeV; Events / 5 GeV",20,0.,100.);
	zptnnlo->Sumw2();
	for(unsigned i=0;i<zptnnlo->GetNbinsX();i++){
		zptnnlo->SetBinContent(i+1,ynlo[i]);
		zptnnlo->SetBinError(i+1,ynloe[i]);
	}
	if(zptnnlo->Integral()>0) zptnnlo->Scale(1./zptnnlo->Integral());
	zptw_signal->Scale(1./zptw_signal->Integral());
	zptw_dy->Scale(1./zptw_dy->Integral());
	TH1D* zptw_sigratio = getDataMC(zptw_signal,zptnnlo);
	TH1D* zptw_mcratio = getDataMC(zptw_dy,zptnnlo);
	drawPlot(zptw_signal,zptnnlo,zptw_sigratio,"Pythia MC","FEWZ 3.1 (NNLO)","Z pt","GeV");
	drawPlot(zptw_dy,zptnnlo,zptw_mcratio,"Madgraph MC","FEWZ 3.1 (NNLO)","Z pt","GeV");
	std::cout << "### Madgraph vs. FEWZ 3.1 ###" << std::endl;
	for(unsigned i=1;i<=zptw_mcratio->GetNbinsX();i++){
		std::cout << "i = " << i << ", mc/nnlo = " << zptw_mcratio->GetBinContent(i) << std::endl;
	}
	std::cout << "### Pythia vs. FEWZ 3.1 ###" << std::endl;
	for(unsigned i=1;i<=zptw_sigratio->GetNbinsX();i++){
		std::cout << "i = " << i << ", mc/nnlo = " << zptw_sigratio->GetBinContent(i) << std::endl;
	}*/

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

void DrawComparison(TString a, TString b, double a_scale, double b_scale, double a_syst, double b_syst, TFile* file, TString title, TString unit){
	TH1D* signalhist = getHisto(a,a_scale,0,file,a_syst);
	TH1D* backgroundhist = getHisto(b,b_scale,2345,file,b_syst);
	signalhist->Scale(1./signalhist->Integral());
	backgroundhist->Scale(1./backgroundhist->Integral());
	TH1D* sb_ratio = getDataMC(signalhist,backgroundhist);
	drawPlot(signalhist,backgroundhist,sb_ratio,"Custom MC","Official MC",title,unit);
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
	signal->SetLineColor(kRed);
	signal->SetLineWidth(3);
	if(signaltop)stack->Add(signal);
	stack->Draw("Histsame");
	total->Draw("E2same");
	if(!signaltop){
		signal->SetLineColor(kRed);
		signal->SetLineWidth(3);
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
			histos.at(i)->SetLineColor(kRed);
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
