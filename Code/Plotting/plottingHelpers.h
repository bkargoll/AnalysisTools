/////////////////////////////////////////////////
//
// plotting helper functions
// to be used in plotting macro
//
/////////////////////////////////////////////////

#include "Riostream.h"

// container to hold all information for one sample
struct sample{
	bool isCombinedSample = false;
	std::vector<sample> sampleList;
	int id = -1;
	double xsec = -1;
	int nEvt = -1;
	int color = -1;
	TString legName = "";
	double mcScale = -1;
	double syst = 0;
};
// create sample struct from raw numbers
sample createSample(TString legName, double xsec, int nEvt, int color){
	sample s;
	s.legName = legName;
	s.xsec = xsec;
	s.nEvt = nEvt;
	s.color = color;
	s.mcScale = lumi * s.xsec / s.nEvt;
	//todo: define syst
	return s;
}
// create struct using dataMCType and SkimSummary info
sample createSample(int dataMcType, TString legName, double xsec, int color, std::map<int,int> nEventsMap){
	sample s;
	s.id = dataMcType;
	s.legName = legName;
	s.xsec = xsec;
	s.nEvt = nEventsMap[s.id];
	s.color = color;
	s.mcScale = lumi * s.xsec / s.nEvt;
	//todo: define syst
	return s;
}
// create slim structed without scaling information
sample createSample(TString legName, int color){
	sample s;
	s.legName = legName;
	s.color = color;
	//todo: define syst
	return s;
}
// combine sample structs into one
sample combineSamples(){
	sample s;
	s.isCombinedSample = true;

}

// function to read number of events from SkimSummary
std::map<int, int> readSkimSummary(TString skimsummaryFile){
	ifstream input_file;
	char *file_=(char*)skimsummaryFile.Data();
	input_file.open(file_, std::ios::in);
	if (!(input_file)){
		std::cout << "\nERROR: Opening SkimSummary file "<< skimsummaryFile <<" fhas failed.\n" << std::endl;
		return false;
	}
	std::cout << "\nOpened SkimSummary SkimEff file: "<< skimsummaryFile <<".\n" << std::endl;
	std::string s;
	int a=0;
	std::map<int, int> nEvtMap;
	while(getline(input_file, s)){
		a++;
		if(a>1000) break;
		std::stringstream line(s);
		TString tmp;
		int id;
		float nevents;
		float neventserr;
		float nevents_sel;
		float neventserr_sel;
		float noweight;
		float noweight_sel;
		line >> tmp >> id >> tmp >> nevents >> tmp >> neventserr >> tmp >> nevents_sel >> tmp >> neventserr_sel >> tmp >> noweight >> tmp >> noweight_sel;
		// only consider main data-mc-types, thus removing all but 2 last digits
		nEvtMap[id%100] += nevents;
	}
	return nEvtMap;
}

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
