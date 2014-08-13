#include "Plots.h"

#include "THStack.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TView.h"
#include "TMath.h"
#include <cmath>
#include <iostream>

int Plots::plotLabel = 0;
TString Plots::File_;
std::vector<TString> Plots::HistogramNames_;

Plots::Plots() :
		doscale(false),
		verbose(false),
		dooneprofile(false) {

}

Plots::~Plots() {
}

void Plots::Set_Plot_Type(TString style, TString label) {
	style.ToLower();
	label.ToLower();
	if (label.Contains("internal"))
		plotLabel = cmsInternal;
	if (label.Contains("private"))
		plotLabel = cmsPrivate;
	if (label.Contains("preliminary"))
		plotLabel = cmsPreliminary;
	if (label.Contains("public"))
		plotLabel = cmsPublic;
	if (style.Contains("style1")) {
		CMSStyle1();
	} else if (style.Contains("style2")) {
		CMSStyle2();
	} else {
		CMSStyle2();
	}
}

void Plots::Plot1D(std::vector<TH1D> histo, std::vector<int> colour, std::vector<TString> legend) {
	std::vector<std::vector<TH1D> > histos;
	histos.push_back(histo);
	Plot1D(histos, colour, legend);
}

void Plots::SaveHistograms(TString File, std::vector<TString> HistogramNames) {
	File_ = File;
	HistogramNames_ = HistogramNames;
}

void Plots::Plot1D(std::vector<std::vector<TH1D> > histo, std::vector<int> colour, std::vector<TString> legend) {

	std::cout << "Plots::Plot1D" << std::endl;
	TCanvas c("c", "c", 200, 10, 750, 750);
	c.Update();
	c.Clear();
	c.Update();
	TLatex latex;
	latex.SetTextSize(0.03);
	latex.SetNDC();

	for (int j = 0; j < histo.size(); j++) {
		if (verbose)
			std::cout << "Plots::Plot1D " << histo.size() << " j= " << j << " color.size()= " << colour.size() << " " << histo.at(j).size() << std::endl;
		if (histo.at(j).size() > 0 && histo.at(j).size() == colour.size()) {
			TLegend leg(0.25, 0.75, 0.9, 0.925);
			leg.SetBorderSize(0);
			leg.SetFillStyle(4000);
			leg.SetFillColor(0);
			leg.SetTextSize(0.03);
			leg.SetMargin(0.15);
			leg.SetNColumns(3);
			leg.SetColumnSeparation(0.05);
			TH1D Total("Total", "Total", histo.at(j).at(0).GetNbinsX(), histo.at(j).at(0).GetXaxis()->GetXmin(), histo.at(j).at(0).GetXaxis()->GetXmax());
			Total.SetFillStyle(3005);
			Total.SetFillColor(1);
			Total.SetLineColor(18);
			Total.SetMarkerColor(1);
			Total.SetMarkerSize(0.001);
			Total.Sumw2();
			TString N = "MCStack";
			N += j;
			N += histo.at(j).at(0).GetName();
			THStack MCHistoStack(N, histo.at(j).at(0).GetName());
			double MC_Integral(0);
			double Data_Integral(0);
			bool cflag = true;
			TString l = "";
			histo.at(j).at(0).SetLineColor(colour.at(0));
			histo.at(j).at(0).SetMarkerColor(colour.at(0));
			if (verbose)
				std::cout << histo.at(j).at(0).GetTitle() << std::endl;
			for (int i = histo.at(j).size() - 1; i >= 0; i--) {
				histo.at(j).at(i).SetLineColor(colour.at(i));
				histo.at(j).at(i).SetFillColor(colour.at(i));
				if (i + 1 < histo.at(j).size()) {
					if (legend.at(i) != l) {
						l = legend.at(i);
						cflag = true;
					}
				}
				if (cflag && histo.at(j).at(i).Integral() > 0) {
					histo.at(j).at(i).SetLineColor(1);
					cflag = false;
				}
			}
			double theIntegral(0);
			for (int i = 0; i < histo.at(j).size(); i++) {
				theIntegral += histo.at(j).at(i).Integral();
				if (i == 0) {
					Data_Integral += histo.at(j).at(0).Integral();
					leg.AddEntry(&histo.at(j).at(i), legend.at(i), "pe");
					theIntegral = 0;
				} else {
					if (colour.size() <= i + 1 && theIntegral >= 0) {
						histo.at(j).at(i).SetLineColor(1);
						leg.AddEntry(&histo.at(j).at(i), legend.at(i), "F");
						theIntegral = 0;
					} else if (colour.size() > i + 1) {
						if (colour.at(i) != colour.at(i + 1) && theIntegral >= 0) {
							histo.at(j).at(i).SetLineColor(1);
							leg.AddEntry(&histo.at(j).at(i), legend.at(i), "F");
							theIntegral = 0;
						}
					}
					if (i != 0) {
						if (histo.at(j).at(i).Integral() > 0) {
							MCHistoStack.Add(&histo.at(j).at(i));
							Total.Add(&histo.at(j).at(i));
						}
					}
					MC_Integral += histo.at(j).at(i).Integral();
					if (verbose)
						std::cout << histo.at(j).at(i).GetTitle() << " " << histo.at(j).at(i).Integral() << " i= " << i << std::endl;
				}
			}
			if (doscale) {
				for (int i = 0; i < histo.at(j).size(); i++) {
					if (i != 0 && Data_Integral > 0 && MC_Integral > 0)
						histo.at(j).at(i).Scale(Data_Integral / MC_Integral);
				}
				if (Data_Integral > 0 && MC_Integral > 0)
					Total.Scale(Data_Integral / MC_Integral);
			}
			if (verbose)
				std::cout << "A" << std::endl;
			for (int l = 0; l < 2; l++) {
				c.Clear();
				c.SetLogy(0);
				if (l == 1) {
					c.Clear();
					c.SetLogy(1);
				}
				if (Data_Integral > 0) {
					double max = histo.at(j).at(0).GetBinContent(histo.at(j).at(0).GetMaximumBin());
					if (l == 1) {
						histo.at(j).at(0).SetMinimum(0.01);
						histo.at(j).at(0).SetMaximum(100 * max);
					} else {
						histo.at(j).at(0).SetMinimum(0);
						histo.at(j).at(0).SetMaximum(1.7 * max);
					}
					histo.at(j).at(0).Draw("E");
					MCHistoStack.Draw("Histsame");
					Total.Draw("E2same");
					histo.at(j).at(0).Draw("Esame");
					histo.at(j).at(0).Draw("axissame");
					TString hTitle = histo.at(j).at(0).GetName();
					if (hTitle.Contains("KSTest")) {
						double KS = histo.at(j).at(0).KolmogorovTest(&Total);
						hTitle = histo.at(j).at(0).GetTitle();
						hTitle += " KS=";
						hTitle += KS;
						latex.DrawLatex(0.225, 0.95, hTitle);
					}
					if (verbose)
						std::cout << "Data: " << histo.at(j).at(0).Integral() << " MC: " << Total.Integral() << std::endl;
				} else {
					MCHistoStack.Draw("Hist");
					double max = Total.GetBinContent(Total.GetMaximumBin());
					if (l == 1) {
						histo.at(j).at(0).SetMinimum(0.01);
						if (max != 0) {
							histo.at(j).at(0).SetMaximum(100 * max);
						} else {
							histo.at(j).at(0).SetMaximum(100);
						}
					} else {
						histo.at(j).at(0).SetMinimum(0);
						histo.at(j).at(0).SetMaximum(1.7 * max);
					}
					histo.at(j).at(0).Draw("E");
					MCHistoStack.Draw("Histsame");
					MCHistoStack.Draw("axissame");
					Total.Draw("E2same");
				}
				leg.Draw();
				c.Update();
				TString name = histo.at(j).at(0).GetName();
				if (l == 1) {
					name += "_log";
				}
				name += "_index_";
				name += j;
				TString EPSName = "EPS/";
				EPSName += name;
				EPSName += ".eps";
				if (l == 0 || !name.Contains("Nminus") || name.Contains("Accumdist")) {
					c.Print(EPSName);
				}
				if (l == 0) {
					for (int a = 0; a < HistogramNames_.size(); a++) {
						if (histo.at(j).at(0).GetName() == HistogramNames_.at(a)) {
							TFile f(File_ + "_EXTRA.root", "RECREATE");
							TString n = HistogramNames_.at(a);
							n.ReplaceAll("Data", "");
							int s = histo.at(j).size() - 1;
							histo.at(j).at(0).Write((histo.at(j).at(0)).GetName());
							histo.at(j).at(s).Write(n + "sig");
							for (int t = 2; t < s; t++) {
								histo.at(j).at(2).Add(&histo.at(j).at(t));
							}
							histo.at(j).at(2).Write(n + "background");
							f.Close();
						}
					}
				}
			}
		}
	}
	std::cout << "Plot1D done " << std::endl;
}

void Plots::Plot1DSignificance(std::vector<TH1D> histo, bool gt, bool lt, std::vector<int> colour, std::vector<TString> legend) {
	// plots significance: assumes the histograms have been normilized
	std::cout << "Plots::Plot1DSignificance" << std::endl;
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.20);
	gStyle->SetPadLeftMargin(0.175);
	gStyle->SetPadRightMargin(0.05);
	TCanvas c("c", "c", 200, 10, 750, 750);
	c.Update();
	c.Clear();
	c.Update();

	histo.at(0).Reset();
	if (histo.size() > 0) {
		for (int b = 1; b <= histo.at(0).GetNbinsX(); b++) {
			if (gt) {
				double sig = histo.at(histo.size() - 1).Integral(b, histo.at(0).GetNbinsX() + 1);
				double bkg(0);
				for (int i = 1; i < histo.size() - 1; i++) {
					bkg += histo.at(i).Integral(b, histo.at(0).GetNbinsX() + 1);
				}
				if (/*1+sig/bkg >0 && */(sig + bkg) != 0) {
					histo.at(0).SetBinContent(b, sig / sqrt(sig + bkg)); //(b,sqrt(2)*sqrt((sig+bkg)*TMath::Log(1+sig/bkg)-sig));//(b,sig/sqrt(sig+bkg));
				} else {
					histo.at(0).SetBinContent(b, 0);
				}
			}
			if (lt) {
				double sig = histo.at(histo.size() - 1).Integral(0, b);
				double bkg(0);
				for (int i = 1; i < histo.size() - 1; i++) {
					bkg += histo.at(i).Integral(0, b);
				}
				if (/*1+sig/bkg >0 && */sig + bkg != 0) {
					histo.at(0).SetBinContent(b, sig / sqrt(sig + bkg)); //(b,sqrt(2)*sqrt((sig+bkg)*TMath::Log(1+sig/bkg)-sig));
				} else {
					histo.at(0).SetBinContent(b, 0);
				}
			}
		}

		histo.at(0).SetYTitle("Signal/#sqrt{Signal+Bkg}"); //("#sqrt{2((S+B)ln(1+S/B)-S)}");//("Signal/#sqrt{Signal+Bkg}");
		histo.at(0).Draw("hist");
		histo.at(0).Draw("axissame");

		TString name = histo.at(0).GetName();
		name += "_sig";
		if (lt)
			name += "lt";
		if (gt)
			name += "gt";
		name += "_index_";
		name += 0;
		TString EPSName = "EPS/";
		EPSName += name;
		EPSName += ".eps";

		if (!name.Contains("Accumdist")) {
			c.Print(EPSName);
			Int_t bmax = histo.at(0).GetMaximumBin();
			if (gt)
				bmax += 1;
			if (verbose)
				std::cout << "Opt. Sig. " << EPSName << " Sig: " << histo.at(0).GetMaximum() << " bin: " << bmax << std::endl;
		}
	}

	std::cout << "Plots::Plot1DSignificance done" << std::endl;
}

void Plots::Plot1Dsigtobkg(std::vector<TH1D> histo, bool gt, bool lt, std::vector<int> colour, std::vector<TString> legend) {
	// plots significance: assumes the histograms have been normilized
	std::cout << "Plots::Plot1Dsigtobkg" << std::endl;
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.20);
	gStyle->SetPadLeftMargin(0.175);
	gStyle->SetPadRightMargin(0.05);
	TCanvas c("c", "c", 200, 10, 750, 750);
	c.Update();
	c.Clear();
	c.Update();

	histo.at(0).Reset();
	if (histo.size() > 2) {
		for (int b = 1; b <= histo.at(0).GetNbinsX(); b++) {
			if (gt) {
				double sig = histo.at(histo.size() - 1).Integral(b, histo.at(0).GetNbinsX() + 1);
				double bkg(0);
				for (int i = 1; i < histo.size() - 1; i++) {
					bkg += histo.at(i).Integral(b, histo.at(0).GetNbinsX() + 1);
				}
				if ((bkg + sig) != 0) {
					histo.at(0).SetBinContent(b, sig / (sig + bkg));
				} else {
					histo.at(0).SetBinContent(b, 0);
				}
			}
			if (lt) {
				double sig = histo.at(histo.size() - 1).Integral(0, b);
				double bkg(0);
				for (int i = 1; i < histo.size() - 1; i++) {
					bkg += histo.at(i).Integral(0, b);
				}
				if (sig + bkg != 0) {
					histo.at(0).SetBinContent(b, sig / (sig + bkg));
				} else {
					histo.at(0).SetBinContent(b, 0);
				}
			}
		}

		for (int b = 1; b <= histo.at(0).GetNbinsX(); b++) {
			double total = histo.at(histo.size() - 1).Integral();
			if (gt) {
				double sig = histo.at(histo.size() - 1).Integral(b, histo.at(0).GetNbinsX() + 1);
				if (total != 0) {
					histo.at(1).SetBinContent(b, sig / total);
				} else {
					histo.at(1).SetBinContent(b, 0);
				}
			}
			if (lt) {
				double sig = histo.at(histo.size() - 1).Integral(0, b);
				if (total != 0) {
					histo.at(1).SetBinContent(b, sig / total);
				} else {
					histo.at(1).SetBinContent(b, 0);
				}
			}
		}

		histo.at(0).SetYTitle("Purity (Black)/Efficiency (Red)");
		histo.at(0).SetMinimum(0.0);
		histo.at(0).SetMaximum(1.0);
		histo.at(0).Draw("hist");
		histo.at(0).Draw("axissame");
		histo.at(1).SetLineColor(2);
		histo.at(1).Draw("histsame");

		TString name = histo.at(0).GetName();
		name += "_sigtobkg";
		if (lt)
			name += "lt";
		if (gt)
			name += "gt";
		name += "_index_";
		name += 0;
		TString EPSName = "EPS/";
		EPSName += name;
		EPSName += ".eps";
		if (!name.Contains("Accumdist")) {
			c.Print(EPSName);
		}
	}
	std::cout << "Plots::Plot1Dsigtobkg done" << std::endl;
}

void Plots::Plot1D_DataMC_Compare(std::vector<TH1D> histo, std::vector<int> colour, std::vector<TString> legend) {
	std::cout << "Plots::Plot1D_DataMC_Compare Start" << std::endl;
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadBottomMargin(0.20);
	gStyle->SetPadLeftMargin(0.175);
	gStyle->SetPadRightMargin(0.05);
	TCanvas c("c", "c", 200, 10, 750, 750);
	c.Clear();

	if (histo.size() > 4) {
		for (int i = 0; i < histo.size(); i++) {
			if (i > 1) {
				histo.at(1).Add(&histo.at(i), 1);
			}
		}
		if (doscale) {
			if (verbose)
				std::cout << "Comparison Integral Scale:" << histo.at(0).Integral() << " " << histo.at(1).Integral() << std::endl;
			histo.at(1).Scale(histo.at(0).Integral() / histo.at(1).Integral());
		}
		for (int b = 1; b <= histo.at(0).GetNbinsX(); b++) {
			double data = histo.at(0).GetBinContent(b);
			double dataerror = histo.at(0).GetBinError(b);
			double mc = histo.at(1).GetBinContent(b);
			double mcerror = histo.at(1).GetBinError(b);
			if (mc != 0) {
				histo.at(0).SetBinContent(b, data / mc);
				histo.at(0).SetBinError(b, dataerror / mc);
				histo.at(2).SetBinContent(b, 1);
				histo.at(2).SetBinError(b, mcerror / mc);

			} else {
				histo.at(0).SetBinContent(b, 0);
				histo.at(0).SetBinError(b, 0);
				histo.at(2).SetBinContent(b, 0);
				histo.at(2).SetBinError(b, 2);
			}
			histo.at(1).SetBinContent(b, data - mc);
			histo.at(1).SetBinError(b, dataerror);
			histo.at(3).SetBinContent(b, 0);
			histo.at(3).SetBinError(b, mcerror);

		}

		histo.at(0).SetYTitle("Data/MC");
		histo.at(0).SetMinimum(0.0);
		histo.at(0).SetMaximum(2.0);
		histo.at(0).Draw("E");
		histo.at(0).Draw("axissame");
		histo.at(2).SetFillColor(5);
		histo.at(2).SetLineColor(5);
		histo.at(2).SetMarkerColor(5);
		histo.at(2).SetMarkerSize(0.0001);
		histo.at(2).Draw("E2same");
		histo.at(0).Draw("Esame");
		TString name = histo.at(0).GetName();
		name += "_sig";
		name += "_index_";
		name += 0;
		TString EPSName = "EPS/";
		EPSName += name;
		EPSName += ".eps";
		if (!name.Contains("Accumdist")) {
			c.Print(EPSName);
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		histo.at(1).SetYTitle("Data-MC");
		histo.at(1).Draw("E");
		histo.at(1).Draw("axissame");
		histo.at(3).SetFillColor(5);
		histo.at(3).SetLineColor(5);
		histo.at(3).SetMarkerColor(5);
		histo.at(3).SetMarkerSize(0.0001);
		histo.at(3).Draw("E2same");
		histo.at(1).Draw("Esame");

		name = histo.at(0).GetName();
		name += "_sigtobkg";
		name += "_index_";
		name += 0;
		EPSName = "EPS/";
		EPSName += name;
		EPSName += ".eps";
		if (!name.Contains("Accumdist")) {
			c.Print(EPSName);
		}

	}
	std::cout << "Plots::Plot1D_DataMC_Compare Finished" << std::endl;
}

void Plots::Plot2D(std::vector<TH2D> histo, std::vector<int> colour, std::vector<TString> legend) {
	std::cout << "Plots::Plot2D" << std::endl;
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadBottomMargin(0.22);
	gStyle->SetPadLeftMargin(0.180);
	gStyle->SetPadRightMargin(0.175);
//	gStyle->SetPalette(53,0);
	if (verbose)
		std::cout << "Starting Plot2D" << histo.size() << " " << colour.size() << " " << legend.size() << std::endl;
	if (histo.size() > 0) {
		std::vector<int> types;
		TString LegType = "";
		double integral = 0;
		for (int i = 0; i < histo.size(); i++) {
			if (verbose) std::cout << "Legend " << legend.at(i) << " Integal " << histo.at(i).Integral() << std::endl;
			if (LegType != legend.at(i)) {
				LegType = legend.at(i);
				integral = 0;
				for (int j = i; j < histo.size(); j++) {
					if (LegType == legend.at(j))
						integral += histo.at(j).Integral();
				}
				if (integral != 0) {
					types.push_back((i!=0) ? i+1: i); // +1 as later on the TotalMC histo is added
					if (verbose) std::cout << "adding one type to types list" << std::endl;
				}
				integral = 0;
			}
		}

		LegType = "";
		std::vector<TH2D>::iterator it;
		it = histo.begin();
		TString tempname = histo.at(0).GetName();
		tempname += "TotalMC";
		histo.insert(it + 1, *((TH2D*) histo.at(0).Clone(tempname.Data())));
		types.insert(types.begin()+1, 1);
		std::vector<TString>::iterator itleg = legend.begin();
		legend.insert(itleg + 1, "Total MC");
		std::vector<int>::iterator itcol = colour.begin();
		colour.insert(itcol + 1, 55);
		histo.at(1).Reset();
		for (int i = 2; i < histo.size(); i++) {
			histo.at(1).Add(&histo.at(i));
		}
		if (verbose)
			std::cout << "n(types) = " << types.size() << std::endl;
		LegType = "null";

		for (int i = 0; i < histo.size(); i++) {
			if (verbose)
				std::cout << "histo i" << i << std::endl;
			if (LegType != legend.at(i)) {
				LegType = legend.at(i);
				for (int j = i + 1; j < histo.size(); j++) {
					if (verbose)
						std::cout << "histo j" << j << std::endl;
					if (legend.at(j) == legend.at(i)){
						histo.at(i).Add(&histo.at(j));
						if (verbose) std::cout << "Added histo " << j << " to histo " << i << std::endl;
					}
				}
				if (verbose)
					std::cout << "legend " << i << " " << histo.at(i).Integral() << std::endl;
				if (histo.at(i).Integral() != 0) {
					if (verbose)
						std::cout << "legend " << i << histo.at(i).Integral() << std::endl;
				}
			}
		}

		// loop to allow splitting of TH2D histograms on multiple canvases
		int nPlotted(0);
		while (nPlotted < types.size()) {
			if (verbose) std::cout << "nPlotted = " << nPlotted << ", ntypes = " << types.size() << std::endl;
			// split canvas for TH2D histograms
			TCanvas c("c", "c", 200, 10, 750, 750);
			c.Clear();

			if (types.size() - nPlotted == 1) {
				if (verbose) std::cout << "Canvas: Divide(1,1)"<< std::endl;
				c.Divide(1, 1);
			} else if (types.size() - nPlotted == 2 || types.size() - nPlotted == 3 || types.size() - nPlotted == 4) {
				if (verbose) std::cout << "Canvas: Divide(2,2)"<< std::endl;
				c.Divide(2, 2);
			} else {
				if (verbose) std::cout << "Canvas: Divide(2,3)"<< std::endl;
				// plot only up to 6 TH2D in one canvas
				c.Divide(2, 3);
			}


			for (int i = nPlotted; (i < nPlotted+6) && (i < types.size());i++) {
				int type = types.at(i);
				if (verbose) std::cout << "loop: " << i << ", type = " << types.at(i) << " legend.at(type) = " << legend.at(type) << std::endl;
				c.cd(i - nPlotted + 1);
				if (verbose) std::cout << "Draw histogram " << type << " into subcanvas " << i - nPlotted +1 << " with Integral = " << histo.at(type).Integral() << std::endl;
				histo.at(type).Draw("colz");
				c.Update();
				TLatex l;
				l.SetNDC();
				l.DrawLatex(0.25, 0.90, legend.at(type));
				}

			TString name = histo.at(0).GetName();
			TString EPSName = "EPS/";
			c.Update();
			name = histo.at(0).GetName();
			name += "_index_";
			name += 0;
			EPSName = "EPS/";
			EPSName += name;
			EPSName += nPlotted;
			EPSName += ".eps";
			c.Print(EPSName);

			nPlotted += 6;
		}
	}
}

void Plots::Plot3D(std::vector<TH3F> histo, std::vector<int> colour, std::vector<TString> legend) {
	std::cout << "Plots::Plot3D" << std::endl;
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadBottomMargin(0.22);
	gStyle->SetPadLeftMargin(0.180);
	gStyle->SetPadRightMargin(0.175);
	gStyle->SetTitleXOffset(2.4);
	gStyle->SetTitleYOffset(2.4);

	TCanvas c("c", "c", 200, 10, 750, 750);
	if (histo.size() > 0) {
		if (verbose)
			std::cout << "Starting Plot3D" << histo.size() << " " << colour.size() << " " << legend.size() << std::endl;
		c.Clear();
		c.Divide(2, 2);
		int type = 0;
		int index = 0;
		std::cout << "A" << std::endl;
		TLegend leg(0.1, 0.4, 0.9, 0.9);
		leg.SetBorderSize(0);
		leg.SetFillStyle(4000);
		leg.SetFillColor(0);
		leg.SetTextSize(0.03);
		leg.SetMargin(0.15);
		leg.SetNColumns(1);
		leg.SetColumnSeparation(0.05);
		for (unsigned int p = 0; p < 4; p++) {
			for (int i = 0; i < histo.size(); i++) {
				histo.at(i).SetMarkerSize(0.25);
				histo.at(i).SetMarkerColor(colour.at(i));
				histo.at(i).GetXaxis()->SetTitleOffset(2.0);
				histo.at(i).GetYaxis()->SetTitleOffset(2.0);
				histo.at(i).GetZaxis()->SetTitleOffset(2.0);
				if (p == 0)
					gPad->SetPhi(30);
				if (p == 1)
					gPad->SetPhi(150);
				if (p == 2)
					gPad->SetPhi(240);
				if (colour.at(i) == 10)
					histo.at(i).SetMarkerColor(16);
				if (colour.at(i) == 414)
					histo.at(i).SetMarkerColor(3);
				if (p != 3) {
					c.cd(p + 1);

					if (i == 0) {
						histo.at(i).Draw();
					} else {
						histo.at(i).Draw("same");
					}
				} else {
					leg.AddEntry(&histo.at(i), legend.at(i), "pe");
				}
			}
		}
		c.cd(4);
		leg.Draw();
		TString name = histo.at(0).GetName();
		TString EPSName = "EPS/";
		c.Update();
		name = histo.at(0).GetName();
		name += "_index_";
		name += 0;
		EPSName = "EPS/";
		EPSName += name;
		EPSName += ".eps";
		c.Print(EPSName);
		if (verbose)
			std::cout << "Histo 3D " << EPSName << std::endl;
	}
}

void Plots::CMSStyle1() {
	std::cout << "Configuring Plots::CMSStyle1()" << std::endl;
	// use plain black on white colors
	Int_t icol = 0; // WHITE
	gStyle->SetFrameBorderMode(icol);
	gStyle->SetFrameFillColor(icol);
	gStyle->SetCanvasBorderMode(icol);
	gStyle->SetCanvasColor(icol);
	gStyle->SetPadBorderMode(icol);
	gStyle->SetPadColor(icol);
	gStyle->SetStatColor(icol);
	//gStyle->SetFillColor(icol); // don't use: white fill color for *all* objects
	// set the paper & margin sizes
	gStyle->SetPaperSize(20, 26);

	// set margin sizes
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadBottomMargin(0.16);
	gStyle->SetPadLeftMargin(0.16);

	// set title offsets (for axis label)
	gStyle->SetTitleXOffset(1.4);
	gStyle->SetTitleYOffset(1.4);

	// use large fonts
	//Int_t font=72; // Helvetica italics
	Int_t font = 42; // Helvetica
	Double_t tsize = 0.05;
	gStyle->SetTextFont(font);

	gStyle->SetTextSize(tsize);
	gStyle->SetLabelFont(font, "x");
	gStyle->SetTitleFont(font, "x");
	gStyle->SetLabelFont(font, "y");
	gStyle->SetTitleFont(font, "y");
	gStyle->SetLabelFont(font, "z");
	gStyle->SetTitleFont(font, "z");

	gStyle->SetLabelSize(tsize, "x");
	gStyle->SetTitleSize(tsize, "x");
	gStyle->SetLabelSize(tsize, "y");
	gStyle->SetTitleSize(tsize, "y");
	gStyle->SetLabelSize(tsize, "z");
	gStyle->SetTitleSize(tsize, "z");

	// use bold lines and markers
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerSize(1.2);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLineStyleString(2, "[12 12]"); // postscript dashes

	// get rid of X error bars
	//gStyle->SetErrorX(0.001);
	// get rid of error bar caps
	gStyle->SetEndErrorSize(0.);

	// do not display any of the standard histogram decorations
	gStyle->SetOptTitle(0);
	//gStyle->SetOptStat(1111);
	gStyle->SetOptStat(0);
	//gStyle->SetOptFit(1111);
	gStyle->SetOptFit(0);
	// put tick marks on top and RHS of plots
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPalette(1);
	gROOT->ForceStyle(true);

}

void Plots::CMSStyle2() {
	std::cout << "Configuring Plots::CMSStyle2()" << std::endl;
	// use plain black on white colors
	Int_t icol = 0; // WHITE
	gStyle->SetFrameBorderMode(icol);
	gStyle->SetFrameFillColor(icol);
	gStyle->SetCanvasBorderMode(icol);
	gStyle->SetCanvasColor(icol);
	gStyle->SetPadBorderMode(icol);
	gStyle->SetPadColor(icol);
	gStyle->SetStatColor(icol);
	//gStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

	// set the paper & margin sizes
	gStyle->SetPaperSize(20, 26);

	// set margin sizes
	gStyle->SetPadTopMargin(0.05);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadBottomMargin(0.16);
	gStyle->SetPadLeftMargin(0.16);

	// set title offsets (for axis label)
	gStyle->SetTitleXOffset(1.4);
	gStyle->SetTitleYOffset(1.4);

	// use large fonts
	//Int_t font=72; // Helvetica italics
	Int_t font = 42; // Helvetica
	Double_t tsize = 0.05;
	gStyle->SetTextFont(font);

	gStyle->SetTextSize(tsize);
	gStyle->SetLabelFont(font, "x");
	gStyle->SetTitleFont(font, "x");
	gStyle->SetLabelFont(font, "y");
	gStyle->SetTitleFont(font, "y");
	gStyle->SetLabelFont(font, "z");
	gStyle->SetTitleFont(font, "z");

	gStyle->SetLabelSize(tsize, "x");
	gStyle->SetTitleSize(tsize, "x");
	gStyle->SetLabelSize(tsize, "y");
	gStyle->SetTitleSize(tsize, "y");
	gStyle->SetLabelSize(tsize, "z");
	gStyle->SetTitleSize(tsize, "z");

	// use bold lines and markers
	gStyle->SetMarkerStyle(20);
	gStyle->SetMarkerSize(1.2);
	gStyle->SetHistLineWidth(2);
	gStyle->SetLineStyleString(2, "[12 12]"); // postscript dashes

	// get rid of X error bars
	//gStyle->SetErrorX(0.001);
	// get rid of error bar caps
	gStyle->SetEndErrorSize(0.);

	// do not display any of the standard histogram decorations
	gStyle->SetOptTitle(0);
	//gStyle->SetOptStat(1111);
	gStyle->SetOptStat(0);
	//gStyle->SetOptFit(1111);
	gStyle->SetOptFit(0);

	// put tick marks on top and RHS of plots
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetPalette(1);
	gROOT->ForceStyle(true);

}

void Plots::CMSLabel(Double_t x, Double_t y, Color_t color) {
	TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize);
	l.SetNDC();
	l.SetTextFont(72);
	l.SetTextColor(color);

	double delx = 0.115 * 696 * gPad->GetWh() / (472 * gPad->GetWw());

	l.DrawLatex(x, y, "CMS");
	TLatex p;
	p.SetNDC();
	p.SetTextFont(42);
	p.SetTextColor(color);
	if (plotLabel == cmsPrivate) {
		p.DrawLatex(x + delx, y, "Private");
	}
	if (plotLabel == cmsPreliminary) {
		p.DrawLatex(x + delx, y, "Preliminary");
	}
}

