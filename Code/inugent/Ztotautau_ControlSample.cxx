#include "Ztotautau_ControlSample.h"
#include "TLorentzVector.h"
#include <cstdlib>
#include "HistoConfig.h"
#include <iostream>

#include "Tools.h"
#include "PDG_Var.h"
#include "TauDataFormat/TauNtuple/interface/DataMCType.h"
#include "TauSolver.h"
#include "SimpleFits/FitSoftware/interface/PDGInfo.h"

Ztotautau_ControlSample::Ztotautau_ControlSample(TString Name_, TString id_):
  Selection(Name_,id_)
  ,channel(muontag)
  ,muoneta(2.1)
  ,jeteta(2.3)
  ,TauTrackPtThreshold(5.0)
{
  if(Get_Name().Contains("muontag")) channel=muontag;
  if(Get_Name().Contains("electrontag")) channel=muontag;  // not implemented yet
  if(Get_Name().Contains("rhotag")) channel=muontag;       // not implemented yet
  if(Get_Name().Contains("threepiontag")) channel=muontag; // not implemented yet
  //verbose=true;
}

Ztotautau_ControlSample::~Ztotautau_ControlSample(){
  for(int j=0; j<Npassed.size(); j++){
    std::cout << "Ztotautau_ControlSample::~Ztotautau_ControlSample Selection Summary before: " 
	 << Npassed.at(j).GetBinContent(1)     << " +/- " << Npassed.at(j).GetBinError(1)     << " after: "
	 << Npassed.at(j).GetBinContent(NCuts) << " +/- " << Npassed.at(j).GetBinError(NCuts) << std::endl;
  }
  std::cout << "Ztotautau_ControlSample::~Ztotautau_ControlSample()" << std::endl;
}

void  Ztotautau_ControlSample::Configure(){
  // Setup Cut Values
  for(int i=0; i<NCuts;i++){
    cut.push_back(0);
    value.push_back(0);
    pass.push_back(false);
    if(i==TriggerOk)          cut.at(TriggerOk)=1;
    if(i==hasTag)             cut.at(hasTag)=1;
    if(i==TagPtmin)           cut.at(TagPtmin)=25;
    if(i==TagPtmax)           cut.at(TagPtmax)=38;
    if(i==TagIso)             cut.at(TagIso)=0.2;
    if(i==NJets)              cut.at(NJets)=1;
    if(i==JetPt)              cut.at(JetPt)=10;
    if(i==deltaPhi)           cut.at(deltaPhi)=TMath::Pi()*7.0/8.0;
    if(i==MET)                cut.at(MET)=40;
    if(i==MT)                 cut.at(MT)=30;
    if(i==etaq)               cut.at(etaq)=0.1;
    if(i==JetTrackPtMax)      cut.at(JetTrackPtMax)=TauTrackPtThreshold;
    if(i==ZMassmin)           cut.at(ZMassmin)=20;
    if(i==ZMassmax)           cut.at(ZMassmax)=80;
    if(i==HT)                 cut.at(HT)=200;
    if(i==charge)             cut.at(charge)=0.0;
    if(i==MaxTracksinJet)     cut.at(MaxTracksinJet)=5.0;
    if(i==MinTracksinJet)     cut.at(MinTracksinJet)=0.0;
  }

  TString hlabel;
  TString htitle;
  for(unsigned int i=0; i<NCuts; i++){
    title.push_back("");
    distindx.push_back(false);
    dist.push_back(std::vector<float>());
    Accumdist.push_back(std::vector<TH1D>());
    Nminus1dist.push_back(std::vector<TH1D>());


    TString c="_Cut_";
    if(i<10)c+="0";
    c+=i;
  
    if(i==TriggerOk){
      title.at(i)="Trigger ";
      hlabel="Trigger ";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TriggerOk_",htitle,2,-0.5,1.5,hlabel,"Events"));
    }
    //-----------
   else if(i==hasTag){
      title.at(i)="has Tag";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="has Tag (bool)";
      if(channel==muontag){
	title.at(i)="Number of Muons $=$";
	title.at(i)+=cut.at(hasTag);
	htitle=title.at(i);
	htitle.ReplaceAll("$","");
	htitle.ReplaceAll("\\","#");
	hlabel="N_{Muons}";
      }
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_hasTag_",htitle,6,-0.5,5.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_hasTag_",htitle,6,-0.5,5.5,hlabel,"Events"));
    }
   else if(i==TagPtmin){
      title.at(i)="$P_{T}^{Tag}>$";
      title.at(i)+=cut.at(TagPtmin);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="P_{T}^{Tag} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagPtmin_",htitle,50,0,100,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagPtmin_",htitle,50,0,100,hlabel,"Events"));
    }
   else if(i==TagPtmax){
     title.at(i)="$P_{T}^{Tag}<$";
     title.at(i)+=cut.at(TagPtmax);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="P_{T}^{Tag} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagPtmax_",htitle,50,0,100,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagPtmax_",htitle,50,0,100,hlabel,"Events"));
   }
   else if(i==TagIso){
      title.at(i)="Relative Isolation (Tag) $<$";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(TagIso));
      title.at(i)+=buffer;
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="Relative Isolation (Tag)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_TagIso_",htitle,40,0,1,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_TagIso_",htitle,40,0,1,hlabel,"Events"));
    }
   else if(i==NJets){
     title.at(i)="$Number of Jet>=$";
     title.at(i)+=cut.at(NJets);
     title.at(i)+="";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="N_{Jets}";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_NJets_",htitle,21,-0.0,20.5,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_NJets_",htitle,21,-0.5,20.5,hlabel,"Events"));
   }
   else if(i==JetPt){
     title.at(i)="$P_{T}^{Jet}>$";
     title.at(i)+=cut.at(JetPt);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="P_{T}^{Jet} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_JetPt_",htitle,20,0,200,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_JetPt_",htitle,20,0,200,hlabel,"Events"));
   }
   else if(i==MET){
      title.at(i)="$E_{T}^{Miss} < $";
      title.at(i)+=cut.at(MET);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="E_{T}^{Miss} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MET_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MET_",htitle,40,0,200,hlabel,"Events"));
    } 
   else if(i==MT){
     title.at(i)="$M_{T} < $";
     title.at(i)+=cut.at(MT);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="M_{T} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MT_",htitle,40,0,200,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MT_",htitle,40,0,200,hlabel,"Events"));
   }
   else if(i==HT){
     title.at(i)="$H_{T} < $";
     title.at(i)+=cut.at(HT);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="H_{T} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_HT_",htitle,50,0,500,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_HT_",htitle,50,0,500,hlabel,"Events"));
   }
   else if(i==etaq){
     title.at(i)="$q_{\\mu}|\\eta|$";
     char buffer[50];
     sprintf(buffer,"%5.2f",cut.at(etaq));
     title.at(i)+=buffer;
     title.at(i)+="";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="q_{#mu}|#eta|";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_etaq_",htitle,40,-2.0,2.0,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_etaq_",htitle,40,-2.0,2.0,hlabel,"Events"));
   }
   else if(i==JetTrackPtMax){
     title.at(i)="$P_{T}^{max\\ Track\\ from\\ Jet} $";
     char buffer[50];
     sprintf(buffer,"%5.2f",cut.at(JetTrackPtMax));
     title.at(i)+=buffer;
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="P_{T}^{Max Track from Jet}";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_JetTrackPtMax_",htitle,20,0,50,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_JetTrackPtMax_",htitle,20,0,50,hlabel,"Events"));
   }
   else if(i==deltaPhi){
      title.at(i)="$\\Delta\\phi(Tag,Jet) < $";
      char buffer[50];
      sprintf(buffer,"%5.2f",cut.at(deltaPhi));
      title.at(i)+=buffer;
      title.at(i)+="(rad)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="#Delta#phi (rad)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_deltaPhi_",htitle,32,0,TMath::Pi(),hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_deltaPhi_",htitle,32,0,TMath::Pi(),hlabel,"Events"));
    } 
   else if(i==ZMassmin){
      title.at(i)="$M_{Z}< $";
      title.at(i)+=cut.at(ZMassmin);
      title.at(i)+="(GeV)";
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="M_{Z} (GeV)";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassmin_",htitle,40,0,200,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassmin_",htitle,40,0,200,hlabel,"Events"));
    } 
   else if(i==ZMassmax){
     title.at(i)="$M_{Z}< $";
     title.at(i)+=cut.at(ZMassmax);
     title.at(i)+="(GeV)";
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="M_{Z} (GeV)";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_ZMassmax_",htitle,40,0,200,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_ZMassmax_",htitle,40,0,200,hlabel,"Events"));
   }
   else if(i==charge){
      title.at(i)="$C_{\\mu}\\times\\sum_{i=0}^{3}C_{i}|_{P_{T}^{max}}$";
      title.at(i)+=cut.at(charge);
      htitle=title.at(i);
      htitle.ReplaceAll("$","");
      htitle.ReplaceAll("\\","#");
      hlabel="C_{#mu}#times #sum_{i=0}^{3}C_{i}|_{P_{T}^{max}}";
      Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_charge_",htitle,21,-10.5,10.5,hlabel,"Events"));
      Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_charge_",htitle,21,-10.5,10.5,hlabel,"Events"));
    } 
   else if(i==MaxTracksinJet){
     title.at(i)="$N_{Tracks}^{Jet}<$";
     title.at(i)+=cut.at(MaxTracksinJet);
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="N_{Tracks}^{Jet}";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MaxTracksinJet_",htitle,21,-0.5,20.5,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MaxTracksinJet_",htitle,21,-0.5,20.5,hlabel,"Events"));
     distindx.at(i)=true;
     Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_MaxTracksinJet_","P_{T,Track} (N-1 Distribution)",100,0,200,"P_{T,Track} (GeV)","Events");
     Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accum1dist_MaxTracksinJet_","P_{T,Track} (Accumulative Distribution)",100,0,200,"P_{T,Track} (GeV)","Events");
   }
   else if(i==MinTracksinJet){
     title.at(i)="$N_{Tracks}^{Jet}>$";
     title.at(i)+=cut.at(MinTracksinJet);
     htitle=title.at(i);
     htitle.ReplaceAll("$","");
     htitle.ReplaceAll("\\","#");
     hlabel="N_{Tracks}^{Jet}";
     Nminus1.push_back(HConfig.GetTH1D(Name+c+"_Nminus1_MinTracksinJet_",htitle,21,-0.5,20.5,hlabel,"Events"));
     Nminus0.push_back(HConfig.GetTH1D(Name+c+"_Nminus0_MinTracksinJet_",htitle,21,-0.5,20.5,hlabel,"Events"));
     distindx.at(i)=true;
     Nminus1dist.at(i)=HConfig.GetTH1D(Name+c+"_Nminus1dist_MinTracksinJet_","P_{T,Track} (N-1 Distribution)",100,0,200,"P_{T,Track} (GeV)","Events");
     Accumdist.at(i)=HConfig.GetTH1D(Name+c+"_Accum1dist_MinTracksinJet_","P_{T,Track} (Accumulative Distribution)",100,0,200,"P_{T,Track} (GeV)","Events");
   }

    //-----------
  }
  // Setup NPassed Histogams
  Npassed=HConfig.GetTH1D(Name+"_NPass","Cut Flow",NCuts+1,-1,NCuts,"Number of Accumulative Cuts Passed","Events");

  MinTauTrackPt=HConfig.GetTH1D(Name+"_MinTauTrackPt_","P_{T,Track}",100,0,25,"P_{T,Track}^{Min} (GeV)","Events");
  MaxTauTrackPt=HConfig.GetTH1D(Name+"_MaxTauTrackPt_","P_{T,Track}",100,0.0,50,"P_{T,Track}^{Max} (GeV)","Events");
  TauTrackdphitheta=HConfig.GetTH1D(Name+"_TauTrackdphitheta_","Track anagle",100,0,4,"#delta(#phi,#theta) (rad)","Events");


  // Setup Extra Histograms
  NVtx=HConfig.GetTH1D(Name+"_NVtx","NVtx",26,-0.5,25.5,"Number of Vertex","Events");
  NGoodVtx=HConfig.GetTH1D(Name+"_NGoodVtx","NGoodVtx",26,-0.05,25.5,"Number of Good Vertex","Events");
  NTrackperVtx=HConfig.GetTH1D(Name+"_NTracksperVtx","NTracksperVtx",151,-0.5,150.5,"Number of Track per Vertex","Events");
  TagEtaPT=HConfig.GetTH2D(Name+"_TagEtaPT","TagEtaPT",25,0,2.5,50,0,50,"|#eta|","P_{T}^{Tag}");

  TauCandFound =HConfig.GetTH1D(Name+"_TauCandFound","TauCandFound",2,-0.5,1.5,"Found #tau Candidate (bool)","Events");
  TauCandEtaPhi =HConfig.GetTH2D(Name+"_TauCandPtEta","TauCandPtEta",25,0,2.5,40,0,200,"|#eta|","#phi (rad)");

  TauCandPhi =HConfig.GetTH1D(Name+"_TauCandPhi","TauCandPhi",32,-TMath::Pi(),TMath::Pi(),"#phi_{KF-#tau} (rad)","Events");
  TauCandPhiRes =HConfig.GetTH1D(Name+"_TauCandPhiRes","TauCandPhiRes",20,-0.02,0.02,"#sigma(#phi_{KF-#tau}) (rad)","Events");
  TauCandEta =HConfig.GetTH1D(Name+"_TauCandEta","TauCandEta",25,-2.5,2.5,"#eta_{KF-#tau}","Events");
  TauCandEtaRes =HConfig.GetTH1D(Name+"_TauCandEtaRes","TauCandEtaRes",50,-0.05,0.05,"#sigma(#eta_{KF-#tau}) ","Events");
  TauCandE =HConfig.GetTH1D(Name+"_TauCandE","TauCandE",40,0,200,"E_{KF-#tau} (GeV)","Events");
  TauCandERes =HConfig.GetTH1D(Name+"_TauCandERes","TauCandERes",20,-50,50,"#sigma(E_{KF-#tau}) (GeV)","Events");
  TauCandP =HConfig.GetTH1D(Name+"_TauCandP","TauCandP",40,0,200,"E_{KF-#tau} (GeV)","Events");
  TauCandPRes =HConfig.GetTH1D(Name+"_TauCandPRes","TauCandPRes",20,-50,50,"#sigma(E_{KF-#tau}) (GeV)","Events");
  TauCandPT =HConfig.GetTH1D(Name+"_TauCandPT","TauCandPT",40,0,200,"E_{KF-#tau} (GeV)","Events");
  TauCandPTRes =HConfig.GetTH1D(Name+"_TauCandPTRes","TauCandPTRes",20,-50,50,"#sigma(E_{KF-#tau}) (GeV)","Events");


  TauCandMass =HConfig.GetTH1D(Name+"_TauCandMass","TauCandMass",100,0,2,"M(KF-#tau) (GeV)","Events");
  TauCandPiMass =HConfig.GetTH1D(Name+"_TauCandPiMass","TauCandPiMass",100,0,1,"M(#pi) (GeV)","Events");
  TauCandNuMass =HConfig.GetTH1D(Name+"_TauCandNuMass","TauCandNuMass",110,-0.1,1,"M(#nu) (GeV)","Events");

  TauSolutionResult = HConfig.GetTH1D(Name+"_TauSolutionResult","TauSolutionResult",3,-1.5,1.5,"Solution: -:R<0:+","Events");
  EstimatedTauE =HConfig.GetTH1D(Name+"_EstimatedTauCandE","EstimatedTauCandE",40,0,200,"E_{Est-#tau} (GeV)","Events");
  EstimatedTauPhi =HConfig.GetTH1D(Name+"_EstimatedTauCandPhi","EstimatedTauCandPhi",32,-TMath::Pi(),TMath::Pi(),"#phi_{Est-#tau} (rad)","Events");
  EstimatedTauEta  =HConfig.GetTH1D(Name+"_EstimatedTauCandEta","EstimatedTauCandEta",25,-2.5,2.5,"#eta_{Est-#tau}","Events");

  EstimatedTauERes =HConfig.GetTH1D(Name+"_TauEstimatedERes","TauEstimatedERes",40,-50,50,"#sigma(E_{Est-#tau}) (GeV)","Events");
  EstimatedTauPhiRes =HConfig.GetTH1D(Name+"_TauEstimatedPhiRes","TauEstimatedPhiRes",100,-0.2,0.2,"#sigma(#phi_{Est-#tau}) (rad)","Events");
  EstimatedTauEtaRes =HConfig.GetTH1D(Name+"_TauEstimatedEtaRes","TauEstimatedEtaRes",100,-0.2,0.2,"#sigma(#eta_{Est-#tau}) ","Events");

  EstimatedTauDirPhiRes=HConfig.GetTH1D(Name+"_EstimatedTauCandDirPhiRes","EstimatedTauCandDirPhiRes",100,-0.2,0.2,"#sigma(#phi_{Direction-#tau}) (rad)","Events");
  EstimatedTauDirEtaRes=HConfig.GetTH1D(Name+"_EstimatedTauCandDirEtaRes","EstimatedTauCandDirEtaRes",100,-0.2,0.2,"#sigma(#eta_{Direction-#tau}) (rad)","Events");

  KFTau_Fit_chiprob=HConfig.GetTH1D(Name+"_KFTau_Fit_prob","KFTau_Fit_prob",25,0,1,"Kinematic Fit Probability","Events");
  KFTau_Fit_a1mass=HConfig.GetTH1D(Name+"_KFTau_Fit_a1mass","KFTau_Fit_a1mass",25,0,2.5,"Kinematic Fit a_{1} Mass","Events");
  KFTau_Fit_chi2=HConfig.GetTH1D(Name+"_KFTau_Fit_chi2","KFTau_Fit_chi",25,0,25,"Kinematic Fit #chi^{2}","Events");
  KFTau_Fit_ndf=HConfig.GetTH1D(Name+"_KFTau_Fit_ndf","KFTau_Fit_ndf",25,0,25,"Kinematic Fit ndf","Events");
  KFTau_Fit_ambiguity=HConfig.GetTH1D(Name+"_KFTau_Fit_ambiguity","KFTau_Fit_ambiguity",5,-2.5,2.5,"Kinematic Fit Ambiguity","Events");
  KFTau_Fit_csum=HConfig.GetTH1D(Name+"_KFTau_Fit_csum","KFTau_Fit_csum",100,0,1.0,"Kinematic Fit csum","Events");
  KFTau_Fit_iterations=HConfig.GetTH1D(Name+"_KFTau_Fit_iterations","KFTau_Fit_iterations",25,0,25,"Kinematic Fit Iterations","Events");
  KFTau_Fit_TauEnergyFraction=HConfig.GetTH1D(Name+"_KFTau_Fit_TauEnergyFraction","KFTau_Fit_TauEnergyFraction",50,0,5.0,"Transverse Energy Fration,","Events");
  KFTau_Fit_PV_PV_significance=HConfig.GetTH1D(Name+"_ KFTau_Fit_PV_PV_significance"," KFTau_Fit_PV_PV_significance",50,0,10,"Kinematic Fit PV-PV_{Rotated} Significance","Events");
  KFTau_Fit_SV_PV_significance=HConfig.GetTH1D(Name+"_ KFTau_Fit_SV_PV_significance","KFTau_Fit_SV_PV_significance",50,0,10,"Kinematic Fit SV-PV_{Rotated} Significance","Events");




  PFTauCandPhi =HConfig.GetTH1D(Name+"_PFTauCandPhi","PFTauCandPhi",32,-TMath::Pi(),TMath::Pi(),"#phi_{KF-#tau} (rad)","Events");
  PFMTauCandPhi =HConfig.GetTH1D(Name+"_PFMTauCandPhi","PFMTauCandPhi",32,-TMath::Pi(),TMath::Pi(),"#phi_{KF-#tau} (rad)","Events");
  KFTauCandPhi =HConfig.GetTH1D(Name+"_KFTauCandPhi","KFTauCandPhi",32,-TMath::Pi(),TMath::Pi(),"#phi_{KF-#tau} (rad)","Events");
  KFSTauCandPhi =HConfig.GetTH1D(Name+"_KFSTauCandPhi","KFSTauCandPhi",32,-TMath::Pi(),TMath::Pi(),"#phi_{KF-#tau} (rad)","Events");
  KFTauVCandPhi =HConfig.GetTH1D(Name+"_KFTauVCandPhi","KFTauVCandPhi",32,-TMath::Pi(),TMath::Pi(),"#phi_{KF-#tau} (rad)","Events");
  KFSTauVCandPhi =HConfig.GetTH1D(Name+"_KFSTauVCandPhi","KFSTauVCandPhi",32,-TMath::Pi(),TMath::Pi(),"#phi_{KF-#tau} (rad)","Events");

  Selection::ConfigureHistograms();
  HConfig.GetHistoInfo(types,CrossSectionandAcceptance,legend,colour);
  for(int i=0;i<CrossSectionandAcceptance.size();i++){
    std::cout << i << " CrossSectionandAcceptance " << CrossSectionandAcceptance.at(i) << std::endl;
  }
}




void  Ztotautau_ControlSample::Store_ExtraDist(){
 Extradist1d.push_back(&NVtx);
 Extradist1d.push_back(&NGoodVtx);
 Extradist1d.push_back(&NTrackperVtx);
 Extradist2d.push_back(&TagEtaPT);

 Extradist1d.push_back(&MaxTauTrackPt);
 Extradist1d.push_back(&MinTauTrackPt);
 Extradist1d.push_back(&TauTrackdphitheta);


 Extradist1d.push_back(&TauCandFound);
 Extradist2d.push_back(&TauCandEtaPhi); 

 Extradist1d.push_back(&TauCandPhi);
 Extradist1d.push_back(&TauCandPhiRes);
 Extradist1d.push_back(&TauCandEta);
 Extradist1d.push_back(&TauCandEtaRes);
 Extradist1d.push_back(&TauCandE);
 Extradist1d.push_back(&TauCandERes);
 Extradist1d.push_back(&TauCandP);
 Extradist1d.push_back(&TauCandPRes);
 Extradist1d.push_back(&TauCandPT);
 Extradist1d.push_back(&TauCandPTRes);

 Extradist1d.push_back(&TauCandMass);
 Extradist1d.push_back(&TauCandPiMass);
 Extradist1d.push_back(&TauCandNuMass);

 Extradist1d.push_back(&TauSolutionResult);
 Extradist1d.push_back(&EstimatedTauE);
 Extradist1d.push_back(&EstimatedTauPhi);
 Extradist1d.push_back(&EstimatedTauEta);

 Extradist1d.push_back(&EstimatedTauERes);
 Extradist1d.push_back(&EstimatedTauPhiRes);
 Extradist1d.push_back(&EstimatedTauEtaRes);

 Extradist1d.push_back(&EstimatedTauDirEtaRes);
 Extradist1d.push_back(&EstimatedTauDirPhiRes);

 Extradist1d.push_back(&KFTau_Fit_chiprob); 
 Extradist1d.push_back(&KFTau_Fit_a1mass);
 Extradist1d.push_back(&KFTau_Fit_chi2);
 Extradist1d.push_back(&KFTau_Fit_ndf);
 Extradist1d.push_back(&KFTau_Fit_ambiguity);
 Extradist1d.push_back(&KFTau_Fit_csum);
 Extradist1d.push_back(&KFTau_Fit_iterations);
 Extradist1d.push_back(&KFTau_Fit_SV_PV_significance);

 Extradist1d.push_back(&PFTauCandPhi);
 Extradist1d.push_back(&PFMTauCandPhi);
 Extradist1d.push_back(&KFTauCandPhi);
 Extradist1d.push_back(&KFSTauCandPhi);
 Extradist1d.push_back(&KFTauVCandPhi);
 Extradist1d.push_back(&KFSTauVCandPhi);


}

void  Ztotautau_ControlSample::doEvent(){
  unsigned int t;
  int id(Ntp->GetMCID());
  if(!HConfig.GetHisto(Ntp->isData(),id,t)){ std::cout << "failed to find id" <<std::endl; return;}

  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() A" << std::endl;

  value.at(TriggerOk)=0;
  if(Ntp->TriggerAccept("IsoMu24"))value.at(TriggerOk)+=1;
  pass.at(TriggerOk)=value.at(TriggerOk)==cut.at(TriggerOk);
  
  // Apply Selection
  unsigned mu_idx(999),nmus(0);
  double mu_pt(0);
  if(channel==muontag){
    for(unsigned int i=0;i<Ntp->NMuons();i++){
      if(Ntp->isGoodMuon(i) && fabs(Ntp->Muon_p4(i).Eta())<muoneta){
	if(mu_pt<Ntp->Muon_p4(i).Pt()){mu_idx=i;mu_pt=Ntp->Muon_p4(i).Pt();}
	nmus++;
      }
    }
    if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() C" << std::endl;
    value.at(hasTag)=nmus;
    pass.at(hasTag)=(value.at(hasTag)==cut.at(hasTag));
   
    value.at(TagPtmin)=mu_pt;
    pass.at(TagPtmin)=(value.at(TagPtmin)>cut.at(TagPtmin));
    value.at(TagPtmax)=mu_pt;
    pass.at(TagPtmax)=(value.at(TagPtmax)<cut.at(TagPtmax));

    if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() C2 - " << mu_idx << " " << Ntp->NMuons() << std::endl;
    if(mu_idx!=999){value.at(TagIso) = (Ntp->Muon_emEt05(mu_idx) + Ntp->Muon_hadEt05(mu_idx) + Ntp->Muon_sumPt05(mu_idx))/Ntp->Muon_p4(mu_idx).Pt();}
    else{value.at(TagIso)=999;}
    pass.at(TagIso)=(value.at(TagIso)<=cut.at(TagIso));
    
  }
  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() d" << std::endl;
  unsigned int jet_idx(999),njets(0);
  std::vector<unsigned int> jet_list;
  double jet_pt(0);
  if(mu_idx!=999){
    for(int i=0;i<Ntp->NPFJets();i++){
      if(Tools::dr(Ntp->Muon_p4(mu_idx),Ntp->PFJet_p4(i))>0.4){
	if(jeteta>fabs(Ntp->PFJet_p4(i).Eta())){
	  double TrackMaxPt=0;
	  std::vector<int> PFJet_Track_idx=Ntp->PFJet_Track_idx(i);
	  for(int j=0; j<PFJet_Track_idx.size();j++){
	    if(PFJet_Track_idx.at(j)<Ntp->NTracks() && PFJet_Track_idx.at(j)>=0){
	      if(Ntp->Track_p4(PFJet_Track_idx.at(j)).Pt()>TrackMaxPt) TrackMaxPt=Ntp->Track_p4(PFJet_Track_idx.at(j)).Pt();
	    }
	  }
	  if(TrackMaxPt>TauTrackPtThreshold){
	    if(jet_pt<Ntp->PFJet_p4(i).Pt()){jet_idx=i;jet_pt=Ntp->PFJet_p4(i).Pt();}
	    njets++;
	  }
	}
      }
    }
  }
  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() e" << std::endl;
  value.at(NJets)=njets;
  pass.at(NJets)=(value.at(NJets)>=cut.at(NJets));

  value.at(JetPt)=jet_pt;
  pass.at(JetPt)=(value.at(JetPt)>=cut.at(JetPt));


  if(jet_idx!=999){
    std::vector<int> PFJet_Track_idx=Ntp->PFJet_Track_idx(jet_idx);
    value.at(JetTrackPtMax)=0;
    for(int i=0; i<PFJet_Track_idx.size();i++){
      if(PFJet_Track_idx.at(i)<Ntp->NTracks() && PFJet_Track_idx.at(i)>=0){
	if(Ntp->Track_p4(PFJet_Track_idx.at(i)).Pt()>value.at(JetTrackPtMax)) value.at(JetTrackPtMax)=Ntp->Track_p4(PFJet_Track_idx.at(i)).Pt();
      }
    }
    pass.at(JetTrackPtMax)=true;//value.at(JetTrackPtMax)>cut.at(JetTrackPtMax);
  }
  else{
    value.at(JetTrackPtMax)=0;
    pass.at(JetTrackPtMax)=false;
  }
  pass.at(JetTrackPtMax)=true;

  if(mu_idx!=999 && jet_idx!=999){
    value.at(deltaPhi)=fabs(Tools::DeltaPhi(Ntp->Muon_p4(mu_idx),Ntp->PFJet_p4(jet_idx)));
  }
  else { 
    value.at(deltaPhi)=0;
  }
  pass.at(deltaPhi)=(fabs(value.at(deltaPhi))>=cut.at(deltaPhi));


  value.at(MET)=Ntp->MET_CorrMVA_et();
  pass.at(MET)=(value.at(MET)<cut.at(MET));

  if(mu_idx!=999){
    value.at(MT)=sqrt(2*(Ntp->MET_CorrMVA_et())*Ntp->Muon_p4(mu_idx).Pt()*fabs(1-cos(Ntp->Muon_p4(mu_idx).Phi()-Ntp->MET_CorrMVA_phi())));
  }
  else{
    value.at(MT)=999;
  }
  pass.at(MT)=(value.at(MT)<=cut.at(MT));

  TLorentzVector Z_lv(0,0,0,0);
  if(jet_idx!=999 && mu_idx!=999){
    std::vector<int> PFJet_Track_idx=Ntp->PFJet_Track_idx(jet_idx);
    TLorentzVector track_p4(1,1,1,1);
    for(int i=0; i<PFJet_Track_idx.size();i++){
      if(PFJet_Track_idx.at(i)<Ntp->NTracks() && PFJet_Track_idx.at(i)>=0){
        if(track_p4.Pt()<Ntp->Track_p4(PFJet_Track_idx.at(i)).Pt())track_p4=Ntp->Track_p4(PFJet_Track_idx.at(i));
      }
    }
    Z_lv+=Ntp->Muon_p4(mu_idx);
    Z_lv+=track_p4;
  }
  value.at(ZMassmin)=Z_lv.M();
  pass.at(ZMassmin)=(value.at(ZMassmin)>cut.at(ZMassmin));
  value.at(ZMassmax)=Z_lv.M();
  pass.at(ZMassmax)=(value.at(ZMassmax)<cut.at(ZMassmax));

  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() g" << std::endl;

  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() h" << std::endl;

  if(mu_idx!=999){
    value.at(etaq)=Ntp->Muon_Charge(mu_idx)*fabs(Ntp->Muon_p4(mu_idx).Eta());
    pass.at(etaq)=value.at(etaq)>cut.at(etaq);
  }
  else{
    value.at(etaq)=0;
    pass.at(etaq)=false;
  }
  pass.at(etaq)=true;

  if(jet_idx!=999 && mu_idx!=999){
    value.at(HT)=0;
    value.at(HT)+=Ntp->PFJet_p4(jet_idx).Pt();
    value.at(HT)+=Ntp->Muon_p4(mu_idx).Pt();
    value.at(HT)+=Ntp->MET_CorrMVA_et();
    for(int i=0;i<Ntp->NPFJets();i++){
      if(i!=jet_idx){
	if(Ntp->isGoodJet(i))value.at(HT)+=Ntp->PFJet_p4(jet_idx).Pt();
      }
    }
  }
  pass.at(HT)=value.at(HT)<cut.at(HT);
  ///////////////////////////////////////
  std::vector<unsigned int> TauTracks;
  if(jet_idx!=999){
    std::vector<int> PFJet_Track_idx=Ntp->PFJet_Track_idx(jet_idx);
    double ptmin(9999), ptmax(0);
    unsigned int ptmaxidx(0);
    for(int i=0; i<PFJet_Track_idx.size();i++){
      if(ptmax<Ntp->Track_p4(PFJet_Track_idx.at(i)).Pt()){ptmax=Ntp->Track_p4(PFJet_Track_idx.at(i)).Pt();ptmaxidx=PFJet_Track_idx.at(i);}
    }
    dist.at(MaxTracksinJet).push_back(ptmax);
    TauTracks.push_back(ptmaxidx);
    for(int i=0; i<PFJet_Track_idx.size();i++){
      if(PFJet_Track_idx.at(i)!=ptmaxidx){
        if(sqrt(pow(Ntp->Track_p4(PFJet_Track_idx.at(i)).DeltaPhi(Ntp->Track_p4(ptmaxidx)),2.0)+
                pow(Ntp->Track_p4(PFJet_Track_idx.at(i)).Theta()-Ntp->Track_p4(ptmaxidx).Theta(),2.0))<0.2){
          TauTracks.push_back(PFJet_Track_idx.at(i));
	  if(ptmin>Ntp->Track_p4(PFJet_Track_idx.at(i)).Pt()) ptmin=Ntp->Track_p4(PFJet_Track_idx.at(i)).Pt();
        }
      }
    }
    dist.at(MinTracksinJet).push_back(ptmin);

    value.at(MaxTracksinJet)=TauTracks.size();
    pass.at(MaxTracksinJet)=true;//cut.at(MaxTracksinJet)>value.at(MaxTracksinJet);                                                                                                                                                          
    value.at(MinTracksinJet)=TauTracks.size();
    pass.at(MinTracksinJet)=true;//cut.at(MinTracksinJet)<value.at(MinTracksinJet);
  }

  value.at(charge)=0;
  if(jet_idx!=999 && mu_idx!=999){
    std::vector<unsigned int> PFJet_Track_idx=TauTracks;//Ntp->PFJet_Track_idx(jet_idx);                                                                                                                                                     
    std::vector<TLorentzVector> track_p4;
    std::vector<float> track_charge;
    for(int i=0; i<PFJet_Track_idx.size();i++){
      if(PFJet_Track_idx.at(i)<Ntp->NTracks() && PFJet_Track_idx.at(i)>=0){
        track_p4.push_back(Ntp->Track_p4(PFJet_Track_idx.at(i)));
        track_charge.push_back(Ntp->Track_charge(PFJet_Track_idx.at(i)));
      }
    }
    for(int i=0; i<track_p4.size();i++){
      double JetCharge(0);
      int ntracks(0);
      for(int j=0; j<track_p4.size();j++){
        if(track_p4.at(i).Pt()<=track_p4.at(j).Pt()){
          ntracks++;
          JetCharge+=track_charge.at(j);
        }
      }
      if((ntracks==3 && track_p4.size()>=3) || (ntracks==track_p4.size() && track_p4.size()<3)){
        value.at(charge)=JetCharge+Ntp->Muon_Charge(mu_idx);
        break;
      }
    }
    pass.at(charge)=(fabs(value.at(charge)-cut.at(charge))<0.5);
  }
  else{
    value.at(charge)=-10.0;
    pass.at(charge)=false;
  }

  /*  if(jet_idx!=999){
    value.at(MaxTracksinJet)=Ntp->PFJet_Track_idx(jet_idx).size();
    pass.at(MaxTracksinJet)=cut.at(MaxTracksinJet)>value.at(MaxTracksinJet); 
    value.at(MinTracksinJet)=Ntp->PFJet_Track_idx(jet_idx).size();
    pass.at(MinTracksinJet)=cut.at(MinTracksinJet)<value.at(MinTracksinJet);
    std::vector<int> PFJet_Track_idx=Ntp->PFJet_Track_idx(jet_idx);
    for(int i=0; i<PFJet_Track_idx.size();i++){
      dist.at(MinTracksinJet).push_back(Ntp->Track_p4(PFJet_Track_idx.at(i)).Pt());
    }
  }

  // Momentum weighted charge is slow only do if nminus1
  // must be last cut
  value.at(charge)=0;
  if(jet_idx!=999 && mu_idx!=999){
    std::vector<int> PFJet_Track_idx=Ntp->PFJet_Track_idx(jet_idx);
    std::vector<TLorentzVector> track_p4;
    std::vector<float> track_charge;
    for(int i=0; i<PFJet_Track_idx.size();i++){ 
      if(PFJet_Track_idx.at(i)<Ntp->NTracks() && PFJet_Track_idx.at(i)>=0){
	track_p4.push_back(Ntp->Track_p4(PFJet_Track_idx.at(i)));
	track_charge.push_back(Ntp->Track_charge(PFJet_Track_idx.at(i)));
      }
    }
    for(int i=0; i<track_p4.size();i++){
      double JetCharge(0);
      int ntracks(0);
      for(int j=0; j<track_p4.size();j++){
	if(track_p4.at(i).Pt()<=track_p4.at(j).Pt()){
	  ntracks++;
	  JetCharge+=track_charge.at(j);
	}
      }
      if((ntracks==3 && track_p4.size()>=3) || (ntracks==1 && track_p4.size()<3)){
	value.at(charge)=JetCharge+Ntp->Muon_Charge(mu_idx);
	break;
      }
    }
    pass.at(charge)=(fabs(value.at(charge)-cut.at(charge))<0.5);
  }  
  else{
    value.at(charge)=-10.0;
    pass.at(charge)=false;
  }

  */
  if(verbose)std::cout << "void  Ztotautau_ControlSample::doEvent() i" << std::endl;

  //////////////////////////////////////////////////////////////////////////////////
  // QCD Control sample
  if(!pass.at(charge) && fabs(value.at(charge))==1){
    if(Ntp->isData()){
      if(!HConfig.GetHisto(!Ntp->isData(),DataMCType::QCD,t)){ std::cout << "failed to find id "<< DataMCType::QCD <<std::endl; return;}
      pass.at(charge)=true;
    }
  }

  double wobs(1),w(1);
  if(!Ntp->isData()){
    w*=Ntp->PUWeightFineBins();
  }
  else{w=1;}
  if(verbose)std::cout << "w=" << w << " " << wobs << " " << w*wobs << std::endl;
  bool status=AnalysisCuts(t,w,wobs); 
 
  ///////////////////////////////////////////////////////////
  // Add plots

  if(!Ntp->isData()){
    for(unsigned int tau_idx=0;tau_idx<Ntp->NMCTaus();tau_idx++){
      double ptmin(9999);
      double ptmax(0);
      unsigned int ptmaxidx(0);
      if((Ntp->MCTau_JAK(0)==2 || Ntp->MCTau_JAK(1)==2) && Ntp->MCTau_p4(1).Pt()>20 && Ntp->MCTau_p4(0).Pt()>20){
 
       for(unsigned int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++){
          if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::pi_plus) || abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::K_plus)){
            if(ptmin>Ntp->MCTauandProd_p4(tau_idx,i).Pt()) ptmin=Ntp->MCTauandProd_p4(tau_idx,i).Pt();
            if(ptmax<Ntp->MCTauandProd_p4(tau_idx,i).Pt()){ ptmax=Ntp->MCTauandProd_p4(tau_idx,i).Pt(); ptmaxidx=i;}
          }
        }
        for(unsigned int i=0; i<Ntp->NMCTauDecayProducts(tau_idx);i++){
          if(abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::pi_plus) || abs(Ntp->MCTauandProd_pdgid(tau_idx,i))==abs(PDGInfo::K_plus)){
            if(i!=ptmaxidx && Ntp->MCTauandProd_p4(tau_idx,i).DeltaR(Ntp->MCTauandProd_p4(tau_idx,ptmaxidx))<0.5){
              TauTrackdphitheta.at(t).Fill(sqrt(pow(Ntp->MCTauandProd_p4(tau_idx,i).DeltaPhi(Ntp->MCTauandProd_p4(tau_idx,ptmaxidx)),2.0)+
						pow(Ntp->MCTauandProd_p4(tau_idx,i).Theta()-Ntp->MCTauandProd_p4(tau_idx,ptmaxidx).Theta(),2.0)),w);
            }
          }
        }
        MinTauTrackPt.at(t).Fill(ptmin,w);
	MaxTauTrackPt.at(t).Fill(ptmax,w);
      }
    }
  }



  if(status){
    if(verbose)std::cout<<"MC type: " << Ntp->GetMCID() <<std::endl;
    NVtx.at(t).Fill(Ntp->NVtx(),w);
    unsigned int nGoodVtx=0;
    for(unsigned int i=0;i<Ntp->NVtx();i++){
      NTrackperVtx.at(t).Fill(Ntp->Vtx_Track_idx(i).size(),w);
      if(Ntp->isVtxGood(i))nGoodVtx++;
    }
    NGoodVtx.at(t).Fill(nGoodVtx,w);;
    if(mu_idx!=999) TagEtaPT.at(t).Fill(fabs(Ntp->Muon_p4(mu_idx).Eta()),Ntp->Muon_p4(mu_idx).Pt(),w);

 
    /////////////////////////////////////////
    //
    // Analyze Tau
    //
    if(jet_idx!=999){
      double mydr=0.5;
      for(unsigned i=0;i<Ntp->NKFTau();i++){
	if(Tools::dr(Ntp->KFTau_TauVis_p4(i),Ntp->PFJet_p4(jet_idx))<mydr){
	  KFTauCandPhi.at(t).Fill(Ntp->KFTau_TauFit_p4(i).Phi(),w);
	  if(Ntp->isGoodKFTau(i))KFSTauCandPhi.at(t).Fill(Ntp->KFTau_TauFit_p4(i).Phi(),w);
	  KFTauVCandPhi.at(t).Fill(Ntp->KFTau_TauVis_p4(i).Phi(),w);
          if(Ntp->isGoodKFTau(i))KFSTauVCandPhi.at(t).Fill(Ntp->KFTau_TauVis_p4(i).Phi(),w);

	}
      }
      for(unsigned i=0;i<Ntp->NPFTaus();i++){
	if(Tools::dr(Ntp->PFTau_p4(i),Ntp->PFJet_p4(jet_idx))<mydr){
	  PFTauCandPhi.at(t).Fill(Ntp->PFTau_p4(i).Phi(),w);
	  if(Ntp->PFTau_hpsDecayMode(i) == 10)PFMTauCandPhi.at(t).Fill(Ntp->PFTau_p4(i).Phi(),w);
	}
      }
    }

    unsigned int tau_idx(999);
    std::cout << "jet_idx " << jet_idx << " Ntp->NKFTau() " << Ntp->NKFTau() << std::endl;
      if(jet_idx!=999){
      double mydr=0.5;
      for(unsigned i=0;i<Ntp->NKFTau();i++){
	std::cout << "tau " <<  i;
        if(Ntp->isGoodKFTau(i)){
	  std::cout << " is good" << std::endl;
          if(Tools::dr(Ntp->KFTau_TauVis_p4(i),Ntp->PFJet_p4(jet_idx))<mydr){
            tau_idx=i;
	    mydr=Tools::dr(Ntp->KFTau_TauVis_p4(i),Ntp->PFJet_p4(jet_idx));
          }
        }
	else std::cout << "failed" << std::endl;
      }
      
      if(tau_idx!=999 ){
    

	TauCandFound.at(t).Fill(1,w);
        TauCandEtaPhi.at(t).Fill(fabs(Ntp->KFTau_TauFit_p4(tau_idx).Eta()),Ntp->KFTau_TauFit_p4(tau_idx).Pt(),w);
	std::cout << "Ntp->KFTau_TauFit_p4(tau_idx).Phi() " << Ntp->KFTau_TauFit_p4(tau_idx).Phi() << std::endl;
        TauCandPhi.at(t).Fill(Ntp->KFTau_TauFit_p4(tau_idx).Phi(),w);
        TauCandEta.at(t).Fill(Ntp->KFTau_TauFit_p4(tau_idx).Eta(),w);
        TauCandE.at(t).Fill(Ntp->KFTau_TauFit_p4(tau_idx).E(),w);
	TauCandP.at(t).Fill(Ntp->KFTau_TauFit_p4(tau_idx).P(),w);
	TauCandPT.at(t).Fill(Ntp->KFTau_TauFit_p4(tau_idx).Pt(),w);
	
	KFTau_Fit_chiprob.at(t).Fill(Ntp->KFTau_Fit_Chi2Prob(tau_idx),w);
	KFTau_Fit_a1mass.at(t).Fill(Ntp->KFTau_Fit_RefitVisibleMass(tau_idx),w);
	KFTau_Fit_chi2.at(t).Fill(Ntp->KFTau_Fit_chi2(tau_idx),w);
	KFTau_Fit_ndf.at(t).Fill(Ntp->KFTau_Fit_ndf(tau_idx),w);
	//KFTau_Fit_ambiguity.at(t).Fill(Ntp->KFTau_Fit_ambiguity(tau_idx),w);
	KFTau_Fit_csum.at(t).Fill(Ntp->KFTau_Fit_csum(tau_idx),w);
	KFTau_Fit_iterations.at(t).Fill(Ntp->KFTau_Fit_iterations(tau_idx),w);
	KFTau_Fit_SV_PV_significance.at(t).Fill(Ntp-> KFTau_Fit_SV_PV_significance(tau_idx),w);
	std::cout << "Ntp-> KFTau_Fit_SV_PV_significance(tau_idx) " << Ntp-> KFTau_Fit_SV_PV_significance(tau_idx) << std::endl;
      }
      else{
	TauCandFound.at(t).Fill(0.0,w);
      }
    }
  }
  if(id==DataMCType::Signal){
    unsigned int mcBoson_idx,mctau_idx;
    if(Ntp->hasSignalTauDecay(PDGInfo::Z0,mcBoson_idx,TauDecay::JAK_A1_3PI,mctau_idx)){
      TLorentzVector MCTau_LV(0,0,0,0),X_LV(0,0,0,0);
      for(unsigned int i=0; i<Ntp->NMCTauDecayProducts(mctau_idx);i++){
	if(abs(Ntp->MCTauandProd_pdgid(mctau_idx,i))==abs(PDGInfo::tau_minus) && Ntp->MCTau_JAK(mctau_idx)==TauDecay::JAK_A1_3PI){
	  MCTau_LV=Ntp->MCTauandProd_p4(mctau_idx,i);
	}
	else if(abs(Ntp->MCTauandProd_pdgid(mctau_idx,i))==abs(PDGInfo::pi_plus) ||
		abs(Ntp->MCTauandProd_pdgid(mctau_idx,i))==abs(PDGInfo::pi0)
		){
	  X_LV+=Ntp->MCTauandProd_p4(mctau_idx,i);
	  if(verbose)std::cout << "MC pi: " << Ntp->MCTauandProd_pdgid(mctau_idx,i) << " 4-vec: " << Ntp->MCTauandProd_p4(mctau_idx,i).E()
			       << " " << Ntp->MCTauandProd_p4(mctau_idx,i).Phi() 
			       << " " << Ntp->MCTauandProd_p4(mctau_idx,i).Theta() 
			       << " " << Ntp->MCTauandProd_p4(mctau_idx,i).M() << std::endl;
	}
      
      }
      unsigned int tau_idx(999);
      double mydr=4.0;
      for(unsigned i=0;i<Ntp->NKFTau();i++){
	if(Ntp->isGoodKFTau(i)){
	  if(Tools::dr(Ntp->KFTau_TauVis_p4(i),MCTau_LV)<mydr){
	    tau_idx=i;
	    mydr=Tools::dr(Ntp->KFTau_TauVis_p4(i),MCTau_LV);
	  }
	}
      }
      if(tau_idx!=999 ){
        TLorentzVector TauVisInput(0,0,0,0);
        TLorentzVector TauSolution(0,0,0,0);
	TLorentzVector TauSolution1;
	TLorentzVector TauSolution2;

      }
    }
  }
}
      
