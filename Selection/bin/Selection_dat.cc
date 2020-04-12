#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <iomanip>
#include <ctime>
#include <map>
#include <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TClass.h"
#include "TApplication.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TH2.h"
#include <TClonesArray.h>

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "CondFormats/BTauObjects/interface/BTagCalibration.h"
#include "CondTools/BTau/interface/BTagCalibrationReader.h"
//#include "../BtagUnc.hh"

//#include "BaconAna/DataFormats/interface/BaconAnaDefs.hh"
//#include "BaconAna/DataFormats/interface/TEventInfo.hh"
//#include "BaconAna/DataFormats/interface/TGenEventInfo.hh"
//#include "BaconAna/DataFormats/interface/TGenParticle.hh"
//#include "BaconAna/DataFormats/interface/TMuon.hh"
//#include "BaconAna/DataFormats/interface/TElectron.hh"
//#include "BaconAna/DataFormats/interface/TVertex.hh"
//#include "BaconAna/DataFormats/interface/TJet.hh"
//#include "BaconAna/DataFormats/interface/TAddJet.hh"
//#include "BaconAna/Utils/interface/TTrigger.hh"
//#include "BaconAna/DataFormats/interface/TLHEWeight.hh"

//#include "WVJJAna/Selection/interface/ScaleFactors.hh"
#include "WVJJAna/Selection/interface/NanoAOD_Data.hh"
#include "WVJJAna/Selection/interface/Utils.hh"
#include "WVJJAna/Selection/interface/METzCalculator.h"

int main (int ac, char** av) {

  std::string inputFile = av[1];
  std::string outputFile = av[2];
  int isMC = atoi(av[3]);
  int era = atoi(av[4]);

  const float MUON_MASS = 0.1056583745;
  const float ELE_MASS  = 0.000511;
  const float W_MASS = 80.385;
  //const float Z_MASS = 91.1876;

  //lepton cuts
  const float LEP_PT_VETO_CUT = 20;
  const float EL_PT_CUT = 35;
  const float EL_ETA_CUT = 2.5;
  const float MU_PT_CUT = 35;
  const float MU_ETA_CUT = 2.4;

  //ak8 jet cuts
  const float AK8_MIN_PT = 200;
  const float AK8_MAX_ETA = 2.4;
  const float AK8_MIN_SDM = 40;
  const float AK8_MAX_SDM = 150;

  //ak4 jet cuts
  //const float AK4_PT_VETO_CUT = 20;
  const float AK4_ETA_CUT = 2.4;
  const float AK4_PT_CUT = 30;
  const float AK4_JJ_MIN_M = 40.0;
  const float AK4_JJ_MAX_M = 150.0;
  const float VBF_MJJ_CUT= 500;

  //cleaning cuts
  const float AK8_LEP_DR_CUT = 1.0;
  const float AK4_AK8_DR_CUT = 0.8;
  const float AK4_DR_CUT = 0.3;

  //2016 csv tag thresholds
  //const float CSV_LOOSE_2016 = 0.5426;
  //const float CSV_MEDIUM_2016 = 0.8484;
  //const float CSV_TIGHT_2016 = 0.9535;

  //2017 csv tag thresholds
  //const float CSV_LOOSE_2017 = 0.5803;
  //const float CSV_MEDIUM_2017 = 0.8838;
  //const float CSV_TIGHT_2017 = 0.9693;
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X

  //ScaleFactors scaleFactor(era);

  std::vector<TLorentzVector> tightMuon;
  std::vector<TLorentzVector> tightEle;
  std::vector<TLorentzVector> goodAK4Jets;

  //std::string iHLTFile="${CMSSW_BASE}/src/BaconAna/DataFormats/data/HLTFile_25ns";
  //const std::string cmssw_base = getenv("CMSSW_BASE");
  //std::string cmssw_base_env = "${CMSSW_BASE}";
  //size_t start_pos = iHLTFile.find(cmssw_base_env);
  //if(start_pos != std::string::npos) {
  //  iHLTFile.replace(start_pos, cmssw_base_env.length(), cmssw_base);
  //}

  //const baconhep::TTrigger triggerMenu(iHLTFile);

  //PUPPI CORRECTIONS
  //2016?
  //TFile* file = TFile::Open("data/puppiCorr.root","READ");
  //TF1* puppisd_corrGEN      = (TF1*)file->Get("puppiJECcorr_gen");
  //TF1* puppisd_corrRECO_cen = (TF1*)file->Get("puppiJECcorr_reco_0eta1v3");
  //TF1* puppisd_corrRECO_for = (TF1*)file->Get("puppiJECcorr_reco_1v3eta2v5");

  //B TAGGING
  //BTagCalibration calib("csvv2", scaleFactor.GetBtagCSV().Data());
  //BTagCalibrationReader bTagReader(BTagEntry::OP_LOOSE,  // working point: can be OP_LOOSE, OP_MEDIUM, OP_TIGHT 
  //                                 "central",             // label for the central value (see the scale factor file)
  //                                 {"up","down"});        // vector of labels for systematics
  //bTagReader.load(calib, BTagEntry::FLAV_B, "comb");      // use the "comb" measurements for b-jets
  //bTagReader.load(calib, BTagEntry::FLAV_C, "comb");      // use the "comb" measurements for c-jets
  //bTagReader.load(calib, BTagEntry::FLAV_UDSG, "incl");   // use the "incl" measurements for light jets

  //JET CORRECTIONS
  JetCorrectorParameters paramAK4chs("data/Autumn18_V19_MC/Autumn18_V19_MC_Uncertainty_AK4PF.txt");
  JetCorrectorParameters paramAK8puppi("data/Autumn18_V19_MC/Autumn18_V19_MC_Uncertainty_AK8PFPuppi.txt");
  JetCorrectionUncertainty *fJetUnc_AK4chs = new JetCorrectionUncertainty(paramAK4chs);
  JetCorrectionUncertainty *fJetUnc_AK8puppi = new JetCorrectionUncertainty(paramAK8puppi);

  //
  //
  //   INPUT/OUTPUT
  //
  //

  TFile *of = new TFile(outputFile.c_str(),"RECREATE");
  TTree *ot = new TTree("Events","Events");
  WVJJData* WVJJTree = new WVJJData(ot);  
  TH1F *totalEvents = new TH1F("TotalEvents","TotalEvents",2,-1,1);

  //baconhep::TEventInfo *info   = new baconhep::TEventInfo();
  //baconhep::TGenEventInfo *gen = new baconhep::TGenEventInfo();
  //TClonesArray *vertexArr   = new TClonesArray("baconhep::TVertex");
  //TClonesArray *muonArr     = new TClonesArray("baconhep::TMuon");
  //TClonesArray *electronArr = new TClonesArray("baconhep::TElectron");
  //TClonesArray *AK4Arr      = new TClonesArray("baconhep::TJet");
  //TClonesArray *PuppiAK8Arr = new TClonesArray("baconhep::TJet");
  //TClonesArray *PuppiAK8AddArr = new TClonesArray("baconhep::TAddJet");
  //TClonesArray *lheWgtArr   = new TClonesArray("baconhep::TLHEWeight");

  TFile *f=0;
  TTree *t=0;

  std::ifstream ifs;
  ifs.open(inputFile.data());
  assert(ifs.is_open());
  std::string line;
  while (getline(ifs,line)) {
    std::stringstream ss(line);
    std::string filetoopen;
    ss >> filetoopen;

    f = TFile::Open(TString("root://cmseos.fnal.gov/")+TString(filetoopen),"read");
    t = (TTree *)f->Get("Events");
    //r = (TTree *)f->Get("Runs");
    if (t==NULL) continue;

    NanoAOD_Data NanoReader = NanoAOD_Data(t);

    for (uint i=0; i < t->GetEntries(); i++) {
      WVJJTree->clearVars();
      NanoReader.GetEntry(i);

      if (era==2018) {
	if(! ( NanoReader.HLT_IsoMu24 || NanoReader.HLT_IsoMu27 || NanoReader.HLT_IsoMu30 || NanoReader.HLT_Mu50 ||
	       NanoReader.HLT_Ele27_WPTight_Gsf || NanoReader.HLT_Ele28_WPTight_Gsf || NanoReader.HLT_Ele32_WPTight_Gsf ||
	       NanoReader.HLT_Ele35_WPTight_Gsf || NanoReader.HLT_Ele38_WPTight_Gsf || NanoReader.HLT_Ele40_WPTight_Gsf) )
	  continue;
	
      }
      /*
	if (era==2016) {
	if(!(triggerMenu.pass("HLT_IsoMu24_v*",info->triggerBits) || triggerMenu.pass("HLT_IsoTkMu24_v*",info->triggerBits) ||  
	     triggerMenu.pass("HLT_Ele27_WPTight_Gsf_v*",info->triggerBits))) continue;
      }
      if (era==2017) {
	if(!(triggerMenu.pass("HLT_IsoMu27_v*",info->triggerBits) || triggerMenu.pass("HLT_Ele32_WPTight_Gsf_v*",info->triggerBits) || 
	     triggerMenu.pass("HLT_Ele35_WPTight_Gsf_v*", info->triggerBits))) continue;
      }
      */
      tightMuon.clear();
      tightEle.clear();
      
      WVJJTree->run = NanoReader.run;
      WVJJTree->evt = NanoReader.event;
      WVJJTree->ls = NanoReader.luminosityBlock;
      
      //vertexArr->Clear();
      //vertexBr->GetEntry(i);
      
      WVJJTree->nPV = NanoReader.PV_npvsGood;

      WVJJTree->puWeight = 1.0;//scaleFactor.GetPUWeight(info->nPUmean, 0);

      // LEPTON SELECTION
      
      //muonArr->Clear();
      //muonBr->GetEntry(i);
      
      int nTightEle=0;
      int nTightMu=0;
      int nVetoEle=0;
      int nVetoMu=0;
      
      for (uint j=0; j < NanoReader.nMuon; j++) {
	//const baconhep::TMuon *mu = (baconhep::TMuon*)((*muonArr)[j]);
	
	if ( abs(NanoReader.Muon_eta[j]) > MU_ETA_CUT ) continue;
	if ( NanoReader.Muon_pt[j] < LEP_PT_VETO_CUT ) continue;
	
	if (!NanoReader.Muon_looseId[j]) continue;
	nVetoMu++;
	
	if (!NanoReader.Muon_tightId[j]) continue;
	if ( NanoReader.Muon_pt[j] < MU_PT_CUT ) continue;
	nTightMu++;
	tightMuon.push_back(TLorentzVector(0,0,0,0));
	tightMuon.back().SetPtEtaPhiM(NanoReader.Muon_pt[j], NanoReader.Muon_eta[j], 
				      NanoReader.Muon_phi[j], MUON_MASS);
	
	if ( NanoReader.Muon_pt[j] > WVJJTree->lep1_pt ) {
	  
	  WVJJTree->lep2_pt = WVJJTree->lep1_pt;
	  WVJJTree->lep2_eta = WVJJTree->lep1_eta;
	  WVJJTree->lep2_phi = WVJJTree->lep1_phi;
	  WVJJTree->lep2_m = WVJJTree->lep1_m;
	  WVJJTree->lep2_iso = WVJJTree->lep1_iso;
	  WVJJTree->lep2_q = WVJJTree->lep1_q;
	  
	  WVJJTree->lep1_pt = NanoReader.Muon_pt[j];
	  WVJJTree->lep1_eta = NanoReader.Muon_eta[j];
	  WVJJTree->lep1_phi = NanoReader.Muon_phi[j];
	  WVJJTree->lep1_m = MUON_MASS;
	  WVJJTree->lep1_iso = NanoReader.Muon_pfRelIso04_all[j];
	  WVJJTree->lep1_q = NanoReader.Muon_charge[j];
	  
	}
	else if ( NanoReader.Muon_pt[j] > WVJJTree->lep2_pt ) {

	  WVJJTree->lep2_pt = NanoReader.Muon_pt[j];
	  WVJJTree->lep2_eta = NanoReader.Muon_eta[j];
	  WVJJTree->lep2_phi = NanoReader.Muon_phi[j];
	  WVJJTree->lep2_m = MUON_MASS;
	  WVJJTree->lep2_iso = NanoReader.Muon_pfRelIso04_all[j];
	  WVJJTree->lep2_q = NanoReader.Muon_charge[j];

	}
      }
      
      //electronArr->Clear();
      //electronBr->GetEntry(i);
      
      for (uint j=0; j < NanoReader.nElectron; j++) {
	//const baconhep::TElectron *ele = (baconhep::TElectron*)((*electronArr)[j]);
	
	if ( abs(NanoReader.Electron_eta[j]) > EL_ETA_CUT ) continue;
	if ( NanoReader.Electron_pt[j] < LEP_PT_VETO_CUT ) continue;
	
	//cut-based ID Fall17 V2 (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
	if (NanoReader.Electron_cutBased[j]<2) continue;
	nVetoEle++;
	
	if (NanoReader.Electron_cutBased[j]<4) continue;
	nTightEle++;
	
	tightEle.push_back(TLorentzVector(0,0,0,0));
	tightEle.back().SetPtEtaPhiM(NanoReader.Electron_pt[j],NanoReader.Electron_eta[j],
				     NanoReader.Electron_phi[j],ELE_MASS);
	
	if ( NanoReader.Electron_pt[j] < EL_PT_CUT ) continue;
	
	//don't try to select electrons unless we don't already
	//have muons
	if (WVJJTree->lep1_m == MUON_MASS) continue;
	
	if ( NanoReader.Electron_pt[j] > WVJJTree->lep1_pt ) {
	  
	  WVJJTree->lep2_pt = WVJJTree->lep1_pt;
	  WVJJTree->lep2_eta = WVJJTree->lep1_eta;
	  WVJJTree->lep2_phi = WVJJTree->lep1_phi;
	  WVJJTree->lep2_m = WVJJTree->lep1_m;
	  WVJJTree->lep2_iso = WVJJTree->lep1_iso;
	  WVJJTree->lep2_q = WVJJTree->lep1_q;
	  
	  WVJJTree->lep1_pt = NanoReader.Electron_pt[j];
	  WVJJTree->lep1_eta = NanoReader.Electron_eta[j];
	  WVJJTree->lep1_phi = NanoReader.Electron_phi[j];
	  WVJJTree->lep1_m = ELE_MASS;
	  WVJJTree->lep1_iso = NanoReader.Electron_pfRelIso03_all[j];
	  WVJJTree->lep1_q = NanoReader.Electron_charge[j];
	  
	}
	else if ( NanoReader.Electron_pt[j] > WVJJTree->lep2_pt ) {
	  
	  WVJJTree->lep2_pt = NanoReader.Electron_pt[j];
	  WVJJTree->lep2_eta = NanoReader.Electron_eta[j];
	  WVJJTree->lep2_phi = NanoReader.Electron_phi[j];
	  WVJJTree->lep2_m = ELE_MASS;
	  WVJJTree->lep2_iso = NanoReader.Electron_pfRelIso03_all[j];
	  WVJJTree->lep2_q = NanoReader.Electron_charge[j];

	}
      }
      
      //check conditions
      if(!(WVJJTree->lep1_pt>0)) continue;
      if ((nTightMu+nTightEle)==0) continue; //no leptons with required ID
      if((nVetoEle+nVetoMu)>2) continue;
      if(nTightMu>0 && nVetoEle>0) continue;
      if(nTightEle>0 && nVetoMu>0) continue;
      if(nTightMu==1 && nVetoMu>1) continue;
      if(nTightEle==1 && nVetoEle>1) continue;
      
      //if (WVJJTree->lep1_pt < 0) continue; // no lepton canddiates
      //if (nVetoLeps>2) continue; // too many leptons
      //if ((nTightMu==1||nTightEle==1) && nVetoLeps>1) continue;
      
      //muon scale variations
      //if (WVJJTree->lep1_m == MU_MASS) {}
      
      //electron scale variations
      if (WVJJTree->lep1_m == ELE_MASS) {
	WVJJTree->lep1_pt_scaleUp = 1.01 * WVJJTree->lep1_pt;
	WVJJTree->lep1_pt_scaleDn = 0.99 * WVJJTree->lep1_pt;
	if (WVJJTree->lep2_pt>0) {
	  WVJJTree->lep2_pt_scaleUp = 1.01 * WVJJTree->lep2_pt;
	  WVJJTree->lep2_pt_scaleDn = 0.99 * WVJJTree->lep2_pt;
	}
      }
      
      if (WVJJTree->lep1_pt > 0 && WVJJTree->lep2_pt > 0) {
	
	TLorentzVector lep1(0,0,0,0);
	lep1.SetPtEtaPhiM( WVJJTree->lep1_pt, WVJJTree->lep1_eta, WVJJTree->lep1_phi, WVJJTree->lep1_m );
	
	TLorentzVector lep2(0,0,0,0);
	lep2.SetPtEtaPhiM( WVJJTree->lep2_pt, WVJJTree->lep2_eta, WVJJTree->lep2_phi, WVJJTree->lep2_m );
	
	TLorentzVector dilep = lep1+lep2;
	
	WVJJTree->dilep_m = dilep.M();
	WVJJTree->dilep_pt = dilep.Pt();
	WVJJTree->dilep_eta = dilep.Eta();
	WVJJTree->dilep_phi = dilep.Phi();
	
	//dilepton scale variations
	
      }
      
      //lepton ID/iso/trigger efficiencies
      //if (WVJJTree->lep1_m == ELE_MASS) {
      //WVJJTree->lep1_idEffWeight = scaleFactor.GetLeptonWeights(WVJJTree->lep1_pt, WVJJTree->lep1_eta, WVJJTree->lep1_m == ELE_MASS ? 11 : 13);
	//WVJJTree->lep1_idEffWeight = GetSFs_Lepton(WVJJTree->lep1_pt, WVJJTree->lep1_eta, hIDIsoEle);
	//WVJJTree->lep1_idEffWeight *= GetSFs_Lepton(WVJJTree->lep1_pt,WVJJTree->lep1_eta, hGSFCorrEle);
	//WVJJTree->lep1_idEffWeight *= GetSFs_Lepton(WVJJTree->lep1_pt,WVJJTree->lep1_eta, hTriggerEle);
	
      //if (WVJJTree->lep2_pt>0) {
	//WVJJTree->lep2_idEffWeight = scaleFactor.GetLeptonWeights(WVJJTree->lep2_pt, WVJJTree->lep2_eta, WVJJTree->lep1_m == ELE_MASS ? 11 : 13);
	  //WVJJTree->lep2_idEffWeight = GetSFs_Lepton(WVJJTree->lep2_pt, WVJJTree->lep2_eta, hIDIsoEle);
	  //WVJJTree->lep2_idEffWeight *= GetSFs_Lepton(WVJJTree->lep2_pt,WVJJTree->lep2_eta, hGSFCorrEle);
	  //WVJJTree->lep2_idEffWeight *= GetSFs_Lepton(WVJJTree->lep2_pt,WVJJTree->lep2_eta, hTriggerEle);
	  //do we even want the trigger eff for the subleading lepton?
      //}

      //else if (WVJJTree->lep1_m == MUON_MASS) {
      //	WVJJTree->lep1_idEffWeight = ScaleFactors.GetLeptonWeights(WVJJTree->lep1_pt, WVJJTree->lep1_eta, 13);
      //	//WVJJTree->lep1_idEffWeight = 1.0; // not implemented yet
      //	//WVJJTree->lep2_idEffWeight = 1.0; // not implemented yet
      //
      //	if (WVJJTree->lep2_pt>0) {
      //	  WVJJTree->lep2_idEffWeight = ScaleFactors.GetLeptonWeights(WVJJTree->lep2_pt, WVJJTree->lep2_eta, 13);
      //	}
      //}
      
      // MET
      
      if (NanoReader.PuppiMET_pt < 0) continue;
      WVJJTree->MET_2017raw = NanoReader.PuppiMET_pt;
      WVJJTree->MET_phi = NanoReader.PuppiMET_phi;
      
      TLorentzVector tempMet(0,0,0,0);
      tempMet.SetPxPyPzE(WVJJTree->MET_2017raw*TMath::Cos(WVJJTree->MET_phi), 
			 WVJJTree->MET_2017raw*TMath::Sin(WVJJTree->MET_phi), 
			 0.0, WVJJTree->MET_2017raw);
      
      TLorentzVector tempMet_Up = tempMet;
      TLorentzVector tempMet_Dn = tempMet;
      
      //AK4Arr->Clear();
      //AK4Br->GetEntry(i);
      
      //2017 only correction for ECAL noise
      for (uint j=0; j<NanoReader.nJet; j++) {
	//const baconhep::TJet *ak4jet = (baconhep::TJet*)((*AK4Arr)[j]);
	//if(ak4jet->pt < AK4_PT_VETO_CUT) continue;
	
	float jecUnc = GetJECunc(NanoReader.Jet_pt[j], NanoReader.Jet_eta[j], fJetUnc_AK4chs);
	
	//HARDCODED
	//if (era==2017 && ak4jet->ptRaw < 50 && abs(ak4jet->eta)>2.65 && abs(ak4jet->eta)<3.139) {
	//  TLorentzVector tempRawJet(0,0,0,0);
	//  tempRawJet.SetPtEtaPhiM(ak4jet->ptRaw, ak4jet->eta, ak4jet->phi, ak4jet->mass);
	//  
	//  tempMet+=tempRawJet;
	//  tempMet_Up+=tempRawJet;
	//  tempMet_Dn+=tempRawJet;
	//}
	
	//else {
	TLorentzVector tempJet(0,0,0,0);
	tempJet.SetPtEtaPhiM(NanoReader.Jet_pt[j], NanoReader.Jet_eta[j],
			     NanoReader.Jet_phi[j], NanoReader.Jet_mass[j]);
	
	TLorentzVector tempJetVar(0,0,0,0);
	tempJetVar.SetPtEtaPhiM(NanoReader.Jet_pt[j]*(1.0+jecUnc), NanoReader.Jet_eta[j],
				NanoReader.Jet_phi[j], NanoReader.Jet_mass[j]);
	tempMet_Up += tempJetVar - tempJet;
	
	tempJetVar.SetPtEtaPhiM(NanoReader.Jet_pt[j]*(1.0-jecUnc), NanoReader.Jet_eta[j], 
				NanoReader.Jet_phi[j], NanoReader.Jet_mass[j]);
	tempMet_Dn += tempJetVar - tempJet;
	//}
      }
      
      WVJJTree->MET = tempMet.Pt();
      WVJJTree->MET_scaleUp = tempMet_Up.Pt();
      WVJJTree->MET_scaleDn = tempMet_Dn.Pt();
      
      if (WVJJTree->lep2_pt<0) {
	
	TLorentzVector lep1(0,0,0,0);
	lep1.SetPtEtaPhiM( WVJJTree->lep1_pt, WVJJTree->lep1_eta, WVJJTree->lep1_phi, WVJJTree->lep1_m );
	
	METzCalculator NeutrinoPz_type0;
	NeutrinoPz_type0.SetMET(tempMet);
	NeutrinoPz_type0.SetLepton(lep1);
	NeutrinoPz_type0.SetLeptonType(WVJJTree->lep1_m == ELE_MASS ? "el" : "mu");
	
	WVJJTree->neu_pz_type0 = NeutrinoPz_type0.Calculate();
	
	NeutrinoPz_type0.SetMET(tempMet_Up);
	WVJJTree->neu_pz_type0_scaleUp = NeutrinoPz_type0.Calculate();
	
	NeutrinoPz_type0.SetMET(tempMet_Dn);
	WVJJTree->neu_pz_type0_scaleDn = NeutrinoPz_type0.Calculate();
      
	TLorentzVector neutrino(0,0,0,0);
	neutrino.SetPxPyPzE(tempMet.Px(), tempMet.Py(), WVJJTree->neu_pz_type0, 
			    sqrt(tempMet.Px()*tempMet.Px()+tempMet.Py()*tempMet.Py()+ WVJJTree->neu_pz_type0* WVJJTree->neu_pz_type0));
	
	TLorentzVector bosLep = lep1+neutrino;
	
	WVJJTree->dilep_m = bosLep.M();
	WVJJTree->dilep_pt = bosLep.Pt();
	WVJJTree->dilep_eta = bosLep.Eta();
	WVJJTree->dilep_phi = bosLep.Phi();
	
	//dilepton scale variations
	
      }
      
      // AK8
      
      //PuppiAK8Arr->Clear();
      //puppiAK8Br->GetEntry(i);
      //PuppiAK8AddArr->Clear();
      //puppiAK8AddBr->GetEntry(i);
      
      float dmW = 3000.0;
      int nGoodFatJet=0;
      
      for (uint j=0; j<NanoReader.nFatJet; j++) {
	//const baconhep::TJet *ak8jet = (baconhep::TJet*)((*PuppiAK8Arr)[j]);
	//const baconhep::TAddJet *ak8addjet = (baconhep::TAddJet*)((*PuppiAK8AddArr)[j]);
	
	if ( NanoReader.FatJet_pt[j]< AK8_MIN_PT ||  fabs(NanoReader.FatJet_eta[j]) > AK8_MAX_ETA ) continue;
	
	float jecUnc = GetJECunc(NanoReader.FatJet_pt[j], NanoReader.FatJet_eta[j], fJetUnc_AK8puppi);
	
	if ( NanoReader.FatJet_msoftdrop[j] < AK8_MIN_SDM || NanoReader.FatJet_msoftdrop[j] > AK8_MAX_SDM ) continue;
	if ( fabs(NanoReader.FatJet_msoftdrop[j] - W_MASS) > dmW ) continue;
	
	bool isClean=true;
	//lepton cleaning
	for ( std::size_t k=0; k<tightEle.size(); k++) {
	  if (deltaR(tightEle.at(k).Eta(), tightEle.at(k).Phi(),
		     NanoReader.FatJet_eta[j], NanoReader.FatJet_phi[j]) < AK8_LEP_DR_CUT)
	    isClean = false;
	}
	for ( std::size_t k=0; k<tightMuon.size(); k++) {
	  if (deltaR(tightMuon.at(k).Eta(), tightMuon.at(k).Phi(),
		     NanoReader.FatJet_eta[j], NanoReader.FatJet_phi[j]) < AK8_LEP_DR_CUT)
	    isClean = false;
	}
	if ( isClean == false ) continue;
	
	WVJJTree->bos_PuppiAK8_m_sd0 = NanoReader.FatJet_msoftdrop[j];
	WVJJTree->bos_PuppiAK8_m_sd0_corr = NanoReader.FatJet_msoftdrop[j];
	WVJJTree->bos_PuppiAK8_tau2tau1 = NanoReader.FatJet_tau2[j]/NanoReader.FatJet_tau1[j];
	WVJJTree->bos_PuppiAK8_pt = NanoReader.FatJet_pt[j];
	WVJJTree->bos_PuppiAK8_eta = NanoReader.FatJet_eta[j];
	WVJJTree->bos_PuppiAK8_phi = NanoReader.FatJet_phi[j];
	
	WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleUp =  WVJJTree->bos_PuppiAK8_m_sd0_corr*(1.0+jecUnc);
	WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleDn =  WVJJTree->bos_PuppiAK8_m_sd0_corr*(1.0-jecUnc);
	WVJJTree->bos_PuppiAK8_pt_scaleUp =  WVJJTree->bos_PuppiAK8_pt*(1.0+jecUnc);
	WVJJTree->bos_PuppiAK8_pt_scaleDn =  WVJJTree->bos_PuppiAK8_pt*(1.0-jecUnc);

	dmW = fabs(NanoReader.FatJet_msoftdrop[j] - W_MASS);
	nGoodFatJet++;
      }
      
      goodAK4Jets.clear();
      //AK4Arr->Clear();
      //AK4Br->GetEntry(i);

      for (uint j=0; j<NanoReader.nJet; j++) {

        float jecUnc = GetJECunc(NanoReader.Jet_pt[j], NanoReader.Jet_eta[j], fJetUnc_AK4chs);
	
	//jet energy scale variations
	if ( NanoReader.Jet_pt[j] < AK4_PT_CUT && NanoReader.Jet_pt[j]*(1.0+jecUnc) < AK4_PT_CUT && 
	     NanoReader.Jet_pt[j]*(1.0-jecUnc) < AK4_PT_CUT) continue;
	//jet ID??

	//if (era==2017 && ak4jet->ptRaw < 50 && abs(ak4jet->eta)>2.65 && abs(ak4jet->eta)<3.139) continue;

	//https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
	if (NanoReader.Jet_eta[j]<2.4 && NanoReader.Jet_pt[j]>30) {
	  if (NanoReader.Jet_btagDeepB[j] > 0.1241) WVJJTree->nBtag_loose++;
	  if (NanoReader.Jet_btagDeepB[j] > 0.4184) WVJJTree->nBtag_medium++;
	  if (NanoReader.Jet_btagDeepB[j] > 0.7527) WVJJTree->nBtag_tight++;
	}

	//if (abs(ak4jet->eta)<2.4 && ak4jet->pt>30) {
	//  if (era==2016) {
	//    if (ak4jet->csv > CSV_LOOSE_2016) WVJJTree->nBtag_loose++;
	//    if (ak4jet->csv > CSV_MEDIUM_2016) WVJJTree->nBtag_medium++;
	//    if (ak4jet->csv > CSV_TIGHT_2016) WVJJTree->nBtag_tight++;
	//  }
	//  else if (era==2017) {
	//    if (ak4jet->csv > CSV_LOOSE_2017) WVJJTree->nBtag_loose++;
	//    if (ak4jet->csv > CSV_MEDIUM_2017) WVJJTree->nBtag_medium++;
	//    if (ak4jet->csv > CSV_TIGHT_2017) WVJJTree->nBtag_tight++;
	//  }
	//}
	
	bool isClean=true;
	// object cleaning
	
	if (nGoodFatJet>0) {
	  if (deltaR(WVJJTree->bos_PuppiAK8_eta, WVJJTree->bos_PuppiAK8_phi,
		     NanoReader.Jet_eta[j], NanoReader.Jet_phi[j]) < AK4_AK8_DR_CUT) {
	    isClean = false;
	  }
	}
	
	for ( std::size_t k=0; k<goodAK4Jets.size(); k++) {
	  if (deltaR(goodAK4Jets.at(k).Eta(), goodAK4Jets.at(k).Phi(),
		     NanoReader.Jet_eta[j], NanoReader.Jet_phi[j]) < AK4_DR_CUT) {
	    isClean = false;
	  }
	}
	for ( std::size_t k=0; k<tightEle.size(); k++) {
	  if (deltaR(tightEle.at(k).Eta(), tightEle.at(k).Phi(),
		     NanoReader.Jet_eta[j],   NanoReader.Jet_phi[j]) < AK4_DR_CUT) {
	    isClean = false;
	  }
	}
	for ( std::size_t k=0; k<tightMuon.size(); k++) {
	  if (deltaR(tightMuon.at(k).Eta(), tightMuon.at(k).Phi(),
		     NanoReader.Jet_eta[j],   NanoReader.Jet_phi[j]) < AK4_DR_CUT) {
	    isClean = false;
	  }
	}
	
	if ( isClean == false ) continue;
	
	goodAK4Jets.push_back(TLorentzVector(0,0,0,0));
	goodAK4Jets.back().SetPtEtaPhiM(NanoReader.Jet_pt[j], NanoReader.Jet_eta[j],
					NanoReader.Jet_phi[j], NanoReader.Jet_mass[j]);
	
      }
      
      int nGoodDijet=0;
      
      uint sel1=1000, sel2=1000;
      if (nGoodFatJet==0) {
	TLorentzVector tmpV1, tmpV2;
	dmW=3000.0;
	for (uint j=0; j<goodAK4Jets.size(); j++) {
	  if ( fabs(goodAK4Jets.at(j).Eta()) < AK4_ETA_CUT ) continue;
	  for(uint k=j+1; k<goodAK4Jets.size(); k++) {
	    if ( fabs(goodAK4Jets.at(k).Eta()) < AK4_ETA_CUT ) continue;
	    TLorentzVector tmpV=goodAK4Jets.at(j)+goodAK4Jets.at(k);
	    
	    if (tmpV.M()<AK4_JJ_MIN_M || tmpV.M()>AK4_JJ_MAX_M) continue;

	    if (fabs(tmpV.M()-W_MASS)>dmW) continue;
	      
	    WVJJTree->bos_j1_AK4_pt =  goodAK4Jets.at(j).Pt();
	    WVJJTree->bos_j1_AK4_eta = goodAK4Jets.at(j).Eta();
	    WVJJTree->bos_j1_AK4_phi = goodAK4Jets.at(j).Phi();
	    WVJJTree->bos_j1_AK4_m =   goodAK4Jets.at(j).M();
	    
	    WVJJTree->bos_j2_AK4_pt =  goodAK4Jets.at(k).Pt();
	    WVJJTree->bos_j2_AK4_eta = goodAK4Jets.at(k).Eta();
	    WVJJTree->bos_j2_AK4_phi = goodAK4Jets.at(k).Phi();
	    WVJJTree->bos_j2_AK4_m =   goodAK4Jets.at(k).M();
	    
	    WVJJTree->bos_AK4AK4_pt =  tmpV.Pt();
	    WVJJTree->bos_AK4AK4_eta = tmpV.Eta();
	    WVJJTree->bos_AK4AK4_phi = tmpV.Phi();
	    WVJJTree->bos_AK4AK4_m =   tmpV.M();
	    
	    sel1=j; sel2=k;
	    dmW=fabs(tmpV.M()-W_MASS);
	    nGoodDijet=1;
	      
	  }
	}
	
	if (nGoodDijet==0) continue;
	
	float jecUnc1 = GetJECunc(WVJJTree->bos_j1_AK4_pt, WVJJTree->bos_j1_AK4_eta, fJetUnc_AK4chs);
	float jecUnc2 = GetJECunc(WVJJTree->bos_j2_AK4_pt, WVJJTree->bos_j2_AK4_eta, fJetUnc_AK4chs);
	
	WVJJTree->bos_j1_AK4_pt_scaleUp = WVJJTree->bos_j1_AK4_pt*(1.0+jecUnc1);
	WVJJTree->bos_j1_AK4_pt_scaleDn = WVJJTree->bos_j1_AK4_pt*(1.0-jecUnc1);
	WVJJTree->bos_j1_AK4_m_scaleUp = WVJJTree->bos_j1_AK4_m*(1.0+jecUnc1);
	WVJJTree->bos_j1_AK4_m_scaleDn = WVJJTree->bos_j1_AK4_m*(1.0-jecUnc1);
	
	WVJJTree->bos_j2_AK4_pt_scaleUp = WVJJTree->bos_j2_AK4_pt*(1.0+jecUnc2);
	WVJJTree->bos_j2_AK4_pt_scaleDn = WVJJTree->bos_j2_AK4_pt*(1.0-jecUnc2);
	WVJJTree->bos_j2_AK4_m_scaleUp = WVJJTree->bos_j2_AK4_m*(1.0+jecUnc2);
	WVJJTree->bos_j2_AK4_m_scaleDn = WVJJTree->bos_j2_AK4_m*(1.0-jecUnc2);
	
	TLorentzVector tempBos1(0,0,0,0);
	TLorentzVector tempBos2(0,0,0,0);
	
	tempBos1.SetPtEtaPhiM(WVJJTree->bos_j1_AK4_pt_scaleUp, WVJJTree->bos_j1_AK4_eta,
			      WVJJTree->bos_j1_AK4_phi, WVJJTree->bos_j1_AK4_m_scaleUp);
	tempBos2.SetPtEtaPhiM(WVJJTree->bos_j2_AK4_pt_scaleUp, WVJJTree->bos_j2_AK4_eta,
			      WVJJTree->bos_j2_AK4_phi, WVJJTree->bos_j2_AK4_m_scaleUp);
	
	WVJJTree->bos_AK4AK4_pt_scaleUp = (tempBos1+tempBos2).Pt();
	WVJJTree->bos_AK4AK4_m_scaleUp = (tempBos1+tempBos2).M();
	
	tempBos1.SetPtEtaPhiM(WVJJTree->bos_j1_AK4_pt_scaleDn, WVJJTree->bos_j1_AK4_eta,
			      WVJJTree->bos_j1_AK4_phi, WVJJTree->bos_j1_AK4_m_scaleDn);
	tempBos2.SetPtEtaPhiM(WVJJTree->bos_j2_AK4_pt_scaleDn, WVJJTree->bos_j2_AK4_eta,
			      WVJJTree->bos_j2_AK4_phi, WVJJTree->bos_j2_AK4_m_scaleDn);
	
	WVJJTree->bos_AK4AK4_pt_scaleDn = (tempBos1+tempBos2).Pt();
	WVJJTree->bos_AK4AK4_m_scaleDn = (tempBos1+tempBos2).M();
	
      } //if (nGoodFatJet==0)
      
      //check we have a hadronic boson candidate
      if ( nGoodFatJet == 0 && nGoodDijet == 0 ) continue;
      
      float tmpMassMax = 0.0;
      int vbf1=-1, vbf2=-1;
      
      for (uint j=0; j<goodAK4Jets.size(); j++) {
	if (j==sel1 || j==sel2) continue;
	for(uint k=j+1; k<goodAK4Jets.size(); k++) {
	  if (k==sel1 || k==sel2) continue;
	  TLorentzVector tempVBF = goodAK4Jets.at(j) + goodAK4Jets.at(k);
	  //require 2 jets be in opposite hemispheres
	  if ( goodAK4Jets.at(j).Eta()*goodAK4Jets.at(k).Eta() > 0 ) continue; 
	  if ( tempVBF.M() < VBF_MJJ_CUT ) continue;
	  if ( tempVBF.M() < tmpMassMax ) continue;
	  tmpMassMax = tempVBF.M();
	  vbf1=j; vbf2=k;
	}
      }
      
      if (vbf1==-1 && vbf2==-1) continue;
      
      TLorentzVector tempVBF = goodAK4Jets.at(vbf1) + goodAK4Jets.at(vbf2);
      
      WVJJTree->vbf1_AK4_pt = goodAK4Jets.at(vbf1).Pt();
      WVJJTree->vbf1_AK4_eta = goodAK4Jets.at(vbf1).Eta();
      WVJJTree->vbf1_AK4_phi = goodAK4Jets.at(vbf1).Phi();
      WVJJTree->vbf1_AK4_m = goodAK4Jets.at(vbf1).M();
      
      WVJJTree->vbf2_AK4_pt = goodAK4Jets.at(vbf2).Pt();
      WVJJTree->vbf2_AK4_eta = goodAK4Jets.at(vbf2).Eta();
      WVJJTree->vbf2_AK4_phi = goodAK4Jets.at(vbf2).Phi();
      WVJJTree->vbf2_AK4_m = goodAK4Jets.at(vbf2).M();
      
      WVJJTree->vbf_pt = tempVBF.Pt();
      WVJJTree->vbf_eta = tempVBF.Eta();
      WVJJTree->vbf_phi = tempVBF.Phi();
      WVJJTree->vbf_m = tempVBF.M();
      
      WVJJTree->vbf_deta = abs( WVJJTree->vbf2_AK4_eta - WVJJTree->vbf1_AK4_eta );
      
      TLorentzVector tempVBF1(0,0,0,0);
      TLorentzVector tempVBF2(0,0,0,0);
      
      float jecUnc1 = GetJECunc(WVJJTree->vbf1_AK4_pt, WVJJTree->vbf1_AK4_eta, fJetUnc_AK4chs);
      float jecUnc2 = GetJECunc(WVJJTree->vbf2_AK4_pt, WVJJTree->vbf2_AK4_eta, fJetUnc_AK4chs);
      
      WVJJTree->vbf1_AK4_pt_scaleUp = WVJJTree->vbf1_AK4_pt*(1.0+jecUnc1);
      WVJJTree->vbf1_AK4_pt_scaleDn = WVJJTree->vbf1_AK4_pt*(1.0-jecUnc1);
      WVJJTree->vbf1_AK4_m_scaleUp = WVJJTree->vbf1_AK4_m*(1.0+jecUnc1);
      WVJJTree->vbf1_AK4_m_scaleDn = WVJJTree->vbf1_AK4_m*(1.0-jecUnc1);
      
      WVJJTree->vbf2_AK4_pt_scaleUp = WVJJTree->vbf2_AK4_pt*(1.0+jecUnc2);
      WVJJTree->vbf2_AK4_pt_scaleDn = WVJJTree->vbf2_AK4_pt*(1.0-jecUnc2);
      WVJJTree->vbf2_AK4_m_scaleUp = WVJJTree->vbf2_AK4_m*(1.0+jecUnc2);
      WVJJTree->vbf2_AK4_m_scaleDn = WVJJTree->vbf2_AK4_m*(1.0-jecUnc2);
      
      tempVBF1.SetPtEtaPhiM(WVJJTree->vbf1_AK4_pt_scaleUp, WVJJTree->vbf1_AK4_eta,
			    WVJJTree->vbf1_AK4_phi, WVJJTree->vbf1_AK4_m_scaleUp);
      tempVBF2.SetPtEtaPhiM(WVJJTree->vbf2_AK4_pt_scaleUp, WVJJTree->vbf2_AK4_eta,
			    WVJJTree->vbf2_AK4_phi, WVJJTree->vbf2_AK4_m_scaleUp);
      
      WVJJTree->vbf_pt_scaleUp = (tempVBF1+tempVBF2).Pt();
      WVJJTree->vbf_m_scaleUp = (tempVBF1+tempVBF2).M();
      
      tempVBF1.SetPtEtaPhiM(WVJJTree->vbf1_AK4_pt_scaleDn, WVJJTree->vbf1_AK4_eta,
			    WVJJTree->vbf1_AK4_phi, WVJJTree->vbf1_AK4_m_scaleDn);
      tempVBF2.SetPtEtaPhiM(WVJJTree->vbf2_AK4_pt_scaleDn, WVJJTree->vbf2_AK4_eta,
			    WVJJTree->vbf2_AK4_phi, WVJJTree->vbf2_AK4_m_scaleDn);
      
      WVJJTree->vbf_pt_scaleDn = (tempVBF1+tempVBF2).Pt();
      WVJJTree->vbf_m_scaleDn = (tempVBF1+tempVBF2).M();
      
      TLorentzVector bosHad(0,0,0,0), bosHad_up(0,0,0,0), bosHad_dn(0,0,0,0);
      TLorentzVector bosLep(0,0,0,0), bosLep_up(0,0,0,0), bosLep_dn(0,0,0,0);
      
      //boosted event
      if (WVJJTree->bos_PuppiAK8_m_sd0_corr > 0) {
	bosHad.SetPtEtaPhiM(WVJJTree->bos_PuppiAK8_pt, WVJJTree->bos_PuppiAK8_eta,
			    WVJJTree->bos_PuppiAK8_phi, WVJJTree->bos_PuppiAK8_m_sd0_corr);
	bosHad_up.SetPtEtaPhiM(WVJJTree->bos_PuppiAK8_pt_scaleUp, WVJJTree->bos_PuppiAK8_eta,
			       WVJJTree->bos_PuppiAK8_phi, WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleUp);
	bosHad_dn.SetPtEtaPhiM(WVJJTree->bos_PuppiAK8_pt_scaleDn, WVJJTree->bos_PuppiAK8_eta,
			       WVJJTree->bos_PuppiAK8_phi, WVJJTree->bos_PuppiAK8_m_sd0_corr_scaleDn);
      }
    //resolved event
      else if (WVJJTree->bos_AK4AK4_m > 0) {
	bosHad.SetPtEtaPhiM(WVJJTree->bos_AK4AK4_pt, WVJJTree->bos_AK4AK4_eta, 
			    WVJJTree->bos_AK4AK4_phi, WVJJTree->bos_AK4AK4_m);
	bosHad_up.SetPtEtaPhiM(WVJJTree->bos_AK4AK4_pt_scaleUp, WVJJTree->bos_AK4AK4_eta, 
			       WVJJTree->bos_AK4AK4_phi, WVJJTree->bos_AK4AK4_m_scaleUp);
	bosHad_dn.SetPtEtaPhiM(WVJJTree->bos_AK4AK4_pt_scaleDn, WVJJTree->bos_AK4AK4_eta, 
			       WVJJTree->bos_AK4AK4_phi, WVJJTree->bos_AK4AK4_m_scaleDn);
      }
      
      bosLep.SetPtEtaPhiM(WVJJTree->dilep_pt, WVJJTree->dilep_eta, WVJJTree->dilep_phi, WVJJTree->dilep_m);
      //lepton variations not implemented
      bosLep_up.SetPtEtaPhiM(WVJJTree->dilep_pt, WVJJTree->dilep_eta, WVJJTree->dilep_phi, WVJJTree->dilep_m);
      bosLep_dn.SetPtEtaPhiM(WVJJTree->dilep_pt, WVJJTree->dilep_eta, WVJJTree->dilep_phi, WVJJTree->dilep_m);
      
      TLorentzVector diBos = bosHad+bosLep;
      TLorentzVector diBos_up = bosHad_up+bosLep_up;
      TLorentzVector diBos_dn = bosHad_dn+bosLep_dn;
      
      WVJJTree->dibos_m = diBos.M();
      WVJJTree->dibos_pt = diBos.Pt();
      WVJJTree->dibos_eta = diBos.Eta();
      WVJJTree->dibos_phi = diBos.Phi();
      
      WVJJTree->dibos_m_scaleUp = diBos_up.M();
      WVJJTree->dibos_m_scaleDn = diBos_dn.M();
      WVJJTree->dibos_pt_scaleUp = diBos_up.Pt();
      WVJJTree->dibos_pt_scaleDn = diBos_dn.Pt();
      
      if (WVJJTree->lep2_pt < 0) {
	WVJJTree->bosCent = std::min( std::min(bosHad.Eta(), bosLep.Eta()) - std::min(WVJJTree->vbf1_AK4_eta, WVJJTree->vbf2_AK4_eta), 
				      std::max(WVJJTree->vbf1_AK4_eta, WVJJTree->vbf2_AK4_eta) - std::max(bosHad.Eta(), bosLep.Eta()) );
      }
      
      WVJJTree->zeppLep = bosLep.Eta() - 0.5*(WVJJTree->vbf1_AK4_eta + WVJJTree->vbf2_AK4_eta);
      WVJJTree->zeppHad = bosHad.Eta() - 0.5*(WVJJTree->vbf1_AK4_eta + WVJJTree->vbf2_AK4_eta);
      
      //if(lheBr) {
      //	lheWgtArr->Clear();
      //	for (int j=0; j<lheWgtArr->GetEntries(); j++) {
      //	  const baconhep::TLHEWeight *lhe = (baconhep::TLHEWeight*)((*lheWgtArr)[j]);
      //	  WVJJTree->LHEWeight[i] = lhe->weight;
      //	}
      //}
      //
      ////need to add photon part
      //for (int j=0; j<AK4Arr->GetEntries(); j++) {
      //	const baconhep::TJet *ak4jet = (baconhep::TJet*)((*AK4Arr)[j]);
      //	if (abs(ak4jet->eta) > 2.0 && abs(ak4jet->eta)<3.0 && ak4jet->pt>30 && ak4jet->pt<500) {
      //	  //float pfRate = hL1PF_jetpt->GetBinContent(hL1PF_jetpt->FindBin(ak4jet->eta, ak4jet->pt));
      //	  WVJJTree->L1PFWeight *= scaleFactor.GetL1PFWeightJet(ak4jet->pt, ak4jet->eta);
      //	}
      //}
      
      ot->Fill();
    }
    delete t; delete f;
    t=0; f=0;
  }
  of->Write();
  of->Close();
  
  delete f;
  f=0; t=0;
  
  return 0;

}