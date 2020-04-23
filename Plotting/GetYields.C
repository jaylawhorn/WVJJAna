#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h> 
#include <TMath.h> 
#include <TSystem.h> 
#include <TFile.h> 
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#endif

struct Sample {
  string sampname;
  string filename;
  float  xsec;
  float  xsecCorr;
  double ntot;
  double wjets_el_wv;
  double wjets_el_wjj;
  double wjets_el_zv;
  double wjets_el_zjj;
  double wjets_mu_wv;
  double wjets_mu_wjj;
  double wjets_mu_zv;
  double wjets_mu_zjj;
  double ttbar_el_wv;
  double ttbar_el_wjj;
  double ttbar_el_zv;
  double ttbar_el_zjj;
  double ttbar_mu_wv;
  double ttbar_mu_wjj;
  double ttbar_mu_zv;
  double ttbar_mu_zjj;
  double loose_el_wv;
  double loose_el_wjj;
  double loose_el_zv;
  double loose_el_zjj;
  double loose_mu_wv;
  double loose_mu_wjj;
  double loose_mu_zv;
  double loose_mu_zjj;
};

Sample makeSample(string sname, string fname, float xsec, float xsecCorr) {

  Sample ret;
  ret.sampname=sname;
  ret.filename=fname;
  ret.xsec=xsec;
  ret.xsecCorr=xsecCorr;
  ret.ntot=0;
  ret.wjets_el_wv=0;
  ret.wjets_el_wjj=0;
  ret.wjets_el_zv=0;
  ret.wjets_el_zjj=0;
  ret.wjets_mu_wv=0;
  ret.wjets_mu_wjj=0;
  ret.wjets_mu_zv=0;
  ret.wjets_mu_zjj=0;
  ret.ttbar_el_wv=0;
  ret.ttbar_el_wjj=0;
  ret.ttbar_el_zv=0;
  ret.ttbar_el_zjj=0;
  ret.ttbar_mu_wv=0;
  ret.ttbar_mu_wjj=0;
  ret.ttbar_mu_zv=0;
  ret.ttbar_mu_zjj=0;
  ret.loose_el_wv=0;
  ret.loose_el_wjj=0;
  ret.loose_el_zv=0;
  ret.loose_el_zjj=0;
  ret.loose_mu_wv=0;
  ret.loose_mu_wjj=0;
  ret.loose_mu_zv=0;
  ret.loose_mu_zjj=0;

  return ret;
};

void Draw2(vector<Sample> samp1, Float_t lumi, Sample &master);

void GetYields() {

  string sampname;
  string filename;
  float  xsec;
  float  xsecCorr;
  int    fillColor;
  bool   doStack;

  vector<Sample> DataM;
  vector<Sample> DataE;
  vector<Sample> VBF_EWK;
  vector<Sample> VBF_QCD;
  vector<Sample> Top;
  vector<Sample> WJets;
  vector<Sample> DYJets;

  //Float_t lumi=35867.06;
  Float_t lumi=41530.0;
  //Float_t lumi=59000.0;
  //Float_t lumi=1000.0;

  ifstream ifs;
  ifs.open("files2017.txt");
  assert(ifs.is_open());
  string line;

  while(getline(ifs,line)) {

    if(line[0]=='#') continue;
    stringstream ss(line);

    ss >> sampname >> filename >> xsec >> xsecCorr;
    
    if (sampname == "VBF_EWK") {
      VBF_EWK.push_back(makeSample(sampname, filename, xsec, xsecCorr));
    }
    else if (sampname == "VBF_QCD") {
      VBF_QCD.push_back(makeSample(sampname, filename, xsec, xsecCorr));
    }
    else if (sampname == "Top") {
      Top.push_back(makeSample(sampname, filename, xsec, xsecCorr));
    }
    else if (sampname == "WJets") {
      WJets.push_back(makeSample(sampname, filename, xsec, xsecCorr));
    }
    else if (sampname == "DYJets") {
      DYJets.push_back(makeSample(sampname, filename, xsec, xsecCorr));
    }
    else if (sampname == "dataE") {
      DataE.push_back(makeSample(sampname, filename, xsec, xsecCorr));
    }
    else if (sampname == "dataM") {
      DataM.push_back(makeSample(sampname, filename, xsec, xsecCorr));
    }
  }

  Sample mVBF_EWK = makeSample("mVBF_EWK","mVBF_EWK",0.0,0.0);
  Sample mVBF_QCD = makeSample("mVBF_QCD","mVBF_QCD",0.0,0.0);
  Sample mTop = makeSample("mTop","mTop",0.0,0.0);
  Sample mWJets = makeSample("mWJets","mWJets",0.0,0.0);
  Sample mDYJets = makeSample("mDYJets","mDYJets",0.0,0.0);
  Sample mDataM = makeSample("mDataM","mDataM",0.0,0.0);

  Draw2(VBF_EWK,lumi,mVBF_EWK);
  Draw2(VBF_QCD,lumi,mVBF_QCD);
  Draw2(Top,lumi,mTop);
  Draw2(WJets,lumi,mWJets);
  Draw2(DYJets,lumi,mDYJets);
  Draw2(DataM,lumi,mDataM);

  cout << "2017,ele,WJets "  << endl;
  cout << "WV:   " << mDataM.wjets_el_wv << ", " << mWJets.wjets_el_wv+mDYJets.wjets_el_wv << ", " << mTop.wjets_el_wv << ", " << mVBF_QCD.wjets_el_wv << ", " << mVBF_EWK.wjets_el_wv << endl;
  cout << "Wjj   " << mDataM.wjets_el_wjj << ", " << mWJets.wjets_el_wjj+mDYJets.wjets_el_wjj << ", " << mTop.wjets_el_wjj << ", " << mVBF_QCD.wjets_el_wjj << ", " << mVBF_EWK.wjets_el_wjj << endl;
  cout << "ZV:   " << mDataM.wjets_el_zv << ", " << mWJets.wjets_el_zv+mDYJets.wjets_el_zv << ", " << mTop.wjets_el_zv << ", " << mVBF_QCD.wjets_el_zv << ", " << mVBF_EWK.wjets_el_zv << endl;
  cout << "Zjj:  " << mDataM.wjets_el_zjj << ", " << mWJets.wjets_el_zjj+mDYJets.wjets_el_zjj << ", " << mTop.wjets_el_zjj << ", " << mVBF_QCD.wjets_el_zjj << ", " << mVBF_EWK.wjets_el_zjj << endl;

  cout << "2017,ele,Loose "  << endl;
  cout << "WV:   " << mDataM.loose_el_wv << ", " << mWJets.loose_el_wv+mDYJets.loose_el_wv << ", " << mTop.loose_el_wv << ", " << mVBF_QCD.loose_el_wv << ", " << mVBF_EWK.loose_el_wv << endl;
  cout << "Wjj   " << mDataM.loose_el_wjj << ", " << mWJets.loose_el_wjj+mDYJets.loose_el_wjj << ", " << mTop.loose_el_wjj << ", " << mVBF_QCD.loose_el_wjj << ", " << mVBF_EWK.loose_el_wjj << endl;
  cout << "ZV:   " << mDataM.loose_el_zv << ", " << mWJets.loose_el_zv+mDYJets.loose_el_zv << ", " << mTop.loose_el_zv << ", " << mVBF_QCD.loose_el_zv << ", " << mVBF_EWK.loose_el_zv << endl;
  cout << "Zjj:  " << mDataM.loose_el_zjj << ", " << mWJets.loose_el_zjj+mDYJets.loose_el_zjj << ", " << mTop.loose_el_zjj << ", " << mVBF_QCD.loose_el_zjj << ", " << mVBF_EWK.loose_el_zjj << endl;

  cout << "2017,ele,TTbar "  << endl;
  cout << "WV:   " << mDataM.ttbar_el_wv << ", " << mWJets.ttbar_el_wv+mDYJets.ttbar_el_wv << ", " << mTop.ttbar_el_wv << ", " << mVBF_QCD.ttbar_el_wv << ", " << mVBF_EWK.ttbar_el_wv << endl;
  cout << "Wjj   " << mDataM.ttbar_el_wjj << ", " << mWJets.ttbar_el_wjj+mDYJets.ttbar_el_wjj << ", " << mTop.ttbar_el_wjj << ", " << mVBF_QCD.ttbar_el_wjj << ", " << mVBF_EWK.ttbar_el_wjj << endl;
  cout << "ZV:   " << mDataM.ttbar_el_zv << ", " << mWJets.ttbar_el_zv+mDYJets.ttbar_el_zv << ", " << mTop.ttbar_el_zv << ", " << mVBF_QCD.ttbar_el_zv << ", " << mVBF_EWK.ttbar_el_zv << endl;
  cout << "Zjj:  " << mDataM.ttbar_el_zjj << ", " << mWJets.ttbar_el_zjj+mDYJets.ttbar_el_zjj << ", " << mTop.ttbar_el_zjj << ", " << mVBF_QCD.ttbar_el_zjj << ", " << mVBF_EWK.ttbar_el_zjj << endl;

  cout << "2017,mu,WJets "  << endl;
  cout << "WV:   " << mDataM.wjets_mu_wv << ", " << mWJets.wjets_mu_wv+mDYJets.wjets_mu_wv << ", " << mTop.wjets_mu_wv << ", " << mVBF_QCD.wjets_mu_wv << ", " << mVBF_EWK.wjets_mu_wv << endl;
  cout << "Wjj   " << mDataM.wjets_mu_wjj << ", " << mWJets.wjets_mu_wjj+mDYJets.wjets_mu_wjj << ", " << mTop.wjets_mu_wjj << ", " << mVBF_QCD.wjets_mu_wjj << ", " << mVBF_EWK.wjets_mu_wjj << endl;
  cout << "ZV:   " << mDataM.wjets_mu_zv << ", " << mWJets.wjets_mu_zv+mDYJets.wjets_mu_zv << ", " << mTop.wjets_mu_zv << ", " << mVBF_QCD.wjets_mu_zv << ", " << mVBF_EWK.wjets_mu_zv << endl;
  cout << "Zjj:  " << mDataM.wjets_mu_zjj << ", " << mWJets.wjets_mu_zjj+mDYJets.wjets_mu_zjj << ", " << mTop.wjets_mu_zjj << ", " << mVBF_QCD.wjets_mu_zjj << ", " << mVBF_EWK.wjets_mu_zjj << endl;

  cout << "2017,mu,Loose "  << endl;
  cout << "WV:   " << mDataM.loose_mu_wv << ", " << mWJets.loose_mu_wv+mDYJets.loose_mu_wv << ", " << mTop.loose_mu_wv << ", " << mVBF_QCD.loose_mu_wv << ", " << mVBF_EWK.loose_mu_wv << endl;
  cout << "Wjj   " << mDataM.loose_mu_wjj << ", " << mWJets.loose_mu_wjj+mDYJets.loose_mu_wjj << ", " << mTop.loose_mu_wjj << ", " << mVBF_QCD.loose_mu_wjj << ", " << mVBF_EWK.loose_mu_wjj << endl;
  cout << "ZV:   " << mDataM.loose_mu_zv << ", " << mWJets.loose_mu_zv+mDYJets.loose_mu_zv << ", " << mTop.loose_mu_zv << ", " << mVBF_QCD.loose_mu_zv << ", " << mVBF_EWK.loose_mu_zv << endl;
  cout << "Zjj:  " << mDataM.loose_mu_zjj << ", " << mWJets.loose_mu_zjj+mDYJets.loose_mu_zjj << ", " << mTop.loose_mu_zjj << ", " << mVBF_QCD.loose_mu_zjj << ", " << mVBF_EWK.loose_mu_zjj << endl;

  cout << "2017,mu,TTbar "  << endl;
  cout << "WV:   " << mDataM.ttbar_mu_wv << ", " << mWJets.ttbar_mu_wv+mDYJets.ttbar_mu_wv << ", " << mTop.ttbar_mu_wv << ", " << mVBF_QCD.ttbar_mu_wv << ", " << mVBF_EWK.ttbar_mu_wv << endl;
  cout << "Wjj   " << mDataM.ttbar_mu_wjj << ", " << mWJets.ttbar_mu_wjj+mDYJets.ttbar_mu_wjj << ", " << mTop.ttbar_mu_wjj << ", " << mVBF_QCD.ttbar_mu_wjj << ", " << mVBF_EWK.ttbar_mu_wjj << endl;
  cout << "ZV:   " << mDataM.ttbar_mu_zv << ", " << mWJets.ttbar_mu_zv+mDYJets.ttbar_mu_zv << ", " << mTop.ttbar_mu_zv << ", " << mVBF_QCD.ttbar_mu_zv << ", " << mVBF_EWK.ttbar_mu_zv << endl;
  cout << "Zjj:  " << mDataM.ttbar_mu_zjj << ", " << mWJets.ttbar_mu_zjj+mDYJets.ttbar_mu_zjj << ", " << mTop.ttbar_mu_zjj << ", " << mVBF_QCD.ttbar_mu_zjj << ", " << mVBF_EWK.ttbar_mu_zjj << endl;


}
void Draw2(vector<Sample> samp1, Float_t lumi, Sample &master) {

  const float MUON_MASS = 0.1056583745;
  const float ELE_MASS  = 0.000511;

  Float_t genWeight=1, puWeight=1, lep1_idEffWeight=1, lep2_idEffWeight=1;
  Float_t L1PFWeight=1;
  Float_t lep1_pt=0, lep1_eta=0, lep1_phi=0, lep1_m=0, lep1_q=0; 
  Float_t lep2_pt=0, lep2_eta=0, lep2_phi=0, lep2_m=0, lep2_q=0;
  Float_t dilep_pt=0, dilep_eta=0, dilep_phi=0, dilep_m=0;
  Float_t MET=0, MET_2017raw=0, pfMET_Corr_phi=0, neu_pz_type0=0;
  Float_t bos_j1_AK4_pt=0, bos_j1_AK4_eta=0, bos_j1_AK4_phi=0, bos_j1_AK4_m=0; 
  Float_t bos_j2_AK4_pt=0, bos_j2_AK4_eta=0, bos_j2_AK4_phi=0, bos_j2_AK4_m=0;
  Float_t bos_AK4AK4_pt=0, bos_AK4AK4_eta=0, bos_AK4AK4_phi=0, bos_AK4AK4_m=0;
  Float_t bos_PuppiAK8_pt=0, bos_PuppiAK8_eta=0, bos_PuppiAK8_phi=0, bos_PuppiAK8_m_sd0_corr=0;
  Float_t bos_PuppiAK8_tau2tau1=0;
  Float_t vbf1_AK4_pt=0, vbf1_AK4_eta=0, vbf1_AK4_phi=0, vbf1_AK4_m=0;
  Float_t vbf2_AK4_pt=0, vbf2_AK4_eta=0, vbf2_AK4_phi=0, vbf2_AK4_m=0;
  Float_t vbf_m=0, vbf_deta=0, dibos_m=0;
  Int_t nBtag_loose=0;
  Int_t nJet30=0, nJet50;
  Float_t bosCent=0, zeppLep=0, zeppHad=0;

  for (uint xx=0; xx<samp1.size(); xx++) {
    cout << samp1.at(xx).filename << endl;

    TFile* infile = new TFile(TString(samp1.at(xx).filename), "READ");   assert(infile);
    TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

    int nTotal=1, nNeg=1;
    if (!(samp1.at(xx).sampname=="dataM" || samp1.at(xx).sampname=="dataE")) {
      TH1F* hTotEvents = (TH1F*) infile->Get("TotalEvents"); assert(hTotEvents);
      samp1.at(xx).ntot=hTotEvents->GetBinContent(2);
      samp1.at(xx).ntot-=2*hTotEvents->GetBinContent(1);
    }
    else {samp1.at(xx).ntot=1;}


    intree->SetBranchAddress("genWeight",&genWeight);
    intree->SetBranchAddress("puWeight",&puWeight);
    intree->SetBranchAddress("L1PFWeight",&L1PFWeight);
    intree->SetBranchAddress("lep1_idEffWeight",&lep1_idEffWeight);
    intree->SetBranchAddress("lep2_idEffWeight",&lep2_idEffWeight);
    
    intree->SetBranchAddress("lep1_pt", &lep1_pt);
    intree->SetBranchAddress("lep1_eta", &lep1_eta);
    intree->SetBranchAddress("lep1_phi", &lep1_phi);
    intree->SetBranchAddress("lep1_m", &lep1_m);
    intree->SetBranchAddress("lep1_q", &lep1_q);
    intree->SetBranchAddress("lep2_pt", &lep2_pt);
    intree->SetBranchAddress("lep2_eta", &lep2_eta);
    intree->SetBranchAddress("lep2_phi", &lep2_phi);
    intree->SetBranchAddress("lep2_m", &lep2_m);
    intree->SetBranchAddress("lep2_q", &lep2_q);
    
    intree->SetBranchAddress("dilep_pt", &dilep_pt);
    intree->SetBranchAddress("dilep_eta", &dilep_eta);
    intree->SetBranchAddress("dilep_phi", &dilep_phi);
    intree->SetBranchAddress("dilep_m", &dilep_m);
    
    intree->SetBranchAddress("MET", &MET);
    intree->SetBranchAddress("MET_phi", &pfMET_Corr_phi);
    intree->SetBranchAddress("neu_pz_type0", &neu_pz_type0);
    
    intree->SetBranchAddress("bos_j1_AK4_pt",  &bos_j1_AK4_pt);
    intree->SetBranchAddress("bos_j1_AK4_eta", &bos_j1_AK4_eta);
    intree->SetBranchAddress("bos_j1_AK4_phi", &bos_j1_AK4_phi);
    intree->SetBranchAddress("bos_j1_AK4_m",   &bos_j1_AK4_m);
    intree->SetBranchAddress("bos_j2_AK4_pt",  &bos_j2_AK4_pt);
    intree->SetBranchAddress("bos_j2_AK4_eta", &bos_j2_AK4_eta);
    intree->SetBranchAddress("bos_j2_AK4_phi", &bos_j2_AK4_phi);
    intree->SetBranchAddress("bos_j2_AK4_m",   &bos_j2_AK4_m);

    intree->SetBranchAddress("bos_AK4AK4_pt",   &bos_AK4AK4_pt);
    intree->SetBranchAddress("bos_AK4AK4_eta",   &bos_AK4AK4_eta);
    intree->SetBranchAddress("bos_AK4AK4_phi",   &bos_AK4AK4_phi);
    intree->SetBranchAddress("bos_AK4AK4_m",   &bos_AK4AK4_m);

    intree->SetBranchAddress("bos_PuppiAK8_pt",  &bos_PuppiAK8_pt);
    intree->SetBranchAddress("bos_PuppiAK8_eta", &bos_PuppiAK8_eta);
    intree->SetBranchAddress("bos_PuppiAK8_phi", &bos_PuppiAK8_phi);
    intree->SetBranchAddress("bos_PuppiAK8_m_sd0_corr", &bos_PuppiAK8_m_sd0_corr);
    intree->SetBranchAddress("bos_PuppiAK8_tau2tau1", &bos_PuppiAK8_tau2tau1);

    intree->SetBranchAddress("vbf1_AK4_pt", &vbf1_AK4_pt);
    intree->SetBranchAddress("vbf1_AK4_eta", &vbf1_AK4_eta);
    intree->SetBranchAddress("vbf1_AK4_phi", &vbf1_AK4_phi);
    intree->SetBranchAddress("vbf1_AK4_m", &vbf1_AK4_m);
    intree->SetBranchAddress("vbf2_AK4_pt", &vbf2_AK4_pt);
    intree->SetBranchAddress("vbf2_AK4_eta", &vbf2_AK4_eta);
    intree->SetBranchAddress("vbf2_AK4_phi", &vbf2_AK4_phi);
    intree->SetBranchAddress("vbf2_AK4_m", &vbf2_AK4_m);

    intree->SetBranchAddress("vbf_m", &vbf_m);

    intree->SetBranchAddress("dibos_m", &dibos_m);

    intree->SetBranchAddress("nJet30", &nJet30);
    intree->SetBranchAddress("nJet50", &nJet50);

    intree->SetBranchAddress("nBtag_loose", &nBtag_loose);

    intree->SetBranchAddress("bosCent", &bosCent);
    intree->SetBranchAddress("zeppLep", &zeppLep);
    intree->SetBranchAddress("zeppHad", &zeppHad);

    for (int i=0; i<intree->GetEntries(); i++) {
      intree->GetEntry(i);

      bool isEle=false, isResolved=false, isZ=false;
      bool isWJets=false, isTTbar=false, isLoose=false;

      if (bos_PuppiAK8_m_sd0_corr > 0 && bos_AK4AK4_m < 0) { isResolved=false; }
      else if (bos_PuppiAK8_m_sd0_corr < 0 && bos_AK4AK4_m > 0) { isResolved=true; }
      else {
	cout << "both or neither of resolved and boosted mass is defined" << endl;
      }
      
      if (lep1_m == ELE_MASS) { isEle=true; }
      else if (lep1_m == MUON_MASS) { isEle=false; }
      else {
	cout << "lepton is not electron or muon! skipping" << endl;
	continue;
      }

      if ( vbf_m < 500) continue;
      if ( fabs(vbf1_AK4_eta - vbf2_AK4_eta)<2.5) continue;
      if ( vbf1_AK4_pt<50 || vbf2_AK4_pt<50 ) continue;

      // resolved
      if (isResolved==true && (bos_j1_AK4_pt<30 || bos_j2_AK4_pt<30)) continue;

      // boosted
      if (isResolved==false && bos_PuppiAK8_pt<200) continue;
      if (isResolved==false && bos_PuppiAK8_tau2tau1>0.55) continue;

      //SR/CR
      if (nBtag_loose>0) {
	if (isResolved==true && (bos_AK4AK4_m>65 && bos_AK4AK4_m<105)) isTTbar=true;
	else if (isResolved==false && (bos_PuppiAK8_m_sd0_corr>65 && bos_PuppiAK8_m_sd0_corr<105)) isTTbar=true;
	else continue;
      }
      else if (nBtag_loose==0) {
	if (isResolved==true && (bos_AK4AK4_m>65 && bos_AK4AK4_m<105)) isLoose=true;
	else if (isResolved==true) isWJets=true;
	if (isResolved==false && (bos_PuppiAK8_m_sd0_corr>65 && bos_PuppiAK8_m_sd0_corr<105)) isLoose=true;
	else if (isResolved==false) isWJets=true;
      }
      
      //leptons
      if (lep2_pt>0) isZ=true;
      if (isEle==true && (lep1_pt<35 || abs(lep1_eta)>2.5 || (abs(lep1_eta)>1.4442 && abs(lep1_eta)<1.566))) continue;
      if (isEle==false && (lep1_pt<35 || abs(lep1_eta)>2.4)) continue;

      //Z(ll)
      if (isZ==true && (dilep_m < 81 || dilep_m > 101)) continue;
      if (isZ==true && isEle==true && (lep2_pt<35 || abs(lep2_eta)>2.5 || (abs(lep2_eta)>1.4442 && abs(lep2_eta)<1.566))) continue;
      if (isZ==true && isEle==false && (lep2_pt<35 || abs(lep2_eta)>2.4)) continue;
      if (isZ==true && (lep1_q*lep2_q)==1) continue;

      //W(lnu)
      if (isZ==false && MET<30) continue;

      float weight=samp1.at(xx).xsec*samp1.at(xx).xsecCorr*lumi*genWeight*puWeight*L1PFWeight;

      if (samp1.at(xx).sampname=="dataM" || samp1.at(xx).sampname=="dataE") weight=1;

      if      (isWJets==true && isEle==true && isResolved==false && isZ==false)  samp1.at(xx).wjets_el_wv+=weight;
      else if (isWJets==true && isEle==true && isResolved==false && isZ==true)   samp1.at(xx).wjets_el_zv+=weight;
      else if (isWJets==true && isEle==true && isResolved==true && isZ==false)   samp1.at(xx).wjets_el_wjj+=weight;
      else if (isWJets==true && isEle==true && isResolved==true && isZ==true)    samp1.at(xx).wjets_el_zjj+=weight;
      else if (isWJets==true && isEle==false && isResolved==false && isZ==false) samp1.at(xx).wjets_mu_wv+=weight;
      else if (isWJets==true && isEle==false && isResolved==false && isZ==true)  samp1.at(xx).wjets_mu_zv+=weight;
      else if (isWJets==true && isEle==false && isResolved==true && isZ==false)  samp1.at(xx).wjets_mu_wjj+=weight;
      else if (isWJets==true && isEle==false && isResolved==true && isZ==true)   samp1.at(xx).wjets_mu_zjj+=weight;

      else if (isTTbar==true && isEle==true && isResolved==false && isZ==false)  samp1.at(xx).ttbar_el_wv+=weight;
      else if (isTTbar==true && isEle==true && isResolved==false && isZ==true)   samp1.at(xx).ttbar_el_zv+=weight;
      else if (isTTbar==true && isEle==true && isResolved==true && isZ==false)   samp1.at(xx).ttbar_el_wjj+=weight;
      else if (isTTbar==true && isEle==true && isResolved==true && isZ==true)    samp1.at(xx).ttbar_el_zjj+=weight;
      else if (isTTbar==true && isEle==false && isResolved==false && isZ==false) samp1.at(xx).ttbar_mu_wv+=weight;
      else if (isTTbar==true && isEle==false && isResolved==false && isZ==true)  samp1.at(xx).ttbar_mu_zv+=weight;
      else if (isTTbar==true && isEle==false && isResolved==true && isZ==false)  samp1.at(xx).ttbar_mu_wjj+=weight;
      else if (isTTbar==true && isEle==false && isResolved==true && isZ==true)   samp1.at(xx).ttbar_mu_zjj+=weight;

      else if (isLoose==true && isEle==true && isResolved==false && isZ==false)  samp1.at(xx).loose_el_wv+=weight;
      else if (isLoose==true && isEle==true && isResolved==false && isZ==true)   samp1.at(xx).loose_el_zv+=weight;
      else if (isLoose==true && isEle==true && isResolved==true && isZ==false)   samp1.at(xx).loose_el_wjj+=weight;
      else if (isLoose==true && isEle==true && isResolved==true && isZ==true)    samp1.at(xx).loose_el_zjj+=weight;
      else if (isLoose==true && isEle==false && isResolved==false && isZ==false) samp1.at(xx).loose_mu_wv+=weight;
      else if (isLoose==true && isEle==false && isResolved==false && isZ==true)  samp1.at(xx).loose_mu_zv+=weight;
      else if (isLoose==true && isEle==false && isResolved==true && isZ==false)  samp1.at(xx).loose_mu_wjj+=weight;
      else if (isLoose==true && isEle==false && isResolved==true && isZ==true)   samp1.at(xx).loose_mu_zjj+=weight;

    }

    master.wjets_el_wv  += samp1.at(xx).wjets_el_wv /samp1.at(xx).ntot;
    master.wjets_el_zv  += samp1.at(xx).wjets_el_zv /samp1.at(xx).ntot;
    master.wjets_el_wjj += samp1.at(xx).wjets_el_wjj/samp1.at(xx).ntot;
    master.wjets_el_zjj += samp1.at(xx).wjets_el_zjj/samp1.at(xx).ntot;
    master.wjets_mu_wv  += samp1.at(xx).wjets_mu_wv /samp1.at(xx).ntot;
    master.wjets_mu_zv  += samp1.at(xx).wjets_mu_zv /samp1.at(xx).ntot;
    master.wjets_mu_wjj += samp1.at(xx).wjets_mu_wjj/samp1.at(xx).ntot;
    master.wjets_mu_zjj += samp1.at(xx).wjets_mu_zjj/samp1.at(xx).ntot;
    master.ttbar_el_wv  += samp1.at(xx).ttbar_el_wv /samp1.at(xx).ntot;
    master.ttbar_el_zv  += samp1.at(xx).ttbar_el_zv /samp1.at(xx).ntot;
    master.ttbar_el_wjj += samp1.at(xx).ttbar_el_wjj/samp1.at(xx).ntot;
    master.ttbar_el_zjj += samp1.at(xx).ttbar_el_zjj/samp1.at(xx).ntot;
    master.ttbar_mu_wv  += samp1.at(xx).ttbar_mu_wv /samp1.at(xx).ntot;
    master.ttbar_mu_zv  += samp1.at(xx).ttbar_mu_zv /samp1.at(xx).ntot;
    master.ttbar_mu_wjj += samp1.at(xx).ttbar_mu_wjj/samp1.at(xx).ntot;
    master.ttbar_mu_zjj += samp1.at(xx).ttbar_mu_zjj/samp1.at(xx).ntot;
    master.loose_el_wv  += samp1.at(xx).loose_el_wv /samp1.at(xx).ntot;
    master.loose_el_zv  += samp1.at(xx).loose_el_zv /samp1.at(xx).ntot;
    master.loose_el_wjj += samp1.at(xx).loose_el_wjj/samp1.at(xx).ntot;
    master.loose_el_zjj += samp1.at(xx).loose_el_zjj/samp1.at(xx).ntot;
    master.loose_mu_wv  += samp1.at(xx).loose_mu_wv /samp1.at(xx).ntot;
    master.loose_mu_zv  += samp1.at(xx).loose_mu_zv /samp1.at(xx).ntot;
    master.loose_mu_wjj += samp1.at(xx).loose_mu_wjj/samp1.at(xx).ntot;
    master.loose_mu_zjj += samp1.at(xx).loose_mu_zjj/samp1.at(xx).ntot;

  }

}
