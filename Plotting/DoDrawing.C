#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h> 
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
  int    nEvents;
  int    nNegEvents;
  int    fillColor;
  bool   doStack;
};

Sample makeSample(string sname, string fname, float xsec, float xsecCorr, 
		  int nEvents, int nNegEvents) {
  Sample ret;
  ret.sampname=sname;
  ret.filename=fname;
  ret.xsec=xsec;
  ret.xsecCorr=xsecCorr;
  ret.nEvents=nEvents;
  ret.nNegEvents=nNegEvents;

  return ret;
};

struct Binning {
  int nBinsMVBF;
  int nBinsDETA;
  int nBinsMVV;
  int nBinsTot;

  float MVBF_LE[15]; //!! careful with hard coding here
  float DETA_LE[15];
  float MVV_LE[15];

};

Binning makeBinning(int nBinsMVBF, int nBinsDETA, int nBinsMVV, 
		    float MVBF_LE[], float DETA_LE[], float MVV_LE[]) {
  
  Binning ret;
  ret.nBinsMVBF=nBinsMVBF;
  ret.nBinsDETA=nBinsDETA;
  ret.nBinsMVV=nBinsMVV;
  ret.nBinsTot=nBinsMVV*nBinsDETA*nBinsMVBF;

  for (int i=0; i<nBinsMVBF+1; i++) ret.MVBF_LE[i]=MVBF_LE[i];

  for (int i=0; i<nBinsDETA+1; i++) ret.DETA_LE[i]=DETA_LE[i];

  for (int i=0; i<nBinsMVV+1; i++) ret.MVV_LE[i]=MVV_LE[i];

  return ret;

}

void Draw2(vector<Sample> samp1, Binning bins, TString outfile, TString var);

void DoDrawing() {

  string sampname;
  string filename;
  float  xsec;
  float  xsecCorr;
  int    nEvents;
  int    nNegEvents;
  int    fillColor;
  bool   doStack;

  vector<Sample> DataM;
  vector<Sample> DataE;
  vector<Sample> VBF_EWK;
  vector<Sample> VBF_QCD;
  vector<Sample> Top;
  vector<Sample> WJets;
  vector<Sample> DYJets;

  ifstream ifs;
  ifs.open("files.txt");
  assert(ifs.is_open());
  string line;

  while(getline(ifs,line)) {

    if(line[0]=='#') continue;
    stringstream ss(line);

    ss >> sampname >> filename >> xsec >> xsecCorr >> nEvents >> nNegEvents;
    
    if (sampname == "VBF_EWK") {
      int index=VBF_EWK.size();
      VBF_EWK.push_back(makeSample(sampname, filename, xsec, xsecCorr, nEvents, nNegEvents));
    }
    else if (sampname == "VBF_QCD") {
      int index=VBF_QCD.size();
      VBF_QCD.push_back(makeSample(sampname, filename, xsec, xsecCorr, nEvents, nNegEvents));
    }
    else if (sampname == "Top") {
      int index=Top.size();
      Top.push_back(makeSample(sampname, filename, xsec, xsecCorr, nEvents, nNegEvents));
    }
    else if (sampname == "WJets") {
      int index=WJets.size();
      WJets.push_back(makeSample(sampname, filename, xsec, xsecCorr, nEvents, nNegEvents));
    }
    else if (sampname == "DYJets") {
      int index=DYJets.size();
      DYJets.push_back(makeSample(sampname, filename, xsec, xsecCorr, nEvents, nNegEvents));
    }
    else if (sampname == "dataE") {
      int index=DataE.size();
      DataE.push_back(makeSample(sampname, filename, xsec, xsecCorr, nEvents, nNegEvents));
    }
    else if (sampname == "dataM") {
      int index=DataM.size();
      DataM.push_back(makeSample(sampname, filename, xsec, xsecCorr, nEvents, nNegEvents));
    }
  }


  const int nBinsMVBF=4;
  const int nBinsDETA=4;
  const int nBinsMVV=7;
  float MVBF_LE[nBinsMVBF+1] = { 500, 600, 800, 1000, 3000 };
  //float MVBF_LE[nBinsMVBF+1] = { 600, 800, 1200, 3000 };
  float DETA_LE[nBinsDETA+1] = { 3.0, 4.0, 5.0, 6.0, 10.0 };
  //float DETA_LE[nBinsDETA+1] = { 4.0, 5.0, 6.0, 10.0 };
  float MVV_LE[nBinsMVV+1] = {150, 300, 450, 600, 1075, 1550, 2025, 2500};
  
  Binning bins=makeBinning(nBinsMVBF, nBinsDETA, nBinsMVV, 
			   MVBF_LE, DETA_LE, MVV_LE);

  Draw2(VBF_EWK,bins,"VBF_EWK_2016loose_drawM.root","nom");
  Draw2(VBF_QCD,bins,"VBF_QCD_2016loose_drawM.root","nom");
  Draw2(Top,bins,"Top_2016loose_drawM.root","nom");
  Draw2(WJets,bins,"WJets_2016loose_drawM.root","nom");
  Draw2(DYJets,bins,"DYJets_2016loose_drawM.root","nom");
  Draw2(DataM,bins,"DataM_2016loose_drawM.root","nom");
  //Draw2(DataE,bins,"DataE_2016loose_drawE.root","nom");


}
void Draw2(vector<Sample> samp1, Binning bins, TString outfile, TString var) {

  TFile *out;
  out = new TFile(outfile,"RECREATE");

  TString histName;

  //Wjj

  histName = Form("%s_mVV_Wjj",samp1.at(0).sampname.c_str());
  TH1D* mVV_Wjj  = new TH1D(histName, histName, bins.nBinsMVV, &bins.MVV_LE[0]);
  mVV_Wjj->Sumw2();
  mVV_Wjj->SetTitle(TString(samp1.at(0).sampname));
  mVV_Wjj->GetXaxis()->SetTitle("m(VV), Wjj");
  
  histName = Form("%s_dETA_Wjj",samp1.at(0).sampname.c_str());
  TH1D* dETA_Wjj = new TH1D(histName, histName, bins.nBinsDETA, &bins.DETA_LE[0]);
  dETA_Wjj->Sumw2();
  dETA_Wjj->SetTitle(TString(samp1.at(0).sampname));
  dETA_Wjj->GetXaxis()->SetTitle("dEta, Wjj");
  
  histName = Form("%s_mVBF_Wjj",samp1.at(0).sampname.c_str());
  TH1D* mVBF_Wjj = new TH1D(histName, histName, bins.nBinsMVBF, &bins.MVBF_LE[0]);
  mVBF_Wjj->Sumw2();
  mVBF_Wjj->SetTitle(TString(samp1.at(0).sampname));
  mVBF_Wjj->GetXaxis()->SetTitle("m(VBF), Wjj");

  histName = Form("%s_mVJ_Wjj",samp1.at(0).sampname.c_str());
  TH1D* mVJ_Wjj = new TH1D(histName, histName, 24, 30, 150);
  mVJ_Wjj->Sumw2();
  mVJ_Wjj->SetTitle(TString(samp1.at(0).sampname));
  mVJ_Wjj->GetXaxis()->SetTitle("m(W) had, Wjj");

  histName = Form("%s_MET_Wjj",samp1.at(0).sampname.c_str());
  TH1D* MET_Wjj = new TH1D(histName, histName, 20, 0, 2000);
  MET_Wjj->Sumw2();
  MET_Wjj->SetTitle(TString(samp1.at(0).sampname));
  MET_Wjj->GetXaxis()->SetTitle("MET, Wjj");

  histName = Form("%s_ptL_Wjj",samp1.at(0).sampname.c_str());
  TH1D* ptL_Wjj = new TH1D(histName, histName, 20, 0, 200);
  ptL_Wjj->Sumw2();
  ptL_Wjj->SetTitle(TString(samp1.at(0).sampname));
  ptL_Wjj->GetXaxis()->SetTitle("pT(lep), Wjj");

  //WV

  histName = Form("%s_mVV_WV",samp1.at(0).sampname.c_str());
  TH1D* mVV_WV = new TH1D(histName, histName, bins.nBinsMVV, &bins.MVV_LE[0]);
  mVV_WV->Sumw2();
  mVV_WV->SetTitle(TString(samp1.at(0).sampname));
  mVV_WV->GetXaxis()->SetTitle("m(VV), WV");

  histName = Form("%s_dETA_WV",samp1.at(0).sampname.c_str());
  TH1D* dETA_WV = new TH1D(histName, histName, bins.nBinsDETA, &bins.DETA_LE[0]);
  dETA_WV->Sumw2();
  dETA_WV->SetTitle(TString(samp1.at(0).sampname));
  dETA_WV->GetXaxis()->SetTitle("dETA, WV");

  histName = Form("%s_mVBF_WV",samp1.at(0).sampname.c_str());
  TH1D* mVBF_WV = new TH1D(histName, histName, bins.nBinsMVBF, &bins.MVBF_LE[0]);
  mVBF_WV->Sumw2();
  mVBF_WV->SetTitle(TString(samp1.at(0).sampname));
  mVBF_WV->GetXaxis()->SetTitle("m(VBF), WV");

  histName = Form("%s_mVJ_WV",samp1.at(0).sampname.c_str());
  TH1D* mVJ_WV = new TH1D(histName, histName, 24, 30, 150);
  mVJ_WV->Sumw2();
  mVJ_WV->SetTitle(TString(samp1.at(0).sampname));
  mVJ_WV->GetXaxis()->SetTitle("m(W) had, WV");

  histName = Form("%s_MET_WV",samp1.at(0).sampname.c_str());
  TH1D* MET_WV = new TH1D(histName, histName, 20, 0, 2000);
  MET_WV->Sumw2();
  MET_WV->SetTitle(TString(samp1.at(0).sampname));
  MET_WV->GetXaxis()->SetTitle("MET, WV");

  histName = Form("%s_ptL_WV",samp1.at(0).sampname.c_str());
  TH1D* ptL_WV = new TH1D(histName, histName, 20, 0, 200);
  ptL_WV->Sumw2();
  ptL_WV->SetTitle(TString(samp1.at(0).sampname));
  ptL_WV->GetXaxis()->SetTitle("pT(lep), WV");

  //Zjj

  histName = Form("%s_mVV_Zjj",samp1.at(0).sampname.c_str());
  TH1D* mVV_Zjj = new TH1D(histName, histName, bins.nBinsMVV, &bins.MVV_LE[0]);
  mVV_Zjj->Sumw2();
  mVV_Zjj->SetTitle(TString(samp1.at(0).sampname));
  mVV_Zjj->GetXaxis()->SetTitle("m(VV), Zjj"); 

  histName = Form("%s_dETA_Zjj",samp1.at(0).sampname.c_str());
  TH1D* dETA_Zjj = new TH1D(histName, histName, bins.nBinsDETA, &bins.DETA_LE[0]);
  dETA_Zjj->Sumw2();
  dETA_Zjj->SetTitle(TString(samp1.at(0).sampname));
  dETA_Zjj->GetXaxis()->SetTitle("dEta, Zjj");

  histName = Form("%s_mVBF_Zjj",samp1.at(0).sampname.c_str());
  TH1D* mVBF_Zjj = new TH1D(histName, histName, bins.nBinsMVBF, &bins.MVBF_LE[0]);
  mVBF_Zjj->Sumw2();
  mVBF_Zjj->SetTitle(TString(samp1.at(0).sampname));
  mVBF_Zjj->GetXaxis()->SetTitle("m(VBF), Zjj");

  histName = Form("%s_mVJ_Zjj",samp1.at(0).sampname.c_str());
  TH1D* mVJ_Zjj = new TH1D(histName, histName, 24, 30, 150);
  mVJ_Zjj->Sumw2();
  mVJ_Zjj->SetTitle(TString(samp1.at(0).sampname));
  mVJ_Zjj->GetXaxis()->SetTitle("m(Z) had, Zjj");

  histName = Form("%s_mVL_Zjj",samp1.at(0).sampname.c_str());
  TH1D* mVL_Zjj = new TH1D(histName, histName, 24, 30, 150);
  mVL_Zjj->Sumw2();
  mVL_Zjj->SetTitle(TString(samp1.at(0).sampname));
  mVL_Zjj->GetXaxis()->SetTitle("m(Z) lep, Zjj");

  //ZV

  histName = Form("%s_mVV_ZV",samp1.at(0).sampname.c_str());
  TH1D* mVV_ZV = new TH1D(histName, histName, bins.nBinsMVV, &bins.MVV_LE[0]);
  mVV_ZV->Sumw2();
  mVV_ZV->SetTitle(TString(samp1.at(0).sampname));
  mVV_ZV->GetXaxis()->SetTitle("m(VV), ZV");

  histName = Form("%s_dETA_ZV",samp1.at(0).sampname.c_str());
  TH1D* dETA_ZV = new TH1D(histName, histName, bins.nBinsDETA, &bins.DETA_LE[0]);
  dETA_ZV->Sumw2();
  dETA_ZV->SetTitle(TString(samp1.at(0).sampname));
  dETA_ZV->GetXaxis()->SetTitle("dEta, ZV");

  histName = Form("%s_mVBF_ZV",samp1.at(0).sampname.c_str());
  TH1D* mVBF_ZV = new TH1D(histName, histName, bins.nBinsMVBF, &bins.MVBF_LE[0]);
  mVBF_ZV->Sumw2();
  mVBF_ZV->SetTitle(TString(samp1.at(0).sampname));
  mVBF_Zjj->GetXaxis()->SetTitle("m(VBF), ZV");

  histName = Form("%s_mVJ_ZV",samp1.at(0).sampname.c_str());
  TH1D* mVJ_ZV = new TH1D(histName, histName, 24, 30, 150);
  mVJ_ZV->Sumw2();
  mVJ_ZV->SetTitle(TString(samp1.at(0).sampname));
  mVJ_ZV->GetXaxis()->SetTitle("m(Z) had, ZV");

  histName = Form("%s_mVL_ZV",samp1.at(0).sampname.c_str());
  TH1D* mVL_ZV = new TH1D(histName, histName, 24, 30, 150);
  mVL_ZV->Sumw2();
  mVL_ZV->SetTitle(TString(samp1.at(0).sampname));
  mVL_ZV->GetXaxis()->SetTitle("m(Z) lep, ZV");

  Float_t lumi=35867.06;
  //Float_t lumi=41530.0;

  Bool_t  isResolved;
  Int_t   type;
  Float_t genWeight=1, pu_Weight=1, btag0Wgt=1, id_eff_Weight=1, trig_eff_Weight=1;
  Float_t l_pt1=0, l_eta1=0, l_phi1=0, l_e1=0;
  Float_t l_pt2=0, l_eta2=0, l_phi2=0, l_e2=0;
  Float_t dilep_pt=0, dilep_eta=0, dilep_phi=0, dilep_m=0;
  Float_t pfMET_Corr=0, pfMET_Corr_phi=0;
  Float_t nu_pz_type0=0;

  Float_t AK4_vjet1_pt=0, AK4_vjet1_eta=0, AK4_vjet1_phi=0, AK4_vjet1_e=0;
  Float_t AK4_vjet2_pt=0, AK4_vjet2_eta=0, AK4_vjet2_phi=0, AK4_vjet2_e=0;
  Float_t ungroomed_PuppiAK8_jet_pt=0, ungroomed_PuppiAK8_jet_eta=0, ungroomed_PuppiAK8_jet_phi=0, ungroomed_PuppiAK8_jet_e=0;
  Float_t PuppiAK8_jet_mass_so_corr=0;

  Float_t vbf_maxpt_j1_pt=0, vbf_maxpt_j1_eta=0, vbf_maxpt_j1_phi=0, vbf_maxpt_j1_e=0, vbf_maxpt_j1_mass=0;
  Float_t vbf_maxpt_j2_pt=0, vbf_maxpt_j2_eta=0, vbf_maxpt_j2_phi=0, vbf_maxpt_j2_e=0, vbf_maxpt_j2_mass=0;

  Int_t   nBTagJet_loose=0;
  Float_t vbf_maxpt_jj_m=0, mass_lvj_type0_PuppiAK8=0, BosonCentrality_type0=0, ZeppenfeldWL_type0=0, ZeppenfeldWH=0;
  Float_t mass_llj_PuppiAK8=0;

  for (uint xx=0; xx<samp1.size(); xx++) {
    cout << samp1.at(xx).filename << endl;

    TFile* infile = new TFile(TString(samp1.at(xx).filename), "READ");   assert(infile);
    TTree* intree = (TTree*) infile->Get("otree"); assert(intree);

    intree->SetBranchAddress("isResolved", &isResolved);
    intree->SetBranchAddress("type", &type);
    intree->SetBranchAddress("genWeight",&genWeight);
    intree->SetBranchAddress("pu_Weight",&pu_Weight);
    intree->SetBranchAddress("btag0Wgt",&btag0Wgt);
    intree->SetBranchAddress("id_eff_Weight",&id_eff_Weight);
    intree->SetBranchAddress("trig_eff_Weight",&trig_eff_Weight);
    
    intree->SetBranchAddress("l_pt1", &l_pt1);
    intree->SetBranchAddress("l_eta1", &l_eta1);
    intree->SetBranchAddress("l_phi1", &l_phi1);
    intree->SetBranchAddress("l_e2", &l_e1);
    intree->SetBranchAddress("l_pt2", &l_pt2);
    intree->SetBranchAddress("l_eta2", &l_eta2);
    intree->SetBranchAddress("l_phi2", &l_phi2);
    intree->SetBranchAddress("l_e2", &l_e2);
    
    intree->SetBranchAddress("dilep_pt", &dilep_pt);
    intree->SetBranchAddress("dilep_eta", &dilep_eta);
    intree->SetBranchAddress("dilep_phi", &dilep_phi);
    intree->SetBranchAddress("dilep_m", &dilep_m);
    
    intree->SetBranchAddress("pfMET_Corr", &pfMET_Corr);
    intree->SetBranchAddress("pfMET_Corr_phi", &pfMET_Corr_phi);
    intree->SetBranchAddress("nu_pz_type0", &nu_pz_type0);
    
    intree->SetBranchAddress("AK8jet_sj1_pt", &AK4_vjet1_pt); //temporary fix
    intree->SetBranchAddress("AK4_vjet1_eta", &AK4_vjet1_eta);
    intree->SetBranchAddress("AK4_vjet1_phi", &AK4_vjet1_phi);
    intree->SetBranchAddress("AK4_vjet1_e", &AK4_vjet1_e);
    intree->SetBranchAddress("AK4_vjet2_pt", &AK4_vjet2_pt);
    intree->SetBranchAddress("AK4_vjet2_eta", &AK4_vjet2_eta);
    intree->SetBranchAddress("AK4_vjet2_phi", &AK4_vjet2_phi);
    intree->SetBranchAddress("AK4_vjet2_e", &AK4_vjet2_e);
    
    intree->SetBranchAddress("ungroomed_PuppiAK8_jet_pt", &ungroomed_PuppiAK8_jet_pt);
    intree->SetBranchAddress("ungroomed_PuppiAK8_jet_eta", &ungroomed_PuppiAK8_jet_eta);
    intree->SetBranchAddress("ungroomed_PuppiAK8_jet_phi", &ungroomed_PuppiAK8_jet_phi);
    intree->SetBranchAddress("ungroomed_PuppiAK8_jet_e", &ungroomed_PuppiAK8_jet_e);
    //missing W(jj) jet scale uncertainties
    if (var=="jesup") {
      intree->SetBranchAddress("ungroomed_PuppiAK8_jet_mass_jes_up", &PuppiAK8_jet_mass_so_corr);
      intree->SetBranchAddress("vbf_maxpt_j1_pt_jes_up", &vbf_maxpt_j1_pt);
      intree->SetBranchAddress("vbf_maxpt_j2_pt_jes_up", &vbf_maxpt_j2_pt);
      intree->SetBranchAddress("vbf_maxpt_jj_m_jes_up", &vbf_maxpt_jj_m);      
      intree->SetBranchAddress("mass_lvj_type0_PuppiAK8_jes_up", &mass_lvj_type0_PuppiAK8);
      intree->SetBranchAddress("mass_llj_PuppiAK8", &mass_llj_PuppiAK8);    
    }
    else if (var=="jesdn"){
      intree->SetBranchAddress("ungroomed_PuppiAK8_jet_mass_jes_dn", &PuppiAK8_jet_mass_so_corr);
      intree->SetBranchAddress("vbf_maxpt_j1_pt_jes_dn", &vbf_maxpt_j1_pt);
      intree->SetBranchAddress("vbf_maxpt_j2_pt_jes_dn", &vbf_maxpt_j2_pt);
      intree->SetBranchAddress("vbf_maxpt_jj_m_jes_dn", &vbf_maxpt_jj_m);      
      intree->SetBranchAddress("mass_lvj_type0_PuppiAK8_jes_dn", &mass_lvj_type0_PuppiAK8);
      intree->SetBranchAddress("mass_llj_PuppiAK8", &mass_llj_PuppiAK8);    
    }
    else {
      intree->SetBranchAddress("PuppiAK8_jet_mass_so_corr", &PuppiAK8_jet_mass_so_corr);
      intree->SetBranchAddress("vbf_maxpt_j1_pt", &vbf_maxpt_j1_pt);
      intree->SetBranchAddress("vbf_maxpt_j2_pt", &vbf_maxpt_j2_pt);
      intree->SetBranchAddress("vbf_maxpt_jj_m", &vbf_maxpt_jj_m);      
      intree->SetBranchAddress("mass_lvj_type0_PuppiAK8", &mass_lvj_type0_PuppiAK8);
      intree->SetBranchAddress("mass_llj_PuppiAK8", &mass_llj_PuppiAK8);    
    }
    
    intree->SetBranchAddress("vbf_maxpt_j1_eta", &vbf_maxpt_j1_eta);
    intree->SetBranchAddress("vbf_maxpt_j1_phi", &vbf_maxpt_j1_phi);
    intree->SetBranchAddress("vbf_maxpt_j1_e", &vbf_maxpt_j1_e);
    intree->SetBranchAddress("vbf_maxpt_j1_mass", &vbf_maxpt_j1_mass);

    intree->SetBranchAddress("vbf_maxpt_j2_eta", &vbf_maxpt_j2_eta);
    intree->SetBranchAddress("vbf_maxpt_j2_phi", &vbf_maxpt_j2_phi);
    intree->SetBranchAddress("vbf_maxpt_j2_e", &vbf_maxpt_j2_e);
    intree->SetBranchAddress("vbf_maxpt_j2_mass", &vbf_maxpt_j2_mass);
    
    intree->SetBranchAddress("nBTagJet_loose", &nBTagJet_loose);

    intree->SetBranchAddress("BosonCentrality_type0", &BosonCentrality_type0);
    intree->SetBranchAddress("ZeppenfeldWL_type0", &ZeppenfeldWL_type0);
    intree->SetBranchAddress("ZeppenfeldWH", &ZeppenfeldWH);    

    for (int i=0; i<intree->GetEntries(); i++) {
      intree->GetEntry(i);
      
      if ( !(nBTagJet_loose==0 &&
	     //((PuppiAK8_jet_mass_so_corr>65) && (PuppiAK8_jet_mass_so_corr<105)) &&
	     ((vbf_maxpt_j1_pt>30) && (vbf_maxpt_j2_pt>30))) ) continue;

      //if ( !(l_pt1>35 && (((type==0)&&(abs(l_eta1)<2.4)) || ((type==1)&&((abs(l_eta1)<2.5)&&!(abs(l_eta1)>1.4442 && abs(l_eta1)<1.566))))) ) continue;
      //if ( !(l_pt1>30 && ((type==1)&&((abs(l_eta1)<2.5)&&!(abs(l_eta1)>1.4442 && abs(l_eta1)<1.566))))) continue;

      if ( !(l_pt1>30 && (((type==0)&&(abs(l_eta1)<2.4)) ))) continue;

      bool isZ=false;

      if ( l_pt2>30 && (((type==0)&&(abs(l_eta2)<2.4)) || ((type==1)&&((abs(l_eta2)<2.5)&&!(abs(l_eta2)>1.4442 && abs(l_eta2)<1.566)))) ) isZ=true;
      //&& 
      //   dilep_m>125 && dilep_m<150 ) isZ=true;
      
      if (!isZ && !( l_pt2<0 && (((type==0)&&(pfMET_Corr>30)) || ((type==1)&&(pfMET_Corr>50))) ) ) continue;

      float weight=(samp1.at(xx).xsec*samp1.at(xx).xsecCorr*lumi*genWeight*trig_eff_Weight*id_eff_Weight*pu_Weight*btag0Wgt)/(1.0*(samp1.at(xx).nEvents-2*samp1.at(xx).nNegEvents));

      if (samp1.at(xx).sampname=="dataM" || samp1.at(xx).sampname=="dataE") weight=1.0;

      float deta=abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta);

      if (isResolved && !isZ) {
	//mVV_Wjj->Fill(mass_lvj_type0_PuppiAK8, weight);
	dETA_Wjj->Fill(deta, weight);
	mVBF_Wjj->Fill(vbf_maxpt_jj_m, weight);
	mVJ_Wjj->Fill(PuppiAK8_jet_mass_so_corr,weight);
	MET_Wjj->Fill(pfMET_Corr,weight);
	ptL_Wjj->Fill(l_pt1,weight);
      }
      else if (!isResolved && !isZ) {
	//mVV_WV->Fill(mass_lvj_type0_PuppiAK8, weight);
	dETA_WV->Fill(deta, weight);
	mVBF_WV->Fill(vbf_maxpt_jj_m, weight);
	mVJ_WV->Fill(PuppiAK8_jet_mass_so_corr,weight);
	MET_WV->Fill(pfMET_Corr,weight);
	ptL_WV->Fill(l_pt1,weight);
      }
      else if (isResolved) {
	//mVV_Zjj->Fill(mass_llj_PuppiAK8, weight);
	dETA_Zjj->Fill(deta, weight);
	mVBF_Zjj->Fill(vbf_maxpt_jj_m, weight);
	mVL_Zjj->Fill(dilep_m,weight);
	mVJ_Zjj->Fill(PuppiAK8_jet_mass_so_corr,weight);
      }
      else {
	//mVV_ZV->Fill(mass_llj_PuppiAK8, weight);
	dETA_ZV->Fill(deta, weight);
	mVBF_ZV->Fill(vbf_maxpt_jj_m, weight);
	mVL_ZV->Fill(dilep_m,weight);
	mVJ_ZV->Fill(PuppiAK8_jet_mass_so_corr,weight);

      }
    }
  }
  
  out->Write();
  out->Close();

  
}
