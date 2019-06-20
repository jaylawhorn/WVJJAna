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

Sample make_Sample(string sname, string fname, float xsec, float xsecCorr, 
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

void Draw2(vector<Sample> samp1, TString outfile);

void Draw() {

  string sampname;
  string filename;
  float  xsec;
  float  xsecCorr;
  int    nEvents;
  int    nNegEvents;
  int    fillColor;
  bool   doStack;

  vector<Sample> DataM;
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
      VBF_EWK.push_back(make_Sample(sampname, filename, xsec, xsecCorr, nEvents, nNegEvents));
    }
    else if (sampname == "VBF_QCD") {
      int index=VBF_QCD.size();
      VBF_QCD.push_back(make_Sample(sampname, filename, xsec, xsecCorr, nEvents, nNegEvents));
    }
    else if (sampname == "Top") {
      int index=Top.size();
      Top.push_back(make_Sample(sampname, filename, xsec, xsecCorr, nEvents, nNegEvents));
    }
    else if (sampname == "WJets") {
      int index=WJets.size();
      WJets.push_back(make_Sample(sampname, filename, xsec, xsecCorr, nEvents, nNegEvents));
    }
    else if (sampname == "DYJets") {
      int index=DYJets.size();
      DYJets.push_back(make_Sample(sampname, filename, xsec, xsecCorr, nEvents, nNegEvents));
    }
    else if (sampname == "dataM") {
      int index=DataM.size();
      DataM.push_back(make_Sample(sampname, filename, xsec, xsecCorr, nEvents, nNegEvents));
    }
  }

  cout << DataM.size() << endl;

  //Draw2(VBF_EWK,"VBF_EWK.root");
  //Draw2(VBF_QCD,"VBF_QCD.root");
  //Draw2(Top,"Top.root");
  //Draw2(WJets,"WJets.root");
  //Draw2(DYJets,"DYJets.root");
  Draw2(DataM,"DataM.root");

}

void Draw2(vector<Sample> samp1, TString outfile) {

  TFile *out = new TFile(outfile,"RECREATE");

  TH1D *hDetajj = new TH1D("hDetajj","hDetajj",50,0,10); hDetajj->Sumw2();
  TH1D *hMjj    = new TH1D("hMjj","hMjj",50,400,2400); hMjj->Sumw2();
  TH2D *hDetaVMjj = new TH2D("hDetaVMjj", "hDetaVMjj", 50, 0, 10, 50, 400, 2400); hDetaVMjj->Sumw2();
  
  int NBINS = 7;
  double MINRange = 150;
  double MAXRange = 2500;
  double massLEdges[8] = {150, 300, 450, 600, 1075, 1550, 2025, 2500};

  //TH1D *hWejj[9];
  //hWejj[0] = new TH1D("hWejj_00","hWejj_00",NBINS, massLEdges); hWejj[0]->Sumw2();
  //hWejj[1] = new TH1D("hWejj_10","hWejj_10",NBINS, massLEdges); hWejj[1]->Sumw2();
  //hWejj[2] = new TH1D("hWejj_20","hWejj_20",NBINS, massLEdges); hWejj[2]->Sumw2();
  //hWejj[3] = new TH1D("hWejj_01","hWejj_01",NBINS, massLEdges); hWejj[3]->Sumw2();
  //hWejj[4] = new TH1D("hWejj_11","hWejj_11",NBINS, massLEdges); hWejj[4]->Sumw2();
  //hWejj[5] = new TH1D("hWejj_21","hWejj_21",NBINS, massLEdges); hWejj[5]->Sumw2();
  //hWejj[6] = new TH1D("hWejj_02","hWejj_02",NBINS, massLEdges); hWejj[6]->Sumw2();
  //hWejj[7] = new TH1D("hWejj_12","hWejj_12",NBINS, massLEdges); hWejj[7]->Sumw2();
  //hWejj[8] = new TH1D("hWejj_22","hWejj_22",NBINS, massLEdges); hWejj[8]->Sumw2();

  TH1D *hWmjj[9];
  hWmjj[0] = new TH1D(Form("hWmjj_00_%s",samp1.at(0).sampname.c_str()),Form("hWmjj_00_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hWmjj[0]->Sumw2();
  hWmjj[1] = new TH1D(Form("hWmjj_10_%s",samp1.at(0).sampname.c_str()),Form("hWmjj_10_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hWmjj[1]->Sumw2();
  hWmjj[2] = new TH1D(Form("hWmjj_20_%s",samp1.at(0).sampname.c_str()),Form("hWmjj_20_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hWmjj[2]->Sumw2();
  hWmjj[3] = new TH1D(Form("hWmjj_01_%s",samp1.at(0).sampname.c_str()),Form("hWmjj_01_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hWmjj[3]->Sumw2();
  hWmjj[4] = new TH1D(Form("hWmjj_11_%s",samp1.at(0).sampname.c_str()),Form("hWmjj_11_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hWmjj[4]->Sumw2();
  hWmjj[5] = new TH1D(Form("hWmjj_21_%s",samp1.at(0).sampname.c_str()),Form("hWmjj_21_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hWmjj[5]->Sumw2();
  hWmjj[6] = new TH1D(Form("hWmjj_02_%s",samp1.at(0).sampname.c_str()),Form("hWmjj_02_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hWmjj[6]->Sumw2();
  hWmjj[7] = new TH1D(Form("hWmjj_12_%s",samp1.at(0).sampname.c_str()),Form("hWmjj_12_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hWmjj[7]->Sumw2();
  hWmjj[8] = new TH1D(Form("hWmjj_22_%s",samp1.at(0).sampname.c_str()),Form("hWmjj_22_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hWmjj[8]->Sumw2();

  //TH1D *hZejj[9];
  //hZejj[0] = new TH1D("hZejj_00","hZejj_00",NBINS, massLEdges); hZejj[0]->Sumw2();
  //hZejj[1] = new TH1D("hZejj_10","hZejj_10",NBINS, massLEdges); hZejj[1]->Sumw2();
  //hZejj[2] = new TH1D("hZejj_20","hZejj_20",NBINS, massLEdges); hZejj[2]->Sumw2();
  //hZejj[3] = new TH1D("hZejj_01","hZejj_01",NBINS, massLEdges); hZejj[3]->Sumw2();
  //hZejj[4] = new TH1D("hZejj_11","hZejj_11",NBINS, massLEdges); hZejj[4]->Sumw2();
  //hZejj[5] = new TH1D("hZejj_21","hZejj_21",NBINS, massLEdges); hZejj[5]->Sumw2();
  //hZejj[6] = new TH1D("hZejj_02","hZejj_02",NBINS, massLEdges); hZejj[6]->Sumw2();
  //hZejj[7] = new TH1D("hZejj_12","hZejj_12",NBINS, massLEdges); hZejj[7]->Sumw2();
  //hZejj[8] = new TH1D("hZejj_22","hZejj_22",NBINS, massLEdges); hZejj[8]->Sumw2();

  TH1D *hZmjj[9];
  hZmjj[0] = new TH1D(Form("hZmjj_00_%s",samp1.at(0).sampname.c_str()),Form("hZmjj_00_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hZmjj[0]->Sumw2();
  hZmjj[1] = new TH1D(Form("hZmjj_10_%s",samp1.at(0).sampname.c_str()),Form("hZmjj_10_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hZmjj[1]->Sumw2();
  hZmjj[2] = new TH1D(Form("hZmjj_20_%s",samp1.at(0).sampname.c_str()),Form("hZmjj_20_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hZmjj[2]->Sumw2();
  hZmjj[3] = new TH1D(Form("hZmjj_01_%s",samp1.at(0).sampname.c_str()),Form("hZmjj_01_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hZmjj[3]->Sumw2();
  hZmjj[4] = new TH1D(Form("hZmjj_11_%s",samp1.at(0).sampname.c_str()),Form("hZmjj_11_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hZmjj[4]->Sumw2();
  hZmjj[5] = new TH1D(Form("hZmjj_21_%s",samp1.at(0).sampname.c_str()),Form("hZmjj_21_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hZmjj[5]->Sumw2();
  hZmjj[6] = new TH1D(Form("hZmjj_02_%s",samp1.at(0).sampname.c_str()),Form("hZmjj_02_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hZmjj[6]->Sumw2();
  hZmjj[7] = new TH1D(Form("hZmjj_12_%s",samp1.at(0).sampname.c_str()),Form("hZmjj_12_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hZmjj[7]->Sumw2();
  hZmjj[8] = new TH1D(Form("hZmjj_22_%s",samp1.at(0).sampname.c_str()),Form("hZmjj_22_%s",samp1.at(0).sampname.c_str()),NBINS, massLEdges); hZmjj[8]->Sumw2();

  for (int i=0; i<9; i++) {
    hWmjj[i]->SetTitle(TString(samp1.at(0).sampname));
    hZmjj[i]->SetTitle(TString(samp1.at(0).sampname));

    hWmjj[i]->GetXaxis()->SetTitle("M_{lvjj}");
    hZmjj[i]->GetXaxis()->SetTitle("M_{lljj}");
  }

  Float_t lumi=35867.06;

  Bool_t  isResolved;
  Int_t   type;
  Float_t genWeight=1, pu_Weight=1, btag0Wgt=1, id_eff_Weight=1, trig_eff_Weight=1;
  Float_t l_pt1=0, l_eta1=0, l_phi1=0, l_e1=0;
  Float_t l_pt2=0, l_eta2=0, l_phi2=0, l_e2=0;
  Float_t l_pt1_Up=0, l_pt1_Down=0, l_pt2_Up=0, l_pt2_Down=0;
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
  Float_t mass_lvj_type0_PuppiAK8_jes_up=0, mass_lvj_type0_PuppiAK8_jes_dn=0;
  Float_t vbf_maxpt_jj_m_jes_up=0, vbf_maxpt_jj_m_jes_dn=0;
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
    
    //intree->SetBranchAddress("l_pt1_Up", &l_pt1_Up);
    //intree->SetBranchAddress("l_pt1_Down", &l_pt1_Down);
    //intree->SetBranchAddress("l_pt2_Up", &l_pt2_Up);
    //intree->SetBranchAddress("l_pt2_Down", &l_pt2_Down);
    
    intree->SetBranchAddress("dilep_pt", &dilep_pt);
    intree->SetBranchAddress("dilep_eta", &dilep_eta);
    intree->SetBranchAddress("dilep_phi", &dilep_phi);
    intree->SetBranchAddress("dilep_m", &dilep_m);
    
    intree->SetBranchAddress("pfMET_Corr", &pfMET_Corr);
    intree->SetBranchAddress("pfMET_Corr_phi", &pfMET_Corr_phi);
    intree->SetBranchAddress("nu_pz_type0", &nu_pz_type0);
    
    //intree->SetBranchAddress("AK4_vjet1_pt", &AK4_vjet1_pt);
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
    intree->SetBranchAddress("PuppiAK8_jet_mass_so_corr", &PuppiAK8_jet_mass_so_corr);
    
    intree->SetBranchAddress("vbf_maxpt_j1_pt", &vbf_maxpt_j1_pt);
    intree->SetBranchAddress("vbf_maxpt_j1_eta", &vbf_maxpt_j1_eta);
    intree->SetBranchAddress("vbf_maxpt_j1_phi", &vbf_maxpt_j1_phi);
    intree->SetBranchAddress("vbf_maxpt_j1_e", &vbf_maxpt_j1_e);
    intree->SetBranchAddress("vbf_maxpt_j1_mass", &vbf_maxpt_j1_mass);
    intree->SetBranchAddress("vbf_maxpt_j2_pt", &vbf_maxpt_j2_pt);
    intree->SetBranchAddress("vbf_maxpt_j2_eta", &vbf_maxpt_j2_eta);
    intree->SetBranchAddress("vbf_maxpt_j2_phi", &vbf_maxpt_j2_phi);
    intree->SetBranchAddress("vbf_maxpt_j2_e", &vbf_maxpt_j2_e);
    intree->SetBranchAddress("vbf_maxpt_j2_mass", &vbf_maxpt_j2_mass);
    
    intree->SetBranchAddress("nBTagJet_loose", &nBTagJet_loose);
    
    intree->SetBranchAddress("vbf_maxpt_jj_m", &vbf_maxpt_jj_m);
    intree->SetBranchAddress("vbf_maxpt_jj_m_jes_up", &vbf_maxpt_jj_m_jes_up);
    intree->SetBranchAddress("vbf_maxpt_jj_m_jes_dn", &vbf_maxpt_jj_m_jes_dn);
    intree->SetBranchAddress("mass_lvj_type0_PuppiAK8", &mass_lvj_type0_PuppiAK8);
    intree->SetBranchAddress("mass_lvj_type0_PuppiAK8_jes_up", &mass_lvj_type0_PuppiAK8_jes_up);
    intree->SetBranchAddress("mass_lvj_type0_PuppiAK8_jes_dn", &mass_lvj_type0_PuppiAK8_jes_dn);
    intree->SetBranchAddress("mass_llj_PuppiAK8", &mass_llj_PuppiAK8);    

    intree->SetBranchAddress("BosonCentrality_type0", &BosonCentrality_type0);
    intree->SetBranchAddress("ZeppenfeldWL_type0", &ZeppenfeldWL_type0);
    intree->SetBranchAddress("ZeppenfeldWH", &ZeppenfeldWH);
    
    
    for (int i=0; i<intree->GetEntries(); i++) {
      intree->GetEntry(i);
      
      if ( !(nBTagJet_loose && isResolved &&
	     ((PuppiAK8_jet_mass_so_corr>65) && (PuppiAK8_jet_mass_so_corr<105)) &&
	     ((vbf_maxpt_j1_pt>30) && (vbf_maxpt_j2_pt>30))) ) continue;
      
      //if ( !(l_pt1>50 && (((type==0)&&(abs(l_eta1)<2.4)) || ((type==1)&&((abs(l_eta1)<2.5)&&!(abs(l_eta1)>1.4442 && abs(l_eta1)<1.566))))) ) continue;
      if ( !(l_pt1>50 && abs(l_eta1)<2.4) ) continue;
      
      float weight=(samp1.at(xx).xsec*samp1.at(xx).xsecCorr*lumi*genWeight*trig_eff_Weight*id_eff_Weight*pu_Weight*btag0Wgt)/(1.0*(samp1.at(xx).nEvents-2*samp1.at(xx).nNegEvents));

      if (samp1.at(xx).sampname=="dataM") weight=1.0;

      hDetajj->Fill( abs(vbf_maxpt_j1_eta - vbf_maxpt_j2_eta), weight);
      hMjj->Fill( vbf_maxpt_jj_m, weight);
      
      hDetaVMjj->Fill(abs(vbf_maxpt_j1_eta - vbf_maxpt_j2_eta), vbf_maxpt_jj_m, weight);

      bool isW=true;
      if (l_pt2>30 && abs(l_eta2)<2.4 && dilep_m>76 && dilep_m<107) isW=false;
      if ( !((l_pt2<0 && pfMET_Corr>50)||isW==false) ) continue;

      if (abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta)>3.0 && abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta)<4.0) {
	if (vbf_maxpt_jj_m>500 && vbf_maxpt_jj_m<600) {
	  if (isW) hWmjj[0]->Fill(mass_lvj_type0_PuppiAK8, weight);
	  else hZmjj[0]->Fill(mass_llj_PuppiAK8, weight);
	}
	else if (vbf_maxpt_jj_m>600 && vbf_maxpt_jj_m<1000) {
	  if (isW) hWmjj[1]->Fill(mass_lvj_type0_PuppiAK8, weight);
	  else hZmjj[1]->Fill(mass_llj_PuppiAK8, weight);
	}
	else if (vbf_maxpt_jj_m>1000) {
	  if (isW) hWmjj[2]->Fill(mass_lvj_type0_PuppiAK8, weight);
	  else hZmjj[2]->Fill(mass_llj_PuppiAK8, weight);
	}
      }
      else if (abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta)>4.0 && abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta)<5.0) {
	if (vbf_maxpt_jj_m>500 && vbf_maxpt_jj_m<600) {
	  if (isW) hWmjj[3]->Fill(mass_lvj_type0_PuppiAK8, weight);
	  else hZmjj[3]->Fill(mass_llj_PuppiAK8, weight);
	}
	else if (vbf_maxpt_jj_m>600 && vbf_maxpt_jj_m<1000) {
	  if (isW) hWmjj[4]->Fill(mass_lvj_type0_PuppiAK8, weight);
	  else hZmjj[4]->Fill(mass_llj_PuppiAK8, weight);
	}
	else if (vbf_maxpt_jj_m>1000) {
	  if (isW) hWmjj[5]->Fill(mass_lvj_type0_PuppiAK8, weight);
	  else hZmjj[5]->Fill(mass_llj_PuppiAK8, weight);
	}
      }

      if (abs(vbf_maxpt_j2_eta-vbf_maxpt_j1_eta)>5.0) {
	if (vbf_maxpt_jj_m>500 && vbf_maxpt_jj_m<600) {
	  if (isW) hWmjj[6]->Fill(mass_lvj_type0_PuppiAK8, weight);
	  else hZmjj[6]->Fill(mass_llj_PuppiAK8, weight);
	}
	else if (vbf_maxpt_jj_m>600 && vbf_maxpt_jj_m<1000) {
	  if (isW) hWmjj[7]->Fill(mass_lvj_type0_PuppiAK8, weight);
	  else hZmjj[7]->Fill(mass_llj_PuppiAK8, weight);
	}
	else if (vbf_maxpt_jj_m>1000) {
	  if (isW) hWmjj[8]->Fill(mass_lvj_type0_PuppiAK8, weight);
	  else hZmjj[8]->Fill(mass_llj_PuppiAK8, weight);
	}
      }
    }
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  gStyle->SetOptStat(0);

  hDetajj->SetTitle(TString(samp1.at(0).sampname));
  hDetajj->GetXaxis()->SetTitle("dEta_{jj}");
  hDetajj->Draw("histe");
  c1->SaveAs(Form("detajj_%s.png",samp1.at(0).sampname.c_str()));

  hMjj->SetTitle(TString(samp1.at(0).sampname));
  hMjj->GetXaxis()->SetTitle("m_{jj}");  
  hMjj->Draw("histe");
  c1->SaveAs(Form("mjj_%s.png",samp1.at(0).sampname.c_str()));
  
  hDetaVMjj->SetTitle(TString(samp1.at(0).sampname));
  hDetaVMjj->GetXaxis()->SetTitle("dEta_{jj}");
  hDetaVMjj->GetYaxis()->SetTitle("m_{jj}");
  hDetaVMjj->Draw("colz");
  c1->SaveAs(Form("2D_deta_mjj_%s.png",samp1.at(0).sampname.c_str()));

  out->Write();
  out->Close();

}

