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
  int    fillColor;
  bool   doStack;
};

Sample makeSample(string sname, string fname, float xsec, float xsecCorr) {

  Sample ret;
  ret.sampname=sname;
  ret.filename=fname;
  ret.xsec=xsec;
  ret.xsecCorr=xsecCorr;

  return ret;
};

struct Binning {
  int nBinsMVBF;
  int nBinsDETA;
  int nBinsMVV;
  int nBinsTot;

  float MVBF_LE[40]; //!! careful with hard coding here
  float DETA_LE[40];
  float MVV_LE[40];

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
  ifs.open("files2018.txt");
  assert(ifs.is_open());
  string line;

  while(getline(ifs,line)) {

    if(line[0]=='#') continue;
    stringstream ss(line);

    ss >> sampname >> filename >> xsec >> xsecCorr;
    
    if (sampname == "VBF_EWK") {
      int index=VBF_EWK.size();
      VBF_EWK.push_back(makeSample(sampname, filename, xsec, xsecCorr));
    }
    else if (sampname == "VBF_QCD") {
      int index=VBF_QCD.size();
      VBF_QCD.push_back(makeSample(sampname, filename, xsec, xsecCorr));
    }
    else if (sampname == "Top") {
      int index=Top.size();
      Top.push_back(makeSample(sampname, filename, xsec, xsecCorr));
    }
    else if (sampname == "WJets") {
      int index=WJets.size();
      WJets.push_back(makeSample(sampname, filename, xsec, xsecCorr));
    }
    else if (sampname == "DYJets") {
      int index=DYJets.size();
      DYJets.push_back(makeSample(sampname, filename, xsec, xsecCorr));
    }
    else if (sampname == "dataE") {
      int index=DataE.size();
      DataE.push_back(makeSample(sampname, filename, xsec, xsecCorr));
    }
    else if (sampname == "dataM") {
      int index=DataM.size();
      DataM.push_back(makeSample(sampname, filename, xsec, xsecCorr));
    }
  }

  const int nBinsMVBF=25;
  //const int nBinsMVBF=10;
  const int nBinsDETA=15;
  const int nBinsMVV=30;
  float MVBF_LE[nBinsMVBF+1] = { 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 
				 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 
				 2500, 2600, 2700, 2800, 2900, 3000 };

//  float MVBF_LE[nBinsMVBF+1] = { 500,700,900,1100,1300, 
//				 1500,1700,1900,2100,2300,
//				 2500};

  float DETA_LE[nBinsDETA+1] = { 2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
				 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0 };

  float MVV_LE[nBinsMVV+1] ={ 0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 
			      1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900,
			      2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000 };
  
  Binning bins=makeBinning(nBinsMVBF, nBinsDETA, nBinsMVV, 
			   MVBF_LE, DETA_LE, MVV_LE);

  Draw2(VBF_EWK,bins,"VBF_EWK_2018_vjets_el.root","nom");
  Draw2(VBF_QCD,bins,"VBF_QCD_2018_vjets_el.root","nom");
  Draw2(Top,bins,"Top_2018_vjets_el.root","nom");
  Draw2(WJets,bins,"WJets_2018_vjets_el.root","nom");
  Draw2(DYJets,bins,"DYJets_2018_vjets_el.root","nom");
  Draw2(DataM,bins,"DataM_2018_vjets_el.root","nom");

  //Draw2(WJets,bins,"w_mlm_jetbin.root","nom");
  //Draw2(DYJets,bins,"w_htbinned.root","nom");
  //Draw2(Top,bins,"w_amc_jetbin.root","nom");


}
void Draw2(vector<Sample> samp1, Binning bins, TString outfile, TString var) {

  TFile *out;
  out = new TFile(outfile,"RECREATE");

  TString histName;

  histName = Form("%s_nJet30_Wjj",samp1.at(0).sampname.c_str());
  TH1D* nJet30_Wjj = new TH1D(histName, histName, 8, 2, 10);
  nJet30_Wjj->Sumw2();
  nJet30_Wjj->SetTitle(TString(samp1.at(0).sampname));
  nJet30_Wjj->GetXaxis()->SetTitle("nJet(pT>30), Wjj");

  histName = Form("%s_nJet50_Wjj",samp1.at(0).sampname.c_str());
  TH1D* nJet50_Wjj = new TH1D(histName, histName, 8, 2, 10);
  nJet50_Wjj->Sumw2();
  nJet50_Wjj->SetTitle(TString(samp1.at(0).sampname));
  nJet50_Wjj->GetXaxis()->SetTitle("nJet(pT>50), Wjj");

  histName = Form("%s_nJet30_WV",samp1.at(0).sampname.c_str());
  TH1D* nJet30_WV = new TH1D(histName, histName, 8, 2, 10);
  nJet30_WV->Sumw2();
  nJet30_WV->SetTitle(TString(samp1.at(0).sampname));
  nJet30_WV->GetXaxis()->SetTitle("nJet(pT>30), WV");

  histName = Form("%s_nJet50_WV",samp1.at(0).sampname.c_str());
  TH1D* nJet50_WV = new TH1D(histName, histName, 8, 2, 10);
  nJet50_WV->Sumw2();
  nJet50_WV->SetTitle(TString(samp1.at(0).sampname));
  nJet50_WV->GetXaxis()->SetTitle("nJet(pT>50), WV");

  histName = Form("%s_nJet30_Zjj",samp1.at(0).sampname.c_str());
  TH1D* nJet30_Zjj = new TH1D(histName, histName, 8, 2, 10);
  nJet30_Zjj->Sumw2();
  nJet30_Zjj->SetTitle(TString(samp1.at(0).sampname));
  nJet30_Zjj->GetXaxis()->SetTitle("nJet(pT>30), Zjj");

  histName = Form("%s_nJet50_Zjj",samp1.at(0).sampname.c_str());
  TH1D* nJet50_Zjj = new TH1D(histName, histName, 8, 2, 10);
  nJet50_Zjj->Sumw2();
  nJet50_Zjj->SetTitle(TString(samp1.at(0).sampname));
  nJet50_Zjj->GetXaxis()->SetTitle("nJet(pT>50), Zjj");

  histName = Form("%s_nJet30_ZV",samp1.at(0).sampname.c_str());
  TH1D* nJet30_ZV = new TH1D(histName, histName, 8, 2, 10);
  nJet30_ZV->Sumw2();
  nJet30_ZV->SetTitle(TString(samp1.at(0).sampname));
  nJet30_ZV->GetXaxis()->SetTitle("nJet(pT>30), ZV");

  histName = Form("%s_nJet50_ZV",samp1.at(0).sampname.c_str());
  TH1D* nJet50_ZV = new TH1D(histName, histName, 8, 2, 10);
  nJet50_ZV->Sumw2();
  nJet50_ZV->SetTitle(TString(samp1.at(0).sampname));
  nJet50_ZV->GetXaxis()->SetTitle("nJet(pT>50), ZV");

  histName = Form("%s_ETA_lep_Wjj",samp1.at(0).sampname.c_str());
  TH1D* ETA_lep_Wjj = new TH1D(histName, histName, 48, -3.5, 3.5);
  ETA_lep_Wjj->Sumw2();
  ETA_lep_Wjj->SetTitle(TString(samp1.at(0).sampname));
  ETA_lep_Wjj->GetXaxis()->SetTitle("Eta(lep), Wjj");

  histName = Form("%s_ETA_bos_Wjj",samp1.at(0).sampname.c_str());
  TH1D* ETA_bos_Wjj = new TH1D(histName, histName, 20, -5, 5);
  ETA_bos_Wjj->Sumw2();
  ETA_bos_Wjj->SetTitle(TString(samp1.at(0).sampname));
  ETA_bos_Wjj->GetXaxis()->SetTitle("Eta(jj), Wjj");

  histName = Form("%s_ETA_lep_WV",samp1.at(0).sampname.c_str());
  TH1D* ETA_lep_WV = new TH1D(histName, histName, 48, -3.5, 3.5);
  ETA_lep_WV->Sumw2();
  ETA_lep_WV->SetTitle(TString(samp1.at(0).sampname));
  ETA_lep_WV->GetXaxis()->SetTitle("Eta(lep), WV");

  histName = Form("%s_ETA_bos_WV",samp1.at(0).sampname.c_str());
  TH1D* ETA_bos_WV = new TH1D(histName, histName, 16, -4, 4);
  ETA_bos_WV->Sumw2();
  ETA_bos_WV->SetTitle(TString(samp1.at(0).sampname));
  ETA_bos_WV->GetXaxis()->SetTitle("Eta(V), WV");

  histName = Form("%s_ETA_lep_ZV",samp1.at(0).sampname.c_str());
  TH1D* ETA_lep_ZV = new TH1D(histName, histName, 12, -3.5, 3.5);
  ETA_lep_ZV->Sumw2();
  ETA_lep_ZV->SetTitle(TString(samp1.at(0).sampname));
  ETA_lep_ZV->GetXaxis()->SetTitle("Eta(lep), ZV");

  histName = Form("%s_ETA_lep_Zjj",samp1.at(0).sampname.c_str());
  TH1D* ETA_lep_Zjj = new TH1D(histName, histName, 12, -3.5, 3.5);
  ETA_lep_Zjj->Sumw2();
  ETA_lep_Zjj->SetTitle(TString(samp1.at(0).sampname));
  ETA_lep_Zjj->GetXaxis()->SetTitle("Eta(lep), Zjj");

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

  histName = Form("%s_ETA1_Wjj",samp1.at(0).sampname.c_str());
  TH1D* ETA1_Wjj = new TH1D(histName, histName, 40, -5, 5);
  ETA1_Wjj->Sumw2();
  ETA1_Wjj->SetTitle(TString(samp1.at(0).sampname));
  ETA1_Wjj->GetXaxis()->SetTitle("Eta(j1), Wjj");

  histName = Form("%s_ETA2_Wjj",samp1.at(0).sampname.c_str());
  TH1D* ETA2_Wjj = new TH1D(histName, histName, 40, -5, 5);
  ETA2_Wjj->Sumw2();
  ETA2_Wjj->SetTitle(TString(samp1.at(0).sampname));
  ETA2_Wjj->GetXaxis()->SetTitle("Eta(j2), Wjj");

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
  TH1D* MET_Wjj = new TH1D(histName, histName, 40, 0, 1000);
  MET_Wjj->Sumw2();
  MET_Wjj->SetTitle(TString(samp1.at(0).sampname));
  MET_Wjj->GetXaxis()->SetTitle("MET, Wjj");

  histName = Form("%s_ptL_Wjj",samp1.at(0).sampname.c_str());
  TH1D* ptL_Wjj = new TH1D(histName, histName, 40, 0, 200);
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

  histName = Form("%s_ETA1_WV",samp1.at(0).sampname.c_str());
  TH1D* ETA1_WV = new TH1D(histName, histName, 40, -5, 5);
  ETA1_WV->Sumw2();
  ETA1_WV->SetTitle(TString(samp1.at(0).sampname));
  ETA1_WV->GetXaxis()->SetTitle("Eta(j1), WV");

  histName = Form("%s_ETA2_WV",samp1.at(0).sampname.c_str());
  TH1D* ETA2_WV = new TH1D(histName, histName, 40, -5, 5);
  ETA2_WV->Sumw2();
  ETA2_WV->SetTitle(TString(samp1.at(0).sampname));
  ETA2_WV->GetXaxis()->SetTitle("Eta(j2), WV");

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
  TH1D* MET_WV = new TH1D(histName, histName, 40, 0, 1000);
  MET_WV->Sumw2();
  MET_WV->SetTitle(TString(samp1.at(0).sampname));
  MET_WV->GetXaxis()->SetTitle("MET, WV");

  histName = Form("%s_ptL_WV",samp1.at(0).sampname.c_str());
  TH1D* ptL_WV = new TH1D(histName, histName, 40, 0, 200);
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

  histName = Form("%s_ETA1_Zjj",samp1.at(0).sampname.c_str());
  TH1D* ETA1_Zjj = new TH1D(histName, histName, 20, -5, 5);
  ETA1_Zjj->Sumw2();
  ETA1_Zjj->SetTitle(TString(samp1.at(0).sampname));
  ETA1_Zjj->GetXaxis()->SetTitle("Eta(j1), Zjj");

  histName = Form("%s_ETA2_Zjj",samp1.at(0).sampname.c_str());
  TH1D* ETA2_Zjj = new TH1D(histName, histName, 20, -5, 5);
  ETA2_Zjj->Sumw2();
  ETA2_Zjj->SetTitle(TString(samp1.at(0).sampname));
  ETA2_Zjj->GetXaxis()->SetTitle("Eta(j2), Zjj");

  histName = Form("%s_mVBF_Zjj",samp1.at(0).sampname.c_str());
  TH1D* mVBF_Zjj = new TH1D(histName, histName, bins.nBinsMVBF, &bins.MVBF_LE[0]);
  mVBF_Zjj->Sumw2();
  mVBF_Zjj->SetTitle(TString(samp1.at(0).sampname));
  mVBF_Zjj->GetXaxis()->SetTitle("m(VBF), Zjj");

  histName = Form("%s_mVJ_Zjj",samp1.at(0).sampname.c_str());
  TH1D* mVJ_Zjj = new TH1D(histName, histName, 12, 30, 150);
  mVJ_Zjj->Sumw2();
  mVJ_Zjj->SetTitle(TString(samp1.at(0).sampname));
  mVJ_Zjj->GetXaxis()->SetTitle("m(Z) had, Zjj");

  histName = Form("%s_mVL_Zjj",samp1.at(0).sampname.c_str());
  TH1D* mVL_Zjj = new TH1D(histName, histName, 12, 60, 120);
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

  histName = Form("%s_ETA1_ZV",samp1.at(0).sampname.c_str());
  TH1D* ETA1_ZV = new TH1D(histName, histName, 20, -5, 5);
  ETA1_ZV->Sumw2();
  ETA1_ZV->SetTitle(TString(samp1.at(0).sampname));
  ETA1_ZV->GetXaxis()->SetTitle("Eta(j1), ZV");

  histName = Form("%s_ETA2_ZV",samp1.at(0).sampname.c_str());
  TH1D* ETA2_ZV = new TH1D(histName, histName, 20, -5, 5);
  ETA2_ZV->Sumw2();
  ETA2_ZV->SetTitle(TString(samp1.at(0).sampname));
  ETA2_ZV->GetXaxis()->SetTitle("Eta(j2), ZV");

  histName = Form("%s_mVBF_ZV",samp1.at(0).sampname.c_str());
  TH1D* mVBF_ZV = new TH1D(histName, histName, bins.nBinsMVBF, &bins.MVBF_LE[0]);
  mVBF_ZV->Sumw2();
  mVBF_ZV->SetTitle(TString(samp1.at(0).sampname));
  mVBF_ZV->GetXaxis()->SetTitle("m(VBF), ZV");

  histName = Form("%s_mVJ_ZV",samp1.at(0).sampname.c_str());
  TH1D* mVJ_ZV = new TH1D(histName, histName, 12, 30, 150);
  mVJ_ZV->Sumw2();
  mVJ_ZV->SetTitle(TString(samp1.at(0).sampname));
  mVJ_ZV->GetXaxis()->SetTitle("m(Z) had, ZV");

  histName = Form("%s_mVL_ZV",samp1.at(0).sampname.c_str());
  TH1D* mVL_ZV = new TH1D(histName, histName, 12, 60, 120);
  mVL_ZV->Sumw2();
  mVL_ZV->SetTitle(TString(samp1.at(0).sampname));
  mVL_ZV->GetXaxis()->SetTitle("m(Z) lep, ZV");

  //Float_t lumi=35867.06;
  //Float_t lumi=41530.0;
  Float_t lumi=59000.0;
  //Float_t lumi=1000.0;
  
  const float MUON_MASS = 0.1056583745;
  const float ELE_MASS  = 0.000511;

  //Float_t genWeight=1, pu_Weight=1, btag0Wgt=1, id_eff_Weight=1, trig_eff_Weight=1;
  Float_t btagWeight=1;
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
  bool trigger_1Mu=false;
  bool trigger_2Mu=false;
  bool trigger_1El=false;
  bool trigger_2El=false;

  for (uint xx=0; xx<samp1.size(); xx++) {
    cout << samp1.at(xx).filename << endl;

    TFile* infile = new TFile(TString(samp1.at(xx).filename), "READ");   assert(infile);
    TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

    double nTotal=1, nNeg=1;
    if (!(samp1.at(xx).sampname=="dataM" || samp1.at(xx).sampname=="dataE")) {

      TH1F* hTotEvents = (TH1F*) infile->Get("TotalEvents"); assert(hTotEvents);
      nTotal = hTotEvents->GetBinContent(2);
      nNeg = hTotEvents->GetBinContent(1);
    }
    
    intree->SetBranchAddress("genWeight",&genWeight);
    intree->SetBranchAddress("puWeight",&puWeight);
    intree->SetBranchAddress("L1PFWeight",&L1PFWeight);
    intree->SetBranchAddress("btagWeight",&btagWeight);
    intree->SetBranchAddress("lep1_idEffWeight",&lep1_idEffWeight);
    intree->SetBranchAddress("lep2_idEffWeight",&lep2_idEffWeight);

    intree->SetBranchAddress("trigger_1Mu",&trigger_1Mu);
    intree->SetBranchAddress("trigger_2Mu",&trigger_2Mu);
    intree->SetBranchAddress("trigger_1El",&trigger_1El);
    intree->SetBranchAddress("trigger_2El",&trigger_2El);

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
    //intree->SetBranchAddress("MET_2017raw", &MET_2017raw);
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
    //intree->SetBranchAddress("vbf_deta", &vbf_deta);

    intree->SetBranchAddress("dibos_m", &dibos_m);

    intree->SetBranchAddress("nJet30", &nJet30);
    intree->SetBranchAddress("nJet50", &nJet50);

    intree->SetBranchAddress("nBtag_loose", &nBtag_loose);

    intree->SetBranchAddress("bosCent", &bosCent);
    intree->SetBranchAddress("zeppLep", &zeppLep);
    intree->SetBranchAddress("zeppHad", &zeppHad);

    cout << intree->GetEntries() << endl;

    for (int i=0; i<intree->GetEntries(); i++) {
    //for (int i=0; i<10000; i++) {
      intree->GetEntry(i);

      bool isEle=false, isResolved=false, isZ=false;

      if (bos_PuppiAK8_m_sd0_corr > 0 && bos_AK4AK4_m < 0) { isResolved=false; }
      else if (bos_PuppiAK8_m_sd0_corr < 0 && bos_AK4AK4_m > 0) { isResolved=true; }
      else {
	//cout << "both or neither of resolved and boosted mass is defined" << endl;
	continue;
      }
      
      if (lep1_m == ELE_MASS) { isEle=true; }
      else if (lep1_m == MUON_MASS) { isEle=false; }
      else {
	cout << "lepton is not electron or muon! skipping" << endl;
	continue;
      }

      //if (fabs(vbf1_AK4_eta)>2.65 && fabs(vbf1_AK4_eta)<3.139) continue;
      //if (fabs(vbf2_AK4_eta)>2.65 && fabs(vbf2_AK4_eta)<3.139) continue;

      if ( vbf_m < 500) continue;
      if ( fabs(vbf1_AK4_eta - vbf2_AK4_eta)<2.5) continue;
      //if ( dibos_m < 500) continue;

      //loose
      //if ( !(nBtag_loose==0 && vbf1_AK4_pt>50 && vbf2_AK4_pt>50) ) continue;
      //if (isResolved==true && (bos_AK4AK4_m<65 ||bos_AK4AK4_m>105)) continue;      
      //if (isResolved==false && (bos_PuppiAK8_m_sd0_corr<65 ||bos_PuppiAK8_m_sd0_corr>105)) continue;

      //ttbar
      //if ( !(nBtag_loose>0 && vbf1_AK4_pt>50 && vbf2_AK4_pt>50) ) continue;
      //if (isResolved==true && (bos_AK4AK4_m<65 ||bos_AK4AK4_m>105)) continue;      
      //if (isResolved==false && (bos_PuppiAK8_m_sd0_corr<65 ||bos_PuppiAK8_m_sd0_corr>105)) continue;

      //wjets
      if ( !(nBtag_loose==0 && vbf1_AK4_pt>50 && vbf2_AK4_pt>50) ) continue;
      if (isResolved==true && (bos_AK4AK4_m>65 &&bos_AK4AK4_m<105)) continue;
      if (isResolved==false && (bos_PuppiAK8_m_sd0_corr>65 &&bos_PuppiAK8_m_sd0_corr<105)) continue;

      //if (isEle==true || !(trigger_2Mu || trigger_1Mu)) continue;
      if (isEle==false || !(trigger_2El || trigger_1El)) continue;

      //remove HEM15/16
      //if (isEle==true && lep1_eta<-1.3 && lep1_phi>-1.57 && lep1_phi<-0.87) continue;

      //if (isResolved==true && (bos_j1_AK4_pt<30 || bos_j2_AK4_pt<30)) continue;
      if (isResolved==true && (bos_j1_AK4_pt<50 || bos_j2_AK4_pt<50)) continue;

      if (isResolved==false && bos_PuppiAK8_pt<200) continue;
      if (isResolved==false && bos_PuppiAK8_tau2tau1>0.55) continue;

      if (lep2_pt>0) isZ=true;
      if (isEle==true && (lep1_pt<35 || abs(lep1_eta)>2.5 || (abs(lep1_eta)>1.4442 && abs(lep1_eta)<1.566))) continue;
      if (isEle==false && (lep1_pt<35 || abs(lep1_eta)>2.4)) continue;

      if (isZ==true && (dilep_m < 81 || dilep_m > 101)) continue;
      if (isZ==true && isEle==true && (lep2_pt<20 || abs(lep2_eta)>2.5 || (abs(lep2_eta)>1.4442 && abs(lep2_eta)<1.566))) continue;
      if (isZ==true && isEle==false && (lep2_pt<20 || abs(lep2_eta)>2.4)) continue;
      if (isZ==true && (lep1_q*lep2_q)==1) continue;

      if (isZ==false && MET<30) continue;
      //if (isZ==false && bosCent<0) continue;

      //float weight=(samp1.at(xx).xsec*samp1.at(xx).xsecCorr*lumi*genWeight*puWeight*L1PFWeight)/(1.0*(nTotal-2*nNeg));
      double weight=(samp1.at(xx).xsec*lumi*genWeight)/(1.0*(nTotal));

      //if (isResolved==false && lep1_pt > 59.9 && lep1_pt<65.1) std::cout << weight << std::endl;

      if (samp1.at(xx).sampname=="dataM" || samp1.at(xx).sampname=="dataE") weight=1.0;
      else if (weight>100) {
	std::cout << "skipping large weight : " << weight << std::endl;
	continue;
      }

      if (isResolved==true && isZ==false) {

	nJet30_Wjj->Fill(nJet30,weight);
	nJet50_Wjj->Fill(nJet50,weight);

	ETA_bos_Wjj->Fill(bos_AK4AK4_eta,weight);
	ETA_lep_Wjj->Fill(lep1_eta,weight);

	mVV_Wjj->Fill(dibos_m, weight);
	dETA_Wjj->Fill(abs(vbf1_AK4_eta - vbf2_AK4_eta), weight);
	ETA1_Wjj->Fill(vbf1_AK4_eta, weight);
	ETA2_Wjj->Fill(vbf2_AK4_eta, weight);
	mVBF_Wjj->Fill(vbf_m, weight);
	mVJ_Wjj->Fill(bos_AK4AK4_m, weight);
	MET_Wjj->Fill(MET, weight);
	ptL_Wjj->Fill(lep1_pt, weight);
      } 
      else if (isResolved==false && isZ==false) {

	nJet30_WV->Fill(nJet30,weight);
	nJet50_WV->Fill(nJet50,weight);

	ETA_bos_WV->Fill(bos_PuppiAK8_eta,weight);
	ETA_lep_WV->Fill(lep1_eta,weight);

	mVV_WV->Fill(dibos_m, weight);
	dETA_WV->Fill(abs(vbf1_AK4_eta - vbf2_AK4_eta), weight);
	ETA1_WV->Fill(vbf1_AK4_eta, weight);
	ETA2_WV->Fill(vbf2_AK4_eta, weight);
	mVBF_WV->Fill(vbf_m, weight);
	mVJ_WV->Fill(bos_PuppiAK8_m_sd0_corr, weight);
	MET_WV->Fill(MET, weight);
	ptL_WV->Fill(lep1_pt,weight);
      }
      else if (isResolved==true && isZ==true) {

	nJet30_Zjj->Fill(nJet30,weight);
	nJet50_Zjj->Fill(nJet50,weight);

	ETA_lep_Zjj->Fill(lep1_eta,weight);
	ETA_lep_Zjj->Fill(lep2_eta,weight);

	mVV_Zjj->Fill(dibos_m, weight);
	dETA_Zjj->Fill(abs(vbf1_AK4_eta - vbf2_AK4_eta), weight);
	ETA1_Zjj->Fill(vbf1_AK4_eta, weight);
	ETA2_Zjj->Fill(vbf2_AK4_eta, weight);
	mVBF_Zjj->Fill(vbf_m, weight);
	mVL_Zjj->Fill(dilep_m,weight);
	mVJ_Zjj->Fill(bos_AK4AK4_m,weight);
      }
      else if (isResolved==false && isZ==true) {

	nJet30_ZV->Fill(nJet30,weight);
	nJet50_ZV->Fill(nJet50,weight);

	ETA_lep_ZV->Fill(lep1_eta,weight);
	ETA_lep_ZV->Fill(lep2_eta,weight);

	mVV_ZV->Fill(dibos_m, weight);
	dETA_ZV->Fill(abs(vbf1_AK4_eta - vbf2_AK4_eta), weight);
	ETA1_ZV->Fill(vbf1_AK4_eta, weight);
	ETA2_ZV->Fill(vbf2_AK4_eta, weight);
	mVBF_ZV->Fill(vbf_m, weight);
	mVL_ZV->Fill(dilep_m, weight);
	mVJ_ZV->Fill(bos_PuppiAK8_m_sd0_corr, weight);
      }
    }
  }
  
  out->Write();
  out->Close();

  
}
