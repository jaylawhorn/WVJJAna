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

void MakeCards() {

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
  ifs.open("files2017.txt");
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
  const int nBinsDETA=3; 
  const int nBinsMVV=4;
  float MVBF_LE[nBinsMVBF+1] = { 500, 600, 800, 1000, 3000 };
  //float MVBF_LE[nBinsMVBF+1] = { 600, 800, 1200, 3000 };
  // float DETA_LE[nBinsDETA+1] = { 3.0, 4.0, 5.0, 6.0, 10.0 };
  float DETA_LE[nBinsDETA+1] = { 4.0, 5.0, 6.0, 10.0 };
  //float MVV_LE[nBinsMVV+1] = {150, 300, 450, 600, 1075, 1550, 2025, 2500};
  float MVV_LE[nBinsMVV+1] = {600, 1075, 1550, 2025, 2500};
  
  Binning bins=makeBinning(nBinsMVBF, nBinsDETA, nBinsMVV, 
			   MVBF_LE, DETA_LE, MVV_LE);

  Draw2(VBF_EWK,bins,"VBF_EWK_2017_makeCards_Dec5.root","nom");
  //Draw2(VBF_EWK,bins,"VBF_EWK_2017_makeCards_Dec5.root","jesup");
  //Draw2(VBF_EWK,bins,"VBF_EWK_2017_makeCards_Dec5.root","jesdn");
  Draw2(VBF_QCD,bins,"VBF_QCD_2017_makeCards_Dec5.root","nom");
  //Draw2(VBF_QCD,bins,"VBF_QCD_2017_makeCards_Dec5.root","jesup");
  //Draw2(VBF_QCD,bins,"VBF_QCD_2017_makeCards_Dec5.root","jesdn");
  Draw2(Top,bins,"Top_2017_makeCards_Dec5.root","nom");
  //Draw2(Top,bins,"Top_2017_makeCards_Dec5.root","jesup");
  //Draw2(Top,bins,"Top_2017_makeCards_Dec5.root","jesdn");
  Draw2(WJets,bins,"WJets_2017_makeCards_Dec5.root","nom");
  //Draw2(WJets,bins,"WJets_2017_makeCards_Dec5.root","jesup");
  //Draw2(WJets,bins,"WJets_2017_makeCards_Dec5.root","jesdn");
  Draw2(DYJets,bins,"DYJets_2017_makeCards_Dec5.root","nom");
  //Draw2(DYJets,bins,"DYJets_2017_makeCards_Dec5.root","jesup");
  //Draw2(DYJets,bins,"DYJets_2017_makeCards_Dec5.root","jesdn");
  Draw2(DataM,bins,"DataM_2017_makeCards_Dec5.root","nom");
  //Draw2(DataE,bins,"DataE_2017_makeCards_Dec5.root","nom");


}
void Draw2(vector<Sample> samp1, Binning bins, TString outfile, TString var) {

  TFile *out;
  if (var=="nom") out = new TFile(outfile,"RECREATE");
  else out = new TFile(outfile,"UPDATE");

  TString histName;

  if (var=="nom") histName = Form("%s_mWjj",samp1.at(0).sampname.c_str());
  else if (var=="jesup") histName = Form("%s_mWjj_jesUp",samp1.at(0).sampname.c_str());
  else if (var=="jesdn") histName = Form("%s_mWjj_jesDown",samp1.at(0).sampname.c_str());
  TH1D* mWjj = new TH1D(histName, histName, bins.nBinsMVBF*bins.nBinsDETA*bins.nBinsMVV, 0, bins.nBinsMVBF*bins.nBinsDETA*bins.nBinsMVV);
  mWjj->Sumw2();
  mWjj->SetTitle(TString(samp1.at(0).sampname));

  if (var=="nom") histName = Form("%s_mWV",samp1.at(0).sampname.c_str());
  else if (var=="jesup") histName = Form("%s_mWV_jesUp",samp1.at(0).sampname.c_str());
  else if (var=="jesdn") histName = Form("%s_mWV_jesDown",samp1.at(0).sampname.c_str());
  TH1D* mWV = new TH1D(histName, histName, bins.nBinsMVBF*bins.nBinsDETA*bins.nBinsMVV, 0, bins.nBinsMVBF*bins.nBinsDETA*bins.nBinsMVV);
  mWV->Sumw2();
  mWV->SetTitle(TString(samp1.at(0).sampname));

  if (var=="nom") histName = Form("%s_mZjj",samp1.at(0).sampname.c_str());
  else if (var=="jesup") histName = Form("%s_mZjj_jesUp",samp1.at(0).sampname.c_str());
  else if (var=="jesdn") histName = Form("%s_mZjj_jesDown",samp1.at(0).sampname.c_str());
  TH1D* mZjj = new TH1D(histName, histName, bins.nBinsMVBF*bins.nBinsDETA*bins.nBinsMVV, 0, bins.nBinsMVBF*bins.nBinsDETA*bins.nBinsMVV);
  mZjj->Sumw2();
  mZjj->SetTitle(TString(samp1.at(0).sampname));

  if (var=="nom") histName = Form("%s_mZV",samp1.at(0).sampname.c_str());
  else if (var=="jesup") histName = Form("%s_mZV_jesUp",samp1.at(0).sampname.c_str());
  else if (var=="jesdn") histName = Form("%s_mZV_jesDown",samp1.at(0).sampname.c_str());
  TH1D* mZV = new TH1D(histName, histName, bins.nBinsMVBF*bins.nBinsDETA*bins.nBinsMVV, 0, bins.nBinsMVBF*bins.nBinsDETA*bins.nBinsMVV);
  mZV->Sumw2();
  mZV->SetTitle(TString(samp1.at(0).sampname));

  Float_t lumi=35867.06;
  //Float_t lumi=41530.0;

  const float MUON_MASS = 0.1056583745;
  const float ELE_MASS  = 0.000511;

  Float_t genWeight=1, puWeight=1, lep1_idEffWeight=1, lep2_idEffWeight=1;
  Float_t L1PFWeight=1;
  Float_t lep1_pt=0, lep1_eta=0, lep1_phi=0, lep1_m=0;
  Float_t lep2_pt=0, lep2_eta=0, lep2_phi=0, lep2_m=0;
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
  Float_t bosCent=0, zeppLep=0, zeppHad=0;

  for (uint xx=0; xx<samp1.size(); xx++) {
    cout << samp1.at(xx).filename << endl;

    TFile* infile = new TFile(TString(samp1.at(xx).filename), "READ");   assert(infile);
    TTree* intree = (TTree*) infile->Get("Events"); assert(intree);

    int nTotal=1, nNeg=1;
    if (!(samp1.at(xx).sampname=="dataM" || samp1.at(xx).sampname=="dataE")) {

      TH1F* hTotEvents = (TH1F*) infile->Get("TotalEvents"); assert(hTotEvents);
      nTotal = hTotEvents->GetBinContent(2);
      nNeg = hTotEvents->GetBinContent(1);
    }

    intree->SetBranchAddress("genWeight",&genWeight);
    intree->SetBranchAddress("puWeight",&puWeight);
    intree->SetBranchAddress("L1PFWeight",&L1PFWeight);
    intree->SetBranchAddress("lep1_idEffWeight",&lep1_idEffWeight);
    intree->SetBranchAddress("lep2_idEffWeight",&lep2_idEffWeight);

    intree->SetBranchAddress("lep1_pt", &lep1_pt);
    intree->SetBranchAddress("lep1_eta", &lep1_eta);
    intree->SetBranchAddress("lep1_phi", &lep1_phi);
    intree->SetBranchAddress("lep1_m", &lep1_m);
    intree->SetBranchAddress("lep2_pt", &lep2_pt);
    intree->SetBranchAddress("lep2_eta", &lep2_eta);
    intree->SetBranchAddress("lep2_phi", &lep2_phi);
    intree->SetBranchAddress("lep2_m", &lep2_m);

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
    //intree->SetBranchAddress("bos_PuppiAK8_m_sd0_corr", &bos_PuppiAK8_m_sd0_corr);
    intree->SetBranchAddress("bos_PuppiAK8_tau2tau1", &bos_PuppiAK8_tau2tau1);

    //intree->SetBranchAddress("vbf1_AK4_pt", &vbf1_AK4_pt);
    intree->SetBranchAddress("vbf1_AK4_eta", &vbf1_AK4_eta);
    intree->SetBranchAddress("vbf1_AK4_phi", &vbf1_AK4_phi);
    intree->SetBranchAddress("vbf1_AK4_m", &vbf1_AK4_m);
    //intree->SetBranchAddress("vbf2_AK4_pt", &vbf2_AK4_pt);
    intree->SetBranchAddress("vbf2_AK4_eta", &vbf2_AK4_eta);
    intree->SetBranchAddress("vbf2_AK4_phi", &vbf2_AK4_phi);
    intree->SetBranchAddress("vbf2_AK4_m", &vbf2_AK4_m);

    //intree->SetBranchAddress("vbf_m", &vbf_m);
    //intree->SetBranchAddress("vbf_deta", &vbf_deta);                                                                            

    //intree->SetBranchAddress("dibos_m", &dibos_m);

    intree->SetBranchAddress("nBtag_loose", &nBtag_loose);

    intree->SetBranchAddress("bosCent", &bosCent);
    intree->SetBranchAddress("zeppLep", &zeppLep);
    intree->SetBranchAddress("zeppHad", &zeppHad);

    if (var=="nom") {
      intree->SetBranchAddress("bos_PuppiAK8_m_sd0_corr", &bos_PuppiAK8_m_sd0_corr);
      intree->SetBranchAddress("vbf1_AK4_pt", &vbf1_AK4_pt);
      intree->SetBranchAddress("vbf2_AK4_pt", &vbf2_AK4_pt);
      intree->SetBranchAddress("vbf_m", &vbf_m);
      intree->SetBranchAddress("dibos_m", &dibos_m);
    }
    else if (var=="jesdn"){
      intree->SetBranchAddress("bos_PuppiAK8_m_sd0_corr_scaleDn", &bos_PuppiAK8_m_sd0_corr);
      intree->SetBranchAddress("vbf1_AK4_pt_scaleDn", &vbf1_AK4_pt);
      intree->SetBranchAddress("vbf2_AK4_pt_scaleDn", &vbf2_AK4_pt);
      intree->SetBranchAddress("vbf_m_scaleDn", &vbf_m);
      intree->SetBranchAddress("dibos_m_scaleDn", &dibos_m);
    }
    else {
      intree->SetBranchAddress("bos_PuppiAK8_m_sd0_corr_scaleUp", &bos_PuppiAK8_m_sd0_corr);
      intree->SetBranchAddress("vbf1_AK4_pt_scaleUp", &vbf1_AK4_pt);
      intree->SetBranchAddress("vbf2_AK4_pt_scaleUp", &vbf2_AK4_pt);
      intree->SetBranchAddress("vbf_m_scaleUp", &vbf_m);
      intree->SetBranchAddress("dibos_m_scaleUp", &dibos_m);
    }

    for (int i=0; i<intree->GetEntries(); i++) {
      intree->GetEntry(i);

      if ( !(nBtag_loose==0 && vbf1_AK4_pt>50 && vbf2_AK4_pt>50) ) continue;

      if ( vbf_m < 500) continue;
      if ( fabs(vbf1_AK4_eta - vbf2_AK4_eta)<4.0) continue;
      if ( dibos_m < 600) continue;

      bool isEle=false, isResolved=false, isZ=false;

      if (lep1_m == ELE_MASS) { isEle=true; }
      else if (lep1_m == MUON_MASS) { isEle=false; }
      else {
        cout << "lepton is not electron or muon! skipping" << endl;
        continue;
      }

      if (bos_PuppiAK8_m_sd0_corr > 0 && bos_AK4AK4_m < 0) { isResolved=false; }
      else if (bos_PuppiAK8_m_sd0_corr < 0 && bos_AK4AK4_m > 0) { isResolved=true; }
      else {
        cout << "both or neither of resolved and boosted mass is defined" << endl;
      }

      if (isResolved==true && (bos_AK4AK4_m<65 ||bos_AK4AK4_m>105)) continue;
      if (isResolved==true && (bos_j1_AK4_pt<30 || bos_j2_AK4_pt<30)) continue;
      if (isResolved==true && bos_AK4AK4_pt>200) continue;

      if (isResolved==false && (bos_PuppiAK8_m_sd0_corr<65 ||bos_PuppiAK8_m_sd0_corr>105)) continue;
      if (isResolved==false && bos_PuppiAK8_pt<200) continue;
      if (isResolved==false && bos_PuppiAK8_tau2tau1<0.55) continue;

      if (lep2_pt>0) isZ=true;

      if (isEle==true && (lep1_pt<35 || abs(lep1_eta)>2.5 || (abs(lep1_eta)>1.4442 && abs(lep1_eta)<1.566))) continue;
      if (isEle==false && (lep1_pt<35 || abs(lep1_eta)>2.4)) continue;

      if (isZ==true && isEle==true && (lep2_pt<30 || abs(lep2_eta)>2.5 || (abs(lep2_eta)>1.4442 && abs(lep2_eta)<1.566))) continue;
      if (isZ==true && isEle==false && (lep2_pt<30 || abs(lep2_eta)>2.4)) continue;

      if (isZ==false && isEle==true && MET<80) continue;
      if (isZ==false && isEle==false && MET<50) continue;

      float weight=(samp1.at(xx).xsec*samp1.at(xx).xsecCorr*lumi*genWeight*puWeight)/(1.0*(nTotal-2*nNeg));
      if (samp1.at(xx).sampname=="dataM" || samp1.at(xx).sampname=="dataE") weight=1.0;

      float deta=abs(vbf1_AK4_eta - vbf2_AK4_eta);

      for (int i=0; i<bins.nBinsMVBF; i++) {
	if (vbf_m < bins.MVBF_LE[i] || vbf_m > bins.MVBF_LE[i+1]) continue;
	
	for (int j=0; j<bins.nBinsDETA; j++) {

	  if (deta > bins.DETA_LE[j] && deta < bins.DETA_LE[j+1]) {

	    for (int k=0; k<bins.nBinsMVV; k++) {

	      if (!isZ && dibos_m > bins.MVV_LE[k] && dibos_m < bins.MVV_LE[k+1]) {

		float indexVal = k + j*bins.nBinsMVV + i*bins.nBinsDETA*bins.nBinsMVV + 0.5;

		if (isResolved) mWjj->Fill(indexVal, weight);
		else mWV->Fill(indexVal, weight);
	      }
	      else if (isZ && dibos_m > bins.MVV_LE[k] && dibos_m < bins.MVV_LE[k+1]) {

		float indexVal = k + j*bins.nBinsMVV + i*bins.nBinsDETA*bins.nBinsMVV + 0.5;		

		if (isResolved)  mZjj->Fill(indexVal, weight);
		else mZV->Fill(indexVal, weight);

	      }
	    }
	  }
	}
      } 
    }
  }
  
  out->Write();
  out->Close();

  
}
