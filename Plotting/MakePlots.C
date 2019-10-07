#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
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

void MakePlots() {

  TFile *fVBFEWK = new TFile("VBF_EWK_2016loose_drawM.root","READ");
  TFile *fVBFQCD = new TFile("VBF_QCD_2016loose_drawM.root","READ");
  TFile *fTop    = new TFile("Top_2016loose_drawM.root",    "READ");
  TFile *fWJets  = new TFile("WJets_2016loose_drawM.root",  "READ");
  TFile *fDYJets = new TFile("DYJets_2016loose_drawM.root", "READ");
  TFile *fDataM  = new TFile("DataM_2016loose_drawM.root",  "READ");

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  gStyle->SetOptStat(0);
  c1->SetLogy();

  TLegend *l = new TLegend(0.75,0.75,0.9,0.9);

  TString jj[4] = { "Wjj","WV","Zjj","ZV" };
  TString ii[7] = { "mVV","dETA","mVBF","mVJ","MET","ptL","mVL"};

  for (int y=0; y<4; y++) {
    for (int x=0; x<7; x++) {

      if (x==0) continue;
      if (y>1 && (x==4 || x==5)) continue; 
      if (y<2 && x==6) continue;

      TH1D *hVBFEWK = (TH1D*) fVBFEWK->Get(Form("VBF_EWK_%s_%s",ii[x].Data(),jj[y].Data()));
      TH1D *hVBFQCD = (TH1D*) fVBFQCD->Get(Form("VBF_QCD_%s_%s",ii[x].Data(),jj[y].Data()));
      TH1D *hTop    = (TH1D*) fTop->Get(Form("Top_%s_%s",ii[x].Data(),jj[y].Data()));
      TH1D *hWJets  = (TH1D*) fWJets->Get(Form("WJets_%s_%s",ii[x].Data(),jj[y].Data()));
      TH1D *hDYJets = (TH1D*) fDYJets->Get(Form("DYJets_%s_%s",ii[x].Data(),jj[y].Data()));
      TH1D *hDataM  = (TH1D*) fDataM->Get(Form("dataM_%s_%s",ii[x].Data(),jj[y].Data()));

      hVBFEWK->SetFillColor(840);
      hVBFEWK->SetLineColor(840);
      hVBFQCD->SetFillColor(400);
      hVBFQCD->SetLineColor(400);
      hTop->SetFillColor(592);
      hTop->SetLineColor(592);
      hWJets->SetFillColor(924);
      hWJets->SetLineColor(924);
      hDYJets->SetFillColor(4);
      hDYJets->SetLineColor(4);
      
      hDataM->SetLineColor(kBlack);
      hDataM->SetMarkerColor(kBlack);
      hDataM->SetMarkerStyle(20);
      hDataM->SetMarkerSize(1.2);
      
      hDYJets->Add(hVBFEWK);
      hVBFQCD->Add(hDYJets);
      hWJets->Add(hVBFQCD);
      hTop->Add(hWJets);
      
      hTop->SetTitle("");
      
      float minaxis = hTop->GetMinimum() < hDataM->GetMinimum() ? hTop->GetMinimum() : hDataM->GetMinimum();
      if (minaxis<0.01) minaxis=1;
      float maxaxis = hTop->GetMaximum() > hDataM->GetMaximum() ? hTop->GetMaximum() : hDataM->GetMaximum();
      hTop->GetYaxis()->SetRangeUser(0.5*minaxis,50*maxaxis);
      hTop->Draw("hist");
      hWJets->Draw("histsame");
      hVBFQCD->Draw("histsame");
      hDYJets->Draw("histsame");
      hVBFEWK->Draw("histsame");
      hDataM->Draw("same");
      
      l->Clear();
      l->AddEntry(hDataM, "Data","pl");
      l->AddEntry(hTop,   "Top","f");
      l->AddEntry(hWJets, "WJets","f");
      l->AddEntry(hVBFQCD, "VV QCD","f");
      l->AddEntry(hDYJets, "DYJets","f");
      l->AddEntry(hVBFEWK, "VBS EWK","f");
      l->Draw();
    
      c1->SaveAs(Form("%s_%s.png",ii[x].Data(),jj[y].Data()));
    }
  }
  
}
