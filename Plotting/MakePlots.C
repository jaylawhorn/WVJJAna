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

  TFile *fVBFEWK = new TFile("VBF_EWK.root","READ");
  TFile *fVBFQCD = new TFile("VBF_QCD.root","READ");
  TFile *fTop    = new TFile("Top.root",    "READ");
  TFile *fWJets  = new TFile("WJets.root",  "READ");
  TFile *fDYJets = new TFile("DYJets.root", "READ");
  TFile *fDataM  = new TFile("DataM.root",  "READ");

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  gStyle->SetOptStat(0);
  c1->SetLogy();

  TLegend *l = new TLegend(0.6,0.6,0.85,0.85);

  TString ii[9] = { "00", "01", "02", "10", "11", "12", "20", "21", "22" };

  for (int x=0; x<9; x++) {

    TH1D *hVBFEWK = (TH1D*) fVBFEWK->Get(Form("hZmjj_%s_VBF_EWK",ii[x].Data()));
    TH1D *hVBFQCD = (TH1D*) fVBFQCD->Get(Form("hZmjj_%s_VBF_QCD",ii[x].Data()));
    TH1D *hTop    = (TH1D*) fTop->Get(Form("hZmjj_%s_Top",ii[x].Data()));
    TH1D *hWJets  = (TH1D*) fWJets->Get(Form("hZmjj_%s_WJets",ii[x].Data()));
    TH1D *hDYJets = (TH1D*) fDYJets->Get(Form("hZmjj_%s_DYJets",ii[x].Data()));
    TH1D *hDataM  = (TH1D*) fDataM->Get(Form("hZmjj_%s_dataM",ii[x].Data()));
    
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
    
    hTop->GetYaxis()->SetRangeUser(0.1,10000);
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
    
    c1->SaveAs(Form("Zmjj_%s.png",ii[x].Data()));
  }

}
