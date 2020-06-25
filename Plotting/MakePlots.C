#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLine.h>
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

  TFile *fVBFEWK = new TFile("VBF_EWK_2018_vjets_mu.root","READ");
  TFile *fVBFQCD = new TFile("VBF_QCD_2018_vjets_mu.root","READ");
  TFile *fTop    =     new TFile("Top_2018_vjets_mu.root","READ");
  TFile *fWJets  =   new TFile("WJets_2018_vjets_mu.root","READ");
  TFile *fDYJets =  new TFile("DYJets_2018_vjets_mu.root","READ");
  TFile *fDataM  =   new TFile("DataM_2018_vjets_mu.root","READ");

  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  gStyle->SetOptStat(0);
  c1->Divide(1,2,0,0);
  c1->cd(1)->SetPad(0,0.3,1.0,1.0);
  c1->cd(1)->SetTopMargin(0.1);
  c1->cd(1)->SetBottomMargin(0.01);
  c1->cd(1)->SetLeftMargin(0.18);  
  c1->cd(1)->SetRightMargin(0.07);  
  c1->cd(1)->SetTickx(1);
  c1->cd(1)->SetTicky(1);  
  c1->cd(2)->SetPad(0,0,1.0,0.3);
  c1->cd(2)->SetTopMargin(0.05);
  c1->cd(2)->SetBottomMargin(0.45);
  c1->cd(2)->SetLeftMargin(0.18);
  c1->cd(2)->SetRightMargin(0.07);
  c1->cd(2)->SetTickx(1);
  c1->cd(2)->SetTicky(1);
  gStyle->SetTitleOffset(1.400,"Y");

  c1->cd(1)->SetLogy();
  TLegend *l = new TLegend(0.6,0.68,0.9,0.85);
  l->SetBorderSize(0);

  TString jj[4] = { "Wjj","WV","Zjj","ZV" };
  TString ii[12] = { "mVV","dETA","mVBF","mVJ","MET","ptL","mVL","ETA1","ETA2","ETA_lep","ETA_bos","nJet50"};

  for (int y=0; y<4; y++) {
    for (int x=0; x<12; x++) {

      //if (x==0) continue;
      if (y>1 && x==4) continue; 
      if (y>1 && x==5) continue; 
      if (y<2 && x==6) continue;
      if (y>1 && x>9) continue;

      TH1D *hVBFEWK = (TH1D*) fVBFEWK->Get(Form("VBF_EWK_%s_%s",ii[x].Data(),jj[y].Data()));
      TH1D *hVBFQCD = (TH1D*) fVBFQCD->Get(Form("VBF_QCD_%s_%s",ii[x].Data(),jj[y].Data()));
      TH1D *hTop    = (TH1D*) fTop->Get(Form("Top_%s_%s",ii[x].Data(),jj[y].Data()));
      TH1D *hWJets  = (TH1D*) fWJets->Get(Form("WJets_%s_%s",ii[x].Data(),jj[y].Data()));
      TH1D *hDYJets = (TH1D*) fDYJets->Get(Form("DYJets_%s_%s",ii[x].Data(),jj[y].Data()));
      TH1D *hDataM  = (TH1D*) fDataM->Get(Form("dataM_%s_%s",ii[x].Data(),jj[y].Data()));
      TH1D *hDiff   = (TH1D*) hDataM->Clone("diff");

      //hVBFEWK->SetFillColor(840);
      //hVBFEWK->SetLineColor(840);
      hVBFEWK->SetLineColor(kBlack);
      hVBFEWK->SetLineStyle(4);
      hVBFEWK->SetLineWidth(3);
      hVBFQCD->SetFillColor(TColor::GetColor(250,202,255));
      hVBFQCD->SetLineColor(TColor::GetColor(250,202,255));
      hTop->SetFillColor(592);
      hTop->SetLineColor(592);
      hWJets->SetFillColor(TColor::GetColor(222,90,106));
      hWJets->SetLineColor(TColor::GetColor(222,90,106));
      hDYJets->SetFillColor(TColor::GetColor(222,90,106));
      hDYJets->SetLineColor(TColor::GetColor(222,90,106));
      
      hDataM->SetLineColor(kBlack);
      hDataM->SetMarkerColor(kBlack);
      hDataM->SetMarkerStyle(20);
      hDataM->SetMarkerSize(1.2);

      if (x==0) {
	std::cout << "V+jets: " << hWJets->Integral()+hDYJets->Integral() << std::endl;
	std::cout << "Top:    " << hTop->Integral() << std::endl;
	std::cout << "QCD-VV: " << hVBFQCD->Integral() << std::endl;
	std::cout << "VBS:    " << hVBFEWK->Integral() << std::endl;
      }

      //hVBFQCD->Add(hVBFEWK);      
      hDYJets->Add(hVBFQCD);
      hWJets->Add(hDYJets);
      hTop->Add(hWJets);

      if (x==0) {
	std::cout << "*************" << std::endl;
	std::cout << "Data:   " << hDataM->Integral() << std::endl;
	std::cout << "Back:   " << hTop->Integral() << std::endl;
	std::cout << "Signal: " << hVBFEWK->Integral() << std::endl;
      }

      hVBFEWK->Scale(10.0);
      
      hDiff->Divide(hTop);
      
      hTop->SetTitle("");
      
      float minaxis = hTop->GetMinimum() < hDataM->GetMinimum() ? hTop->GetMinimum() : hDataM->GetMinimum();
      if (minaxis<0.01) minaxis=1;
      float maxaxis = hTop->GetMaximum() > hDataM->GetMaximum() ? hTop->GetMaximum() : hDataM->GetMaximum();
      hTop->GetYaxis()->SetRangeUser(0.1,50*maxaxis);
      c1->cd(1);
      hTop->Draw("hist");
      hWJets->Draw("histsame");
      hDYJets->Draw("histsame");
      hVBFQCD->Draw("histsame");
      hVBFEWK->Draw("histsame");
      hDataM->Draw("same");
      c1->RedrawAxis();
      
      l->Clear();
      l->AddEntry(hDataM, "Data","pl");
      l->AddEntry(hTop,   "Top","f");
      l->AddEntry(hDYJets, "VJets","f");
      l->AddEntry(hVBFQCD, "VV QCD","f");
      l->AddEntry(hVBFEWK, "VBS EWK (x10)","l");
      l->Draw();

      c1->cd(2);

      hDiff->SetMarkerStyle(8);
      hDiff->SetMarkerSize(1.25);
      if (jj[y]=="Zjj" || jj[y]=="ZV") {
	hDiff->GetYaxis()->SetRangeUser(0.0,3.0);
      } else {
	hDiff->GetYaxis()->SetRangeUser(0.0,3.0);
      }
      hDiff->GetXaxis()->SetTitleOffset(0.9);
      hDiff->GetXaxis()->SetTitleSize(0.15);
      hDiff->GetXaxis()->SetLabelSize(0.15);
      hDiff->GetYaxis()->SetTitleSize(0.1);
      hDiff->GetYaxis()->SetTitleOffset(0.5);
      hDiff->GetYaxis()->CenterTitle(true);
      hDiff->GetYaxis()->SetLabelSize(0.1);
      hDiff->GetYaxis()->SetTitle("Data/MC");
      hDiff->SetTitle("");

      //hWJets->Draw("axis");
      //hDiff->GetYaxis()->SetRangeUser(0,2);

      hDiff->GetXaxis()->SetNdivisions(505);
      hDiff->GetYaxis()->SetNdivisions(505);

      hDiff->Draw("pe");

      float minX = hDiff->GetBinLowEdge(1);
      float maxX = hDiff->GetBinLowEdge(hDiff->GetNbinsX())+hDiff->GetBinWidth(hDiff->GetNbinsX());
      
      TLine lup(minX, 1, maxX, 1);
      //TLine ldn(minX, -5, maxX, -5);
      lup.Draw("same");
      //ldn.Draw("same");

      c1->SaveAs(Form("%s_%s_2018_vjets_mu_LO_NJET.png",ii[x].Data(),jj[y].Data()));
    }
  }
  
}
