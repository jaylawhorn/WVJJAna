{

  TFile *f = new TFile("comb_0.root","UPDATE");

  TH1D *data_obs_mWjj = (TH1D*) DYJets_mWjj->Clone("data_obs_mWjj");

  data_obs_mWjj->Add(Top_mWjj);
  data_obs_mWjj->Add(VBF_EWK_mWjj);
  data_obs_mWjj->Add(VBF_QCD_mWjj);
  data_obs_mWjj->Add(WJets_mWjj);

  TH1D *data_obs_mWV = (TH1D*) DYJets_mWV->Clone("data_obs_mWV");

  data_obs_mWV->Add(Top_mWV);
  data_obs_mWV->Add(VBF_EWK_mWV);
  data_obs_mWV->Add(VBF_QCD_mWV);
  data_obs_mWV->Add(WJets_mWV);

  TH1D *data_obs_mZjj = (TH1D*) DYJets_mZjj->Clone("data_obs_mZjj");

  data_obs_mZjj->Add(Top_mZjj);
  data_obs_mZjj->Add(VBF_EWK_mZjj);
  data_obs_mZjj->Add(VBF_QCD_mZjj);
  data_obs_mZjj->Add(WJets_mZjj);

  TH1D *data_obs_mZV = (TH1D*) DYJets_mZV->Clone("data_obs_mZV");

  data_obs_mZV->Add(Top_mZV);
  data_obs_mZV->Add(VBF_EWK_mZV);
  data_obs_mZV->Add(VBF_QCD_mZV);
  data_obs_mZV->Add(WJets_mZV);

  f->Write();
  

}
