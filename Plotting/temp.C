{

  TFile *f = new TFile("comb_XX.root","UPDATE");

  TH1D *data_obs_mWjj = (TH1D*) dataM_mWjj->Clone("data_obs_mWjj");

  data_obs_mWjj->Add(dataE_mWjj);

  TH1D *data_obs_mWV = (TH1D*) dataM_mWV->Clone("data_obs_mWV");

  data_obs_mWV->Add(dataE_mWV);

  TH1D *data_obs_mZjj = (TH1D*) dataM_mZjj->Clone("data_obs_mZjj");

  data_obs_mZjj->Add(dataE_mZjj);

  TH1D *data_obs_mZV = (TH1D*) dataM_mZV->Clone("data_obs_mZV");

  data_obs_mZV->Add(dataE_mZV);

  f->Write();
  

}
