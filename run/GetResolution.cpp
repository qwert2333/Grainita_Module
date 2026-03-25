vector<double> getRes(double En, TString s_scheme){

  double totE;

  TString s_En; s_En.Form("%d", (int)En);
  //TFile *rfile = new TFile("simdir/sim_pi-_"+s_En+"GeV_3mmGS.root", "read");
  TFile *rfile = new TFile("GrainitaCalo_"+s_scheme+"_gamma_"+s_En+"GeV.root", "read");
  cout<<"Run En point: "<<s_En<<endl;
  TTree *clutree = (TTree*)rfile->Get("eventTree");

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();

  TH1F* h_Egam = new TH1F("h_Egam", "h_Egam", 100, En*0.8, En*1.05);
  clutree->Draw("EdepCrystal/1000.>>h_Egam");

  vector<double> results;
  if(h_Egam->GetEntries()>50){
    //gStyle->SetOptFit(1111);
    h_Egam->Fit("gaus","Q");
    double Pmean = h_Egam->GetFunction("gaus")->GetParameter(1);
    double Pres = h_Egam->GetFunction("gaus")->GetParameter(2);
    h_Egam->Fit("gaus","Q","",Pmean-1.5*Pres, Pmean+1.5*Pres);
    Pmean = h_Egam->GetFunction("gaus")->GetParameter(1);
    Pres = h_Egam->GetFunction("gaus")->GetParameter(2);
    h_Egam->Fit("gaus","Q","",Pmean-1.5*Pres, Pmean+1.5*Pres);
    Pmean = h_Egam->GetFunction("gaus")->GetParameter(1);
    Pres = h_Egam->GetFunction("gaus")->GetParameter(2);
    double Pmean_err = h_Egam->GetFunction("gaus")->GetParError(1);
    double Pres_err = h_Egam->GetFunction("gaus")->GetParError(2);
    results.push_back(Pmean);
    results.push_back(Pres);
    results.push_back(Pmean_err);
    results.push_back(Pres_err);

    double chi2 = h_Egam->GetFunction("gaus")->GetChisquare();
    int ndof = h_Egam->GetFunction("gaus")->GetNDF();
    TLatex latex;
    latex.SetTextSize(0.04);
    latex.SetTextAlign(13);  //align at top
    latex.SetTextFont(42);
    latex.DrawLatexNDC(0.3, 0.8, "En = "+s_En+" GeV");
    latex.DrawLatexNDC(0.3, 0.7, Form("sigma/mean = %.2f/%.2f = %.2f", Pres, Pmean, Pres/Pmean));
    latex.DrawLatexNDC(0.3, 0.6, Form("chi2/ndf = %.3f/%d", chi2, ndof));
  }
  else{
    results.push_back(0);
    results.push_back(0);
    results.push_back(0);
    results.push_back(0);
  }

  //h_Egam->Draw();

  return results;
}


void GetResolution(){
  const int En_pointsize = 2; 
  TString s_scheme[5] = {"Pitch3mm_fiber66", "Pitch5mm_fiber40", "Pitch7mm_fiber28", "Pitch10mm_fiber20", "Pitch15mm_fiber13"};
  double En_point[En_pointsize] = {1., 10.};
  double En_over1_sq[En_pointsize] = {0.};

  double ELin[5][En_pointsize] = {0.};
  double ERes[5][En_pointsize] = {0.};
  double ERes_sq[5][En_pointsize] = {0.};

  for(int i=0; i<5; i++){
    for(int j=0; j<En_pointsize; j++){
      vector<double> tmp_result = getRes(En_point[j], s_scheme[i]);     
      ELin[i][j] = tmp_result[0];
      ERes[i][j] = tmp_result[1]/tmp_result[0];
      ERes_sq[i][j] = ERes[i][j]*ERes[i][j];
    }
  }


  for(int j=0; j<En_pointsize; j++){
    En_over1_sq[j] = 1./En_point[j]/En_point[j];
  }

  TCanvas *cv1 = new TCanvas("cv1", "cv1", 1000, 800);
  TLegend *l1 = new TLegend(0.2, 0.6, 0.5, 0.85);
  TGraph *gr_Eres[5];
  int color[5] = {1, 2, 4, 6, 8};
  for(int i=0; i<5; i++){
    gr_Eres[i] = new TGraph(En_pointsize, En_over1_sq, ERes_sq[i]);
    gr_Eres[i]->GetXaxis()->SetTitle("(1/E)^{2} / GeV^{-2}");
    gr_Eres[i]->GetYaxis()->SetTitle("(#sigma_{E} / E)^{2}");
    gr_Eres[i]->GetYaxis()->SetRangeUser(0, 2e-4);
    gr_Eres[i]->SetLineColor(color[i]);
    gr_Eres[i]->SetLineWidth(2);
    gr_Eres[i]->SetMarkerSize(1);
    if(i==0) gr_Eres[i]->Draw("ALP");
    else gr_Eres[i]->Draw("LP");

    TString s_pitch = s_scheme[i](0, s_scheme[i].First('_'));
    l1->AddEntry(gr_Eres[i], s_pitch, "l");

    gr_Eres[i]->Fit("pol1", "Q", "", En_over1_sq[0], En_over1_sq[1]);
    cout << s_pitch << " Stochastic term "<< sqrt(gr_Eres[i]->GetFunction("pol1")->GetParameter(1))*100 << "%, Constant term "<< sqrt(gr_Eres[i]->GetFunction("pol1")->GetParameter(0))*100<<"%"<<endl;
  }
  l1->SetBorderSize(0);
  l1->Draw();
  cv1->Draw();

}
