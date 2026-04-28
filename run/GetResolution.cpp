TRandom3 rndm(2077);
double response_map0[7][7] = { 1.92843, 3.07318, -6.09594, -7.20984, -4.89420, 2.93879, 3.34451,
0.78368, 2.67172, -6.05233, -5.98395, -0.56035, 7.52742, 3.75024,
2.64146, 1.26341, -4.37351, -7.44680, 0.64540, 1.20680, 4.89643,
1.87612, -2.61481, -7.82113, -5.65042, -2.84322, -3.10660, -0.05807,
1.69076, -1.83813, -5.50838, -5.34737, -1.53522, 1.40096, 3.20557,
7.67758, 2.89708, -4.82156, -0.85286, 2.44410, 2.02479, 8.05585,
5.23757, 2.79756, 0.96678, -3.64029, 0.22140, 4.70631, 6.38108,
};
double response_map1[7][7] = { 6.96517, 3.81123, 2.76829, 2.33175, 0.44024, 6.32259, 6.04933,
10.11910, 3.15215, -2.85353, -0.57657, 2.80487, 4.54815, 5.77608,
6.41911, 1.59115, -2.10792, -1.50952, -1.18852, 1.28123, 7.28820,
-3.84496, -0.66577, -5.11809, -6.39939, -6.39484, -1.13020, -2.88659,
-0.59270, -2.99046, -4.92565, -5.34925, -7.45341, -0.86279, -3.89184,
1.56377, -1.35745, -2.78436, -6.73437, -4.96666, 0.42888, 2.71463,
0.56692, -0.42992, 0.03094, -3.27039, -2.73893, 3.12867, 2.92165,
};
double sigma_nonuniform = 4.25;

double response_map_gen[7][7] = {0.};
for(int i=0; i<7; i++){
for(int j=0; j<7; j++){
  response_map_gen[i][j] = rndm.Gaus(0, sigma_nonuniform);
  cout<<response_map_gen[i][j]<<"  ";
}
cout<<endl;}


vector<double> getRes_smearOptical(double En, TString s_scheme){

  double totE;
  vector<int>* cellID = nullptr;
  vector<double>* vecEdep = nullptr;

  TString s_En; s_En.Form("%d", (int)En);
  TFile *rfile = new TFile("GrainitaCalo_"+s_scheme+"_gamma_"+s_En+"GeV.root", "read");
  cout<<"Run En point: "<<s_En<<endl;
  TTree *clutree = (TTree*)rfile->Get("eventTree");
  clutree->SetBranchAddress("EdepCrystal", &totE);
  clutree->SetBranchAddress("vecCellID", &cellID);
  clutree->SetBranchAddress("vecEdep", &vecEdep);

  double lightYield = 10;
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();

  //TH1F* h_Egam = new TH1F("h_Egam", "h_Egam", 100, En*1000*lightYield*0.8, En*1000*lightYield*1.02);
  TH1F* h_Egam = new TH1F("h_Egam", "h_Egam", 100, En*0.92, En*1.02);
  //clutree->Draw("EdepCrystal/1000.>>h_Egam");
  for(int ievt=0; ievt<clutree->GetEntries(); ievt++){
    clutree->GetEntry(ievt);
    totE = 0.;
    for(int ii=0; ii<cellID->size(); ii++){
      int id = cellID->at(ii); 
      int x = id%1000 - 500 + 84;    // Transeverse
      int y = id/1000 - 500 + 84;    // Transverse

      //double response = response_map0[x%7][y%7]; 
      //double response = response_map1[x%7][y%7]; 
      double response = response_map_gen[x%7][y%7]; 
      //double response = rndm.Gaus(0, sigma_nonuniform); 
      totE += (1+response/100.) * vecEdep->at(ii);
    }
    //double Nph = rndm.Gaus(totE*lightYield, sqrt(totE*lightYield));
    //h_Egam->Fill(Nph);
    h_Egam->Fill(totE/1000.);
  }

  vector<double> results;
  if(h_Egam->GetEntries()>50){
    gStyle->SetOptFit(1111);
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

    //TF1* fCB = new TF1("fCB", "crystalball", En*0.92, En*1.);
    //fCB->SetParameters(
    //  h_Egam->GetMaximum(),   // N
    //  h_Egam->GetMean(),      // mean
    //  h_Egam->GetRMS(),       // sigma
    //  1.5,               // alpha
    //  2.0                // n
    //);
    //h_Egam->Fit(fCB, "R");
    //double Pmean = fCB->GetParameter(1);
    //double Pres = fCB->GetParameter(2);
    //double Pmean_err = fCB->GetParError(1);
    //double Pres_err = fCB->GetParError(2);

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
    latex.DrawLatexNDC(0.2, 0.8, "En = "+s_En+" GeV");
    latex.DrawLatexNDC(0.2, 0.7, Form("sigma/mean = %.2f/%.2f = %.2f", Pres, Pmean, Pres/Pmean));
    latex.DrawLatexNDC(0.2, 0.6, Form("chi2/ndf = %.3f/%d", chi2, ndof));
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

vector<double> getRes(double En, TString s_scheme){

  double totE;

  TString s_En; s_En.Form("%d", (int)En);
  //TFile *rfile = new TFile("simdir/sim_pi-_"+s_En+"GeV_3mmGS.root", "read");
  TFile *rfile = new TFile("GrainitaCalo_"+s_scheme+"_gamma_"+s_En+"GeV.root", "read");
  cout<<"Run En point: "<<s_En<<endl;
  TTree *clutree = (TTree*)rfile->Get("eventTree");

  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
  c1->cd();

  TH1F* h_Egam = new TH1F("h_Egam", "h_Egam", 100, En*0.9, En*1.02);
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
  const int En_pointsize = 4; 
  TString s_scheme[5] = {"Pitch3mm_fiber66", "Pitch5mm_fiber40", "Pitch7mm_fiber24", "Pitch10mm_fiber20", "Pitch15mm_fiber13"};
  double En_point[En_pointsize] = {1., 5., 10., 25.};
  double En_over1_sq[En_pointsize] = {0.};

  double ELin[5][En_pointsize] = {0.};
  double ERes[5][En_pointsize] = {0.};
  double ERes_sq[5][En_pointsize] = {0.};

  //for(int i=0; i<5; i++){
  int i=2;
    for(int j=0; j<En_pointsize; j++){
      vector<double> tmp_result = getRes_smearOptical(En_point[j], s_scheme[i]);     
      ELin[i][j] = tmp_result[0];
      ERes[i][j] = tmp_result[1]/tmp_result[0];
      ERes_sq[i][j] = ERes[i][j]*ERes[i][j];
    }
  //}


  for(int j=0; j<En_pointsize; j++){
    En_over1_sq[j] = 1./En_point[j];
  }

  TCanvas *cv1 = new TCanvas("cv1", "cv1", 1000, 800);
  TLegend *l1 = new TLegend(0.2, 0.6, 0.5, 0.85);
  TGraph *gr_Eres[5];
  int color[5] = {1, 2, 4, 6, 8};
  //for(int i=0; i<5; i++){
    gr_Eres[i] = new TGraph(En_pointsize, En_over1_sq, ERes_sq[i]);
    gr_Eres[i]->GetXaxis()->SetTitle("(1/E) / GeV^{-1}");
    gr_Eres[i]->GetYaxis()->SetTitle("(#sigma_{E} / E)^{2}");
    //gr_Eres[i]->GetYaxis()->SetRangeUser(0, 2e-4);
    gr_Eres[i]->SetLineColor(color[i]);
    gr_Eres[i]->SetLineWidth(2);
    gr_Eres[i]->SetMarkerSize(10);
    //if(i==0) gr_Eres[i]->Draw("ALP");
    //else gr_Eres[i]->Draw("LP");
    gr_Eres[i]->Draw("ALP");

    TString s_pitch = s_scheme[i](0, s_scheme[i].First('_'));
    l1->AddEntry(gr_Eres[i], s_pitch, "l");

    gr_Eres[i]->Fit("pol1", "Q", "", En_over1_sq[0], En_over1_sq[1]);
    cout << s_pitch << " Stochastic term "<< sqrt(gr_Eres[i]->GetFunction("pol1")->GetParameter(1))*100 << "%, Constant term "<< sqrt(gr_Eres[i]->GetFunction("pol1")->GetParameter(0))*100<<"%"<<endl;
  //}
  l1->SetBorderSize(0);
  l1->Draw();
  cv1->Draw();

}


void DrawEnHist(){
  gStyle->SetLegendTextSize(0.04);
  gStyle->SetLegendFont(42);

  int pitch[5] = {4, 1, 7, 10, 15};
  TString s_scheme[5] = {"Box4/GrainitaCalo_Pitch7mm_fiber14", "Box2/GrainitaCalo_Pitch7mm_fiber28", "Box1/GrainitaCalo_Pitch7mm_fiber56", "Pitch7mm_fiber56", "Pitch10mm_fiber40"};
  TFile *rfile[5];
  TTree *rtree[5];
  TH1D *h_EnDep[4][5];
  //TFile *wfile = new TFile("Plots.root", "recreate");
  for(int i=0; i<3; i++){
    rfile[i] = new TFile(s_scheme[i]+"_gamma_GunPlane30cm_1GeV.root", "read");
    rtree[i] = (TTree*)rfile[i]->Get("eventTree");
    h_EnDep[0][i] = new TH1D(Form("h_EnCrys_pitch%dmm", pitch[i]), "", 100, 850, 1050);
    h_EnDep[1][i] = new TH1D(Form("h_EnFiber_pitch%dmm", pitch[i]), "", 100, 0, 100);
    h_EnDep[2][i] = new TH1D(Form("h_EnFrame_pitch%dmm", pitch[i]), "", 100, 0, 50);
    h_EnDep[3][i] = new TH1D(Form("h_EnSum_pitch%dmm", pitch[i]), "", 100, 700, 1100);

    rtree[i]->Draw(Form("EdepCrystal>>h_EnCrys_pitch%dmm", pitch[i] ) );
    rtree[i]->Draw(Form("EdepFiberCore+EdepFiberClad>>h_EnFiber_pitch%dmm", pitch[i] ) );
    rtree[i]->Draw(Form("EdepCarbonFrame>>h_EnFrame_pitch%dmm", pitch[i] ) );
    rtree[i]->Draw(Form("EdepCrystal+EdepFiberCore+EdepFiberClad+EdepCarbonFrame>>h_EnSum_pitch%dmm", pitch[i] ) );
  }

  TCanvas *cv1 = new TCanvas("cv1", "cv1", 1000, 800);
  TLegend *l1 = new TLegend(0.2, 0.7, 0.65, 0.85);
  int color[5] = {1, 2, 4, 6, 8};
  double pitch_double[5] = {0.};
  for(int i=0; i<3; i++){
    h_EnDep[2][i]->SetStats(0);
    h_EnDep[2][i]->GetXaxis()->SetTitle("E_{Carbon Frame} / MeV");
    h_EnDep[2][i]->SetLineColor(color[i]);
    h_EnDep[2][i]->SetLineWidth(4);
    h_EnDep[2][i]->SetMarkerSize(0);
    h_EnDep[2][i]->Draw("same");
  }
  l1->AddEntry(h_EnDep[2][0], Form("4#times4 Boxes frame, Mean %.2f", h_EnDep[2][0]->GetMean()), "l");
  l1->AddEntry(h_EnDep[2][1], Form("2#times2 Boxes frame, Mean %.2f", h_EnDep[2][1]->GetMean()), "l");
  l1->AddEntry(h_EnDep[2][2], Form("No inside frame, Mean %.2f",      h_EnDep[2][2]->GetMean()), "l");
  l1->SetBorderSize(0);
  l1->Draw();
  cv1->Draw();

}
