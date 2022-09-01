//#############################################################################
void SetStyleVariables(TStyle *t2kStyle){

  // use plain black on white colors
  t2kStyle->SetFrameBorderMode(0);
  t2kStyle->SetCanvasBorderMode(0);
  t2kStyle->SetPadBorderMode(0);
  t2kStyle->SetPadColor(0);
  t2kStyle->SetCanvasColor(0);
  t2kStyle->SetStatColor(0);
  //t2kStyle->SetFillColor(0);
  t2kStyle->SetLegendBorderSize(1);

  // set the paper & margin sizes
  t2kStyle->SetPaperSize(20,26);
  t2kStyle->SetPadTopMargin(0.1);
  t2kStyle->SetPadRightMargin(0.15); //0.05 
  t2kStyle->SetPadBottomMargin(0.15);
  t2kStyle->SetPadLeftMargin(0.15);

  // use large Times-Roman fonts
  t2kStyle->SetTextFont(132);
  t2kStyle->SetTextSize(0.08);
  t2kStyle->SetLabelFont(132,"x");
  t2kStyle->SetLabelFont(132,"y");
  t2kStyle->SetLabelFont(132,"z");
  t2kStyle->SetLabelSize(0.05,"x");
  t2kStyle->SetTitleSize(0.06,"x");
  t2kStyle->SetLabelSize(0.05,"y");
  t2kStyle->SetTitleSize(0.06,"y");
  t2kStyle->SetLabelSize(0.05,"z");
  t2kStyle->SetTitleSize(0.06,"z");
  t2kStyle->SetLabelFont(132,"t");
  t2kStyle->SetTitleFont(132,"x");
  t2kStyle->SetTitleFont(132,"y");
  t2kStyle->SetTitleFont(132,"z");
  t2kStyle->SetTitleFont(132,"t");
  t2kStyle->SetTitleFillColor(0);
  t2kStyle->SetTitleX(0.25);
  t2kStyle->SetTitleFontSize(0.08);
  t2kStyle->SetTitleFont(132,"pad");

  //t2kStyle->SetPadGridX(true);
  //t2kStyle->SetPadGridY(true);

  // use bold lines and markers
  //  t2kStyle->SetMarkerStyle(20);
  t2kStyle->SetHistLineWidth(1.85);
  t2kStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  //  t2kStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  t2kStyle->SetOptTitle(0);
  t2kStyle->SetOptStat(0);
  t2kStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  //t2kStyle->SetPadTickX(1);
  t2kStyle->SetPadTickY(1);

  t2kStyle->SetPalette(1,0);  // use the nice red->blue palette
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  //Double_t green[NRGBs] = { 0.00, 0.00, 0.00, 0.00, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                   NCont);
  t2kStyle->SetNumberContours(NCont);

  // End of definition of t2kStyle
}

void analysis_getGraph(const char* input, TH1D* hGraph, double sigma, int tkiVar, double costhmu, double costhpi, double pmu, double ppi){
  TFile *f = new TFile(input,"OPEN");
  TTree *t = (TTree*)f->Get("TKITree");
  Int_t target;
  Float_t dptt, pN, dat, dphit, dpt, W_p_pi, W_nu_mu, costheta_adler, phi_adler;
  Float_t pimom, picostheta, mumom, mucostheta;
  Float_t weight;
  Float_t Q2;
  t->SetBranchAddress("dptt",&dptt);
  t->SetBranchAddress("pN",&pN);
  t->SetBranchAddress("dat",&dat);
  t->SetBranchAddress("dphit",&dphit);
  t->SetBranchAddress("dpt",&dpt);
  t->SetBranchAddress("W_p_pi",&W_p_pi);
  t->SetBranchAddress("W_nu_mu",&W_nu_mu);
  //t->SetBranchAddress("costheta_adler",&costheta_adler);
  //t->SetBranchAddress("phi_adler",&phi_adler);
  t->SetBranchAddress("target",&target);
  t->SetBranchAddress("weight",&weight);
  t->SetBranchAddress("Q2",&Q2);
  t->SetBranchAddress("pimom",&pimom);
  t->SetBranchAddress("picostheta",&picostheta);
  t->SetBranchAddress("mumom",&mumom);
  t->SetBranchAddress("mucostheta",&mucostheta);

  //TH1D* hdptt = new TH1D("hdptt", "Delta pTT", 99, -300, 300);
  //TH1D* hpN = new TH1D("hpN", "pN", 100, 0, 2000);
  //TH1D* hdat = new TH1D("hdat", "Delta AlphaT", 100, 0, 180);

  double nEvt_H = 0;
  std::cout<<input<<", nEntry = "<<t->GetEntries()<<std::endl;
  for (int i = 0; i < t->GetEntries(); i++){
        t->GetEntry(i);
        if (pimom<ppi || picostheta<costhpi || mumom<pmu || mucostheta<costhmu) continue;
        if (tkiVar==0) {hGraph->Fill(-Q2/1.e6,weight);}
        else if  (tkiVar==1) {
          if (mucostheta<0||mucostheta>0.8) continue;
          hGraph->Fill(mumom/1.e3,weight);
        }
        else if  (tkiVar==2) {
          if (mucostheta<0.8||mucostheta>0.85) continue;
          hGraph->Fill(mumom/1.e3,weight);
        }
        else if  (tkiVar==3) {
          if (mucostheta<0.85||mucostheta>0.9) continue;
          hGraph->Fill(mumom/1.e3,weight);
        }
        else if  (tkiVar==4) {
          if (mucostheta<0.9||mucostheta>1.) continue;
          hGraph->Fill(mumom/1.e3,weight);
        }
        else if  (tkiVar==5) {
          hGraph->Fill(pimom/1.e3,weight);
        }
        else if  (tkiVar==6) {
          hGraph->Fill(acos(picostheta),weight);
        }
        
  }
  //std::cout<<"nEvt_H = "<<nEvt_H<<std::endl;
  //if (tkiVar==2) for (int i=1;i<=hGraph->GetXaxis()->GetNbins();i++) hGraph->Fill(hGraph->GetBinCenter(i),(nEvt_H+0.)/hGraph->GetXaxis()->GetNbins());
  for (int j=1;j<=hGraph->GetXaxis()->GetNbins();j++) hGraph->SetBinContent(j,hGraph->GetBinContent(j)*sigma/hGraph->GetBinWidth(j)/1.e-38);
}

//#############################################################################


void makePlots_xsec_chi2(std::string inFile, int tkiVar, bool perdof){
  std::string fileList[] = {"1pi_NEUT.root"
                            ,"1pi_GENIE_R3.00.06_G18_01a.root"
                           ,"1pi_GENIE_R3.00.06_G18_10b.root"
               
                           ,"1pi_GiBUU.root"
                           
                           };
  const int nFiles = sizeof(fileList)/sizeof(fileList[0]);

  const int nDpttbins = 5;
  const double Dpttbins[nDpttbins+1] = {-700,-300,-100,100,300,700};
  const int nPnbins = 4;
  const double Pnbins[nPnbins+1] = {0,120,240,600,1500};
  const int nDatbins = 3;
  const double Datbins[nDatbins+1] = {0,60,120,180};
  double Bins[100];
  std::string tkiVarS;
  std::string ylabel;
  std::string xlabel;
  std::string chi2ylabel;
  int nBins = -1;
  double costhmu,  costhpi,  pmu,  ppi;
  double norm_costh=1;
  double xmin, xmax;
  //double chi2_dptt[nFiles] = {11.4,4.77,11.3,10.2,7.95,22.0};
  //double chi2_pN[nFiles] = {9.79,2.11,10.5,14.3,12.0,46.2};
  //double chi2_dat[nFiles] = {1.63,2.87,2.35,5.81,3.13,3.37};
  double chi2_tki[nFiles];
  TFile* fIn = new TFile(inFile.c_str());
  TH1D* h_unfolded_xsec;
  switch (tkiVar) {
      case 0:
        h_unfolded_xsec = (TH1D*)fIn->Get("Q2");
        ylabel = "#frac{d#sigma}{dQ^{2}}(10^{-38}cm^{2}/(GeV/c)^{2}/Nucleon)";
        xlabel = "Q^{2} (GeV/c)^{2}";
        costhmu=0.2;costhpi=0.2;pmu=200;ppi=200;
        xmin=0; xmax=2;
        break;
      case 1:
        h_unfolded_xsec = (TH1D*)fIn->Get("p_mu_theta_mu_0");
        ylabel = "#frac{d^{2}#sigma}{dp_{#mu}dcos#theta}(10^{-38}cm^{2}/(GeV/c)/Nucleon)";
        xlabel = "p_{#mu} (GeV/c)";
        costhmu=0.;costhpi=-1;pmu=0;ppi=0;
        norm_costh=0.8;
        xmin=0; xmax=2;
        break;
      case 2:
        h_unfolded_xsec = (TH1D*)fIn->Get("p_mu_theta_mu_1");
        ylabel = "#frac{d^{2}#sigma}{dp_{#mu}dcos#theta}(10^{-38}cm^{2}/(GeV/c)/Nucleon)";
        xlabel = "p_{#mu} (GeV/c)";
        costhmu=0.8;costhpi=-1;pmu=0;ppi=0;
        norm_costh=0.05;
        xmin=0; xmax=2;
        break;
      case 3:
        h_unfolded_xsec = (TH1D*)fIn->Get("p_mu_theta_mu_2");
        ylabel = "#frac{d^{2}#sigma}{dp_{#mu}dcos#theta}(10^{-38}cm^{2}/(GeV/c)/Nucleon)";
        xlabel = "p_{#mu} (GeV/c)";
        costhmu=0.85;costhpi=-1;pmu=0;ppi=0;
        norm_costh=0.05;
        xmin=0; xmax=2;
        break;
      case 4:
        h_unfolded_xsec = (TH1D*)fIn->Get("p_mu_theta_mu_3");
        ylabel = "#frac{d^{2}#sigma}{dp_{#mu}dcos#theta}(10^{-38}cm^{2}/(GeV/c)/Nucleon)";
        xlabel = "p_{#mu} (GeV/c)";
        costhmu=0.9;costhpi=-1;pmu=0;ppi=0;
        norm_costh=0.1;
        xmin=0; xmax=2;
        break;
      case 5:
        h_unfolded_xsec = (TH1D*)fIn->Get("Momentum_pion");
        ylabel = "#frac{d#sigma}{dp_{#pi}}(10^{-38}cm^{2}/(GeV/c)/Nucleon)";
        xlabel = "p_{#pi} (GeV/c)";
        costhmu=0.2;costhpi=0.2;pmu=200;ppi=0;
        xmin=0; xmax=3;
        break;
      case 6:
        h_unfolded_xsec = (TH1D*)fIn->Get("Theta_pion");
        ylabel = "#frac{d#sigma}{d#theta_{#pi}}(10^{-38}cm^{2}/(rad)/Nucleon)";
        xlabel = "#theta_{#pi} (rad)";
        costhmu=0.2;costhpi=0.;pmu=200;ppi=200;
        xmin=0; xmax=3.14;
        break;
      default:
        printf("***Warning: Not a valid TKI variable***\n");
        break;
  }

  //TFile* fIn = new TFile(inFile.c_str());
  //TH1D* sel_best_fit_results = (TH1D*)fIn -> Get("sel_best_fit");

  //TH1D* h_unfolded_xsec = new TH1D("","",nBins,Bins);
  TH1D* h_Model_xsec[nFiles];
  for (int i=0;i<nFiles;i++) {
    h_Model_xsec[i] = (TH1D*)h_unfolded_xsec->Clone();
    h_Model_xsec[i]->Reset();
  }

  //for (int i=0;i<nBins;i++){
    //h_unfolded_xsec->SetBinContent(i+1,sel_best_fit_results->GetBinContent(i+1));
    //h_unfolded_xsec->SetBinError(i+1,sel_best_fit_results->GetBinError(i+1));
  //}


  h_unfolded_xsec->SetLineColor(kRed);
  h_unfolded_xsec->SetMarkerStyle(kFullCircle);
  h_unfolded_xsec->GetXaxis()->SetTitle(xlabel.c_str());
  h_unfolded_xsec->GetYaxis()->SetTitle(ylabel.c_str());
  h_unfolded_xsec->GetYaxis()->SetRangeUser(0,h_unfolded_xsec->GetMaximum()*1.5);
  h_unfolded_xsec->GetXaxis()->SetRangeUser(xmin,xmax);
  h_unfolded_xsec->Draw("E0 X0)");

  /*TMatrixDSym* xsec_cov = (TMatrixDSym*)fIn -> Get("xsec_cov");
  TMatrixDSym xsec_cov_new;
  xsec_cov_new.ResizeTo(nBins, nBins);
  xsec_cov_new.Zero();
  for (int i=0;i<nBins;i++) for (int j=0;j<nBins;j++) {
    xsec_cov_new(i, j) += (*xsec_cov)[i][j]*1E84*(Bins[i+1]-Bins[i])*(Bins[j+1]-Bins[j]);
    //std::cout<<"xsec_cov_new("<<i<<","<<j<<")="<<xsec_cov_new(i, j)<<std::endl;
  }
  double det = 0;
  double total_add = 0;
  TDecompLU inv_test;
  TMatrixD inv_matrix(xsec_cov_new);
  while (!inv_test.InvertLU(inv_matrix, 1E-48, &det))
    {
        std::cerr <<  "In AnaFitParameters::SetCovarianceMatrix():\n"
                  << "Covariance matrix is non invertable. Determinant is " << det
                  << std::endl;
        for(int i = 0; i < xsec_cov_new.GetNrows(); ++i)
            xsec_cov_new(i,i) += 1e-6;
        total_add += 1e-6;
    }
  std::cout << "Added " << total_add << " to force positive-definite."
              << std::endl;
  xsec_cov_new.SetMatrixArray(inv_matrix.GetMatrixArray());
  std::cout << "Covariance matrix inverted successfully." << std::endl;*/

  double norm[nFiles] = {1,1,1,1e-39};
  for (int i=0;i<nFiles;i++) {
     analysis_getGraph(fileList[i].c_str(),h_Model_xsec[i],norm[i]/norm_costh, tkiVar, costhmu,  costhpi,  pmu,  ppi);
     h_Model_xsec[i]->SetLineColor(i+1);
     h_Model_xsec[i]->SetLineWidth(2);
     h_Model_xsec[i]->SetLineStyle(i+1);
     //if (i==0||i==1||i==3) 
     h_Model_xsec[i]->Draw("same hist");
     /*double chi2_model = 0;
     for (int j=1;j<=nBins;j++)  for (int k=1;k<=nBins;k++)  {
       double x = (sel_best_fit_results->GetBinContent(j)- h_Model_xsec[i]->GetBinContent(j))*1E42*(Bins[j]-Bins[j-1]);
       double y = (sel_best_fit_results->GetBinContent(k)- h_Model_xsec[i]->GetBinContent(k))*1E42*(Bins[k]-Bins[k-1]);
       if(!std::isnan(x)&&!std::isnan(y)) chi2_model += x*y*xsec_cov_new(j-1,k-1);
     }
     chi2_tki[i]=chi2_model;*/
  }

  
  const std::string xsecModel[nFiles] =
                            {"NEUT 5.4.0 (RFG)",
                              "GENIE R3.00.06_G18_01a (RS, RFG, hA)",
                             "GENIE R3.00.06_G18_10b (BS, LFG, hN)",
                             
                             "GiBUU 2019",
                             
                            };

  auto legend2 = new TLegend(0.27,0.59,0.83,0.88);
  legend2->AddEntry(h_unfolded_xsec,"T2K Data","lep");
  for (int i=0;i<nFiles;i++) {
    std::string legendString = Form("%s",xsecModel[i].c_str());
    //if (i==0||i==1||i==3) 
    legend2->AddEntry(h_Model_xsec[i],legendString.c_str(),"l");
  }
  legend2->SetBorderSize(0);
  legend2->SetFillStyle(0);
  if (tkiVar==0) legend2->Draw();
return;
  /*TCanvas *c2 = new TCanvas();
  TH1D* hXsecChi2 = new TH1D("","",nFiles,0,nFiles);
  for (int i=1;i<=nFiles;i++) {
    hXsecChi2->GetXaxis()->SetBinLabel(i,xsecModel[i-1].c_str());
    hXsecChi2->SetBinContent(i,chi2_tki[i-1]);
  }
  hXsecChi2->GetYaxis()->SetTitle(Form("#chi^{2}_{xsec}%s",chi2ylabel.c_str()));
  hXsecChi2->SetLineWidth(2);
  hXsecChi2->SetMinimum(0);
  hXsecChi2->Draw();

  TCanvas *c3 = new TCanvas();
  TH1D* hXsecChi2PerDof = (TH1D*)hXsecChi2->Clone();
  hXsecChi2PerDof->Scale(1./nBins);
  hXsecChi2PerDof->SetMinimum(0);
  hXsecChi2PerDof->GetYaxis()->SetTitle(Form("#chi^{2}_{xsec}/ndof%s",chi2ylabel.c_str()));
  hXsecChi2PerDof->Draw();*/
}

void makePlots_xsec_nuisance(){
  //**** Set Style for Plots ****
  TStyle *t2kstyle = new TStyle("T2K","T2K approved plots style");
  SetStyleVariables(t2kstyle);
  gROOT->SetStyle("T2K");

  TCanvas *c1 = new TCanvas();
  makePlots_xsec_chi2("Q2.root",0,false);

  TCanvas *c2 = new TCanvas();
  makePlots_xsec_chi2("PmuThetamu.root",1,false);
  TCanvas *c3 = new TCanvas();
  makePlots_xsec_chi2("PmuThetamu.root",2,false);
  TCanvas *c4 = new TCanvas();
  makePlots_xsec_chi2("PmuThetamu.root",3,false);
  TCanvas *c5 = new TCanvas();
  makePlots_xsec_chi2("PmuThetamu.root",4,false);
  TCanvas *c6 = new TCanvas();
  makePlots_xsec_chi2("MomentumPion.root",5,false);
  TCanvas *c7 = new TCanvas();
  makePlots_xsec_chi2("Thetapion.root",6,false);
return;
 

}
