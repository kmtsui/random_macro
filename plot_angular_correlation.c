int counter;

TH1D* plot(std::string config, double sys_lo, double sys_hi)
{
    gStyle->SetOptStat(0);

    //std::string config = "fullring";
    //std::string config = "noring1";

    int nmPMT_min = 100;
    int step = 100;
    const int nstep = 40;
    TH1D* hist_resol_bias = new TH1D(config.c_str(),config.c_str(),nstep,nmPMT_min-step/2,nmPMT_min+(nstep-0.5)*step);
    TH1D* hist_resol = new TH1D("","",nstep,nmPMT_min-step/2,nmPMT_min+(nstep-0.5)*step);
    TH1D* hist_resol_bias_1ns = new TH1D("","",nstep,nmPMT_min-step/2,nmPMT_min+(nstep-0.5)*step);
    TH1D* hist_resol_bias_led = new TH1D("","",nstep,nmPMT_min-step/2,nmPMT_min+(nstep-0.5)*step);
    TH1D* hist_globalcc = new TH1D("","",nstep,nmPMT_min-step/2,nmPMT_min+(nstep-0.5)*step);
    double alpha_nom[nstep];
    for (int i=nstep-1;i>=0;i--)
    {
        int nmPMT = nmPMT_min+i*step;
        TFile* f = new TFile(Form("diffuser4_400nm_nominal_mPMT_%i_%s.root",nmPMT,config.c_str()),"OPEN");
        TH1D* hist_alpha_result = (TH1D*)f->Get("hist_alpha_result");
        TH1D* hist_alpha_error_final = (TH1D*)f->Get("hist_alpha_error_final");
        TVectorT<double>* res_globalcc = (TVectorT<double>*)f->Get("res_globalcc");
        TMatrixTSym<double>* res_cor_matrix = (TMatrixTSym<double>*)f->Get("res_cor_matrix");
        double value = (*res_cor_matrix)[1][4];
        // hist_globalcc->SetBinContent(i+1,(*res_globalcc)[1]);
        // alpha_nom[i] = 10800;//hist_alpha_result->GetBinContent(1);//10800;//
        // hist_resol_bias->SetBinContent(i+1,hist_alpha_result->GetBinContent(1)/alpha_nom[i]);
        // hist_resol_bias->SetBinError(i+1,hist_alpha_error_final->GetBinContent(1)/alpha_nom[i]);

        hist_resol->SetBinContent(i+1,fabs(value)+sys_lo);
    }

    hist_resol->GetXaxis()->SetTitle("#mPMT");
    //hist_resol->GetYaxis()->SetTitle("Estimated stat. + sys. uncert. in L_{#alpha}");
    hist_resol->GetYaxis()->SetTitle("Largest correlation between L_{#alpha} and A(#theta)");
    hist_resol->GetYaxis()->SetTitleOffset(1.5);
    hist_resol->GetYaxis()->SetRangeUser(0.45,0.65);
    
    hist_resol->SetLineColor(counter+1);
    hist_resol->SetLineStyle(counter+1);
    hist_resol->SetLineWidth(3);
    counter++;

    return hist_resol;
}

void plot_angular_correlation()
{
    gStyle->SetOptStat(0);
    
    counter = 0;
    
    // TH1D* h_full = plot("fullring",0.017,0.009);
    // TH1D* h_half = plot("halfring",0.024,0.011);
    // TH1D* h_noout = plot("noring1",0.046,0.012);
    // TH1D* h_nomid = plot("noring2",0.024,0.010);
    TH1D* h_full = plot("fullring_pol",0.00,0.005);
    TH1D* h_half = plot("halfring_pol",0.00,0.0075);
    TH1D* h_noout = plot("noring1_pol",0.05,0.012);
    TH1D* h_nomid = plot("noring2_pol",0.0,0.006);

    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    h_full->Draw("hist same");
    h_half->Draw("hist same");
    h_noout->Draw("hist same");
    h_nomid->Draw("hist same");
    TLegend* legend = new TLegend(0.45,0.65,0.75,0.9);
    legend->AddEntry(h_full,"Full ring","l");
    legend->AddEntry(h_half,"Half ring","l");
    legend->AddEntry(h_noout,"No outer","l");
    legend->AddEntry(h_nomid,"No middle","l");
    legend->Draw("same");
    c1->SaveAs("Lalpha_angular_correlation.pdf");
}