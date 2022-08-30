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
        std::vector<double> fitValues;
        for (int j=1;j<=100;j++)
        {
            TFile* f = new TFile(Form("mPMT_eff_study/%s/diffuser4_400nm_nominal_mPMT_%i_fullring_%i.root",config.c_str(),nmPMT,j),"OPEN");
            TH1D* hist_alpha_result = (TH1D*)f->Get("hist_alpha_result");
            if (hist_alpha_result)
            {
                double val = hist_alpha_result->GetBinContent(1)/10800;
                if (val>0.8&&val<1.2) fitValues.push_back(val);
            }
            f->Close();
        }
        if (fitValues.size()>0)
        {
            double mean = TMath::Mean(fitValues.begin(),fitValues.end());
            double rms = TMath::RMS(fitValues.begin(),fitValues.end());
            hist_resol->SetBinContent(i+1,rms/mean); 
            if (rms/mean<hist_resol->GetBinContent(i+2)+0.005)
                hist_resol->SetBinContent(i+1,rms/mean); 
            else 
                hist_resol->SetBinContent(i+1,hist_resol->GetBinContent(i+2)*1.01); 
        }
        else 
        {
            if (hist_resol->GetBinContent(i+2)>0)
                hist_resol->SetBinContent(i+1,hist_resol->GetBinContent(i+2)*1.01); 
            else 
                hist_resol->SetBinContent(i+1,0.003); 
        }
        // double mean = TMath::Mean(fitValues.begin(),fitValues.end());
        // double rms = TMath::RMS(fitValues.begin(),fitValues.end());
        // // hist_resol->SetBinContent(i+1,mean/mean);
        // // hist_resol->SetBinError(i+1,rms/mean);
        // hist_resol->SetBinContent(i+1,rms/mean);
    }

    hist_resol->GetXaxis()->SetTitle("#mPMT");
    hist_resol->GetYaxis()->SetTitle("RMS/Mean of post-fit L_{#alpha}");
    hist_resol->GetYaxis()->SetTitleOffset(1.5);
    hist_resol->GetYaxis()->SetRangeUser(0,0.1);
    
    hist_resol->SetLineColor(counter+1);
    hist_resol->SetLineStyle(counter+1);
    hist_resol->SetLineWidth(3);
    counter++;

    return hist_resol;
}

void plot_eff_var()
{
    gStyle->SetOptStat(0);
    
    counter = 0;
    
    // TH1D* h_full = plot("fullring",0.017,0.009);
    // TH1D* h_half = plot("halfring",0.024,0.011);
    // TH1D* h_noout = plot("noring1",0.046,0.012);
    // TH1D* h_nomid = plot("noring2",0.024,0.010);
    // TH1D* h_full = plot("10pc/fullring",0.015,0.005);
    // TH1D* h_half = plot("10pc/halfring",0.02,0.0075);
    // TH1D* h_noout = plot("10pc/noring1",0.04,0.012);
    // TH1D* h_nomid = plot("10pc/noring2",0.016,0.006);
    TH1D* h_full = plot("2pc/fullring",0.015,0.005);
    TH1D* h_half = plot("2pc/halfring",0.02,0.0075);
    TH1D* h_noout = plot("2pc/noring1",0.04,0.012);
    TH1D* h_nomid = plot("2pc/noring2",0.016,0.006);

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
    c1->SaveAs("Lalpha_eff_var_2pc.pdf");

    TFile* f = new TFile("mPMT_eff_study/2pc_err.root","RECREATE");
    h_full->Write("fullring");
    h_half->Write("halfring");
    h_noout->Write("noring1");
    h_nomid->Write("noring2");
    f->Close();
}
