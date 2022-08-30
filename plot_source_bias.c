int counter;

TH1D* plot(std::string config)
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
        hist_globalcc->SetBinContent(i+1,(*res_globalcc)[1]);
        alpha_nom[i] = 10800;//hist_alpha_result->GetBinContent(1);//10800;//
        hist_resol_bias->SetBinContent(i+1,hist_alpha_result->GetBinContent(1)/alpha_nom[i]);
        hist_resol_bias->SetBinError(i+1,hist_alpha_error_final->GetBinContent(1)/alpha_nom[i]);

        //hist_resol->SetBinContent(i+1,hist_alpha_error_final->GetBinContent(1)/alpha_nom[i]);
    }


    for (int i=0;i<nstep;i++)
    {
        int nmPMT = nmPMT_min+i*step;
        TFile* f = new TFile(Form("diffuser4_400nm_nominal_mPMT_%i_%s_ledprofile.root",nmPMT,config.c_str()),"OPEN");
        TH1D* hist_alpha_result = (TH1D*)f->Get("hist_alpha_result");
        TH1D* hist_alpha_error_final = (TH1D*)f->Get("hist_alpha_error_final");
        if (!hist_alpha_result) continue;
        hist_resol->SetBinContent(i+1,fabs(hist_alpha_result->GetBinContent(1)/alpha_nom[i]-1));
        //hist_resol->SetBinError(i+1,hist_alpha_error_final->GetBinContent(1)/alpha_nom[i]);
    }

    hist_resol->GetXaxis()->SetTitle("#mPMT");
    hist_resol->GetYaxis()->SetTitle("Bias in L_{#alpha}");
    hist_resol->GetYaxis()->SetTitleOffset(1.5);
    hist_resol->GetYaxis()->SetRangeUser(0,0.1);
    
    hist_resol->SetLineColor(counter+1);
    hist_resol->SetLineStyle(counter+1);
    hist_resol->SetLineWidth(3);
    counter++;

    return hist_resol;
}

void plot_source_bias()
{
    gStyle->SetOptStat(0);
    
    counter = 0;
    
    TH1D* h_full = plot("fullring");
    TH1D* h_half = plot("halfring");
    TH1D* h_noout = plot("noring1");
    TH1D* h_nomid = plot("noring2");

    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    h_full->Draw("same");
    h_half->Draw("same");
    h_noout->Draw("same");
    h_nomid->Draw("same");
    h_full->Fit("pol1","W");
    h_half->Fit("pol1","W");
    h_noout->Fit("pol1","W");
    h_nomid->Fit("pol1","W");
    TLegend* legend = new TLegend(0.45,0.65,0.75,0.9);
    legend->AddEntry(h_full,"Full ring","l");
    legend->AddEntry(h_half,"Half ring","l");
    legend->AddEntry(h_noout,"No outer","l");
    legend->AddEntry(h_nomid,"No middle","l");
    legend->Draw("same");
    c1->SaveAs("Lalpha_source_bias.pdf");
}