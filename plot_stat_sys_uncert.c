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
        hist_resol->SetBinContent(i+1,hist_alpha_result->GetBinContent(1)/alpha_nom[i]-1);
        hist_resol->SetBinError(i+1,hist_alpha_error_final->GetBinContent(1)/alpha_nom[i]);
    }

    // for (int i=nstep;i>0;i--)
    // {
    //     double val;
    //     if (i==nstep) val = fabs(hist_resol->GetBinContent(i));
    //     else val = std::max(fabs(hist_resol->GetBinContent(i)),fabs(hist_resol->GetBinContent(i+1)));
    //     //val = std::max(val, hist_resol->GetBinError(i));
    //     hist_resol->SetBinContent(i,val);
    // }
    double prev = 0;
    int np = 0;
    TGraph* g = new TGraph(2);
    g->SetPoint(0,hist_resol->GetBinCenter(1),sys_lo);
    g->SetPoint(1,hist_resol->GetBinCenter(nstep),sys_hi);
    // for (int i=nstep;i>0;i--)
    // {
    //     //double val = hist_resol->GetBinContent(i)*hist_resol->GetBinContent(i) + hist_resol->GetBinError(i)*hist_resol->GetBinError(i);
    //     double val = std::max(fabs(hist_resol->GetBinContent(i)), fabs(hist_resol->GetBinError(i)));
    //     //val = sqrt(val);
    //     if (val>prev)
    //     {
    //         g->SetPoint(np,hist_resol->GetBinCenter(i),val);
    //         np++;
    //         prev = val;
    //     }
    //     else if (i==1)
    //         g->SetPoint(np,hist_resol->GetBinCenter(i),prev);
    //     //hist_resol->SetBinContent(i,val);
    // }

    // eff_var
    TFile* feff = new TFile("../mPMT_eff_study/10pc_err.root");
    TH1D* hist_eff_var = (TH1D*)feff->Get(config.c_str());

    for (int i=nstep;i>0;i--)
    {
        double val = g->Eval(hist_resol->GetBinCenter(i),0,"");
        val*=val;
        val+=hist_resol->GetBinError(i)*hist_resol->GetBinError(i) + 0.01*0.01;
        val+=hist_eff_var->GetBinContent(i)*hist_eff_var->GetBinContent(i);
        //val = 0.01*0.01;
        val = sqrt(val);
        hist_resol->SetBinContent(i,val);
    }

    hist_resol->GetXaxis()->SetTitle("#mPMT");
    hist_resol->GetYaxis()->SetTitle("Estimated stat. + sys. uncert. in L_{#alpha}");
    //hist_resol->GetYaxis()->SetTitle("Sys. uncert. in L_{#alpha} from late light");
    hist_resol->GetYaxis()->SetTitleOffset(1.5);
    hist_resol->GetYaxis()->SetRangeUser(0,0.12);
    
    hist_resol->SetLineColor(counter+1);
    hist_resol->SetLineStyle(counter+1);
    hist_resol->SetLineWidth(3);
    counter++;

    return hist_resol;
}

void plot_stat_sys_uncert()
{
    gStyle->SetOptStat(0);
    
    counter = 0;
    
    // TH1D* h_full = plot("fullring",0.017,0.009);
    // TH1D* h_half = plot("halfring",0.024,0.011);
    // TH1D* h_noout = plot("noring1",0.046,0.012);
    // TH1D* h_nomid = plot("noring2",0.024,0.010);
    TH1D* h_full = plot("fullring",0.015,0.005);
    TH1D* h_half = plot("halfring",0.02,0.0075);
    TH1D* h_noout = plot("noring1",0.04,0.012);
    TH1D* h_nomid = plot("noring2",0.016,0.006);

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
    c1->SaveAs("Lalpha_stat_sys_uncert_eff_var_10pc.pdf");
}