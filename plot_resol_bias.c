void plot(std::string config)
{
    gStyle->SetOptStat(0);

    //std::string config = "fullring";
    //std::string config = "noring1";

    int nmPMT_min = 1000;
    int step = 100;
    const int nstep = 31;
    TH1D* hist_resol_bias = new TH1D(config.c_str(),config.c_str(),nstep,nmPMT_min-step/2,nmPMT_min+(nstep-0.5)*step);
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
        alpha_nom[i] = hist_alpha_result->GetBinContent(1);//10800;//
        hist_resol_bias->SetBinContent(i+1,hist_alpha_result->GetBinContent(1)/alpha_nom[i]);
        hist_resol_bias->SetBinError(i+1,hist_alpha_error_final->GetBinContent(1)/alpha_nom[i]);
    }
    for (int i=0;i<nstep;i++)
    {
        int nmPMT = nmPMT_min+i*step;
        TFile* f = new TFile(Form("diffuser4_400nm_nominal_mPMT_%i_%s_1ns.root",nmPMT,config.c_str()),"OPEN");
        TH1D* hist_alpha_result = (TH1D*)f->Get("hist_alpha_result");
        TH1D* hist_alpha_error_final = (TH1D*)f->Get("hist_alpha_error_final");
        hist_resol_bias_1ns->SetBinContent(i+1,hist_alpha_result->GetBinContent(1)/alpha_nom[i]);
        hist_resol_bias_1ns->SetBinError(i+1,hist_alpha_error_final->GetBinContent(1)/alpha_nom[i]);
    }
    for (int i=0;i<nstep;i++)
    {
        int nmPMT = nmPMT_min+i*step;
        TFile* f = new TFile(Form("diffuser4_400nm_nominal_mPMT_%i_%s_led_bias5pc.root",nmPMT,config.c_str()),"OPEN");
        TH1D* hist_alpha_result = (TH1D*)f->Get("hist_alpha_result");
        TH1D* hist_alpha_error_final = (TH1D*)f->Get("hist_alpha_error_final");
        if (!hist_alpha_result) continue;
        hist_resol_bias_led->SetBinContent(i+1,hist_alpha_result->GetBinContent(1)/alpha_nom[i]);
        hist_resol_bias_led->SetBinError(i+1,hist_alpha_error_final->GetBinContent(1)/alpha_nom[i]);
    }
    TCanvas* c1 = new TCanvas();
    hist_resol_bias->GetXaxis()->SetTitle("# mPMTs");
    hist_resol_bias->GetYaxis()->SetTitle("Ratio to nominal");
    hist_resol_bias->GetYaxis()->SetRangeUser(0.9,1.1);
    hist_resol_bias->Draw();
    hist_resol_bias_1ns->SetLineColor(kRed);
    hist_resol_bias_1ns->Draw("same");
    hist_resol_bias_led->SetLineColor(kViolet);
    hist_resol_bias_led->Draw("same");
    hist_resol_bias_led->Fit("pol1");
    c1->SetGridy();
    hist_globalcc->Draw();
    c1->SaveAs(Form("resol_bias_%s.pdf",config.c_str()));

    gStyle->SetPaintTextFormat("4.1f");
    c1->SetGridy(false);
    TFile* f = new TFile(Form("diffuser4_400nm_nominal_mPMT_%i_%s.root",4000,config.c_str()),"OPEN");
    TMatrixDSym* res_cor_matrix =  (TMatrixDSym*)f -> Get("res_cor_matrix");
    res_cor_matrix->Draw("colz text");
    TH2D *h_postfit_cor = (TH2D*) c1->GetPrimitive("TMatrixDBase");
    h_postfit_cor->GetZaxis()->SetRangeUser(-1,1);
    h_postfit_cor->GetXaxis()->SetBinLabel(2,"L_{#alpha}");
    h_postfit_cor->GetXaxis()->SetBinLabel(3,"A_{mPMT}");
    h_postfit_cor->GetXaxis()->SetRangeUser(1,42);
    h_postfit_cor->GetYaxis()->SetBinLabel(2,"L_{#alpha}");
    h_postfit_cor->GetYaxis()->SetBinLabel(3,"A_{mPMT}");
    h_postfit_cor->GetYaxis()->SetRangeUser(1,42);
    TH1D* hist_alpha_cor = new TH1D("","",40,0,40);
    for (int i=1;i<=40;i++)
    {
        hist_alpha_cor->SetBinContent(i,-h_postfit_cor->GetBinContent(2,i+2));
    }
    hist_alpha_cor->Draw();
    c1->SaveAs(Form("%s_cor.pdf",config.c_str()));
}

void plot_resol_bias()
{
    plot("fullring");
    plot("noring1");
    plot("noring2");
    plot("halfring");
}