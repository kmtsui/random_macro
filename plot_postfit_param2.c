TH1D* plot_mPMT(std::string filename, double truth_alpha){
    gStyle->SetOptStat(0);

    TFile* f = new TFile(filename.c_str());
    TH1D* hist_mPMT_angular_result = (TH1D*)f->Get("hist_mPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final = (TH1D*)f->Get("hist_mPMT_angular_error_final");
    TH1D* hist_angular = new TH1D("","",40,0.,1);
    for (int i=1;i<=hist_mPMT_angular_result->GetNbinsX();i++)
    {
        hist_angular->SetBinContent(i,hist_mPMT_angular_result->GetBinContent(i));
        hist_angular->SetBinError(i,hist_mPMT_angular_error_final->GetBinContent(i));
    }
    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX());
    int ndof = ((TTree*)f->Get("PMTTree"))->GetEntries()-41;
    double chi2perdof = chi2/ndof;
    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    hist_angular->GetXaxis()->SetTitle("cos#theta");
    hist_angular->GetYaxis()->SetTitle("A_{mPMT}");
    hist_angular->Draw("E0");
    double fitted_alpha = ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);
    TLatex latex;
    latex.SetTextSize(0.045);
    double miny = hist_angular->GetMinimum();
    double maxy = hist_angular->GetMaximum();
    double diff = (maxy-miny)/10;
    latex.DrawLatex(0.55,miny+diff*3,Form("#chi^{2}/dof = %4.0f/%i",chi2,ndof));
    latex.DrawLatex(0.55,miny+diff*2,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha));
    latex.DrawLatex(0.55,miny+diff*1,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha));
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    c1->SetGridy(false);
    c1->SetCanvasSize(3000, 1500);
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
    figname += "_cor";
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    return hist_angular;
}

TH1D* plot_mPMT_source(std::string filename, double truth_alpha){
    gStyle->SetOptStat(0);

    TFile* f = new TFile(filename.c_str());
    TH1D* hist_mPMT_angular_result = (TH1D*)f->Get("hist_mPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final = (TH1D*)f->Get("hist_mPMT_angular_error_final");
    TH1D* hist_angular = new TH1D("","",20,0.5,1);
    if (hist_mPMT_angular_result)
        for (int i=1;i<=hist_mPMT_angular_result->GetNbinsX();i++)
        {
            hist_angular->SetBinContent(i,hist_mPMT_angular_result->GetBinContent(i));
            hist_angular->SetBinError(i,hist_mPMT_angular_error_final->GetBinContent(i));
        }
    int ndof = ((TTree*)f->Get("PMTTree"))->GetEntries()-21;
    TH1D* hist_LED_cosths_result = (TH1D*)f->Get("hist_LED_cosths_result");
    TH1D* hist_LED_cosths_prior = (TH1D*)f->Get("hist_LED_cosths_prior");
    TH1D* hist_LED_cosths_error_final = (TH1D*)f->Get("hist_LED_cosths_error_final");
    TH1D* hist_LED_cosths_error_prior = (TH1D*)f->Get("hist_LED_cosths_error_prior");
    TH1D* hist_LED_cosths_prefit = new TH1D("","",480,-40,0);
    TH1D* hist_LED_cosths_postfit = new TH1D("","",480,-40,0);
    int nledpar = 0;
    for (int i=1;i<=hist_LED_cosths_prefit->GetNbinsX();i++)
    {
        hist_LED_cosths_prefit->SetBinContent(i,hist_LED_cosths_prior->GetBinContent(i));
        hist_LED_cosths_prefit->SetBinError(i,hist_LED_cosths_error_prior->GetBinContent(i));
        hist_LED_cosths_postfit->SetBinContent(i,hist_LED_cosths_result->GetBinContent(i));
        hist_LED_cosths_postfit->SetBinError(i,hist_LED_cosths_error_final->GetBinContent(i));
        if (hist_LED_cosths_error_final->GetBinContent(i)>1.e-6) 
        {
            ndof++;
            nledpar++;
        }
    }
    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX()-1);
    double chi2perdof = chi2/ndof;
    TH1D* chi2_syst_periter = (TH1D*)f->Get("chi2_syst_periter");
    double chi2syst = chi2_syst_periter->GetBinContent(chi2_syst_periter->GetNbinsX()-1);
    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    hist_angular->GetXaxis()->SetTitle("cos#theta");
    hist_angular->GetYaxis()->SetTitle("A_{mPMT}");
    hist_angular->Draw("E0");
    double fitted_alpha = ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);
    TLatex latex;
    latex.SetTextSize(0.045);
    double miny = hist_angular->GetMinimum()*1.1;
    double maxy = hist_angular->GetMaximum();
    double diff = (maxy-miny)/10;
    latex.DrawLatex(0.75,miny+diff*2,Form("#chi^{2}/dof = %4.0f/%i",chi2,ndof));
    latex.DrawLatex(0.75,miny+diff,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha));
    latex.DrawLatex(0.75,miny,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha));
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    c1->SetGridy(false);
    TMatrixDSym* res_cor_matrix =  (TMatrixDSym*)f -> Get("res_cor_matrix");
    res_cor_matrix->Draw("colz text");
    TH2D *h_postfit_cor = (TH2D*) c1->GetPrimitive("TMatrixDBase");
    h_postfit_cor->GetZaxis()->SetRangeUser(-1,1);
    h_postfit_cor->GetXaxis()->SetBinLabel(2,"L_{#alpha}");
    h_postfit_cor->GetXaxis()->SetBinLabel(3,"A_{mPMT}");
    h_postfit_cor->GetXaxis()->SetRangeUser(1,22);
    h_postfit_cor->GetYaxis()->SetBinLabel(2,"L_{#alpha}");
    h_postfit_cor->GetYaxis()->SetBinLabel(3,"A_{mPMT}");
    h_postfit_cor->GetYaxis()->SetRangeUser(1,22);
    //figname += "_cor";
    c1->SaveAs(Form("%s_cor.pdf",figname.c_str()));

    hist_LED_cosths_prefit->GetXaxis()->SetTitle("-#theta_{s} (deg)");
    hist_LED_cosths_prefit->GetYaxis()->SetTitle("I(#theta_{s},#phi_{s})");
    hist_LED_cosths_prefit->SetLineColor(kRed);
    hist_LED_cosths_prefit->Draw("E0");
    hist_LED_cosths_postfit->SetLineColor(kBlue);
    hist_LED_cosths_postfit->Draw("E0 same");
    TLegend* legend = new TLegend(0.5,0.2,0.8,0.4);
    legend->AddEntry(hist_LED_cosths_prefit,"Pre-fit","lep");
    legend->AddEntry(hist_LED_cosths_postfit,"Post-fit","lep");
    legend->AddEntry((TObject*)0,Form("#chi^{2}_{syst}/npar=%2.0f/%i",chi2syst,nledpar),"");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs(Form("%s_intensity.pdf",figname.c_str()));

    TH1D* hist_LED_cosths_residual = (TH1D*)hist_LED_cosths_postfit->Clone();
    TH1D* hist_LED_cosths_residual_distribution = new TH1D("","",100,-0.1,0.1);
    for (int i=1;i<=hist_LED_cosths_residual->GetNbinsX();i++)
    {
        double val = (hist_LED_cosths_postfit->GetBinContent(i)-hist_LED_cosths_prefit->GetBinContent(i))/hist_LED_cosths_prefit->GetBinContent(i);
        double err = hist_LED_cosths_postfit->GetBinError(i)/hist_LED_cosths_prefit->GetBinContent(i);
        hist_LED_cosths_residual->SetBinContent(i,val);
        hist_LED_cosths_residual->SetBinError(i,err);
        if (err>0)
            hist_LED_cosths_residual_distribution->Fill(val);
    }
    hist_LED_cosths_residual->Draw();
    c1->SaveAs(Form("%s_intensity_residual.pdf",figname.c_str()));
    hist_LED_cosths_residual_distribution->Draw();
    c1->SaveAs(Form("%s_intensity_residual_distribution.pdf",figname.c_str()));

    return hist_angular;
}

void plot_mPMT_3angular(std::string filename, double truth_alpha){
    gStyle->SetOptStat(0);

    TFile* f = new TFile(filename.c_str());
    TH1D* hist_mPMT_angular_result = (TH1D*)f->Get("hist_mPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final = (TH1D*)f->Get("hist_mPMT_angular_error_final");
    TH1D* hist_angular = new TH1D("","",40,0.5,1);
    for (int i=1;i<=hist_mPMT_angular_result->GetNbinsX();i++)
    {
        hist_angular->SetBinContent(i,hist_mPMT_angular_result->GetBinContent(i));
        hist_angular->SetBinError(i,hist_mPMT_angular_error_final->GetBinContent(i));
    }
    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX()-1);
    int ndof = ((TTree*)f->Get("PMTTree"))->GetEntries()-41-12;
    double chi2perdof = chi2/ndof;
    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    hist_angular->GetXaxis()->SetTitle("cos#theta");
    hist_angular->GetYaxis()->SetTitle("A_{#theta}");
    hist_angular->Draw("E0");
    double fitted_alpha = ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);
    TLatex latex;
    latex.SetTextSize(0.045);
    double miny = hist_angular->GetMinimum()*1.2;
    latex.DrawLatex(0.75,miny+0.02,Form("#chi^{2}/dof = %4.0f/%i",chi2,ndof));
    latex.DrawLatex(0.75,miny+0.01,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha));
    latex.DrawLatex(0.75,miny,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha));
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    std::string anotherfilename = filename;
    anotherfilename.erase(anotherfilename.find("_3angular"),9);
    TH1D* another_angular = plot_mPMT(anotherfilename,truth_alpha);
    c1->cd();
    another_angular->SetLineColor(kRed);
    another_angular->Draw("E0 same");
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    TH1D* hist_phim_result = (TH1D*)f->Get("hist_phim_result");
    TH1D* hist_phim_error_final = (TH1D*)f->Get("hist_phim_error_final");
    TH1D* hist_mPMT_id_result = (TH1D*)f->Get("hist_mPMT_id_result");
    TH1D* hist_mPMT_id_error_final = (TH1D*)f->Get("hist_mPMT_id_error_final");
    double phim_bins[] = {0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0,TMath::Pi()};
    TH1D* hist_phim = new TH1D("","",11,phim_bins);
    for (int i=1;i<=11;i++)
    {
        hist_phim->SetBinContent(i,hist_phim_result->GetBinContent(i));
        hist_phim->SetBinError(i,hist_phim_error_final->GetBinContent(i));
    }
    hist_phim->GetXaxis()->SetTitle("#phi_{m}");
    hist_phim->GetYaxis()->SetTitle("A_{#phi_{m}}");
    hist_phim->GetYaxis()->SetRangeUser(0.99,1.01);
    hist_phim->Draw("E0");
    latex.DrawLatex(1.5,0.996,Form("A_{Ring}(1) = 1 (fixed)"));
    latex.DrawLatex(1.5,0.995,Form("A_{Ring}(2) = %4.4f +/- %4.4f",hist_mPMT_id_result->GetBinContent(2),hist_mPMT_id_error_final->GetBinContent(2)));
    latex.DrawLatex(1.5,0.994,Form("A_{Ring}(3) = %4.4f +/- %4.4f",hist_mPMT_id_result->GetBinContent(3),hist_mPMT_id_error_final->GetBinContent(3)));
    c1->SaveAs(Form("%s_phim.pdf",figname.c_str()));

}

TH1D* plot_BnLPMT(std::string filename, double truth_alpha){
    gStyle->SetOptStat(0);

    TFile* f = new TFile(filename.c_str());
    TH1D* hist_mPMT_angular_result = (TH1D*)f->Get("hist_BnLPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final = (TH1D*)f->Get("hist_BnLPMT_angular_error_final");
    TH1D* hist_angular = new TH1D("","",40,0.,1);
    int nPar = 1;
    double miny = 1.; double maxy=0;
    for (int i=1;i<=hist_mPMT_angular_result->GetNbinsX();i++)
    {
        if (hist_mPMT_angular_error_final->GetBinContent(i)<1.e-7) continue;
        hist_angular->SetBinContent(i,hist_mPMT_angular_result->GetBinContent(i));
        hist_angular->SetBinError(i,hist_mPMT_angular_error_final->GetBinContent(i));
        nPar++;
        if (miny>hist_mPMT_angular_result->GetBinContent(i)) miny = hist_mPMT_angular_result->GetBinContent(i);
        if (maxy<hist_mPMT_angular_result->GetBinContent(i)) maxy = hist_mPMT_angular_result->GetBinContent(i);
    }
    hist_angular->GetYaxis()->SetRangeUser(miny*0.95,maxy*1.05);
    miny *= 1.1;
    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX());
    int ndof = ((TTree*)f->Get("PMTTree"))->GetEntries()-nPar;
    double chi2perdof = chi2/ndof;
    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    hist_angular->GetXaxis()->SetTitle("cos#theta");
    hist_angular->GetYaxis()->SetTitle("A_{BnL}");
    hist_angular->GetYaxis()->SetTitleOffset(1.1);
    hist_angular->Draw("E0");
    double fitted_alpha = ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);
    TLatex latex;
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.55,miny+0.02,Form("#chi^{2}/dof = %4.0f/%i",chi2,ndof));
    latex.DrawLatex(0.55,miny+0.01,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha));
    latex.DrawLatex(0.55,miny,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha));
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    c1->SetGridy(false);
    c1->SetCanvasSize(3000, 1500);
    TMatrixDSym* res_cor_matrix =  (TMatrixDSym*)f -> Get("res_cor_matrix");
    res_cor_matrix->Draw("colz text");
    TH2D *h_postfit_cor = (TH2D*) c1->GetPrimitive("TMatrixDBase");
    h_postfit_cor->GetZaxis()->SetRangeUser(-1,1);
    h_postfit_cor->GetXaxis()->SetBinLabel(2,"L_{#alpha}");
    h_postfit_cor->GetXaxis()->SetBinLabel(3,"A_{BnL}");
    h_postfit_cor->GetXaxis()->SetRangeUser(1,42);
    h_postfit_cor->GetYaxis()->SetBinLabel(2,"L_{#alpha}");
    h_postfit_cor->GetYaxis()->SetBinLabel(3,"A_{BnL}");
    h_postfit_cor->GetYaxis()->SetRangeUser(1,42);
    figname += "_cor";
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    return hist_angular;
}

TH1D* plot_BnLPMT_source(std::string filename, double truth_alpha){
    gStyle->SetOptStat(0);

    TFile* f = new TFile(filename.c_str());
    TH1D* hist_mPMT_angular_result = (TH1D*)f->Get("hist_BnLPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final = (TH1D*)f->Get("hist_BnLPMT_angular_error_final");
    TH1D* hist_angular = new TH1D("","",20,0.5,1);
    int nPar = 1;
    double miny = 1.; double maxy=0;
    for (int i=1;i<=hist_mPMT_angular_result->GetNbinsX();i++)
    {
        if (hist_mPMT_angular_error_final->GetBinContent(i)<1.e-7) continue;
        hist_angular->SetBinContent(i,hist_mPMT_angular_result->GetBinContent(i));
        hist_angular->SetBinError(i,hist_mPMT_angular_error_final->GetBinContent(i));
        nPar++;
        if (miny>hist_mPMT_angular_result->GetBinContent(i)) miny = hist_mPMT_angular_result->GetBinContent(i);
        if (maxy<hist_mPMT_angular_result->GetBinContent(i)) maxy = hist_mPMT_angular_result->GetBinContent(i);
    }
    hist_angular->GetYaxis()->SetRangeUser(miny*0.95,maxy*1.05);
    miny *= 1.1;
    int ndof = ((TTree*)f->Get("PMTTree"))->GetEntries()-nPar;
    TH1D* hist_LED_cosths_result = (TH1D*)f->Get("hist_LED_cosths_result");
    TH1D* hist_LED_cosths_prior = (TH1D*)f->Get("hist_LED_cosths_prior");
    TH1D* hist_LED_cosths_error_final = (TH1D*)f->Get("hist_LED_cosths_error_final");
    TH1D* hist_LED_cosths_error_prior = (TH1D*)f->Get("hist_LED_cosths_error_prior");
    TH1D* hist_LED_cosths_prefit = new TH1D("","",480,-40,0);
    TH1D* hist_LED_cosths_postfit = new TH1D("","",480,-40,0);
    for (int i=1;i<=hist_LED_cosths_prefit->GetNbinsX();i++)
    {
        hist_LED_cosths_prefit->SetBinContent(i,hist_LED_cosths_prior->GetBinContent(i));
        hist_LED_cosths_prefit->SetBinError(i,hist_LED_cosths_error_prior->GetBinContent(i));
        hist_LED_cosths_postfit->SetBinContent(i,hist_LED_cosths_result->GetBinContent(i));
        hist_LED_cosths_postfit->SetBinError(i,hist_LED_cosths_error_final->GetBinContent(i));
        if (hist_LED_cosths_error_final->GetBinContent(i)>1.e-6) ndof++;
    }
    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX()-1);
    double chi2perdof = chi2/ndof;
    TH1D* chi2_syst_periter = (TH1D*)f->Get("chi2_syst_periter");
    double chi2syst = chi2_syst_periter->GetBinContent(chi2_syst_periter->GetNbinsX()-1);
    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    hist_angular->GetXaxis()->SetTitle("cos#theta");
    hist_angular->GetYaxis()->SetTitle("A_{BnL}");
    hist_angular->Draw("E0");
    double fitted_alpha = ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);
    TLatex latex;
    latex.SetTextSize(0.045);
    double diff = (maxy-miny)/10;
    latex.DrawLatex(0.75,miny+diff*2,Form("#chi^{2}/dof = %4.0f/%i",chi2,ndof));
    latex.DrawLatex(0.75,miny+diff,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha));
    latex.DrawLatex(0.75,miny,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha));
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    c1->SetGridy(false);
    TMatrixDSym* res_cor_matrix =  (TMatrixDSym*)f -> Get("res_cor_matrix");
    res_cor_matrix->Draw("colz text");
    TH2D *h_postfit_cor = (TH2D*) c1->GetPrimitive("TMatrixDBase");
    h_postfit_cor->GetZaxis()->SetRangeUser(-1,1);
    h_postfit_cor->GetXaxis()->SetBinLabel(2,"L_{#alpha}");
    h_postfit_cor->GetXaxis()->SetBinLabel(3,"A_{BnL}");
    h_postfit_cor->GetXaxis()->SetRangeUser(1,22);
    h_postfit_cor->GetYaxis()->SetBinLabel(2,"L_{#alpha}");
    h_postfit_cor->GetYaxis()->SetBinLabel(3,"A_{BnL}");
    h_postfit_cor->GetYaxis()->SetRangeUser(1,22);
    //figname += "_cor";
    c1->SaveAs(Form("%s_cor.pdf",figname.c_str()));

    hist_LED_cosths_prefit->GetXaxis()->SetTitle("-#theta_{s} (deg)");
    hist_LED_cosths_prefit->GetYaxis()->SetTitle("I(#theta_{s},#phi_{s})");
    hist_LED_cosths_prefit->SetLineColor(kRed);
    hist_LED_cosths_prefit->Draw("E0");
    hist_LED_cosths_postfit->SetLineColor(kBlue);
    hist_LED_cosths_postfit->Draw("E0 same");
    TLegend* legend = new TLegend(0.5,0.2,0.8,0.4);
    legend->AddEntry(hist_LED_cosths_prefit,"Pre-fit","lep");
    legend->AddEntry(hist_LED_cosths_postfit,"Post-fit","lep");
    legend->AddEntry((TObject*)0,Form("#chi^{2}_{syst}=%2.0f",chi2syst),"");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs(Form("%s_intensity.pdf",figname.c_str()));

    return hist_angular;
}

void plot_postfit_combined(std::string filename, double truth_alpha)
{
    gStyle->SetOptStat(0);

    TFile* f = new TFile(filename.c_str());
    TH1D* hist_mPMT_angular_result = (TH1D*)f->Get("hist_mPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final = (TH1D*)f->Get("hist_mPMT_angular_error_final");
    TH1D* hist_mPMT_angular = new TH1D("","",40,0,1);
    int nPar = 1;
    double miny = 1.; double maxy=0;
    for (int i=1;i<=hist_mPMT_angular_result->GetNbinsX();i++)
    {
        if (hist_mPMT_angular_result->GetBinContent(i)<1.e-7) continue;
        hist_mPMT_angular->SetBinContent(i,hist_mPMT_angular_result->GetBinContent(i));
        hist_mPMT_angular->SetBinError(i,hist_mPMT_angular_error_final->GetBinContent(i));
        nPar++;
        if (miny>hist_mPMT_angular->GetBinContent(i)) miny = hist_mPMT_angular->GetBinContent(i);
        if (maxy<hist_mPMT_angular->GetBinContent(i)) maxy = hist_mPMT_angular->GetBinContent(i);
    }
    TH1D* hist_BnLPMT_angular_result = (TH1D*)f->Get("hist_BnLPMT_angular_result");
    TH1D* hist_BnLPMT_angular_error_final = (TH1D*)f->Get("hist_BnLPMT_angular_error_final");
    TH1D* hist_BnLPMT_angular = new TH1D("","",40,0,1);
    for (int i=1;i<=hist_BnLPMT_angular_result->GetNbinsX();i++)
    {
        if (hist_BnLPMT_angular_error_final->GetBinContent(i)<1.e-7) continue;
        hist_BnLPMT_angular->SetBinContent(i,hist_BnLPMT_angular_result->GetBinContent(i));
        hist_BnLPMT_angular->SetBinError(i,hist_BnLPMT_angular_error_final->GetBinContent(i));
        nPar++;
        if (miny>hist_BnLPMT_angular->GetBinContent(i)) miny = hist_BnLPMT_angular->GetBinContent(i);
        if (maxy<hist_BnLPMT_angular->GetBinContent(i)) maxy = hist_BnLPMT_angular->GetBinContent(i);
    }

    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX());
    int ndof = ((TTree*)f->Get("PMTTree"))->GetEntries()-nPar;
    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    hist_mPMT_angular->GetXaxis()->SetTitle("cos#theta");
    hist_mPMT_angular->GetYaxis()->SetTitle("A_{X}");
    hist_mPMT_angular->GetYaxis()->SetRangeUser(miny*0.95,maxy*1.05);
    hist_mPMT_angular->SetLineColor(kBlue);
    hist_mPMT_angular->Draw("E0");
    hist_BnLPMT_angular->SetLineColor(kRed);
    hist_BnLPMT_angular->Draw("E0 same");
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    TLegend* legend = new TLegend(0.2,0.65,0.6,0.95);
    legend->AddEntry(hist_mPMT_angular,"A_{mPMT}","lep");
    legend->AddEntry(hist_BnLPMT_angular,"A_{BnL}","lep");
    legend->AddEntry((TObject*)0,Form("#chi^{2}/dof = %4.0f/%i",chi2,ndof),"");
    double fitted_alpha = ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);
    legend->AddEntry((TObject*)0,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha),"");
    legend->AddEntry((TObject*)0,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha),"");
    legend->Draw();
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    gStyle->SetPaintTextFormat("4.1f");
    c1->SetGridy(false);
    //c1->SetCanvasSize(3000, 1500);
    TMatrixDSym* res_cor_matrix =  (TMatrixDSym*)f -> Get("res_cor_matrix");
    res_cor_matrix->Draw("colz");
    TH2D *h_postfit_cor = (TH2D*) c1->GetPrimitive("TMatrixDBase");
    h_postfit_cor->GetZaxis()->SetRangeUser(-1,1);
    h_postfit_cor->GetXaxis()->SetBinLabel(2,"L_{#alpha}");
    h_postfit_cor->GetXaxis()->SetBinLabel(3,"A_{BnL}");
    h_postfit_cor->GetXaxis()->SetBinLabel(43,"A_{mPMT}");
    h_postfit_cor->GetXaxis()->SetRangeUser(1,82);
    h_postfit_cor->GetYaxis()->SetBinLabel(2,"L_{#alpha}");
    h_postfit_cor->GetYaxis()->SetBinLabel(3,"A_{BnL}");
    h_postfit_cor->GetYaxis()->SetBinLabel(43,"A_{mPMT}");
    h_postfit_cor->GetYaxis()->SetRangeUser(1,82);
    figname += "_cor";
    c1->SaveAs(Form("%s.pdf",figname.c_str()));
}

void plot_postfit_param(){
    gStyle->SetPaintTextFormat("4.2f");

    //plot_mPMT_3angular("diffusr4_400nm_nominal_mPMT_3angular.root",10648);
    //plot_mPMT("diffusr4_400nm_nominal_mPMT_ledprofile.root",10648);
    //plot_mPMT("diffusr4_400nm_nominal_mPMT_ledprofile_thetas_fix.root",10648);
    //plot_mPMT("diffusr4_400nm_nominal_mPMT_ledprofile_thetas_phis_fix.root",10648);

    //plot_mPMT_source("diffusr4_400nm_nominal_mPMT_ledprofile_thetas_phis_1percent.root",10648);
    //plot_mPMT_source("diffusr4_400nm_nominal_mPMT_ledprofile_thetas_phis_5percent.root",10648);
    //plot_mPMT_source("diffusr4_400nm_nominal_mPMT_ledprofile_thetas_phis_10percent.root",10648);
    // plot_mPMT_source("diffusr4_400nm_nominal_mPMT_pol_ledprofile_fullring.root",10648);
    // plot_mPMT_source("diffusr4_400nm_nominal_mPMT_pol_ledprofile_noring1.root",10648);
    // plot_mPMT_source("diffusr4_400nm_nominal_mPMT_pol_ledprofile_noring2.root",10648);

    //plot_mPMT("diffuser4_400nm_nominal_mPMT_simplest.root",10648);
    //plot_mPMT("diffuser4_350nm_nominal_mPMT_simplest.root",6351);
    //plot_BnLPMT("diffuser4_400nm_nominal_BnL_simplest.root",10648);
    //plot_BnLPMT("diffuser4_350nm_nominal_BnL_simplest.root",6351);
    plot_postfit_combined("diffuser4_400nm_nominal_combined_simplest.root",10648);
}