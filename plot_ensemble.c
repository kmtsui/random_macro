void make_ensemble(int pmt_type=0)
{
    gStyle->SetOptStat(1);

    double truth_alpha = 10648;
    double truth_beta = -0.15;
    TH1D* hist_alpha = new TH1D("","",100,truth_alpha*0.8,truth_alpha*1.2);
    TH1D* hist_slope = new TH1D("","",100,-0.5,0.5);
    TH1D* hist_chi2 = pmt_type ==0 ? new TH1D("","",100,0,50000) : new TH1D("","",100,0,50000);

    int nFiles = 1000;
    for (int i=0;i<nFiles;i++)
    {
        TFile* f = pmt_type ==0 ? new TFile(Form("BnL/diffuser4_20_1600bins/diffuser4_20_400nm_toy%i.root",i),"OPEN") : new TFile(Form("mPMT/diffuser4_20_400bins/diffuser4_20_400nm_toy%i.root",i),"OPEN");
        TH1D* hist = (TH1D*)f->Get("hist_alphaZ_result");
        hist_alpha->Fill(hist->GetBinContent(1));
        hist_slope->Fill(hist->GetBinContent(2));
        TVectorT<double>*	chi2_tuple_postfit = (TVectorT<double>*)f->Get("chi2_tuple_postfit");
        hist_chi2->Fill((*chi2_tuple_postfit)[3]);
        f->Close();
    }

    TCanvas* c1 = new TCanvas();
    hist_alpha->GetXaxis()->SetTitle("Best-fit L_{#alpha} (cm)");
    hist_alpha->GetYaxis()->SetTitle("Count");
    //hist_alpha->SetTitle(Form("#sigma_{offset}=%1.1fns,#sigma_{smear}=%1.1fns",offset,smear));
    hist_alpha->Draw();
    TLine* l = new TLine(truth_alpha,0,truth_alpha,hist_alpha->GetMaximum());
    l->SetLineStyle(2);
    l->SetLineWidth(3);
    l->SetLineColor(kRed);
    l->Draw("same");
    if (pmt_type==0) c1->SaveAs("BnL_alpha_ensemble_joint.pdf");
    else c1->SaveAs("mPMT_alpha_ensemble_joint.pdf");

    hist_slope->Draw();
    hist_slope->GetXaxis()->SetTitle("Best-fit #beta");
    hist_slope->GetYaxis()->SetTitle("Count");
    l = new TLine(truth_beta,0,truth_beta,hist_slope->GetMaximum());
    l->SetLineStyle(2);
    l->SetLineWidth(3);
    l->SetLineColor(kRed);
    l->Draw("same");
    if (pmt_type==0) c1->SaveAs("BnL_beta_ensemble_joint.pdf");
    else c1->SaveAs("mPMT_beta_ensemble_joint.pdf");

    TCanvas* c2 = new TCanvas();
    hist_slope->GetXaxis()->SetTitle("Post-fit #chi^{2}");
    hist_slope->GetYaxis()->SetTitle("Count");
    hist_chi2->Draw();
    if (pmt_type==0) c2->SaveAs("BnL_chi2_ensemble.pdf");
    else c2->SaveAs("mPMT_chi2_ensemble.pdf");
}

TH1D* make_mPMT_ensemble(std::string dirname)
{
    gStyle->SetOptStat(0);

    double truth_alpha = 10648;
    double truth_beta = -0.15;

    double start = 100, end = 2000, step = 100;
    TH1D* h_alpha = new TH1D("","",(end-start)/step+1,start-step/2,end+step/2);

    TCanvas* c1 = new TCanvas();
    int nFiles = dirname=="mPMT" ? 1000 : 100;
    for (int i=0;i<h_alpha->GetNbinsX();i++)
    {
        TH1D* hist_alpha = new TH1D("","",100,truth_alpha*0.8,truth_alpha*1.2);
        int nmpmts = start+i*step;

        for (int j=0;j<nFiles;j++)
        {
            TFile* f = new TFile(Form("%s/%imPMTs/diffuser4_400nm_toy%i.root",dirname.c_str(),nmpmts,j),"OPEN");
            if (!f->IsOpen()) continue;
            TH1D* hist = (TH1D*)f->Get("hist_alphaZ_result");
            if (dirname=="mPMT")
                if (nmpmts==800)
                    if (hist->GetBinContent(1)>12000) continue;
            hist_alpha->Fill(hist->GetBinContent(1));
            f->Close();
        }

        h_alpha->SetBinContent(i+1,hist_alpha->GetMean());
        h_alpha->SetBinError(i+1,hist_alpha->GetRMS());

        hist_alpha->Draw();
        c1->SaveAs(Form("mPMT_%i.pdf",nmpmts));
    }

    h_alpha->Draw();
    c1->SaveAs("mPMT_number_alpha.pdf");

    for (int i=1;i<=h_alpha->GetNbinsX();i++)
    {
        h_alpha->SetBinContent(i,h_alpha->GetBinError(i)/h_alpha->GetBinContent(i));
    }
    h_alpha->GetXaxis()->SetTitle("#mPMTs");
    h_alpha->GetYaxis()->SetTitle("Estimated total uncertainty in L_{#alpha}");
    h_alpha->GetYaxis()->SetTitleOffset(1.3);
    h_alpha->GetYaxis()->SetRangeUser(0,0.06);
    h_alpha->SetLineWidth(2);
    h_alpha->Draw("hist");
    c1->SaveAs("mPMT_number_alpha_uncertainty.pdf");

    return h_alpha;

}

void plot_ensemble()
{
    // make_ensemble(0);
    // make_ensemble(1);

    TH1D* h_full = make_mPMT_ensemble("mPMT");
    TH1D* h_noRing1 = make_mPMT_ensemble("noRing1");
    TH1D* h_noRing2 = make_mPMT_ensemble("noRing2");

    TCanvas* c1 = new TCanvas();
    h_full->SetLineColor(kBlue);
    h_noRing1->SetLineColor(kRed);
    h_noRing2->SetLineColor(kGreen);
    h_full->SetLineStyle(1);
    h_noRing1->SetLineStyle(2);
    h_noRing2->SetLineStyle(3);
    h_noRing1->Draw("hist");
    h_noRing2->Draw("hist same");
    h_full->Draw("hist same");
    TLegend* l = new TLegend(0.5,0.5,0.8,0.8);
    l->AddEntry(h_full,"Full","l");
    l->AddEntry(h_noRing2,"No middle","l");
    l->AddEntry(h_noRing1,"No outer","l");
    l->Draw();
    c1->SaveAs("mPMT_number_ring_alpha_uncertainty.pdf");
}