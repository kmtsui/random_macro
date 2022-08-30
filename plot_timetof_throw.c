void make_plot(int pmt_type=0)
{
    gStyle->SetOptStat(1);

    double truht_alpha = 10648;
    TH1D* hist_alpha = new TH1D("","",100,truht_alpha*0.97,truht_alpha*1.03);
    TH1D* hist_chi2 = pmt_type ==0 ? new TH1D("","",100,5000,15000) : new TH1D("","",100,10000,25000);

    int nFiles = 1000;
    for (int i=0;i<nFiles;i++)
    {
        TFile* f = pmt_type ==0 ? new TFile(Form("BnL/diffuser4_400nm_nominal_BnL_toy%i.root",i),"OPEN") : new TFile(Form("mPMT/diffuser4_400nm_nominal_mPMT_toy%i.root",i),"OPEN");
        TH1D* hist = (TH1D*)f->Get("hist_alpha_result");
        hist_alpha->Fill(hist->GetBinContent(1));
        TVectorT<double>*	chi2_tuple_postfit = (TVectorT<double>*)f->Get("chi2_tuple_postfit");
        hist_chi2->Fill((*chi2_tuple_postfit)[3]);
        f->Close();
    }

    TCanvas* c1 = new TCanvas();
    hist_alpha->GetXaxis()->SetTitle("Best-fit L_{#alpha} (cm)");
    hist_alpha->GetYaxis()->SetTitle("Count");
    //hist_alpha->SetTitle(Form("#sigma_{offset}=%1.1fns,#sigma_{smear}=%1.1fns",offset,smear));
    hist_alpha->Draw();
    TLine* l = new TLine(truht_alpha,0,truht_alpha,hist_alpha->GetMaximum());
    l->SetLineStyle(2);
    l->SetLineWidth(3);
    l->SetLineColor(kRed);
    l->Draw("same");
    if (pmt_type==0) c1->SaveAs("BnL_timetof_offset_throw.pdf");
    else c1->SaveAs("mPMT_timetof_offset_throw.pdf");

    TCanvas* c2 = new TCanvas();
    hist_chi2->Draw();
    if (pmt_type==0) c2->SaveAs("BnL_timetof_offset_throw_chi2.pdf");
    else c2->SaveAs("mPMT_timetof_offset_throw_chi2.pdf");
}

void plot_timetof_throw()
{
    make_plot(0);
    make_plot(1);
}