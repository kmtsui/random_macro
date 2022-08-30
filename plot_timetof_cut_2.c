void plot_timetof_cut(bool pmt_type = 0)
{
    gStyle->SetOptStat(0);
    
    double truth_alpha = 10648;
    double start = pmt_type==0 ? 10 : 1.0;
    double stepsize = pmt_type==0 ? 1 : 0.1;
    int nStep = pmt_type==0 ? 11: 21;
    TH1D* hTimetof = new TH1D("","",nStep,start-0.5*stepsize,start+(nStep-0.5)*stepsize);

    for (int i=0;i<nStep;i++)
    {
        double cut = start+i*stepsize-950;
        TFile* f = pmt_type==0 ? new TFile(Form("diffuser4_400nm_nominal_BnL_timetof_%1.1f.root",cut)): new TFile(Form("diffuser4_400nm_nominal_mPMT_timetof_%1.1f.root",cut));
        TH1D* hist_alpha_result = (TH1D*)f->Get("hist_alpha_result");
        TH1D* hist_alpha_error_final = (TH1D*)f->Get("hist_alpha_error_final");
        hTimetof->SetBinContent(i+1,hist_alpha_result->GetBinContent(1));
        hTimetof->SetBinError(i+1,hist_alpha_error_final->GetBinContent(1));
        f->Close();
    }

    TCanvas* c1 = new TCanvas();
    hTimetof->GetXaxis()->SetTitle("t_{D} cut (ns)");
    hTimetof->GetYaxis()->SetTitle("Best-fit L_{#alpha} (cm)");
    hTimetof->GetYaxis()->SetTitleOffset(1.3);
    hTimetof->Draw("E0");
    TLine* l = new TLine(start-0.5*stepsize,truth_alpha,start+(nStep-0.5)*stepsize,truth_alpha);
    l->SetLineStyle(2);
    l->SetLineWidth(3);
    l->SetLineColor(kRed);
    l->Draw("same");
    if (pmt_type==0) c1->SaveAs("BnL_timetof_cut.pdf");
    else c1->SaveAs("mPMT_timetof_cut.pdf");
}