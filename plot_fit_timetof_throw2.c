void plot_fit_timetof_throw()
{
    double truht_alpha = 10650;
    TH1D* hist_alpha = new TH1D("","",40,/*truht_alpha*0.94,truht_alpha*1.06*/10000,11000);

    int nFiles = 100;
    for (int i=1;i<=nFiles;i++)
    {
        TFile* f = new TFile(Form("time_shift_smear/offset_1.0_smear_1.0/fitoutput_diffuser4_400nm_nominal_mPMT_20angular_%i.root",i),"OPEN");
        //TFile* f = new TFile(Form("time_shift_smear/fitoutput_diffuser4_400nm_nominal_mPMT_220angular_%i.root",i),"OPEN");
        //TFile* f = new TFile(Form("timetof_throw/fitoutput_diffuser4_400nm_nominal_mPMT_220angular_scatterfit_%i.root",i),"OPEN");
        TH1D* hist = (TH1D*)f->Get("hist_alpha_result");
        hist_alpha->Fill(hist->GetBinContent(1));
        f->Close();
    }

    hist_alpha->GetXaxis()->SetTitle("Fitted #alpha (cm)");
    hist_alpha->GetYaxis()->SetTitle("Count");
    hist_alpha->Draw();
}