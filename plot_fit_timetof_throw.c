void make_plot(double offset, double smear)
{
    double truht_alpha = 10650;
    TH1D* hist_alpha = new TH1D("","",40,/*truht_alpha*0.94,truht_alpha*1.06*/truht_alpha*0.9,truht_alpha*1.1);

    int nFiles = 100;
    for (int i=1;i<=nFiles;i++)
    {
        //TFile* f = new TFile(Form("time_throw/offset_%1.1f_smear_%1.1f/fitoutput_diffuser4_400nm_nominal_mPMT_%i.root",offset,smear,i),"OPEN");
        //TFile* f = new TFile(Form("time_throw/BnL_offset_%1.1f_smear_%1.1f/fitoutput_diffuser4_350nm_nominal_BnL_%i.root",offset,smear,i),"OPEN");
        TFile* f = new TFile(Form("time_throw/BnL_offset_%1.1f_smear_%1.1f/fitoutput_diffuser4_400nm_nominal_BnL_%i.root",offset,smear,i),"OPEN");
        //TFile* f = new TFile(Form("time_shift_smear/fitoutput_diffuser4_400nm_nominal_mPMT_220angular_%i.root",i),"OPEN");
        //TFile* f = new TFile(Form("timetof_throw/fitoutput_diffuser4_400nm_nominal_mPMT_220angular_scatterfit_%i.root",i),"OPEN");
        TH1D* hist = (TH1D*)f->Get("hist_alpha_result");
        hist_alpha->Fill(hist->GetBinContent(1));
        f->Close();
    }

    TCanvas* c1 = new TCanvas();
    hist_alpha->GetXaxis()->SetTitle("Best-fit L_{#alpha} (cm)");
    hist_alpha->GetYaxis()->SetTitle("Count");
    hist_alpha->SetTitle(Form("#sigma_{offset}=%1.1fns,#sigma_{smear}=%1.1fns",offset,smear));
    hist_alpha->Draw();
    c1->SaveAs(Form("BnL_alpha_offset_%1.1f_smearing_%1.1f.pdf",offset,smear));
}

void plot_fit_timetof_throw()
{
    // make_plot(0.2,0.2);
    // make_plot(0.2,0.5);
    // make_plot(0.2,1.0);
    // make_plot(0.5,0.2);
    // make_plot(0.5,0.5);
    // make_plot(0.5,1.0);
    // make_plot(1.0,0.2);
    // make_plot(1.0,0.5);
    make_plot(1.0,1.0);
    make_plot(2.0,2.0);
    make_plot(3.0,3.0);
    make_plot(4.0,4.0);
    make_plot(5.0,5.0);
}