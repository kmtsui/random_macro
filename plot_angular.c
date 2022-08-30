void plot_angular()
{
    gStyle->SetOptStat(0);

    TFile* f = new TFile("test4.10.root");
    TH1D* hist_mPMT_angular_result = (TH1D*)f->Get("hist_mPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final = (TH1D*)f->Get("hist_mPMT_angular_error_final");
    TH1D* hist_angular = new TH1D("angular_response","angular_response",20,0.5,1);
    for (int i=1;i<=hist_mPMT_angular_result->GetNbinsX();i++)
    {
        hist_angular->SetBinContent(i,hist_mPMT_angular_result->GetBinContent(i));
        hist_angular->SetBinError(i,hist_mPMT_angular_error_final->GetBinContent(i));
    }
    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    hist_angular->GetXaxis()->SetTitle("cos#theta");
    hist_angular->Draw("E0");

    TCanvas* c2 = new TCanvas();
    TH1D* hist_LED_cosths_result = (TH1D*)f->Get("hist_LED_cosths_result");
    TH1D* hist_LED_cosths_error_final = (TH1D*)f->Get("hist_LED_cosths_error_final");
    TH1D* hist_LED_cosths_prior = (TH1D*)f->Get("hist_LED_cosths_prior");
    TH1D* hist_LED = new TH1D("source_intensity","source_intensity",480,-40,0);
    TH1D* hist_LED_prior = (TH1D*)hist_LED->Clone();
    for (int i=1;i<=hist_LED->GetNbinsX();i++)
    {
        hist_LED->SetBinContent(i,hist_LED_cosths_result->GetBinContent(i));
        hist_LED->SetBinError(i,hist_LED_cosths_error_final->GetBinContent(i));
        hist_LED_prior->SetBinContent(i,hist_LED_cosths_prior->GetBinContent(i));
    }
    hist_LED->GetXaxis()->SetTitle("#theta_{s} (deg)");
    hist_LED->SetLineColor(kRed);
    hist_LED->Draw("E0");
    hist_LED_prior->Draw("same");
}