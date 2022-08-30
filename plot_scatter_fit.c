double calc_scatter_length(double ray)
{
    return 13186.7/0.75*ray;
}

void plot_scatter_fit()
{
    const int nFiles = 5;

    std::string filename[] =
    {
        "fitoutput_diffuser_4_400nm_abs_0.3_ray_3.653.root",
        "fitoutput_diffuser_4_400nm_abs_0.8_ray_0.88132.root",
        "fitoutput_diffuser4_400nm_nominal_scatterfit.root",
        "fitoutput_diffuser_4_400nm_abs_1.8_ray_0.7034.root",
        "fitoutput_diffuser_4_400nm_abs_2.3_ray_0.6796.root",
    };

    double scatter_length[] = {calc_scatter_length(3.653),calc_scatter_length(0.88132),calc_scatter_length(0.75),calc_scatter_length(0.7034),calc_scatter_length(0.6796)};

    int nBins = 60;

    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas();
    TH1D* hist_scatter[nFiles];
    TLegend* legend = new TLegend(0.1,0.65,0.5,0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    for (int i=0;i<nFiles;i++)
    {
        hist_scatter[i] = new TH1D("","",nBins,0,nBins);
        TFile* f = new TFile(filename[i].c_str(),"OPEN");
        TH1D* hist_scatter_result = (TH1D*)f->Get("hist_scatter_result");
        TH1D* hist_scatter_error_final = (TH1D*)f->Get("hist_scatter_error_final");

        for (int j=1;j<=nBins;j++)
        {
            if (hist_scatter_error_final->GetBinContent(j)>0) 
            {
                hist_scatter[i]->SetBinContent(j,hist_scatter_result->GetBinContent(j));
                hist_scatter[i]->SetBinError(j,hist_scatter_error_final->GetBinContent(j));
            }
        }

        hist_scatter[i]->SetLineWidth(3);
        hist_scatter[i]->SetLineColor(i+1);
        if (i==4) hist_scatter[i]->SetLineColor(i+2);
        hist_scatter[i]->Draw("same");
        legend->AddEntry(hist_scatter[i],Form("Scattering length = %1.0f cm", scatter_length[i]));
    }
    legend->Draw("same");
}