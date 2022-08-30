int counter;

TH1D* get_plot(std::string fname)
{
    TFile* f = new TFile(fname.c_str(),"OPEN");
    TMatrixDSym* res_cor_matrix =  (TMatrixDSym*)f -> Get("res_cor_matrix");
    int nBins = 40;
    TH1D* hist_alpha_cor = new TH1D("","",nBins,0,1);
    for (int i=1;i<=nBins;i++)
    {
        hist_alpha_cor->SetBinContent(i,-(*res_cor_matrix)[1][i+1]);
    }
    hist_alpha_cor->GetXaxis()->SetTitle("cos#theta");
    hist_alpha_cor->GetYaxis()->SetTitle("|Cor(L_{#alpha},A_{#theta})|");
    hist_alpha_cor->GetYaxis()->SetRangeUser(0,1);
    
    hist_alpha_cor->SetLineColor(counter+1);
    hist_alpha_cor->SetLineStyle(counter+1);
    hist_alpha_cor->SetLineWidth(2);
    counter++;

    return hist_alpha_cor;
}

void plot_correlation()
{
    gStyle->SetOptStat(0);
    
    counter = 0;
    
    TH1D* hist_fullring = get_plot("diffuser4_400nm_nominal_mPMT_4000_fullring.root");
    TH1D* hist_halfring = get_plot("diffuser4_400nm_nominal_mPMT_4000_halfring.root");
    TH1D* hist_noring1 = get_plot("diffuser4_400nm_nominal_mPMT_4000_noring1.root");
    TH1D* hist_noring2 = get_plot("diffuser4_400nm_nominal_mPMT_4000_noring2.root");

    TCanvas* c1 = new TCanvas();
    hist_fullring->Draw("same");
    hist_halfring->Draw("same");
    hist_noring1->Draw("same");
    hist_noring2->Draw("same");
    TLegend* legend = new TLegend(0.15,0.6,0.45,0.85);
    legend->AddEntry(hist_fullring,"Full ring","l");
    legend->AddEntry(hist_halfring,"Half ring","l");
    legend->AddEntry(hist_noring1,"No outer","l");
    legend->AddEntry(hist_noring2,"No middle","l");
    legend->Draw("same");
    c1->SaveAs("Lalpha_cor.pdf");

}