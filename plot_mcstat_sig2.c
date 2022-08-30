void plot_mcstat_sig2(std::string filename)
{
    TFile* f = new TFile(filename.c_str());
    TTree* t = (TTree*)f->Get("PMTTree");
    double indirectPEerr2, nPE_pred;
    t->SetBranchAddress("indirectPEerr2",&indirectPEerr2);
    t->SetBranchAddress("nPE_pred",&nPE_pred);

    TH1D* hSig2 = new TH1D("","",100,0,0.001);
    for (int i=0;i<t->GetEntries();i++)
    {
        t->GetEntry(i);
        hSig2->Fill(indirectPEerr2/nPE_pred/nPE_pred);
    }
    TCanvas* c1 = new TCanvas();
    hSig2->GetXaxis()->SetTitle("#sigma_{j}^{2}");
    hSig2->Draw();
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    c1->SaveAs(Form("%s_sig2.pdf",figname.c_str()));
}

void plot_mcstat_sig2()
{
    plot_mcstat_sig2("diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root");
    plot_mcstat_sig2("diffuser4_400nm_nominal_BnL_pol_indirect.root");
}