void make_zpos_eff_prior()
{
    TFile* f = new TFile("zpos_eff_input.root","RECREATE");
    int nBins = 40;
    double effVar = 0.02;
    TH1D* h_prior = new TH1D("","",nBins,0,nBins);
    TMatrixDSym cov_matrix(nBins);
    for (int i=1; i<=nBins; i++)
    {
        double val = 1. - effVar + 2.*effVar*(i-1.)/(nBins-1.);
        h_prior->SetBinContent(i,val);
    }
    h_prior->Write("prior_2pcVar");

    effVar = 0.05;
    for (int i=1; i<=nBins; i++)
    {
        double val = 1. - effVar + 2.*effVar*(i-1.)/(nBins-1.);
        h_prior->SetBinContent(i,val);
    }
    h_prior->Write("prior_5pcVar");

    for (int j=1;j<=25;j++)
    {
        double err = j*0.001;
        double errpc = err*100;
        for (int i=1; i<=nBins; i++)
        {
            cov_matrix[i-1][i-1] = err*err;
        }
        cov_matrix.Write(Form("cov_%1.1fpcAbs",errpc));
    }


    f->Close();
}