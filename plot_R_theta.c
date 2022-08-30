void plot_R_theta_mPMT(){
    TChain* fChain_pmtType1_digitized = new TChain("hitRate_pmtType1");
    fChain_pmtType1_digitized->Add("../diffuser4_400nm_nominal/absorption_diffuser*.root");

    TFile* f = new TFile(fChain_pmtType1_digitized->GetFile()->GetName());
    TTree* t = (TTree*)f->Get("pmt_type1");
    double R, costh, cosths, omega, phim;
    t->SetBranchAddress("R",&R);
    t->SetBranchAddress("costh",&costh);
    t->SetBranchAddress("cosths",&cosths);
    t->SetBranchAddress("omega",&omega);
    t->SetBranchAddress("phim",&phim);
    const int nPMTs = t->GetEntries();
    std::cout<<"nPMTs = "<<nPMTs<<std::endl;
    std::vector<double> cosths_array, costh_array, R_array;
    std::vector<double> peak_count, tail_count, ratio_count;

    bool mPMT = true;

    int nPMT_total = t->GetEntries();
    int nPMTpermPMT = 19;
    if (mPMT) nPMT_total = nPMT_total/nPMTpermPMT;

    int nPMT_active = 2000;
    double PMT_frac = (nPMT_active+0.)/(nPMT_total);
    int PMT_count = 0;
    std::vector<bool> pmt_mask;
    for (int i=0;i<nPMT_total;i++)
    {
        if ((PMT_count+0.)/(i+1.)<PMT_frac && PMT_count<nPMT_active) 
        {
            if (!mPMT) pmt_mask.push_back(0);
            else 
            {    
                for (int j=i*nPMTpermPMT;j<(i+1)*nPMTpermPMT;j++) {
                    pmt_mask.push_back(0);
                }
            }
            PMT_count++;
        } 
        else 
        {
            if (!mPMT) pmt_mask.push_back(1);
            else
            {
                for (int j=i*nPMTpermPMT;j<(i+1)*nPMTpermPMT;j++) {
                    pmt_mask.push_back(1);
                }
            }
        }
    }


    TH2D* h_R_theta = new TH2D("","",20,0,1,20,1000,9000);

    for (int i=0;i<nPMTs;i++){
        t->GetEntry(i);
        if (!pmt_mask[i]) if (cosths>0.767) h_R_theta->Fill(costh,R);
    }

    gStyle->SetOptStat(0);

    TCanvas* c1 = new TCanvas();
    h_R_theta->GetXaxis()->SetTitle("cos#theta");
    h_R_theta->GetYaxis()->SetTitle("R (cm)");
    h_R_theta->Draw("colz");
    c1->SaveAs("R_theta_diffuser4_mPMT.pdf");

}

void plot_R_theta_BnL(){
    TChain* fChain_pmtType1_digitized = new TChain("hitRate_pmtType1");
    fChain_pmtType1_digitized->Add("../diffuser4_400nm_nominal/absorption_diffuser*.root");

    TFile* f = new TFile(fChain_pmtType1_digitized->GetFile()->GetName());
    TTree* t = (TTree*)f->Get("pmt_type0");
    double R, costh, cosths, omega, phim;
    t->SetBranchAddress("R",&R);
    t->SetBranchAddress("costh",&costh);
    t->SetBranchAddress("cosths",&cosths);
    t->SetBranchAddress("omega",&omega);
    t->SetBranchAddress("phim",&phim);
    const int nPMTs = t->GetEntries();
    std::cout<<"nPMTs = "<<nPMTs<<std::endl;
    std::vector<double> cosths_array, costh_array, R_array;
    std::vector<double> peak_count, tail_count, ratio_count;

    TH2D* h_R_theta = new TH2D("","",20,0,1,20,1000,9000);

    for (int i=0;i<nPMTs;i++){
        t->GetEntry(i);
        if (cosths>0.767) h_R_theta->Fill(costh,R);
    }

    gStyle->SetOptStat(0);

    TCanvas* c1 = new TCanvas();
    h_R_theta->GetXaxis()->SetTitle("cos#theta");
    h_R_theta->GetYaxis()->SetTitle("R (cm)");
    h_R_theta->Draw("colz");
    c1->SaveAs("R_theta_diffuser4_BnLPMT.pdf");

}

void plot_R_theta(){
    plot_R_theta_BnL();
}