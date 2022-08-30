TH2D* get_cosths_tof_hist(char* filename, double norm, double offset=0)
{
    TChain* fChain_digitized = new TChain("hitRate_pmtType1");
    fChain_digitized->Add(filename);
    
    TFile* f = new TFile(fChain_digitized->GetFile()->GetName());
    TTree* t = (TTree*)f->Get("pmt_type1");
    double R, costh, cosths, omega, phim;
    t->SetBranchAddress("R",&R);
    t->SetBranchAddress("costh",&costh);
    t->SetBranchAddress("cosths",&cosths);
    t->SetBranchAddress("omega",&omega);
    t->SetBranchAddress("phim",&phim);
    const int nPMTs = t->GetEntries();
    std::cout<<"nPMTs = "<<nPMTs<<std::endl;
    double cosths_array[nPMTs], costh_array[nPMTs], R_array[nPMTs];

    for (int i=0;i<nPMTs;i++){
        t->GetEntry(i);
        cosths_array[i]=cosths;
        costh_array[i]=costh;
        R_array[i]=R;
    }

    double nPE, timetof;
    int PMT_id;
    fChain_digitized->SetBranchAddress("nPE",&nPE);
    fChain_digitized->SetBranchAddress("timetof",&timetof);
    fChain_digitized->SetBranchAddress("PMT_id",&PMT_id);
    TH2D* h_Timetof = new TH2D("","",100,-5,5,5,0.767,1.0);
    for (unsigned long int i=0;i<fChain_digitized->GetEntries();i++){
        if (i%10000==0) std::cout<<".";
        fChain_digitized->GetEntry(i);
        if (costh_array[PMT_id]>0.5)
            h_Timetof->Fill(timetof+offset,cosths_array[PMT_id],nPE);
    }

    h_Timetof->Scale(1./norm);

    return h_Timetof;
}

void plot_scatter_slide(){
    char* nominal_file = (char*)"../diffuser4_350nm_nominal/absorption_diffuser4_*.root";
    char* no_scatter_file = (char*)"../WCSIM_output/out_diffuser_4_350nm_abs_0.096_ray_1000000.root";

    TH2D* h_nominal = get_cosths_tof_hist(nominal_file,1,950);
    //TH2D* h_nominal = get_cosths_tof_hist(no_scatter_file,1);

    TCanvas* c1 = new TCanvas();
    for (int i=1;i<=h_nominal->GetNbinsY();i++)
    {
        TH1D* hProj_dist = h_nominal->ProjectionX(Form("_px%i",i),i,i,"");
        hProj_dist->SetLineColor(i);
        //hProj_dist->SetLineWidth(3);
        hProj_dist->Scale(1/hProj_dist->GetMaximum());
        hProj_dist->Draw("hist same");
        //c1->SaveAs(Form("cosths_bin%i.pdf",i));
    }
    c1->SaveAs("cosths_bin_all.pdf");


}