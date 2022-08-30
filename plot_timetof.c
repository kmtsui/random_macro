void plot_timetof(){
    // TChain* fChain_pmtType1_digitized = new TChain("hitRate_pmtType1");
    // fChain_pmtType1_digitized->Add("../production/out_diffuser4_400nm_nominal_*.root");

    // TChain* fChain_pmtType0_digitized = new TChain("hitRate_pmtType0");
    // fChain_pmtType0_digitized->Add("../production/out_diffuser4_400nm_nominal_*.root");

    // TChain* fChain_pmtType1_digitized = new TChain("hitRate_pmtType1");
    // fChain_pmtType1_digitized->Add("../collimator32_400nm_nominal/out_collimator32_nominal_*.root");

    // TChain* fChain_pmtType0_digitized = new TChain("hitRate_pmtType0");
    // fChain_pmtType0_digitized->Add("../collimator32_400nm_nominal/out_collimator32_nominal_*.root");

    TChain* fChain_pmtType1_digitized = new TChain("hitRate_pmtType1");
    fChain_pmtType1_digitized->Add("../production/diffuser4_400nm_nominal_raw.root");

    TChain* fChain_pmtType0_digitized = new TChain("hitRate_pmtType1");
    fChain_pmtType0_digitized->Add("../production/diffuser4_400nm_nominal_raw.root");

    // TChain* fChain_pmtType1_raw = new TChain("hitRate_pmtType1");
    // fChain_pmtType1_raw->Add("../production/out_diffuser4_400nm_nominal_*.root");

    // TChain* fChain_pmtType0_raw = new TChain("hitRate_pmtType0");
    // fChain_pmtType0_raw->Add("../production/out_diffuser4_400nm_nominal_*.root");

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

    for (int i=0;i<nPMTs;i++){
        t->GetEntry(i);
        cosths_array.push_back(cosths);
        costh_array.push_back(costh);
        R_array.push_back(R);
    }

    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas();

    double nPE, timetof, nPE_digi, timetof_digi;
    int PMT_id;

    // fChain_pmtType1_raw->SetBranchAddress("nPE",&nPE);
    // fChain_pmtType1_raw->SetBranchAddress("timetof",&timetof);
    // fChain_pmtType1_raw->SetBranchAddress("PMT_id",&PMT_id);
    // TH1D* h_Timetof_pmtType1_raw = new TH1D("","",100,-0.5,0.5);
    // for (unsigned long int i=0;i<fChain_pmtType1_raw->GetEntries();i++){
    //     fChain_pmtType1_raw->GetEntry(i);
    //     if (R_array[PMT_id]>5000)
    //     h_Timetof_pmtType1_raw->Fill(timetof);
    // }
    // h_Timetof_pmtType1_raw->GetXaxis()->SetTitle("t_{c} (ns)");
    // h_Timetof_pmtType1_raw->Draw();
    // c1->SaveAs("plots/timetof_raw_mPMT.pdf");

    // fChain_pmtType0_raw->SetBranchAddress("nPE",&nPE);
    // fChain_pmtType0_raw->SetBranchAddress("timetof",&timetof);
    // fChain_pmtType0_raw->SetBranchAddress("PMT_id",&PMT_id);
    // TH1D* h_Timetof_pmtType0_raw = new TH1D("","",100,-1.5,1.5);
    // for (unsigned long int i=0;i<fChain_pmtType0_raw->GetEntries();i++){
    //     fChain_pmtType0_raw->GetEntry(i);
    //     h_Timetof_pmtType0_raw->Fill(timetof);
    // }
    // h_Timetof_pmtType0_raw->GetXaxis()->SetTitle("t_{c} (ns)");
    // h_Timetof_pmtType0_raw->Draw();
    // c1->SaveAs("plots/timetof_raw_BnLPMT.pdf");

    // fChain_pmtType1_digitized->SetBranchAddress("nPE",&nPE);
    // fChain_pmtType1_digitized->SetBranchAddress("timetof",&timetof);
    // fChain_pmtType1_digitized->SetBranchAddress("nPE_digi",&nPE_digi);
    // fChain_pmtType1_digitized->SetBranchAddress("timetof_digi",&timetof_digi);
    // fChain_pmtType1_digitized->SetBranchAddress("PMT_id",&PMT_id);
    // TH1D* h_Timetof_pmtType1_digitized_direct = new TH1D("","",100,-5,5);
    // TH1D* h_Timetof_pmtType1_digitized_indirect = new TH1D("","",100,-5,5);
    // for (unsigned long int i=0;i<fChain_pmtType1_digitized->GetEntries();i++){
    //     fChain_pmtType1_digitized->GetEntry(i);
    //     if (timetof<0.6) h_Timetof_pmtType1_digitized_direct->Fill(timetof_digi+1,nPE_digi);
    //     else h_Timetof_pmtType1_digitized_indirect->Fill(timetof_digi+1,nPE_digi);
    // }
    // h_Timetof_pmtType1_digitized_direct->Add(h_Timetof_pmtType1_digitized_indirect);
    // h_Timetof_pmtType1_digitized_direct->GetXaxis()->SetTitle("t_{D} (ns)");
    // h_Timetof_pmtType1_digitized_direct->Draw("hist");
    // h_Timetof_pmtType1_digitized_indirect->SetLineColor(kRed);
    // h_Timetof_pmtType1_digitized_indirect->Draw("hist same");
    // TLegend* legend = new TLegend(0.2,0.5,0.4,0.65);
    // legend->AddEntry(h_Timetof_pmtType1_digitized_direct,"All","l");
    // legend->AddEntry(h_Timetof_pmtType1_digitized_indirect,"Indirect","l");
    // legend->SetBorderSize(0);
    // legend->Draw();
    // c1->SaveAs("plots/timetof_digitized_direct_indirect_mPMT.pdf");

    fChain_pmtType0_digitized->SetBranchAddress("nPE",&nPE);
    fChain_pmtType0_digitized->SetBranchAddress("timetof",&timetof);
    fChain_pmtType0_digitized->SetBranchAddress("nPE_digi",&nPE_digi);
    fChain_pmtType0_digitized->SetBranchAddress("timetof_digi",&timetof_digi);
    fChain_pmtType0_digitized->SetBranchAddress("PMT_id",&PMT_id);
    TH1D* h_Timetof_pmtType0_digitized_direct = new TH1D("","",100,-5,5);
    TH1D* h_Timetof_pmtType0_digitized_indirect = new TH1D("","",100,-5,5);
    for (unsigned long int i=0;i<fChain_pmtType0_digitized->GetEntries();i++){
        fChain_pmtType0_digitized->GetEntry(i);
        if (timetof<0.3)  h_Timetof_pmtType0_digitized_direct->Fill(timetof_digi,nPE_digi);
        else h_Timetof_pmtType0_digitized_indirect->Fill(timetof_digi,nPE_digi);
    }
    h_Timetof_pmtType0_digitized_direct->Add(h_Timetof_pmtType0_digitized_indirect);
    h_Timetof_pmtType0_digitized_direct->GetXaxis()->SetTitle("t_{D} (ns)");
    h_Timetof_pmtType0_digitized_direct->Draw("hist");
    h_Timetof_pmtType0_digitized_indirect->SetLineColor(kRed);
    h_Timetof_pmtType0_digitized_indirect->Draw("hist same");
    //TLegend* legend = new TLegend(0.5,0.5,0.7,0.65);
    TLegend* legend = new TLegend(0.65,0.5,0.85,0.65);
    legend->AddEntry(h_Timetof_pmtType0_digitized_direct,"All","l");
    legend->AddEntry(h_Timetof_pmtType0_digitized_indirect,"t_{c}>0.3ns","l");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs("plots/timetof_digitized_direct_indirect_mPMT.pdf");
}