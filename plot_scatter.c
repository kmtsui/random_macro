void plot_scatter(){
    TChain* fChain_digitized = new TChain("hitRate_pmtType1");
    fChain_digitized->Add("../diffuser4_350nm_nominal/absorption_diffuser*.root");

    TChain* fChain_raw = new TChain("hitRate_pmtType1");
    fChain_raw->Add("../diffuser4_350nm_x2absorption/raw_*.root");

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
    double peak_count[nPMTs], tail_count[nPMTs], ratio_count[nPMTs];

    for (int i=0;i<nPMTs;i++){
        t->GetEntry(i);
        cosths_array[i]=cosths;
        costh_array[i]=costh;
        R_array[i]=R;
        peak_count[i]=0;
        tail_count[i]=0;
        //ratio_count[i]=0;
    }

    double nPE, timetof;
    int PMT_id;
    fChain_digitized->SetBranchAddress("nPE",&nPE);
    fChain_digitized->SetBranchAddress("timetof",&timetof);
    fChain_digitized->SetBranchAddress("PMT_id",&PMT_id);
    TH1D* h_Timetof = new TH1D("","",100,-955,-945);
    TH1D* h_R_digitized_3ns = new TH1D("","",10,1000,9000);
    TH1D* h_R_digitized_5ns = (TH1D*)h_R_digitized_3ns->Clone();
    TH1D* h_R_raw = (TH1D*)h_R_digitized_3ns->Clone();
    for (unsigned long int i=0;i<fChain_digitized->GetEntries();i++){
        fChain_digitized->GetEntry(i);
        if (timetof<-947) peak_count[PMT_id]+=nPE;
        if (timetof>-947&&timetof<-945) tail_count[PMT_id]+=nPE;
        if (costh_array[PMT_id]>0.5 && cosths_array[PMT_id]>0.767 )
        //if (cosths_array[PMT_id]<0.766 )
        {
            h_Timetof->Fill(timetof,nPE);
            if (timetof<-947) h_R_digitized_3ns->Fill(R_array[PMT_id],nPE);
            if (timetof<-945) h_R_digitized_5ns->Fill(R_array[PMT_id],nPE);
        }
    }

    TCanvas* c1 = new TCanvas();
    TGraph* gr = new TGraph();
    TH1D* h_R_peak = new TH1D("","",20,1000,9000);
    TH1D* h_costh_peak = new TH1D("","",20,0.5,1);
    TH1D* h_cosths_peak = new TH1D("","",10,0.767,1);
    TH1D* h_R_tail = (TH1D*)h_R_peak->Clone();
    TH1D* h_costh_tail = (TH1D*)h_costh_peak->Clone();
    TH1D* h_cosths_tail = (TH1D*)h_cosths_peak->Clone();
    int nPoints = 0;
    for (int i=0;i<nPMTs;i++){
        if (cosths_array[i]<0.767) continue;
        if (costh_array[i]<0.5) continue;
        if (peak_count[i]<=0) continue;
        double ratio = tail_count[i]/peak_count[i];
        gr->SetPoint(nPoints,cosths_array[i],ratio);
        nPoints++;

        h_R_peak->Fill(R_array[i],peak_count[i]);
        h_R_tail->Fill(R_array[i],tail_count[i]);
        h_costh_peak->Fill(costh_array[i],peak_count[i]);
        h_costh_tail->Fill(costh_array[i],tail_count[i]);
        h_cosths_peak->Fill(cosths_array[i],peak_count[i]);
        h_cosths_tail->Fill(cosths_array[i],tail_count[i]);
    }
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->Draw("AP");
    c1->SaveAs("test.pdf");
    
    gStyle->SetOptStat(0);

    h_R_tail->Divide(h_R_peak);
    h_R_tail->GetXaxis()->SetTitle("R (cm)");
    h_R_tail->GetYaxis()->SetTitle("Ratio");
    h_R_tail->Draw();
    c1->SaveAs("h_R_ratio.pdf");

    h_costh_tail->Divide(h_costh_peak);
    h_costh_tail->GetXaxis()->SetTitle("cos#theta");
    h_costh_tail->GetYaxis()->SetTitle("Ratio");
    h_costh_tail->Draw();
    c1->SaveAs("h_costh_ratio.pdf");

    h_cosths_tail->Divide(h_cosths_peak);
    h_cosths_tail->GetXaxis()->SetTitle("cos#theta_{s}");
    h_cosths_tail->GetYaxis()->SetTitle("Ratio");
    h_cosths_tail->Draw();
    c1->SaveAs("h_cosths_ratio.pdf");

    h_Timetof->GetXaxis()->SetTitle("timetof (ns)");
    h_Timetof->GetYaxis()->SetTitle("Count");
    h_Timetof->Draw("hist");
    c1->SaveAs("h_timetof.pdf");

    // fChain_raw->SetBranchAddress("nPE",&nPE);
    // fChain_raw->SetBranchAddress("timetof",&timetof);
    // fChain_raw->SetBranchAddress("PMT_id",&PMT_id);
    // for (unsigned long int i=0;i<fChain_raw->GetEntries();i++){
    //     fChain_raw->GetEntry(i);
    //     if (costh_array[PMT_id]>0.5 && cosths_array[PMT_id]>0.767 )
    //     {
    //         if (timetof<0.3) h_R_raw->Fill(R_array[PMT_id],nPE);
    //     }
    // }
    // TGraph* gr = new TGraph();
    // int nPoints = 0;
    // for (int i=1;i<=h_R_digitized_3ns->GetNbinsX();i++){
    //     if (h_R_digitized_3ns->GetBinContent(i)<=0) continue;
    //     double x = h_R_digitized_3ns->GetBinContent(i)/h_R_digitized_5ns->GetBinContent(i);
    //     double y = h_R_raw->GetBinContent(i)/h_R_digitized_3ns->GetBinContent(i);
    //     gr->SetPoint(nPoints,x,y);
    //     nPoints++;
    // }
    // gr->SetMarkerColor(4);
    // gr->SetMarkerStyle(21);
    // gr->Draw("AP");
}