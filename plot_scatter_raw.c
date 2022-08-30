void plot_scatter_raw(){
    TChain* fChain_digitized = new TChain("hitRate_pmtType1");
    fChain_digitized->Add("/hepstore/kmtsui/hyperk/out_diffuser_4_350nm_nominal_trackscattering_digitized.root");

    TChain* fChain_raw = new TChain("hitRate_pmtType1");
    fChain_raw->Add("/hepstore/kmtsui/hyperk/out_diffuser_4_350nm_nominal_trackscattering.root");

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
    //double cosths_array[nPMTs], costh_array[nPMTs];//, R_array[nPMTs];
    //double peak_count[nPMTs], direct_photons[nPMTs];
    //double tail_count[nPMTs];//, scatter_photons[nPMTs];
    std::vector<double> cosths_array, costh_array, peak_count, direct_photons, tail_count, scatter_photons;

    for (int i=0;i<nPMTs;i++){
        t->GetEntry(i);
        //cosths_array[i]=cosths;
        //costh_array[i]=costh;
        //R_array[i]=R;
        //peak_count[i]=0;
        //tail_count[i]=0;
        //direct_photons[i]=0;
        //scatter_photons[i]=0;
        cosths_array.push_back(cosths);
        costh_array.push_back(costh);
        peak_count.push_back(0);
        direct_photons.push_back(0);
        tail_count.push_back(0);
        scatter_photons.push_back(0);
    }

    double nPE, timetof;
    int PMT_id;
    fChain_digitized->SetBranchAddress("nPE",&nPE);
    fChain_digitized->SetBranchAddress("timetof",&timetof);
    fChain_digitized->SetBranchAddress("PMT_id",&PMT_id);
    for (unsigned long int i=0;i<fChain_digitized->GetEntries();i++){
        fChain_digitized->GetEntry(i);
        if (timetof<3) peak_count[PMT_id]+=nPE;
        if (timetof>3&&timetof<5) tail_count[PMT_id]+=nPE;
    }

    TH1D* h_Time_direct  = new TH1D("","",1000,-1,1);
    TH1D* h_Time_scatter = new TH1D("","",1000,-1,1);

    double nPE_scatter;
    fChain_raw->SetBranchAddress("nPE",&nPE);
    fChain_raw->SetBranchAddress("nPE_scatter",&nPE_scatter);
    fChain_raw->SetBranchAddress("timetof",&timetof);
    fChain_raw->SetBranchAddress("PMT_id",&PMT_id);
    for (unsigned long int i=0;i<fChain_raw->GetEntries();i++){
        fChain_raw->GetEntry(i);
        if (timetof<0.3) direct_photons[PMT_id]+=nPE-nPE_scatter;
        if (timetof<0.3) scatter_photons[PMT_id]+=nPE_scatter;
        if (costh_array[PMT_id]>0.5 && cosths_array[PMT_id]>0.767)
        {
            h_Time_direct->Fill(timetof,nPE-nPE_scatter);
            h_Time_scatter->Fill(timetof,nPE_scatter);
        }
    }

    // TH1D* h_cosths_direct = new TH1D("","",20,0.767,1);
    // TH1D* h_cosths_peak   = new TH1D("","",20,0.767,1);
    // TH1D* h_cosths_tail   = new TH1D("","",20,0.767,1);
    // TH1D* h_ratio_direct = new TH1D("","",10,0.,0.05);
    // TH1D* h_ratio_peak   = new TH1D("","",10,0.,0.05);
    // for (int i=0;i<nPMTs;i++){
    //     if (costh_array[i]>0.5 && cosths_array[i]>0.767) {
    //         h_cosths_direct->Fill(cosths_array[i],direct_photons[i]);
    //         h_cosths_peak->Fill(cosths_array[i],peak_count[i]);
    //         h_cosths_tail->Fill(cosths_array[i],tail_count[i]);

    //         if (peak_count[i]>0){
    //             double ratio = tail_count[i]/direct_photons[i];
    //             h_ratio_direct->Fill(ratio,direct_photons[i]);
    //             h_ratio_peak->Fill(ratio,peak_count[i]);
    //         }
    //     }
    // }
    // TCanvas* c1 = new TCanvas();
    // h_cosths_direct->Divide(h_cosths_peak);
    // h_cosths_direct->GetXaxis()->SetTitle("cos#theta_{s}");
    // h_cosths_direct->Draw();

    // TCanvas* c2 = new TCanvas();
    // h_cosths_tail->Divide(h_cosths_peak);
    // h_cosths_tail->GetXaxis()->SetTitle("cos#theta_{s}");
    // h_cosths_tail->GetYaxis()->SetTitle("Ratio");
    // h_cosths_tail->Draw();

    // gStyle->SetOptStat(0);
    // TCanvas* c3 = new TCanvas();
    // h_Time_direct->GetXaxis()->SetTitle("timetof (ns)");
    // h_Time_direct->GetYaxis()->SetTitle("Count");
    // h_Time_scatter->SetLineColor(kRed);
    // h_Time_direct->Draw("hist");
    // h_Time_scatter->Draw("hist same");

    TH1D* h_direct_ratio = new TH1D("h_direct_ratio","",10,0.767,1.);
    TH1D* h_tail_ratio = new TH1D("h_tail_ratio","",10,0.767,1.);
    TH1D* h_tail_count = new TH1D("h_tail_count","",10,0.767,1.);
    for (int i=0;i<nPMTs;i++){
        if (costh_array[i]>0.5 && costh_array[i]<1. && cosths_array[i]>0.767 && cosths_array[i]<1. && peak_count[i]>0 ) {
            double x = tail_count[i]/peak_count[i]; 
            double y = scatter_photons[i]/(direct_photons[i]);
            h_direct_ratio->Fill(cosths_array[i],y);
            h_tail_ratio->Fill(cosths_array[i],x);
            h_tail_count->Fill(cosths_array[i],1);
        }
    }

    TCanvas* c1 = new TCanvas();
    h_direct_ratio->Divide(h_tail_count);
    h_direct_ratio->Draw("hist");
    
    TCanvas* c2 = new TCanvas();
    h_tail_ratio->Divide(h_tail_count);
    h_tail_ratio->Draw("hist");

    TCanvas* c3 = new TCanvas();
    h_tail_count->Draw("hist");
    
return;
    TGraph* gr = new TGraph();
    int nPoints = 0;
    for (int i=0;i<nPMTs;i++){
        if (costh_array[i]>0.5 && cosths_array[i]>0.767 && peak_count[i]>0 && tail_count[i]>0) {
            double x = tail_count[i]/peak_count[i];
            double y = direct_photons[i]/peak_count[i];
            //double x = tail_count[i];
            //double y = scatter_photons[i];
            gr->SetPoint(nPoints,x,y);
            nPoints++;
        }
    }
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->Draw("AP");
}