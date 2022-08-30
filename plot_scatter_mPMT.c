void plot_scatter_mPMT(){

    gStyle->SetOptStat(0);
    TFile* f400 = new TFile("diffuser4_400nm_scattering_map_mPMT.root");
    TFile* f350 = new TFile("diffuser4_350nm_scattering_map_mPMT.root");
    TH1D* h400 = (TH1D*)f400->Get("scattering_map_2_200");
    TH1D* h350 = (TH1D*)f350->Get("scattering_map_2_200");
    TCanvas* c = new TCanvas();
    h400->GetYaxis()->SetRangeUser(0,1.6);
    h400->GetXaxis()->SetTitle("PMT_id");
    h400->GetYaxis()->SetTitle("r");
    h400->Draw("hist");
    h350->SetLineColor(kRed);
    h350->Draw("hist same");
    c->SaveAs("test.pdf");

return;
    TChain* fChain_digitized = new TChain("hitRate_pmtType1");
    //fChain_digitized->Add("../WCSIM_output/out_diffuser_4_350nm_nominal_trackscattering.root");
    //fChain_digitized->Add("../WCSIM_output/out_diffuser_4_350nm_abs_1.3_ray_0.5_trackscattering.root");
    //fChain_digitized->Add("../production/out_diffuser4_350nm_nominal_*.root");
    fChain_digitized->Add("../production/out_diffuser4_350nm_nominal_*.root");

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
    std::vector<double> cosths_array, costh_array, R_array;
    //double peak_count[nPMTs], tail_count[nPMTs], ratio_count[nPMTs];

    for (int i=0;i<nPMTs;i++){
         t->GetEntry(i);
         cosths_array.push_back(cosths);
         costh_array.push_back(costh);
         R_array.push_back(R);
    //     R_array[i]=R;
    //     peak_count[i]=0;
    //     tail_count[i]=0;
    //     //ratio_count[i]=0;
    }

    double nPE, timetof, nPE_digi, timetof_digi;
    int PMT_id, nRaySct, nReflec;
    fChain_digitized->SetBranchAddress("nPE",&nPE);
    fChain_digitized->SetBranchAddress("timetof",&timetof);
    fChain_digitized->SetBranchAddress("PMT_id",&PMT_id);
    fChain_digitized->SetBranchAddress("nRaySct",&nRaySct);
    fChain_digitized->SetBranchAddress("nReflec",&nReflec);
    fChain_digitized->SetBranchAddress("nPE_digi",&nPE_digi);
    fChain_digitized->SetBranchAddress("timetof_digi",&timetof_digi);
    TH1D* h_timetof_all_10 = new TH1D("","",nPMTs,0,nPMTs);
    TH1D* h_timetof_10 = new TH1D("","",nPMTs,0,nPMTs);
    TH1D* h_timetof_10_200 = new TH1D("","",nPMTs,0,nPMTs);
    TH1D* h_timetof_all = new TH1D("","",100,-10,20);
    TH1D* h_timetof_reflec = (TH1D*)h_timetof_all->Clone();
    TH1D* h_timetof_raysct = (TH1D*)h_timetof_all->Clone();
    TH1D* h_indirect = new TH1D("","",1000,-1,1);
    for (unsigned long int i=0;i<fChain_digitized->GetEntries();i++){
        fChain_digitized->GetEntry(i);
        //if (cosths_array[PMT_id]<0.767) continue;
        //if (costh_array[PMT_id]<0.5) continue;
        if (timetof_digi<10) h_timetof_all_10->Fill(PMT_id+0.5,nPE_digi);
        h_timetof_all->Fill(timetof_digi,nPE_digi);
        // if (nRaySct+nReflec==0) 
        // //if (timetof<0) 
        // {
        //     continue;
        // }
        if (nReflec>0) h_timetof_reflec->Fill(timetof_digi,nPE_digi);
        if (nRaySct>0) h_timetof_raysct->Fill(timetof_digi,nPE_digi);
        if (nReflec+nRaySct>0) h_indirect->Fill(timetof,nPE);
        // if (nRaySct>0 || (nReflec>0 && timetof>5))
        // {
        //     if (timetof_digi<2) h_timetof_10->Fill(PMT_id+0.5,nPE_digi);
        //     else if (timetof_digi<200) h_timetof_10_200->Fill(PMT_id+0.5,nPE_digi);
        //     // if (timetof<0) h_timetof_10->Fill(PMT_id+0.5);
        //     // else if (timetof<200) h_timetof_10_200->Fill(PMT_id+0.5);
        // }
        if ((timetof>0.6) && timetof_digi<2) h_timetof_10->Fill(PMT_id+0.5,nPE_digi);
        else if (timetof_digi>2&&timetof_digi<200) h_timetof_10_200->Fill(PMT_id+0.5);
    }

    TFile* fout = new TFile("diffuser4_350nm_scattering_map_mPMT.root","RECREATE");
    double sratio;
    TTree* st = new TTree("st","st");
    st->Branch("cosths",&cosths);
    st->Branch("costh",&costh);
    st->Branch("R",&R);
    st->Branch("sratio",&sratio);

    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas();
    TGraph* gr = new TGraph();
    int nPoints = 0;
    TH3D* h_pmt_map = new TH3D("","",10,0.767,1,10,0.5,1,10,1000,9000);
    TH3D* h_scatter_map = (TH3D*)h_pmt_map->Clone();
    for (int i=0;i<nPMTs;i++){
        //if (i>15500)
        if (h_timetof_10_200->GetBinContent(i+1)>0)
        {
            double x = h_timetof_10_200->GetBinContent(i+1);
            double y = h_timetof_10->GetBinContent(i+1);
            //gr->SetPoint(nPoints,h_timetof_10_200->GetBinContent(i+1)/h_timetof_all_10->GetBinContent(i+1),h_timetof_10->GetBinContent(i+1)/h_timetof_all_10->GetBinContent(i+1));
            gr->SetPoint(nPoints,x,y);
            nPoints++;

            double val = y/x;
            h_pmt_map->Fill(cosths_array[i],costh_array[i],R_array[i]);
            h_scatter_map->Fill(cosths_array[i],costh_array[i],R_array[i],val);

            cosths = cosths_array[i];
            costh = costh_array[i];
            R = R_array[i];
            sratio = val;
            st->Fill();
        }
    }
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->GetXaxis()->SetTitle("PE_{indirect,3<t_{D}<200}");
    gr->GetYaxis()->SetTitle("PE_{indirect,t_{D}<3}");
    gr->Draw("AP");
    //TFitResultPtr r = gr->Fit("pol2","S");
    //const ROOT::Fit::FitResult::IModelFunction * func = r->FittedFunction();
    c1->SaveAs("test1.pdf");

    // double x[] = {0};
    // for (int i=1;i<=h_timetof_10->GetNbinsX();i++)
    // {
    //     x[0] = h_timetof_10->GetBinContent(i);
    //     h_timetof_10->SetBinContent(i,0);
    //     if (x[0]<=0) continue;
    //     double y = (*func)(x);
    //     if (y<0) continue;
    //     h_timetof_10->SetBinContent(i,y/x[0]);
    // }
    for (int i=1;i<=h_timetof_10->GetNbinsX();i++)
    {
        double x = h_timetof_10_200->GetBinContent(i);
        double y = h_timetof_10->GetBinContent(i);
        if (x>0&&y>0)
        {
            double val = y/x;
            double err = sqrt(1./x+1./y);
            // double val = y;
            // double err = sqrt(1./y);
            h_timetof_10->SetBinContent(i,val);
            h_timetof_10->SetBinError(i,err);
        }
    }
    //h_timetof_10->Divide(h_timetof_10_200);
    h_timetof_10->Draw("hist");
    c1->SaveAs("test2.pdf");

    h_timetof_all->GetXaxis()->SetTitle("t_{D} (ns)");
    h_timetof_all->Draw("hist");
    h_timetof_reflec->SetLineColor(kRed);
    h_timetof_reflec->Draw("hist same");
    h_timetof_raysct->SetLineColor(kGreen);
    h_timetof_raysct->Draw("hist same");
    c1->SaveAs("test3.pdf");

    h_indirect->GetXaxis()->SetTitle("t_{c} (ns)");
    h_indirect->Draw("hist");
    c1->SaveAs("test4.pdf");

    h_timetof_10->Write("scattering_map_2_200");
    h_scatter_map->Divide(h_pmt_map);
    h_scatter_map->Write("h_scatter_map");
    st->Write();
    fout->Close();
}