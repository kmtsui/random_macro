// define a function with 4 parameters
double fitf(double *x,double *par) {
    double tpeak = par[0];
    double a = par[1];
    double sigma = par[2];

    double x0 = par[3];
    double b = par[4];
    double sigma1 = par[5];
    double sigma2 = par[6];

    double val = a*exp(-(x[0]-tpeak)*(x[0]-tpeak)/sigma/sigma);

    if (x[0]<x0)
        return val+b*exp(-(x[0]-x0)*(x[0]-x0)/sigma1/sigma1);
    else
        return val+b*exp(-(x[0]-x0)/sigma2);
}

TH1D* get_tof_hist(char* filename, double norm, double offset=0)
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
    TH1D* h_Timetof = new TH1D("","",100,-5,5);
    for (unsigned long int i=0;i<fChain_digitized->GetEntries();i++){
        if (i%10000==0) std::cout<<".";
        fChain_digitized->GetEntry(i);
        if (costh_array[PMT_id]>0.5 && cosths_array[PMT_id]>0.9 && cosths_array[PMT_id]<1.0)
            h_Timetof->Fill(timetof+offset,nPE);
    }

    h_Timetof->Scale(1./norm);

    return h_Timetof;
}

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
    TH2D* h_Timetof = new TH2D("","",100,-5,5,22,0.76,1.0);
    for (unsigned long int i=0;i<fChain_digitized->GetEntries();i++){
        if (i%10000==0) std::cout<<".";
        fChain_digitized->GetEntry(i);
        if (costh_array[PMT_id]>0.5)
            h_Timetof->Fill(timetof+offset,cosths_array[PMT_id],nPE);
    }

    h_Timetof->Scale(1./norm);

    return h_Timetof;
}


void plot_scatter_subtraction(){
    char* nominal_file = (char*)"../diffuser4_350nm_nominal/absorption_diffuser4_*.root";
    char* no_scatter_file = (char*)"../WCSIM_output/out_diffuser_4_350nm_abs_0.096_ray_1000000.root";

    TH1D* h_nominal = get_tof_hist(nominal_file,1270*10000*100,950);
    TH1D* h_no_scatter = get_tof_hist(no_scatter_file,2000*100000);

    h_nominal->SetName("h_nominal");
    h_no_scatter->SetName("h_no_scatter");

    for (int i=1;i<=h_nominal->GetNbinsX();i++)
    {
        if (i<50) continue;
        if (h_nominal->GetBinContent(i)>0 && h_nominal->GetBinContent(i)<h_no_scatter->GetBinContent(i))
        {
            double scale_factor = h_no_scatter->GetBinContent(i)/h_nominal->GetBinContent(i);
            //h_nominal->Scale(scale_factor);
        }
    }

    TCanvas* c1 = new TCanvas();
    h_nominal->SetLineColor(kBlue);
    h_nominal->SetLineWidth(3);
    h_no_scatter->SetLineColor(kRed);
    h_no_scatter->SetLineWidth(3);
    h_nominal->Draw("hist same");
    h_no_scatter->Draw("hist same");

    TCanvas* c2 = new TCanvas();
    TH1D* h_diff = (TH1D*)h_nominal->Clone();
    h_diff->Add(h_no_scatter,-1);
    h_diff->Draw("hist");

    TCanvas* c3 = new TCanvas();
    TH1D* h_nominal_fit = (TH1D*)h_nominal->Clone();
    h_nominal_fit->Scale(1270*10000*100);
    h_nominal_fit->Draw("hist");
    TF1 *func = new TF1("fit",fitf,-5,5,7);
    func->SetParNames ("tpeak","a","sigma","x0","b","sigma1","sigma2");
    func->SetParameters(h_nominal_fit->GetBinCenter(h_nominal_fit->GetMaximumBin()),h_nominal_fit->GetMaximum(),0.6,2,h_nominal_fit->GetBinContent(h_nominal_fit->GetNbinsX()),1,0.2);
    //func->FixParameter(0,h_nominal_fit->GetBinCenter(h_nominal_fit->GetMaximumBin()));
    //func->FixParameter(3,2);
    auto fitresult = h_nominal_fit->Fit("fit","S");
    fitresult->Print();
    func->SetLineColor(kRed);
    func->SetLineWidth(3);
    func->Draw("same");
}