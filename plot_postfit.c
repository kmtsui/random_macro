void plot_postfit(){
    TFile* f = new TFile("fitoutput_mPMT_test_superfineangular.root");
    TTree* t = (TTree*)f->Get("PMTTree");

    // Declaration of leaf types
    int sampleId;
    double nPE;
    double R;
    double costh;
    double cosths;
    double costhm;
    double phim;
    double omega;
    int PMT_id;
    int mPMT_id;
    std::vector<double>* weight = 0;
    t->SetBranchAddress("sampleId",&sampleId);
    t->SetBranchAddress("nPE",&nPE);
    t->SetBranchAddress("R",&R);
    t->SetBranchAddress("costh",&costh);
    t->SetBranchAddress("cosths",&cosths);
    t->SetBranchAddress("costhm",&costhm);
    t->SetBranchAddress("phim",&phim);
    t->SetBranchAddress("omega",&omega);
    t->SetBranchAddress("PMT_id",&PMT_id);
    t->SetBranchAddress("mPMT_id",&mPMT_id);
    t->SetBranchAddress("weight",&weight);

    TF1* attenuation = new TF1("attenuation","exp(-x/10650.)",1000,9000);
    //attenuation->Draw();

    const int nPMTs = t->GetEntries();
    double R_array[nPMTs], attenuation_array[nPMTs];
    double R_array_err[nPMTs], attenuation_array_err[nPMTs];

    for (int i=0;i<nPMTs;i++){
        t->GetEntry(i);
        double frac_err = 1./sqrt(nPE);
        double weight_fit = (*weight)[0]*omega*(*weight)[2];
        double val = nPE/weight_fit;
        R_array[i] = R;
        attenuation_array[i] = val;
        R_array_err[i] = 0;
        attenuation_array_err[i] = val*frac_err;
        //if (i%1000==0) std::cout<<R_array[i]<<" "<<attenuation_array[i]<<std::endl;
    }

    TGraphErrors* gr = new TGraphErrors(nPMTs,R_array,attenuation_array,R_array_err,attenuation_array_err);
    gr->SetMarkerColor(4);
    gr->SetMarkerStyle(21);
    gr->Draw("AP");
    attenuation->Draw("same");
}