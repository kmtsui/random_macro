// define a function with 4 parameters
double fitf(double *x,double *par) {
    double tpeak = par[0];
    double a = par[1];
    double sigma1 = par[2];
    double sigma2 = par[3];

    double argSq = (x[0]-tpeak)*(x[0]-tpeak);
    if (x[0]<tpeak)
        return a*exp(-argSq/sigma1/sigma1);
    else
        return a*exp(-argSq/sigma2/sigma2);
}

void plot(std::string filename, double timetof_min, double timetof_dt, const int timetof_nfiles, double truth_alpha){
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);

    double trigger_offset = 950;

    double chi2[timetof_nfiles];

    timetof_min += trigger_offset;

    TCanvas* c1 = new TCanvas();
    TH1D* h_alpha = new TH1D("","",timetof_nfiles,timetof_min-timetof_dt/2.,timetof_min+(timetof_nfiles-0.5)*timetof_dt);
    double R_list[5] = {1000,3000,5000,6500,8000};
    TH1D* h_attenuation[5];
    for (int i=0;i<5;i++)
        h_attenuation[i] = (TH1D*)h_alpha->Clone();
    for (int i=1;i<=h_alpha->GetNbinsX();i++){
        double timetofcut = timetof_min+(i-1)*timetof_dt - trigger_offset;
        TFile f(Form("%s_%1.1f.root",filename.c_str(), timetofcut));
        TH1D* hist_alpha_result = (TH1D*)f.Get("hist_alpha_result");
        TH1D* hist_alpha_error_final = (TH1D*)f.Get("hist_alpha_error_final");
        TH1D* chi2_total_periter = (TH1D*)f.Get("chi2_total_periter");
        chi2[i-1]=chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX()-1);
        h_alpha->SetBinContent(i,hist_alpha_result->GetBinContent(1)-70);
        h_alpha->SetBinError(i,hist_alpha_error_final->GetBinContent(1));
        //h_alpha->SetBinContent(i,best_fit[i-1]);
        //h_alpha->SetBinError(i,fit_err[i-1]);

        for (int j=0;j<5;j++)
            h_attenuation[j]->SetBinContent(i,exp(-R_list[j]/hist_alpha_result->GetBinContent(1))/exp(-R_list[j]/truth_alpha));
    }
    h_alpha->GetXaxis()->SetTitle("t_{D} cut (ns)");
    h_alpha->GetYaxis()->SetTitle("Fitted L_{#alpha} (cm)");
    h_alpha->GetYaxis()->SetTitleOffset(1.5);
    h_alpha->Draw();
    TLine *line1 = new TLine(timetof_min-timetof_dt/2.,truth_alpha,timetof_min+(timetof_nfiles-0.5)*timetof_dt,truth_alpha);
    line1->SetLineWidth(3);
    line1->SetLineStyle(9);
    line1->SetLineColor(kBlue);
    line1->Draw();
    TLatex latex;
    latex.SetTextSize(0.06);
    latex.SetTextColor(kBlue);
    latex.DrawLatex(timetof_min,truth_alpha*1.003,"True L_{#alpha}");
    c1->SaveAs(Form("%s_cut.pdf",filename.c_str()));
return;

    TCanvas* c_atten = new TCanvas();
    TLegend* legend = new TLegend(0.2,0.6,0.5,0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    for (int j=4;j>=0;j--)
    {
        legend->AddEntry(h_attenuation[j],Form("R = %i cm",(int)R_list[j]),"l");
        h_attenuation[j]->GetXaxis()->SetTitle("timetof cut (ns)");
        h_attenuation[j]->GetYaxis()->SetTitle("Ratio to true attenuation factor");
        h_attenuation[j]->SetLineWidth(3);
        h_attenuation[j]->SetLineStyle(j+1);
        h_attenuation[j]->SetLineColor(j+1);
        h_attenuation[j]->Draw("hist same");
    }
    legend->Draw("same");


return;

    TCanvas* c2 = new TCanvas();
    TH1D* h_chi2 = new TH1D("","",timetof_nfiles,timetof_min-timetof_dt/2.,timetof_min+(timetof_nfiles-0.5)*timetof_dt);
    for (int i=1;i<=h_chi2->GetNbinsX();i++){
        h_chi2->SetBinContent(i,chi2[i-1]);
    }
    h_chi2->GetXaxis()->SetTitle("timetof cut (ns)");
    h_chi2->GetYaxis()->SetTitle("Post-fit #chi^{2}");
    h_chi2->Draw();

    TChain* fChain = new TChain("hitRate_pmtType1");
    fChain->Add("../diffuser4_400nm_nominal/absorption_diffuser*.root");
    //fChain->Add("../diffuser4_350nm_x2absorption/absorption_diffuser*.root");
    //fChain->Add("/hepstore/kmtsui/hyperk/out.root");

    std::cout<<"fChain->GetFile()->GetName() = "<<fChain->GetFile()->GetName()<<std::endl;
    TFile* f = new TFile(fChain->GetFile()->GetName());
    TTree* t = (TTree*)f->Get("pmt_type1");
    double R, costh, cosths, omega, phim;
    t->SetBranchAddress("R",&R);
    t->SetBranchAddress("costh",&costh);
    t->SetBranchAddress("cosths",&cosths);
    t->SetBranchAddress("omega",&omega);
    t->SetBranchAddress("phim",&phim);
    int nPMTs = t->GetEntries();
    double cosths_array[nPMTs], costh_array[nPMTs], R_array[nPMTs];
    for (int i=0;i<nPMTs;i++){
        t->GetEntry(i);
        cosths_array[i]=cosths;
        costh_array[i]=costh;
        R_array[i]=R;
    }


    double nPE, timetof;
    int PMT_id;
    fChain->SetBranchAddress("nPE",&nPE);
    fChain->SetBranchAddress("timetof",&timetof);
    fChain->SetBranchAddress("PMT_id",&PMT_id);
    TH1D* h_Timetof = new TH1D("","",1000,-952,-942);
    TH1D* h_R_all = new TH1D("","",100,1000,9000);
    TH1D* h_R_cut = (TH1D*) h_R_all->Clone();
    TH2D* h_R_Timetof = new TH2D("","",100,-952,-942,10,1000,9000);
    for (unsigned long int i=0;i<fChain->GetEntries();i++){
        fChain->GetEntry(i);
        if (costh_array[PMT_id]>0.5 && cosths_array[PMT_id]>0.767 )
        //if (costh_array[PMT_id]>0.5 && cosths_array[PMT_id]<0.766 )
        {
            h_Timetof->Fill(timetof,nPE);
            h_R_Timetof->Fill(timetof,R_array[PMT_id],nPE);
            if (timetof<-946.5) h_R_all->Fill(R_array[PMT_id],nPE);
            if (timetof<-947.5) h_R_cut->Fill(R_array[PMT_id],nPE);
        }
    }
    TCanvas* c3 = new TCanvas();
    h_Timetof->GetXaxis()->SetTitle("timetof (ns)");
    h_Timetof->Draw("hist");
    auto fitresult = h_Timetof->Fit("gaus","S","",-951,-948);
    TF1 *fa = new TF1("fa","[0]*TMath::Exp(-0.5*(x-[1])*(x-[1])/[2]/[2])",-952,-942);
    fa->SetParameters(fitresult->Parameter(0),fitresult->Parameter(1),fitresult->Parameter(2));
    fa->SetLineColor(kRed);
    fa->Draw("same");
    double integral_timewindow = 0;
    double integral_all = 0;
    for (int i = 1; i<= h_Timetof->GetNbinsX();i++) {
        if (h_Timetof->GetBinCenter(i)<-947) {integral_timewindow+=h_Timetof->GetBinContent(i);integral_all=integral_timewindow;}
        else {
            integral_all+=h_Timetof->GetBinContent(i);
            if (integral_all>=2*integral_timewindow) {
                std::cout<<"Cut at timetof = "<<h_Timetof->GetBinCenter(i)<<std::endl;
                break;
            }
        }
    }
    std::cout<<"Direct photon time window integral = "<<h_Timetof->Integral()<<std::endl;

    TCanvas* c4 = new TCanvas();
    h_R_cut->Divide(h_R_all);
    h_R_cut->Draw();

    const int np = 2;
    double x[np] = {0.97,0.99};
    double y[np] = {0.99,1.02};
    TGraph gr(np,x,y);

    TCanvas* c5 = new TCanvas();
    h_R_Timetof->Draw("colz");
    for(int i=1;i<=h_R_Timetof->GetNbinsY();i++){
        TH1D* hProj_dist = h_R_Timetof->ProjectionX(Form("_px%i",i),i,i,"e");
        double peak_count = hProj_dist->Integral(1,50);
        double tail_count = hProj_dist->Integral(51,65);
        std::cout<<"Bin "<<i<<", ratio = "<<1./(1.-tail_count/peak_count)<<std::endl;
        //std::cout<<"Bin "<<i<<", ratio = "<<peak_count/(peak_count+tail_count)<<", correction = "<<1./gr.Eval(peak_count/(peak_count+tail_count))<<std::endl;
    }
}

void plot_timetof_cut(){
    //plot("time_cut/diffuser4_400nm_nominal_mPMT_timetof", -948, 0.1, 21, 10648);
    plot("time_cut/diffuser4_400nm_nominal_BnLPMT_timetof", -948, 0.5, 21, 10648);
}