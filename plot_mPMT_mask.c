void plot_mPMT_mask(){
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(0);

    const int nFiles = 12;
    double nPMT[nFiles] = {500, 1000, 1500, 2000, 3000, 4000, 
                        5000,6000, 7000, 8000, 9000, 9616};
    double best_fit[nFiles], fit_err[nFiles], chi2[nFiles], frac_err[nFiles];

    double binEdges[nFiles+1];
    binEdges[0] = nPMT[0] - (nPMT[1]-nPMT[0])/2.;
    for (int i=1;i<nFiles;i++){
        binEdges[i] = (nPMT[i-1] + nPMT[i])/2;
    }
    binEdges[nFiles] = (nPMT[nFiles-1] - nPMT[nFiles-2])/2.+nPMT[nFiles-1];
    TH1D* hist_alpha = new TH1D("","",nFiles,binEdges);
    
    double truth_alpha = 10650;
    for (int i=0;i<nFiles;i++){
        TFile f(Form("fitoutput_%imPMT.root",(int)nPMT[i]));
        TH1D* hist_alpha_result = (TH1D*)f.Get("hist_alpha_result");
        TH1D* hist_alpha_error_final = (TH1D*)f.Get("hist_alpha_error_final");
        TH1D* chi2_total_periter = (TH1D*)f.Get("chi2_total_periter");
        best_fit[i] = hist_alpha_result->GetBinContent(1);
        fit_err[i] = hist_alpha_error_final->GetBinContent(1);
        chi2[i]=chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX()-1);
        frac_err[i] = fit_err[i]/truth_alpha;

        hist_alpha->SetBinContent(i+1,best_fit[i]);
        hist_alpha->SetBinError(i+1,fit_err[i]);
    }

    TCanvas* c1 = new TCanvas();
    //hist_alpha->GetXaxis()->SetRangeUser(0,10000);
    hist_alpha->GetXaxis()->SetTitle("Number of mPMTs");
    hist_alpha->GetYaxis()->SetTitle("Fitted #alpha (cm)");
    hist_alpha->GetYaxis()->SetTitleOffset(1.5);
    hist_alpha->GetYaxis()->SetRangeUser(truth_alpha*0.9,truth_alpha*1.1);
    hist_alpha->Draw();
    TLine *line1 = new TLine(binEdges[0],truth_alpha,binEdges[nFiles],truth_alpha);
    line1->SetLineWidth(3);
    line1->SetLineStyle(9);
    line1->SetLineColor(kBlue);
    line1->Draw();
    auto tl = new TLegend(0.6,0.6,0.8,0.8);
    tl->AddEntry(hist_alpha,"Fitted #alpha","lep");
    tl->AddEntry(line1,"Truth #alpha","l");
    tl->SetBorderSize(0);
    tl->SetFillStyle(0);
    tl->Draw();
    c1->SaveAs("plots/diffuser4_400nm_nominal_fitted_alpha.pdf");

    TCanvas* c2 = new TCanvas();
    TGraph* gr_frac_err = new TGraph(nFiles,nPMT,frac_err);
    gr_frac_err->SetTitle("");
    gr_frac_err->GetXaxis()->SetTitle("Number of mPMTs");
    gr_frac_err->GetYaxis()->SetTitle("Fraction error of fitted #alpha");
    gr_frac_err->GetYaxis()->SetTitleOffset(1.5);
    gr_frac_err->Draw("A*");
    c2->SaveAs("plots/diffuser4_400nm_nominal_fitted_alpha_error.pdf");



//     /*double cut_value[] = {-948,-947.9,-947.8,-947.7,-947.6,
//                           -947.5,-947.4,-947.3,-947.2,-947.1,
//                           -947.0,-946.9,-946.8,-946.7,-946.6,
//                           -946.5,-946.4,-946.3,-946.2,-946.1,-946.0
//                           };
//     double best_fit[] = {10454.6,10478.1,10496.2,10512,10526.7,
//                          10542.4,10558,10570.4,10584,10598.1,
//                          10610.8,10625.4,10639.5,10654.8,10669.2,
//                          10682.9,10695.6,10709.6,10723.3,10737.7,10751.1
//                          };
//     double fit_err[] = {26.1755,26.1398,26.1172,26.1136,26.1244,
//                         26.5384,26.202,26.4786,26.5114,26.3332,
//                         27.3065,27.1043,26.4991,27.0526,26.6243,
//                         26.6814,26.7345,26.793,27.4835,27.4781,27.442
//                         };
//     double chi2[] = {420.647905178950623,426.606803253819919,431.984865867160124,431.685334245097181,435.734122513476791,
//                      440.89749561223698,442.582374246697441,446.358975709943479,448.264620029241769,451.697346738721421,
//                      452.896750629614985,453.62356646664449,459.70399555233513,461.191868323573715,462.737524831582959,
//                      464.565971258331047,464.42846009481633,470.846379758606076,472.230008616049247,473.718318824266134,476.666893267659816
//                      };*/

//     double timetof_min = -948;
//     double timetof_dt = 0.1;
//     const int timetof_nfiles = 21;
//     double chi2[timetof_nfiles];

//     double truth_alpha = 17858.242;
//     TCanvas* c1 = new TCanvas();
//     TH1D* h_alpha = new TH1D("","",timetof_nfiles,timetof_min-timetof_dt/2.,timetof_min+(timetof_nfiles-0.5)*timetof_dt);
//     for (int i=1;i<=h_alpha->GetNbinsX();i++){
//         double timetofcut = -948+(i-1)*0.1;
//         TFile f(Form("fitoutput_diffuser4_400nm_x2scattering_mPMT_fineAngular_timetof_%1.1f.root",timetofcut));
//         TH1D* hist_alpha_result = (TH1D*)f.Get("hist_alpha_result");
//         TH1D* hist_alpha_error_final = (TH1D*)f.Get("hist_alpha_error_final");
//         TH1D* chi2_total_periter = (TH1D*)f.Get("chi2_total_periter");
//         chi2[i-1]=chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX()-1);
//         h_alpha->SetBinContent(i,hist_alpha_result->GetBinContent(1));
//         h_alpha->SetBinError(i,hist_alpha_error_final->GetBinContent(1));
//         //h_alpha->SetBinContent(i,best_fit[i-1]);
//         //h_alpha->SetBinError(i,fit_err[i-1]);
//     }
//     h_alpha->GetXaxis()->SetTitle("timetof cut (ns)");
//     h_alpha->GetYaxis()->SetTitle("Fitted #alpha (cm)");
//     h_alpha->Draw();
//     TLine *line1 = new TLine(timetof_min-timetof_dt/2.,truth_alpha,timetof_min+(timetof_nfiles-0.5)*timetof_dt,truth_alpha);
//     line1->SetLineWidth(3);
//     line1->SetLineStyle(9);
//     line1->SetLineColor(kBlue);
//     line1->Draw();
//     TLatex latex;
//     latex.SetTextSize(0.06);
//     latex.SetTextColor(kBlue);
//     latex.DrawLatex(-948,truth_alpha*1.001,"Truth #alpha");

//     TCanvas* c2 = new TCanvas();
//     TH1D* h_chi2 = new TH1D("","",timetof_nfiles,timetof_min-timetof_dt/2.,timetof_min+(timetof_nfiles-0.5)*timetof_dt);
//     for (int i=1;i<=h_chi2->GetNbinsX();i++){
//         h_chi2->SetBinContent(i,chi2[i-1]);
//     }
//     h_chi2->GetXaxis()->SetTitle("timetof cut (ns)");
//     h_chi2->GetYaxis()->SetTitle("Post-fit #chi^{2}");
//     h_chi2->Draw();
// return;
//     TChain* fChain = new TChain("hitRate_pmtType1");
//     fChain->Add("../diffuser4_350nm_x2absorption/absorption_diffuser*.root");
//     //fChain->Add("../diffuser4_350nm_x2absorption/absorption_diffuser*.root");
//     //fChain->Add("/hepstore/kmtsui/hyperk/out.root");

//     std::cout<<"fChain->GetFile()->GetName() = "<<fChain->GetFile()->GetName()<<std::endl;
//     TFile* f = new TFile(fChain->GetFile()->GetName());
//     TTree* t = (TTree*)f->Get("pmt_type1");
//     double R, costh, cosths, omega, phim;
//     t->SetBranchAddress("R",&R);
//     t->SetBranchAddress("costh",&costh);
//     t->SetBranchAddress("cosths",&cosths);
//     t->SetBranchAddress("omega",&omega);
//     t->SetBranchAddress("phim",&phim);
//     int nPMTs = t->GetEntries();
//     double cosths_array[nPMTs], costh_array[nPMTs], R_array[nPMTs];
//     for (int i=0;i<nPMTs;i++){
//         t->GetEntry(i);
//         cosths_array[i]=cosths;
//         costh_array[i]=costh;
//         R_array[i]=R;
//     }


//     double nPE, timetof;
//     int PMT_id;
//     fChain->SetBranchAddress("nPE",&nPE);
//     fChain->SetBranchAddress("timetof",&timetof);
//     fChain->SetBranchAddress("PMT_id",&PMT_id);
//     TH1D* h_Timetof = new TH1D("","",1000,-952,-942);
//     TH1D* h_R_all = new TH1D("","",100,1000,9000);
//     TH1D* h_R_cut = (TH1D*) h_R_all->Clone();
//     for (unsigned long int i=0;i<fChain->GetEntries();i++){
//         fChain->GetEntry(i);
//         if (costh_array[PMT_id]>0.5 && cosths_array[PMT_id]>0.767 )
//         //if (costh_array[PMT_id]>0.5 && cosths_array[PMT_id]<0.766 )
//         {
//             h_Timetof->Fill(timetof,nPE);
//             if (timetof<-946.5) h_R_all->Fill(R_array[PMT_id],nPE);
//             if (timetof<-947.5) h_R_cut->Fill(R_array[PMT_id],nPE);
//         }
//     }
//     TCanvas* c3 = new TCanvas();
//     h_Timetof->GetXaxis()->SetTitle("timetof (ns)");
//     h_Timetof->Draw("hist");
//     auto fitresult = h_Timetof->Fit("gaus","S","",-951,-949);
//     TF1 *fa = new TF1("fa","[0]*TMath::Exp(-0.5*(x-[1])*(x-[1])/[2]/[2])",-952,-942);
//     fa->SetParameters(fitresult->Parameter(0),fitresult->Parameter(1),fitresult->Parameter(2));
//     fa->SetLineColor(kRed);
//     fa->Draw("same");
//     double integral_timewindow = 0;
//     double integral_all = 0;
//     for (int i = 1; i<= h_Timetof->GetNbinsX();i++) {
//         if (h_Timetof->GetBinCenter(i)<-947) {integral_timewindow+=h_Timetof->GetBinContent(i);integral_all=integral_timewindow;}
//         else {
//             integral_all+=h_Timetof->GetBinContent(i);
//             if (integral_all>=2*integral_timewindow) {
//                 std::cout<<"Cut at timetof = "<<h_Timetof->GetBinCenter(i)<<std::endl;
//                 break;
//             }
//         }
//     }
//     std::cout<<"Direct photon time window integral = "<<h_Timetof->Integral()<<std::endl;

//     TCanvas* c4 = new TCanvas();
//     h_R_cut->Divide(h_R_all);
//     h_R_cut->Draw();
}