double PoissonLLH (double mc, double w2, double data)
{

        // Standard Poisson LLH.
        double chi2 = 0.0;
        if(mc > 0.0)
        {
            chi2 = 2 * (mc - data);
            if(data > 0.0)
                chi2 += 2 * data * std::log(data / mc);
        }

        return (chi2 >= 0.0) ? chi2 : 0.0;

}

void CalcChi2_mPMT(){
    gStyle->SetOptStat(0);

    TFile* f = new TFile("fitoutput_mPMT_947.root");
    TH1D* evhist_sam1_data = (TH1D*)f->Get("evhist_sam0_data");
    TH1D* evhist_sam1_pred_final = (TH1D*)f->Get("evhist_sam0_pred_final");

    double* exp_w  = evhist_sam1_pred_final->GetArray();
    double* exp_w2 = evhist_sam1_pred_final->GetSumw2()->GetArray();
    double* data   = evhist_sam1_data->GetArray();

    TFile* fPMT = new TFile("../diffuser4_400nm_nominal/absorption_diffuser4_0.root");
    TTree* t = (TTree*)fPMT->Get("pmt_type1");
    double R, costh, cosths, omega, phim;
    int mPMT_id;
    t->SetBranchAddress("R",&R);
    t->SetBranchAddress("costh",&costh);
    t->SetBranchAddress("cosths",&cosths);
    t->SetBranchAddress("omega",&omega);
    t->SetBranchAddress("phim",&phim);
    t->SetBranchAddress("mPMT_id",&mPMT_id);
    int nPMTs = t->GetEntries();
    TH1D* h_R = new TH1D("","",100,1000,9000);
    TH1D* h_R_all = (TH1D*) h_R->Clone();
    TH1D* h_costh = new TH1D("","",100,0.5,1);
    TH1D* h_costh_all = (TH1D*) h_costh->Clone();
    TH2D* h_R_costh = new TH2D("","",20,0.5,1,20,1000,9000);
    TH2D* h_R_costh_all = (TH2D*)h_R_costh->Clone();
    TH1D* h_mPMTid = new TH1D("","",19,0,19);
    TH1D* h_mPMTid_all = (TH1D*) h_mPMTid->Clone();

    double chi2 = 0.0;
    int nbins = evhist_sam1_data->GetNbinsX();
    int nonzerobins = 0;
    TH1D* h_chi2 = new TH1D("","",nbins,0,nbins);
    TH1D* h_chi2_dist = new TH1D("","",300,0,30);
    for(unsigned int i = 1; i <= nbins; ++i)
    {
        double val = PoissonLLH(exp_w[i], exp_w2[i], data[i]);
        chi2 += val;
        if (data[i]>0) {
            nonzerobins++;

            t->GetEntry(i-1);
            h_R_all->Fill(R);
            h_costh_all->Fill(costh);
            h_R_costh_all->Fill(costh,R);
            h_mPMTid_all->Fill(mPMT_id+0.5);

            if (val>5) {
                h_R->Fill(R);
                h_costh->Fill(costh);
                h_R_costh->Fill(costh,R);
                h_mPMTid->Fill(mPMT_id+0.5);
            }
        }
        h_chi2->SetBinContent(i,val);
        h_chi2_dist->Fill(val);

    }

    TCanvas* c1 = new TCanvas();
    h_chi2->Draw("hist");

    TCanvas* c11 = new TCanvas();
    h_chi2_dist->Draw("hist");

    TCanvas* c_R = new TCanvas();
    h_R->Divide(h_R_all);
    h_R->Draw("hist");

    TCanvas* c_costh = new TCanvas();
    h_costh->Divide(h_costh_all);
    h_costh->Draw("hist");

    TCanvas* c_mPMTid = new TCanvas();
    h_mPMTid->Divide(h_mPMTid_all);
    h_mPMTid->Draw("hist");
    
    TCanvas* c_R_costh = new TCanvas();
    h_R_costh->Divide(h_R_costh_all);
    h_R_costh->Draw("colz");

    std::cout<<"chi2 = "<<chi2<<", nonzerobins = "<<nonzerobins<<std::endl;


return;
    TCanvas* c2 = new TCanvas();
    TH1D* hist_mPMT_angular_result = (TH1D*)f->Get("hist_mPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final = (TH1D*)f->Get("hist_mPMT_angular_error_final");
    int nbins_costh = hist_mPMT_angular_result->GetNbinsX();
    double costh_min = 0.5, costh_max = 1.0; 
    TH1D* h_angular = new TH1D("","",nbins_costh,costh_min,costh_max);
    for (int i = 1; i<=nbins_costh; i++ ){
        h_angular->SetBinContent(i,hist_mPMT_angular_result->GetBinContent(i));
        h_angular->SetBinError(i,hist_mPMT_angular_error_final->GetBinContent(i));
    }
    h_angular->GetXaxis()->SetTitle("cos#theta");
    h_angular->GetYaxis()->SetTitle("A(cos#theta)");
    h_angular->Draw();

    gStyle->SetPaintTextFormat("0.2f");
    TCanvas* c3 = new TCanvas();
    TMatrixDSym* res_cor_matrix = (TMatrixDSym*)f->Get("res_cor_matrix");
    res_cor_matrix->Draw("colz text");
    TH2D *h_postfit_cor = (TH2D*) c3->GetPrimitive("TMatrixDBase");
    h_postfit_cor->GetXaxis()->SetTitle("Parameter Bin");
    h_postfit_cor->GetYaxis()->SetTitle("Parameter Bin");
    h_postfit_cor->GetZaxis()->SetRangeUser(-1,1);

}

void CalcChi2_BnLPMT(){
    gStyle->SetOptStat(0);

    TFile* f = new TFile("fitoutput_BnLPMT_test.root");
    TH1D* evhist_sam1_data = (TH1D*)f->Get("evhist_sam0_data");
    TH1D* evhist_sam1_pred_final = (TH1D*)f->Get("evhist_sam0_pred_final");

    double* exp_w  = evhist_sam1_pred_final->GetArray();
    double* exp_w2 = evhist_sam1_pred_final->GetSumw2()->GetArray();
    double* data   = evhist_sam1_data->GetArray();

    double chi2 = 0.0;
    int nbins = evhist_sam1_data->GetNbinsX();
    int nonzerobins = 0;
    TH1D* h_chi2 = new TH1D("","",nbins,0,nbins);
    for(unsigned int i = 1; i <= nbins; ++i)
    {
        double val = PoissonLLH(exp_w[i], exp_w2[i], data[i]);
        chi2 += val;
        if (data[i]>0) nonzerobins++;
        h_chi2->SetBinContent(i,val);
    }

    TCanvas* c1 = new TCanvas();
    h_chi2->Draw("hist");
    
    std::cout<<"chi2 = "<<chi2<<", nonzerobins = "<<nonzerobins<<std::endl;

    TCanvas* c2 = new TCanvas();
    TH1D* hist_mPMT_angular_result = (TH1D*)f->Get("hist_BnLPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final = (TH1D*)f->Get("hist_BnLPMT_angular_error_final");
    int nbins_costh = hist_mPMT_angular_result->GetNbinsX();
    double costh_min = 0.5, costh_max = 1.0; 
    TH1D* h_angular = new TH1D("","",nbins_costh,costh_min,costh_max);
    for (int i = 1; i<=nbins_costh; i++ ){
        h_angular->SetBinContent(i,hist_mPMT_angular_result->GetBinContent(i));
        h_angular->SetBinError(i,hist_mPMT_angular_error_final->GetBinContent(i));
    }
    h_angular->GetXaxis()->SetTitle("cos#theta");
    h_angular->GetYaxis()->SetTitle("A(cos#theta)");
    h_angular->Draw();

    gStyle->SetPaintTextFormat("0.2f");
    TCanvas* c3 = new TCanvas();
    TMatrixDSym* res_cor_matrix = (TMatrixDSym*)f->Get("res_cor_matrix");
    res_cor_matrix->Draw("colz text");
    TH2D *h_postfit_cor = (TH2D*) c3->GetPrimitive("TMatrixDBase");
    h_postfit_cor->GetXaxis()->SetTitle("Parameter Bin");
    h_postfit_cor->GetYaxis()->SetTitle("Parameter Bin");
    h_postfit_cor->GetZaxis()->SetRangeUser(-1,1);

}

void CalcChi2_combined(){
    gStyle->SetOptStat(0);

    TFile* f = new TFile("fitoutput_combined.root");

    TH1D* evhist_sam1_data = (TH1D*)f->Get("evhist_sam0_data");
    TH1D* evhist_sam1_pred_final = (TH1D*)f->Get("evhist_sam0_pred_final");

    double* exp_w  = evhist_sam1_pred_final->GetArray();
    double* exp_w2 = evhist_sam1_pred_final->GetSumw2()->GetArray();
    double* data   = evhist_sam1_data->GetArray();

    double chi2 = 0.0;
    int nbins = evhist_sam1_data->GetNbinsX();
    int nonzerobins = 0;
    TH1D* h_chi2_BnL = new TH1D("","",nbins,0,nbins);
    for(unsigned int i = 1; i <= nbins; ++i)
    {
        double val = PoissonLLH(exp_w[i], exp_w2[i], data[i]);
        chi2 += val;
        if (data[i]>0) nonzerobins++;
        h_chi2_BnL->SetBinContent(i,val);
    }

    TCanvas* c1_BnL = new TCanvas();
    h_chi2_BnL->Draw("hist");
    
    std::cout<<"BnL chi2 = "<<chi2<<", nonzerobins = "<<nonzerobins<<std::endl;

    evhist_sam1_data = (TH1D*)f->Get("evhist_sam1_data");
    evhist_sam1_pred_final = (TH1D*)f->Get("evhist_sam1_pred_final");

    exp_w  = evhist_sam1_pred_final->GetArray();
    exp_w2 = evhist_sam1_pred_final->GetSumw2()->GetArray();
    data   = evhist_sam1_data->GetArray();

    chi2 = 0.0;
    nbins = evhist_sam1_data->GetNbinsX();
    nonzerobins = 0;
    TH1D* h_chi2_mPMT = new TH1D("","",nbins,0,nbins);
    for(unsigned int i = 1; i <= nbins; ++i)
    {
        double val = PoissonLLH(exp_w[i], exp_w2[i], data[i]);
        chi2 += val;
        if (data[i]>0) nonzerobins++;
        h_chi2_mPMT->SetBinContent(i,val);
    }

    TCanvas* c1_mPMT = new TCanvas();
    h_chi2_mPMT->Draw("hist");
    
    std::cout<<"mPMT chi2 = "<<chi2<<", nonzerobins = "<<nonzerobins<<std::endl;

    TCanvas* c2_BnL = new TCanvas();
    TH1D* hist_BnLPMT_angular_result = (TH1D*)f->Get("hist_BnLPMT_angular_result");
    TH1D* hist_BnLPMT_angular_error_final = (TH1D*)f->Get("hist_BnLPMT_angular_error_final");
    int nbins_costh = hist_BnLPMT_angular_result->GetNbinsX();
    double costh_min = 0.5; double costh_max = 1.0; 
    TH1D* h_angular_BnL = new TH1D("","",nbins_costh,costh_min,costh_max);
    for (int i = 1; i<=nbins_costh; i++ ){
        h_angular_BnL->SetBinContent(i,hist_BnLPMT_angular_result->GetBinContent(i));
        h_angular_BnL->SetBinError(i,hist_BnLPMT_angular_error_final->GetBinContent(i));
    }
    h_angular_BnL->GetXaxis()->SetTitle("cos#theta");
    h_angular_BnL->GetYaxis()->SetTitle("A_{j}(cos#theta)");
    h_angular_BnL->Draw();

    TCanvas* c2_mPMT = new TCanvas();
    TH1D* hist_mPMT_angular_result = (TH1D*)f->Get("hist_mPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final = (TH1D*)f->Get("hist_mPMT_angular_error_final");
    nbins_costh = hist_mPMT_angular_result->GetNbinsX();
    costh_min = 0.5; costh_max = 1.0; 
    TH1D* h_angular_mPMT = new TH1D("","",nbins_costh,costh_min,costh_max);
    for (int i = 1; i<=nbins_costh; i++ ){
        h_angular_mPMT->SetBinContent(i,hist_mPMT_angular_result->GetBinContent(i));
        h_angular_mPMT->SetBinError(i,hist_mPMT_angular_error_final->GetBinContent(i));
    }
    h_angular_mPMT->GetXaxis()->SetTitle("cos#theta");
    h_angular_mPMT->GetYaxis()->SetTitle("A_{i}(cos#theta)");
    h_angular_mPMT->Draw();

    gStyle->SetPaintTextFormat("0.2f");
    gStyle->SetTextSize(0.5);
    //gStyle->SetPaperSize(500,500);
    TCanvas* c3 = new TCanvas("","",2000,900);
    TMatrixDSym* res_cor_matrix = (TMatrixDSym*)f->Get("res_cor_matrix");
    res_cor_matrix->Draw("colz text");
    TH2D *h_postfit_cor = (TH2D*) c3->GetPrimitive("TMatrixDBase");
    h_postfit_cor->GetXaxis()->SetTitle("Parameter Bin");
    h_postfit_cor->GetYaxis()->SetTitle("Parameter Bin");
    h_postfit_cor->GetZaxis()->SetRangeUser(-1,1);
    c3->SaveAs("cor_matrix.pdf");

}

void CalcChi2(){
    CalcChi2_combined();
    //CalcChi2_mPMT();
}