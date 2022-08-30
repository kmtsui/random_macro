std::vector<int> pol_orders; // order of polynomial in each piece 
std::vector<double> pol_range; // applicable range for each polynomial
std::vector<double> CalcPol(const double* par, std::vector<double> costh_array)
{
    // Derive the coefficients for each polynomial
    std::vector<std::vector<double>> pol_coeff;
    std::vector<double> pol_p0, pol_p1;
    int par_index = 0;
    for (int i=0;i<pol_orders.size();i++)
    {
        std::vector<double> coeff;
        if (i==0) // for the first polynomial, the coefficients are unconstrained
        {
            for (int j=0;j<=pol_orders[i];j++)
            {
                coeff.push_back(par[j]);
                par_index++;
            } 
        }
        else // for others, we need to match the 0-th and 1-st order derivatives 
        {
            coeff.push_back(pol_p0[i-1]);
            coeff.push_back(pol_p1[i-1]);
            for (int j=2;j<=pol_orders[i];j++)
            {
                coeff.push_back(par[par_index]);
                par_index++;
            } 
        }
        pol_coeff.emplace_back(coeff);
        // std::cout<<"pol "<<i<<": ";
        // for (int j=0;j<pol_coeff[i].size();j++)
        //     std::cout<<pol_coeff[i][j]<<" ";
        // std::cout<<std::endl;
        double p0 = 0;
        double p1 = 0;
        for (int j=0;j<=pol_orders[i];j++) // store the 0-th and 1-st order derivatives at end-point as boundary conditions
        {
            p0 += coeff[j]*TMath::Power(pol_range[i+1]-pol_range[i],j);
            p1 += coeff[j]*j*TMath::Power(pol_range[i+1]-pol_range[i],j-1);
        } 
        pol_p0.push_back(p0);
        pol_p1.push_back(p1);
    }

    std::vector<double> angular_response;
    for (int k=0;k<costh_array.size();k++) // actually calculate the angular response
    {
        double costh = costh_array[k];
        double val = 0;
        for (int i=0;i<pol_orders.size();i++)
        {
            if (costh>=pol_range[i] && costh<pol_range[i+1])
            {
                //std::cout<<"Using pol"<<i<<std::endl;
                for (int j=0;j<=pol_orders[i];j++)
                {
                    val += pol_coeff[i][j]*TMath::Power(costh-pol_range[i],j);
                }
            }
        }
        //std::cout<<"CalcPol:: costh = "<<costh<<", val = "<<val<<std::endl;
        angular_response.push_back(val);
    }
        

    return angular_response;
}

void plot_poly_attenz(std::string file1, std::string file2, int pmttype, std::vector<int> orders, std::vector<double> ranges)
{
    gStyle->SetOptStat(0);

    pol_orders = orders; pol_range = ranges;

    TFile* fold = new TFile(file1.c_str());
    TFile* fnew = new TFile(file2.c_str());
    TH1D* hist_mPMT_angular_result_old = pmttype==0 ? (TH1D*)fold->Get("hist_BnLPMT_angular_result") : (TH1D*)fold->Get("hist_mPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final_old = pmttype==0 ? (TH1D*)fold->Get("hist_BnLPMT_angular_error_final") : (TH1D*)fold->Get("hist_mPMT_angular_error_final");
    TH1D* hist_mPMT_angular_pol_result = pmttype==0 ? (TH1D*)fnew->Get("hist_BnLPMT_angular_pol_result") : (TH1D*)fnew->Get("hist_mPMT_angular_pol_result");
    TH1D* hist_mPMT_angular_error_final_new = pmttype==0 ? (TH1D*)fnew->Get("hist_BnLPMT_angular_pol_error_final") : (TH1D*)fnew->Get("hist_mPMT_angular_pol_error_final");

    double truth_alpha = 10648;
    double fitted_alpha = ((TH1D*)fnew->Get("hist_alphaZ_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)fnew->Get("hist_alphaZ_error_final"))->GetBinContent(1);
    double truth_slope = 0.15;
    double fitted_slope = ((TH1D*)fnew->Get("hist_alphaZ_result"))->GetBinContent(2);
    double error_slope = ((TH1D*)fnew->Get("hist_alphaZ_error_final"))->GetBinContent(2);

    TVectorT<double>* chi2_tuple_postfit = (TVectorT<double>*)fnew->Get("chi2_tuple_postfit");
    double chi2 = (*chi2_tuple_postfit)[3];
    double chi2_stat = (*chi2_tuple_postfit)[0];
    double chi2_syst = (*chi2_tuple_postfit)[1];
    TH1D* evhist_sam0_pred_final = (TH1D*)fnew->Get("evhist_sam0_pred_final");
    int ndof = 0; //((TTree*)fnew->Get("PMTTree"))->GetEntries();
    for (int i=1;i<=evhist_sam0_pred_final->GetNbinsX();i++)
        if (evhist_sam0_pred_final->GetBinContent(i)>0) ndof++;

    TH1D* hist_angular_old = new TH1D("","",40,0,1);
    TH1D* hist_angular_new = new TH1D("","",40,0,1);
    std::vector<double> costh_array;
    for (int i=1;i<=hist_mPMT_angular_result_old->GetNbinsX();i++)
    {
        double val = hist_mPMT_angular_result_old->GetBinContent(i);
        if (hist_mPMT_angular_error_final_old->GetBinContent(i)<1.e-9) continue;
        hist_angular_old->SetBinContent(i,val);
        //hist_angular_old->SetBinError(i,hist_mPMT_angular_error_final_old->GetBinContent(i));
        // if (hist_mPMT_angular_error_final_new->GetBinContent(i)<1.e-9) continue;
        // hist_angular_new->SetBinContent(i,hist_mPMT_angular_result_new->GetBinContent(i));
        //hist_angular_new->SetBinError(i,hist_mPMT_angular_error_final_new->GetBinContent(i));

        costh_array.push_back(hist_angular_old->GetBinCenter(i));
    }

    double poly_params[100];
    for (int i=1;i<=hist_mPMT_angular_pol_result->GetNbinsX();i++)
    {
        poly_params[i-1] = hist_mPMT_angular_pol_result->GetBinContent(i);
    }
    std::vector<double> angular_response = CalcPol(poly_params,costh_array);
    for (int i=1;i<=hist_angular_new->GetNbinsX();i++)
    {
        hist_angular_new->SetBinContent(i,angular_response[i-1]);
    }

    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    hist_angular_old->GetXaxis()->SetTitle("cos#theta");
    hist_angular_old->GetYaxis()->SetTitle("A(cos#theta)");
    hist_angular_old->Draw();
    hist_angular_new->SetLineColor(kRed);
    hist_angular_new->Draw("same");
    TLegend* legend = new TLegend(0.15,0.55,0.55,0.85);
    legend->AddEntry(hist_angular_old,"Nominal MC","l");
    legend->AddEntry(hist_angular_new,"Fake data","l");
    legend->AddEntry((TObject*)0,Form("#chi^{2}_{stat}/bin=%2.0f/%i, #chi^{2}_{syst}=%2.0f",chi2_stat,ndof,chi2_syst),"");
    legend->AddEntry((TObject*)0,Form("Truth L_{#alpha} = %3.0f cm",truth_alpha),"");
    legend->AddEntry((TObject*)0,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha),"");
    legend->AddEntry((TObject*)0,Form("Truth slope = %0.2f /cm",truth_slope),"");
    legend->AddEntry((TObject*)0,Form("Fitted slope = %0.2f +/- %0.2f /cm",fitted_slope,error_slope),"");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();
    c1->SaveAs("test.pdf");

    hist_angular_new->Divide(hist_angular_old);
    hist_angular_new->Draw();
    c1->SaveAs("test1.pdf");

}

void plot_poly(std::string file1, std::string file2, int pmttype, std::vector<int> orders, std::vector<double> ranges)
{
    gStyle->SetOptStat(0);

    pol_orders = orders; pol_range = ranges;

    TFile* fold = new TFile(file1.c_str());
    TFile* fnew = new TFile(file2.c_str());
    TH1D* hist_mPMT_angular_result_old = pmttype==0 ? (TH1D*)fold->Get("hist_BnLPMT_angular_result") : (TH1D*)fold->Get("hist_mPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final_old = pmttype==0 ? (TH1D*)fold->Get("hist_BnLPMT_angular_error_final") : (TH1D*)fold->Get("hist_mPMT_angular_error_final");
    TH1D* hist_mPMT_angular_pol_result = pmttype==0 ? (TH1D*)fnew->Get("hist_BnLPMT_angular_pol_result") : (TH1D*)fnew->Get("hist_mPMT_angular_pol_result");
    TH1D* hist_mPMT_angular_error_final_new = pmttype==0 ? (TH1D*)fnew->Get("hist_BnLPMT_angular_pol_error_final") : (TH1D*)fnew->Get("hist_mPMT_angular_pol_error_final");

    double fitted_alpha = ((TH1D*)fnew->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)fnew->Get("hist_alpha_error_final"))->GetBinContent(1);

    TH1D* chi2_total_periter = (TH1D*)fnew->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX()-1);
    int ndof = ((TTree*)fnew->Get("PMTTree"))->GetEntries()-1;
    ndof -= pmttype==0 ? 0 : 12;

    TH1D* hist_angular_old = new TH1D("","",40,0,1);
    TH1D* hist_angular_new = new TH1D("","",40,0,1);
    std::vector<double> costh_array;
    for (int i=1;i<=hist_mPMT_angular_result_old->GetNbinsX();i++)
    {
        double val = hist_mPMT_angular_result_old->GetBinContent(i);
        if (hist_mPMT_angular_error_final_old->GetBinContent(i)<1.e-9) continue;
        hist_angular_old->SetBinContent(i,val);
        //hist_angular_old->SetBinError(i,hist_mPMT_angular_error_final_old->GetBinContent(i));
        // if (hist_mPMT_angular_error_final_new->GetBinContent(i)<1.e-9) continue;
        // hist_angular_new->SetBinContent(i,hist_mPMT_angular_result_new->GetBinContent(i));
        //hist_angular_new->SetBinError(i,hist_mPMT_angular_error_final_new->GetBinContent(i));
        ndof--;

        costh_array.push_back(hist_angular_old->GetBinCenter(i));
    }

    double poly_params[100];
    for (int i=1;i<=hist_mPMT_angular_pol_result->GetNbinsX();i++)
    {
        poly_params[i-1] = hist_mPMT_angular_pol_result->GetBinContent(i);
    }
    std::vector<double> angular_response = CalcPol(poly_params,costh_array);
    for (int i=1;i<=hist_angular_new->GetNbinsX();i++)
    {
        hist_angular_new->SetBinContent(i,angular_response[i-1]);
    }

    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    hist_angular_old->GetXaxis()->SetTitle("cos#theta");
    hist_angular_old->GetYaxis()->SetTitle("A(cos#theta)");
    hist_angular_old->Draw();
    hist_angular_new->SetLineColor(kRed);
    hist_angular_new->Draw("same");
    TLegend* legend = new TLegend(0.15,0.6,0.55,0.85);
    legend->AddEntry(hist_angular_old,"Nominal MC","l");
    legend->AddEntry(hist_angular_new,"Fake data","l");
    legend->AddEntry((TObject*)0,Form("#chi^{2}/dof=%2.0f/%i",chi2,ndof),"");
    legend->AddEntry((TObject*)0,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha),"");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();
    c1->SaveAs("test.pdf");

    hist_angular_new->Divide(hist_angular_old);
    hist_angular_new->Draw();
    c1->SaveAs("test1.pdf");

}

void plot(std::string file1, std::string file2, int pmttype = 0)
{
    gStyle->SetOptStat(0);

    TFile* fold = new TFile(file1.c_str());
    TFile* fnew = new TFile(file2.c_str());
    TH1D* hist_mPMT_angular_result_old = pmttype==0 ? (TH1D*)fold->Get("hist_BnLPMT_angular_result") : (TH1D*)fold->Get("hist_mPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final_old = pmttype==0 ? (TH1D*)fold->Get("hist_BnLPMT_angular_error_final") : (TH1D*)fold->Get("hist_mPMT_angular_error_final");
    TH1D* hist_mPMT_angular_result_new = pmttype==0 ? (TH1D*)fnew->Get("hist_BnLPMT_angular_result") : (TH1D*)fnew->Get("hist_mPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final_new = pmttype==0 ? (TH1D*)fnew->Get("hist_BnLPMT_angular_error_final") : (TH1D*)fnew->Get("hist_mPMT_angular_error_final");

    double fitted_alpha = ((TH1D*)fnew->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)fnew->Get("hist_alpha_error_final"))->GetBinContent(1);

    TH1D* chi2_total_periter = (TH1D*)fnew->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX()-1);
    int ndof = ((TTree*)fnew->Get("PMTTree"))->GetEntries()-1;
    ndof -= pmttype==0 ? 0 : 12;

    TH1D* hist_angular_old = new TH1D("","",40,0,1);
    TH1D* hist_angular_new = new TH1D("","",40,0,1);
    for (int i=1;i<=hist_mPMT_angular_result_old->GetNbinsX();i++)
    {
        double val = hist_mPMT_angular_result_old->GetBinContent(i);
        if (hist_mPMT_angular_error_final_old->GetBinContent(i)<1.e-9) continue;
        hist_angular_old->SetBinContent(i,val);
        //hist_angular_old->SetBinError(i,hist_mPMT_angular_error_final_old->GetBinContent(i));
        if (hist_mPMT_angular_error_final_new->GetBinContent(i)<1.e-9) continue;
        hist_angular_new->SetBinContent(i,hist_mPMT_angular_result_new->GetBinContent(i));
        //hist_angular_new->SetBinError(i,hist_mPMT_angular_error_final_new->GetBinContent(i));
        ndof--;
    }
    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    hist_angular_old->GetXaxis()->SetTitle("cos#theta");
    hist_angular_old->GetYaxis()->SetTitle("A(cos#theta)");
    hist_angular_old->Draw();
    hist_angular_new->SetLineColor(kRed);
    hist_angular_new->Draw("same");
    TLegend* legend = new TLegend(0.15,0.6,0.55,0.85);
    legend->AddEntry(hist_angular_old,"Nominal MC","l");
    legend->AddEntry(hist_angular_new,"Fake data","l");
    legend->AddEntry((TObject*)0,Form("#chi^{2}/dof=%2.0f/%i",chi2,ndof),"");
    legend->AddEntry((TObject*)0,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha),"");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();
    c1->SaveAs("test.pdf");

    hist_angular_new->Divide(hist_angular_old);
    hist_angular_new->Draw();
    c1->SaveAs("test1.pdf");

}

void plot_attenz(std::string file1, std::string file2, int pmttype = 0)
{
    gStyle->SetOptStat(0);

    TFile* fold = new TFile(file1.c_str());
    TFile* fnew = new TFile(file2.c_str());
    TH1D* hist_mPMT_angular_result_old = pmttype==0 ? (TH1D*)fold->Get("hist_BnLPMT_angular_result") : (TH1D*)fold->Get("hist_mPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final_old = pmttype==0 ? (TH1D*)fold->Get("hist_BnLPMT_angular_error_final") : (TH1D*)fold->Get("hist_mPMT_angular_error_final");
    TH1D* hist_mPMT_angular_result_new = pmttype==0 ? (TH1D*)fnew->Get("hist_BnLPMT_angular_result") : (TH1D*)fnew->Get("hist_mPMT_angular_result");
    TH1D* hist_mPMT_angular_error_final_new = pmttype==0 ? (TH1D*)fnew->Get("hist_BnLPMT_angular_error_final") : (TH1D*)fnew->Get("hist_mPMT_angular_error_final");

    double fitted_alpha = ((TH1D*)fnew->Get("hist_alphaZ_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)fnew->Get("hist_alphaZ_error_final"))->GetBinContent(1);
    double fitted_slope = ((TH1D*)fnew->Get("hist_alphaZ_result"))->GetBinContent(2);
    double error_slope = ((TH1D*)fnew->Get("hist_alphaZ_error_final"))->GetBinContent(2);

    TH1D* chi2_total_periter = (TH1D*)fnew->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX()-1);
    int ndof = ((TTree*)fnew->Get("PMTTree"))->GetEntries()-2;
    ndof -= pmttype==0 ? 0 : 12;

    TH1D* hist_angular_old = new TH1D("","",40,0,1);
    TH1D* hist_angular_new = new TH1D("","",40,0,1);
    for (int i=1;i<=hist_mPMT_angular_result_old->GetNbinsX();i++)
    {
        double val = hist_mPMT_angular_result_old->GetBinContent(i);
        if (hist_mPMT_angular_error_final_old->GetBinContent(i)<1.e-9) continue;
        hist_angular_old->SetBinContent(i,val);
        //hist_angular_old->SetBinError(i,hist_mPMT_angular_error_final_old->GetBinContent(i));
        if (hist_mPMT_angular_error_final_new->GetBinContent(i)<1.e-9) continue;
        hist_angular_new->SetBinContent(i,hist_mPMT_angular_result_new->GetBinContent(i));
        //hist_angular_new->SetBinError(i,hist_mPMT_angular_error_final_new->GetBinContent(i));
        ndof--;
    }
    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    hist_angular_old->GetXaxis()->SetTitle("cos#theta");
    hist_angular_old->GetYaxis()->SetTitle("A(cos#theta)");
    hist_angular_old->Draw();
    hist_angular_new->SetLineColor(kRed);
    hist_angular_new->Draw("same");
    TLegend* legend = new TLegend(0.15,0.6,0.55,0.85);
    legend->AddEntry(hist_angular_old,"Nominal MC","l");
    legend->AddEntry(hist_angular_new,"Fake data","l");
    legend->AddEntry((TObject*)0,Form("#chi^{2}/dof=%2.0f/%i",chi2,ndof),"");
    legend->AddEntry((TObject*)0,Form("Fitted L_{#alpha} (z_{0}) = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha),"");
    legend->AddEntry((TObject*)0,Form("Fitted #beta = %3.3f +/- %3.3f",fitted_slope,error_slope),"");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();
    c1->SaveAs("test.pdf");

    hist_angular_new->Divide(hist_angular_old);
    hist_angular_new->Draw();
    c1->SaveAs("test1.pdf");

}

void plot_angular_compare()
{
    //plot("TN/diffusr4_400nm_nominal_combined_fullcosth_40degcosths_newMC.root","TN/diffusr20_400nm_zdep_BnL_upperbarrel.root",0);
    //plot("TN/diffusr4_400nm_nominal_combined_fullcosth_40degcosths_newMC.root","TN/diffusr4_400nm_nominal_mPMT_z0.15_20deg_10percentprior.root",1);
    //plot_attenz("TN/diffusr4_400nm_nominal_combined_fullcosth_40degcosths_newMC.root","test.root",0);
    //plot("systematic_study/source_4/1pc_mPMT1000_corrected.root","test.root",1);
    plot_poly_attenz("systematic_study/source_4/1pc_mPMT1000_corrected.root","test.root",1,std::vector<int>{3,3,3}, std::vector<double>{0,0.6,0.75,1.01});
}