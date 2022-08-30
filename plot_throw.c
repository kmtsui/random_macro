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

void plot_throw()
{
    int nFiles = 1000;

    pol_orders = std::vector<int>{3,3,3};
    pol_range = std::vector<double>{0,0.6,0.75,1.01};

    double truth_alpha = 10650;
    TH1D* hAlpha = new TH1D("","",100,truth_alpha*0.8,truth_alpha*1.2);
    double truth_alphaZ = 0.15;
    TH1D* hAlphaZ = new TH1D("","",100,truth_alphaZ*0.,truth_alphaZ*2.0);
    TH1D* hChi2 = new TH1D("","",100,200,500);
    TH1D* hChi2_stat = new TH1D("","",100,200,500);
    TH1D* hChi2_syst = new TH1D("","",100,0,100);
    TH1D* hAngular = new TH1D("","",100,0,1);
    std::vector<std::vector<double>> vAngular(100,std::vector<double>{});
    std::vector<double> costh_array;
    for (int i=1;i<=hAngular->GetNbinsX();i++)
        costh_array.push_back(hAngular->GetBinCenter(i));
    
    TH1D* hPths = new TH1D("","",40,-40,0);
    TH1D* hVths = new TH1D("","",40,-40,0);
    TH1D* hEffz = new TH1D("","",40,-3300,3300);
    std::vector<std::vector<double>> vPths(40,std::vector<double>{});
    std::vector<std::vector<double>> vVths(40,std::vector<double>{});
    std::vector<std::vector<double>> vEffz(40,std::vector<double>{});

    for (int i=0;i<nFiles;i++)
    {
        TFile f(Form("all3/throw_all_toy%i.root",i));
        TH1D* hist_alphaZ_result = (TH1D*)f.Get("hist_alphaZ_result");
        hAlpha->Fill(hist_alphaZ_result->GetBinContent(1));
        hAlphaZ->Fill(hist_alphaZ_result->GetBinContent(2));
        TVectorD* chi2_tuple_postfit = (TVectorD*)f.Get("chi2_tuple_postfit");
        hChi2->Fill((*chi2_tuple_postfit)[3]);
        hChi2_stat->Fill((*chi2_tuple_postfit)[0]);
        hChi2_syst->Fill((*chi2_tuple_postfit)[1]);

        TH1D* hist_mPMT_angular_pol_result = (TH1D*)f.Get("hist_mPMT_angular_pol_result");
        double poly_params[100];
        for (int j=1;j<=hist_mPMT_angular_pol_result->GetNbinsX();j++)
        {
            poly_params[j-1] = hist_mPMT_angular_pol_result->GetBinContent(j);
        }
        std::vector<double> angular_response = CalcPol(poly_params,costh_array);
        for (int j=0;j<100;j++)
        {
            vAngular[j].push_back(angular_response[j]);
        }

        TH1D* hist_LED_theta_result = (TH1D*)f.Get("hist_LED_theta_result");
        TH1D* hist_LED_phi_result = (TH1D*)f.Get("hist_LED_phi_result");
        TH1D* hist_zpos_eff_result = (TH1D*)f.Get("hist_zpos_eff_result");
        for (int j=1;j<=40;j++)
        {
            vPths[j-1].push_back(hist_LED_theta_result->GetBinContent(j));
            vVths[j-1].push_back(hist_LED_phi_result->GetBinContent(j));
            vEffz[j-1].push_back(hist_zpos_eff_result->GetBinContent(j));
        }

        f.Close();
    }

    for (int j=0;j<100;j++)
    {
        double mean = TMath::Mean(vAngular[j].begin(),vAngular[j].end());
        double rms = TMath::RMS(vAngular[j].begin(),vAngular[j].end());
        hAngular->SetBinContent(j+1,mean);
        hAngular->SetBinError(j+1,rms);
    }

    TCanvas* c1 = new TCanvas();
    hAlpha->GetXaxis()->SetTitle("Fitted L_{#alpha} (cm)");
    hAlpha->Draw();
    c1->SaveAs("test.pdf");
    hAlphaZ->GetXaxis()->SetTitle("Fitted slope (cm^{-1})");
    hAlphaZ->Draw();
    c1->SaveAs("test1.pdf");

    hChi2->GetXaxis()->SetTitle("Total #chi^2");
    hChi2->Draw();
    c1->SaveAs("test2.pdf");
    hChi2_stat->GetXaxis()->SetTitle("#chi^2_{stat}");
    hChi2_stat->Draw();
    c1->SaveAs("test2.1.pdf");
    hChi2_syst->GetXaxis()->SetTitle("#chi^2_{syst}");
    hChi2_syst->Draw();
    c1->SaveAs("test2.2.pdf");

    gStyle->SetOptStat(0);
    c1->SetGridy();
    hAngular->GetXaxis()->SetTitle("cos#theta");
    hAngular->GetYaxis()->SetTitle("Fitted A(cos#theta)");
    hAngular->Draw("E0");
    c1->SaveAs("test3.pdf");

    TFile f("all4/throw_all.root");
    for (int j=0;j<40;j++)
    {
        double mean = TMath::Mean(vPths[j].begin(),vPths[j].end());
        double rms = TMath::RMS(vPths[j].begin(),vPths[j].end());
        hPths->SetBinContent(j+1,mean);
        hPths->SetBinError(j+1,rms);
    }
    TH1D* hist_LED_theta_prior = (TH1D*)f.Get("hist_LED_theta_prior");
    TH1D* hist_LED_theta_error_prior = (TH1D*)f.Get("hist_LED_theta_error_prior");
    hPths->GetXaxis()->SetTitle("-#theta_{s} (deg)");
    hPths->GetYaxis()->SetTitle("P(#theta_{s})");
    hPths->Draw("E0");
    for (int i=1;i<=hist_LED_theta_prior->GetNbinsX();i++)
    {
        TBox* box = new TBox(hPths->GetBinLowEdge(i),hist_LED_theta_prior->GetBinContent(i)-hist_LED_theta_error_prior->GetBinContent(i),
                             hPths->GetBinLowEdge(i+1),hist_LED_theta_prior->GetBinContent(i)+hist_LED_theta_error_prior->GetBinContent(i));
        box->SetFillColor(856);
        box->Draw("same");
    }
    hPths->Draw("E0 same");
    c1->SaveAs("test4.pdf");

    for (int j=0;j<40;j++)
    {
        double mean = TMath::Mean(vVths[j].begin(),vVths[j].end());
        double rms = TMath::RMS(vVths[j].begin(),vVths[j].end());
        hVths->SetBinContent(j+1,mean);
        hVths->SetBinError(j+1,rms);
    }
    TH1D* hist_LED_phi_prior = (TH1D*)f.Get("hist_LED_phi_prior");
    TH1D* hist_LED_phi_error_prior = (TH1D*)f.Get("hist_LED_phi_error_prior");
    hVths->GetXaxis()->SetTitle("-#theta_{s} (deg)");
    hVths->GetYaxis()->SetTitle("V(#theta_{s})");
    hVths->GetYaxis()->SetTitleOffset(1.2);
    hVths->Draw("E0");
    for (int i=1;i<=hist_LED_phi_prior->GetNbinsX();i++)
    {
        TBox* box = new TBox(hVths->GetBinLowEdge(i),hist_LED_phi_prior->GetBinContent(i)-hist_LED_phi_error_prior->GetBinContent(i),
                             hVths->GetBinLowEdge(i+1),hist_LED_phi_prior->GetBinContent(i)+hist_LED_phi_error_prior->GetBinContent(i));
        box->SetFillColor(856);
        box->Draw("same");
    }
    hVths->Draw("E0 same");
    c1->SaveAs("test5.pdf");

    for (int j=0;j<40;j++)
    {
        double mean = TMath::Mean(vEffz[j].begin(),vEffz[j].end());
        double rms = TMath::RMS(vEffz[j].begin(),vEffz[j].end());
        hEffz->SetBinContent(j+1,mean);
        hEffz->SetBinError(j+1,rms);
    }
    TH1D* hist_zpos_eff_prior = (TH1D*)f.Get("hist_zpos_eff_prior");
    TH1D* hist_zpos_eff_error_prior = (TH1D*)f.Get("hist_zpos_eff_error_prior");
    hEffz->GetXaxis()->SetTitle("z (cm)");
    hEffz->GetYaxis()->SetTitle("#epsilon(z)");
    hEffz->GetYaxis()->SetTitleOffset(1.2);
    hEffz->Draw("E0");
    for (int i=1;i<=hist_zpos_eff_prior->GetNbinsX();i++)
    {
        TBox* box = new TBox(hEffz->GetBinLowEdge(i),hist_zpos_eff_prior->GetBinContent(i)-hist_zpos_eff_error_prior->GetBinContent(i),
                             hEffz->GetBinLowEdge(i+1),hist_zpos_eff_prior->GetBinContent(i)+hist_zpos_eff_error_prior->GetBinContent(i));
        box->SetFillColor(856);
        box->Draw("same");
    }
    hEffz->Draw("E0 same");
    c1->SaveAs("test6.pdf");
}