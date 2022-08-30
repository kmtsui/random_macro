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


int counter;

TH1D* plot(std::string config, double sys_lo, double sys_hi)
{
    gStyle->SetOptStat(0);

    //std::string config = "fullring";
    //std::string config = "noring1";

    std::vector<double> costh_val{0.5,1};

    int nmPMT_min = 100;
    int step = 100;
    const int nstep = 40;
    TH1D* hist_resol_bias = new TH1D(config.c_str(),config.c_str(),nstep,nmPMT_min-step/2,nmPMT_min+(nstep-0.5)*step);
    TH1D* hist_resol = new TH1D("","",nstep,nmPMT_min-step/2,nmPMT_min+(nstep-0.5)*step);
    TH1D* hist_resol_bias_1ns = new TH1D("","",nstep,nmPMT_min-step/2,nmPMT_min+(nstep-0.5)*step);
    TH1D* hist_resol_bias_led = new TH1D("","",nstep,nmPMT_min-step/2,nmPMT_min+(nstep-0.5)*step);
    TH1D* hist_globalcc = new TH1D("","",nstep,nmPMT_min-step/2,nmPMT_min+(nstep-0.5)*step);
    double alpha_nom[nstep];
    for (int i=nstep-1;i>=0;i--)
    {
        int nmPMT = nmPMT_min+i*step;
        std::vector<double> fitValues;
        for (int j=1;j<=100;j++)
        {
            TFile* f = new TFile(Form("mPMT_eff_study/%s/diffuser4_400nm_nominal_mPMT_%i_fullring_%i.root",config.c_str(),nmPMT,j),"OPEN");
            TH1D* hist_alpha_result = (TH1D*)f->Get("hist_alpha_result");
            TH1D* hist_mPMT_angular_pol_result = (TH1D*)f->Get("hist_mPMT_angular_pol_result");
            if (hist_alpha_result)
            {
                const int npars = hist_mPMT_angular_pol_result->GetNbinsX();
                double pol_pars[npars];
                for (int k=1;k<=npars;k++)
                {
                    pol_pars[k-1] = hist_mPMT_angular_pol_result->GetBinContent(k);
                }
                std::vector<double> fit_val = CalcPol(pol_pars,costh_val);

                double val = fit_val[0]/fit_val[1];
                fitValues.push_back(val);
            }
            f->Close();
        }
        if (fitValues.size()>0)
        {
            double mean = TMath::Mean(fitValues.begin(),fitValues.end());
            double rms = TMath::RMS(fitValues.begin(),fitValues.end());
            hist_resol->SetBinContent(i+1,rms/mean); 
            if (rms/mean<0.0005)
                hist_resol->SetBinContent(i+1,0.005); 
            else if (rms/mean/hist_resol->GetBinContent(i+2)<1.5)
                hist_resol->SetBinContent(i+1,rms/mean); 
            else 
            {
                if (hist_resol->GetBinContent(i+2)>0.0005)
                    hist_resol->SetBinContent(i+1,hist_resol->GetBinContent(i+2)*1.1); 
                else 
                    hist_resol->SetBinContent(i+1,0.005); 
            }
        }
        else 
        {
            if (hist_resol->GetBinContent(i+2)>0.0005)
                hist_resol->SetBinContent(i+1,hist_resol->GetBinContent(i+2)*1.01); 
            else 
                hist_resol->SetBinContent(i+1,0.001); 
        }
        // double mean = TMath::Mean(fitValues.begin(),fitValues.end());
        // double rms = TMath::RMS(fitValues.begin(),fitValues.end());
        // // hist_resol->SetBinContent(i+1,mean/mean);
        // // hist_resol->SetBinError(i+1,rms/mean);
        // hist_resol->SetBinContent(i+1,rms/mean);
    }

    hist_resol->GetXaxis()->SetTitle("#mPMT");
    hist_resol->GetYaxis()->SetTitle("RMS/Mean of post-fit A(cos#theta=0.5)/A(cos#theta=1.0)");
    hist_resol->GetYaxis()->SetTitleOffset(1.5);
    hist_resol->GetYaxis()->SetRangeUser(0,0.05);
    
    hist_resol->SetLineColor(counter+1);
    hist_resol->SetLineStyle(counter+1);
    hist_resol->SetLineWidth(3);
    counter++;

    return hist_resol;
}

void plot_angular_eff_var()
{
    gStyle->SetOptStat(0);
    
    int orders[] = {2,3,3};
    double ranges[] = {0.0,0.6,0.75,1.01};
    pol_orders.assign(orders, orders+sizeof(orders)/sizeof(int));
    pol_range.assign(ranges, ranges+sizeof(ranges)/sizeof(double));

    counter = 0;
    
    // TH1D* h_full = plot("fullring",0.017,0.009);
    // TH1D* h_half = plot("halfring",0.024,0.011);
    // TH1D* h_noout = plot("noring1",0.046,0.012);
    // TH1D* h_nomid = plot("noring2",0.024,0.010);
    TH1D* h_full = plot("10pc/fullring",0.015,0.005);
    TH1D* h_half = plot("10pc/halfring",0.02,0.0075);
    TH1D* h_noout = plot("10pc/noring1",0.04,0.012);
    TH1D* h_nomid = plot("10pc/noring2",0.016,0.006);
    // TH1D* h_full = plot("2pc/fullring",0.015,0.005);
    // TH1D* h_half = plot("2pc/halfring",0.02,0.0075);
    // TH1D* h_noout = plot("2pc/noring1",0.04,0.012);
    // TH1D* h_nomid = plot("2pc/noring2",0.016,0.006);

    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    h_full->Draw("hist same");
    h_half->Draw("hist same");
    h_noout->Draw("hist same");
    h_nomid->Draw("hist same");
    TLegend* legend = new TLegend(0.45,0.65,0.75,0.9);
    legend->AddEntry(h_full,"Full ring","l");
    legend->AddEntry(h_half,"Half ring","l");
    legend->AddEntry(h_noout,"No outer","l");
    legend->AddEntry(h_nomid,"No middle","l");
    legend->Draw("same");
    c1->SaveAs("angular_eff_var_10pc.pdf");

    TFile* f = new TFile("mPMT_eff_study/10pc_angular_err.root","RECREATE");
    h_full->Write("fullring");
    h_half->Write("halfring");
    h_noout->Write("noring1");
    h_nomid->Write("noring2");
    f->Close();
}
