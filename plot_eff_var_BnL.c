int counter;

TH1D* plot(std::string config, double sys_lo, double sys_hi)
{
    gStyle->SetOptStat(0);

    //std::string config = "fullring";
    //std::string config = "noring1";

    int nmPMT_min = 100;
    int step = 100;
    const int nstep = 40;
    TH1D* hist_resol_bias = new TH1D(config.c_str(),config.c_str(),nstep,nmPMT_min-step/2,nmPMT_min+(nstep-0.5)*step);
    TH1D* hist_resol = new TH1D("","",100,10800*0.8,10800*1.2);
    TH1D* hist_resol_bias_1ns = new TH1D("","",nstep,nmPMT_min-step/2,nmPMT_min+(nstep-0.5)*step);
    TH1D* hist_resol_bias_led = new TH1D("","",nstep,nmPMT_min-step/2,nmPMT_min+(nstep-0.5)*step);
    TH1D* hist_globalcc = new TH1D("","",nstep,nmPMT_min-step/2,nmPMT_min+(nstep-0.5)*step);
    double alpha_nom[nstep];
    std::vector<double> fitValues;
    for (int i=0;i<100;i++)
    {
        TFile* f = new TFile(Form("%stoy%i.root",config.c_str(),i),"OPEN");
        TH1D* hist_alpha_result = (TH1D*)f->Get("hist_alpha_result");
        if (hist_alpha_result)
        {
            double val = hist_alpha_result->GetBinContent(1)/10800;
            if (val>0.8&&val<1.2) fitValues.push_back(val);
            hist_resol->Fill(hist_alpha_result->GetBinContent(1));
        }
        f->Close();
    }
    if (fitValues.size()>0)
    {
        double mean = TMath::Mean(fitValues.begin(),fitValues.end());
        double rms = TMath::RMS(fitValues.begin(),fitValues.end());
        std::cout<<"mean = "<<mean<<", rms = "<<rms<<", rms/mean = "<<rms/mean<<std::endl;
    }

    hist_resol->GetXaxis()->SetTitle("Post-fit L_{#alpha}");
    hist_resol->GetYaxis()->SetTitle("Entries");
    hist_resol->GetYaxis()->SetTitleOffset(1.5);
    //hist_resol->GetYaxis()->SetRangeUser(0,0.1);
    
    hist_resol->SetLineColor(counter+1);
    hist_resol->SetLineStyle(counter+1);
    hist_resol->SetLineWidth(3);
    counter++;

    return hist_resol;
}

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

TH1D* plot_angular(std::string config, double sys_lo, double sys_hi)
{
    gStyle->SetOptStat(0);

    //std::string config = "fullring";
    //std::string config = "noring1";

    double costh_min = 0, costh_max = 1;
    int nbins = 40;
    TH1D* hist_resol = new TH1D("","",nbins,costh_min,costh_max);
    std::vector<std::vector<double>> fitValues(nbins);
    for (int i=0;i<100;i++)
    {
        TFile* f = new TFile(Form("%stoy%i.root",config.c_str(),i),"OPEN");
        TH1D* hist_alpha_result = (TH1D*)f->Get("hist_alpha_result");
        TH1D* hist_BnLPMT_angular_result = (TH1D*)f->Get("hist_BnLPMT_angular_result");
        if (hist_alpha_result)
        {
            for (int j=1;j<=nbins;j++)
            {
                double val = hist_BnLPMT_angular_result->GetBinContent(j)/hist_BnLPMT_angular_result->GetBinContent(nbins);
                fitValues[j-1].push_back(val);
            }
        }
        f->Close();
    }
    for (int j=0;j<nbins;j++)
    {
        if (fitValues[j].size()>0)
        {
            double mean = TMath::Mean(fitValues[j].begin(),fitValues[j].end());
            double rms = TMath::RMS(fitValues[j].begin(),fitValues[j].end());
            //std::cout<<"mean = "<<mean<<", rms = "<<rms<<", rms/mean = "<<rms/mean<<std::endl;
            hist_resol->SetBinContent(j+1,rms/mean);
            hist_resol->SetBinError(j+1,rms);
        }
    }

    hist_resol->GetXaxis()->SetTitle("cos#theta");
    hist_resol->GetYaxis()->SetTitle("A(cos#theta)");
    hist_resol->GetYaxis()->SetTitleOffset(1.5);
    //hist_resol->GetYaxis()->SetRangeUser(0,0.1);
    
    hist_resol->SetLineColor(counter+1);
    hist_resol->SetLineStyle(counter+1);
    hist_resol->SetLineWidth(3);
    counter++;

    return hist_resol;
}

TH2D* plot_angular_2d(std::string config, double sys_lo, double sys_hi)
{
    gStyle->SetOptStat(0);

    //std::string config = "fullring";
    //std::string config = "noring1";

    double costh_min = 0, costh_max = 1;
    int nbins = 40;
    const int ntoys = 100;
    TH1D* toy_result[ntoys];
    TH1D* toy_mean = new TH1D("","",nbins,0,nbins);
    TH2D* hist_cov = new TH2D("","",nbins,0,nbins,nbins,0,nbins);
    TH2D* hist_cor = new TH2D("","",nbins,0,nbins,nbins,0,nbins);
    std::vector<std::vector<double>> fitValues(nbins);
    for (int i=0;i<ntoys;i++)
    {
        toy_result[i] = new TH1D("","",nbins,0,nbins);
        TFile* f = new TFile(Form("%stoy%i.root",config.c_str(),i),"OPEN");
        TH1D* hist_alpha_result = (TH1D*)f->Get("hist_alpha_result");
        TH1D* hist_BnLPMT_angular_result = (TH1D*)f->Get("hist_BnLPMT_angular_result");
        if (hist_alpha_result)
        {
            toy_result[i]->SetBinContent(1,hist_alpha_result->GetBinContent(1));
            for (int j=1;j<=nbins-1;j++)
            {
                double val = hist_BnLPMT_angular_result->GetBinContent(j)/hist_BnLPMT_angular_result->GetBinContent(nbins);
                toy_result[i]->SetBinContent(j+1,val);
            }
        }
        toy_mean->Add(toy_result[i]);
        f->Close();
    }
    toy_mean->Scale(1./ntoys);

    TMatrixDSym toy_cov;
    TMatrixDSym toy_cor;

    toy_cov.ResizeTo(nbins, nbins);
    toy_cov.Zero();

    toy_cor.ResizeTo(nbins, nbins);
    toy_cor.Zero();

    for(int t=0;t<ntoys;t++)
    {
        for(int i = 0; i < nbins; ++i)
        {
            for(int j = 0; j < nbins; ++j)
            {
                const double x = toy_result[t]->GetBinContent(i + 1) - toy_mean->GetBinContent(i + 1);
                const double y = toy_result[t]->GetBinContent(j + 1) - toy_mean->GetBinContent(j + 1);
                toy_cov(i, j) += x * y / (1.0 * ntoys);
            }
        }
    }

    for(int i = 0; i < nbins; ++i)
    {
        for(int j = 0; j < nbins; ++j)
        {
            const double x = toy_cov(i, i);
            const double y = toy_cov(j, j);
            const double z = toy_cov(i, j);
            toy_cor(i, j) = z / (sqrt(x * y));

            if(std::isnan(toy_cor(i, j)))
                toy_cor(i, j) = 0.0;

            hist_cov->SetBinContent(i+1,j+1,toy_cov(i, j));
            hist_cor->SetBinContent(i+1,j+1,toy_cor(i, j));
        }
    }

    hist_cor->GetZaxis()->SetRangeUser(-1,1);
    return hist_cor;
}

void plot_eff_var_BnL()
{
    gStyle->SetOptStat(0);
    
    counter = 0;
    
    // TH1D* h_full = plot("fullring",0.017,0.009);
    // TH1D* h_half = plot("halfring",0.024,0.011);
    // TH1D* h_noout = plot("noring1",0.046,0.012);
    // TH1D* h_nomid = plot("noring2",0.024,0.010);
    TH1D* h_full = plot("BnL_eff_study/2pc",0.015,0.005);
    TH1D* h_half = plot("BnL_eff_study/10pc",0.02,0.0075);
    // TH1D* h_noout = plot("10pc/noring1",0.04,0.012);
    // TH1D* h_nomid = plot("10pc/noring2",0.016,0.006);
    // TH1D* h_full = plot("2pc/fullring",0.015,0.005);
    // TH1D* h_half = plot("2pc/halfring",0.02,0.0075);
    // TH1D* h_noout = plot("2pc/noring1",0.04,0.012);
    // TH1D* h_nomid = plot("2pc/noring2",0.016,0.006);

    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    h_full->Draw("hist same");
    h_half->Draw("hist same");
    // h_noout->Draw("hist same");
    // h_nomid->Draw("hist same");
    TLegend* legend = new TLegend(0.45,0.65,0.75,0.9);
    legend->AddEntry(h_full,"2 percent","l");
    legend->AddEntry(h_half,"10 percent","l");
    // legend->AddEntry(h_noout,"No outer","l");
    // legend->AddEntry(h_nomid,"No middle","l");
    legend->Draw("same");
    c1->SaveAs("Lalpha_eff_var_BnL.pdf");

    int orders[] = {2,3};
    double ranges[] = {0.0,0.6,1.01};
    pol_orders.assign(orders, orders+sizeof(orders)/sizeof(int));
    pol_range.assign(ranges, ranges+sizeof(ranges)/sizeof(double));

    counter = 0;
    TH1D* h_full_angular = plot_angular("BnL_eff_study/2pc",0.015,0.005);
    TH1D* h_half_angular = plot_angular("BnL_eff_study/10pc",0.02,0.0075);
    TCanvas* c2 = new TCanvas();
    c2->SetGridy();
    h_half_angular->Draw("E0 same");
    h_full_angular->Draw("E0 same");
    TLegend* legend2 = new TLegend(0.25,0.65,0.55,0.9);
    legend2->AddEntry(h_full,"2 percent","l");
    legend2->AddEntry(h_half,"10 percent","l");
    legend2->Draw("same");
    c2->SaveAs("angular_eff_var_BnL.pdf");

    TCanvas* c3 = new TCanvas();
    TH2D* h_half_angular_2d = plot_angular_2d("BnL_eff_study/10pc",0.02,0.0075);
    gStyle->SetPaintTextFormat("0.2f");
    h_half_angular_2d->Draw("colz text");
    c3->SaveAs("angular2d_eff_var_BnL.pdf");
}
