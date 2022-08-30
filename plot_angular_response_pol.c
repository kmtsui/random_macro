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

void plot_angular_response_pol_mPMT()
{
    int orders[] = {2,3,3};
    double ranges[] = {0.5,0.6,0.75,1.0};
    pol_orders.assign(orders, orders+sizeof(orders)/sizeof(int));
    pol_range.assign(ranges, ranges+sizeof(ranges)/sizeof(double));
    double par1[] = {0.20845,0.252388,-0.288586,0.315924,-0.183385};
    //double par2[] = {0.220456,0.264191,-0.326777,0.404403,-0.296335};
    double par2[] = {0.206815,0.234586,-0.0227563,0.104737,0.143297};
    TH1D* h_fit_pol1 = new TH1D("","",1000,0.5,1.0);
    TH1D* h_fit_pol2 = new TH1D("","",1000,0.5,1.0);
    std::vector<double> pltpts;
    for (int i=1;i<=h_fit_pol1->GetNbinsX();i++)
    {
        pltpts.push_back(h_fit_pol1->GetBinCenter(i));
    }
    std::vector<double> valpts1 = CalcPol(par1,pltpts);
    std::vector<double> valpts2 = CalcPol(par2,pltpts);
    for (int i=1;i<=h_fit_pol1->GetNbinsX();i++)
    {
        h_fit_pol1->SetBinContent(i,valpts1[i-1]);
        h_fit_pol2->SetBinContent(i,valpts2[i-1]);
    }
    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas();
    h_fit_pol1->SetLineColor(kBlue);
    h_fit_pol2->SetLineColor(kRed);
    h_fit_pol1->Draw("same");
    h_fit_pol2->Draw("same");
    c1->SaveAs("test.pdf");
    TH1D* h_fit_pol_ratio = (TH1D*)h_fit_pol1->Clone();
    h_fit_pol1->Divide(h_fit_pol2);
    h_fit_pol1->Draw("hist");
    c1->SaveAs("test1.pdf");
}

void plot_angular_response_pol_BnL()
{
    int orders[] = {2,3};
    double ranges[] = {0.5,0.6,1.0};
    pol_orders.assign(orders, orders+sizeof(orders)/sizeof(int));
    pol_range.assign(ranges, ranges+sizeof(ranges)/sizeof(double));

    const int nCurves = 5;
    const int npars = 5;
    double pars[nCurves][npars] = {
        {0.213834,0.218174,0.121592,-0.0341528,0.358366},
        {0.21563,0.205651,0.182075,0.0307789,0.232192},
        {0.215885,0.199067,0.210368,0.0114371,0.27878},
        {0.216972,0.208817,0.163935,0.035006,0.207446},
        {0.213328,0.203111,0.203218,-0.0413365,0.393139},
    };
    //double truth_alpha[nCurves] = {10648,11001,10159,12304,8860};
    double truth_alpha[nCurves] = {6351,6430,6236,7510,5157};
    int wavelength = 350;
    std::string water[] = {"nominal","abwplus20pc","abwminus20pc","rayplus20pc","rayminus20pc"};

    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas();
    TLegend* legend = new TLegend(0.1,0.15,0.6,0.45);
    TH1D* h_fit_pol[nCurves];
    TH1D* h_fit_ratio[nCurves];
    for (int i=0;i<nCurves;i++)
    {
        TFile* f = new TFile(Form("diffusr4_%inm_%s_BnL_indirect.root",wavelength,water[i].c_str()));
        double alpha = ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
        double error = ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);
        TVectorD* res_vector = (TVectorD*)f->Get("res_vector");
        int startingIndex = 2;
        for (int j=0;j<npars;j++)
        {
            pars[i][j]= (*res_vector)[j+startingIndex];
        }
        h_fit_pol[i] = new TH1D("","",1000,0.5,1.0);
        std::vector<double> pltpts;
        for (int j=1;j<=h_fit_pol[i]->GetNbinsX();j++)
        {
            pltpts.push_back(h_fit_pol[i]->GetBinCenter(j));
        }
        std::vector<double> valpts = CalcPol(pars[i],pltpts);
        for (int j=1; j<=h_fit_pol[i]->GetNbinsX();j++)
        {
            h_fit_pol[i]->SetBinContent(j,valpts[j-1]);
        }
        h_fit_pol[i]->SetLineColor(i+1);
        h_fit_pol[i]->SetLineStyle(i+1);
        h_fit_pol[i]->SetLineWidth(3);
        h_fit_pol[i]->GetXaxis()->SetTitle("cos#theta");
        //h_fit_pol[i]->Draw("same");
        h_fit_ratio[i] = (TH1D*)h_fit_pol[i]->Clone();
        h_fit_ratio[i]->Divide(h_fit_pol[0]);
        h_fit_ratio[i]->GetYaxis()->SetRangeUser(0.93,1.05);
        h_fit_ratio[i]->GetYaxis()->SetTitle("A/A_{nominal}");
        h_fit_ratio[i]->Draw("hist same");
        legend->AddEntry(h_fit_ratio[i],Form("%s: %3.0f ; %3.0f#pm%2.0f",water[i].c_str(),truth_alpha[i],alpha,error),"l");
    }
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();
    c1->SaveAs("test.pdf");
}

void plot_angular_response_pol()
{
    plot_angular_response_pol_BnL();
}