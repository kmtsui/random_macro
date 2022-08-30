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

void make_angular_response_prior()
{
    int old = 1;
    int pmttype = 0;
    std::string outname = old ?  ( pmttype == 0 ? "BnL_angular_prior_old.root" : "mPMT_angular_prior_old.root") : ( pmttype == 0 ? "BnL_angular_prior.root" : "mPMT_angular_prior.root");

    TFile* f = old ? new TFile("TN/diffusr4_400nm_nominal_combined_pol_fullcosth_40degcosths.root") : new TFile("TN/diffusr4_400nm_nominal_combined_pol_fullcosth_40degcosths_newMC.root");
    TVectorD* res_vector = (TVectorD*)f->Get("res_vector");
    TMatrixDSym* res_cov_matrix = (TMatrixDSym*)f->Get("res_cov_matrix");

    int startingIndex = pmttype == 0 ? 2 : 7;
    int orders0[] = {2,3}; int orders1[] = {2,3,3};
    double ranges0[] = {0.,0.6,1.0}; double ranges1[] = {0.,0.6,0.75,1.0};
    int npolpars = pmttype == 0 ? 5 : 7;

    if (pmttype==0)
    {
        pol_orders.assign(orders0, orders0+sizeof(orders0)/sizeof(int));
        pol_range.assign(ranges0, ranges0+sizeof(ranges0)/sizeof(double));
    }
    else if (pmttype==1)
    {
        pol_orders.assign(orders1, orders1+sizeof(orders1)/sizeof(int));
        pol_range.assign(ranges1, ranges1+sizeof(ranges1)/sizeof(double));
    }

    double par_pol[100];
    for (int i=0;i<npolpars;i++)
        par_pol[i] = (*res_vector)[i+startingIndex];
    
    TFile* fout = new TFile(outname.c_str(),"RECREATE");

    int nParameters = 40;
    double costh_min = 0.;
    double costh_max = 1.0;
    TH1D* h_fit_pol = new TH1D("","",nParameters,costh_min,costh_max);
    std::vector<double> pltpts;
    for (int i=1;i<=h_fit_pol->GetNbinsX();i++)
    {
        pltpts.push_back(h_fit_pol->GetBinCenter(i));
    }
    std::vector<double> valpts = CalcPol(par_pol,pltpts);
    for (int i=1;i<=h_fit_pol->GetNbinsX();i++)
    {
        h_fit_pol->SetBinContent(i,valpts[i-1]);
    }
    h_fit_pol->Write("prior");

    TMatrixDSym cov_matrix(nParameters);
    for(int r = 0; r < nParameters; ++r)
    {
        for(int c = 0; c < nParameters; ++c)
        {
            cov_matrix[r][c] = r==c? pow(h_fit_pol->GetBinContent(r+1)*0.01,2) : 0;
        }
    }
    cov_matrix.Write("cov_matrix_1percent");
    for(int r = 0; r < nParameters; ++r)
    {
        for(int c = 0; c < nParameters; ++c)
        {
            cov_matrix[r][c] = r==c? pow(h_fit_pol->GetBinContent(r+1)*0.05,2) : 0;
        }
    }
    cov_matrix.Write("cov_matrix_5percent");
    for(int r = 0; r < nParameters; ++r)
    {
        for(int c = 0; c < nParameters; ++c)
        {
            cov_matrix[r][c] = r==c? pow(h_fit_pol->GetBinContent(r+1)*0.10,2) : 0;
        }
    }
    cov_matrix.Write("cov_matrix_10percent");

    fout->Close();
    f->Close();
}