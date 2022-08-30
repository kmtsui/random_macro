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

double BarlowLLH(double mc, double w2, double data)
{
    // Solving for the quadratic equation,
    // beta^2 + (mu * sigma^2 - 1)beta - data * sigma^2) = 0
    // where sigma^2 is the relative variance.
    double rel_var = w2 / (mc * mc);
    double b       = (mc * rel_var) - 1;
    double c       = 4 * data * rel_var;

    double beta   = (-b + std::sqrt(b * b + c)) / 2.0;
    double mc_hat = mc * beta;

    // Calculate the following LLH:
    //-2lnL = 2 * beta*mc - data + data * ln(data / (beta*mc)) + (beta-1)^2 / sigma^2
    // where sigma^2 is the same as above.
    double chi2 = 0.0;
    //if(data <= 0.0)
    //{
    //    chi2 = 2 * mc_hat;
    //    chi2 += (beta - 1) * (beta - 1) / rel_var;
    //}
    if(mc_hat > 0.0)
    {
        chi2 = 2 * (mc_hat - data);
        if(data > 0.0)
            chi2 += 2 * data * std::log(data / mc_hat);

        if (rel_var > 0.0) chi2 += (beta - 1) * (beta - 1) / rel_var;
    }

    return (chi2 >= 0.0) ? chi2 : 0.0;
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

TH1D* get_plot(std::string filename, double truth_alpha, int pmt_type = 0, int startingIndex = 2, bool led = false, bool attez = false){
    gStyle->SetOptStat(0);

    int orders_BnL[] = {3,3};
    double ranges_BnL[] = {0.0,0.7,1.0};
    if (led) ranges_BnL[0]=0.2;
    int orders_mPMT[] = {3,3,3,4};
    double ranges_mPMT[] = {0.,0.3,0.6,0.75,1.0};
    pol_orders.clear();
    pol_range.clear();
    if (pmt_type==0)
    {
        pol_orders.assign(orders_BnL, orders_BnL+sizeof(orders_BnL)/sizeof(int));
        pol_range.assign(ranges_BnL, ranges_BnL+sizeof(ranges_BnL)/sizeof(double));
    }
    else
    {
        pol_orders.assign(orders_mPMT, orders_mPMT+sizeof(orders_mPMT)/sizeof(int));
        pol_range.assign(ranges_mPMT, ranges_mPMT+sizeof(ranges_mPMT)/sizeof(double));
    }

    int npars = pol_orders[0]+1;
    for (int i=1;i<pol_orders.size();i++)
        npars+=pol_orders[i]-1;
    double pars[100];
    std::vector<std::string> par_name;
    for (int i=0;i<pol_orders.size();i++)
    {
        int j=0;
        if (i!=0) j=2;
        for (;j<pol_orders[i]+1;j++)
        {
            par_name.push_back(Form("pol%i_p%i",pol_orders[i],j));
        }
    }

    TFile* f = new TFile(filename.c_str());
    TVectorD* res_vector = (TVectorD*)f->Get("res_vector");
    for (int j=0;j<npars;j++)
    {
        pars[j]= (*res_vector)[j+startingIndex];
    }
    TH1D* h_fit_pol = new TH1D("","",1000,pol_range.front(),pol_range.back());
    std::vector<double> pltpts;
    for (int j=1;j<=h_fit_pol->GetNbinsX();j++)
    {
        pltpts.push_back(h_fit_pol->GetBinCenter(j));
    }
    std::vector<double> valpts = CalcPol(pars,pltpts);
    for (int j=1; j<=h_fit_pol->GetNbinsX();j++)
    {
        h_fit_pol->SetBinContent(j,valpts[j-1]);
    }
    h_fit_pol->SetLineWidth(3);
    h_fit_pol->GetXaxis()->SetTitle("cos#theta");
    h_fit_pol->GetYaxis()->SetTitle("A(cos#theta)");

    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX());
    int ndof = ((TTree*)f->Get("PMTTree"))->GetEntries()-npars-1;
    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    h_fit_pol->Draw("");
    double fitted_alpha = attez ? ((TH1D*)f->Get("hist_alphaZ_result"))->GetBinContent(1) : ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = attez ? ((TH1D*)f->Get("hist_alphaZ_error_final"))->GetBinContent(1) : ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);
    TLatex latex;
    latex.SetTextSize(0.045);
    double miny = h_fit_pol->GetMinimum()*1.1;
    double maxy = h_fit_pol->GetMaximum();
    double diff = (maxy-miny)/10;
    latex.DrawLatex(0.1,miny+diff*7,Form("#chi^{2}/dof = %4.0f/%i",chi2,ndof));
    latex.DrawLatex(0.1,miny+diff*6,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha));
    latex.DrawLatex(0.1,miny+diff*5,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha));
    for (int i=1;i<pol_range.size();i++)
    {
        TLine* l = new TLine(pol_range[i],h_fit_pol->GetMinimum(),pol_range[i],maxy);
        l->SetLineStyle(2);
        l->Draw("same");
    }
    for (int i=0;i<pol_orders.size();i++)
    {
        latex.DrawLatex(pol_range[i]+0.02,maxy,Form("pol%i",pol_orders[i]));
    }
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    c1->SetGridy(false);
    TMatrixDSym* res_cor_matrix =  (TMatrixDSym*)f -> Get("res_cor_matrix");
    res_cor_matrix->Draw("colz text");
    TH2D *h_postfit_cor = (TH2D*) c1->GetPrimitive("TMatrixDBase");
    h_postfit_cor->GetZaxis()->SetRangeUser(-1,1);
    h_postfit_cor->GetXaxis()->SetBinLabel(2,"L_{#alpha}");
    //h_postfit_cor->GetXaxis()->SetBinLabel(3,"A_{mPMT}");
    h_postfit_cor->GetXaxis()->SetRangeUser(1,1+1+npars);
    h_postfit_cor->GetYaxis()->SetBinLabel(2,"L_{#alpha}");
    for (int i=3;i<3+npars;i++)
    {
        h_postfit_cor->GetXaxis()->SetBinLabel(i,par_name[i-3].c_str());
        h_postfit_cor->GetYaxis()->SetBinLabel(i,par_name[i-3].c_str());
    }
    //h_postfit_cor->GetYaxis()->SetBinLabel(3,"A_{mPMT}");
    h_postfit_cor->GetYaxis()->SetRangeUser(1,1+1+npars);
    figname += "_cor";
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    return h_fit_pol;
}

void get_plot_3angular(std::string filename, double truth_alpha)
{
    TH1D* h_mPMT = get_plot(filename,truth_alpha,1);

    TFile* f = new TFile(filename.c_str());

    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX());
    int ndof = ((TTree*)f->Get("PMTTree"))->GetEntries()-1-7-12;
    double chi2perdof = chi2/ndof;
    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    h_mPMT->GetXaxis()->SetTitle("cos#theta");
    h_mPMT->GetYaxis()->SetTitle("A_{#theta}");
    h_mPMT->Draw("");
    double fitted_alpha = ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);
    TLatex latex;
    latex.SetTextSize(0.045);
    double miny = h_mPMT->GetMinimum()*1.2;
    latex.DrawLatex(0.55,miny+0.06,Form("#chi^{2}/dof = %4.0f/%i",chi2,ndof));
    latex.DrawLatex(0.55,miny+0.04,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha));
    latex.DrawLatex(0.55,miny+0.02,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha));
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    std::string anotherfilename = filename;
    anotherfilename.erase(anotherfilename.find("_3angular"),9);
    TH1D* another_angular = get_plot(anotherfilename,truth_alpha,1);
    c1->cd();
    another_angular->SetLineColor(kRed);
    another_angular->SetLineStyle(2);
    another_angular->Draw("same");
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    TH1D* hist_phim_result = (TH1D*)f->Get("hist_phim_result");
    TH1D* hist_phim_error_final = (TH1D*)f->Get("hist_phim_error_final");
    TH1D* hist_mPMT_id_result = (TH1D*)f->Get("hist_mPMT_id_result");
    TH1D* hist_mPMT_id_error_final = (TH1D*)f->Get("hist_mPMT_id_error_final");
    double phim_bins[] = {0,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0,TMath::Pi()};
    TH1D* hist_phim = new TH1D("","",11,phim_bins);
    for (int i=1;i<=11;i++)
    {
        hist_phim->SetBinContent(i,hist_phim_result->GetBinContent(i));
        hist_phim->SetBinError(i,hist_phim_error_final->GetBinContent(i));
    }
    hist_phim->GetXaxis()->SetTitle("#phi_{m}");
    hist_phim->GetYaxis()->SetTitle("A_{#phi_{m}}");
    hist_phim->GetYaxis()->SetRangeUser(0.99,1.01);
    hist_phim->Draw("E0");
    latex.DrawLatex(1.5,0.996,Form("A_{Ring}(1) = 1 (fixed)"));
    latex.DrawLatex(1.5,0.995,Form("A_{Ring}(2) = %4.4f +/- %4.4f",hist_mPMT_id_result->GetBinContent(2),hist_mPMT_id_error_final->GetBinContent(2)));
    latex.DrawLatex(1.5,0.994,Form("A_{Ring}(3) = %4.4f +/- %4.4f",hist_mPMT_id_result->GetBinContent(3),hist_mPMT_id_error_final->GetBinContent(3)));
    c1->SaveAs(Form("%s_phim.pdf",figname.c_str()));
}

void get_plot_combined(std::string filename, double truth_alpha)
{
    TH1D* h_BnL = get_plot(filename,truth_alpha,0,2);
    TH1D* h_mPMT = get_plot(filename,truth_alpha,1,8);

    TFile* f = new TFile(filename.c_str());
    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX());
    int npars_BnL = 6;
    int npars_mPMT = 11;
    int npars = npars_BnL+npars_mPMT;
    int ndof = ((TTree*)f->Get("PMTTree"))->GetEntries()-1;
    double fitted_alpha = ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);

    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    h_mPMT->SetLineColor(kBlue);
    h_mPMT->SetMaximum(h_BnL->GetMaximum()*1.1);
    h_mPMT->Draw();
    h_BnL->SetLineColor(kRed);
    h_BnL->Draw("same");
    TLegend* legend = new TLegend(0.15,0.55,0.5,0.85);
    legend->AddEntry(h_BnL,"A_{BnL}","lep");
    legend->AddEntry(h_mPMT,"A_{mPMT}","lep");
    legend->AddEntry((TObject*)0,Form("#chi^{2}/dof = %4.0f/%i",chi2,ndof),"");
    legend->AddEntry((TObject*)0,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha),"");
    legend->AddEntry((TObject*)0,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha),"");
    legend->Draw();
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    gStyle->SetPaintTextFormat("4.3f");
    c1->SetGridy(false);
    TMatrixDSym* res_cor_matrix =  (TMatrixDSym*)f -> Get("res_cor_matrix");
    res_cor_matrix->Draw("colz text");
    TH2D *h_postfit_cor = (TH2D*) c1->GetPrimitive("TMatrixDBase");
    h_postfit_cor->GetZaxis()->SetRangeUser(-1,1);
    h_postfit_cor->GetXaxis()->SetBinLabel(2,"L_{#alpha}");
    h_postfit_cor->GetXaxis()->SetBinLabel(3,"BnL_pol");
    h_postfit_cor->GetXaxis()->SetBinLabel(3+npars_BnL,"mPMT_pol");
    h_postfit_cor->GetXaxis()->SetRangeUser(1,1+1+npars);
    h_postfit_cor->GetYaxis()->SetBinLabel(2,"L_{#alpha}");
    h_postfit_cor->GetYaxis()->SetBinLabel(3,"BnL_pol");
    h_postfit_cor->GetYaxis()->SetBinLabel(3+npars_BnL,"mPMT_pol");
    h_postfit_cor->GetYaxis()->SetRangeUser(1,1+1+npars);
    figname += "_cor";
    c1->SaveAs(Form("%s.pdf",figname.c_str()));
}

void get_plot_led(std::string filename, double truth_alpha, int pmt_type = 0)
{
    TH1D* h_pol_led = get_plot(filename,truth_alpha,pmt_type);
    std::string filename_nom = filename;
    filename_nom.erase(filename_nom.find("_led"));
    //filename_nom.substr(0,filename_nom.find("_led"));
    filename_nom += ".root";
    TH1D* h_pol_nom = get_plot(filename_nom,truth_alpha,pmt_type);

    TFile* f = new TFile(filename.c_str());
    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX()-1);
    int npars = pmt_type==0 ? 5:7;
    int ndof = ((TTree*)f->Get("PMTTree"))->GetEntries()-npars-1;
    double fitted_alpha = ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);

    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    h_pol_led->SetLineColor(kBlue);
    h_pol_led->Draw();
    h_pol_nom->SetLineColor(kRed);
    h_pol_nom->SetLineStyle(2);
    h_pol_nom->Draw("same");
    TLegend* legend = new TLegend(0.5,0.15,0.9,0.45);
    legend->AddEntry(h_pol_led,"Fake data","lep");
    legend->AddEntry(h_pol_nom,"Nominal MC","lep");
    legend->AddEntry((TObject*)0,Form("#chi^{2}/dof = %4.0f/%i",chi2,ndof),"");
    legend->AddEntry((TObject*)0,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha),"");
    legend->AddEntry((TObject*)0,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha),"");
    legend->Draw();
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    TH1D* hist_LED_cosths_result = (TH1D*)f->Get("hist_LED_cosths_result");
    TH1D* hist_LED_cosths_prior = (TH1D*)f->Get("hist_LED_cosths_prior");
    TH1D* hist_LED_cosths_error_final = (TH1D*)f->Get("hist_LED_cosths_error_final");
    TH1D* hist_LED_cosths_error_prior = (TH1D*)f->Get("hist_LED_cosths_error_prior");
    TH1D* hist_LED_cosths_prefit = new TH1D("","",480,-40,0);
    TH1D* hist_LED_cosths_postfit = new TH1D("","",480,-40,0);
    for (int i=1;i<=hist_LED_cosths_prefit->GetNbinsX();i++)
    {
        hist_LED_cosths_prefit->SetBinContent(i,hist_LED_cosths_prior->GetBinContent(i));
        hist_LED_cosths_prefit->SetBinError(i,hist_LED_cosths_error_prior->GetBinContent(i));
        hist_LED_cosths_postfit->SetBinContent(i,hist_LED_cosths_result->GetBinContent(i));
        hist_LED_cosths_postfit->SetBinError(i,hist_LED_cosths_error_final->GetBinContent(i));
        if (hist_LED_cosths_error_final->GetBinContent(i)>1.e-6) ndof--;
    }
    hist_LED_cosths_prefit->GetXaxis()->SetTitle("-#theta_{s} (deg)");
    hist_LED_cosths_prefit->GetYaxis()->SetTitle("I(#theta_{s},#phi_{s})");
    hist_LED_cosths_prefit->SetLineColor(kRed);
    hist_LED_cosths_prefit->Draw("E0");
    hist_LED_cosths_postfit->SetLineColor(kBlue);
    hist_LED_cosths_postfit->Draw("E0 same");
    legend = new TLegend(0.5,0.2,0.8,0.4);
    legend->AddEntry(hist_LED_cosths_prefit,"Pre-fit","lep");
    legend->AddEntry(hist_LED_cosths_postfit,"Post-fit","lep");
    TH1D* chi2_syst_periter = (TH1D*)f->Get("chi2_syst_periter");
    double chi2syst = chi2_syst_periter->GetBinContent(chi2_syst_periter->GetNbinsX()-1);
    legend->AddEntry((TObject*)0,Form("#chi^{2}_{syst}=%2.0f",chi2syst),"");
    legend->SetBorderSize(0);
    legend->Draw();
    c1->SaveAs(Form("%s_intensity.pdf",figname.c_str()));

}

void get_plot_compare_nominal(std::string filename, std::string filename_nom, double truth_alpha, int pmt_type = 0, int npars=24)
{
    TH1D* h_pol_led = pmt_type ==0 ? get_plot(filename,truth_alpha,pmt_type,2,true) : get_plot(filename,truth_alpha,pmt_type);
    TH1D* h_pol_nom = get_plot(filename_nom,truth_alpha,pmt_type);

    TFile* f = new TFile(filename.c_str());
    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX());
    TH1D* evhist_sam0_pred_final = (TH1D*)f->Get("evhist_sam0_pred_final");
    int ndof = 0;
    for (int i=1;i<=evhist_sam0_pred_final->GetNbinsX();i++)
        if (evhist_sam0_pred_final->GetBinContent(i)>0)
            ndof++;
    ndof -= npars;
    //int ndof = ((TTree*)f->Get("PMTTree"))->GetEntries()-npars;
    double fitted_alpha = ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);

    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    h_pol_led->SetLineColor(kBlue);
    h_pol_nom->SetLineColor(kRed);
    h_pol_nom->SetLineStyle(2);
    h_pol_nom->GetXaxis()->SetRangeUser(h_pol_led->GetXaxis()->GetXmin(),h_pol_led->GetXaxis()->GetXmax());
    h_pol_nom->GetYaxis()->SetRangeUser(h_pol_nom->GetMinimum()*0.8,h_pol_nom->GetMaximum()*1.2);
    h_pol_nom->Draw("same");
    h_pol_led->Draw("same");
    TLegend* legend = new TLegend(0.15,0.55,0.55,0.85);
    legend->AddEntry(h_pol_led,"Fake data","lep");
    legend->AddEntry(h_pol_nom,"Nominal MC","lep");
    legend->AddEntry((TObject*)0,Form("#chi^{2}/dof = %4.0f/%i",chi2,ndof),"");
    legend->AddEntry((TObject*)0,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha),"");
    legend->AddEntry((TObject*)0,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha),"");
    legend->Draw();
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    c1->SaveAs(Form("%s.pdf",figname.c_str()));
}

void get_plot_linz_compare_nominal(std::string filename, std::string filename_nom, double truth_alpha, int pmt_type = 0, int npars=24)
{
    TH1D* h_pol_led = pmt_type ==0 ? get_plot(filename,truth_alpha,pmt_type,2,true) : get_plot(filename,truth_alpha,pmt_type);
    TH1D* h_pol_nom = get_plot(filename_nom,truth_alpha,pmt_type);

    TFile* f = new TFile(filename.c_str());
    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX());
    double chi2_stat = ((TH1D*)f->Get("chi2_stat_periter"))->GetBinContent(chi2_total_periter->GetNbinsX());
    double chi2_syst = ((TH1D*)f->Get("chi2_syst_periter"))->GetBinContent(chi2_total_periter->GetNbinsX());
    TH1D* chi2_syst_periter = (TH1D*)f->Get("chi2_syst_periter");
    TH1D* evhist_sam0_pred_final = (TH1D*)f->Get("evhist_sam0_pred_final");
    int ndof = 0;
    for (int i=1;i<=evhist_sam0_pred_final->GetNbinsX();i++)
        if (evhist_sam0_pred_final->GetBinContent(i)>0)
            ndof++;
    for (int j=1;;j++)
    {
        TH1D* evhist = (TH1D*)f->Get(Form("evhist_sam%i_pred_final",j));
        if (!evhist) break;
        for (int i=1;i<=evhist->GetNbinsX();i++)
            if (evhist->GetBinContent(i)>0)
                ndof++;
    }
    double fitted_alpha = ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);

    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    h_pol_led->SetLineColor(kBlue);
    h_pol_nom->SetLineColor(kRed);
    h_pol_nom->SetLineStyle(2);
    h_pol_nom->GetXaxis()->SetRangeUser(h_pol_led->GetXaxis()->GetXmin(),h_pol_led->GetXaxis()->GetXmax());
    h_pol_nom->GetYaxis()->SetRangeUser(h_pol_nom->GetMinimum()*0.8,h_pol_nom->GetMaximum()*1.2);
    h_pol_nom->Draw("same");
    h_pol_led->Draw("same");
    TLegend* legend = new TLegend(0.15,0.55,0.55,0.85);
    legend->AddEntry(h_pol_led,"Fake data","lep");
    legend->AddEntry(h_pol_nom,"Nominal MC","lep");
    legend->AddEntry((TObject*)0,Form("#chi^{2}_{stat}/bin = %4.0f/%i",chi2,ndof),"");
    legend->AddEntry((TObject*)0,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha),"");
    legend->AddEntry((TObject*)0,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha),"");
    legend->Draw();
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    TH1D* hist_zpos_eff_result = (TH1D*)f->Get("hist_zpos_eff_result");
    TH1D* hist_zpos_eff_prior = (TH1D*)f->Get("hist_zpos_eff_prior");
    TH1D* hist_zpos_eff_error_final = (TH1D*)f->Get("hist_zpos_eff_error_final");
    TH1D* hist_zpos_eff_error_prior = (TH1D*)f->Get("hist_zpos_eff_error_prior");
    TH1D* h_linz = new TH1D("","",hist_zpos_eff_result->GetNbinsX(),-3300,3300);
    for (int i=1;i<=h_linz->GetNbinsX();i++)
    {
        h_linz->SetBinContent(i,hist_zpos_eff_result->GetBinContent(i));
        h_linz->SetBinError(i,hist_zpos_eff_error_final->GetBinContent(i));
    }
    h_linz->GetXaxis()->SetTitle("z (cm)");
    h_linz->GetYaxis()->SetTitle("#epsilon(z)");
    h_linz->GetYaxis()->SetRangeUser(hist_zpos_eff_prior->GetBinContent(1)-2*hist_zpos_eff_error_prior->GetBinContent(1),
                                     hist_zpos_eff_prior->GetBinContent(h_linz->GetNbinsX())+2*hist_zpos_eff_error_prior->GetBinContent(h_linz->GetNbinsX()));
    h_linz->Draw();
    for (int i=1;i<=h_linz->GetNbinsX();i++)
    {
        if (hist_zpos_eff_prior->GetBinContent(i)>0)
        {
            TBox* box = new TBox(h_linz->GetBinLowEdge(i),hist_zpos_eff_prior->GetBinContent(i)-hist_zpos_eff_error_prior->GetBinContent(i),
                                 h_linz->GetBinLowEdge(i+1),hist_zpos_eff_prior->GetBinContent(i)+hist_zpos_eff_error_prior->GetBinContent(i));
            box->SetFillColor(kRed);
            box->Draw("same");
        }
    }
    h_linz->Draw("same");
    hist_zpos_eff_prior->SetLineColor(kRed);
    legend = new TLegend(0.2,0.6,0.4,0.8);
    legend->AddEntry(hist_zpos_eff_prior,"Pre-fit","lep");
    legend->AddEntry(h_linz,"Post-fit","lep");
    legend->AddEntry((TObject*)0,Form("#chi^{2}_{syst}=%2.1f",chi2_syst),"");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();
    c1->SaveAs(Form("%s_effz.pdf",figname.c_str()));
}

void get_plot_attenz_compare_nominal(std::string filename, std::string filename_nom, double truth_alpha, double truth_beta, int pmt_type = 0, int npars=24)
{
    TH1D* h_pol_led = pmt_type ==0 ? get_plot(filename,truth_alpha,pmt_type,3,true,true) : get_plot(filename,truth_alpha,pmt_type,3,false,true);
    TH1D* h_pol_nom = get_plot(filename_nom,truth_alpha,pmt_type);

    TFile* f = new TFile(filename.c_str());
    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX());
    TH1D* evhist_sam0_pred_final = (TH1D*)f->Get("evhist_sam0_pred_final");
    int ndof = 0;
    for (int i=1;i<=evhist_sam0_pred_final->GetNbinsX();i++)
        if (evhist_sam0_pred_final->GetBinContent(i)>0)
            ndof++;
    for (int j=1;;j++)
    {
        TH1D* evhist = (TH1D*)f->Get(Form("evhist_sam%i_pred_final",j));
        if (!evhist) break;
        for (int i=1;i<=evhist->GetNbinsX();i++)
            if (evhist->GetBinContent(i)>0)
                ndof++;
    }
    ndof -= npars;
    //int ndof = ((TTree*)f->Get("PMTTree"))->GetEntries()-npars;
    double fitted_alpha = ((TH1D*)f->Get("hist_alphaZ_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)f->Get("hist_alphaZ_error_final"))->GetBinContent(1);
    double fitted_beta = ((TH1D*)f->Get("hist_alphaZ_result"))->GetBinContent(2);
    double error_beta = ((TH1D*)f->Get("hist_alphaZ_error_final"))->GetBinContent(2);

    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    h_pol_led->SetLineColor(kBlue);
    h_pol_nom->SetLineColor(kRed);
    h_pol_nom->SetLineStyle(2);
    h_pol_nom->SetLineWidth(5);
    h_pol_nom->GetXaxis()->SetRangeUser(h_pol_led->GetXaxis()->GetXmin(),h_pol_led->GetXaxis()->GetXmax());
    h_pol_nom->GetYaxis()->SetRangeUser(h_pol_nom->GetMinimum()*0.8,h_pol_nom->GetMaximum()*1.2);
    h_pol_nom->Draw("same");
    h_pol_led->Draw("same");
    TLegend* legend = new TLegend(0.15,0.55,0.55,0.85);
    legend->AddEntry(h_pol_led,"Fake data","lep");
    legend->AddEntry(h_pol_nom,"Nominal MC","lep");
    legend->AddEntry((TObject*)0,Form("#chi^{2}/dof = %4.0f/%i",chi2,ndof),"");
    legend->AddEntry((TObject*)0,Form("Truth L_{#alpha}(0) = %4.0f cm",truth_alpha),"");
    legend->AddEntry((TObject*)0,Form("Fitted L_{#alpha}(0) = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha),"");
    legend->AddEntry((TObject*)0,Form("Truth #beta = %4.2f",truth_beta),"");
    legend->AddEntry((TObject*)0,Form("Fitted #beta = %3.3f +/- %3.3f",fitted_beta,error_beta),"");
    legend->Draw();
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    c1->SaveAs(Form("%s.pdf",figname.c_str()));
}

void get_plot_compare_nominal_joint(std::string filename, std::string filename_nom_BnL, std::string filename_nom_mPMT, double truth_alpha)
{
    TH1D* h_pol_BnL  = get_plot(filename,truth_alpha,0,2,true);
    TH1D* h_pol_mPMT = get_plot(filename,truth_alpha,1,8);
    TH1D* h_pol_BnL_nom  = get_plot(filename_nom_BnL,truth_alpha,0);
    TH1D* h_pol_mPMT_nom = get_plot(filename_nom_mPMT,truth_alpha,1);

    TFile* f = new TFile(filename.c_str());
    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX());
    TH1D* evhist_sam0_pred_final = (TH1D*)f->Get("evhist_sam0_pred_final");
    TH1D* evhist_sam0_data = (TH1D*)f->Get("evhist_sam0_data");
    TH1D* evhist_err2_sam0_pred_final = (TH1D*)f->Get("evhist_err2_sam0_pred_final");
    int ndof_BnL = 0;
    double chi2_BnL = 0;
    for (int i=1;i<=evhist_sam0_pred_final->GetNbinsX();i++)
    {
        if (evhist_sam0_pred_final->GetBinContent(i)>0)
        {
            ndof_BnL++;
            chi2_BnL += BarlowLLH(evhist_sam0_pred_final->GetBinContent(i),evhist_err2_sam0_pred_final->GetBinContent(i),evhist_sam0_data->GetBinContent(i));
        }
    }
    TH1D* evhist_sam1_pred_final = (TH1D*)f->Get("evhist_sam1_pred_final");
    TH1D* evhist_sam1_data = (TH1D*)f->Get("evhist_sam1_data");
    TH1D* evhist_err2_sam1_pred_final = (TH1D*)f->Get("evhist_err2_sam1_pred_final");
    int ndof_mPMT = 0;
    double chi2_mPMT = 0;
    for (int i=1;i<=evhist_sam1_pred_final->GetNbinsX();i++)
    {
        if (evhist_sam1_pred_final->GetBinContent(i)>0)
        {
            ndof_mPMT++;
            chi2_mPMT += BarlowLLH(evhist_sam1_pred_final->GetBinContent(i),evhist_err2_sam1_pred_final->GetBinContent(i),evhist_sam1_data->GetBinContent(i));
        }
    }
    std::cout<<"ndof_BnL = "<<ndof_BnL<<", ndof_mPMT = "<<ndof_mPMT<<std::endl;
    //ndof -= npars;
    //int ndof = ((TTree*)f->Get("PMTTree"))->GetEntries()-npars;
    double fitted_alpha = ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);

    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    h_pol_BnL->SetLineColor(kBlue);
    h_pol_BnL_nom->SetLineColor(kRed);
    h_pol_BnL_nom->SetLineStyle(2);
    h_pol_BnL_nom->GetXaxis()->SetRangeUser(h_pol_BnL->GetXaxis()->GetXmin(),h_pol_BnL->GetXaxis()->GetXmax());
    h_pol_BnL_nom->Draw("same");
    h_pol_BnL->Draw("same");
    TLegend* legend = new TLegend(0.15,0.55,0.55,0.85);
    legend->AddEntry(h_pol_BnL,"Fake data","lep");
    legend->AddEntry(h_pol_BnL_nom,"Nominal MC","lep");
    legend->AddEntry((TObject*)0,Form("#chi^{2}_{stat}/bin = %4.0f/%i",chi2_BnL,ndof_BnL),"");
    legend->AddEntry((TObject*)0,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha),"");
    legend->AddEntry((TObject*)0,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha),"");
    legend->Draw();
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    c1->SaveAs(Form("%s_BnL.pdf",figname.c_str()));

    TCanvas* c2 = new TCanvas();
    c2->SetGridy();
    h_pol_mPMT->SetLineColor(kBlue);
    h_pol_mPMT_nom->SetLineColor(kRed);
    h_pol_mPMT_nom->SetLineStyle(2);
    h_pol_mPMT_nom->GetXaxis()->SetRangeUser(h_pol_mPMT->GetXaxis()->GetXmin(),h_pol_mPMT->GetXaxis()->GetXmax());
    h_pol_mPMT_nom->Draw("");
    h_pol_mPMT->Draw("same");
    legend = new TLegend(0.15,0.55,0.55,0.85);
    legend->AddEntry(h_pol_mPMT,"Fake data","lep");
    legend->AddEntry(h_pol_mPMT_nom,"Nominal MC","lep");
    legend->AddEntry((TObject*)0,Form("#chi^{2}_{stat}/bin = %4.0f/%i",chi2_mPMT,ndof_mPMT),"");
    legend->AddEntry((TObject*)0,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha),"");
    legend->AddEntry((TObject*)0,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha),"");
    legend->Draw();
    c2->SaveAs(Form("%s_mPMT.pdf",figname.c_str()));
}

void get_plot_chi2(std::string filename)
{
    gStyle->SetOptStat(0);

    TFile* f = new TFile(filename.c_str());
    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX());
    TH1D* evhist_sam0_pred_final = (TH1D*)f->Get("evhist_sam0_pred_final");
    TH1D* evhist_sam0_data = (TH1D*)f->Get("evhist_sam0_data");
    TH1D* evhist_err2_sam0_pred_final = (TH1D*)f->Get("evhist_err2_sam0_pred_final");
    TH1D* hChi2 = new TH1D("","",evhist_sam0_pred_final->GetNbinsX(),0,evhist_sam0_pred_final->GetNbinsX());
    TH2D* hChi2_R = new TH2D("","",100,2600,8100,100,0,20);
    TH2D* hChi2_costh = new TH2D("","",100,0,1,100,0,20);
    TH2D* hChi2_cosths = new TH2D("","",100,0.767,1,100,0,20);
    TH2D* hChi2_zpos = new TH2D("","",100,-3300,3300,100,0,20);
    TH1D* hRatio = new TH1D("","",evhist_sam0_pred_final->GetNbinsX(),0,evhist_sam0_pred_final->GetNbinsX());
    TTree* PMTTree = (TTree*)f->Get("PMTTree");
    double R, costh, cosths, phis, zpos;
    PMTTree->SetBranchAddress("R",&R);
    PMTTree->SetBranchAddress("costh",&costh);
    PMTTree->SetBranchAddress("cosths",&cosths);
    PMTTree->SetBranchAddress("phis",&phis);
    PMTTree->SetBranchAddress("zpos",&zpos);
    for (int i=1;i<=evhist_sam0_pred_final->GetNbinsX();i++)
    {
        PMTTree->GetEntry(i-1);
        if (evhist_sam0_pred_final->GetBinContent(i)>0)
        {
            double val = BarlowLLH(evhist_sam0_pred_final->GetBinContent(i),evhist_err2_sam0_pred_final->GetBinContent(i),evhist_sam0_data->GetBinContent(i));
            hChi2->Fill(i-0.5,val);
            hRatio->Fill(i-0.5,evhist_sam0_pred_final->GetBinContent(i)/evhist_sam0_data->GetBinContent(i));

            hChi2_R->Fill(R,val);
            hChi2_costh->Fill(costh,val);
            hChi2_cosths->Fill(cosths,val);
            hChi2_zpos->Fill(zpos,val);
        }
    }

    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    hChi2->Draw("hist");
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    c1->SaveAs(Form("%s_chi2.pdf",figname.c_str()));
    hRatio->Draw("hist");
    c1->SaveAs(Form("%s_ratio.pdf",figname.c_str()));
    hChi2_R->Draw("colz");
    c1->SaveAs(Form("%s_chi2_R.pdf",figname.c_str()));
    hChi2_costh->Draw("colz");
    c1->SaveAs(Form("%s_chi2_costh.pdf",figname.c_str()));
    hChi2_cosths->Draw("colz");
    c1->SaveAs(Form("%s_chi2_cosths.pdf",figname.c_str()));
    hChi2_zpos->Draw("colz");
    c1->SaveAs(Form("%s_chi2_zpos.pdf",figname.c_str()));
}

void get_plot_theta_phi(std::string filename, std::string var)
{
    gStyle->SetOptStat(0);

    TFile* f = new TFile(filename.c_str());

    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    TH1D* hist_LED_theta_result = (TH1D*)f->Get(Form("hist_LED_%s_result",var.c_str()));
    TH1D* hist_LED_theta_prior = (TH1D*)f->Get(Form("hist_LED_%s_prior",var.c_str()));
    TH1D* hist_LED_theta_error_final = (TH1D*)f->Get(Form("hist_LED_%s_error_final",var.c_str()));
    TH1D* hist_LED_theta_error_prior = (TH1D*)f->Get(Form("hist_LED_%s_error_prior",var.c_str()));
    TH1D* hist_LED_theta_prefit = new TH1D("","",40,-40,0);
    TH1D* hist_LED_theta_postfit = new TH1D("","",40,-40,0);
    double chi2syst = 0;
    for (int i=1;i<=hist_LED_theta_prefit->GetNbinsX();i++)
    {
        hist_LED_theta_prefit->SetBinContent(i,hist_LED_theta_prior->GetBinContent(i));
        hist_LED_theta_prefit->SetBinError(i,hist_LED_theta_error_prior->GetBinContent(i));
        hist_LED_theta_postfit->SetBinContent(i,hist_LED_theta_result->GetBinContent(i));
        hist_LED_theta_postfit->SetBinError(i,hist_LED_theta_error_final->GetBinContent(i));

        double val = hist_LED_theta_result->GetBinContent(i)-hist_LED_theta_prior->GetBinContent(i);
        val *= 1./hist_LED_theta_error_prior->GetBinContent(i);
        chi2syst += val*val;
    }
    hist_LED_theta_postfit->GetXaxis()->SetTitle("-#theta_{s} (deg)");
    if (var == "theta") hist_LED_theta_postfit->GetYaxis()->SetTitle("P(#theta_{s})");
    else if (var == "phi") hist_LED_theta_postfit->GetYaxis()->SetTitle("V(#theta_{s})");
    hist_LED_theta_postfit->GetYaxis()->SetTitleOffset(1.2);
    if (var == "theta") hist_LED_theta_postfit->SetMaximum(1.1);
    else if (var == "phi") hist_LED_theta_postfit->SetMaximum(0.025);
    hist_LED_theta_prefit->SetLineColor(kRed);
    // hist_LED_theta_prefit->SetFillColor(kRed);
    // hist_LED_theta_prefit->Draw("E3");
    hist_LED_theta_postfit->SetLineColor(kBlue);
    hist_LED_theta_postfit->Draw("E0");
    for (int i=1;i<=hist_LED_theta_prefit->GetNbinsX();i++)
    {
        if (hist_LED_theta_prefit->GetBinContent(i)>0)
        {
        TBox* box = new TBox(hist_LED_theta_prefit->GetBinLowEdge(i),hist_LED_theta_prefit->GetBinContent(i)-hist_LED_theta_prefit->GetBinError(i),
                                hist_LED_theta_prefit->GetBinLowEdge(i+1),hist_LED_theta_prefit->GetBinContent(i)+hist_LED_theta_prefit->GetBinError(i));
        box->SetFillColor(kRed);
        box->Draw("same");
        }
    }
    hist_LED_theta_postfit->Draw("E0 same");
    TLegend* legend;
    if (var == "theta") legend = new TLegend(0.5,0.2,0.8,0.4);
    else if (var == "phi") legend = new TLegend(0.5,0.6,0.8,0.8);
    legend->AddEntry(hist_LED_theta_prefit,"Pre-fit","lep");
    legend->AddEntry(hist_LED_theta_postfit,"Post-fit","lep");
    TH1D* chi2_syst_periter = (TH1D*)f->Get("chi2_syst_periter");
    //double chi2syst = chi2_syst_periter->GetBinContent(chi2_syst_periter->GetNbinsX());
    legend->AddEntry((TObject*)0,Form("#chi^{2}_{syst}=%2.1f",chi2syst),"");
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->Draw();
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    if (var == "theta") c1->SaveAs(Form("%s_thetas.pdf",figname.c_str()));
    else if (var == "phi") c1->SaveAs(Form("%s_phis.pdf",figname.c_str()));
}

void plot_postfit_pol()
{
    gStyle->SetPaintTextFormat("4.4f");

    //get_plot("diffuser4_400nm_nominal_mPMT_pol.root",10648,1);
    //get_plot("diffuser4_400nm_nominal_BnL_pol.root",10648,0);
    //get_plot_combined("diffuser4_400nm_nominal_combined_pol.root",10648);
    //get_plot_3angular("diffuser4_400nm_nominal_mPMT_pol_3angular.root",10648);
    //get_plot("diffuser4_400nm_nominal_BnL_pol_indirect.root",10648,0);
    //get_plot("diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root",10648,1);
    //get_plot("diffuser4_350nm_nominal_mPMT_pol_3angular_indirect.root",6351,1);
    //get_plot("diffuser4_350nm_nominal_BnL_pol_indirect.root",6351,0);

    //get_plot_compare_nominal("diffuser4_400nm_nominal_BnL_pol_indirect_led_uniform.root","diffuser4_400nm_nominal_BnL_pol_indirect.root",10648,0,7);
    //get_plot_compare_nominal("diffuser4_400nm_nominal_mPMT_pol_3angular_indirect_led_uniform.root","diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root",10648,1,24);

    //get_plot_compare_nominal("diffuser4_400nm_nominal_BnL_pol_indirect_led_fixedprior.root","diffuser4_400nm_nominal_BnL_pol_indirect.root",10648,0,7);
    //get_plot_compare_nominal("diffuser4_400nm_nominal_mPMT_pol_3angular_indirect_led_fixedprior.root","diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root",10648,1,24);

    // get_plot_compare_nominal("diffuser4_400nm_nominal_BnL_pol_indirect_led_5pc_50pc_prior.root","diffuser4_400nm_nominal_BnL_pol_indirect.root",10648,0,7);
    // get_plot_theta_phi("diffuser4_400nm_nominal_BnL_pol_indirect_led_5pc_50pc_prior.root","theta");
    // get_plot_theta_phi("diffuser4_400nm_nominal_BnL_pol_indirect_led_5pc_50pc_prior.root","phi");

    // get_plot_compare_nominal("diffuser4_400nm_nominal_mPMT_pol_3angular_indirect_led_5pc_50pc_prior.root","diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root",10648,1,24);
    // get_plot_theta_phi("diffuser4_400nm_nominal_mPMT_pol_3angular_indirect_led_5pc_50pc_prior.root","theta");
    // get_plot_theta_phi("diffuser4_400nm_nominal_mPMT_pol_3angular_indirect_led_5pc_50pc_prior.root","phi");

    // get_plot_compare_nominal_joint("diffuser4_400nm_nominal_combined_pol_indirect_binning.root","diffuser4_400nm_nominal_BnL_pol_indirect.root","diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root",10648);
    // get_plot_theta_phi("diffuser4_400nm_nominal_combined_pol_indirect_binning.root","theta");
    // get_plot_theta_phi("diffuser4_400nm_nominal_combined_pol_indirect_binning.root","phi");

    //get_plot_compare_nominal("diffuser4_400nm_z-0.15_BnL.root","diffuser4_400nm_nominal_BnL_pol_indirect.root",10648,0,7);
    //get_plot_compare_nominal("diffuser4_400nm_z-0.15_mPMT.root","diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root",10648,1,24);

    // get_plot_attenz_compare_nominal("diffuser4_400nm_z-0.15_BnL_fixedbeta.root", "diffuser4_400nm_nominal_BnL_pol_indirect.root",10648, -0.15, 0, 7);
    // get_plot_attenz_compare_nominal("diffuser4_400nm_z-0.15_mPMT_fixedbeta.root", "diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root",10648, -0.15, 1, 24);

    // get_plot_attenz_compare_nominal("diffuser4_400nm_z-0.15_BnL_attenz.root", "diffuser4_400nm_nominal_BnL_pol_indirect.root",10648, -0.15, 0, 8);
    // get_plot_attenz_compare_nominal("diffuser4_400nm_z-0.15_mPMT_attenz.root", "diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root",10648, -0.15, 1, 25);

    // get_plot_attenz_compare_nominal("diffuser4_20_400nm_z-0.15_BnL_attenz.root", "diffuser4_400nm_nominal_BnL_pol_indirect.root",10648, -0.15, 0, 8);
    // get_plot_attenz_compare_nominal("diffuser4_20_400nm_z-0.15_mPMT_attenz.root", "diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root",10648, -0.15, 1, 25);

    //get_plot_compare_nominal("diffuser4_400nm_linz5pc_BnL.root","diffuser4_400nm_nominal_BnL_pol_indirect.root",10648,0,7);
    //get_plot_compare_nominal("diffuser4_400nm_linz2pc_mPMT.root","diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root",10648,1,24);

    // get_plot_linz_compare_nominal("diffuser4_400nm_linz5pc_BnL_fit.root","diffuser4_400nm_nominal_BnL_pol_indirect.root",10648,0,7);
    // get_plot_linz_compare_nominal("diffuser4_400nm_linz2pc_mPMT_fit.root","diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root",10648,1,24);
    // get_plot_linz_compare_nominal("diffuser4_20_400nm_linz5pc_BnL_fit.root","diffuser4_400nm_nominal_BnL_pol_indirect.root",10648,0,7);
    // get_plot_linz_compare_nominal("diffuser4_20_400nm_linz2pc_mPMT_fit.root","diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root",10648,1,24);

    //get_plot_chi2("../test.root");
    get_plot_attenz_compare_nominal("diffuser4_400nm_z-0.15_r_2.5pc_mPMT.root", "diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root",10648, -0.15, 1, 25);
    get_plot_attenz_compare_nominal("diffuser4_400nm_z-0.15_r_5pc_mPMT.root", "diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root",10648, -0.15, 1, 25);
    get_plot_attenz_compare_nominal("diffuser4_400nm_z-0.15_r_10pc_mPMT.root", "diffuser4_400nm_nominal_mPMT_pol_3angular_indirect.root",10648, -0.15, 1, 25);
}
