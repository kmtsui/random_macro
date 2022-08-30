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

TH1D* get_plot(std::string filename, double truth_alpha, int pmt_type = 0, int startingIndex = 2){
    gStyle->SetOptStat(0);

    int orders_BnL[] = {2,3};
    double ranges_BnL[] = {0.5,0.6,1.0};
    int orders_mPMT[] = {2,3,3};
    double ranges_mPMT[] = {0.,0.6,0.75,1.0};
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
        npars+=pol_orders[1]-1;
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
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX()-1);
    int ndof = ((TTree*)f->Get("PMTTree"))->GetEntries()-npars-1;
    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    h_fit_pol->Draw("");
    double fitted_alpha = ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);
    TLatex latex;
    latex.SetTextSize(0.045);
    double miny = h_fit_pol->GetMinimum()*1.1;
    double maxy = h_fit_pol->GetMaximum();
    double diff = (maxy-miny)/10;
    latex.DrawLatex(0.76,miny+diff*2,Form("#chi^{2}/dof = %4.0f/%i",chi2,ndof));
    latex.DrawLatex(0.76,miny+diff,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha));
    latex.DrawLatex(0.76,miny,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha));
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
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX()-1);
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
    latex.DrawLatex(0.75,miny+0.02,Form("#chi^{2}/dof = %4.0f/%i",chi2,ndof));
    latex.DrawLatex(0.75,miny+0.01,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha));
    latex.DrawLatex(0.75,miny,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha));
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
    TH1D* h_mPMT = get_plot(filename,truth_alpha,1,7);

    TFile* f = new TFile(filename.c_str());
    TH1D* chi2_total_periter = (TH1D*)f->Get("chi2_total_periter");
    double chi2 = chi2_total_periter->GetBinContent(chi2_total_periter->GetNbinsX()-1);
    int npars_BnL = 5;
    int npars_mPMT = 7;
    int npars = npars_BnL+npars_mPMT;
    int ndof = ((TTree*)f->Get("PMTTree"))->GetEntries()-1;
    double fitted_alpha = ((TH1D*)f->Get("hist_alpha_result"))->GetBinContent(1);
    double error_alpha = ((TH1D*)f->Get("hist_alpha_error_final"))->GetBinContent(1);

    TCanvas* c1 = new TCanvas();
    c1->SetGridy();
    h_mPMT->SetLineColor(kBlue);
    h_mPMT->Draw();
    h_BnL->SetLineColor(kRed);
    h_BnL->Draw("same");
    TLegend* legend = new TLegend(0.5,0.15,0.9,0.45);
    legend->AddEntry(h_BnL,"A_{BnL}","lep");
    legend->AddEntry(h_mPMT,"A_{mPMT}","lep");
    legend->AddEntry((TObject*)0,Form("#chi^{2}/dof = %4.0f/%i",chi2,ndof),"");
    legend->AddEntry((TObject*)0,Form("Truth L_{#alpha} = %4.0f cm",truth_alpha),"");
    legend->AddEntry((TObject*)0,Form("Fitted L_{#alpha} = %3.0f +/- %2.0f cm",fitted_alpha,error_alpha),"");
    legend->Draw();
    std::string figname = filename;
    figname.erase(figname.find(".root"),5);
    c1->SaveAs(Form("%s.pdf",figname.c_str()));

    gStyle->SetPaintTextFormat("4.2f");
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

void plot_postfit_pol()
{
    //plot_angular_response_pol_BnL();

    // get_plot("diffusr4_400nm_nominal_mPMT_pol.root",10648,1);
    // get_plot("diffusr4_350nm_nominal_mPMT_pol.root",6351,1);
    // get_plot("diffusr4_400nm_nominal_BnLPMT_pol.root",10648,0);
    // get_plot("diffusr4_350nm_nominal_BnLPMT_pol.root",6351,0);
    
    //get_plot_3angular("diffusr4_400nm_nominal_mPMT_pol_3angular.root",10648);

    //get_plot_combined("diffusr4_400nm_nominal_combined_pol.root",10648);

    //get_plot_led("diffusr4_400nm_nominal_BnL_pol_ledprofile.root",10648,0);
    //get_plot_led("diffusr4_400nm_nominal_mPMT_pol_ledprofile.root",10648,1);

    // get_plot_led("diffusr4_400nm_nominal_mPMT_pol_3angular_indirect_led1percent.root",10648,1);
    // get_plot_led("diffusr4_400nm_nominal_mPMT_pol_3angular_indirect_led5percent.root",10648,1);
    // get_plot_led("diffusr4_400nm_nominal_mPMT_pol_3angular_indirect_led10percent.root",10648,1);
    // get_plot_led("diffusr4_400nm_nominal_BnL_pol_indirect_led1percent.root",10648,0);
    // get_plot_led("diffusr4_400nm_nominal_BnL_pol_indirect_led5percent.root",10648,0);
    // get_plot_led("diffusr4_400nm_nominal_BnL_pol_indirect_led10percent.root",10648,0);

    TH1D* hist_full = get_plot("diffusr4_400nm_nominal_mPMT_pol_ledprofile_fullring.root",10648,1);
    TH1D* hist_noring1 = get_plot("diffusr4_400nm_nominal_mPMT_pol_ledprofile_noring1.root",10648,1);
    TH1D* hist_noring2 = get_plot("diffusr4_400nm_nominal_mPMT_pol_ledprofile_noring2.root",10648,1);
    TCanvas* c1 = new TCanvas();
    hist_full->SetLineColor(kBlack);
    hist_noring1->SetLineColor(kBlue);
    hist_noring2->SetLineColor(kRed);
    hist_full->Draw("same");
    hist_noring1->Draw("same");
    hist_noring2->Draw("same");
    c1->SaveAs("test.pdf");
}
