const int MAXNPAR = 1000;
const int MAXNSCNDPART = 1000;

//#############################################################################
void SetStyleVariables(TStyle *t2kStyle){

  // use plain black on white colors
  t2kStyle->SetFrameBorderMode(0);
  t2kStyle->SetCanvasBorderMode(0);
  t2kStyle->SetPadBorderMode(0);
  t2kStyle->SetPadColor(0);
  t2kStyle->SetCanvasColor(0);
  t2kStyle->SetStatColor(0);
  //t2kStyle->SetFillColor(0);
  t2kStyle->SetLegendBorderSize(1);

  // set the paper & margin sizes
  t2kStyle->SetPaperSize(20,26);
  t2kStyle->SetPadTopMargin(0.1);
  t2kStyle->SetPadRightMargin(0.05); //0.15
  t2kStyle->SetPadBottomMargin(0.15);
  t2kStyle->SetPadLeftMargin(0.15);

  t2kStyle->SetCanvasDefH(550);

  // use large Times-Roman fonts
  t2kStyle->SetTextFont(132);
  t2kStyle->SetTextSize(0.08);
  t2kStyle->SetLabelFont(132,"x");
  t2kStyle->SetLabelFont(132,"y");
  t2kStyle->SetLabelFont(132,"z");
  t2kStyle->SetLabelSize(0.06,"x");
  t2kStyle->SetTitleSize(0.072,"x");
  t2kStyle->SetTitleOffset(0.88,"x");
  t2kStyle->SetLabelSize(0.05,"y");
  t2kStyle->SetTitleSize(0.08,"y");
  t2kStyle->SetTitleOffset(0.7,"y");
  t2kStyle->SetLabelSize(0.06,"z");
  t2kStyle->SetTitleSize(0.06,"z");
  t2kStyle->SetLabelFont(132,"t");
  t2kStyle->SetTitleFont(132,"x");
  t2kStyle->SetTitleFont(132,"y");
  t2kStyle->SetTitleFont(132,"z");
  t2kStyle->SetTitleFont(132,"t");
  t2kStyle->SetTitleFillColor(0);
  t2kStyle->SetTitleX(0.25);
  t2kStyle->SetTitleFontSize(0.08);
  t2kStyle->SetTitleFont(132,"pad");

  //t2kStyle->SetPadGridX(true);
  //t2kStyle->SetPadGridY(true);

  // use bold lines and markers
  //  t2kStyle->SetMarkerStyle(20);
  t2kStyle->SetHistLineWidth(1.85);
  t2kStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars and y error bar caps
  //  t2kStyle->SetErrorX(0.001);

  // do not display any of the standard histogram decorations
  t2kStyle->SetOptTitle(0);
  t2kStyle->SetOptStat(0);
  t2kStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  //t2kStyle->SetPadTickX(1);
  t2kStyle->SetPadTickY(1);

  t2kStyle->SetPalette(1,0);  // use the nice red->blue palette
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  //Double_t green[NRGBs] = { 0.00, 0.00, 0.00, 0.00, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue,
                                   NCont);
  t2kStyle->SetNumberContours(NCont);

  // End of definition of t2kStyle
}

const int nGraphs=3;

void get_plots(const char* input , TH1D* hGraph[nGraphs]){
//TH1D* get_plots(const char* input){
    TFile* f = new TFile(input);
    TTree* t = (TTree*)f->Get("h1");

    Int_t npar;
    UChar_t ipv[MAXNPAR]; //This is GEANT3 particle code of primary particles.
    Int_t nscndprt;
    Int_t iprtscnd[MAXNSCNDPART]; //This is PDG code of secondary particles
    Float_t pscnd[MAXNSCNDPART][3];
    Int_t lmecscnd[MAXNSCNDPART]; //This is GEANT3 interaction code
    t->SetBranchAddress("npar",&npar);
    t->SetBranchAddress("ipv",&ipv);
    t->SetBranchAddress("nscndprt",&nscndprt);
    t->SetBranchAddress("iprtscnd",&iprtscnd);
    t->SetBranchAddress("pscnd",&pscnd);
    t->SetBranchAddress("lmecscnd",&lmecscnd);

    hGraph[0] = new TH1D("","",10,0,10);//multiplicity
    //std::cout<<"hGraph[0]->GetNbinsX()="<<hGraph[0]->GetNbinsX()<<std::endl;
    hGraph[1] = new TH1D("","",40,0,10);//hGammaEnergy
    hGraph[2] = new TH1D("","",40,0,10);//hTotalGammaEnergy

    for (int i=0;i<t->GetEntries();i++){
        t->GetEntry(i);
        //std::cout<<"Event "<<i<<", number of primaries = "<<npar;
        //for (int j=0;j<npar;j++)
            //std::cout<<", ipv["<<j<<"]="<<(int)ipv[j];
        //std::cout<<std::endl;
        int nGam = 0;
        double totalGammaEnergy=0;
        //std::cout<<"Event "<<i<<", nCapture produces gammas of energy ";
        for (int j=0;j<nscndprt;j++) {
            if (iprtscnd[j]==22&&lmecscnd[j]==18) {
                double energy=sqrt(pscnd[j][0]*pscnd[j][0]+pscnd[j][1]*pscnd[j][1]+pscnd[j][2]*pscnd[j][2]);
                if (energy>2.224&&energy<2.2246) continue;
                totalGammaEnergy+=energy;
                //std::cout<<energy<<" MeV, ";
                nGam++;
                hGraph[1]->Fill(energy);
            }
        }
        //std::cout<<"total gamma energy = "<<totalGammaEnergy<<" MeV "<<std::endl;
        if (totalGammaEnergy>0) {
            hGraph[0]->Fill(nGam);
            hGraph[2]->Fill(totalGammaEnergy);
        }
    }
    //return hGraph;
    /*TCanvas* c1 = new TCanvas();
    hMultiplicity->Draw();
    c1->SaveAs("plots/multiplicity.pdf");
    hGammaEnergy->Draw();
    c1->SaveAs("plots/Egam.pdf");
    hTotalGammaEnergy->Draw();
    c1->SaveAs("plots/Egamtotal.pdf");*/

}

void analysis_skdetsim(){
    TStyle *t2kstyle = new TStyle("T2K","T2K approved plots style");
    SetStyleVariables(t2kstyle);
    gROOT->SetStyle("T2K");

    const int nFiles=5;
    std::string fileList[] = {
        "test_0.root","test_4.root","test_2.root","test_3.root","test_1.root",
    };
    std::string modelName[] = {
        "geant4 default", "geant4 photo-evaporation", "GLG4Sim", "ggarnet", "ANNRI_Gd"
    };
    std::string varname[] = {
        "Gamma multiplicity", "Gamma energy (MeV)", "Total gamma energy (MeV)"
    };
    std::string graphname[] = {
        "plots/multiplicity.pdf", "plots/Egam.pdf", "plots/Egamtotal.pdf"
    };
    int linecolor[]={kBlack, kRed, kBlue, kGreen, kViolet};

    TH1D* hGraph[nFiles][nGraphs];
    /*std::vector<TH1D*> hGraph = std::vector<TH1D*>();
    for (int i=0;i<nFiles;i++){
        TH1D* hGraph_tmp[3];
        hGraph.push_back(hGraph_tmp);
    }*/

    for (int i=0;i<nFiles;i++) {
        get_plots(fileList[i].c_str(),hGraph[i]);
    }

    for (int i=0;i<nGraphs;i++) {
        TCanvas* c1 = new TCanvas();
        double maxy = 0;
        for (int j=0;j<nFiles;j++) {
            hGraph[j][i]->Scale(1/hGraph[j][i]->Integral());
            if (hGraph[j][i]->GetMaximum()>maxy) maxy = hGraph[j][i]->GetMaximum();
        }
        for (int j=nFiles-1;j>=0;j--) {
            hGraph[j][i]->SetMaximum(maxy*1.1);
            hGraph[j][i]->SetLineColor(linecolor[j]);
            hGraph[j][i]->SetLineWidth(3);
            hGraph[j][i]->SetLineStyle(j+1);
            hGraph[j][i]->GetXaxis()->SetTitle(varname[i].c_str());
            //hGraph[j][i]->Scale(1/hGraph[j][i]->Integral());
            hGraph[j][i]->GetYaxis()->SetTitle("pdf");
            hGraph[j][i]->Draw("same");
        }
        TLegend* legend2;
        if (i==2) legend2 = new TLegend(0.3,0.59,0.7,0.88);
        else if (i==1) legend2 = new TLegend(0.42,0.59,0.82,0.88);
        else if (i==0) legend2 = new TLegend(0.2,0.59,0.6,0.82);
        for (int j=0;j<nFiles;j++)
            legend2->AddEntry(hGraph[j][i],modelName[j].c_str(),"l");
        legend2->SetBorderSize(0);
        legend2->SetFillStyle(0);
        legend2->Draw();
        c1->SaveAs(graphname[i].c_str());
    }
}