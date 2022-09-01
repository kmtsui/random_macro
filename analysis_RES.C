//#include "../src/event1.h"

Float_t CalcDeltapTT(const Float_t nuDir[3], const Float_t muDir[3], const Float_t pi_dir[3], const Float_t pi_mom, const Float_t p_dir[3], const Float_t p_mom){
    Float_t transDir[3] = {muDir[1]*nuDir[2]-muDir[2]*nuDir[1],muDir[2]*nuDir[0]-muDir[0]*nuDir[2],muDir[0]*nuDir[1]-muDir[1]*nuDir[0]};
    Float_t deltapTT=((pi_dir[0]*pi_mom+p_dir[0]*p_mom)*transDir[0]+(pi_dir[1]*pi_mom+p_dir[1]*p_mom)*transDir[1]+(pi_dir[2]*pi_mom+p_dir[2]*p_mom)*transDir[2])/sqrt(transDir[0]*transDir[0]+transDir[1]*transDir[1]+transDir[2]*transDir[2]);
    return deltapTT;
};

Float_t CalcW_p_pi(const Float_t nuDir[3], const Float_t muDir[3], const Float_t mu_mom, const Float_t pi_dir[3], const Float_t pi_mom, const Float_t p_dir[3], const Float_t p_mom){
    Double_t ProtonMass = 938.272;
    Double_t PionMass = 139.570;
    Double_t pi_mom_x = pi_mom*pi_dir[0];
    Double_t pi_mom_y = pi_mom*pi_dir[1];
    Double_t pi_mom_z = pi_mom*pi_dir[2];
    Double_t p_mom_x = p_mom*p_dir[0];
    Double_t p_mom_y = p_mom*p_dir[1];
    Double_t p_mom_z = p_mom*p_dir[2];
    Double_t Epi = sqrt(PionMass*PionMass+pi_mom_x*pi_mom_x+pi_mom_y*pi_mom_y+pi_mom_z*pi_mom_z);
    Double_t Ep = sqrt(ProtonMass*ProtonMass+p_mom_x*p_mom_x+p_mom_y*p_mom_y+p_mom_z*p_mom_z);
    Double_t W = sqrt((Epi+Ep)*(Epi+Ep)-(pi_mom_x+p_mom_x)*(pi_mom_x+p_mom_x)-(pi_mom_y+p_mom_y)*(pi_mom_y+p_mom_y)-(pi_mom_z+p_mom_z)*(pi_mom_z+p_mom_z));

    return W;
};

Float_t CalcpW_nu_mu(const Float_t nuDir[3], const Float_t muDir[3], const Float_t mu_mom, const Float_t pi_dir[3], const Float_t pi_mom, const Float_t p_dir[3], const Float_t p_mom){
    Double_t MuonMass =  105.65837;// }//in MeV //google = wiki
    Double_t NeutronMass = 939.565;//}//in GeV //wiki
    Double_t ProtonMass = 938.272;//}//in GeV //google = wiki
    Double_t PionMass = 139.570;//}//in GeV //wiki
    Double_t MA = 11.174860*1000;
    Double_t MAstar = 10.263718*1000;

    TVector3 pp(p_mom*p_dir[0],p_mom*p_dir[1],p_mom*p_dir[2]);
    TVector3 ppi(pi_mom*pi_dir[0],pi_mom*pi_dir[1],pi_mom*pi_dir[2]);
    TVector3 pmu(mu_mom*muDir[0],mu_mom*muDir[1],mu_mom*muDir[2]);
    TVector3 pnu(nuDir[0],nuDir[1],nuDir[2]);

    TVector3 plmu = (pmu*pnu.Unit())*pnu.Unit();
    TVector3 plp = (pp*pnu.Unit())*pnu.Unit();
    TVector3 plpi = (ppi*pnu.Unit())*pnu.Unit();

    TVector3 plbaryon=plp+plpi;

    TVector3 vdPt = pmu.Cross(pnu.Unit())+pp.Cross(pnu.Unit())+ppi.Cross(pnu.Unit());

    const double pT = vdPt.Mag();
    const double Eprim  = sqrt(pmu.Mag()*pmu.Mag()+MuonMass*MuonMass);
    const double kprimL = plmu.Mag();
    const double Epprim = sqrt(pp.Mag()*pp.Mag()+ProtonMass*ProtonMass)+sqrt(ppi.Mag()*ppi.Mag()+PionMass*PionMass);
    const double pprimL = plbaryon.Mag();
    const double factor = MA - Eprim - Epprim + kprimL + pprimL;
    const double pL = -(MAstar*MAstar + pT*pT-factor*factor)/2.0/factor;
    Float_t neutronmomentum = sqrt(pL*pL + pT*pT);

    double Enu = kprimL+pprimL-pL;
    TVector3 pN = pL*pnu.Unit()+vdPt;
    TVector3 pW = pnu*Enu-plmu+pL*pnu.Unit()+pp.Cross(pnu.Unit())+ppi.Cross(pnu.Unit());
    double W = (Enu+MA-Eprim-sqrt(MAstar*MAstar+neutronmomentum*neutronmomentum));
    W=sqrt(W*W-pW.Mag()*pW.Mag());

    return W;
};

Float_t CalcpN(const Float_t nuDir[3], const Float_t muDir[3], const Float_t mu_mom, const Float_t pi_dir[3], const Float_t pi_mom, const Float_t p_dir[3], const Float_t p_mom){
    Double_t MuonMass =  105.65837;// }//in MeV //google = wiki
    Double_t NeutronMass = 939.565;//}//in GeV //wiki
    Double_t ProtonMass = 938.272;//}//in GeV //google = wiki
    Double_t PionMass = 139.570;//}//in GeV //wiki
    Double_t MA = 11.174860*1000;
    Double_t MAstar = 10.263718*1000;

    TVector3 pp(p_mom*p_dir[0],p_mom*p_dir[1],p_mom*p_dir[2]);
    TVector3 ppi(pi_mom*pi_dir[0],pi_mom*pi_dir[1],pi_mom*pi_dir[2]);
    TVector3 pmu(mu_mom*muDir[0],mu_mom*muDir[1],mu_mom*muDir[2]);
    TVector3 pnu(nuDir[0],nuDir[1],nuDir[2]);

    TVector3 plmu = (pmu*pnu.Unit())*pnu.Unit();
    TVector3 plp = (pp*pnu.Unit())*pnu.Unit();
    TVector3 plpi = (ppi*pnu.Unit())*pnu.Unit();

    TVector3 plbaryon=plp+plpi;

    TVector3 vdPt = pmu.Cross(pnu.Unit())+pp.Cross(pnu.Unit())+ppi.Cross(pnu.Unit());

    const double pT = vdPt.Mag();
    const double Eprim  = sqrt(pmu.Mag()*pmu.Mag()+MuonMass*MuonMass);
    const double kprimL = plmu.Mag();
    const double Epprim = sqrt(pp.Mag()*pp.Mag()+ProtonMass*ProtonMass)+sqrt(ppi.Mag()*ppi.Mag()+PionMass*PionMass);
    const double pprimL = plbaryon.Mag();
    const double factor = MA - Eprim - Epprim + kprimL + pprimL;
    const double pL = -(MAstar*MAstar + pT*pT-factor*factor)/2.0/factor;
    Float_t neutronmomentum = sqrt(pL*pL + pT*pT);

    return neutronmomentum;
};

Float_t CalcdeltaAlphaT(const Float_t nuDir[3], const Float_t muDir[3], const Float_t mu_mom, const Float_t pi_dir[3], const Float_t pi_mom, const Float_t p_dir[3], const Float_t p_mom){
    TVector3 pp(p_mom*p_dir[0],p_mom*p_dir[1],p_mom*p_dir[2]);
    TVector3 ppi(pi_mom*pi_dir[0],pi_mom*pi_dir[1],pi_mom*pi_dir[2]);
    TVector3 pmu(mu_mom*muDir[0],mu_mom*muDir[1],mu_mom*muDir[2]);
    TVector3 pnu(nuDir[0],nuDir[1],nuDir[2]);

    TVector3 pTmu = pmu.Cross(pnu.Unit());
    TVector3 vdPt = pmu.Cross(pnu.Unit())+pp.Cross(pnu.Unit())+ppi.Cross(pnu.Unit());

    Float_t deltaAlphaT=-pTmu.Unit()*vdPt.Unit();
    if (deltaAlphaT>-1&&deltaAlphaT<1) return acos(deltaAlphaT)*180./3.14159265;
    else return -999;

};

Float_t CalcdeltaPhiT(const Float_t nuDir[3], const Float_t muDir[3], const Float_t mu_mom, const Float_t pi_dir[3], const Float_t pi_mom, const Float_t p_dir[3], const Float_t p_mom){
    TVector3 pp(p_mom*p_dir[0],p_mom*p_dir[1],p_mom*p_dir[2]);
    TVector3 ppi(pi_mom*pi_dir[0],pi_mom*pi_dir[1],pi_mom*pi_dir[2]);
    TVector3 pmu(mu_mom*muDir[0],mu_mom*muDir[1],mu_mom*muDir[2]);
    TVector3 pnu(nuDir[0],nuDir[1],nuDir[2]);

    TVector3 pTmu = pmu.Cross(pnu.Unit());
    TVector3 vdPt = pp.Cross(pnu.Unit())+ppi.Cross(pnu.Unit());

    Float_t deltaAlphaT=-pTmu.Unit()*vdPt.Unit();
    if (deltaAlphaT>=-1&&deltaAlphaT<=1) return acos(deltaAlphaT)*180./3.14159265;
    else return -999;

};

Float_t CalcDpt(const Float_t nuDir[3], const Float_t muDir[3], const Float_t mu_mom, const Float_t pi_dir[3], const Float_t pi_mom, const Float_t p_dir[3], const Float_t p_mom){
    TVector3 pp(p_mom*p_dir[0],p_mom*p_dir[1],p_mom*p_dir[2]);
    TVector3 ppi(pi_mom*pi_dir[0],pi_mom*pi_dir[1],pi_mom*pi_dir[2]);
    TVector3 pmu(mu_mom*muDir[0],mu_mom*muDir[1],mu_mom*muDir[2]);
    TVector3 pnu(nuDir[0],nuDir[1],nuDir[2]);

    TVector3 vdPt = pmu.Cross(pnu.Unit())+pp.Cross(pnu.Unit())+ppi.Cross(pnu.Unit());

    return vdPt.Mag();

};

void analysis_getGraph(const char* input){
  ROOT::EnableImplicitMT(8); // Enable implicit parallel computing

  TFile *f = new TFile(input);
  TTree *t = (TTree*)f->Get("treeout");
  //create a pointer to event
  event *e   = new event();
  t->SetBranchAddress("e",&e);
  TH1D* xsections = (TH1D*)f->Get("xsections");

  //string out_reweighted(input);
  //out_reweighted = out_reweighted+".reweighted.weights";
  //TFile *ff = new TFile(out_reweighted.c_str());
  //TTree *tf = (TTree*)ff->Get("weights");
  //Double_t reweights;
  //tf->SetBranchAddress("weight",&reweights);

  TH1D* fluxhist;
  //int nfluxbins;
  double fluxbins[1000];
  TFile* fflux = new TFile("/hepstore/kmtsui/T2K/work/xsLLhFitter_super/inputs/nd5_tuned13av7p1_13anom_run1-10b_numode.root");
  fluxhist = (TH1D*)fflux->Get("enu_nd5_tuned13a_numu");
  const int nfluxbins = fluxhist->GetNbinsX();
  for (int i=1;i<=nfluxbins+1;i++){
    fluxbins[i-1] = fluxhist->GetBinLowEdge(i)*1000;
  }
  TH2D* rewhist1 = new TH2D("rewhist_Enu_Q2_C_res_pip","",nfluxbins,fluxbins,100,0,10);
  TH2D* rewhist2 = new TH2D("rewhist_Enu_k_C_res_pip","",nfluxbins,fluxbins,100,0,600);
  double Q2bins[201];
  double Wbins[201];
  double jacbins[201];
  double q0bins[201];
  for (int i=0;i<201;i++){
    Q2bins[i]=i*(10./200.);
    q0bins[i]=i*(2000./200.);
    jacbins[i]=i*(3.e6/200.);
    Wbins[i]=1080.+(1600.-1080.)/200.*i;
  }
  TH3D* rewhist3_11 = new TH3D("rewhist_Enu_Q2_W_C_res_p_pip","",nfluxbins,fluxbins,200,Q2bins,200,Wbins);
  TH3D* rewhist3_12 = new TH3D("rewhist_Enu_Q2_W_C_res_n_pi0","",nfluxbins,fluxbins,200,Q2bins,200,Wbins);
  TH3D* rewhist3_13 = new TH3D("rewhist_Enu_Q2_W_C_res_n_pip","",nfluxbins,fluxbins,200,Q2bins,200,Wbins);
  TH3D* rewhist31_11 = new TH3D("rewhist_Enu_Q2_q0_C_res_p_pip","",nfluxbins,fluxbins,200,Q2bins,200,q0bins);
  TH3D* rewhist31_12 = new TH3D("rewhist_Enu_Q2_q0_C_res_n_pi0","",nfluxbins,fluxbins,200,Q2bins,200,q0bins);
  TH3D* rewhist31_13 = new TH3D("rewhist_Enu_Q2_q0_C_res_n_pip","",nfluxbins,fluxbins,200,Q2bins,200,q0bins);
  TH3D* rewhist3[nfluxbins][3];
  TH3D* rewhist3_ppi[nfluxbins][3];
  const int nQ2bins = 100;
  TH3D* rewhist5[nfluxbins][nQ2bins][3];
  for (int i=0;i<nfluxbins;i++){
    rewhist3[i][0]= new TH3D(Form("rewhist_Enu%i_Q2_q0_W_C_res_p_pip",i),"",100,0,10,100,0,5000,100,1080,1600);
    rewhist3[i][1]= new TH3D(Form("rewhist_Enu%i_Q2_q0_W_C_res_n_pi0",i),"",100,0,10,100,0,5000,100,1080,1600);
    rewhist3[i][2]= new TH3D(Form("rewhist_Enu%i_Q2_q0_W_C_res_n_pip",i),"",100,0,10,100,0,5000,100,1080,1600);
    rewhist3_ppi[i][0]= new TH3D(Form("rewhist_Enu%i_Q2_q0_ppi_C_res_p_pip",i),"",100,0,10,100,0,5000,100,0,1500);
    rewhist3_ppi[i][1]= new TH3D(Form("rewhist_Enu%i_Q2_q0_ppi_C_res_n_pi0",i),"",100,0,10,100,0,5000,100,0,1500);
    rewhist3_ppi[i][2]= new TH3D(Form("rewhist_Enu%i_Q2_q0_ppi_C_res_n_pip",i),"",100,0,10,100,0,5000,100,0,100);
    for (int j=0;j<nQ2bins;j++){
      rewhist5[i][j][0]= new TH3D(Form("rewhist_Enu%i_Q2%i_q0_W_gam_C_res_p_pip",i,j),"",60,0,3000,40,1080,1600,20,1,2);
      rewhist5[i][j][1]= new TH3D(Form("rewhist_Enu%i_Q2%i_q0_W_gam_C_res_n_pi0",i,j),"",60,0,3000,40,1080,1600,20,1,2);
      rewhist5[i][j][2]= new TH3D(Form("rewhist_Enu%i_Q2%i_q0_W_gam_C_res_n_pip",i,j),"",60,0,3000,40,1080,1600,20,1,2);
    }
  }

  TH2D* rewhist11 = new TH2D("rewhist_Enu_Q2_C_res_p_pip","",nfluxbins,fluxbins,nQ2bins,0,10);
  TH2D* rewhist12 = new TH2D("rewhist_Enu_Q2_C_res_n_pi0","",nfluxbins,fluxbins,nQ2bins,0,10);
  TH2D* rewhist13 = new TH2D("rewhist_Enu_Q2_C_res_n_pip","",nfluxbins,fluxbins,nQ2bins,0,10);

  TH2D* rewhistEnuW11 = new TH2D("rewhist_Enu_W_C_res_p_pip","",nfluxbins,fluxbins,200,1080,1600);
  TH2D* rewhistEnuW12 = new TH2D("rewhist_Enu_W_C_res_n_pi0","",nfluxbins,fluxbins,200,1080,1600);
  TH2D* rewhistEnuW13 = new TH2D("rewhist_Enu_W_C_res_n_pip","",nfluxbins,fluxbins,200,1080,1600);

  std::string outname = "RES_";
  outname+=input;
  TFile* out_file = TFile::Open(outname.c_str(), "RECREATE");
  TTree* out_seltree = new TTree("RES","RES");
  Float_t dptt, pN, dat, dphit, dpt, W_p_pi, W_nu_mu;
  Float_t pmom, pcostheta, pimom, picostheta, mumom, mucostheta;
  Float_t W, Q2, costhetamupi, target_mom, weight, Enu;
  Float_t W_truth, Ehad;
  Double_t norm, res_nu_E, res_jacobian, q0;
  Int_t channel, pi_q, res_delta;
  Int_t target, pn, p_out, n_out;
  out_seltree -> Branch("dptt", &dptt, "dptt/F");
  out_seltree -> Branch("pN", &pN, "pN/F");
  out_seltree -> Branch("dat", &dat, "dat/F");
  out_seltree -> Branch("dphit", &dphit, "dphit/F");
  out_seltree -> Branch("dpt", &dpt, "dpt/F");
  out_seltree -> Branch("W_p_pi", &W_p_pi, "W_p_pi/F");
  out_seltree -> Branch("W_nu_mu", &W_nu_mu, "W_nu_mu/F");
  out_seltree -> Branch("channel", &channel, "channel/I");
  out_seltree -> Branch("pi_q", &pi_q, "pi_q/I");
  out_seltree -> Branch("target", &target, "target/I");
  out_seltree -> Branch("pn", &pn, "pn/I");
  out_seltree -> Branch("p_out", &p_out, "p_out/I");
  out_seltree -> Branch("n_out", &n_out, "n_out/I");
  out_seltree -> Branch("res_delta", &res_delta, "res_delta/I");
  out_seltree -> Branch("pmom", &pmom, "pmom/F");
  out_seltree -> Branch("pcostheta", &pcostheta, "pcostheta/F");
  out_seltree -> Branch("pimom", &pimom, "pimom/F");
  out_seltree -> Branch("picostheta", &picostheta, "picostheta/F");
  out_seltree -> Branch("mumom", &mumom, "mumom/F");
  out_seltree -> Branch("mucostheta", &mucostheta, "mucostheta/F");
  out_seltree -> Branch("W", &W, "W/F");
  out_seltree -> Branch("W_truth", &W_truth, "W_truth/F");
  out_seltree -> Branch("Q2", &Q2, "Q2/F");
  out_seltree -> Branch("costhetamupi", &costhetamupi, "costhetamupi/F");
  out_seltree -> Branch("target_mom", &target_mom, "target_mom/F");
  out_seltree -> Branch("weight", &weight, "weight/F");
  out_seltree -> Branch("norm", &norm, "norm/F");
  out_seltree -> Branch("res_nu_E", &res_nu_E, "res_nu_E/D");
  out_seltree -> Branch("res_jacobian", &res_jacobian, "res_jacobian/D");
  out_seltree -> Branch("q0", &q0, "q0/D");
  out_seltree -> Branch("Enu", &Enu, "Enu/F");
  out_seltree -> Branch("Ehad", &Ehad, "Ehad/F");

  const auto nEntries = t->GetEntries();
  t->SetImplicitMT(true);

  for (auto i : ROOT::TSeqUL(nEntries)){
    if (i%10000==0) {
      std::cout << "."  ;
      std::cout.flush();
    }
    t->GetEntry(i);
    p_out=0;n_out=0;
    //tf->GetEntry(i);
    if (e->dyn!=2) continue;
    res_nu_E = e->res_nu.momentum();
    channel = e->dyn;
    norm = e->norm;
    res_delta = e->flag.res_delta;
    pn = e->in[1].charge();
    Enu = e->in[0].momentum();
    target_mom = e->in[1].momentum();
    weight = e->weight;//*reweights;
    res_jacobian = e->res_jacobian;
    q0 = e->q0();
    Float_t nuDir[3], muDir[3], pi_dir[3], pi_mom, p_dir[3], p_mom;
    double Pmu=0, Ppi = 0, Pp=0;
    int nMu=0, nPi=0, nP=0, nOtherMesons=0;
    nuDir[0]=e->in[0].p().x/e->in[0].momentum();
    nuDir[1]=e->in[0].p().y/e->in[0].momentum();
    nuDir[2]=e->in[0].p().z/e->in[0].momentum();
    vect Wp4 = e->in[0].p4()+e->in[1].p4()-e->out[0].p4();
    Ehad = 0;
    W_truth = Wp4*Wp4;
    for (int k = 0; k < e->out.size(); k++) {
      if (k>0) Ehad+= e->out[k].E();
      if (abs(e->out[k].pdg)==13) { 
        nMu++;
        Pmu=e->out[k].momentum();
        muDir[0]=e->out[k].p().x/e->out[k].momentum();
        muDir[1]=e->out[k].p().y/e->out[k].momentum();
        muDir[2]=e->out[k].p().z/e->out[k].momentum();
      } else if (abs(e->out[k].pdg)==211||e->out[k].pdg==111) {
        nPi++;
        Ppi=e->out[k].momentum();
        pi_dir[0]=e->out[k].p().x/e->out[k].momentum();
        pi_dir[1]=e->out[k].p().y/e->out[k].momentum();
        pi_dir[2]=e->out[k].p().z/e->out[k].momentum();
        pi_q = e->out[k].charge();
      } 
      else if (e->out[k].pdg==2212) p_out++;
      else if (e->out[k].pdg==2112) n_out++;
    }
    /*for (int k = 0; k < e->post.size(); k++) {
      if (e->post[k].pdg==13) { 
        nMu++;
        Pmu=e->post[k].momentum();
        muDir[0]=e->post[k].p().x/e->post[k].momentum();
        muDir[1]=e->post[k].p().y/e->post[k].momentum();
        muDir[2]=e->post[k].p().z/e->post[k].momentum();
      } else if (abs(e->post[k].pdg)==211||e->post[k].pdg==111) {
        nPi++;
        Ppi=e->post[k].momentum();
        pi_dir[0]=e->post[k].p().x/e->post[k].momentum();
        pi_dir[1]=e->post[k].p().y/e->post[k].momentum();
        pi_dir[2]=e->post[k].p().z/e->post[k].momentum();
        pi_q = e->post[k].charge();
      }  
      else if (e->post[k].pdg==2212) p_out++;
      else if (e->post[k].pdg==2112) n_out++;
    }*/
    double costheat_mu = nuDir[0]*muDir[0]+nuDir[1]*muDir[1]+nuDir[2]*muDir[2];
    double costheat_pi = nuDir[0]*pi_dir[0]+nuDir[1]*pi_dir[1]+nuDir[2]*pi_dir[2];
    double costheat_p = nuDir[0]*p_dir[0]+nuDir[1]*p_dir[1]+nuDir[2]*p_dir[2];
    costhetamupi = muDir[0]*pi_dir[0]+muDir[1]*pi_dir[1]+muDir[2]*pi_dir[2];
    W = e->W();
    Q2 = e->q2();
    if (nMu==1&&nPi==1) {
      pimom = Ppi; picostheta = costheat_pi; mumom = Pmu; mucostheta = costheat_mu;
      if (e->par.nucleus_p==1) target=1;
      else target=6;
      out_seltree->Fill();
      if (channel==2&&pi_q==1&&target==6){
        rewhist1->Fill(Enu,-Q2/1.e6);
        rewhist2->Fill(Enu,target_mom);
        //rewhist3->Fill(Enu,-Q2/1.e6,mumom);
      } 
      int fluxbin = fluxhist->GetXaxis()->FindBin(Enu/1000.)-1;
      int q2bin = rewhist11->GetYaxis()->FindBin(-Q2/1.e6)-1;
      if (fluxbin>=nfluxbins) fluxbin = nfluxbins-1;
      if (q2bin>=nQ2bins) q2bin = nQ2bins-1;
      if (channel==2&&pn==1&&pi_q==1&&target==6){
        rewhist11->Fill(Enu,-Q2/1.e6);
        rewhist3_11->Fill(Enu,-Q2/1.e6,W);
        rewhist31_11->Fill(Enu,-Q2/1.e6,q0);
        rewhist3[fluxbin][0]->Fill(-Q2/1.e6,q0,W);
        rewhist5[fluxbin][q2bin][0]->Fill(q0,W,Ehad/W);
        rewhist3_ppi[fluxbin][0]->Fill(-Q2/1.e6,q0,pimom);
        rewhistEnuW11->Fill(Enu,W);
      }
      if (channel==2&&pn==0&&pi_q==0&&target==6){
        rewhist12->Fill(Enu,-Q2/1.e6);
        rewhist3_12->Fill(Enu,-Q2/1.e6,W);
        rewhist31_12->Fill(Enu,-Q2/1.e6,q0);
        rewhist3[fluxbin][1]->Fill(-Q2/1.e6,q0,W);
        rewhist5[fluxbin][q2bin][1]->Fill(q0,W,Ehad/W);
        rewhist3_ppi[fluxbin][1]->Fill(-Q2/1.e6,q0,pimom);
        rewhistEnuW12->Fill(Enu,W);
      }
      if (channel==2&&pn==0&&pi_q==1&&target==6){
        rewhist13->Fill(Enu,-Q2/1.e6);
        rewhist3_13->Fill(Enu,-Q2/1.e6,W);
        rewhist31_13->Fill(Enu,-Q2/1.e6,q0);
        rewhist3[fluxbin][2]->Fill(-Q2/1.e6,q0,W);
        rewhist5[fluxbin][q2bin][2]->Fill(q0,W,Ehad/W);
        rewhist3_ppi[fluxbin][2]->Fill(-Q2/1.e6,q0,pimom);
        rewhistEnuW13->Fill(Enu,W);
      }
    }
  }
  rewhist1->Scale(xsections->Integral()/t->GetEntries());
  rewhist2->Scale(xsections->Integral()/t->GetEntries());
  //rewhist3->Scale(xsections->Integral()/t->GetEntries());
  rewhist3_11->Scale(xsections->Integral()/t->GetEntries());
  rewhist3_12->Scale(xsections->Integral()/t->GetEntries());
  rewhist3_13->Scale(xsections->Integral()/t->GetEntries());
  rewhist31_11->Scale(xsections->Integral()/t->GetEntries());
  rewhist31_12->Scale(xsections->Integral()/t->GetEntries());
  rewhist31_13->Scale(xsections->Integral()/t->GetEntries());
  rewhist11->Scale(xsections->Integral()/t->GetEntries());
  rewhist12->Scale(xsections->Integral()/t->GetEntries());
  rewhist13->Scale(xsections->Integral()/t->GetEntries());
  rewhistEnuW11->Scale(xsections->Integral()/t->GetEntries());
  rewhistEnuW12->Scale(xsections->Integral()/t->GetEntries());
  rewhistEnuW13->Scale(xsections->Integral()/t->GetEntries());
  out_file->cd();
  out_seltree->Write();
  xsections->Write();
  rewhist1->Write();
  rewhist2->Write();
  //rewhist3->Write();
  rewhist3_11->Write();
  rewhist3_12->Write();
  rewhist3_13->Write();
  rewhist31_11->Write();
  rewhist31_12->Write();
  rewhist31_13->Write();
  rewhist11->Write();
  rewhist12->Write();
  rewhist13->Write();
  rewhistEnuW11->Write();
  rewhistEnuW12->Write();
  rewhistEnuW13->Write();
  for (int i=0;i<nfluxbins;i++){
    rewhist3[i][0]->Scale(xsections->Integral()/t->GetEntries());
    rewhist3[i][1]->Scale(xsections->Integral()/t->GetEntries());
    rewhist3[i][2]->Scale(xsections->Integral()/t->GetEntries());
    rewhist3[i][0]->Write();
    rewhist3[i][1]->Write();
    rewhist3[i][2]->Write();

    rewhist3_ppi[i][0]->Scale(xsections->Integral()/t->GetEntries());
    rewhist3_ppi[i][1]->Scale(xsections->Integral()/t->GetEntries());
    rewhist3_ppi[i][2]->Scale(xsections->Integral()/t->GetEntries());
    rewhist3_ppi[i][0]->Write();
    rewhist3_ppi[i][1]->Write();
    rewhist3_ppi[i][2]->Write();

    for (int j=0;j<nQ2bins;j++) {
      rewhist5[i][j][0]->Scale(xsections->Integral()/t->GetEntries());
      rewhist5[i][j][1]->Scale(xsections->Integral()/t->GetEntries());
      rewhist5[i][j][2]->Scale(xsections->Integral()/t->GetEntries());
      rewhist5[i][j][0]->Write();
      rewhist5[i][j][1]->Write();
      rewhist5[i][j][2]->Write();
    }
  }
  //hdptt->Draw();
  out_file->Close();
}

void analysis_RES(){

  gStyle->SetOptStat(0);

  //TH1D* hGraph_LFG[3];TH1D* hGraph_RFG[3];TH1D* hGraph_BRRFG[3];TH1D* hGraph_SF[3];
  //analysis_getGraph("out_LFG_C.root");
  //analysis_getGraph("out_RFG_C.root");
  //analysis_getGraph("out_SF_C.root");
  //analysis_getGraph("out_BRRFG_C.root");
  //analysis_getGraph("out_RFG_C_LFGEb.root");
  //analysis_getGraph("out_RFG_C_Eb27.root");
  //analysis_getGraph("out_proton2.root");
  //analysis_getGraph("out_neutron2.root");
  analysis_getGraph("out_RFG_O_numub_noPDD_Eb50.root");
  //analysis_getGraph("out_RFG_C_noPDD_Eb27.root");
  //analysis_getGraph("test.root");
  gApplication->Terminate();

}
