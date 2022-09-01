#include "SFGAnalysisProject_cxx/SFGAnalysisProject_cxxProjectHeaders.h"

// gSystem->Load("SFGAnalysisProject_cxx/SFGAnalysisProject.so");

// typedef TClonesArray TSFGAlgoResContainer;
// using namespace ND::TSFGReconModule;

std::vector<int> FindParentChild(std::vector<TLorentzVector> startingPos, std::vector<TLorentzVector> stoppingPos, std::vector<int> parents)
{
    // for (int i=0;i<parents.size();i++)
    //     std::cout<<parents[i]<<"\t";
    // std::cout<<std::endl;

    bool allFilled = true;
    for (int i=0;i<parents.size();i++)
    {
        if (parents[i]==-1) 
        {
            allFilled = false;
            break;
        }
    }
    if (allFilled) return parents;

    bool update = false;
    std::vector<double> distances(parents.size(),1.e6);
    std::vector<int> parentIdx(parents.size(),-1);
    for (int i=0;i<startingPos.size();i++)
    {
        if (parents[i]!=-1) continue;
        for (int j=0;j<startingPos.size();j++)
        {
            if (parents[j]==-1 || i==j ) continue;
            double dist1 = (startingPos[i]-startingPos[j]).P();
            double dist2 = (startingPos[i]-stoppingPos[j]).P();
            double dist3 = (stoppingPos[i]-startingPos[j]).P();
            double dist4 = (stoppingPos[i]-stoppingPos[j]).P();
            double dist = std::min( std::min(dist1,dist2) , std::min(dist3,dist4) );
            if (dist<distances[i])
            {
                parentIdx[i] = j;
                distances[i] = dist;
            }
        }
    }
    for (int i=0;i<startingPos.size();i++)
    {
        if (distances[i]<100)
        {
            parents[i] = parentIdx[i];
            update = true;
        }
    }
    if (!update)
    {
        double minDist = 1.e6;
        double minIdx = -1;
        for (int i=0;i<parents.size();i++)
        {
            if (minDist>distances[i])
            {
                minDist = distances[i];
                minIdx = i;
            }
        }
        parents[minIdx] = parentIdx[minIdx];
    }

    return FindParentChild(startingPos, stoppingPos, parents);
}

int FindParent(std::vector<TLorentzVector> stoppingPos, TLorentzVector startPos, TVector3 startDir, int TrkId)
{
    int idx = -1;
    double minDist = 1.e6;
    //std::cout<<"Track "<<TrkId<<std::endl;
    for (int i=0;i<stoppingPos.size();i++)
    {
        if (i==TrkId) continue;
        TVector3 relPos = (startPos-stoppingPos[i]).Vect();
        TVector3 relDir = relPos.Unit();
        //if (relDir.Dot(startDir)<0) continue;
        double dist = relPos.Mag();
        if (dist<50)
        {
            if ( dist<minDist )
            {
                idx = i;
                minDist = dist;
            }
        }
        else if (relDir.Dot(startDir)>0)
        {
            if ( dist<minDist )
            {
                idx = i;
                minDist = dist;
            }
        }

        //std::cout<<"Parent "<<i<<": "<<relDir.Dot(startDir)<<" "<<dist<<" "<<idx<<std::endl;
    }

    return idx;
}

void selection(std::string filename)
{
    std::vector<std::string> pidname = {"NotSet","Other","Shower","EM","Electron","Gamma","Hadronic","PiZero","Light","Muon","Pion","Heavy","Proton","Kaon"};

    TFile* f = new TFile(filename.c_str());
    //if (gSystem->Load("SFGAnalysisProject_cxx/SFGAnalysisProject.so")<0)
    //    f->MakeProject("SFGAnalysisProject_cxx","ND::TSFGReconModule","recreate++");
    //gSystem->Load("SFGAnalysisProject_cxx/SFGAnalysisProject_cxx.so");

    double ProtonBraggCutLow = 850.0;
    double PionBraggCutLow = 430.0;
    double PionBraggCutHigh = 1000.0;
    int VertexCut = 2;

    TTree* sfgTree = (TTree*)f->Get("ReconDir/SFG");
    TTree* trajTree = (TTree*)f->Get("TruthDir/Trajectories");

    int NAlgoResults, NVertices, NParticles, NTracks, NShowers, NClusters, NNodes, NHits, NFibers, NTrueHits, EventID;
    TClonesArray* fAlgoResults = new TClonesArray("ND::TSFGReconModule::TSFGAlgoRes");
    TClonesArray* fParticles = new TClonesArray("ND::TSFGReconModule::TSFGParticle");
    TClonesArray* fNodes = new TClonesArray("ND::TSFGReconModule::TSFGNode");
    TClonesArray* fTracks = new TClonesArray("ND::TSFGReconModule::TSFGTrack");
    TClonesArray* fClusters = new TClonesArray("ND::TSFGReconModule::TSFGCluster");
    sfgTree->SetBranchAddress("NAlgoResults",&NAlgoResults);
    sfgTree->SetBranchAddress("AlgoResults",&fAlgoResults);
    sfgTree->SetBranchAddress("NParticles",&NParticles);
    sfgTree->SetBranchAddress("Particles",&fParticles);
    sfgTree->SetBranchAddress("NNodes",&NNodes);
    sfgTree->SetBranchAddress("Nodes",&fNodes);
    sfgTree->SetBranchAddress("NTracks",&NTracks);
    sfgTree->SetBranchAddress("Tracks",&fTracks);
    sfgTree->SetBranchAddress("NClusters",&NClusters);
    sfgTree->SetBranchAddress("Clusters",&fClusters);
    sfgTree->SetBranchAddress("EventID",&EventID);

    int NTraj;
    TClonesArray* fTrajectories = new TClonesArray("ND::TTruthTrajectoriesModule::TTruthTrajectory");
    trajTree->SetBranchAddress("NTraj",&NTraj);
    trajTree->SetBranchAddress("Trajectories",&fTrajectories);

    std::string outname = filename;
    outname.erase(outname.find("anal_000_bsdv01_2.root"));
    outname+="sele_000.root";
    //std::cout<<"outname = "<<outname<<std::endl;

    TFile* fout = new TFile(outname.c_str(),"RECREATE");
    int eventId, trackId, trackPID, truthId, truthPID, trackContain, ProtonBragg, PionBragg;
    double trackQuality, truthMomentum, truthCharge, BraggParameter, trackMomentum, trackLength, stopX, stopY, stopZ, ddedx;
    double lightPID, heavyPID, emPID;
    double MomentumByRangePion, MomentumDeDxPion;
    TTree* sel = new TTree("sel","sel");
    sel->Branch("eventId",&eventId);
    sel->Branch("trackId",&trackId);
    sel->Branch("trackPID",&trackPID);
    sel->Branch("trackQuality",&trackQuality);
    sel->Branch("trackContain",&trackContain);
    sel->Branch("trackMomentum",&trackMomentum);
    sel->Branch("trackLength",&trackLength);
    sel->Branch("stopX",&stopX);
    sel->Branch("stopY",&stopY);
    sel->Branch("stopZ",&stopZ);
    sel->Branch("lightPID",&lightPID);
    sel->Branch("heavyPID",&heavyPID);
    sel->Branch("emPID",&emPID);
    sel->Branch("BraggParameter",&BraggParameter);
    sel->Branch("ProtonBragg",&ProtonBragg);
    sel->Branch("PionBragg",&PionBragg);
    sel->Branch("MomentumByRangePion",&MomentumByRangePion);
    sel->Branch("MomentumDeDxPion",&MomentumDeDxPion);
    sel->Branch("truthId",&truthId);
    sel->Branch("truthPID",&truthPID);
    sel->Branch("truthMomentum",&truthMomentum);
    sel->Branch("truthCharge",&truthCharge);
    int maintrkId, maintrkPID, maintrkContain, maintrkProtonBragg, maintrkPionBragg, nSecondary, nCloseTrks;
    double maintrkQuality, maintrktruthMomentum, maintrktruthCharge, maintrkPIDWeight, maintrkPionWeight, dx, dth, maintrkBraggParameter, maintrkMomentum, maintrkLength;
    double maintrkstopX, maintrkstopY, maintrkstopZ, maintrkddedx, maintrkdedxmax;
    double maintrkLightPID, maintrkHeavyPID, maintrkEMPID;
    double mergeLength, truthStop;
    double maintrkMomentumByRangePion, maintrkMomentumDeDxPion;
    int vtxpid[100];
    double vtxdist[100], vtxtime[100];
    int nCloseClus;
    double clusterdist[1000];
    double clusEDep;
    int nRecoChild, nRecoEM, nRecoLight, nRecoHeavy, dcye;
    double totalEdep, RecoEdep, maintrkBB, mergeEdep, extraEdep;
    double maintrkEdep ,nodeCombineEdep, maintrkShowerEnergy;
    double truthSecondaryEnergy, transverseComponent;
    int mergeLight, mergeHeavy, mergeEM, mergeTruthPID, mergeContain;
    int nTrajPoints;
    int delta_nclust;
    double trajPointDist[100], trajPointBB[100];
    TTree* maintrk = new TTree("maintrk","maintrk");
    maintrk->Branch("eventId",&eventId);
    maintrk->Branch("maintrkId",&maintrkId);
    maintrk->Branch("maintrkPID",&maintrkPID);
    maintrk->Branch("maintrkQuality",&maintrkQuality);
    maintrk->Branch("maintrkPIDWeight",&maintrkPIDWeight);
    maintrk->Branch("maintrkPionWeight",&maintrkPionWeight);
    maintrk->Branch("maintrkContain",&maintrkContain);
    maintrk->Branch("maintrkMomentum",&maintrkMomentum);
    maintrk->Branch("maintrkLength",&maintrkLength);
    maintrk->Branch("maintrkBraggParameter",&maintrkBraggParameter);
    maintrk->Branch("maintrkProtonBragg",&maintrkProtonBragg);
    maintrk->Branch("maintrkPionBragg",&maintrkPionBragg);
    maintrk->Branch("maintrkstopX",&maintrkstopX);
    maintrk->Branch("maintrkstopY",&maintrkstopY);
    maintrk->Branch("maintrkstopZ",&maintrkstopZ);
    maintrk->Branch("maintrkLightPID",&maintrkLightPID);
    maintrk->Branch("maintrkHeavyPID",&maintrkHeavyPID);
    maintrk->Branch("maintrkEMPID",&maintrkEMPID);
    maintrk->Branch("MomentumByRangePion",&maintrkMomentumByRangePion);
    maintrk->Branch("MomentumDeDxPion",&maintrkMomentumDeDxPion);
    maintrk->Branch("truthMomentum",&maintrktruthMomentum);
    maintrk->Branch("truthSecondaryEnergy",&truthSecondaryEnergy);
    maintrk->Branch("truthCharge",&maintrktruthCharge);
    maintrk->Branch("dx",&dx);
    maintrk->Branch("dth",&dth);
    maintrk->Branch("nRecoChild",&nRecoChild);
    maintrk->Branch("nRecoEM",&nRecoEM);
    maintrk->Branch("nRecoLight",&nRecoLight);
    maintrk->Branch("nRecoHeavy",&nRecoHeavy);
    maintrk->Branch("totalEdep",&totalEdep);
    maintrk->Branch("maintrkShowerEnergy",&maintrkShowerEnergy);
    maintrk->Branch("maintrkEdep",&maintrkEdep);
    maintrk->Branch("nodeCombineEdep",&nodeCombineEdep);
    maintrk->Branch("RecoEdep",&RecoEdep);
    maintrk->Branch("extraEdep",&extraEdep);
    maintrk->Branch("nSecondary",&nSecondary);
    maintrk->Branch("dcye",&dcye);
    maintrk->Branch("maintrkddedx",&maintrkddedx);
    maintrk->Branch("maintrkdedxmax",&maintrkdedxmax);
    maintrk->Branch("maintrkBB",&maintrkBB);
    maintrk->Branch("truthStop",&truthStop);
    maintrk->Branch("nCloseTrks",&nCloseTrks);
    maintrk->Branch("mergeEdep",&mergeEdep);
    maintrk->Branch("mergeLight",&mergeLight);
    maintrk->Branch("mergeHeavy",&mergeHeavy);
    maintrk->Branch("mergeEM",&mergeEM);
    maintrk->Branch("mergeContain",&mergeContain);
    maintrk->Branch("mergeLength",&mergeLength);
    maintrk->Branch("mergeTruthPID",&mergeTruthPID);
    maintrk->Branch("NClusters",&NClusters);
    maintrk->Branch("nCloseClus",&nCloseClus);
    maintrk->Branch("delta_nclust",&delta_nclust);
    maintrk->Branch("clusEDep",&clusEDep);
    maintrk->Branch("nTrajPoints",&nTrajPoints);
    maintrk->Branch("transverseComponent",&transverseComponent);
    maintrk->Branch("vtxdist",vtxdist,"vtxdist[nSecondary]/D");
    maintrk->Branch("vtxtime",vtxtime,"vtxtime[nSecondary]/D");
    maintrk->Branch("vtxpid",vtxpid,"vtxpid[nSecondary]/I");
    maintrk->Branch("clusterdist",clusterdist,"clusterdist[NClusters]/D");
    maintrk->Branch("trajPointDist",trajPointDist,"trajPointDist[nTrajPoints]/D");
    maintrk->Branch("trajPointBB",trajPointBB,"trajPointBB[nTrajPoints]/D");

    double delta, delta_mom, delta_pid, delta_clusterdist, delta_clustertime, delta_excess;
    TTree* deltatree = new TTree("deltatree","deltatree");
    deltatree->Branch("delta",&delta);
    deltatree->Branch("delta_mom",&delta_mom);
    deltatree->Branch("delta_pid",&delta_pid);
    deltatree->Branch("delta_clusterdist",&delta_clusterdist);
    deltatree->Branch("delta_clustertime",&delta_clustertime);
    deltatree->Branch("delta_excess",&delta_excess);

    TTree* neutroncluster = new TTree("neutroncluster","neutroncluster");
    double neutrondist, neutronmom, neutrontime, neutronedp, neutrontimediff;
    neutroncluster->Branch("neutrondist",&neutrondist);
    neutroncluster->Branch("neutronmom",&neutronmom);
    neutroncluster->Branch("neutrontime",&neutrontime);
    neutroncluster->Branch("neutrontimediff",&neutrontimediff);
    neutroncluster->Branch("neutronedp",&neutronedp);

    TTree* deltacluster = new TTree("deltacluster","deltacluster");
    double deltacluster_dist, deltacluster_time, deltacluster_edp, deltacluster_pid, deltacluster_mom;
    deltacluster->Branch("deltacluster_dist",&deltacluster_dist);
    deltacluster->Branch("deltacluster_time",&deltacluster_time);
    deltacluster->Branch("deltacluster_edp",&deltacluster_edp);
    deltacluster->Branch("deltacluster_pid",&deltacluster_pid);
    deltacluster->Branch("deltacluster_mom",&deltacluster_mom);

    TTree* showercone = new TTree("showercone","showercone");
    double shower_costh, shower_energy;
    int shower_trkclus;
    showercone->Branch("shower_costh",&shower_costh);
    showercone->Branch("shower_energy",&shower_energy);
    showercone->Branch("shower_trkclus",&shower_trkclus);

    TH1D* hist_dedx = new TH1D("","",100,0,2000);
    TH1D* hist_trkBB = new TH1D("","",100,0,100);

    //TCanvas* c1 = new TCanvas();
    for (int i=0;i<sfgTree->GetEntries();i++)
    {
        fAlgoResults->Clear();
        sfgTree->GetEntry(i);

        fTrajectories->Clear();
        trajTree->GetEntry(i);

        eventId = EventID;
        maintrkId = -1;
        maintrkPID = -1;
        maintrkQuality = 0;
        maintrktruthMomentum = 0;
        maintrktruthCharge = 0;
        dx = 0;
        dth = 0;

        float mainTrkCharge = 0;
        double mainTrkStartZ = 0;
        double mainTrkChargeTmp=0;
        double mainTrkIdTmp;

        if (NAlgoResults==0 || NParticles==0) 
        {
            ND::TTruthTrajectoriesModule::TTruthTrajectory* trj = (ND::TTruthTrajectoriesModule::TTruthTrajectory*)fTrajectories->At(0);
            maintrktruthMomentum = trj->InitMomentum.P();
            maintrk->Fill();
            continue;
        }
        //std::cout<<"NAlgoResults = "<<NAlgoResults<<", NParticles = "<<NParticles<<std::endl;
        ND::TSFGReconModule::TSFGAlgoRes* result = (ND::TSFGReconModule::TSFGAlgoRes*)fAlgoResults->At(0);
        std::vector<int> Particles = result->Particles;
        //std::cout<<"\n\n\nIn result "<<i<<", NParticles = "<<Particles.size()<<std::endl;
        std::vector<TLorentzVector> startingPos;
        std::vector<TLorentzVector> stoppingPos;
        std::vector<TVector3> startingDir;
        std::vector<double> trklength;
        std::vector<double> trkEdep;
        TLorentzVector maintrkStartPos, maintrkStopPos;
        std::vector<int> trkPID;
        std::vector<int> trkTruthPID;
        std::vector<int> trkContain;
        std::vector<double> trkddedx(Particles.size(),0);
        //TH1D* hist[Particles.size()];
        for (int j=0;j<Particles.size();j++)
        {
            trackId = j;
            trackContain = 0;
            trackPID = -1;
            lightPID = 0;
            heavyPID = 0;
            emPID = 0;
            float maxweight = 0;
            float pionweight = 0;
            ND::TSFGReconModule::TSFGParticle* par = (ND::TSFGReconModule::TSFGParticle*)fParticles->At(Particles[j]);
            //std::cout<<"\nParticle "<<j<<std::endl;
            //std::cout<<"Position = " << par->Position[0] << " " << par->Position[1] << " " << par->Position[2] << " " << par->Position[3] << std::endl;
            //std::cout<<"PosVariance = " << par->PosVariance[0] << " " << par->PosVariance[1] << " " << par->PosVariance[2] << " " << par->PosVariance[3] << std::endl;
            //std::cout<<"Direction = " << par->Direction[0] << " " << par->Direction[1] << " " << par->Direction[2] << std::endl;
            //std::cout<<"DirVariance = " << par->DirVariance[0] << " " << par->DirVariance[1] << " " << par->DirVariance[2] << std::endl;
            //std::cout<<"Momentum = " <<par->Momentum<<" Charge = " <<par->Charge<<" Quality = " <<par->Quality <<std::endl;
            trackMomentum = par->Momentum;
            TLorentzVector Position = par->Position;
            startingPos.push_back(Position);
            TLorentzVector PosVariance = par->PosVariance;
            TVector3 Direction = par->Direction;
            TVector3 DirVariance = par->DirVariance;
            startingDir.push_back(Direction);
            trackQuality = par->Quality;
            MomentumByRangePion = par->MomentumByRangePion;
            MomentumDeDxPion = par->MomentumDeDxPion;
            vector<int>    PID = par->PID;              //
            vector<float>  PID_weight = par->PID_weight;
            //std::cout<<"PID\t weight\n";
            for (int k=0;k<PID.size();k++)
            {
                //std::cout<<pidname[PID[k]]<<"\t" <<PID_weight[k]<< std::endl;
                if (PID[k]==8||PID[k]==9||PID[k]==10)
                {
                    if (PID[k]==10)
                        trackContain = 1;
                    if (pionweight<PID_weight[k])
                        pionweight = PID_weight[k];
                    if (lightPID<PID_weight[k])
                        lightPID=PID_weight[k];
                }
                if (PID[k]==11||PID[k]==12)
                {
                    if (heavyPID<PID_weight[k])
                        heavyPID=PID_weight[k];
                }
                if (PID[k]==3)
                {
                    if (emPID<PID_weight[k])
                        emPID=PID_weight[k];
                }
                if (maxweight<PID_weight[k])
                {
                    trackPID = PID[k];
                    maxweight = PID_weight[k];
                }
            }
            trkPID.push_back(trackPID);
            trkContain.push_back(trackContain);

            std::vector<int> Nodes = par->Nodes;
            std::vector<float> Nodededx;
            int count = Nodes.size();
            int avgW = 5;
            std::vector<float> NodededxAvg(count,0);
            //hist[j] = new TH1D("","",count,0,count);
            for (int k=0;k<count;k++)
            {
                ND::TSFGReconModule::TSFGNode* node = (ND::TSFGReconModule::TSFGNode*)fNodes->At(Nodes[k]);
                Nodededx.push_back(node->EDeposit);
                for (int l=k; l<std::min(k+avgW,count);l++)
                    NodededxAvg[l] += node->EDeposit/avgW;
                if (k>=avgW&&k<count-avgW)
                {
                    double localvar = fabs((NodededxAvg[k]-NodededxAvg[k-1])/NodededxAvg[k-1]);
                    if (trkddedx[j]<localvar)
                    {
                        trkddedx[j] = localvar;
                        ddedx = localvar;
                    }
                }
                //hist[j]->SetBinContent(k,node->EDeposit);
                if (k==count-1)
                {
                    stoppingPos.push_back(node->Position);
                    stopX = node->Position.X();
                    stopY = node->Position.Y();
                    stopZ = node->Position.Z();
                }
            }
            int node_num = 5;
            if(count < node_num) node_num = count;
            std::vector<float> node_dedx;
            for(int k = 0; k < node_num; k++){
                if((count-k-1) < VertexCut) continue;
                node_dedx.push_back(Nodededx[count-k-1]);
            }
            float dedx_sum = 0;
            for(int k = 0; k < node_dedx.size(); k++) dedx_sum += node_dedx[k];
            if(node_dedx.size() > 0) dedx_sum /= node_dedx.size();
            BraggParameter = dedx_sum;
            if (BraggParameter>=ProtonBraggCutLow) ProtonBragg=1;
            else ProtonBragg=0;
            if (BraggParameter>=PionBraggCutLow&&BraggParameter<=PionBraggCutHigh) PionBragg=1;
            else PionBragg=0;

            std::vector<int> Tracks = par->Tracks;
            ND::TSFGReconModule::TSFGTrack* trk = (ND::TSFGReconModule::TSFGTrack*)fTracks->At(Tracks[0]);
            trackLength = trk->Length;
            trklength.push_back(trackLength);
            trkEdep.push_back(trk->EDeposit);

            std::vector<int> Truth_TrajIds = par->Truth_TrajIds;
            std::vector<float> Truth_ChargeShare = par->Truth_ChargeShare;
            truthId = -1;
            truthPID = 0;
            float maxcharge = 0;
            float totalcharge = 0;
            truthMomentum = 0;
            truthCharge = 0;
            bool isMainTrk = false;
            //std::cout<<"Truth particles\n";
            for (int k=0;k<Truth_TrajIds.size();k++)
            {
                totalcharge += Truth_ChargeShare[k];
                for (int l=0;l<NTraj;l++)
                {
                    ND::TTruthTrajectoriesModule::TTruthTrajectory* trj = (ND::TTruthTrajectoriesModule::TTruthTrajectory*)fTrajectories->At(l);
                    if (trj->ID==Truth_TrajIds[k])
                    {
                        //std::cout<<trj->ID<<"\t"<<trj->PDG<<"\t"<<trj->InitMomentum.P()<<"\t"<<Truth_ChargeShare[k] <<std::endl;
                        if (maxcharge<Truth_ChargeShare[k])
                        {
                            truthId = trj->ID;
                            truthPID = trj->PDG;
                            maxcharge = Truth_ChargeShare[k];
                            truthMomentum = trj->InitMomentum.P();
                        }

                        if (trj->ID==1)
                        {
                            TLorentzVector InitPosition = trj->InitPosition;
                            TVector3 InitDirection = trj->InitMomentum.Vect().Unit();
                            if 
                            ( 
                                ( (Position-InitPosition).P()<200 && Direction.Dot(InitDirection)>0.1 && Position.Z()<mainTrkStartZ ) ||
                                ( (stoppingPos[j]-InitPosition).P()<200 && -Direction.Dot(InitDirection)>0.1 && stoppingPos[j].Z()<mainTrkStartZ ) 
                                // fabs((Position[0]-InitPosition[0])/PosVariance[0])<1. &&
                                // fabs((Position[1]-InitPosition[1])/PosVariance[1])<1. &&
                                // fabs((Position[2]-InitPosition[2])/PosVariance[2])<1. &&
                                // fabs((Direction[0]-InitDirection[0])/DirVariance[0])<1. &&
                                // fabs((Direction[1]-InitDirection[1])/DirVariance[1])<1. &&
                                // fabs((Direction[2]-InitDirection[2])/DirVariance[2])<1. &&
                                //mainTrkCharge<Truth_ChargeShare[k]
                            )
                            {
                                dx = (Position-InitPosition).P();
                                dth = Direction.Dot(InitDirection);
                                maintrkId = trackId;
                                maintrkPID = trackPID;
                                maintrkQuality = trackQuality;
                                mainTrkCharge = Truth_ChargeShare[k];
                                maintrktruthMomentum = trj->InitMomentum.P();
                                maintrkPIDWeight = maxweight;
                                maintrkPionWeight = pionweight;
                                maintrkContain = trackContain;
                                maintrkBraggParameter = BraggParameter;
                                maintrkProtonBragg = ProtonBragg;
                                maintrkPionBragg = PionBragg;
                                maintrkMomentum = trackMomentum;
                                maintrkLength = trackLength;
                                maintrkstopX = stopX;
                                maintrkstopY = stopY;
                                maintrkstopZ = stopZ;
                                mainTrkStartZ = std::min(Position.Z(),stoppingPos[j].Z());
                                maintrkStopPos = Position.Z()<stoppingPos[j].Z() ? stoppingPos[j] : Position;
                                maintrkLightPID = lightPID;
                                maintrkHeavyPID = heavyPID;
                                maintrkEMPID = emPID;
                                maintrkddedx = ddedx;
                                maintrkMomentumByRangePion = MomentumByRangePion;
                                maintrkMomentumDeDxPion = MomentumDeDxPion;
                                isMainTrk = true;
                            }
                        }

                        break;
                    }
                }
            }
            trkTruthPID.push_back(truthPID);

            truthCharge = maxcharge/totalcharge;

            if (isMainTrk) maintrktruthCharge = mainTrkCharge/totalcharge;

            sel->Fill();
        }

        mergeLight=0; mergeHeavy=0; mergeEM=0;
        mergeLength = maintrkLength;
        mergeEdep = trkEdep[maintrkId];
        extraEdep = 0;
        nSecondary = Particles.size();
        nCloseTrks = 0;
        nRecoChild = 0; nRecoLight = 0; nRecoEM = 0; nRecoHeavy = 0;
        dcye = 0;
        for (int j=0;j<nSecondary;j++)
        {
            if (maintrkId==j || maintrkId<0) vtxdist[j]=-1;
            else 
            {
                vtxdist[j] = std::min( (maintrkStopPos-startingPos[j]).P() , (maintrkStopPos-stoppingPos[j]).P() );
                vtxtime[j] = std::min( startingPos[j].T() - maintrkStopPos.T() , stoppingPos[j].T() - maintrkStopPos.T() );
                vtxpid[j] = trkTruthPID[j];
                //if (vtxtime[j]>30) dcye = 1;
                if (vtxdist[j]<20) 
                {
                    nCloseTrks++;
                    if (trkPID[j]>=8&&trkPID[j]<=10) mergeLight = 1;
                    else if (trkPID[j]==3) mergeEM = 1;
                    else mergeHeavy = 1;

                    extraEdep += trkEdep[j];
                    mergeEdep += trkEdep[j];

                    mergeLength += trklength[j];
                    mergeTruthPID = trkTruthPID[j];
                    mergeContain = trkContain[j];

                    if (vtxtime[j]>30) dcye = 1;

                    // nRecoChild++;
                    // if (trkPID[j]>=8&&trkPID[j]<=10) nRecoLight++;
                    // else if (trkPID[j]==3) nRecoEM++;
                    // else nRecoHeavy++;
                } 
                // else if ( maintrkId == FindParent(stoppingPos, startingPos[j], startingDir[j], j) )
                // {
                    // nRecoChild++;
                    // if (trkPID[j]>=8&&trkPID[j]<=10) nRecoLight++;
                    // else if (trkPID[j]==3) nRecoEM++;
                    // else nRecoHeavy++;
                // }
            }
        }
        std::vector<int> Clusters=result->Clusters;
        NClusters = Clusters.size();
        nCloseClus = 0;
        clusEDep = 0;
        totalEdep = 0;
        RecoEdep = 0;
        std::vector<TLorentzVector> clusterPosition;
        std::vector<double> clusterEDeposit;
        for (int j=0;j<NClusters;j++)
        {
            ND::TSFGReconModule::TSFGCluster* clu = (ND::TSFGReconModule::TSFGCluster*)fClusters->At(Clusters[j]);
            TLorentzVector cluPosition = clu->Position;
            if (maintrkId<0) clusterdist[j]=-1;
            else 
            {
                clusterdist[j] = (maintrkStopPos-cluPosition).P();
                if (clusterdist[j]<200) 
                {
                    nCloseClus++;
                    clusEDep += clu->EDeposit;
                }
                //std::cout<<"clusterdist["<<j<<"]="<<clusterdist[j]<<std::endl;
            }
            clusterPosition.push_back(clu->Position);
            clusterEDeposit.push_back(clu->EDeposit);
        }

        maintrkBB = -1;
        maintrkdedxmax=0;
        ND::TTruthTrajectoriesModule::TTruthTrajectory* ptrj = (ND::TTruthTrajectoriesModule::TTruthTrajectory*)fTrajectories->At(0);
        TLorentzVector FinalPosition = ptrj->FinalPosition;
        truthStop = -1;
        if (maintrkId>=0)
        {
            maintrkShowerEnergy = 0;
            transverseComponent = 0;
            maintrkEdep = trkEdep[maintrkId];
            // hist[maintrkId]->Draw();
            // if (eventId/100000==1000)
            //     c1->SaveAs(Form("SFGevent3D/event_%i_dedx.pdf",eventId));
            //std::cout<<"Event "<<eventId<<std::endl;
            std::vector<int> parents(nSecondary,-1);
            parents[maintrkId] = maintrkId;
            std::vector<int> parentsNew = FindParentChild(startingPos, stoppingPos, parents);
            for (int j=0;j<nSecondary;j++)
            {
                totalEdep += trkEdep[j];
                if (j==maintrkId) 
                {
                    RecoEdep += trkEdep[j];
                    continue;
                }
                else
                {
                    shower_trkclus = 1;
                    shower_costh = (startingPos[j]-startingPos[maintrkId]).Vect().Unit().Dot(startingDir[maintrkId]);
                    shower_energy = trkEdep[j];
                    showercone->Fill();

                    if (shower_costh>0.6) maintrkShowerEnergy+=shower_energy;
                    else if ((startingPos[j]-stoppingPos[maintrkId]).P()<20) maintrkShowerEnergy+=shower_energy;
                }
                if ( maintrkId == FindParent(stoppingPos, startingPos[j], startingDir[j], j) )
                //if (parentsNew[j]==maintrkId)
                {
                    nRecoChild++;
                    RecoEdep += trkEdep[j];
                    if (trkPID[j]>=8&&trkPID[j]<=10) nRecoLight++;
                    else if (trkPID[j]==3) nRecoEM++;
                    else nRecoHeavy++;
                }
            }

            std::vector<TLorentzVector> childPos;
            std::vector<TLorentzVector> childEndPos;
            std::vector<TLorentzVector> childMom;
            std::vector<int> childPid;
            std::vector<double> childDelta;
            std::vector<double> childExcess;
            std::vector<ND::TTruthTrajectoriesModule::TTruthTrajectoryPoint> PrimaryPoints;
            truthSecondaryEnergy = 0;
            for (int l=0;l<NTraj;l++)
            {
                ND::TTruthTrajectoriesModule::TTruthTrajectory* trj = (ND::TTruthTrajectoriesModule::TTruthTrajectory*)fTrajectories->At(l);
                if (trj->ParentID==1)
                {
                    childPos.push_back(trj->InitPosition);
                    childEndPos.push_back(trj->FinalPosition);
                    childMom.push_back(trj->InitMomentum);
                    childPid.push_back(trj->PDG);
                    childDelta.push_back(0);
                    childExcess.push_back(0);

                    truthSecondaryEnergy += trj->InitMomentum.E()-trj->Mass;
                }
                // if (trj->ID==1)
                //     PrimaryPoints = trj->Points;
                if (trj->PDG==2112)
                {
                    for (int c=0;c<clusterPosition.size();c++)
                    {
                        if ((trj->FinalPosition-clusterPosition[c]).P()<50)
                        {
                            neutrondist = (trj->FinalPosition-trj->InitPosition).P();
                            neutrontime =  (trj->FinalPosition-trj->InitPosition).T();
                            neutrontimediff = fabs(trj->FinalPosition.T()-clusterPosition[c].T());
                            neutronmom = trj->InitMomentum.P();
                            neutronedp = clusterEDeposit[c];
                            neutroncluster->Fill();
                            break;
                        }
                    }
                }
            }
            ND::TTruthTrajectoriesModule::TTruthTrajectory* trjP = (ND::TTruthTrajectoriesModule::TTruthTrajectory*)fTrajectories->At(0);
            nTrajPoints = trjP->Points.size();
            std::vector<TLorentzVector> PrimaryPointsLV;
            for (int l=0;l<nTrajPoints;l++)
            {
                PrimaryPointsLV.push_back(TLorentzVector(trjP->Points[l].PositionX,trjP->Points[l].PositionY,trjP->Points[l].PositionZ,trjP->Points[l].PositionT));
                trajPointDist[l]=1.e6;
                trajPointBB[l]=0.;
            }

            ND::TSFGReconModule::TSFGParticle* par = (ND::TSFGReconModule::TSFGParticle*)fParticles->At(Particles[maintrkId]);
            std::vector<int> Nodes = par->Nodes;
            int count = Nodes.size();
            std::vector<double> dedx_sum(count,0);
            std::vector<double> dedx_max(childPos.size(),-1);
            std::vector<double> trkBB_max(childPos.size(),-1);
            std::vector<std::vector<double>> dedx_slide(count,std::vector<double>());
            std::vector<double> dedx_ma(count,0);
            std::vector<double> dedx_rms(count,0);
            std::vector<double> dedx_BB(count,0);
            std::vector<double> dedx_node;
            TLorentzVector deltapos;
            // TCanvas c1;
            // TH1D hist_dedx_local("","",count,0,count);
            // TH1D hist_bb_u("","",count,0,count);
            // TH1D hist_bb_l("","",count,0,count);
            double dedxsum = 0;
            int period = 5;
            int nsig = 1;
            nodeCombineEdep = trkEdep[maintrkId];
            std::vector<bool> useTrk(Particles.size(),false);
            useTrk[maintrkId] = true;
            std::vector<bool> useClus(clusterPosition.size(),false);
            TVector3 prevDir = ((ND::TSFGReconModule::TSFGNode*)fNodes->At(Nodes[0]))->Direction;
            for (int k=0;k<count;k++) 
            {
                ND::TSFGReconModule::TSFGNode* node = (ND::TSFGReconModule::TSFGNode*)fNodes->At(Nodes[k]);
                TLorentzVector nodePosition = node->Position;

                TVector3 curDir = node->Direction;
                double costh = prevDir.Dot(curDir);
                if (costh<0.9)
                {
                    transverseComponent += (1.-(k+0.)/count)*sqrt(1-costh*costh);
                }
                prevDir = curDir;

                for (int l = k; l < std::min(count,k+period); l++ )
                {
                    dedx_sum[l]+=node->EDeposit/period;
                    dedx_slide[l].push_back(node->EDeposit);
                }

                dedx_node.push_back(node->EDeposit);
                dedx_ma[k] = TMath::Mean(dedx_slide[k].begin(),dedx_slide[k].end());
                dedx_rms[k] = TMath::RMS(dedx_slide[k].begin(),dedx_slide[k].end());
                dedx_BB[k] = (node->EDeposit-dedx_ma[k-1])/dedx_rms[k-1];

                for (int l=1;l<nTrajPoints;l++)
                {
                    double val = (nodePosition-PrimaryPointsLV[l]).P();
                    if (trajPointDist[l]>val)
                        trajPointDist[l] = val;
                    if (val<100)
                    {
                        if (trajPointBB[l]<dedx_BB[k])
                            trajPointBB[l]=dedx_BB[k];
                        // if (k>=count-5)
                        //     trajPointBB[l] = 0;
                    }
                }

                if ((FinalPosition-nodePosition).P()<20)
                    truthStop = (k+0.)/count;
                if (k==count-1&&truthStop<0)
                    truthStop = (FinalPosition-nodePosition).P();

                if (k>=period)
                    // if (node->EDeposit>dedx_ma[k]+nsig*dedx_rms[k])
                    //     maintrkBB = (k+1.)/count;
                {
                    if (maintrkBB<dedx_BB[k])
                    {
                        maintrkBB = dedx_BB[k];
                        deltapos = nodePosition;
                    }
                }

                // hist_dedx_local.SetBinContent(k+1,node->EDeposit);
                // hist_bb_u.SetBinContent(k+2,dedx_ma[k]+nsig*dedx_rms[k]);
                // hist_bb_l.SetBinContent(k+2,dedx_ma[k]);

                dedxsum += node->EDeposit;
                if (maintrkdedxmax<node->EDeposit)
                    maintrkdedxmax=node->EDeposit;

                for (int l=0;l<childPos.size();l++)
                {
                    if ((nodePosition-childPos[l]).P()<50)
                    {
                        if (dedx_max[l]<dedx_sum[k])
                            dedx_max[l]=dedx_sum[k];
                        if (trkBB_max[l]<dedx_BB[k])
                            trkBB_max[l]=dedx_BB[k];

                        if (childDelta[l]<dedx_BB[k])
                            childDelta[l]=dedx_BB[k];

                        if (childExcess[l]<node->EDeposit-dedx_ma[k-1])
                            childExcess[l] = node->EDeposit-dedx_ma[k-1];
                    }
                }

                for (int c=0;c<clusterPosition.size();c++)
                {
                    if (!useClus[c] && (nodePosition-clusterPosition[c]).P()<50 )
                    {
                        nodeCombineEdep += clusterEDeposit[c];
                        useClus[c] = true;
                    }
                }

                for (int l=0;l<startingPos.size();l++)
                {
                    if (!useTrk[l])
                    {
                        if ((nodePosition-startingPos[l]).P()<20)
                        {
                            nodeCombineEdep +=  trkEdep[l];
                            useTrk[l] = true;
                        }
                        else if ((startingPos[l]-nodePosition).Vect().Unit().Dot(startingDir[l])>0.8)
                        {
                            nodeCombineEdep +=  trkEdep[l];
                            useTrk[l] = true;
                        }
                    }
                }
            }

            delta_nclust = 0;
            for (int c=0;c<clusterPosition.size();c++)
            {
                if ((clusterPosition[c]-deltapos).P()<1000&&clusterPosition[c].T()-deltapos.T()<10)
                {
                    delta_nclust++;
                }

                deltacluster_dist = (clusterPosition[c]-deltapos).P();
                deltacluster_time = clusterPosition[c].T()-deltapos.T();
                deltacluster_edp = clusterEDeposit[c];
                deltacluster_pid = 0;
                deltacluster_mom = -1;
                for (int l=0;l<childPos.size();l++)
                {
                    if ((clusterPosition[c]-childEndPos[l]).P()<50 && (childPos[l]-deltapos).P()<50)
                    {
                        deltacluster_pid = childPid[l];
                        deltacluster_mom = childMom[l].P();
                        break;
                    }
                }
                deltacluster->Fill();

                shower_trkclus = 0;
                shower_costh = (clusterPosition[c]-startingPos[maintrkId]).Vect().Unit().Dot(startingDir[maintrkId]);
                shower_energy = clusterEDeposit[c];
                showercone->Fill();

                if (shower_costh>0.6) maintrkShowerEnergy+=shower_energy;
            }

            for (int l=0;l<childPos.size();l++)
            {
                if (childDelta[l]>0)
                {
                    delta = childDelta[l];
                    delta_mom = childMom[l].P();
                    delta_pid = childPid[l];
                    delta_excess = childExcess[l];
                    delta_clusterdist = -1;
                    delta_clustertime = -1;
                    if (delta_pid==2112)
                    {
                        bool withCluster = false;
                        for (int c=0;c<clusterPosition.size();c++)
                        {
                            if ((clusterPosition[c]-childEndPos[l]).P()<50)
                            {
                                withCluster = true;
                                break;
                            }
                        }
                        if (withCluster)
                        {
                            delta_clusterdist = (childPos[l]-childEndPos[l]).P();
                            delta_clustertime = childEndPos[l].T()-childPos[l].T();
                        }
                    }
                    deltatree->Fill();
                }
            }

            // maintrkBB = TMath::RMS(dedx_node.begin(),dedx_node.end())/TMath::Mean(dedx_node.begin(),dedx_node.end());
            // hist_dedx_local.Draw();
            // hist_bb_u.SetLineColor(kRed);
            // hist_bb_u.Draw("same");
            // hist_bb_l.SetLineColor(kGreen);
            // hist_bb_l.Draw("same");
            // c1.SaveAs(Form("SFGevent3D/event_%i_dedx.pdf",eventId));
            if (count>5)
                maintrkdedxmax = maintrkdedxmax/(dedxsum/(count-5.));
            else 
                maintrkdedxmax = 0;
                
            for (int l=0;l<childPos.size();l++)
                if (dedx_max[l]>0)
                {
                    hist_dedx->Fill(dedx_max[l]);
                    hist_trkBB->Fill(trkBB_max[l]);
                }
            for (int l=1;l<nTrajPoints;l++)
                if (trajPointBB[l]>0)
                    hist_trkBB->Fill(trajPointBB[l]);

        }

        maintrk->Fill();

        // for (int j=0;j<NAlgoResults;j++)
        // {
        //     ND::TSFGReconModule::TSFGAlgoRes* result = (ND::TSFGReconModule::TSFGAlgoRes*)fAlgoResults->At(j);
        //     std::cout<<"Name = "<<result->AlgorithmName<<std::endl;
        // }
    }

    fout->cd();
    sel->Write();
    maintrk->Write();
    deltatree->Write();
    neutroncluster->Write();
    deltacluster->Write();
    showercone->Write();
    hist_dedx->Write("hist_dedx");
    hist_trkBB->Write("hist_trkBB");
    fout->Close();
    f->Clear();
}

void test()
{
    gSystem->Load("SFGAnalysisProject_cxx/SFGAnalysisProject_cxx.so");

    TSystemDirectory dir(gSystem->pwd(),gSystem->pwd());
    TList *files = dir.GetListOfFiles();
    for (int i=0;i<files->GetSize();i++)
    {
        std::string fname = (std::string)files->At(i)->GetName();
        if (fname.find("anal_000_bsdv01_2.root")!=std::string::npos)
        //if (fname.find("oa_pg_pip_00001234-1000_zthr6a35oiuz_anal_000_bsdv01_2.root")!=std::string::npos)
        {
            std::cout<<"Running file: "<<fname<<std::endl;
            selection(fname);
        }
    }
}