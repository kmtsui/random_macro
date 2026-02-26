#define HEMI_CUDA_DISABLE 1
#define NOSKLIBRARIES

R__LOAD_LIBRARY($WCSIM_BUILD_DIR/lib/libWCSimRoot.so)
R__LOAD_LIBRARY(libfiTQun.so)

// run_fitqun_macro.C
//
// Example ROOT macro to demonstrate how to run fiTQun from a shared library.
//
// Usage in ROOT:
// root[] .L run_fitqun_macro.C
// root[] run_fitqun_macro("path/to/your/wcsim_file.root")

#include <iostream>
#include <vector>
#include <string>

// It's crucial to include the headers that define the classes 
// and common blocks you want to interact with.
#include "fiTQun.h"
#include "fiTQun_shared.h"
#include "spliTChan.h"
#include "WCSimWrap.h"

// spliTChanOutC.h defines the MAXSE macro used in fitqunoutC.h,
// so it must be included first.
// #include "spliTChanOutC.h" 
// #include "fitqunoutC.h" // To access the output common blocks

// Manually define the struct layout to match the library (from fitqunoutC.h)
// This avoids including the header which causes symbol redefinition issues in ROOT macros.
#define MAXNPEAK 10
#define MAXPX 7

extern "C" {
    struct fitqun1r_common {
        int    fqnse;
        int    fqitwnd[MAXNPEAK];
        int    fqipeak[MAXNPEAK];
        int    fqnhitpmt[MAXNPEAK];
        float  fqtotq[MAXNPEAK];
        float  fq0rtotmu[MAXNPEAK];
        float  fq0rnll[MAXNPEAK];
        int    fqn50[MAXNPEAK];
        float  fqq50[MAXNPEAK];
        int    fq1rpcflg[MAXNPEAK][MAXPX];
        float  fq1rmom[MAXNPEAK][MAXPX];
        float  fq1rt0[MAXNPEAK][MAXPX];
        float  fq1rtotmu[MAXNPEAK][MAXPX];
        float  fq1rnll[MAXNPEAK][MAXPX];
        float  fq1rpos[MAXNPEAK][MAXPX][3];
        float  fq1rdir[MAXNPEAK][MAXPX][3];
        float  fq1rdconv[MAXNPEAK][MAXPX];
        float  fq1reloss[MAXNPEAK][MAXPX];
    };
    extern struct fitqun1r_common fitqun1r_;
}

#include "TRuntimeParameters.hxx"

using namespace fiTQun_parameters;

std::vector<std::vector<double>> setup_eventDiplayXY(const char* inFileName)
{

    TFile* f = new TFile(inFileName, "read");
    if (!f || !f->IsOpen()) {
        std::cerr << "Error: Could not open file: " << inFileName << std::endl;
        exit(9);
    }

    // Geometry tree - only need 1 "event"
    WCSimRootGeom *geo = 0;
    TTree *geotree = (TTree*)f->Get("wcsimGeoT");
    geotree->SetBranchAddress("wcsimrootgeom", &geo);
    std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
    if (geotree->GetEntries() == 0) {
        exit(9);
    }
    geotree->GetEntry(0);
    int nPMTs_type0=geo->GetWCNumPMT();
    std::cout << "geo has " << nPMTs_type0 << " PMTs" << std::endl;
    std::vector<TVector3> pmt_posT(nPMTs_type0);
    double max_r = 0, max_z = 0;
    for (int i=0;i<nPMTs_type0;i++) 
    {
        WCSimRootPMT pmt;
        pmt = geo->GetPMT(i);
        std::vector<double> pos(3);
        for(int j=0;j<3;j++){
            pos[j] = pmt.GetPosition(j);
        }

        TVector3 pmtpos(pos[0],pos[1],pos[2]);
        pmt_posT[i] = pmtpos;
        
        // y-axis is vertical
        if (max_z<fabs(pos[1])) max_z=fabs(pos[1]);
        if (max_r<sqrt(pos[0]*pos[0]+pos[2]*pos[2]))
            if (fabs(pmt.GetOrientation(1))>0.5) max_r = sqrt(pos[0]*pos[0]+pos[2]*pos[2]);
    }

    double barrelCut = max_z-10;
    TH2D* hist_event_display = new TH2D("Charges","Charges",250,-TMath::Pi()*max_r,TMath::Pi()*max_r,250,-max_z-2*max_r,max_z+2*max_r);
    std::vector<std::vector<double>> eventDiplayXY;
    for (int i=0;i<nPMTs_type0;i++)
    {
        // rotation for event display
        double y = -pmt_posT.at(i).x();
        double x =  pmt_posT.at(i).z();
        double z =  pmt_posT.at(i).y();
        std::vector<double> pmtXY;
        if (fabs(z)<barrelCut) // barrel
        {
            double th = atan2(y,x);
            pmtXY.push_back(-max_r*th);
            pmtXY.push_back(z);
        }
        else if (z>barrelCut) //top
        {
            pmtXY.push_back(-y);
            pmtXY.push_back(max_z+max_r-x);
        }
        else //bot
        {
            pmtXY.push_back(-y);
            pmtXY.push_back(-max_z-max_r+x);
        }
        eventDiplayXY.push_back(pmtXY);
    }

    return eventDiplayXY;
}

void run_fitqun_macro(const char* inFileName) {
    // 1. Load the shared library
    //    This makes all the classes and functions from your library available in the ROOT session.
    std::cout << "Successfully loaded libfiTQun.so" << std::endl;

    // 2. Initialize Parameters
    //    In runfiTQun.cc, this is done by parsing command-line options. Here, we do it manually.
    //    This loads fiTQun.parameters.dat and spliTChan.parameters.dat.
    //    The TRuntimeParameters object is a singleton that holds all parameters.
    //    fiTQun_shared and other classes will fetch their parameters from it.
    TRuntimeParameters::Get().ReadParamOverrideFile("/eos/home-k/katsui/fiTQun/ParameterOverrideFiles/nuPRISMBeamTest_16cShort_mPMT.parameters.dat");
    // TRuntimeParameters::Get().ReadParamOverrideFile("spliTChan.parameters.dat");
    std::cout << "Successfully loaded parameter files." << std::endl;

    // 3. Initialize WCSimWrap
    //    This singleton wraps the WCSim ROOT file and provides access to event data.
    //    The constructor loads the file and reads the first event's geometry.
    WCSimWrap::Get((char*)inFileName);
    if (WCSimWrap::Get()->NEvt() == 0) {
        std::cerr << "Error: Could not open WCSim file or file contains no events: " << inFileName << std::endl;
        return;
    }
    std::cout << "Successfully initialized WCSimWrap with file: " << inFileName << std::endl;
    std::cout << "Number of events in file: " << WCSimWrap::Get()->NEvt() << std::endl;

    std::vector<std::vector<double>> eventDiplayXY = setup_eventDiplayXY(inFileName);

    // 4. Instantiate fiTQun
    //    The constructor requires the number of PMTs, which it gets from the WCSimWrap singleton.
    int writeHistos = TRuntimeParameters::Get().GetWriteHistos();
    fiTQun* thefit = new fiTQun(WCSimWrap::Get()->NPMT());
    spliTChan * theSubEvents = new spliTChan(writeHistos, WCSimWrap::Get()->NPMT());
    
    int fQuiet = TRuntimeParameters::Get().GetParameterI("fiTQun.fQuiet");
    fiTQun_shared::Get()->Cantyoubequiet(fQuiet);

    // Which fits should be run?
    
    bool flgSubEvt = TRuntimeParameters::Get().GetParameterI("fiTQun.DoSubEvent");
    
    bool flg1Rfit = TRuntimeParameters::Get().GetParameterI("fiTQun.DoSingleRingFits");
    
    int flgRingType[nPID];
    for (int i=0;i<nPID;i++) flgRingType[i] = 0;
    flgRingType[ie] = TRuntimeParameters::Get().GetParameterI("fiTQun.DoElectron1RFit");
    flgRingType[imu] = TRuntimeParameters::Get().GetParameterI("fiTQun.DoMuon1RFit");
    

    // std::cout << "flgSubEvt = " << flgSubEvt << std::endl;
    // std::cout << "flg1Rfit = " << flg1Rfit << std::endl;
    // std::cout << "flgMSFit = " << doMSFit << std::endl;
    
    // Options for all fits
    bool fFitCprof = TRuntimeParameters::Get().GetParameterI("fiTQun.UseFitCProfile");
    if (fFitCprof) std::cout << "Using fit Cherenkov profile!" << std::endl;
    bool fFittprof = TRuntimeParameters::Get().GetParameterI("fiTQun.UseFittProfile");
    if (fFittprof) std::cout << "Using fit time profile!" << std::endl;
    
    int iDetGeom=TRuntimeParameters::Get().GetParameterI("fiTQun.DetGeomType");
    
    thefit->ReadSharedParams(fiTQun_shared::Get()->GetScatflg(),flgRingType,fFitCprof,fFittprof,iDetGeom);
    TString strDetName = "WCSim";

    thefit->SetQEEffCorr(TRuntimeParameters::Get().GetParameterD(Form("fiTQun.QEEffCorr%s",strDetName.Data())));

    TRuntimeParameters::Get().PrintListOfParameters();

    double cqtot = fiTQun_shared::Get()->GetTotqCoef();


    std::cout << "Successfully instantiated fiTQun." << std::endl;
    
    // 5. Event Loop
    int n_events_to_process = 1;
    if (WCSimWrap::Get()->NEvt() < n_events_to_process) {
        n_events_to_process = WCSimWrap::Get()->NEvt();
    }

    for (int i_event = 0; i_event < n_events_to_process; ++i_event) {
        std::cout << "\n=================================================" << std::endl;
        std::cout << "Processing Event: " << i_event << std::endl;
        
        // Load the event from the WCSim file
        if (WCSimWrap::Get()->LoadEntry(i_event) < 0) {
            std::cerr << "Error loading event " << i_event << std::endl;
            continue;
        }

        // Initialize fiTQun's internal state for the current event.
        // This is a CRUCIAL step. It copies the TQ (time/charge) info 
        // from WCSimWrap into the data structures (common blocks) used by the fitters.
        thefit->InitEvent(i_event);
        thefit->SetTimeWindow(0,-9999,9999);

        // 6. Run a Fit
        //    Let's run a single-ring muon fit as an example.
        std::cout << "Running single-ring muon fit..." << std::endl;

        double snglTrkParams[fiTQun_shared::nSnglTrkParams] = {0.0, 150, -42.47625, 0.0, 1.571, -1.571, 340.0, 0.0};
        int pc_flg;
        int seed_flg = 1; // Vertex prefit
        
        // This is the main fitting call for a single ring muon hypothesis.
        thefit->Resetmu();
        thefit->Do1RFit(ie, snglTrkParams, pc_flg, seed_flg);
        thefit->Resetmu();
        seed_flg = 2;
        thefit->Do1RFit(imu, snglTrkParams, pc_flg, seed_flg);

        std::cout << "Fit finished." << std::endl;

        // =======================================================================
        // Example: Calculate expected PMT charge for specific track parameters
        // =======================================================================
        // Parameters: {x(cm), y(cm), z(cm), t(ns), theta(rad), phi(rad), p(MeV/c), dconv(cm)}
        // for WCTE, (x,y,z)-->(x,-z,y)
        double manualTrackParams[8] = {0.0, 150, -42.47625, 0.0, 1.571, -1.571, 340.0, 0.0};
        
        int nPMTs = WCSimWrap::Get()->NPMT();
        double* expectedCharges = new double[nPMTs];
        
        // Calculate expected charges for Electron hypothesis (ie)
        // This fills expectedCharges with direct + scattered light
        thefit->Get1Rmudist(imu, manualTrackParams, expectedCharges);
        
        std::cout << "Manual Calculation - Expected charge for PMT 0: " << expectedCharges[0] << std::endl;
        std::cout << "Observed charge for PMT 0: " << fiTQun_shared::chrg[0] << std::endl;

        TH1D* hist_pmtQ_fiTQun = new TH1D("pmtQ_fiTQun","pmtQ_fiTQun",nPMTs,0,nPMTs);
        TH1D* hist_pmtQ_data = new TH1D("pmtQ_data","pmtQ_data",nPMTs,0,nPMTs);
        double max_r = 307.5926/2;
        double max_z = 271.4235/2;
        TH2D* hist_event_display_fiTQun = new TH2D("Charges_fiTQun","Charges_fiTQun",250,-TMath::Pi()*max_r,TMath::Pi()*max_r,250,-max_z-2*max_r,max_z+2*max_r);
        TH2D* hist_event_display_data = new TH2D("Charges_data","Charges_data",250,-TMath::Pi()*max_r,TMath::Pi()*max_r,250,-max_z-2*max_r,max_z+2*max_r);
        for (int iPMT = 0; iPMT < nPMTs; ++iPMT) {
            hist_pmtQ_fiTQun->SetBinContent(iPMT+1, expectedCharges[iPMT]);
            hist_pmtQ_data->SetBinContent(iPMT+1, fiTQun_shared::chrg[iPMT]);

            hist_event_display_fiTQun->Fill(eventDiplayXY[iPMT][0],eventDiplayXY[iPMT][1],expectedCharges[iPMT]);
            hist_event_display_data->Fill(eventDiplayXY[iPMT][0],eventDiplayXY[iPMT][1],fiTQun_shared::chrg[iPMT]);
        
        }
        hist_pmtQ_fiTQun->Draw();
        hist_pmtQ_data->SetLineColor(kRed);
        hist_pmtQ_data->Draw("same");
        gPad->Update();
        gPad->SaveAs("test_1s.pdf");
        hist_event_display_fiTQun->Draw("colz");
        gPad->Update();
        gPad->SaveAs("test_2s.pdf");
        hist_event_display_data->Draw("colz");
        gPad->Update();
        gPad->SaveAs("test_3s.pdf");


        
        delete[] expectedCharges;
        // =======================================================================

        // 7. Access Results
        //    The results are stored in the Fortran-style common blocks defined in fitqunoutC.h.
        //    We can access them directly. `fitqun1r_` holds the results for single-ring fits.
        //    The results for the electron fit are at index `ie`.
        //
        //    Note: The results are stored per sub-event. Since we are not running the sub-event
        //    finder, the results will be in the first slot (index 0).
        int sub_event_index = 0;
        float momentum = fitqun1r_.fq1rmom[sub_event_index][imu];
        float nll = fitqun1r_.fq1rnll[sub_event_index][imu];
        
        std::cout << "--- Fit Results (Muon Hypothesis) ---" << std::endl;
        std::cout << "Reconstructed Momentum: " << momentum << " MeV/c" << std::endl;
        std::cout << "Negative Log-Likelihood: " << nll << std::endl;
        std::cout << "Reconstructed Vertex (x, y, z, t, dirx, diry, dirz): " 
                  << fitqun1r_.fq1rpos[sub_event_index][imu][0] << ", "
                  << fitqun1r_.fq1rpos[sub_event_index][imu][1] << ", "
                  << fitqun1r_.fq1rpos[sub_event_index][imu][2] << ", "
                  << fitqun1r_.fq1rt0[sub_event_index][imu] << ", "
                  << fitqun1r_.fq1rdir[sub_event_index][imu][0] << ", "
                  << fitqun1r_.fq1rdir[sub_event_index][imu][1] << ", "
                  << fitqun1r_.fq1rdir[sub_event_index][imu][2] << ", "
                  << std::endl;
    }

    delete thefit;
    std::cout << "\nMacro finished." << std::endl;
}
