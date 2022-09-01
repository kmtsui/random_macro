/// Modified from TMVAClassificationApplication.C
/// This macro provides a simple example on how to use the trained classifiers
/// within an analysis module

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;

void TMVA_DNN_Reweight( TString myMethodList = "" )
{

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   Use["DNN_CPU"] = 0;         // CUDA-accelerated DNN training.
   Use["DNN_GPU"] = 1;         // Multi-core accelerated DNN.

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t Enu, mumom, mucostheta, Q2, pimom, picostheta, costhetamupi;
   reader->AddVariable( "Enu", &Enu );
   reader->AddVariable( "mumom", &mumom );
   reader->AddVariable( "mucostheta", &mucostheta );
   reader->AddVariable( "Q2", &Q2 );
   reader->AddVariable( "pimom", &pimom );
   reader->AddVariable( "picostheta", &picostheta );
   reader->AddVariable( "costhetamupi", &costhetamupi );


   // Book the MVA methods

   TString dir    = "dataset/weights/";
   TString prefix = "TMVAClassification";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile );
      }
   }



   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //
   TFile *input1(0), *input2(0);
   TString fname_Eb0 = "/hepstore/kmtsui/T2K/NuWro/nuwro/output/RES_out_RFG_C_noPDD.root";
   TString fname_Eb27 = "/hepstore/kmtsui/T2K/NuWro/nuwro/output/RES_out_RFG_C_noPDD_Eb27.root";
   
   input1 = TFile::Open( fname_Eb0 ); 
   input2 = TFile::Open( fname_Eb27 ); 
   
   if (!input1 || !input2) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input files: " << input1->GetName() << ", "<< input2->GetName() << std::endl;

   // Event loop

   // Prepare the event tree
   // - Here the variable names have to corresponds to your tree
   // - You can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   std::cout << "--- Select signal sample" << std::endl;
   TTree* t_Eb0 = (TTree*)input1->Get("RES");
   t_Eb0->SetBranchAddress( "Enu", &Enu );
   t_Eb0->SetBranchAddress( "mumom", &mumom );
   t_Eb0->SetBranchAddress( "mucostheta", &mucostheta );
   t_Eb0->SetBranchAddress( "Q2", &Q2 );
   t_Eb0->SetBranchAddress( "pimom", &pimom );
   t_Eb0->SetBranchAddress( "picostheta", &picostheta );
   t_Eb0->SetBranchAddress( "costhetamupi", &costhetamupi );


   TH1D* h_weight = new TH1D("h_weight","",100,0,10);
   TH1D* h_Q2_Eb0 = new TH1D("h_Q2_Eb0","",20,0,2);
   TH1D* h_Q2_Eb27 = new TH1D("h_Q2_Eb27","",20,0,2);
   TH1D* h_Q2_Eb27_rw = new TH1D("h_Q2_Eb27_rw","",20,0,2);

   std::cout << "--- Processing: " << t_Eb0->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<t_Eb0->GetEntries();ievt++) {

      if (ievt>1000000) break;

      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      t_Eb0->GetEntry(ievt);

      double val = reader->EvaluateMVA("DNN_GPU method");
      double weight = (1.-val)/val;

      h_weight->Fill(weight);
      h_Q2_Eb0->Fill(-Q2/1.e6);
      h_Q2_Eb27_rw->Fill(-Q2/1.e6,weight);
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   TCanvas* c1 = new TCanvas();
   h_weight->Draw();

   TCanvas* c2 = new TCanvas();
   h_Q2_Eb0->Draw();
   h_Q2_Eb27_rw->SetLineColor(kRed);
   h_Q2_Eb27_rw->Draw("same");

   // Write histograms

   TFile *target  = new TFile( "TMVApp.root","RECREATE" );
   h_weight->Write();
   h_Q2_Eb0->Write();
   h_Q2_Eb27_rw->Write();

   target->Close();

   std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;

   delete reader;

   std::cout << "==> TMVAClassificationApplication is done!" << std::endl << std::endl;
}

int main( int argc, char** argv )
{
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   TMVA_DNN_Reweight(methodList);
   return 0;
}
