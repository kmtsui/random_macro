/// Modified from TMVAClassification.C
/// This macro provides examples for the training and testing of the
/// TMVA DNN classifiers.
///

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int TMVA_DNN_Train( TString myMethodList = "" )
{
   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc

   // Methods to be processed can be given as an argument; use format:
   //
   //     mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

#ifdef R__HAS_TMVAGPU
   Use["DNN_GPU"]         = 1; // CUDA-accelerated DNN training.
#else
   Use["DNN_GPU"]         = 0;
#endif

#ifdef R__HAS_TMVACPU
   Use["DNN_CPU"]         = 0; // Multi-core accelerated DNN.
#else
   Use["DNN_CPU"]         = 0;
#endif

   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return 1;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // Here the preparation phase begins

   // Read training and test data
   // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
   TFile *input1(0), *input2(0);
   TString fname_Eb0 = "/hepstore/kmtsui/T2K/NuWro/nuwro/output/RES_out_RFG_C_noPDD.root";
   TString fname_Eb27 = "/hepstore/kmtsui/T2K/NuWro/nuwro/output/RES_out_RFG_C_noPDD_Eb27.root";
   
   input1 = TFile::Open( fname_Eb0 ); 
   input2 = TFile::Open( fname_Eb27 ); 
   
   if (!input1 || !input2) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassification       : Using input files: " << input1->GetName() << ", "<< input2->GetName() << std::endl;

   // Register the training and test trees

   TTree *signalTree     = (TTree*)input1->Get("RES");
   TTree *background     = (TTree*)input2->Get("RES");

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   TMVA::DataLoader *dataloader=new TMVA::DataLoader("dataset");
   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
   dataloader->AddVariable( "Enu", 'F' );
   dataloader->AddVariable( "mumom", 'F' );
   dataloader->AddVariable( "mucostheta", 'F' );
   dataloader->AddVariable( "Q2", 'F' );
   dataloader->AddVariable( "pimom", 'F' );
   dataloader->AddVariable( "picostheta", 'F' );
   dataloader->AddVariable( "costhetamupi", 'F' );

   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables

   //dataloader->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
   //dataloader->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );


   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;

   // You can add an arbitrary number of signal or background trees
   dataloader->AddSignalTree    ( signalTree,     signalWeight );
   dataloader->AddBackgroundTree( background, backgroundWeight );

   // To give different trees for training and testing, do as follows:
   //
   //     dataloader->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
   //     dataloader->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );

   // Use the following code instead of the above two or four lines to add signal and background
   // training and test events "by hand"
   // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
   //      variable definition, but simply compute the expression before adding the event
   // ```cpp
   // // --- begin ----------------------------------------------------------
   // std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
   // Float_t  treevars[4], weight;
   //
   // // Signal
   // for (UInt_t ivar=0; ivar<4; ivar++) signalTree->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   // for (UInt_t i=0; i<signalTree->GetEntries(); i++) {
   //    signalTree->GetEntry(i);
   //    for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //    // add training and test events; here: first half is training, second is testing
   //    // note that the weight can also be event-wise
   //    if (i < signalTree->GetEntries()/2.0) dataloader->AddSignalTrainingEvent( vars, signalWeight );
   //    else                              dataloader->AddSignalTestEvent    ( vars, signalWeight );
   // }
   //
   // // Background (has event weights)
   // background->SetBranchAddress( "weight", &weight );
   // for (UInt_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   // for (UInt_t i=0; i<background->GetEntries(); i++) {
   //    background->GetEntry(i);
   //    for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //    // add training and test events; here: first half is training, second is testing
   //    // note that the weight can also be event-wise
   //    if (i < background->GetEntries()/2) dataloader->AddBackgroundTrainingEvent( vars, backgroundWeight*weight );
   //    else                                dataloader->AddBackgroundTestEvent    ( vars, backgroundWeight*weight );
   // }
   // // --- end ------------------------------------------------------------
   // ```
   // End of tree registration

   // Set individual event weights (the variables must exist in the original TTree)
   // -  for signal    : `dataloader->SetSignalWeightExpression    ("weight1*weight2");`
   // -  for background: `dataloader->SetBackgroundWeightExpression("weight1*weight2");`
   //dataloader->SetBackgroundWeightExpression( "weight" );

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

   // Tell the dataloader how to use the training and testing events
   //
   // If no numbers of events are given, half of the events in the tree are used
   // for training, and the other half for testing:
   //
   //    dataloader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
   //
   // To also specify the number of testing events, use:
   //
   //    dataloader->PrepareTrainingAndTestTree( mycut,
   //         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
   //dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
   //                                       "nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V" );
   dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
                                           "nTrain_Signal=100000:nTest_Signal=100000:nTrain_Background=100000:nTest_Background=100000:SplitMode=Random:NormMode=NumEvents:!V" );

   // ### Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation

   // Multi-architecture DNN implementation.
   if (Use["DNN_CPU"] or Use["DNN_GPU"]) {
      // General layout.
      TString layoutString ("Layout=RELU|100,RELU|100,RELU|100,SIGMOID");

      // Training strategies.
      TString training0("Repetitions=1,"
                        "ConvergenceSteps=200,BatchSize=256,TestRepetitions=10,"
                        "Regularization=L2,"
                        "DropConfig=0.2+0.2+0.2, Multithreading=True");
      TString training1("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                        "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
      TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
                        "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                        "WeightDecay=1e-4,Regularization=L2,"
                        "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
      TString trainingStrategyString ("TrainingStrategy=");
      trainingStrategyString += training0;// + "|" + training1 + "|" + training2;

      // General Options.
      TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                          "WeightInitialization=XAVIERUNIFORM:CreateMVAPdfs=True");
      dnnOptions.Append (":"); dnnOptions.Append (layoutString);
      dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);

      // Cuda implementation.
      if (Use["DNN_GPU"]) {
         TString gpuOptions = dnnOptions + ":Architecture=GPU";
         factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_GPU", gpuOptions);
      }
      // Multi-core CPU implementation.
      if (Use["DNN_CPU"]) {
         TString cpuOptions = dnnOptions + ":Architecture=CPU";
         factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU", cpuOptions);
      }
   }


   // For an example of the category classifier usage, see: TMVAClassificationCategory
   //
   // --------------------------------------------------------------------------------------------------
   //  Now you can optimize the setting (configuration) of the MVAs using the set of training events
   // STILL EXPERIMENTAL and only implemented for BDT's !
   //
   //     factory->OptimizeAllMethods("SigEffAtBkg0.01","Scan");
   //     factory->OptimizeAllMethods("ROCIntegral","FitGA");
   //
   // --------------------------------------------------------------------------------------------------

   // Now you can tell the factory to train, test, and evaluate the MVAs
   //
   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;
   delete dataloader;
   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

   return 0;
}

int main( int argc, char** argv )
{
   // Select methods (don't look at this code - not of interest)
   TString methodList;
   for (int i=1; i<argc; i++) {
      TString regMethod(argv[i]);
      if(regMethod=="-b" || regMethod=="--batch") continue;
      if (!methodList.IsNull()) methodList += TString(",");
      methodList += regMethod;
   }
   return TMVA_DNN_Train(methodList);
}
