#include "PDAnalysis/EleMVACalibrationProducer/interface/EleMVACalibrationProducer.h"

// #include <vector>
// #include <string>
#include <math.h>
// #include <iostream>
// #include <sstream>
// 
// #include "TFile.h"
// #include "TTree.h"
// #include "TBranch.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
// #include "TMVA/Tools.h"
#include "TMVA/PyMethodBase.h"
// #include "TH1.h"
#include "TF1.h"
// #include "TString.h"
#include "TMath.h"
#include "TEfficiency.h"
// #include "TGraph.h"
// #include "TGraphErrors.h"
// #include "TGraphAsymmErrors.h"
// #include "TCanvas.h"
#include "Math/MinimizerOptions.h"
#include "TFitResult.h"
#include "TPaveStats.h"
// 

void EleMVACalibrationProducer::beginJob()
{
//   treeName="EleMVAsecondTree";
  initTree();
  
  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(10); //ksiourmen
  gStyle->SetOptFit(1); //pcev
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls( 15000 );
  
  getUserParameter( "verbose", verbose );
  getUserParameter("writeOutput", writeOutput);
  getUserParameter("treeListName", treeListName);
//   getUserParameter("sampleName", sampleName);
  getUserParameter("inputTreeName", inputTreeName);
  getUserParameter("mvaMethod", mvaMethod);
  getUserParameter("useTightSelection", useTightSelection);
  getUserParameter("weightControlPlots", weightControlPlots);
  getUserParameter("useSyst", useSyst);
  getUserParameter("nBinsCal", nBinsCal);
  getUserParameter("systematics2017", systematics2017);
  getUserParameter("systematics2018", systematics2018);
  getUserParameter("mvaInputPath", mvaInputPath);
  getUserParameter("eleIdWP", eleIdWP);
  getUserParameter("elePtWP", elePtWP);
  getUserParameter("eleDzWP", eleDzWP);
  getUserParameter("eleDxyWP", eleDxyWP);
  getUserParameter("eleDRBWP", eleDRBWP);
  getUserParameter("mvaValueMinCutOff", mvaValueMinCutOff);
  getUserParameter("mvaValueMaxCutOff", mvaValueMaxCutOff);
  
//   sampleName = inputFileName;
  
//   do
//   {
//     sampleName = sampleName.substr(sampleName.find("/") + 1);
//   } 
//   while(sampleName.find("/") != std::string::npos);
//   
//   sampleName = sampleName.substr(0, sampleName.find("."));

  std::cout << "I N F O : EleMVACalibrationProducer::beginJob() - Parameters:" << std::endl;
  std::cout << "          treeListName = " << treeListName << std::endl;
  std::cout << "          inputTreeName = " << inputTreeName << std::endl;
//   std::cout << "          sampleName = " << sampleName << std::endl;
  std::cout << "          mvaMethod = " << mvaMethod << std::endl;
  std::cout << "          mvaInputPath = " << mvaInputPath << std::endl;
  std::cout << "          useTightSelection = " << useTightSelection << std::endl;
  std::cout << "          nBinsCal = " << nBinsCal << std::endl;
  std::cout << "          mvaValueMinCutOff = " << mvaValueMinCutOff << std::endl;
  std::cout << "          mvaValueMaxCutOff = " << mvaValueMaxCutOff << std::endl;
  
//   auto *inputFile = new TFile(inputFileName.c_str());
//   auto *inputTree = (TTree*)inputFile->Get(inputTreeName.c_str());
//   if(inputFile->IsZombie())
//   {
//     return;
//   }

  if(treeListName.find("Bs") != std::string::npos)
  {
    process = "BsJPsiPhi";
    minMass = 5.20;
    maxMass = 5.65;
    dirPath = "./Bs";
    if(treeListName.find("DG0") != std::string::npos)
    {
      process += "DG0";
    }
  }
  if(treeListName.find("Bu") != std::string::npos)
  {
    process = "BuJPsiK";
    minMass = 5.10;
    maxMass = 5.65;
    x1MassSB = 5.44;
    x2MassSB = 5.64;
    dirPath = "./Bu";
  }
  
  if(treeListName.find("MC") != std::string::npos)
  {
    process = process + "MC";
    dirPath += "MC";
  }
  if(treeListName.find("Data") != std::string::npos)
  {
    process = process + "Data";
    dirPath += "Data";
    isData = true;
  }
  
  if(treeListName.find("2017") != std::string::npos)
  {
    process = process + "2017";
    dirPath += "2017";
    if(useSyst)
    {
      systematics = systematics2017;
    }
  }
  
  if(treeListName.find("2018") != std::string::npos)
  {
    process = process + "2018";
    dirPath += "2018";
    if(useSyst)
    {
      systematics = systematics2018;
    }
  }
  
  gSystem->mkdir(dirPath.c_str());
  
  if(verbose)
  {
    std::cout << "EleMVACalibrationProducer::beginJob() - FILES OPENED.\n";
  }
  
  // PER-EVENT VARIABLES
  evtW = {-1., -1.}; // per-event mistag rate
  totalP = 0.; // total tagging power
  totalPBinned = 0.;
  totalPBinnedErr = 0.;
  
  binSizeCal = (1.-0.)/nBinsCal; // (1-0) -> mva score range
  
  //COMPUTE MVA
  TMVA::PyMethodBase::PyInitialize();
  reader = new TMVA::Reader("!Color:Silent");
  double mvaValue = -1.;
  float eleDxyDiff;
  float eleDzDiff;
  reader->AddVariable("elePt", &elePt);
  reader->AddVariable("eleEta", &eleEta);
  reader->AddVariable("eleDxy", &eleDxy);
  reader->AddVariable("eleExy", &eleExy);
  reader->AddVariable("eleDz", &eleDz);
  reader->AddVariable("eleEz", &eleEz);
//   reader->AddVariable("eleIDNIV2RawVal", &eleIDNIV2RawVal);
  reader->AddVariable("eleIDNIV2Val", &eleIDNIV2Val);
  reader->AddVariable("eleDRB", &eleDRB);
  reader->AddVariable("elePFIsoScaled", &elePFIsoScaled);
  reader->AddVariable("eleConeCleanPt", &eleConeCleanPt);
  reader->AddVariable("eleConeCleanPtRel", &eleConeCleanPtRel);
  reader->AddVariable("eleConeCleanDR", &eleConeCleanDR);
  reader->AddVariable("eleConeCleanEnergyRatio", &eleConeCleanEnergyRatio);
  reader->AddVariable("eleConeCleanQ", &eleConeCleanQ);
//   reader->AddVariable("eleConeCleanAvgDxy", &eleConeCleanAvgDxy);
//   reader->AddVariable("eleDxy-eleConeCleanAvgDxy", &eleDxyDiff);
//   reader->AddVariable("eleDz-eleConeCleanAvgDz", &eleDzDiff);
//   reader->AddVariable("eleConeCleanStdDevDxy", &eleConeCleanStdDevDxy);
//   reader->AddVariable("eleConeCleanStdDevDz", &eleConeCleanStdDevDz);
  
//   eleDxyDiff = eleDxy-eleConeCleanAvgDxy;
//   eleDzDiff = eleDz-eleConeCleanAvgDz;
  
  reader->BookMVA(mvaMethod, mvaInputPath + "TMVAClassification_" + mvaMethod + ".weights.xml");
  
  vEleIdcuts = {-1., -0.9999, -0.9995, -0.999, -0.995, -0.99, -0.95, -0.9, -0.8, -0.7, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 0.85, 0.9, 0.95, 0.99, 1.};
  
  for(unsigned int iCut = 0; iCut < vEleIdcuts.size(); iCut++)
  {
    vEleIdcutsLog.push_back(TMath::Log(vEleIdcuts[iCut]+2));
  }
}



// void EleMVACalibrationProducer::setupNtuple()
// {
//   initTree();
// }



// bool EleMVACalibrationProducer::getEntry(int ientry)
// {
//   currentTree->GetEntry(ientry);
//   return true;
// }



void EleMVACalibrationProducer::book()
{
  autoSavedObject = hNGenB = new TH1D( "hNGenB", "hNGenB", 10, -0.5, 9.5 );
  autoSavedObject = hIsTagRightVsTagTruth = Create2DHistogram<TH2D>("hIsTagRightVsTagTruth", "IsTagRight vs. tag truth", 5, -2.5, 2.5, 5, -2.5, 2.5, "TagTruth", "IsTagRight");
  
  std::vector<double> vMvaScoreBinEdges;
  for(int i = 0; i <= nBinsCal; i++)
  {
    vMvaScoreBinEdges.push_back(i * binSizeCal);
  }
  
  vhMassCalRightTag = CreateVectorOf1DHistograms<TH1D>("hMassCalRightTag", "Mass calibration - right tag", nMassBins, minMass, maxMass, "M [GeV]", "N. events", "mvaScore", vMvaScoreBinEdges);
  autoSavedObject = reinterpret_cast<std::vector<TObject*>* >(vhMassCalRightTag);
  vhMassCalWrongTag = CreateVectorOf1DHistograms<TH1D>("hMassCalWrongTag", "Mass calibration - wrong tag", nMassBins, minMass, maxMass, "M [GeV]", "N. events", "mvaScore", vMvaScoreBinEdges);
  autoSavedObject = reinterpret_cast<std::vector<TObject*>* >(vhMassCalWrongTag);
  
  
//   TH1F *hMassCalRT[nBinsCal];
//   TH1F *hMassCalWT[nBinsCal];
//   double *wCalc = new double[nBinsCal]; // measured mistag rate
//   double *wCalcEdgeL = new double[nBinsCal];
//   double *wCalcEdgeH = new double[nBinsCal];
  
  for(int i = 0; i < nBinsCal; i++)
  {
    vWCalc.push_back(0.);
    vWCalcLowEdge.push_back(0.);
    vWCalcHighEdge.push_back(0.);
  }
  
  //HISTOGRAMS BOOKING
  nBinsMva = 100;
  autoSavedObject = hMva          = Create1DHistogram<TH1D>("hMva",          "mva value",    nBinsMva, 0.0, 1.0, "MVA output", "N. events");
  autoSavedObject = hMvaRightTag  = Create1DHistogram<TH1D>("hMvaRightTag",  "mva value - right tag events", nBinsMva, 0.0, 1.0, "MVA output", "N. events");
  autoSavedObject = hMvaWrongTag  = Create1DHistogram<TH1D>("hMvaWrongTag",  "mva value - wrong tag events", nBinsMva, 0.0, 1.0, "MVA output", "N. events");
  autoSavedObject = hMassTot      = Create1DHistogram<TH1D>("hMassTot",      "B mass - all events", nMassBins, minMass, maxMass, "M [GeV]", "N. events");
  autoSavedObject = hMassRightTag = Create1DHistogram<TH1D>("hMassRightTag", "B mass - right tag events", nMassBins, minMass, maxMass, "M [GeV]", "N. events");
  autoSavedObject = hMassWrongTag = Create1DHistogram<TH1D>("hMassWrongTag", "B mass - wrong tag events", nMassBins, minMass, maxMass, "M [GeV]", "N. events");
  autoSavedObject = hMassNoTag    = Create1DHistogram<TH1D>("hMassNoTag",    "B mass - no tag events", nMassBins, minMass, maxMass, "M [GeV]", "N. events");

  for(unsigned int i = 0; i < vEleIdcuts.size(); i++)
  {
    std::ostringstream os;
    float value = vEleIdcuts.at(i);
    os << value;
    std::string sValue = os.str();
    std::string hName = "hMassRightTag_eleIDNIV2Val_gt_" + sValue;
    std::string hTitle = "B mass - right tag events - eleIDNIV2Val > " + sValue;
    TH1D* hTemp = Create1DHistogram<TH1D>(hName.c_str(), hTitle.c_str(), nMassBins, minMass, maxMass, "M [GeV]", "N. events");
    vhMassRightTag.push_back(hTemp);
    hName = "hMassWrongTag_eleIDNIV2Val_gt_" + sValue;
    hTitle = "B mass - wrong tag events - eleIDNIV2Val > " + sValue;
    hTemp = Create1DHistogram<TH1D>(hName.c_str(), hTitle.c_str(), nMassBins, minMass, maxMass, "M [GeV]", "N. events");
    vhMassWrongTag.push_back(hTemp);
  }
  autoSavedObject = reinterpret_cast<std::vector<TObject*>* >(&vhMassRightTag);
  autoSavedObject = reinterpret_cast<std::vector<TObject*>* >(&vhMassWrongTag);

  autoSavedObject = hElePt                      = Create1DHistogram<TH1D>("hElePt",                   "electron p_{T}",                    100,   0. ,  50. , "p_{T} [GeV]",      "N. electrons");
  autoSavedObject = hEleEta                     = Create1DHistogram<TH1D>("hEleEta",                  "electron #eta",                     100,  -3. ,   3. , "#eta",              "N. electrons");
  autoSavedObject = hEleDxy                     = Create1DHistogram<TH1D>("hEleDxy",                  "electron d_{xy}",                   100,  -0.1,   0.1, "d_{xy} [cm]",       "N. electrons");
  autoSavedObject = hEleExy                     = Create1DHistogram<TH1D>("hEleExy",                  "electron d_{xy} error",             100,   0. ,   0.1, "e_{xy} [cm]",       "N. electrons");
  autoSavedObject = hEleDz                      = Create1DHistogram<TH1D>("hEleDz",                   "electron d_{z}",                    100,  -0.2,   0.2, "d_{z} [cm]",        "N. electrons");
  autoSavedObject = hEleEz                      = Create1DHistogram<TH1D>("hEleEz",                   "electron d_{z} error",              100,   0. ,   0.1, "e_{z} [cm]",        "N. electrons");
  autoSavedObject = hEleIDNIV2Val               = Create1DHistogram<TH1D>("hEleIDNIV2Val",            "electron ID non-isoV2 val.",        100,  -1. ,   1. , "EleIDNIV2Val",      "N. electrons");
  autoSavedObject = hEleDRB                     = Create1DHistogram<TH1D>("hEleDRB",                  "electron #Delta(R) from B",         100,   0. ,   6. , "#Delta(R)",         "N. electrons");
  autoSavedObject = hElePFIsoScaled             = Create1DHistogram<TH1D>("hElePFIsoScaled",          "electron PF isolation scaled",      100,   0. ,  10. , "PF Isolation",     "N. electrons");
  autoSavedObject = hEleConeCleanPt             = Create1DHistogram<TH1D>("hEleConeCleanPt",          "p_{T} of clean electron cone",      100,   0. , 100. , "p_{T} [GeV]",     "N. electrons");
  autoSavedObject = hEleConeCleanPtRel          = Create1DHistogram<TH1D>("hEleConeCleanPtRel",       "p_{T,rel} of clean electron cone",  100, -50. ,  50. , "p_{T,rel}  [GeV]", "N. electrons");
  autoSavedObject = hEleConeCleanDR             = Create1DHistogram<TH1D>("hEleConeCleanDR",          "#Delta(R) of clean electron cone",  100,   0. ,   0.3, "#Delta(R)",         "N. electrons");
  autoSavedObject = hEleConeCleanEnergyRatio    = Create1DHistogram<TH1D>("hEleConeCleanEnergyRatio", "E ratio of clean electron cone",    100,   0. ,   4. , "E(e)/E(cone)",      "N. electrons");
  autoSavedObject = hEleConeCleanQ              = Create1DHistogram<TH1D>("hEleConeCleanQ",           "Weighted Q of clean electron cone", 100,  -1. ,   1. , "weighted Q(cone)",  "N. electrons");
  
  autoSavedObject = gNRightTagVsEleIdCut = CreateGraph<TGraph>("gNRightTagVsEleIdCut", "N. of right tags vs. eleId cut", "EleID cut", "N. of right tags", vEleIdcutsLog.size());
  autoSavedObject = gNWrongTagVsEleIdCut = CreateGraph<TGraph>("gNWrongTagVsEleIdCut", "N. of wrong tags vs. eleId cut", "EleID cut", "N. of wrong tags", vEleIdcutsLog.size());
  autoSavedObject = gNNoTagVsEleIdCut = CreateGraph<TGraph>("gNNoTagVsEleIdCut", "N. of no tags vs. eleId cut", "EleID cut", "N. of no tags", vEleIdcutsLog.size());
  autoSavedObject = gBaseEffVsEleIdCut = CreateGraph<TGraph>("gBaseEffVsEleIdCut", "Base efficiency vs. eleId cut", "EleID cut", "Base eff.", vEleIdcutsLog.size());
  autoSavedObject = gBaseRightTagEffVsEleIdCut = CreateGraph<TGraph>("gBaseRightTagEffVsEleIdCut", "Base right tag efficiency vs. eleId cut", "EleID cut", "Base right tag eff.", vEleIdcutsLog.size());
  autoSavedObject = gBaseWrongTagEffVsEleIdCut = CreateGraph<TGraph>("gBaseWrongTagEffVsEleIdCut", "Base wrong tag efficiency vs. eleId cut", "EleID cut", "Base wrong tag eff.", vEleIdcutsLog.size());
  autoSavedObject = gBaseMistagVsEleIdCut = CreateGraph<TGraph>("gBaseMistagVsEleIdCut", "Base mistag vs. eleId cut", "EleID cut", "Base mistag", vEleIdcutsLog.size());
  autoSavedObject = gBaseDilutionVsEleIdCut = CreateGraph<TGraph>("gBaseDilutionVsEleIdCut", "Base dilution vs. eleId cut", "EleID cut", "Base dilution", vEleIdcutsLog.size());
  autoSavedObject = gBasePowerVsEleIdCut = CreateGraph<TGraph>("gBasePowerVsEleIdCut", "Base power vs. eleId cut", "EleID cut", "Base power", vEleIdcutsLog.size());

  autoSavedObject = hCalFitStatusRightTagVsMvaScore = Create1DHistogram<TH1D>("hCalFitStatusRightTagVsMvaScore", "Status of right tag calibration fits", nBinsCal, 0., 1., "mistag calc.", "Fit status");
  autoSavedObject = hCalFitStatusWrongTagVsMvaScore = Create1DHistogram<TH1D>("hCalFitStatusWrongTagVsMvaScore", "Status of wrong tag calibration fits", nBinsCal, 0., 1., "mistag calc.", "Fit status");
  
  autoSavedObject = hCalFitChi2NDOFRightTagVsMvaScore = Create1DHistogram<TH1D>("hCalFitChi2NDOFRightTagVsMvaScore", "#chi^2/NDOF of right tag calibration fits", nBinsCal, 0., 1., "mistag calc.", "#chi^2/NDOF");
  autoSavedObject = hCalFitChi2NDOFWrongTagVsMvaScore = Create1DHistogram<TH1D>("hCalFitChi2NDOFWrongTagVsMvaScore", "#chi^2/NDOF of wrong tag calibration fits", nBinsCal, 0., 1., "mistag calc.", "#chi^2/NDOF");
    
  autoSavedObject = gNRightTagVsBaseEff = CreateGraph<TGraph>("gNRightTagVsBaseEff", "N. of right tags vs. base eff.", "Base eff.", "N. of right tags", vEleIdcutsLog.size());
  autoSavedObject = gNWrongTagVsBaseEff = CreateGraph<TGraph>("gNWrongTagVsBaseEff", "N. of wrong tags vs. base eff.", "Base eff.", "N. of wrong tags", vEleIdcutsLog.size());
  autoSavedObject = gNNoTagVsBaseEff = CreateGraph<TGraph>("gNNoTagVsBaseEff", "N. of no tags vs. base eff.", "Base eff.", "N. of no tags", vEleIdcutsLog.size());
  autoSavedObject = gBaseEffVsBaseEff = CreateGraph<TGraph>("gBaseEffVsBaseEff", "Base efficiency vs. base eff.", "Base eff.", "Base eff.", vEleIdcutsLog.size());
  autoSavedObject = gBaseRightTagEffVsBaseEff = CreateGraph<TGraph>("gBaseRightTagEffVsBaseEff", "Base right tag efficiency vs. base eff.", "Base eff.", "Base right tag eff.", vEleIdcutsLog.size());
  autoSavedObject = gBaseWrongTagEffVsBaseEff = CreateGraph<TGraph>("gBaseWrongTagEffVsBaseEff", "Base wrong tag efficiency vs. base eff.", "Base eff.", "Base wrong tag eff.", vEleIdcutsLog.size());
  autoSavedObject = gBaseMistagVsBaseEff = CreateGraph<TGraph>("gBaseMistagVsBaseEff", "Base mistag vs. base eff.", "Base eff.", "Base mistag", vEleIdcutsLog.size());
  autoSavedObject = gBaseDilutionVsBaseEff = CreateGraph<TGraph>("gBaseDilutionVsBaseEff", "Base dilution vs. base eff.", "Base eff.", "Base dilution", vEleIdcutsLog.size());
  autoSavedObject = gBasePowerVsBaseEff = CreateGraph<TGraph>("gBasePowerVsBaseEff", "Base power vs. base eff.", "Base eff.", "Base power", vEleIdcutsLog.size());

  autoSavedObject = gBaseWrongTagEffVsBaseRightTagEff = CreateGraph<TGraph>("gBaseWrongTagEffVsBaseRightTagEff", "Base wrong tag eff. vs. base right tag eff.", "Base right tag eff.", "Base wrong tag eff.", vEleIdcutsLog.size());
  
  if(verbose)
  {
    std::cout << "EleMVACalibrationProducer::book() - BOOKING COMPLETED.\n";
  }
}



bool EleMVACalibrationProducer::analyze(int entry, int event_file, int event_tot)
{
  
  if(entry%100000==0)
  {
    std::cout << "----- at event " << entry << std::endl;
  }
  
  //EVENT SELECTION
  if(!JPsiTrkTrkHltBit)
  {
    return false;
  }
  if(useTightSelection && !tightB)
  {
    return false;
  }
  
  hMassTot->Fill(BMass, evtWeight);
  
  //ELECTRON SELECTION
  if(!eleSelected)
  {
    hMassNoTag->Fill(BMass, evtWeight);
    return false;
  }
  if(eleIDNIV2Val <= eleIdWP)
  {
    hMassNoTag->Fill(BMass, evtWeight);
    return false;
  }
  if(elePt <= elePtWP)
  {
    hMassNoTag->Fill(BMass, evtWeight);
    return false;
  }
  if(fabs(eleDz) >= eleDzWP)
  {
    hMassNoTag->Fill(BMass, evtWeight);
    return false;
  }
  if(fabs(eleDxy) >= eleDxyWP)
  {
    hMassNoTag->Fill(BMass, evtWeight);
    return false;
  }
  if(eleDRB <= eleDRBWP)
  {
    hMassNoTag->Fill(BMass, evtWeight);
    return false;
  }

  if(eleDxy<=-999. || eleDxy>=999.)
  {
    hMassNoTag->Fill(BMass, evtWeight);
    return false;
  }
  
  if(eleConeCleanAvgDxy<=-999. || eleConeCleanAvgDxy>=999)
  {
    hMassNoTag->Fill(BMass, evtWeight);
    return false;
  }

  if(eleConeCleanStdDevDxy<=0. || eleConeCleanStdDevDxy>=999)
  {
    hMassNoTag->Fill(BMass, evtWeight);
    return false;
  }

  //TAGGING
  double mvaValue = reader->EvaluateMVA(mvaMethod);
  
  if(mvaValue < mvaValueMinCutOff)
  {
    hMassNoTag->Fill(BMass, evtWeight);
    return false;
  }
  if(mvaValue > mvaValueMaxCutOff)
  {
    hMassNoTag->Fill(BMass, evtWeight);
    return false;
  }

  hMva->Fill(mvaValue, evtWeight);
  evtW[0] = 1 - mvaValue;
  totalP += pow(1. - 2. * evtW[0], 2) * evtWeight;
  
  int evtTag = -1 * eleCharge;
  bool isTagRight = TMath::Sign(1, BidGen) == evtTag;
  if(tagTruth)
    tagTruth = 0;
  else
    tagTruth = 1;

  hIsTagRightVsTagTruth->Fill(tagTruth, isTagRight);
  if(isTagRight != tagTruth)
  {
    std::cout << "E R R O R! EleMVACalibrationProducer::analyze(): isTagRight != tagTruth!\n";
//     exit(1);
  }
  
  for(int j = 0; j < nBinsCal; j++)
  {
    if((evtW[0] >= (double)j * binSizeCal) && (evtW[0] < ((double)j * binSizeCal + binSizeCal)))
    {
      if(isTagRight)
      {
        vhMassCalRightTag->at(j)->Fill(BMass, evtWeight);
      }
      else
      {
        vhMassCalWrongTag->at(j)->Fill(BMass, evtWeight);
      }
      vWCalc[j] += evtW[0] * evtWeight;
      break;
    }
  }
  
//   int iEleIdCut = -999.;
//   for(int i = vEleIdcuts.size() - 1; i > -1; i--)
//   {
//     if(eleIDNIV2Val >vEleIdcuts.at(i))
//     {
//       iEleIdCut = i;
//       break;
//     }
//   }
  
  if(isTagRight)
  {
    hMvaRightTag->Fill(mvaValue, evtWeight);
    hMassRightTag->Fill(BMass, evtWeight);
    for(unsigned int i = 0; i < vEleIdcuts.size(); i++)
    {
      if(eleIDNIV2Val >vEleIdcuts.at(i))
      {
        vhMassRightTag[i]->Fill(BMass, evtWeight);
      }
    }
  }
  else
  {
    hMvaWrongTag->Fill(mvaValue, evtWeight);
    hMassWrongTag->Fill(BMass, evtWeight);
    for(unsigned int i = 0; i < vEleIdcuts.size(); i++)
    {
      if(eleIDNIV2Val >= vEleIdcuts.at(i))
      {
        vhMassWrongTag[i]->Fill(BMass, evtWeight);
      }
    }
  }
  
  hNGenB->Fill(nGenB);
  
  // Fill control histograms
  double evtWeightCtrl = (weightControlPlots? evtWeight : 1);
  hElePt->Fill(elePt, evtWeightCtrl);
  hEleEta->Fill(eleEta, evtWeightCtrl);
  hEleDxy->Fill(eleDxy, evtWeightCtrl);
  hEleExy->Fill(eleExy, evtWeightCtrl);
  hEleDz->Fill(eleDz, evtWeightCtrl);
  hEleEz->Fill(eleEz, evtWeightCtrl);
  hEleIDNIV2Val->Fill(eleIDNIV2Val, evtWeightCtrl);
  hEleDRB->Fill(eleDRB, evtWeightCtrl);
  hElePFIsoScaled->Fill(elePFIsoScaled, evtWeightCtrl);
  hEleConeCleanPt->Fill(eleConeCleanPt, evtWeightCtrl);
  hEleConeCleanPtRel->Fill(eleConeCleanPtRel, evtWeightCtrl);
  hEleConeCleanDR->Fill(eleConeCleanDR, evtWeightCtrl);
  hEleConeCleanEnergyRatio->Fill(eleConeCleanEnergyRatio, evtWeightCtrl);
  hEleConeCleanQ->Fill(eleConeCleanQ, evtWeightCtrl);
  
  return true;
}



void EleMVACalibrationProducer::endJob()
{
  // PERFORMANCE OUTPUT
  double nRightTag = hMassRightTag->Integral(); //Integral() takes in consideration event weights
  double nWrongTag = hMassWrongTag->Integral();
  double nNoTag = hMassNoTag->Integral();
  double nTot = hMassTot->Integral();
  
  if(isData)
  { 
    std::cout << "Real data: count events with fit...\n";
    //for data fit mass
    nTot = CountEventsWithFit(hMassTot).first; //Fit of the total histogram need to be called first
    nRightTag = CountEventsWithFit(hMassRightTag).first;
    nWrongTag = CountEventsWithFit(hMassWrongTag).first;
    nNoTag = CountEventsWithFit(hMassNoTag).first;
  }
  
  double effBase = (double)(nRightTag + nWrongTag) / (nRightTag + nWrongTag + nNoTag);
  double wBase = (double)nWrongTag / (nRightTag + nWrongTag);
  double pBase = effBase * pow(1. - 2. * wBase, 2);
  
  std::cout << "nTot = " << nTot << std::endl;
  std::cout << "nRightTag = " << nRightTag << std::endl;
  std::cout << "nWrongTag = " << nWrongTag << std::endl;
  std::cout << "nNoTag = " << nNoTag << std::endl;
  
  std::cout << "Cutting on eleId value:\n";
  
  for(unsigned int i = 0; i < vEleIdcuts.size(); i++)
  {
    double cutValue = vEleIdcuts[i];
    double cutValueLog = vEleIdcutsLog[i];
    double rightTagIntegral = vhMassRightTag[i]->Integral();
    double wrongTagIntegral = vhMassWrongTag[i]->Integral();
    std::cout << "Ele Id cut value = " << cutValue << " - nRightTag = " << rightTagIntegral << " - nWrongTag = " << wrongTagIntegral << std::endl;
    
    double nNoTagVsEleIdCut = nTot - rightTagIntegral - wrongTagIntegral;
    double baseEffVsEleIdCut = (double)(rightTagIntegral + wrongTagIntegral) / nTot;
    double baseRightTagEffVsEleIdCut = (double)(rightTagIntegral) / nTot;
    double baseWrongTagEffVsEleIdCut = (double)(wrongTagIntegral) / nTot;
    double baseMistagVsEleIdCut = (double)(wrongTagIntegral) / (rightTagIntegral + wrongTagIntegral);
    double baseDilutionVsEleIdCut = (1. - 2. * baseMistagVsEleIdCut);
    double basePowerVsEleIdCut = baseEffVsEleIdCut * pow(1. - 2. * baseMistagVsEleIdCut, 2);
    if(baseEffVsEleIdCut == 0)
    {
      baseMistagVsEleIdCut = 0;
      basePowerVsEleIdCut = 0;
    }
    
    gNRightTagVsEleIdCut->SetPoint(i, cutValueLog, rightTagIntegral);
    gNWrongTagVsEleIdCut->SetPoint(i, cutValueLog, wrongTagIntegral);
    gNNoTagVsEleIdCut->SetPoint(i, cutValueLog, nNoTagVsEleIdCut);
    gBaseEffVsEleIdCut->SetPoint(i, cutValueLog, baseEffVsEleIdCut);
    gBaseRightTagEffVsEleIdCut->SetPoint(i, cutValueLog, baseRightTagEffVsEleIdCut);
    gBaseWrongTagEffVsEleIdCut->SetPoint(i, cutValueLog, baseWrongTagEffVsEleIdCut);
    gBaseMistagVsEleIdCut->SetPoint(i, cutValueLog, baseMistagVsEleIdCut);
    gBaseDilutionVsEleIdCut->SetPoint(i, cutValueLog, baseDilutionVsEleIdCut);
    gBasePowerVsEleIdCut->SetPoint(i, cutValueLog, basePowerVsEleIdCut);
    
    gNRightTagVsBaseEff->SetPoint(i, baseEffVsEleIdCut, rightTagIntegral);
    gNWrongTagVsBaseEff->SetPoint(i, baseEffVsEleIdCut, wrongTagIntegral);
    gNNoTagVsBaseEff->SetPoint(i, baseEffVsEleIdCut, nNoTagVsEleIdCut);
    gBaseEffVsBaseEff->SetPoint(i, baseEffVsEleIdCut, baseEffVsEleIdCut);
    gBaseRightTagEffVsBaseEff->SetPoint(i, baseEffVsEleIdCut, baseRightTagEffVsEleIdCut);
    gBaseWrongTagEffVsBaseEff->SetPoint(i, baseEffVsEleIdCut, baseWrongTagEffVsEleIdCut);
    gBaseMistagVsBaseEff->SetPoint(i, baseEffVsEleIdCut, baseMistagVsEleIdCut);
    gBaseDilutionVsBaseEff->SetPoint(i, baseEffVsEleIdCut, baseDilutionVsEleIdCut);
    gBasePowerVsBaseEff->SetPoint(i, baseEffVsEleIdCut, basePowerVsEleIdCut);
    
    gBaseWrongTagEffVsBaseRightTagEff->SetPoint(i, baseRightTagEffVsEleIdCut, baseWrongTagEffVsEleIdCut);
  }
  
  autoSavedObject = CreateCanvas("cIsTagRightVsTagTruth", "colz", false, false, false, hIsTagRightVsTagTruth);
  
  autoSavedObject = CreateCanvas("cNRightTagVsEleIdCut", 21, 1, false, false, gNRightTagVsEleIdCut);
  autoSavedObject = CreateCanvas("cNWrongTagVsEleIdCut", 21, 1, false, false,  gNWrongTagVsEleIdCut);
  autoSavedObject = CreateCanvas("cNNoTagVsEleIdCut", 21, 1, false, false,  gNNoTagVsEleIdCut);
  autoSavedObject = CreateCanvas("cBaseEffVsEleIdCut", 21, 1, false, false,  gBaseEffVsEleIdCut);
  autoSavedObject = CreateCanvas("cBaseRightTagEffVsEleIdCut", 21, 1, false, false,  gBaseRightTagEffVsEleIdCut);
  autoSavedObject = CreateCanvas("cBaseWrongTagEffVsEleIdCut", 21, 1, false, false,  gBaseWrongTagEffVsEleIdCut);
  autoSavedObject = CreateCanvas("cBaseMistagVsEleIdCut", 21, 1, false, false,  gBaseMistagVsEleIdCut);
  autoSavedObject = CreateCanvas("cBaseDilutionVsEleIdCut", 21, 1, false, false,  gBaseDilutionVsEleIdCut);
  autoSavedObject = CreateCanvas("cBasePowerVsEleIdCut", 21, 1, false, false,  gBasePowerVsEleIdCut);
  
  autoSavedObject = CreateCanvas("cNRightTagVsBaseEff", 21, 1, false, false,  gNRightTagVsBaseEff);
  autoSavedObject = CreateCanvas("cNWrongTagVsBaseEff", 21, 1, false, false,  gNWrongTagVsBaseEff);
  autoSavedObject = CreateCanvas("cNNoTagVsBaseEff", 21, 1, false, false,  gNNoTagVsBaseEff);
  autoSavedObject = CreateCanvas("cBaseEffVsBaseEff", 21, 1, false, false,  gBaseEffVsBaseEff);
  autoSavedObject = CreateCanvas("cBaseRightTagEffVsBaseEff", 21, 1, false, false,  gBaseRightTagEffVsBaseEff);
  autoSavedObject = CreateCanvas("cBaseWrongTagEffVsBaseEff", 21, 1, false, false,  gBaseWrongTagEffVsBaseEff);
  autoSavedObject = CreateCanvas("cBaseMistagVsBaseEff", 21, 1, false, false,  gBaseMistagVsBaseEff);
  autoSavedObject = CreateCanvas("cBaseDilutionVsBaseEff", 21, 1, false, false,  gBaseDilutionVsBaseEff);
  autoSavedObject = CreateCanvas("cBasePowerVsBaseEff", 21, 1, false, false,  gBasePowerVsBaseEff);

  autoSavedObject = CreateCanvas("cBaseWrongTagEffVsBaseRightTagEff", 21, 1, false, false,  gBaseWrongTagEffVsBaseRightTagEff);
  
  std::cout << std::endl;
  std::cout << "Base efficiency = " << 100 * effBase << "%\n";
  std::cout << "Base mistag = " << 100 * wBase << "%\n";
  std::cout << "Base power = " << 100 * pBase << "%\n";
  
  totalP /= (double)(nRightTag + nWrongTag + nNoTag);
  std::cout << std::endl;
  std::cout << "Per-event-mistag power (not calibrated) = " << 100. * totalP << "% (+" << 100 * (totalP - pBase) / pBase << "%)\n";
  std::cout << std::endl;
  
  // CALIBRATION
  for(int j = 0; j < nBinsCal; j++)
  {
    vWCalc[j] /= (vhMassCalRightTag->at(j)->Integral() + vhMassCalWrongTag->at(j)->Integral());
  }
  
  std::vector<double> vX;
  std::vector<double> vY;
  std::vector<double> vEXL;
  std::vector<double> vEXH;  
  std::vector<double> vEYL;
  std::vector<double> vEYH;

  int minEntries = 0;
  if(isData)
  {
    minEntries = 20;
  }
  int rebinThreshold = 1500;
  
  for(int j = 0; j < nBinsCal; j++)
  {
    std::pair<double, double> calRightTag; // .first = nEvt; .second = sigma(nEvt)
    std::pair<double, double> calWrongTag;
    double wMeas;
    double wMeasErr;
    double wMeasErrLow;
    double wMeasErrHigh;
    
    calRightTag.first = vhMassCalRightTag->at(j)->Integral();
    calWrongTag.first = vhMassCalWrongTag->at(j)->Integral();
    calRightTag.second = sqrt(calRightTag.first);
    calWrongTag.second = sqrt(calWrongTag.first);
    
    if(calRightTag.first <= minEntries && calWrongTag.first <= minEntries )
    {
      continue;
    }
    
    if(isData)
    {
      if(calRightTag.first < rebinThreshold) 
      {
        vhMassCalRightTag->at(j)->Rebin();
      }
      if(calWrongTag.first < rebinThreshold)
      {
        vhMassCalWrongTag->at(j)->Rebin();
      }
      if(calRightTag.first >= minEntries)
      {
        calRightTag = CountEventsWithFit(vhMassCalRightTag->at(j));
        hCalFitStatusRightTagVsMvaScore->SetBinContent(j, fitResult->Status());
        hCalFitChi2NDOFRightTagVsMvaScore->SetBinContent(j, fitResult->Chi2() / fitResult->Ndf());
        int nPars = fitResult->NTotalParameters();
        for(auto iPar = 0; iPar < nPars; iPar++)
        {
          std::string parName = fitResult->ParName(iPar);
          for(unsigned int iHisto = 0; iHisto < vhCalParsRightTagVsMvaScore.size(); iHisto++)
          {
            if(strstr(vhCalParsRightTagVsMvaScore[iHisto]->GetName(), parName.c_str()) != 0)
            {
              vhCalParsRightTagVsMvaScore[iHisto]->SetBinContent(j, fitResult->GetParams()[iPar]);
              vhCalParsRightTagVsMvaScore[iHisto]->SetBinError(j, fitResult->GetErrors()[iPar]);
            }
          }
        }
      }
      if(calWrongTag.first >= minEntries)
      {
        calWrongTag = CountEventsWithFit(vhMassCalWrongTag->at(j));
        hCalFitStatusWrongTagVsMvaScore->SetBinContent(j, fitResult->Status());
        hCalFitChi2NDOFWrongTagVsMvaScore->SetBinContent(j, fitResult->Chi2() / fitResult->Ndf());
        int nPars = fitResult->NTotalParameters();
        for(auto iPar = 0; iPar < nPars; iPar++)
        {
          std::string parName = fitResult->ParName(iPar);
          for(unsigned int iHisto = 0; iHisto < vhCalParsWrongTagVsMvaScore.size(); iHisto++)
          {
            if(strstr(vhCalParsWrongTagVsMvaScore[iHisto]->GetName(), parName.c_str()) != 0)
            {
              vhCalParsWrongTagVsMvaScore[iHisto]->SetBinContent(j, fitResult->GetParams()[iPar]);
              vhCalParsWrongTagVsMvaScore[iHisto]->SetBinError(j, fitResult->GetErrors()[iPar]);
            }
          }
        }
      }
      if(calRightTag.second < 1)
      {
        calRightTag.second = 2 * sqrt(calRightTag.first);
      }
      if(calWrongTag.second < 1)
      {
        calWrongTag.second = 2 * sqrt(calWrongTag.first);
      }
      wMeas = calWrongTag.first / (calWrongTag.first + calRightTag.first);
      wMeasErr = sqrt(pow(calWrongTag.first, 2) * pow(calRightTag.second, 2) +
                      pow(calRightTag.first, 2) * pow(calWrongTag.second, 2)) / 
                      pow(calRightTag.first + calWrongTag.first, 2);
      
      wMeasErr += wMeasErr * systematics;
      
      wMeasErrHigh = wMeasErr;
      wMeasErrLow = wMeasErr;
      if(wMeas + wMeasErrHigh > 1.)
      {
        wMeasErrHigh = 1. - wMeas;
      }
      if(wMeas - wMeasErrLow < 0.) 
      {
        wMeasErrLow = wMeas - 0.;
      }
    }
    else
    {
      wMeas = calWrongTag.first / (calWrongTag.first + calRightTag.first);
      wMeasErrHigh = TEfficiency::AgrestiCoull(calWrongTag.first + calRightTag.first, calWrongTag.first, 0.6827, true) - wMeas;
      wMeasErrLow = wMeas - TEfficiency::AgrestiCoull(calWrongTag.first + calRightTag.first, calWrongTag.first, 0.6827, false);
      wMeasErr = std::max(wMeasErrHigh, wMeasErrLow);
    }
    
    vX.push_back(vWCalc[j]);
    vEXL.push_back(vWCalc[j] - ((double)j * binSizeCal));
    vEXH.push_back(((double)j * binSizeCal + binSizeCal) - vWCalc[j]);
    vY.push_back(wMeas); // measured mistag
    vEYL.push_back(wMeasErrLow);
    vEYH.push_back(wMeasErrHigh);
    
    totalPBinned += (calRightTag.first + calWrongTag.first) * pow(1. - 2. * wMeas, 2);
    // totalPBinnedErr
    std::cout << "BIN " << j << " -- vWCalc = " << vWCalc[j];
    std::cout << " -- nRightTag " << (int)calRightTag.first << " +- " << (int)calRightTag.second << " -- nWrongTag " << (int)calWrongTag.first << " +- " << (int)calWrongTag.second;
    std::cout << " -- wMeas " << wMeas <<" +- " << wMeasErr << std::endl << std::endl;
  }
  
  autoSavedObject = reinterpret_cast<std::vector<TObject*>* >(&vhCalParsRightTagVsMvaScore);
  autoSavedObject = reinterpret_cast<std::vector<TObject*>* >(&vhCalParsWrongTagVsMvaScore);

  
  totalPBinned /= (double)nTot;
  
  std::cout << std::endl;
  std::cout << "Per-event-mistag power (calibrated bins) = " << 100. * totalPBinned << "% +- " << 100. * totalPBinnedErr << " (+" << 100 * (totalPBinned - pBase) / pBase << "%)\n\n\n";
  
  auto *c1 = new TCanvas("c1", "c1", 1000, 1600);
  TPad *pad1 = new TPad("pad1", "", 0.0, 0.3, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "", 0.0, 0.0, 1.0, 0.3);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  gPad->SetGrid();
  
  auto *gCal = new TGraphAsymmErrors(vX.size(), &vX[0], &vY[0], 0, 0, &vEYL[0], &vEYH[0]);
  auto *gCalErr = new TGraphAsymmErrors(vX.size(), &vX[0], &vY[0], &vEXL[0], &vEXH[0], &vEYL[0], &vEYH[0]);
  auto *fCal = new TF1("osElectronCal", "[0]+[1]*x", 0., 1.);
  
  TFitResultPtr fitresultCal = gCal->Fit("osElectronCal","S");
  fCal = gCal->GetFunction("osElectronCal");
  
  double q = fCal->GetParameter(0);
  double m = fCal->GetParameter(1);
  
  std::cout << std::endl;
  std::cout << "q = " << q << " +- " << fCal->GetParError(0) << " [" << abs(q) / fCal->GetParError(0) << " s.d.]\n";
  std::cout << "m = " << m << " +- " << fCal->GetParError(1) << " [" << abs(m - 1) / fCal->GetParError(1) << " s.d.]\n";
  
  std::vector<double> wResY;
  std::vector<double> wResEY;
  std::vector<double> wResEYH;
  std::vector<double> wResEYL;
  std::vector<double> wRatioY;
  std::vector<double> wRatioEY;
  std::vector<double> wRatioEYH;
  std::vector<double> wRatioEYL;
  
  for(unsigned int j = 0; j < vX.size(); j++)
  {
    double dev = vY[j] - fCal->Eval(vX[j]);
    double sigma;
    if(dev >= 0) 
    {
      sigma = vEYL[j];
    }
    else
    {
      sigma = vEYH[j];
    }
    
    wResEYH.push_back(vEYH[j] / sigma);
    wResEYL.push_back(vEYL[j] / sigma);
    wResY.push_back(dev / sigma);
    wResEY.push_back(1.);
  }
  
  auto *gCalRes = new TGraphAsymmErrors(vX.size(), &vX[0], &wResY[0], &vEXL[0], &vEXH[0], &wResEYL[0], &wResEYH[0]);
  
  gCal->SetMarkerStyle(20);
  gCal->SetMarkerSize(1);
  gCal->SetMaximum(1.);
  gCal->SetMinimum(0.);
  gCal->GetXaxis()->SetLimits(0., 1.);
  gCal->SetTitle(("calibration " + process).c_str());
  gCal->GetXaxis()->SetTitle("mistag calc.");
  gCal->GetYaxis()->SetTitle("mistag meas.");
  gCal->Draw("AP");
  gCalErr->Draw("EZ same");
  fCal->Draw("same");
  gPad->Modified(); 
  gPad->Update();
  TPaveStats *st = (TPaveStats*)gCal->FindObject("stats");
  st->SetX1NDC(0.65);
  st->SetX2NDC(0.95);
  st->SetY1NDC(0.15);
  st->SetY2NDC(0.30);
  gPad->Modified(); 
  gPad->Update();
  
  pad2->cd();
  gPad->SetGrid();
  gCalRes->SetMarkerStyle(20);
  gCalRes->SetMarkerSize(1);
  gCalRes->GetXaxis()->SetLimits(0.0, 1.02);
  gCalRes->Draw("APZ");
  gCalRes->SetTitle("");
  gCalRes->GetYaxis()->SetTitle("# s.d.");
  auto *y0_ = new TF1("", "0.", 0., 1.02);
  y0_->SetLineColor(kBlack);
  y0_->Draw("SAME");
  
  if(writeOutput)
  {
    c1->Print(("calibration_" + process + ".pdf").c_str());
  }
  
  auto *c2 = new TCanvas();
  gPad->SetGrid();
  hMva->SetMarkerStyle(20);
  hMva->SetMarkerSize(.75);
  hMva->SetTitle(("dnnDistribution " + process).c_str());
  hMva->GetXaxis()->SetTitle("dnn score (right tag probability)");
  hMva->GetXaxis()->SetNdivisions(10 + 100 * (int)(nMassBins / 10), kFALSE);
  gStyle->SetOptStat(0);
  hMva->Draw("HIST PL");
  if(writeOutput)
  {
    c2->Print(("dnnDistribution_" + process + ".pdf").c_str());
  }
  gStyle->SetOptStat(10);
  
  //FUNCTIONS
  if(writeOutput)
  {
    auto *fo = new TFile(("OSElectronTaggerCalibration" + process + ".root").c_str(), "RECREATE");
    fo->cd();
    fCal->Write();
    fitresultCal->Write();
    fo->Close();
    delete fo;
  }
}



std::pair<double, double> EleMVACalibrationProducer::CountEventsWithFit(TH1 *hist)
{
  std::string name = hist->GetName();
  std::string title = hist->GetTitle();
  // cout<<" ---  now fitting "<<title<<endl;
  std::cout << "EleMVACalibrationProducer::CountEventsWithFit(\"" << name << "\")...\n";
  
  double nEntries = hist->GetEntries();
  
  bool isTot = name == "hMassTot" ? true : false;
  bool lowStat = nEntries <= 0 ? true : false;
  bool highStat = nEntries > 0 ? true : false;
  
  std::cout << "isTot: " << (isTot?"true":"false") << std::endl;
  std::cout << "lowStat: " << (lowStat?"true":"false") << std::endl;
  std::cout << "highStat: " << (highStat?"true":"false") << std::endl;
  
  TRandom3 *r3 = new TRandom3();
  double mean = 5.3663;
  double sigma = 0.015;
  if(process.find("BsJPsiPhi") != std::string::npos)
  {
    mean = 5.3663;
    std::cout << "Process is BsJPsiPhi. mean = " << mean << std::endl;
  }
  if(process.find("BuJPsiK") != std::string::npos)
  {
    mean = 5.2793;
    std::cout << "Process is BuJPsiK. mean = " << mean << std::endl;
  }
  if(process.find("BdJPsiKx") != std::string::npos)
  {
    mean = 5.2796;
    std::cout << "Process is BdJPsiKx. mean = " << mean << std::endl;
  }
  
  TString sgnDef = "[1]*TMath::Gaus(x, [0], [2], true)";
  sgnDef +=       "+[3]*[1]*TMath::Gaus(x, [0], [4], true)";
  
  TString bkgDef = "[5]";
  if(!lowStat)  bkgDef += "+[6]*x";
  if(highStat) bkgDef += "+[7]*TMath::Erfc([8]*(x-[9]))";
  
  bkgDef = "(" + bkgDef + ">= 0 ? " + bkgDef + " : 0 )";
  bkgDef = "(" + bkgDef + ">= 0 ? " + bkgDef + " : 0 )";
  
  TString funcDef = sgnDef + "+" + bkgDef;
  
  std::cout << "Function used for fit: \"" << funcDef << "\"\n";
  
  
  std::cout << "Fit limits: minMass = " << minMass << ", maxMass = " << maxMass << std::endl;
  
  TF1 *func = new TF1("func", funcDef, minMass, maxMass);
  TF1 *sgn = new TF1("sgn", sgnDef, minMass, maxMass);
  TF1 *bkg = new TF1("bkg", bkgDef, minMass, maxMass);
  
  func->SetParName(0, "mean");
  func->SetParName(1, "norm1");
  func->SetParName(2, "sigma1");
  func->SetParName(3, "frac");
  func->SetParName(4, "sigma2");
  
  sgn->SetParName(0, "mean");
  sgn->SetParName(1, "norm1");
  sgn->SetParName(2, "sigma1");
  sgn->SetParName(3, "frac");
  sgn->SetParName(4, "sigma2");
  
  func->SetParName(5, "bkgConst");
  bkg->SetParName(5, "bkgConst");
  
  if(!lowStat)
  {
    func->SetParName(6, "bkgLinSlope");
    bkg->SetParName(6, "bkgLinSlope");
  }

  if(highStat)
  {
    func->SetParName(7, "erfcNorm");
    func->SetParName(8, "erfcWidth");
    func->SetParName(9, "erfcShift");
    
    bkg->SetParName(7, "erfcNorm");
    bkg->SetParName(8, "erfcWidth");
    bkg->SetParName(9, "erfcShift");
  }
  
  
  // Book fit stability histograms on-the-fly
  int nPars = func->GetNpar();
  
  for(auto iPar = 0; iPar < nPars; iPar++)
  {
    std::string parName = func->GetParName(iPar);
    bool found = false;
    for(unsigned int iHisto = 0; iHisto < vhCalParsRightTagVsMvaScore.size(); iHisto++)
    {
      if(strstr(vhCalParsRightTagVsMvaScore[iHisto]->GetName(), parName.c_str()) != 0)
      {
        found = true;
        break;
      }
    }
    if(!found)
    {
//       if(verbose)
        std::cout << "I N F O : Histogram for parameter " << parName << " not found. Creating it...\n";
      std::string histoName = "h" + parName + "RightTagVsMvaScore";
      std::string histoTitle = parName + " right tag vs. Mva score";
      TH1D* hTemp = Create1DHistogram<TH1D>(histoName.c_str(), histoTitle.c_str(), nBinsCal, 0., 1., "mistag calc.", parName.c_str());
      vhCalParsRightTagVsMvaScore.push_back(hTemp);
      
      histoName = "h" + parName + "WrongTagVsMvaScore";
      histoTitle = parName + " wrong tag vs. Mva score";
      hTemp = Create1DHistogram<TH1D>(histoName.c_str(), histoTitle.c_str(), nBinsCal, 0., 1., "mistag calc.", parName.c_str());
      vhCalParsWrongTagVsMvaScore.push_back(hTemp);
    }
  }
  
  //SIGNAL
  double limit = hist->GetEntries() * hist->GetBinWidth(1);

  double rnd1 = r3->Gaus(1.,0.01);
  double rnd3 = r3->Gaus(1.,0.01);
  double rnd2 = r3->Gaus(1.,0.0001);
  double rnd4 = r3->Gaus(1.,0.001);
  
  func->SetParameter(0, mean);
  func->SetParameter(1, limit * rnd1);
  func->SetParameter(3, 1 * rnd3);
  func->SetParameter(2, sigma * rnd2);
  func->SetParameter(4, sigma * rnd4);
  func->SetParLimits(1, 0, 4 * limit);
  func->SetParLimits(3, 0, 5);
  func->SetParLimits(2, 0.0001, 0.5);
  func->SetParLimits(4, 0.0001, 0.5);

  sgn->SetParameter(0, mean);
  sgn->SetParameter(1, limit * rnd1);
  sgn->SetParameter(3, 1 * rnd3);
  sgn->SetParameter(2, sigma * rnd2);
  sgn->SetParameter(4, sigma * rnd4);
  sgn->SetParLimits(1, 0, 4 * limit);
  sgn->SetParLimits(3, 0, 5);
  sgn->SetParLimits(2, 0.0001, 0.5);
  sgn->SetParLimits(4, 0.0001, 0.5);
  
  //BKG    
  TAxis *xaxis = hist->GetXaxis();
  int binx1 = xaxis->FindBin(x1MassSB);
  int binx2 = xaxis->FindBin(x2MassSB);
  double y2 = hist->GetBinContent(binx2);
  double y1 = hist->GetBinContent(binx1);
  double x2 = hist->GetBinCenter(binx2);
  double x1 = hist->GetBinCenter(binx1);
  std::cout << "x1 = " << x1 << ", y1 = " << y1 << std::endl;
  std::cout << "x2 = " << x2 << ", y2 = " << y2 << std::endl;
  double m = (y2 - y1) / (x2 - x1);
  std::cout << "Starting value for m = " << m << std::endl;
  if(m > 0)
  {
    std::cout << "Starting value for m is positive. Setting m = -1. \n";
    m = -1;
  }
  if(nEntries < 100)
  {
    std::cout << "Histogram has < 100 entries. Setting m = -1. \n";
    m = -1;
  }
  
  func->SetParameter(6, m);  
  bkg->SetParameter(6, m);
  
  if(lowStat)
  {
    func->SetParameter(5, nEntries / 50.);
    bkg->SetParameter(5, nEntries / 50.);
  }
  else
  {
    func->SetParameter(5, nEntries / 50. - m * x2);
    bkg->SetParameter(5, nEntries / 50. - m * x2);  
  }

  if(lowStat)
  {
//     func->SetParameter(5, 1);
    func->SetParLimits(5, 0, 1e3);
    
    bkg->SetParLimits(5, 0, 1e3);    
  }
  if(!lowStat)
  {
    func->SetParLimits(6, -1e5, 0);
    
    bkg->SetParLimits(6, -1e5, 0);
  }
  
  if(highStat)
  {
    func->SetParameter(7, hist->GetBinContent(2) / 2.);
    func->SetParameter(8, 10);
    func->SetParameter(9, 5);
    func->SetParLimits(7, 0, hist->GetBinContent(2) * 2.);

    bkg->SetParameter(7, hist->GetBinContent(2) / 2.);
    bkg->SetParameter(8, 10);
    bkg->SetParameter(9, 5);
    bkg->SetParLimits(7, 0, hist->GetBinContent(2) * 2.);
  }
  
  //FIXING PARAMETERS
  if(!isTot)
  {
    std::cout << "Fixing parameters to fit to tot statistics.\n";
    func->FixParameter(0, meanTotMass);
    func->FixParameter(2, sigma1TotMass);
    func->FixParameter(4, sigma2TotMass);
    func->FixParameter(3, fracTotMass);
    
    func->FixParameter(8, erfcWidthTotMass);
    func->FixParameter(9, erfcShiftTotMass);
    
    sgn->FixParameter(0, meanTotMass);
    sgn->FixParameter(2, sigma1TotMass);
    sgn->FixParameter(4, sigma2TotMass);
    sgn->FixParameter(3, fracTotMass);
    
    bkg->FixParameter(8, erfcWidthTotMass);
    bkg->FixParameter(9, erfcShiftTotMass);    
  }
  
  if(process.find("MC") != std::string::npos)
  {
    std::cout << "MC: Fixing more parameters...\n";
    func->FixParameter(5, 0);
    func->FixParameter(6, 0);
    func->FixParameter(7, 0);
    func->FixParameter(8, 0);
    func->FixParameter(9, 0);
    
    bkg->FixParameter(5, 0);
    bkg->FixParameter(6, 0);
    bkg->FixParameter(7, 0);
    bkg->FixParameter(8, 0);
    bkg->FixParameter(9, 0);
  }
  
  func->SetNpx(2000);
  sgn->SetNpx(2000);
  bkg->SetNpx(2000);
  
  std::cout << "Before fit, background values:\n";
  for(auto x = minMass; x < maxMass; x+=hist->GetBinWidth(1))
  {
    std::cout << "(" << x << ", " << (*bkg)(x) << ")" << " ";
  }
  std::cout << std::endl;
  
  auto c5 = new TCanvas();
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(.75);
  fitResult = hist->Fit("func", "LRS");
  int fitstatus = fitResult;
  int covstatus = fitResult->CovMatrixStatus();
  if(fitstatus != 0)
  {
    std::cout << "STATUS of " << name << " --> " << fitstatus << std::endl;
  }
  if(covstatus != 3)
  {
    std::cout << "COV STATUS of " << name << " --> " << covstatus << std::endl;
  }
  TF1 *fit = hist->GetFunction("func");
  if(isTot)
  {
    std::cout << "Tot fit. Saving parameters...\n";
    meanTotMass = fit->GetParameter("mean");
    sigma1TotMass = fit->GetParameter("sigma1");
    sigma2TotMass = fit->GetParameter("sigma2");
    fracTotMass = fit->GetParameter("frac");
    erfcWidthTotMass = fit->GetParameter("erfcWidth");
    erfcShiftTotMass = fit->GetParameter("erfcShift");
    
    std::cout << "meanTotMass = " << meanTotMass << std::endl;
    std::cout << "sigma1TotMass = " << sigma1TotMass << std::endl;
    std::cout << "sigma2TotMass = " << sigma2TotMass << std::endl;
    std::cout << "erfcWidthTotMass = " << erfcWidthTotMass << std::endl;
    std::cout << "erfcShiftTotMass = " << erfcShiftTotMass << std::endl;
  }
  
  hist->Draw("PE");
  hist->SetMinimum(0);
  
  // PLOTTING
  TF1 *f1 = new TF1("f1", "[0]*TMath::Gaus(x, [1], [2], true)", minMass, maxMass);
  TF1 *f2 = new TF1("f2", "[0]*TMath::Gaus(x, [1], [2], true)", minMass, maxMass);
  TF1 *f4;
  TF1 *f5 = new TF1("f5", "[0]+[1]*x", minMass, maxMass);
  
  f1->SetParameters(fit->GetParameter("norm1"), fit->GetParameter("mean"), fit->GetParameter("sigma1"));
  f2->SetParameters(fit->GetParameter("norm1")*fit->GetParameter("frac"), fit->GetParameter("mean"), fit->GetParameter("sigma2"));
  
  if(lowStat)
  {
    f4 = new TF1("f4", "[0]", minMass, maxMass);
    f4->SetParameter(0, fit->GetParameter(5));
  }
  
  if(highStat)
  {
    f4 = new TF1("f4", "[0]+[1]*x+[2]*TMath::Erfc([3]*(x-[4]))", minMass, maxMass);
    f4->SetParameters(fit->GetParameter(5), fit->GetParameter(6),
                      fit->GetParameter(7), fit->GetParameter(8),
                      fit->GetParameter(9));
    f5->SetParameters(fit->GetParameter(5), fit->GetParameter(6));
  }
  
  if(!lowStat && !highStat)
  {
    f4 = new TF1("f4","[0]+[1]*x", minMass, maxMass);
    f4->SetParameters(fit->GetParameter(5), fit->GetParameter(6));
  }
  
  f1->SetLineColor(kBlue);
  f2->SetLineColor(kViolet);
  f4->SetLineColor(kGreen);
  if(highStat)
  {
    f5->SetLineColor(kOrange);
  }
  f1->SetLineStyle(2);
  f2->SetLineStyle(2);
  f4->SetLineStyle(2);
  if(highStat)
  {
    f5->SetLineStyle(2);
  }
  
  f1->Draw("same");
  f2->Draw("same");
  f4->Draw("same");
  if(highStat)
  {
    f5->Draw("same");
  }
  
  if(writeOutput)
  {
    c5->Print((dirPath + "/" + name + ".pdf").c_str());
  }
  
  double nEvt = fit->GetParameter(1);
  nEvt += fit->GetParameter(3)*fit->GetParameter(1);
  std::cout << "nEvt = " << fit->GetParameter(1) << " + (" << fit->GetParameter(3) << " * " << fit->GetParameter(1) << ") = " << nEvt << std::endl;
  nEvt /= hist->GetBinWidth(1);
  std::cout << "After division by bin width: nEvt = " << nEvt << std::endl;
  
  TMatrixDSym cov = fitResult->GetCovarianceMatrix();
//   // This is for N = P1 + P3
//   double errN = sqrt(cov(1,1) + cov(3,3) + 2 * cov(1,3)) / hist->GetBinWidth(1);
  // This is for N = P1 + P1*P3
  double errN = sqrt( (1+fit->GetParameter(3)) * (1+fit->GetParameter(3)) * cov(1,1) +        // |d(N)/d(P1)|^2 * cov(11)
                      fit->GetParameter(1) * fit->GetParameter(1) * cov(3,3) +                // |d(N)/d(P3)|^2 * cov(33)
                      2 * fabs((1+fit->GetParameter(3)) * fit->GetParameter(1)) * cov(1,3) )  // 2*|d(N)/d(P1)*d(N)/d(P3)| * cov(33)
                / hist->GetBinWidth(1);
  std::cout << "errN = " << errN << std::endl;  
  return std::make_pair(nEvt, errN);  
}
