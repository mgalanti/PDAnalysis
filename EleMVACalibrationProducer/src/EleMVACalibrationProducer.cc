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
  getUserParameter("inputFileName", inputFileName);
  getUserParameter("inputTreeName", inputTreeName);
  getUserParameter("mvaMethod", mvaMethod);
  getUserParameter("useTightSelection", useTightSelection);
  getUserParameter("useSyst", useSyst);
  getUserParameter("nBinsCal", nBinsCal);
  getUserParameter("systematics2017", systematics2017);
  getUserParameter("systematics2018", systematics2018);
  getUserParameter("mvaInputPath", mvaInputPath);
  getUserParameter("eleIdWP", eleIdWP);
  getUserParameter("elePtWP", elePtWP);
  getUserParameter("eleDzWP", eleDzWP);
  
  sampleName = inputFileName;
  
  do
  {
    sampleName = sampleName.substr(sampleName.find("/") + 1);
  } 
  while(sampleName.find("/") != std::string::npos);
  
  sampleName = sampleName.substr(0, sampleName.find("."));

  
  std::cout << "I N F O : EleMVACalibrationProducer::beginJob() - Parameters:" << std::endl;
  std::cout << "          inputFileName = " << inputFileName << std::endl;
  std::cout << "          sampleName = " << sampleName << std::endl;
  std::cout << "          mvaMethod = " << mvaMethod << std::endl;
  std::cout << "          mvaInputPath = " << mvaInputPath << std::endl;
  std::cout << "          useTightSelection = " << useTightSelection << std::endl;
  std::cout << "          nBinsCal = " << nBinsCal << std::endl;
  
  
//   auto *inputFile = new TFile(inputFileName.c_str());
//   auto *inputTree = (TTree*)inputFile->Get(inputTreeName.c_str());
//   if(inputFile->IsZombie())
//   {
//     return;
//   }

  
  if(inputFileName.find("Bs") != std::string::npos)
  {
    process = "BsJPsiPhi";
    minMass = 5.20;
    maxMass = 5.65;
    dirPath = "./Bs";
    if(inputFileName.find("DG0") != std::string::npos)
    {
      process += "DG0";
    }
  }
  if(inputFileName.find("Bu") != std::string::npos)
  {
    process = "BuJPsiK";
    minMass = 5.10;
    maxMass = 5.65;
    x1MassSB = 5.44;
    x2MassSB = 5.64;
    dirPath = "./Bu";
  }
  
  if(inputFileName.find("MC") != std::string::npos)
  {
    process = process + "MC";
    dirPath += "MC";
  }
  if(inputFileName.find("Data") != std::string::npos)
  {
    process = process + "Data";
    dirPath += "Data";
    isData = true;
  }
  
  if(inputFileName.find("2017") != std::string::npos)
  {
    process = process + "2017";
    dirPath += "2017";
    if(useSyst)
    {
      systematics = systematics2017;
    }
  }
  
  if(inputFileName.find("2018") != std::string::npos)
  {
    process = process + "2018";
    dirPath += "2018";
    if(useSyst)
    {
      systematics = systematics2018;
    }
  }
  
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
  reader->AddVariable("elePt", &elePt);
  reader->AddVariable("eleEta", &eleEta);
  reader->AddVariable("eleDxy", &eleDxy);
  reader->AddVariable("eleExy", &eleExy);
  reader->AddVariable("eleDz", &eleDz);
  reader->AddVariable("eleEz", &eleEz);
  reader->AddVariable("eleIDNIV2Val", &eleIDNIV2Val);
  reader->AddVariable("eleDRB", &eleDRB);
  reader->AddVariable("elePFIsoScaled", &elePFIsoScaled);
  reader->AddVariable("eleConeCleanPt", &eleConeCleanPt);
  reader->AddVariable("eleConeCleanPtRel", &eleConeCleanPtRel);
  reader->AddVariable("eleConeCleanDR", &eleConeCleanDR);
  reader->AddVariable("eleConeCleanEnergyRatio", &eleConeCleanEnergyRatio);
  reader->AddVariable("eleConeCleanQ", &eleConeCleanQ);
  reader->BookMVA(mvaMethod, mvaInputPath + "TMVAClassification_" + mvaMethod + ".weights.xml");
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
  if(useTightSelection && !tightEvent)
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
  
  //TAGGING
  double mvaValue = reader->EvaluateMVA(mvaMethod);
  hMva->Fill(mvaValue, evtWeight);
  evtW[0] = 1 - mvaValue;
  totalP += pow(1. - 2. * evtW[0], 2) * evtWeight;
  
  // FIXME: Uncomment these as soon as the electron charge is in the ntuples
//   int evtTag = -1 * eleCharge;
//   bool isTagRight = TMath::Sign(1, BidGen) == evtTag;
//   if(isTagRight != tagTruth)
//   {
//     std::cout << "W A R N I N G! EleMVACalibrationProducer::analyze(): isTagRight != tagTruth!\n";
//   }
  
  // FIXME: Temporary solution
  bool isTagRight = tagTruth;
  
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
  
  if(isTagRight)
  {
    hMvaRightTag->Fill(mvaValue, evtWeight);
    hMassRightTag->Fill(BMass, evtWeight);
  }
  else
  {
    hMvaWrongTag->Fill(mvaValue, evtWeight);
    hMassWrongTag->Fill(BMass, evtWeight);
  }
  
  hNGenB->Fill(nGenB);
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
      }
      if(calWrongTag.first >= minEntries)
      {
        calWrongTag = CountEventsWithFit(vhMassCalWrongTag->at(j));
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
    c1->Print(("calibration" + process + ".pdf").c_str());
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
    auto *fo = new TFile(("OSMuonTaggerCalibration" + process + ".root").c_str(), "RECREATE");
    fo->cd();
    fCal->Write();
    fitresultCal->Write();
    fo->Close();
    delete fo;
  }
}



std::pair<double, double> EleMVACalibrationProducer::CountEventsWithFit(TH1 *hist)
{
  std::string title = hist->GetTitle();
  // cout<<" ---  now fitting "<<title<<endl;
  
  bool isTot = title == "hMassTot" ? true : false;
  bool lowStat = hist->GetEntries() <= 250 ? true : false;
  bool highStat = hist->GetEntries() > 1500 ? true : false;
  
  TRandom3 *r3 = new TRandom3();
  double mean = 5.3663;
  double sigma = 0.015;
  if(process.find("BsJPsiPhi") != std::string::npos)
  {
    mean = 5.3663;
  }
  if(process.find("BuJPsiK") != std::string::npos)
  {
    mean = 5.2793;
  }
  if(process.find("BdJPsiKx") != std::string::npos)
  {
    mean = 5.2796;
  }
  
  TString sgnDef = "[1]*TMath::Gaus(x, [0], [2], true)";
  sgnDef +=       "+[3]*TMath::Gaus(x, [0], [4], true)";
  
  TString bkgDef = "[5]";
  if(!lowStat)  bkgDef += "+[6]*x";
  if(highStat) bkgDef += "+[7]*TMath::Erfc([8]*(x-[9]))";
  
  bkgDef = "(" + bkgDef + ">= 0 ? " + bkgDef + " : 0 )";
  
  TString funcDef = sgnDef + "+" + bkgDef;
  
  TF1 *func = new TF1("func", funcDef, minMass, maxMass);
  
  func->SetParName(0, "mean");
  func->SetParName(1, "A1");
  func->SetParName(2, "sigma1");
  func->SetParName(3, "A2");
  func->SetParName(4, "sigma2");
  
  //SIGNAL
  double limit = hist->GetEntries() * hist->GetBinWidth(1);
  
  func->SetParameter(0, mean);
  func->SetParameter(1, limit / 2 * r3->Gaus(1.,0.01));
  func->SetParameter(3, limit / 2 * r3->Gaus(1.,0.01));
  func->SetParameter(2, sigma * r3->Gaus(1.,0.01));
  func->SetParameter(4, sigma * r3->Gaus(1.,0.01));
  func->SetParLimits(1, 0, 2 * limit);
  func->SetParLimits(3, 0, 2 * limit);
  func->SetParLimits(2, 0.001, 0.1);
  func->SetParLimits(4, 0.001, 0.1);
  
  //BKG    
  TAxis *xaxis = hist->GetXaxis();
  int binx1 = xaxis->FindBin(x1MassSB);
  int binx2 = xaxis->FindBin(x2MassSB);
  double y2 = hist->GetBinContent(binx2);
  double y1 = hist->GetBinContent(binx1);
  double x2 = hist->GetBinCenter(binx2);
  double x1 = hist->GetBinCenter(binx1);
  double m = (y2 - y1) / (x2 - x1);
  if(m > 0)
  {
    m = -1;
  }
  
  func->SetParameter(5, 10);
  func->SetParameter(6, m);
  
  if(lowStat)
  {
    func->SetParameter(5, 1);
    func->SetParLimits(5, 0, 1e3);
  }
  
  if(highStat)
  {
    func->SetParameter(7, hist->GetBinContent(2) / 2);
    func->SetParameter(8, 10);
    func->SetParameter(9, 5);
    func->SetParLimits(7, 0, hist->GetBinContent(2) * 1.5);
  }
  
  //FIXING PARAMETERS
  if(!isTot)
  {
    func->FixParameter(0, meanTotMass);
    func->FixParameter(2, sigma1TotMass);
    func->FixParameter(4, sigma2TotMass);
  }
  
  if(process.find("MC") != std::string::npos)
  {
    func->FixParameter(5, 0);
    func->FixParameter(6, 0);
    func->FixParameter(7, 0);
    func->FixParameter(8, 0);
    func->FixParameter(9, 0);
  }
  
  func->SetNpx(2000);
  
  auto c5 = new TCanvas();
  hist->SetMarkerStyle(20);
  hist->SetMarkerSize(.75);
  TFitResultPtr r = hist->Fit("func", "LRSQ");
  int fitstatus = r;
  int covstatus = r->CovMatrixStatus();
  if(fitstatus != 0)
  {
    std::cout << "STATUS of " << title << " --> " << fitstatus << std::endl;
  }
  if(covstatus != 3)
  {
    std::cout << "COV STATUS of " << title << " --> " << covstatus << std::endl;
  }
  TF1 *fit = hist->GetFunction("func");
  if(isTot)
  {
    meanTotMass = fit->GetParameter("mean");
    sigma1TotMass = fit->GetParameter("sigma1");
    sigma2TotMass = fit->GetParameter("sigma2");
  }
  
  hist->Draw("PE");
  hist->SetMinimum(0);
  
  // PLOTTING
  TF1 *f1 = new TF1("f1", "[0]*TMath::Gaus(x, [1], [2], true)", minMass, maxMass);
  TF1 *f2 = new TF1("f2", "[0]*TMath::Gaus(x, [1], [2], true)", minMass, maxMass);
  TF1 *f4;
  TF1 *f5 = new TF1("f5", "[0]+[1]*x", minMass, maxMass);
  
  f1->SetParameters(fit->GetParameter("A1"), fit->GetParameter("mean"), fit->GetParameter("sigma1"));
  f2->SetParameters(fit->GetParameter("A2"), fit->GetParameter("mean"), fit->GetParameter("sigma2"));
  
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
    c5->Print((dirPath + "/" + title + ".pdf").c_str());
  }
  
  double nEvt = fit->GetParameter(1);
  nEvt += fit->GetParameter(3);
  nEvt /= hist->GetBinWidth(1);
  
  TMatrixDSym cov = r->GetCovarianceMatrix();
  double errN = sqrt(cov(1,1) + cov(3,3) + 2 * cov(1,3)) / hist->GetBinWidth(1);
  
  return std::make_pair(nEvt, errN);  
}
