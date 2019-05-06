#include <iostream>
#include <sstream>
#include <string>
#include <bitset>
#include <math.h>

#include "ElectronEfficiencyAnalyzer.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "Math/LorentzVector.h"

using namespace std;



ElectronEfficiencyAnalyzer::~ElectronEfficiencyAnalyzer() {
}



void ElectronEfficiencyAnalyzer::beginJob() {

  MGBaseAnalyzer::beginJob();
  MGSelector::SetSelectionString(evtSelection);


  // user parameters are retrieved as strings by using their names;
  // numeric parameters ( int, float or whatever ) can be directly set
  // by passing the corresponding variable,
  // e.g. getUserParameter( "name", x )

  getUserParameter( "verbose", verbose );

  // Get the ID flags for the electrons using the enums...
  HZZV1IDEleBit = PDEnumString::mvaEleID_Spring16_HZZ_V1_wpLoose;   
  HZZV2IDEleBit = PDEnumString::mvaEleID_Fall17_iso_V2_wpHZZ;   
  MVANILIDEleBit = PDEnumString::mvaEleID_Fall17_noIso_V2_wpLoose;

  // user parameters are set as names associated to a string, 
  // default values can be set in the analyzer class contructor
    
//   std::cout << "ElectronEfficiencyAnalyzer::beginJob(): sampleName   = " << sampleName << std::endl;
//   std::cout << "                                        evtSelection = " << evtSelection << std::endl;

  setUserParameter("verbose", "f");
//   std::cout << "use_info is " << use_info << std::endl;
  
  std::pair<int, std::string> elePair;
  
  elePair = std::make_pair(0, "");
  genEleAllMap.insert(elePair);
  
  elePair = std::make_pair(0, "");
  genEleFromBMap.insert(elePair);
  elePair = std::make_pair(1, "FromB");
  genEleFromBMap.insert(elePair);
//   elePair = std::make_pair(2, "NotFromB");
//   genEleFromBMap.insert(elePair);
  
  elePair = std::make_pair(0, "");
  genEleMatchedMap.insert(elePair);
  elePair = std::make_pair(10, "Matched");
  genEleMatchedMap.insert(elePair);
//   elePair = std::make_pair(20, "NotMatched");
//   genEleMatchedMap.insert(elePair);
  
  elePair = std::make_pair(0, "");
  genEleBsChargeCorrMap.insert(elePair);
  elePair = std::make_pair(100, "NegCorr");
  genEleBsChargeCorrMap.insert(elePair);
  elePair = std::make_pair(200, "ZeroCorr");
  genEleBsChargeCorrMap.insert(elePair);
  elePair = std::make_pair(300, "PosCorr");
  genEleBsChargeCorrMap.insert(elePair);
    
  for(auto itAll = genEleAllMap.begin(); itAll != genEleAllMap.end(); itAll++)
  {
    int iAll = itAll->first;
    std::string strAll = itAll->second;
    for(auto itFromB = genEleFromBMap.begin(); itFromB != genEleFromBMap.end(); itFromB++)
    {
      int iFromB = itFromB->first;
      std::string strFromB = itFromB->second;
      
      int iGraph = iAll + iFromB;
      std::string strGraph = strAll + strFromB;
      if(iGraph == 0)
        strGraph = "All";
      elePair = std::make_pair(iGraph, strGraph);
      genEleGraphMap.insert(elePair);
      
      for(auto itMatched = genEleMatchedMap.begin(); itMatched != genEleMatchedMap.end(); itMatched++)
      {
        int iMatched = itMatched->first;
        std::string strMatched = itMatched->second;
        
        for(auto itBsChargeCorr = genEleBsChargeCorrMap.begin(); itBsChargeCorr != genEleBsChargeCorrMap.end(); itBsChargeCorr++)
        {
          int iBsChargeCorr = itBsChargeCorr->first;
          std::string strBsChargeCorr = itBsChargeCorr->second;
          
          int iHist = iAll + iFromB + iMatched + iBsChargeCorr;
          std::string strHist = strAll + strFromB + strMatched + strBsChargeCorr;
          if(iHist == 0)
            strHist = "All";
          elePair = std::make_pair(iHist, strHist);
          genEleHistMap.insert(elePair);
        }
      }
    }
  }
  elePair = std::make_pair(0, "");
  eleAllMap.insert(elePair);
  
  elePair = std::make_pair(0, "");
  eleGenMap.insert(elePair);
  elePair = std::make_pair(1, "True");
  eleGenMap.insert(elePair);
//   elePair = std::make_pair(2, "Fake");
//   eleGenMap.insert(elePair);
  elePair = std::make_pair(3, "TrueFromB");
  eleGenMap.insert(elePair);
//   elePair = std::make_pair(4, "TrueNotFromB");
//   eleGenMap.insert(elePair);
  
  elePair = std::make_pair(0, "");
  elePVMap.insert(elePair);
  elePair = std::make_pair(10, "SamePV");
  elePVMap.insert(elePair);
//   elePair = std::make_pair(20, "OtherPV");
//   elePVMap.insert(elePair);
  
  elePair = std::make_pair(0, "");
  eleSelMap.insert(elePair);
  elePair = std::make_pair(100, "HZZV1Sel");
  eleSelMap.insert(elePair);
//   elePair = std::make_pair(200, "NotHZZV1Sel");
//   eleSelMap.insert(elePair);
  elePair = std::make_pair(300, "HZZV2Sel");
  eleSelMap.insert(elePair);
//   elePair = std::make_pair(400, "NotHZZV2Sel");
//   eleSelMap.insert(elePair);
  elePair = std::make_pair(500, "MVANILSel");
  eleSelMap.insert(elePair);
//   elePair = std::make_pair(600, "NotMVANILSel");
//   eleSelMap.insert(elePair);
  
  elePair = std::make_pair(0, "");
  eleBsChargeCorrMap.insert(elePair);
  elePair = std::make_pair(1000, "NegCorr");
  eleBsChargeCorrMap.insert(elePair);
  elePair = std::make_pair(2000, "ZeroCorr");
  eleBsChargeCorrMap.insert(elePair);
  elePair = std::make_pair(3000, "PosCorr");
  eleBsChargeCorrMap.insert(elePair);

  for(auto itAll = eleAllMap.begin(); itAll != eleAllMap.end(); itAll++)
  {
    int iAll = itAll->first;
    std::string strAll = itAll->second;
    for(auto itGen = eleGenMap.begin(); itGen != eleGenMap.end(); itGen++)
    {
      int iGen = itGen->first;
      std::string strGen = itGen->second;
      for(auto itPV = elePVMap.begin(); itPV != elePVMap.end(); itPV++)
      {
        int iPV = itPV->first;
        std::string strPV = itPV->second;
                
        for(auto itSel = eleSelMap.begin(); itSel != eleSelMap.end(); itSel++)
        {
          int iSel = itSel->first;
          std::string strSel = itSel->second;
          
          for(auto itBsChargeCorr = eleBsChargeCorrMap.begin(); itBsChargeCorr != eleBsChargeCorrMap.end(); itBsChargeCorr++)
          {
            int iBsChargeCorr = itBsChargeCorr->first;
            std::string strBsChargeCorr = itBsChargeCorr->second;
            
            int iHist = iAll + iGen + iPV + iSel + iBsChargeCorr;
            std::string strHist = strAll + strGen + strPV + strSel + strBsChargeCorr;
            if(iHist == 0)
              strHist = "All";
            elePair = std::make_pair(iHist, strHist);
            eleHistMap.insert(elePair);
            
            int iGraph = iAll + iGen + iPV + iSel + iBsChargeCorr;
            std::string strGraph = strAll + strGen + strPV + strSel + strBsChargeCorr;
            if(iGraph == 0)
              strGraph = "All";
            elePair = std::make_pair(iGraph, strGraph);
            eleGraphMap.insert(elePair);
          }
        }
      }
    }
  }
  
  return;
}


void ElectronEfficiencyAnalyzer::book() 
{

  // putting "autoSavedObject" in front of the histo creation 
  // it's automatically marked for saving on file; the option 
  // is uneffective when not using the full utility
  const char* nGenEleAxisName = "N. of gen. electrons";
  const char* genElectronAxisName = "N. gen. electrons";
  
  const char* recoEfficiencyAxisName = "Reco efficiency";
  const char* efficiencyAxisName = "Efficiency";

  const char* nEleAxisName = "N. of electrons";
  const char* ptAxisName = "p_{T} [GeV]";
  const char* phiAxisName = "#phi";
  const char* etaAxisName = "#eta";
  const char* energyAxisName = "E [GeV]";
  const char* chargeAxisName = "Charge";
//   const char* muonAxisName = "N. muons";
  const char* electronAxisName = "N. electrons";
  const char* eventAxisName = "N. of events";
//   const char* isolationAxisName = "Isolation";
//   const char* effectiveAreaAxisName = "Effective area";
  const char* idAxisName = "ID";
//   const char* consChaAxisName = "Charge consistency";
//   const char* ebeeGapAxisName = "electron in EB-EE gap";
  const char* dBAxisName = "dB [cm]";
  const char* bsChargeCorrAxisName = "Charge correlation w.r.t B_{s}";
  const char* recoChargeCorrAxisName = "Charge correlation w.r.t reco e";
  const char* genChargeCorrAxisName = "Charge correlation w.r.t gen e";
  const char* HZZV1IDAxisName = "HZZ V1 ID flag";
  const char* HZZV2IDAxisName = "HZZ V2 ID flag";
  const char* MVANILIDAxisName = "MVA 17V2 Non-Iso Loose ID flag";
  const char* MVAOutputAxisName = "MVA output";
  const char* MVACatAxisName = "MVA category";
  const char* distanceAxisName = "distance [cm]";
  
  vPtBinEdges.push_back(0);
  vPtBinEdges.push_back(2);
  vPtBinEdges.push_back(5);
  vPtBinEdges.push_back(10);
  vPtBinEdges.push_back(20);
  vPtBinEdges.push_back(40);
  
  vMVACatBinEdges.push_back(-0.5);
  vMVACatBinEdges.push_back(0.5);
  vMVACatBinEdges.push_back(1.5);
  vMVACatBinEdges.push_back(2.5);
  vMVACatBinEdges.push_back(3.5);
  vMVACatBinEdges.push_back(4.5);
  vMVACatBinEdges.push_back(5.5);
  vMVACatBinEdges.push_back(6.5);
  vMVACatBinEdges.push_back(7.5);
  
  std::cout << "Booking Gen histograms...\n";
  std::string histFullName;
  std::string histFullTitle;
  TH1D* hist;
  TH2D* hist2D;
  std::vector<TH1D*>* vHist;
  
  for(auto itGenEleHistMap = genEleHistMap.begin(); itGenEleHistMap != genEleHistMap.end(); itGenEleHistMap++)
  {
    int iHist = itGenEleHistMap->first;
    std::string nameHist = itGenEleHistMap->second;
    std::cout << "iHist = " << iHist << ", nameHist = " << nameHist << std::endl;
 
    // cd's to the directory named nameHist in the output root file
    autoSavedObject = gDirectory->mkdir(("Hist_gen_" + nameHist).c_str());    
    
    histFullName = "hN" + nameHist + "GenEle";
    histFullTitle = "N. of " + nameHist + " Gen Electrons";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 11, -0.5, 10.5, nGenEleAxisName, eventAxisName);
    hNGenEle.insert(std::make_pair(iHist, hist));
    autoSavedObject = hNGenEle[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "GenEleCharge";
    histFullTitle = nameHist + " Gen Electrons charge";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 3, -1.5, 1.5, chargeAxisName, genElectronAxisName);
    hGenEleCharge.insert(std::make_pair(iHist, hist));
    autoSavedObject = hGenEleCharge[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";

    histFullName = "h" + nameHist + "GenElePt";
    histFullTitle = nameHist + " Gen Electrons p_{T}";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 100, 0.0, 100.0, ptAxisName, genElectronAxisName);
    hGenElePt.insert(std::make_pair(iHist, hist));
    autoSavedObject = hGenElePt[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "GenElePhi";
    histFullTitle = nameHist + " Gen Electrons #phi";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 100, -constants::PI, constants::PI, phiAxisName, genElectronAxisName);
    hGenElePhi.insert(std::make_pair(iHist, hist));
    autoSavedObject = hGenElePhi[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "GenEleEta";
    histFullTitle = nameHist + " Gen Electrons #eta";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 100, -3.0, 3.0, etaAxisName, genElectronAxisName);
    hGenEleEta.insert(std::make_pair(iHist, hist));
    autoSavedObject = hGenEleEta[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "GenElePtVsEta";
    histFullTitle = nameHist + " Gen Electrons p_{T} vs. #eta";
    hist2D = Create2DHistogram<TH2D>(histFullName.c_str(), histFullTitle.c_str(), 100, -3, 3, 100, 0, 100, etaAxisName, ptAxisName);
    hGenElePtVsEta.insert(std::make_pair(iHist, hist2D));
    autoSavedObject = hGenElePtVsEta[iHist];
    std::cout << hist2D->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "GenEleE";
    histFullTitle = nameHist + " Gen Electrons energy";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 200, 0., 200.0, energyAxisName, genElectronAxisName);
    hGenEleE.insert(std::make_pair(iHist, hist));
    autoSavedObject = hGenEleE[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";

    histFullName = "h" + nameHist + "GenEleDb";
    histFullTitle = nameHist + " Gen Electrons dB w.r.t. Gen PV";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 100, 0.0, 2.0, dBAxisName, genElectronAxisName);
    hGenEleDb.insert(std::make_pair(iHist, hist));
    autoSavedObject = hGenEleDb[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";

    histFullName = "h" + nameHist + "GenEleBsChargeCorr";
    histFullTitle = nameHist + " Gen Electrons charge correlation w.r.t. Bs";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 5, -2.5, 2.5, bsChargeCorrAxisName, genElectronAxisName);
    hGenEleBsChargeCorr.insert(std::make_pair(iHist, hist));
    autoSavedObject = hGenEleBsChargeCorr[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";

    histFullName = "h" + nameHist + "GenEleRecoChargeCorr";
    histFullTitle = nameHist + " Gen Electrons charge correlation w.r.t. matched reco electron";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 3, -1.5, 1.5, recoChargeCorrAxisName, genElectronAxisName);
    hGenEleRecoChargeCorr.insert(std::make_pair(iHist, hist));
    autoSavedObject = hGenEleRecoChargeCorr[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";

    // cd's back to the root directory of the output root file
    autoSavedObject = gDirectory->GetDirectory("..");
  }

  std::cout << "Booking Gen graphs...\n";
  std::string graphFullName;
  std::string graphFullTitle;
  TGraphAsymmErrors* graph;

  for(auto itGenEleGraphMap = genEleGraphMap.begin(); itGenEleGraphMap != genEleGraphMap.end(); itGenEleGraphMap++)
  {
    int iGraph = itGenEleGraphMap->first;
    std::string nameGraph = itGenEleGraphMap->second;
    std::cout << "iGraph = " << iGraph << ", nameGraph = " << nameGraph << std::endl;
    
    graphFullName = "g" + nameGraph + "GenEleRecoEffVsCharge";
    graphFullTitle = "Reco efficiency of " + nameGraph + " Gen Electrons vs. charge";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), chargeAxisName, recoEfficiencyAxisName);
    gGenEleRecoEffVsCharge.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gGenEleRecoEffVsCharge[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";
  
    graphFullName = "g" + nameGraph + "GenEleRecoEffVsPt";
    graphFullTitle = "Reco efficiency of " + nameGraph + " Gen Electrons vs. gen p_{T}";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), ptAxisName, recoEfficiencyAxisName);
    gGenEleRecoEffVsPt.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gGenEleRecoEffVsPt[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";
  
    graphFullName = "g" + nameGraph + "GenEleRecoEffVsPhi";
    graphFullTitle = "Reco efficiency of " + nameGraph + " Gen Electrons vs. gen #phi";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), phiAxisName, recoEfficiencyAxisName);
    gGenEleRecoEffVsPhi.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gGenEleRecoEffVsPhi[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";
    
    graphFullName = "g" + nameGraph + "GenEleRecoEffVsEta";
    graphFullTitle = "Reco efficiency of " + nameGraph + " Gen Electrons vs. gen #eta";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), etaAxisName, recoEfficiencyAxisName);
    gGenEleRecoEffVsEta.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gGenEleRecoEffVsEta[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";
    
    graphFullName = "g" + nameGraph + "GenEleRecoEffVsE";
    graphFullTitle = "Reco efficiency of " + nameGraph + " Gen Electrons vs. gen energy";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), energyAxisName, recoEfficiencyAxisName);
    gGenEleRecoEffVsE.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gGenEleRecoEffVsE[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";

    graphFullName = "g" + nameGraph + "GenEleRecoEffVsDb";
    graphFullTitle = "Reco efficiency of " + nameGraph + " Gen Electrons vs. gen dB";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), dBAxisName, recoEfficiencyAxisName);
    gGenEleRecoEffVsDb.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gGenEleRecoEffVsDb[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";

    graphFullName = "g" + nameGraph + "GenEleRecoEffVsBsChargeCorr";
    graphFullTitle = "Reco efficiency of " + nameGraph + " Gen Electrons vs. charge correlation w.r.t. B_{s}";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), bsChargeCorrAxisName, recoEfficiencyAxisName);
    gGenEleRecoEffVsBsChargeCorr.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gGenEleRecoEffVsBsChargeCorr[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";

    graphFullName = "g" + nameGraph + "GenEleRecoEffVsRecoChargeCorr";
    graphFullTitle = "Reco efficiency of " + nameGraph + " Gen Electrons vs. charge correlation w.r.t. reco e";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), recoChargeCorrAxisName, recoEfficiencyAxisName);
    gGenEleRecoEffVsRecoChargeCorr.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gGenEleRecoEffVsRecoChargeCorr[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";
  }
  
  std::cout << "Booking Reco histograms...\n";

  for(auto itEleHistMap = eleHistMap.begin(); itEleHistMap != eleHistMap.end(); itEleHistMap++)
  {
    int iHist = itEleHistMap->first;
    std::string nameHist = itEleHistMap->second;
    std::cout << "iHist = " << iHist << ", nameHist = " << nameHist << std::endl;

    // cd's to the directory named nameHist in the output root file
    autoSavedObject = gDirectory->mkdir(("Hist_" + nameHist).c_str());    
    
    histFullName = "hN" + nameHist + "Ele";
    histFullTitle = "N. of " + nameHist + " electrons";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 11, -0.5, 10.5, nEleAxisName, eventAxisName);
    hNEle.insert(std::make_pair(iHist, hist));
    autoSavedObject = hNEle[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "EleCharge";
    histFullTitle = nameHist + " Electrons charge";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 3, -1.5, 1.5, chargeAxisName, electronAxisName);
    hEleCharge.insert(std::make_pair(iHist, hist));
    autoSavedObject = hEleCharge[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";

    histFullName = "h" + nameHist + "ElePt";
    histFullTitle = nameHist + " Electrons p_{T}";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 100, 0.0, 100.0, ptAxisName, electronAxisName);
    hElePt.insert(std::make_pair(iHist, hist));
    autoSavedObject = hElePt[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "ElePhi";
    histFullTitle = nameHist + " Electrons #phi";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 100, -constants::PI, constants::PI, phiAxisName, electronAxisName);
    hElePhi.insert(std::make_pair(iHist, hist));
    autoSavedObject = hElePhi[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "EleEta";
    histFullTitle = nameHist + " Electrons #eta";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 100, -3.0, 3.0, etaAxisName, electronAxisName);
    hEleEta.insert(std::make_pair(iHist, hist));
    autoSavedObject = hEleEta[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "ElePtVsEta";
    histFullTitle = nameHist + " Electrons p_{T} vs. #eta";
    hist2D = Create2DHistogram<TH2D>(histFullName.c_str(), histFullTitle.c_str(), 100, -3, 3, 100, 0, 100, etaAxisName, ptAxisName);
    hElePtVsEta.insert(std::make_pair(iHist, hist2D));
    autoSavedObject = hElePtVsEta[iHist];
    std::cout << hist2D->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "EleE";
    histFullTitle = nameHist + " Electrons energy";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 200, 0., 200.0, energyAxisName, electronAxisName);
    hEleE.insert(std::make_pair(iHist, hist));
    autoSavedObject = hEleE[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "EleIDs";
    histFullTitle = nameHist + " Electrons IDs";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 1000, 0.0, 4294967295, idAxisName, electronAxisName);
    hEleIDs.insert(std::make_pair(iHist, hist));
    autoSavedObject = hEleIDs[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "EleDb";
    histFullTitle = nameHist + " Electrons dB w.r.t PV";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 100, 0.0, 2.0, dBAxisName, electronAxisName);
    hEleDb.insert(std::make_pair(iHist, hist));
    autoSavedObject = hEleDb[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "EleBsChargeCorr";
    histFullTitle = nameHist + " Electrons charge correlation w.r.t B_{s}";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 5, -2.5, 2.5, bsChargeCorrAxisName, electronAxisName);
    hEleBsChargeCorr.insert(std::make_pair(iHist, hist));
    autoSavedObject = hEleBsChargeCorr[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "EleGenChargeCorr";
    histFullTitle = nameHist + " Electrons charge correlation w.r.t gen e";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 3, -1.5, 1.5, genChargeCorrAxisName, electronAxisName);
    hEleGenChargeCorr.insert(std::make_pair(iHist, hist));
    autoSavedObject = hEleGenChargeCorr[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "ElePassesHZZV1ID";
    histFullTitle = nameHist + " Electrons passing the HZZ V1 ID selection";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 2, -0.5, 1.5, HZZV1IDAxisName, electronAxisName);
    hElePassesHZZV1ID.insert(std::make_pair(iHist, hist));
    autoSavedObject = hElePassesHZZV1ID[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "ElePassesHZZV2ID";
    histFullTitle = nameHist + " Electrons passing the HZZ V2 ID selection";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 2, -0.5, 1.5, HZZV2IDAxisName, electronAxisName);
    hElePassesHZZV2ID.insert(std::make_pair(iHist, hist));
    autoSavedObject = hElePassesHZZV2ID[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "ElePassesMVANILID";
    histFullTitle = nameHist + " Electrons passing the MVA 17V2 Non-Iso loose ID selection";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 2, -0.5, 1.5, MVANILIDAxisName, electronAxisName);
    hElePassesMVANILID.insert(std::make_pair(iHist, hist));
    autoSavedObject = hElePassesMVANILID[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "ElePVZDistanceFromBestPV";
    histFullTitle = "Distance along z of " + nameHist + " Electrons PV from best B PV";
    hist = Create1DHistogram<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 100, -20., 20., distanceAxisName, electronAxisName);
    hElePVZDistanceFromBestPV.insert(std::make_pair(iHist, hist));
    autoSavedObject = hElePVZDistanceFromBestPV[iHist];
    std::cout << hist->ClassName() << "* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "EleEta";
    histFullTitle = nameHist + " Electrons #eta";
    vHist = CreateVectorOf1DHistograms<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 100, -3.0, 3.0, etaAxisName, electronAxisName, "pT", vPtBinEdges);
    vhEleEtaVsPt.insert(std::make_pair(iHist, vHist));
    autoSavedObject = reinterpret_cast<std::vector<TObject*>* >(vhEleEtaVsPt[iHist]);
    std::cout << "std::vector<" << hist->ClassName() << "*>* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "EleHZZMVAOutput";
    histFullTitle = nameHist + " Electrons HZZ MVA output";
    vHist = CreateVectorOf1DHistograms<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 100, -1.0, 1.0, MVAOutputAxisName, electronAxisName, "pT", vPtBinEdges);
    vhEleHZZMVAOutputVsPt.insert(std::make_pair(iHist, vHist));
    autoSavedObject = reinterpret_cast<std::vector<TObject*>* >(vhEleHZZMVAOutputVsPt[iHist]);
    std::cout << "std::vector<" << hist->ClassName() << "*>* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "EleHZZMVAOutput";
    histFullTitle = nameHist + " Electrons HZZ MVA output";
    vHist = CreateVectorOf1DHistograms<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 100, -1.0, 1.0, MVAOutputAxisName, electronAxisName, "HZZ MVA category", vMVACatBinEdges);
    vhEleHZZMVAOutputVsHZZMVACat.insert(std::make_pair(iHist, vHist));
    autoSavedObject = reinterpret_cast<std::vector<TObject*>* >(vhEleHZZMVAOutputVsHZZMVACat[iHist]);
    std::cout << "std::vector<" << hist->ClassName() << "*>* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "EleHZZMVACat";
    histFullTitle = nameHist + " Electrons HZZ MVA category";
    vHist = CreateVectorOf1DHistograms<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 8, -0.5, 7.5, MVACatAxisName, electronAxisName, "pT", vPtBinEdges);
    vhEleHZZMVACatVsPt.insert(std::make_pair(iHist, vHist));
    autoSavedObject = reinterpret_cast<std::vector<TObject*>* >(vhEleHZZMVACatVsPt[iHist]);
    std::cout << "std::vector<" << hist->ClassName() << "*>* " << histFullName << " booked.\n";
    
    histFullName = "h" + nameHist + "EleNIMVAOutput";
    histFullTitle = nameHist + " Electrons Fall17V2 non-iso MVA output";
    vHist = CreateVectorOf1DHistograms<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 100, -2.0, 2.0, MVAOutputAxisName, electronAxisName, "pT", vPtBinEdges);
    vhEleNIMVAOutputVsPt.insert(std::make_pair(iHist, vHist));
    autoSavedObject = reinterpret_cast<std::vector<TObject*>* >(vhEleNIMVAOutputVsPt[iHist]);
    std::cout << "std::vector<" << hist->ClassName() << "*>* " << histFullName << " booked.\n";

    histFullName = "h" + nameHist + "EleNIMVAOutput";
    histFullTitle = nameHist + " Electrons Fall17V2 non-iso MVA output";
    vHist = CreateVectorOf1DHistograms<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 100, -2.0, 2.0, MVAOutputAxisName, electronAxisName, "Fall17V2 non-iso MVA category", vMVACatBinEdges);
    vhEleNIMVAOutputVsNIMVACat.insert(std::make_pair(iHist, vHist));
    autoSavedObject = reinterpret_cast<std::vector<TObject*>* >(vhEleNIMVAOutputVsNIMVACat[iHist]);
    std::cout << "std::vector<" << hist->ClassName() << "*>* " << histFullName << " booked.\n";

    histFullName = "h" + nameHist + "EleNIMVACat";
    histFullTitle = nameHist + " Electrons Fall17V2 non-iso MVA category";
    vHist = CreateVectorOf1DHistograms<TH1D>(histFullName.c_str(), histFullTitle.c_str(), 8, -0.5, 7.5, MVACatAxisName, electronAxisName, "pT", vPtBinEdges);
    vhEleNIMVACatVsPt.insert(std::make_pair(iHist, vHist));
    autoSavedObject = reinterpret_cast<std::vector<TObject*>* >(vhEleNIMVACatVsPt[iHist]);
    std::cout << "std::vector<" << hist->ClassName() << "*>* " << histFullName << " booked.\n";
    
    // cd's back to the root directory of the output root file
    autoSavedObject = gDirectory->GetDirectory("..");
  }

  std::cout << "Booking Reco graphs...\n";

  for(auto itEleGraphMap = eleGraphMap.begin(); itEleGraphMap != eleGraphMap.end(); itEleGraphMap++)
  {
    int iGraph = itEleGraphMap->first;
    std::string nameGraph = itEleGraphMap->second;
    
    int nDigits = 1;
    if (iGraph != 0)
      nDigits = (int)pow(10,floor(log10(iGraph)));

    std::string nameDenominator = "All";
    if(iGraph != 0)
      nameDenominator = eleGraphMap[iGraph%nDigits];
    std::cout << "iGraph = " << iGraph << ", nameGraph = " << nameGraph << ", nameDenominator = " << nameDenominator << std::endl;
    
    graphFullName = "g" + nameGraph + "EleEffVsCharge";
    graphFullTitle = "Efficiency of " + nameGraph + " over " + nameDenominator + " Electrons vs. charge";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), chargeAxisName, efficiencyAxisName);
    gEleEffVsCharge.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gEleEffVsCharge[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";
  
    graphFullName = "g" + nameGraph + "EleEffVsPt";
    graphFullTitle = "Efficiency of " + nameGraph + " over " + nameDenominator + " Electrons vs. reco p_{T}";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), ptAxisName, efficiencyAxisName);
    gEleEffVsPt.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gEleEffVsPt[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";
  
    graphFullName = "g" + nameGraph + "EleEffVsPhi";
    graphFullTitle = "Efficiency of " + nameGraph + " over " + nameDenominator + " Electrons vs. reco #phi";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), phiAxisName, efficiencyAxisName);
    gEleEffVsPhi.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gEleEffVsPhi[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";
    
    graphFullName = "g" + nameGraph + "EleEffVsEta";
    graphFullTitle = "Efficiency of " + nameGraph + " over " + nameDenominator + " Electrons vs. reco #eta";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), etaAxisName, efficiencyAxisName);
    gEleEffVsEta.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gEleEffVsEta[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";
    
    graphFullName = "g" + nameGraph + "EleEffVsE";
    graphFullTitle = "Efficiency of " + nameGraph + " over " + nameDenominator + " Electrons vs. reco energy";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), energyAxisName, efficiencyAxisName);
    gEleEffVsE.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gEleEffVsE[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";

    graphFullName = "g" + nameGraph + "EleEffVsDb";
    graphFullTitle = "Efficiency of " + nameGraph + " over " + nameDenominator + " Electrons vs. reco dB";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), dBAxisName, efficiencyAxisName);
    gEleEffVsDb.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gEleEffVsDb[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";

    graphFullName = "g" + nameGraph + "EleEffVsBsChargeCorr";
    graphFullTitle = "Efficiency of " + nameGraph + " over " + nameDenominator + " Electrons vs. charge correlation w.r.t. B_{s}";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), bsChargeCorrAxisName, efficiencyAxisName);
    gEleEffVsBsChargeCorr.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gEleEffVsBsChargeCorr[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";

    graphFullName = "g" + nameGraph + "EleEffVsGenChargeCorr";
    graphFullTitle = "Efficiency of " + nameGraph + " over " + nameDenominator + " Electrons vs. charge correlation w.r.t. gen e";
    graph = CreateGraph<TGraphAsymmErrors>(graphFullName.c_str(), graphFullTitle.c_str(), genChargeCorrAxisName, efficiencyAxisName);
    gEleEffVsGenChargeCorr.insert(std::make_pair(iGraph, graph));
    autoSavedObject = gEleEffVsGenChargeCorr[iGraph];
    std::cout << graph->ClassName() << "* " << graphFullName << " booked.\n";
  }
  
  autoSavedObject = gEleEffVsPtTest = CreateGraph<TGraphAsymmErrors>("gEleEffVsPtTest", "Efficiency of the HZZ V1 selection vs. pt for SamePV electrons", ptAxisName, "Efficiency");

  return;
}


void ElectronEfficiencyAnalyzer::reset() {
// automatic reset
  autoReset();
  return;
}


bool ElectronEfficiencyAnalyzer::analyze( int entry, int event_file, int event_tot ) {

  if ( verbose ) {
    cout << " +++++++++++++++++++++++++++ " << endl;
    cout << "entry: "
         << entry << " " << event_file << " " << event_tot << endl;
    cout << "run: " <<   runNumber << " , "
         << "evt: " << eventNumber << endl;
  }
  else {
//    if ( !( event_file % 10000 ) || !( event_tot % 10000 ) )
    if ( !( event_tot % 10000 ) && event_tot )
      cout << event_file << " " << event_tot << endl;
  }
  
  std::vector<std::pair<int,int> > selectedObjects;
  std::vector<std::pair<int,int> > selectedObjectsLoose;
  std::vector<std::pair<int,int> > selectedObjectsTight;
  
  bool evtSelected = SelectBsToJPsiPhiEvent(evtSelection.c_str(), selectedObjects);
  
  bool evtSelectedLoose = SelectBsToJPsiPhiEvent("eleTagLooseV0", selectedObjectsLoose);  
  bool evtSelectedTight = SelectBsToJPsiPhiEvent("eleTagTightV0", selectedObjectsTight);
  
//   if(evtSelectedLoose)
//     std::cout << "Event " << eventNumber << " is" << (evtSelectedLoose?" ":" NOT ") << "selected by the loose selection!\n";

  int iSelObject = 0;
  int iBestPV = -1;
  int iBestBs = -1;
  for (auto itSelObjects = selectedObjects.begin(); itSelObjects != selectedObjects.end(); itSelObjects++)
  {
//     std::cout << "Event selection: selected object #" << iSelObject << ": type = " << itSelObjects->first << ", index = " << itSelObjects->second << std::endl;
    if(itSelObjects->first == PDEnumString::recPV)
      iBestPV = itSelObjects->second;
    if(itSelObjects->first == PDEnumString::recSvt)
      iBestBs = itSelObjects->second;
    iSelObject++;
  }
  
  int iGenBs = -1;
  if (iBestBs > -1)
    iGenBs = GetClosestGenNoLL(iBestBs);
  int idGenBs = 0;
  if(iGenBs > 0)
  {
    idGenBs = genId->at(iGenBs);
  }
  
  std::cout << "iGenBs = " << iGenBs << ", idGenBs = " << idGenBs << std::endl;

  int iSelObjectLoose = 0;
  for (auto itSelObjects = selectedObjectsLoose.begin(); itSelObjects != selectedObjectsLoose.end(); itSelObjects++)
  {
//     std::cout << "Loose selection: selected object #" << iSelObjectLoose << ": type = " << itSelObjects->first << ", index = " << itSelObjects->second << std::endl;
    iSelObjectLoose++;
  }

  int iSelObjectTight = 0;
//   if(evtSelectedTight)
//     std::cout << "Event " << eventNumber << " is" << (evtSelectedTight?" ":" NOT ") << "selected by the tight selection!\n";
  for (auto itSelObjects = selectedObjectsTight.begin(); itSelObjects != selectedObjectsTight.end(); itSelObjects++)
  {
//     std::cout << "Tight selection: selected object #" << iSelObjectTight << ": type = " << itSelObjects->first << ", index = " << itSelObjects->second << std::endl;
    iSelObjectTight++;
  }
  
  std::vector<int> allGenElectrons = GetAllGenElectrons();
//   unsigned int iGenEle = 0;
  int nGenElectrons = allGenElectrons.size();
//   for(auto itGenEle = allGenElectrons.begin(); itGenEle != allGenElectrons.end(); itGenEle++)
//   {
//     std::cout << "Gen electron n. " << iGenEle << " has iGen = " << *itGenEle << ", genId = " << genId->at(*itGenEle) << std::endl;
//     int iBMot = RecursiveLookForMotherIds(*itGenEle, listBMesonsAndBaryons, false);
//     if(iBMot > -1)
//     {
//       std::cout << "    The electrons comes from a B meson or baryon!!!\n";
//       std::cout << "        iBMot = " << iBMot << ", idBMot = " << genId->at(iBMot) << std::endl;
//     }
//     PrintMotherChain(*itGenEle);
// //     RecursivePrintMothers(*itGenEle);
//     std::cout << "\n";
//     iGenEle++;
//   }
  
//   std::cout << "Matching of gen electrons from B to reco electrons:\n";
//   std::vector<int> allGenElectronsFromB = GetAllGenElectronsFromB();
//   std::map<int, int> mGenEleFromBMatch;
//   for(unsigned int i = 0; i < allGenElectronsFromB.size(); i++)
//   {
//     int iGenEle = allGenElectronsFromB[i];
//     int iMatchedEle = GetClosestRecoElectron(genPt->at(iGenEle), genEta->at(iGenEle), genPhi->at(iGenEle));
//     mGenEleFromBMatch.insert(std::make_pair(iGenEle, iMatchedEle));
//     std::cout << "Electron with iGen = " << iGenEle << " is matched to reco electron " << iMatchedEle << std::endl; 
//   }
  
  

  // Main event selection
  if(!evtSelected)
    return false;
  
//   std::cout << "iBestPV = " << iBestPV << std::endl;

  // Gen plots
  int iGenElectron;
  float ptGenEle;
  float phiGenEle;
  float etaGenEle;
  float energyGenEle;
  int chargeGenEle;
  float dBGenEle;
  int bsChargeCorrGenEle;
  int recoChargeCorrGenEle;
  std::map<int, int> genEleN;
  for(auto itHistMap = genEleHistMap.begin(); itHistMap != genEleHistMap.end(); itHistMap++)
  {
    std::pair<int, int> genEleNPair = std::make_pair(itHistMap->first, 0);
    genEleN.insert(genEleNPair);
  }
  
  for(int i = 0; i < nGenElectrons; i++)
  {
    iGenElectron = allGenElectrons[i];

    std::vector<int> genEleFromBType;
    genEleFromBType.push_back(0);
    std::vector<int> genEleMatchedType;
    genEleMatchedType.push_back(0);
    std::vector<int> genEleBsChargeCorrType;
    genEleBsChargeCorrType.push_back(0);
    
    
    ptGenEle = genPt->at(iGenElectron);
    phiGenEle = genPhi->at(iGenElectron);
    etaGenEle = genEta->at(iGenElectron);
    energyGenEle = genE->at(iGenElectron);
    chargeGenEle = genCharge->at(iGenElectron);
    // FIXME: fill this variable with the correct value...
    dBGenEle = 0;
    bsChargeCorrGenEle = GetGenLepBsChargeCorrelation(iGenElectron, iGenBs);
    
    
    int iBMot = RecursiveLookForMotherIds(iGenElectron, listBMesonsAndBaryons, false);
    if(iBMot > -1)
    {
      genEleFromBType.push_back(1);
    }
//     else
//     {
//       genEleFromBType.push_back(2);
//     }
    
    if(iBMot > 1 && bsChargeCorrGenEle < 0)
    {
      std::cout << "POSSIBLE PATHOLOGIC CASE: electron from B that has the wrong charge correlation w.r.t. the Bs\n";
      std::cout << "iGenElectron = " << iGenElectron << std::endl;
      std::cout << "Mothers:\n";
      PrintMotherChain(iGenElectron);
      std::cout << std::endl;
    }
    
    int iMatchedEle = GetClosestRecoElectron(genPt->at(iGenElectron), genEta->at(iGenElectron), genPhi->at(iGenElectron));
    recoChargeCorrGenEle = 0;
    if(iMatchedEle > -1)
    {
      genEleMatchedType.push_back(10);
      recoChargeCorrGenEle = genCharge->at(iGenElectron) * eleCharge->at(iMatchedEle);
    }
//     else
//     {
//       genEleMatchedType.push_back(20);
//     }
    
    if(bsChargeCorrGenEle == -1)
    {
      genEleBsChargeCorrType.push_back(100);
    }
    else if(bsChargeCorrGenEle == 0)
    {
      genEleBsChargeCorrType.push_back(200);
    }
    else if(bsChargeCorrGenEle == 1)
    {
      genEleBsChargeCorrType.push_back(300);
    }
    else
    {
      std::cout << "E R R O R ! ElectronEfficiencyAnalyzer::analyze(...): Gen ele charge correlation value " << bsChargeCorrGenEle << " not supported!\n";
      std::cout << "            Exiting...\n";
      exit(1);
    }

    for(auto itFromB = genEleFromBType.begin(); itFromB != genEleFromBType.end(); itFromB++)
    {
      int iFromB = *itFromB;
      for(auto itMatched = genEleMatchedType.begin(); itMatched != genEleMatchedType.end(); itMatched++)
      {
        int iMatched = *itMatched;
        for(auto itBsChargeCorr = genEleBsChargeCorrType.begin(); itBsChargeCorr != genEleBsChargeCorrType.end(); itBsChargeCorr++)
        {
          int iBsChargeCorr = *itBsChargeCorr;
          int iHist = iFromB + iMatched + iBsChargeCorr;
          
          genEleN[iHist]++;
          hGenEleCharge[iHist]        ->Fill(chargeGenEle);   
          hGenElePt[iHist]            ->Fill(ptGenEle);   
          hGenElePhi[iHist]           ->Fill(phiGenEle);
          hGenEleEta[iHist]           ->Fill(etaGenEle);
          hGenElePtVsEta[iHist]       ->Fill(etaGenEle, ptGenEle);  
          hGenEleE[iHist]             ->Fill(energyGenEle);   
          hGenEleDb[iHist]            ->Fill(dBGenEle);   
          hGenEleBsChargeCorr[iHist]  ->Fill(bsChargeCorrGenEle);
          hGenEleRecoChargeCorr[iHist]->Fill(recoChargeCorrGenEle);
        }
      }
    }
  }
  for(auto itGenEleType = genEleHistMap.begin(); itGenEleType != genEleHistMap.end(); itGenEleType++)
  {
    hNGenEle[itGenEleType->first]->Fill(genEleN[itGenEleType->first]);
  }
  
  // Reco plots
  int iElectron;
  float ptEle;
  float phiEle;
  float etaEle;
  float energyEle;
  int chargeEle;
//   float chaIsoEle;
//   float neuIsoEle;
//   float phoIsoEle;
//   float pchIsoEle;
//   float aEffEle;
  int idsEle;
//   bool consChaEle;
//   bool ebeeGapEle;
  float dBEle;
  int bsChargeCorrEle;
  int genChargeCorrEle;
  bool HZZV1IDEle;
  bool HZZV2IDEle;
  bool MVANILIDEle;
  int gsfPVtxEle;
  float HZZMVAOutputEle = -9999;
  int HZZMVACatEle = -9999;
  float NIMVAOutputEle = -9999;
  int NIMVACatEle = -9999;
  float IMVAOutputEle = -9999;
  int IMVACatEle = -9999;
  std::map<int, int> eleN;
  for(auto itHistMap = eleHistMap.begin(); itHistMap != eleHistMap.end(); itHistMap++)
  {
    std::pair<int, int> eleNPair = std::make_pair(itHistMap->first, 0);
    eleN.insert(eleNPair);
  }
//   int eleN[eleHistMap.size()] = {0};
//   int eleHistType = 0;
//   int eleGenType = 0;
//   int elePVType = 0;
//   int eleSelType = 0;

  for(iElectron = 0; iElectron < nElectrons; ++iElectron)
  {
    std::vector<int> eleGenType;
    eleGenType.push_back(0);
    std::vector<int> elePVType;
    elePVType.push_back(0);
    std::vector<int> eleSelType;
    eleSelType.push_back(0);
    std::vector<int> eleBsChargeCorrType;
    eleBsChargeCorrType.push_back(0);
  
//     std::cout << "Electron " << iElectron << " is associated to PV " << eleGsfPVtx->at(iElectron) << std::endl;
    ptEle = elePt->at(iElectron);
    phiEle = elePhi->at(iElectron);
    etaEle = eleEta->at(iElectron);
    chargeEle = eleCharge->at(iElectron);
    energyEle = eleE->at(iElectron);
//     chaIsoEle = eleChaIso->at(iElectron);
//     neuIsoEle = eleNeuIso->at(iElectron);
//     phoIsoEle = elePhoIso->at(iElectron);
//     pchIsoEle = elePCHIso->at(iElectron);
//     aEffEle = eleAEff->at(iElectron);
    idsEle = eleIDs->at(iElectron);
//     consChaEle = eleConsCha->at(iElectron);
//     ebeeGapEle = eleEBEEGap->at(iElectron);
    dBEle = eleDb->at(iElectron);
    
    int iMatchedGenP = GetClosestGen(ptEle, etaEle, phiEle);
    int idMatchedGenP = -9999;
    if (iMatchedGenP >= 0)
      idMatchedGenP = genId->at(iMatchedGenP);
    
    int iMatchedGenPNoLL = GetClosestGenNoLL(ptEle, etaEle, phiEle);
    int idMatchedGenPNoLL = -9999;
    if (iMatchedGenPNoLL >= 0)
      idMatchedGenPNoLL = genId->at(iMatchedGenPNoLL);
    // MG: By default I use the LL version of MC matching.
    //     However, one may want to know when there is a doubt 
    //     that an electron is real or fake. 
    //     Let's print it out...
    if (iMatchedGenP != iMatchedGenPNoLL && (abs(idMatchedGenP) == 11 || abs(idMatchedGenPNoLL) == 11))
    {
      std::cout << "W A R N I N G! LL and NoLL version of the electron MC matching differ on whether the electron is a fake or not!\n";
      std::cout << "               runNumber = " << runNumber << ", eventNumber = " << eventNumber << ", iElectron = " << iElectron << std::endl; 
      float dRMatched = -9999;
      float dPtMatched = -9999;
      if(iMatchedGenP >= 0)
      {
        dRMatched  = deltaR(etaEle, phiEle, genEta->at(iMatchedGenP), genPhi->at(iMatchedGenP));
        dPtMatched = fabs(genPt->at(iMatchedGenP) - ptEle)/genPt->at(iMatchedGenP);
      }
      std::cout << "               Gen particle matched to electron has i = " << iMatchedGenP << ", id = " << idMatchedGenP << std::endl;
      std::cout << "                   dRMatched = " << dRMatched << ", dPtMatched = " << dPtMatched << std::endl;

      float dRMatchedNoLL = -9999;
      float dPtMatchedNoLL = -9999;
      if(iMatchedGenPNoLL >= 0)
      {
        dRMatchedNoLL  = deltaR(etaEle, phiEle, genEta->at(iMatchedGenPNoLL), genPhi->at(iMatchedGenPNoLL));
        dPtMatchedNoLL = fabs(genPt->at(iMatchedGenPNoLL) - ptEle)/genPt->at(iMatchedGenPNoLL);
      }
      std::cout << "               Gen particle matched to electron has i = " << iMatchedGenPNoLL << ", id = " << idMatchedGenPNoLL << " (no LL requirement)" << std::endl;
      std::cout << "                   dRMatchedNoLL = " << dRMatchedNoLL << ", dPtMatchedNoLL = " << dPtMatchedNoLL << std::endl << std::endl;
    }
    std::cout << "I N F O : Reco electron " << iElectron << " is matched to gen particle at index " << iMatchedGenP << " with Id " << idMatchedGenP << std::endl;

    genChargeCorrEle = 0;
    if(abs(idMatchedGenP) == 11)
    {
      eleGenType.push_back(1);
      genChargeCorrEle = chargeEle * genCharge->at(iMatchedGenP);
      // Sub-case: true electron coming (or not) from a B hadron
      if(RecursiveLookForMotherIds(iMatchedGenP, listBMesonsAndBaryons, false) > -1)
      {
        eleGenType.push_back(3);
      }
//       else
//       {
//         eleGenType.push_back(4);
//       }
    }
//     else
//     {
//       eleGenType.push_back(2);
//     }
    
    gsfPVtxEle = eleGsfPVtx->at(iElectron) - nPVertices; 
//     std::cout << "gsfPVtxEle = " << gsfPVtxEle << std::endl;
    
    double eleDzWrtBestPV = dZ(iElectron, iBestPV);
//     std::cout << "eleDzWrtBestPV = " << eleDzWrtBestPV << std::endl;
    
    
//     if(gsfPVtxEle == iBestPV)
    if(fabs(eleDzWrtBestPV) < 1.)
    {
      elePVType.push_back(10);
    }
//     else
//     {
//       elePVType.push_back(20);
//     }
    
    //     std::cout << "This event has " << useObjType->size() << " userInfo entries\n";
    
    for (int iUserInfo = 0; iUserInfo < nUserInfo; iUserInfo++)
    {
//       std::cout << "Looping on userInfo object #" << iUserInfo << std::endl;
      if(useObjType->at(iUserInfo) == 3 && useObjIndex->at(iUserInfo) == iElectron)
      {
//         std::cout << "Found information for this electron...\n";
        if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Spring16HZZV1Values)
        {
//           std::cout << "Found HZZMVAOutputEle value!\n";
          HZZMVAOutputEle = useInfoValue->at(iUserInfo);
        }
        if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Spring16HZZV1Categories)
        {
//           std::cout << "Found HZZMVACatEle value!\n";
          HZZMVACatEle = useInfoValue->at(iUserInfo);
        }
        if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17NoIsoV2Values)
        {
//           std::cout << "Found NIMVAOutputEle value!\n";
          NIMVAOutputEle = useInfoValue->at(iUserInfo);
        }
        if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17NoIsoV2Categories)
        {
//           std::cout << "Found NIMVACatEle value!\n";
          NIMVACatEle = useInfoValue->at(iUserInfo);
        }
      }
    }
    
//     std::cout << "DEBUG: idsEle         = " << std::bitset<32>(idsEle) << std::endl;
//     std::cout << "DEBUG: eleIDBit       = " << std::bitset<32>(eleIDBit) << std::endl;
//     std::cout << "DEBUG: HZZIDEleBit    = " << std::bitset<32>(HZZIDEleBit) << std::endl;
//     std::cout << "DEBUG: MVANILIDEleBit = " << std::bitset<32>(MVANILIDEleBit) << std::endl;
    
    HZZV1IDEle = idsEle & HZZV1IDEleBit;
    HZZV2IDEle = idsEle & HZZV2IDEleBit;
    MVANILIDEle = idsEle & MVANILIDEleBit;
//     bool elePassesID;
//     if (eleID == "All")
//       elePassesID = true;
//     else
//       elePassesID = idsEle & eleIDBit;
//     if (eleSelection == "HZZ")
//       eleID = HZZIDEle;
//     else if (eleSelection == "MVANIL")
//       eleID = MVANILIDEle;
//     else if (eleSelection == "Open")
//       eleID = true;
//     else
//       eleID = false;

//     std::cout << "DEBUG: elePassesID = " << elePassesID << std::endl;
//     std::cout << "DEBUG: HZZIDEle    = " << HZZIDEle << std::endl;
//     std::cout << "DEBUG: MVANILIDEle = " << MVANILIDEle << std::endl;
    
//     if(elePassesID)
//     {
//       eleSelType.push_back(100);
//     }
//     else
//     {
//       eleSelType.push_back(200);
//     }

    if(HZZV1IDEle)
    {
      eleSelType.push_back(100);
    }
//     else
//     {
//       eleSelType.push_back(200);
//     }
    if(HZZV2IDEle)
    {
      eleSelType.push_back(300);
    }
//     else
//     {
//       eleSelType.push_back(400);
//     }
    if(MVANILIDEle)
    {
      eleSelType.push_back(500);
    }
//     else
//     {
//       eleSelType.push_back(600);
//     }

    // Finds the charge correlation between the gen particle matched to the electron and the Bs
    // Returns +/-1 only if the matched gen particle is a lepton (11, 13)
    bsChargeCorrEle = GetGenLepBsChargeCorrelation(iMatchedGenP, iGenBs);
    // If the above method returns 0, try with the reco electron charge...
    if(!bsChargeCorrEle)
    {
      bsChargeCorrEle = GetBsChargeCorrelation(chargeEle, iGenBs);
      // In this case, correlation value goes to the +/-2 bins
      bsChargeCorrEle*=2;
    }
    if(bsChargeCorrEle > 0)
    {
      eleBsChargeCorrType.push_back(3000);
    }
    else if(bsChargeCorrEle == 0)
    {
      eleBsChargeCorrType.push_back(2000);
    }
    else // if(bsChargeCorrEle < 0)
    {
      eleBsChargeCorrType.push_back(1000);
    }
//     else
//     {
//       std::cout << "E R R O R ! ElectronEfficiencyAnalyzer::analyze(...): Ele charge correlation value " << bsChargeCorrEle << " not supported!\n";
//       std::cout << "            Exiting...\n";
//       exit(1);
//     }
    
    for(auto itGen = eleGenType.begin(); itGen != eleGenType.end(); itGen++)
    {
      int iGen = *itGen;
      for(auto itPV = elePVType.begin(); itPV != elePVType.end(); itPV++)
      {
        int iPV = *itPV;
        for(auto itSel = eleSelType.begin(); itSel != eleSelType.end(); itSel++)
        {
          int iSel = *itSel;
          for(auto itBsChargeCorr = eleBsChargeCorrType.begin(); itBsChargeCorr != eleBsChargeCorrType.end(); itBsChargeCorr++)
          {
            int iBsChargeCorr = *itBsChargeCorr;
            int iHist = iGen + iPV + iSel + iBsChargeCorr;
            eleN[iHist]++;
            hEleCharge[iHist]        ->Fill(chargeEle);   
            hElePt[iHist]            ->Fill(ptEle);   
            hElePhi[iHist]           ->Fill(phiEle);
            hEleEta[iHist]           ->Fill(etaEle);
            hElePtVsEta[iHist]       ->Fill(etaEle, ptEle);   
            hEleE[iHist]             ->Fill(energyEle);
            hEleIDs[iHist]           ->Fill(idsEle);
            hEleDb[iHist]            ->Fill(dBEle);
            hEleBsChargeCorr[iHist]  ->Fill(bsChargeCorrEle);
            hEleGenChargeCorr[iHist] ->Fill(genChargeCorrEle);
            hElePassesHZZV1ID[iHist] ->Fill(HZZV1IDEle);
            hElePassesHZZV2ID[iHist] ->Fill(HZZV2IDEle);
            hElePassesMVANILID[iHist]->Fill(MVANILIDEle);
            hElePVZDistanceFromBestPV[iHist]->Fill(pvtZ->at(gsfPVtxEle) - pvtZ->at(iBestPV));
            Fill(vhEleEtaVsPt[iHist], etaEle, ptEle, vPtBinEdges);
            Fill(vhEleHZZMVAOutputVsPt[iHist], HZZMVAOutputEle, ptEle, vPtBinEdges);
            Fill(vhEleHZZMVAOutputVsHZZMVACat[iHist], HZZMVAOutputEle, HZZMVACatEle, vMVACatBinEdges);
            Fill(vhEleHZZMVACatVsPt[iHist], HZZMVACatEle, ptEle, vPtBinEdges);
            Fill(vhEleNIMVAOutputVsPt[iHist], NIMVAOutputEle, ptEle, vPtBinEdges);
            Fill(vhEleNIMVAOutputVsNIMVACat[iHist], NIMVAOutputEle, NIMVACatEle, vMVACatBinEdges);
            Fill(vhEleNIMVACatVsPt[iHist], NIMVACatEle, ptEle, vPtBinEdges);    
          }
        }
      }
    }
  }
  for(auto itEleType = eleHistMap.begin(); itEleType != eleHistMap.end(); itEleType++)
  {
    hNEle[itEleType->first]->Fill(eleN[itEleType->first]);
  }
  return true;
}


void ElectronEfficiencyAnalyzer::endJob() 
{
  
  std::string cFullName;
  TCanvas* cTemp;
  std::vector<TCanvas*> * vcTemp;
  
  for(auto itGenEleHistMap = genEleHistMap.begin(); itGenEleHistMap != genEleHistMap.end(); itGenEleHistMap++)
  {
    int iCanvas = itGenEleHistMap->first;
    std::string nameCanvas = itGenEleHistMap->second;
    
    std::cout << "Producing canvases for iCanvas = " << iCanvas << " and nameCanvas = " << nameCanvas << std::endl;

    autoSavedObject = gDirectory->mkdir(("Canvas_gen_" + nameCanvas).c_str());
    
    cFullName = "cN" + nameCanvas + "GenEle";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, true, hNGenEle[iCanvas]);
    cNGenEle.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cNGenEle[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "GenEleCharge";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, true, hGenEleCharge[iCanvas]);
    cGenEleCharge.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cGenEleCharge[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "GenElePt";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, true, hGenElePt[iCanvas]);
    cGenElePt.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cGenElePt[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "GenElePhi";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, false, hGenElePhi[iCanvas]);
    cGenElePhi.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cGenElePhi[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "GenEleEta";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, false, hGenEleEta[iCanvas]);
    cGenEleEta.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cGenEleEta[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "GenElePtVsEta";
    cTemp = CreateCanvas(cFullName.c_str(), "colz", false, false, true, hGenElePtVsEta[iCanvas]);
    cGenElePtVsEta.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cGenElePtVsEta[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "GenEleE";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, true, hGenEleE[iCanvas]);
    cGenEleE.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cGenEleE[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "GenEleDb";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, true, hGenEleDb[iCanvas]);
    cGenEleDb.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cGenEleDb[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "GenEleBsChargeCorr";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, false, hGenEleBsChargeCorr[iCanvas]);
    cGenEleBsChargeCorr.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cGenEleBsChargeCorr[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "GenEleRecoChargeCorr";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, false, hGenEleRecoChargeCorr[iCanvas]);
    cGenEleRecoChargeCorr.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cGenEleRecoChargeCorr[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    autoSavedObject = gDirectory->GetDirectory("..");
  }

  for(auto itGenEleGraphMap = genEleGraphMap.begin(); itGenEleGraphMap != genEleGraphMap.end(); itGenEleGraphMap++)
  {
    int iCanvas = itGenEleGraphMap->first;
    std::string nameCanvas = itGenEleGraphMap->second;

    autoSavedObject = gDirectory->mkdir(("Canvas_eff_gen_" + nameCanvas).c_str());
    
    std::cout << "Calculating efficiencies and producing canvases for iCanvas = " << iCanvas << " and nameCanvas = " << nameCanvas << std::endl;
    FeldmanCousinsDivide(hGenEleCharge[iCanvas+10], hGenEleCharge[iCanvas], gGenEleRecoEffVsCharge[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "GenEleRecoEffVsCharge").c_str(), 21, 1, false, false, gGenEleRecoEffVsCharge[iCanvas]);
  
    FeldmanCousinsDivide(hGenElePt[iCanvas+10], hGenElePt[iCanvas], gGenEleRecoEffVsPt[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "GenEleRecoEffVsPt").c_str(), 21, 1, false, false, gGenEleRecoEffVsPt[iCanvas]);
    
    FeldmanCousinsDivide(hGenElePhi[iCanvas+10], hGenElePhi[iCanvas], gGenEleRecoEffVsPhi[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "GenEleRecoEffVsPhi").c_str(), 21, 1, false, false, gGenEleRecoEffVsPhi[iCanvas]);
    
    FeldmanCousinsDivide(hGenEleEta[iCanvas+10], hGenEleEta[iCanvas], gGenEleRecoEffVsEta[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "GenEleRecoEffVsEta").c_str(), 21, 1, false, false, gGenEleRecoEffVsEta[iCanvas]);
    
    FeldmanCousinsDivide(hGenEleE[iCanvas+10], hGenEleE[iCanvas], gGenEleRecoEffVsE[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "GenEleRecoEffVsE").c_str(), 21, 1, false, false, gGenEleRecoEffVsE[iCanvas]);
    
    FeldmanCousinsDivide(hGenEleDb[iCanvas+10], hGenEleDb[iCanvas], gGenEleRecoEffVsDb[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "GenEleRecoEffVsDb").c_str(), 21, 1, false, false, gGenEleRecoEffVsDb[iCanvas]);
    
    FeldmanCousinsDivide(hGenEleBsChargeCorr[iCanvas+10], hGenEleBsChargeCorr[iCanvas], gGenEleRecoEffVsBsChargeCorr[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "GenEleRecoEffVsBsChargeCorr").c_str(), 21, 1, false, false, gGenEleRecoEffVsBsChargeCorr[iCanvas]);
    
    FeldmanCousinsDivide(hGenEleRecoChargeCorr[iCanvas+10], hGenEleRecoChargeCorr[iCanvas], gGenEleRecoEffVsRecoChargeCorr[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "GenEleRecoEffVsRecoChargeCorr").c_str(), 21, 1, false, false, gGenEleRecoEffVsRecoChargeCorr[iCanvas]);
    
    autoSavedObject = gDirectory->GetDirectory("..");
  }
  
  for(auto itEleHistMap = eleHistMap.begin(); itEleHistMap != eleHistMap.end(); itEleHistMap++)
  {
    int iCanvas = itEleHistMap->first;
    std::string nameCanvas = itEleHistMap->second;
    
    std::cout << "Producing canvases for iCanvas = " << iCanvas << " and nameCanvas = " << nameCanvas << std::endl;

    autoSavedObject = gDirectory->mkdir(("Canvas_" + nameCanvas).c_str());
    
    cFullName = "cN" + nameCanvas + "Ele";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, true, hNEle[iCanvas]);
    cNEle.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cNEle[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "EleCharge";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, true, hEleCharge[iCanvas]);
    cEleCharge.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cEleCharge[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "ElePt";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, true, hElePt[iCanvas]);
    cElePt.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cElePt[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "ElePhi";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, false, hElePhi[iCanvas]);
    cElePhi.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cElePhi[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "EleEta";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, false, hEleEta[iCanvas]);
    cEleEta.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cEleEta[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "ElePtVsEta";
    cTemp = CreateCanvas(cFullName.c_str(), "colz", false, false, true, hElePtVsEta[iCanvas]);
    cElePtVsEta.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cElePtVsEta[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "EleE";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, true, hEleE[iCanvas]);
    cEleE.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cEleE[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "EleIDs";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, false, hEleIDs[iCanvas]);
    cEleIDs.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cEleIDs[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "EleDb";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, true, hEleDb[iCanvas]);
    cEleDb.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cEleDb[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "EleBsChargeCorr";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, false, hEleBsChargeCorr[iCanvas]);
    cEleBsChargeCorr.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cEleBsChargeCorr[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "EleGenChargeCorr";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, false, hEleGenChargeCorr[iCanvas]);
    cEleGenChargeCorr.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cEleGenChargeCorr[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "ElePassesHZZV1ID";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, true, hElePassesHZZV1ID[iCanvas]);
    cElePassesHZZV1ID.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cElePassesHZZV1ID[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "ElePassesHZZV2ID";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, true, hElePassesHZZV2ID[iCanvas]);
    cElePassesHZZV2ID.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cElePassesHZZV2ID[iCanvas];
    std::cout << "TCanvas " << cFullName << " produced.\n";

    cFullName = "c" + nameCanvas + "ElePVDistanceFromBestPV";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, false, hElePVZDistanceFromBestPV[iCanvas]);
    cElePVDistanceFromBestPV.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cElePVDistanceFromBestPV[iCanvas];
    std::cout << "TCanvas* " << cFullName << " produced.\n";
    
    cFullName = "c" + nameCanvas + "ElePassesMVANILID";
    cTemp = CreateCanvas(cFullName.c_str(), 0, 21, 1, false, true, hElePassesMVANILID[iCanvas]);
    cElePassesMVANILID.insert(std::make_pair(iCanvas, cTemp));
    autoSavedObject = cElePassesMVANILID[iCanvas];
    std::cout << "TCanvas* " << cFullName << " produced.\n";
    
    vcTemp = CreateCanvases(0, 21, 1, false, false, *vhEleEtaVsPt[iCanvas]);
    vcEleEtaVsPt.insert(std::make_pair(iCanvas, vcTemp));
    autoSavedObject = reinterpret_cast<std::vector<TObject*>*>(vcEleEtaVsPt[iCanvas]);
    std::cout << "std::vector<TCanvas*>* " << cFullName << " produced.\n";
    
    vcTemp = CreateCanvases(0, 21, 1, false, false, *vhEleHZZMVAOutputVsPt[iCanvas]);
    vcEleHZZMVAOutputVsPt.insert(std::make_pair(iCanvas, vcTemp));
    autoSavedObject = reinterpret_cast<std::vector<TObject*>*>(vcEleHZZMVAOutputVsPt[iCanvas]);
    std::cout << "std::vector<TCanvas*>* " << cFullName << " produced.\n";
    
    vcTemp = CreateCanvases(0, 21, 1, false, false, *vhEleHZZMVAOutputVsHZZMVACat[iCanvas]);
    vcEleHZZMVAOutputVsHZZMVACat.insert(std::make_pair(iCanvas, vcTemp));
    autoSavedObject = reinterpret_cast<std::vector<TObject*>*>(vcEleHZZMVAOutputVsHZZMVACat[iCanvas]);
    std::cout << "std::vector<TCanvas*>* " << cFullName << " produced.\n";

    vcTemp = CreateCanvases(0, 21, 1, false, false, *vhEleHZZMVACatVsPt[iCanvas]);
    vcEleHZZMVACatVsPt.insert(std::make_pair(iCanvas, vcTemp));
    autoSavedObject = reinterpret_cast<std::vector<TObject*>*>(vcEleHZZMVACatVsPt[iCanvas]);
    std::cout << "std::vector<TCanvas*>* " << cFullName << " produced.\n";
    
    vcTemp = CreateCanvases(0, 21, 1, false, false, *vhEleNIMVAOutputVsPt[iCanvas]);
    vcEleNIMVAOutputVsPt.insert(std::make_pair(iCanvas, vcTemp));
    autoSavedObject = reinterpret_cast<std::vector<TObject*>*>(vcEleNIMVAOutputVsPt[iCanvas]);
    std::cout << "std::vector<TCanvas*>* " << cFullName << " produced.\n";
    
    vcTemp = CreateCanvases(0, 21, 1, false, false, *vhEleNIMVAOutputVsNIMVACat[iCanvas]);
    vcEleNIMVAOutputVsNIMVACat.insert(std::make_pair(iCanvas, vcTemp));
    autoSavedObject = reinterpret_cast<std::vector<TObject*>*>(vcEleNIMVAOutputVsNIMVACat[iCanvas]);
    std::cout << "std::vector<TCanvas*>* " << cFullName << " produced.\n";

    vcTemp = CreateCanvases(0, 21, 1, false, false, *vhEleNIMVACatVsPt[iCanvas]);
    vcEleNIMVACatVsPt.insert(std::make_pair(iCanvas, vcTemp));
    autoSavedObject = reinterpret_cast<std::vector<TObject*>*>(vcEleNIMVACatVsPt[iCanvas]);
    std::cout << "std::vector<TCanvas*>* " << cFullName << " produced.\n";
    
    autoSavedObject = gDirectory->GetDirectory("..");
  }


  for(auto itEleGraphMap = eleGraphMap.begin(); itEleGraphMap != eleGraphMap.end(); itEleGraphMap++)
  {
    int iCanvas = itEleGraphMap->first;
    if(iCanvas == 0)
      continue;
    int nDigits = (int)pow(10,floor(log10(iCanvas)));
    std::string nameCanvas = itEleGraphMap->second;

    autoSavedObject = gDirectory->mkdir(("Canvas_eff_" + nameCanvas).c_str());
    
    
    std::cout << "Calculating efficiencies and producing canvases for iCanvas = " << iCanvas << " and nameCanvas = " << nameCanvas << std::endl;
    FeldmanCousinsDivide(hEleCharge[iCanvas], hEleCharge[iCanvas%nDigits], gEleEffVsCharge[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "EleEffVsCharge").c_str(), 21, 1, false, false, gEleEffVsCharge[iCanvas]);
  
    FeldmanCousinsDivide(hElePt[iCanvas], hElePt[iCanvas%nDigits], gEleEffVsPt[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "EleEffVsPt").c_str(), 21, 1, false, false, gEleEffVsPt[iCanvas]);
    
    FeldmanCousinsDivide(hElePhi[iCanvas], hElePhi[iCanvas%nDigits], gEleEffVsPhi[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "EleEffVsPhi").c_str(), 21, 1, false, false, gEleEffVsPhi[iCanvas]);
    
    FeldmanCousinsDivide(hEleEta[iCanvas], hEleEta[iCanvas%nDigits], gEleEffVsEta[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "EleEffVsEta").c_str(), 21, 1, false, false, gEleEffVsEta[iCanvas]);
    
    FeldmanCousinsDivide(hEleE[iCanvas], hEleE[iCanvas%nDigits], gEleEffVsE[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "EleEffVsE").c_str(), 21, 1, false, false, gEleEffVsE[iCanvas]);
    
    FeldmanCousinsDivide(hEleDb[iCanvas], hEleDb[iCanvas%nDigits], gEleEffVsDb[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "EleEffVsDb").c_str(), 21, 1, false, false, gEleEffVsDb[iCanvas]);

    FeldmanCousinsDivide(hEleBsChargeCorr[iCanvas], hEleBsChargeCorr[iCanvas%nDigits], gEleEffVsBsChargeCorr[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "EleEffVsBsChargeCorr").c_str(), 21, 1, false, false, gEleEffVsBsChargeCorr[iCanvas]);

    FeldmanCousinsDivide(hEleGenChargeCorr[iCanvas], hEleGenChargeCorr[iCanvas%nDigits], gEleEffVsGenChargeCorr[iCanvas]);
    autoSavedObject = CreateCanvas(("c" + nameCanvas + "EleEffVsGenChargeCorr").c_str(), 21, 1, false, false, gEleEffVsGenChargeCorr[iCanvas]);

    autoSavedObject = gDirectory->GetDirectory("..");    
  }

    
  FeldmanCousinsDivide(hElePt[1+10+100], hElePt[1+10], gEleEffVsPtTest);
  autoSavedObject = CreateCanvas("cEleEffVsPtTest", 21, 1, false, false, gEleEffVsPtTest);
  
  return;
}
