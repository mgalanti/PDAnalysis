#include <iostream>
#include <sstream>
#include <string>
#include <bitset>
#include <math.h>

#include "TSystem.h"

#include "EleMVASecondNtupleProducer.h"

#include "PDAnalysis/EleMVASecondNtupleProducer/interface/EleMVASecondNtupleWriter.h"

#include "TDirectory.h"
#include "TBranch.h"
#include "TTree.h"
#include "TCanvas.h"
#include "Math/LorentzVector.h"

using namespace std;



EleMVASecondNtupleProducer::~EleMVASecondNtupleProducer() {
}



void EleMVASecondNtupleProducer::beginJob() {

  MGBaseAnalyzer::beginJob();
  MGSelector::SetSelectionString(evtSelection);


  // user parameters are retrieved as strings by using their names;
  // numeric parameters ( int, float or whatever ) can be directly set
  // by passing the corresponding variable,
  // e.g. getUserParameter( "name", x )

  getUserParameter("tightSelection", tightSelection);
  getUserParameter("nConeIterations", nConeIterations);
  getUserParameter("coneTolerance", coneTolerance);
  
  // user parameters are set as names associated to a string, 
  // default values can be set in the analyzer class contructor
  
  tWriter = new EleMVASecondNtupleWriter; // second ntuple
  
  std::string secondNtupleBaseName = getUserParameter("secondNtupleBaseName");
  
  // Set the second ntuple file name
  int tries = 0;
  secondNtupleFileName = secondNtupleBaseName + "__" + sampleName + "__" + evtSelection + "__" + std::to_string(tries++) + ".root";
  std::cout << "EleMVASecondNtupleProducer::beginJob(): setting variables:\n";
  std::cout << "                                        secondNtupleFileName = " << secondNtupleFileName << std::endl;
  bool secondNtupleOrigNameNotOk = false;
  while(!gSystem->AccessPathName(gSystem->ExpandPathName(secondNtupleFileName.c_str())))
  {
    secondNtupleOrigNameNotOk = true;
    if(tries == 1)
    {
      std::cout << "W A R N I N G! Second ntuple file \"" << secondNtupleFileName << "\" already exists!\n";
      std::cout << "               Trying a different name...\n";
    }
    secondNtupleFileName = secondNtupleBaseName + "__" + sampleName + "__" + evtSelection + "__" + std::to_string(tries++) + ".root";
    if(tries > 9999)
      break;
  }
  if(!gSystem->AccessPathName(gSystem->ExpandPathName(secondNtupleFileName.c_str())))
  {
    std::cout << "E R R O R! Output file \"" << secondNtupleFileName << "\" already exists!\n";
    std::cout << "           Exiting...\n";
    exit(1);
  }
  if(secondNtupleOrigNameNotOk)
    std::cout << "                                        New secondNtupleFileName = " << secondNtupleFileName << std::endl;
  
  // Actually open the second ntuple file for writing
  tWriter->open(secondNtupleFileName, "RECREATE");
  
  // Set the generator-level B meson mass according to the particle looked for in the selection
  trueBMass = 0;
  if(evtSelection.substr(0,2).compare("Bs") == 0)
    trueBMass = constants::BsMass;
  else if(evtSelection.substr(0,2).compare("Bu") == 0)
    trueBMass = constants::BuMass;
  else
  {
    std::cout << "E R R O R ! Cannot determine the selected meson from the selection string!\n";
    std::cout << "            The selection string is: \"" << evtSelection << "\"\n";
    std::cout << "            Please fix the configuration file.\n";
    std::cout << "            Exiting...\n";
    exit(1);
  }
  std::cout << "EleMVASecondNtupleProducer::beginJob(): I N F O. trueBMass = " << trueBMass << std::endl;

  return;
}


void EleMVASecondNtupleProducer::book() 
{

  // putting "autoSavedObject" in front of the histo creation 
  // it's automatically marked for saving on file; the option 
  // is uneffective when not using the full utility
  
  float minMass = 5.0;
  float maxMass = 5.5;
  int nMassBins = 250;

  float minDr = 0;
  float maxDr = 5.0;
  int nDrBins = 100;

  float minPtRel = 0;
  float maxPtRel = 50.0;
  int nPtRelBins = 100;
  
  autoSavedObject = hEleBTrkDistance = Create2DHistogram<TH2D>("hEleBTrkDistance", "Distance in #Delta(R) and #Delta(p_{T})/p_{T} between electron and B tracks", 100, 0., 5., 100, 0., 2., "#Delta(R)", "2#times|p_{T,ele}-p_{T,Trk}|/(p_{T,ele}+p_{T,Trk})");
  
  autoSavedObject = hWeightedBMass = Create1DHistogram<TH1D>("hWeightedBMass", "B mass (entries weighted by number of gen B in the event)", nMassBins, minMass, maxMass, "M(B) [GeV]", "Event weights");
  
  autoSavedObject = hEleConeDistance = Create1DHistogram<TH1D>("hEleConeDistance", "Distance between electron and cone axis", nDrBins, minDr, maxDr, "#Delta(R)", "N. Electrons");
  
  autoSavedObject = hEleConePtRel = Create1DHistogram<TH1D>("hEleConePtRel", "PtRel of electron with respect to cone", nPtRelBins, minPtRel, maxPtRel, "PtRel", "N. Electrons");
  
  autoSavedObject = hEleConeCleanDistance = Create1DHistogram<TH1D>("hEleConeCleanDistance", "Distance between electron and coneclean axis", nDrBins, minDr, maxDr, "#Delta(R)", "N. Electrons");
  
  autoSavedObject = hEleConeCleanPtRel = Create1DHistogram<TH1D>("hEleConeCleanPtRel", "PtRel of electron with respect to coneclean", nPtRelBins, minPtRel, maxPtRel, "PtRel", "N. Electrons");
  
  return;
}


void EleMVASecondNtupleProducer::reset() {
// automatic reset
  autoReset();
  return;
}


bool EleMVASecondNtupleProducer::analyze(int entry, int event_file, int event_tot) 
{
  if (verbose) 
  {
    cout << " +++++++++++++++++++++++++++ " << endl;
    cout << "entry: "
         << entry << " " << event_file << " " << event_tot << endl;
    cout << "run: " <<   runNumber << " , "
         << "evt: " << eventNumber << endl;
  }
  else 
  {
//    if (!(event_file % 10000) || !(event_tot % 10000))
    if (!( event_tot % 10000 ) && event_tot)
      cout << event_file << " " << event_tot << endl;
  }
  if(entry == 0)
    checkBranches();
  setCurrentEvt(event_tot);
  
  // Main per-event analysis code goes here
  
  tWriter->Reset();

  // Check if the event passes the default selection defined in the cfg
  bool evtSelected = SelectEvent();
  
//   if(!evtSelected)
//     return false;
 
  // Check if the event passes the tight selection
  std::vector<std::pair<int, int> > selectedObjectsTight;
  bool tightEvent = SelectEvent(tightSelection, selectedObjectsTight);
  
  bool tightB = false;
  for (auto itSelObjectsTight = selectedObjectsTight.begin(); itSelObjectsTight != selectedObjectsTight.end(); itSelObjectsTight++)
  {
    if(itSelObjectsTight->first == PDEnumString::recSvt)
      tightB = true;
  }
  
  // If the event passes the tight selection, let's use those candidates instead of the loose ones
  if(tightEvent)
    selectedObjects = selectedObjectsTight;

  if(verbose && tightB && !tightEvent)
  {
    std::cout << "I N F O : signal B is tight but event is not tight!\n";
  }

  if(!tightB && tightEvent)
  {
    std::cout << "E R R O R ! signal B is not tight but event is tight!\n";
    exit(1);
  }
  
  int iSelObject = 0;
  int iBestPV = -1;
  int iBestB = -1;
  int iBestEle = -1;
  for (auto itSelObjects = selectedObjects.begin(); itSelObjects != selectedObjects.end(); itSelObjects++)
  {
//     std::cout << "Event selection: selected object #" << iSelObject << ": type = " << itSelObjects->first << ", index = " << itSelObjects->second << std::endl;
    if(itSelObjects->first == PDEnumString::recPV)
      iBestPV = itSelObjects->second;
    if(itSelObjects->first == PDEnumString::recSvt)
      iBestB = itSelObjects->second;
    if(itSelObjects->first == PDEnumString::recElectron)
      iBestEle = itSelObjects->second;
    iSelObject++;
  }
  
  // Reject the event if the signal side has not been selected
  if(iBestPV == -1 || iBestB == -1)
    return false;
  
//   // Depending on the selection string provided in the configuration, it is not granted 
//   // that we have all the needed objects at this point, even if the event passed the selection.
//   // Thus, let's check explicitly.  
//   if(iBestPV == -1 || iBestB == -1 || iBestEle == -1)
//   {
//     std::cout << "E R R O R ! Event passed the selection \"" << evtSelection << "\" but not all the signal-side objects are available!\n";
//     std::cout << "            iBestPV = " << iBestPV << ", iBestB = " << iBestB << ", iBestEle = " << iBestEle << std::endl;
//     std::cout << "            Please fix the configuration file or the MGSelector::SelectEvent() code to avoid this error.";
//     std::cout << "            Exiting...\n";
//     exit(1);
//   }

  // HLT bit values
  bool jPsiMuHltBit = hlt(PDEnumString::HLT_Dimuon0_Jpsi3p5_Muon2_v) || hlt(PDEnumString::HLT_Dimuon0_Jpsi_Muon_v);
  bool jPsiTrkTrkHltBit = hlt(PDEnumString::HLT_DoubleMu4_JpsiTrkTrk_Displaced_v);
  bool jPsiTrkHltBit = hlt(PDEnumString::HLT_DoubleMu4_JpsiTrk_Displaced_v);
  
  std::vector<int> tracksFromB = tracksFromSV(iBestB);
  TLorentzVector pB = GetTLorentzVectorFromJPsiX(iBestB);
  int iBestJPsi = (subVtxFromSV(iBestB)).at(0);
  std::vector<int> tracksFromJPsi = tracksFromSV(iBestJPsi);

  // Truth information for B
  int iGenB = -1;
  int nGenB = -1;
  int idGenB = 0;
  int genBMixStatus = -1;
  int evtWeight = 1;

  // Other variables needed for the ntuple production
  std::vector<int> allBHadrons;
  std::vector<int> longLivedBHadrons;

  // Use generator information if available
  if(has_gen)
  {
    // Truth information for B
    allBHadrons = GetAllBHadrons(true);  // true removes particles with mix status = 2
    longLivedBHadrons = GetAllLongLivedBHadrons(true); // true removes particles with mix status = 2
    iGenB = GetClosestGenInList(pB.Pt(), pB.Eta(), pB.Phi(), longLivedBHadrons, 0.4, 0.4);
    nGenB = longLivedBHadrons.size();
    // Do not consider events where the signal-side is not matched to a gen B
    if(iGenB < 0)
      return false;

    idGenB = genId->at(iGenB);

    // Do not consider events where the signal B is matched to a hadron different from the one looked for
    if(selectionSubStrings[0].compare("BsToJPsiPhi") == 0 && abs(idGenB) != 531)
    {
      std::cout << "I N F O: Event rejected - looking for BsToJPsiPhi but signal B is matched to a hadron with id = " << idGenB << ".\n";
      return false;
    }
    if(selectionSubStrings[0].compare("BuToJPsiK") == 0 && abs(idGenB) != 521)
    {
      std::cout << "I N F O: Event rejected - looking for BuToJPsiK but signal B is matched to a hadron with id = " << idGenB << ".\n";
      return false;
    }
    // Check whether the matched B has mixed or will mix
    genBMixStatus = GetMixStatus(iGenB);
    
    if(genBMixStatus == 2)
    {
      std::cout << "I N F O: Event rejected - mix status of generated B is 2.\n";
      return false;
    }
    
    // If the matched B hadron comes from mixing, then go back to the pre-mixing one
    if(genBMixStatus == 1)
    {
      iGenB = RecursiveLookForMotherIds(iGenB, {-idGenB});
      idGenB = genId->at(iGenB);
    }

    // Change event weight according to how many "interesting" B hadrons are there
    for(auto iGen: longLivedBHadrons)
    {
      if(iGen == iGenB)
        continue;
      if(abs(genId->at(iGen)) == abs(idGenB))
        evtWeight = 2;
    }
    
  }
  else  // If no gen information is available, only the BuToJPsiK channel can be useful
  {
    // Set id randomly for Bs - this should not be used
    if(selectionSubStrings[0].compare("BsToJPsiPhi") == 0)
    {
      idGenB = ((double)rand() / (RAND_MAX)) < 0.5 ? +531 : -531;
      std::cout << "EleMVASecondNtupleProducer::analyze():  W A R N I N G ! Trying to save tagging information for BsToJPsiPhi but\n";
      std::cout << "                                                        no generator information is available in the input files!\n";
      std::cout << "                                                        idGenB is filled randomly to " << idGenB << "!\n";
    }
    // Bu is a self-tagging channel - OK!
    if(selectionSubStrings[0].compare("BuToJPsiK") == 0)
    {
      for(auto iTrk: tracksFromB)
      {
        if(iTrk == tracksFromJPsi[0] || iTrk == tracksFromJPsi[1])
          continue;
        idGenB = trkCharge->at(iTrk) > 0 ? +521 : -521;
      }
      
    }
  }
  
  // Write general event variables to secondary ntuple
  (tWriter->evtNumber) = event_tot;
  (tWriter->evtWeight) = evtWeight;
  (tWriter->tightEvent) = tightEvent;

  (tWriter->nGenB) = nGenB;
  
  (tWriter->JPsiMuHltBit) = jPsiMuHltBit;
  (tWriter->JPsiTrkTrkHltBit) = jPsiTrkTrkHltBit;
  (tWriter->JPsiTrkHltBit) = jPsiTrkHltBit;  
  
  (tWriter->iPV) = iBestPV;
  
  (tWriter->tightB) = tightB;
    
  // Write signal-side variables to secondary ntuple
  (tWriter->BPt) = pB.Pt();
  (tWriter->BEta) = pB.Eta();
  (tWriter->BPhi) = pB.Phi();
  (tWriter->BMass) = svtMass->at(iBestB);
  
  hWeightedBMass->Fill(svtMass->at(iBestB), evtWeight);
  
  (tWriter->BLxy) = GetCt2D(pB, iBestB) / (trueBMass / pB.Pt());
  (tWriter->BCt2DBS) = GetCt2D(pB, iBestB, trueBMass);
  
  (tWriter->BCt2DPV) = GetCt2DPV(pB, iBestB, iBestPV, trueBMass);
  (tWriter->BCt2DPVErr) = GetCt2DPVErr(pB, iBestB, iBestPV, trueBMass);
  (tWriter->BCt2DPVSigmaUnit) = GetCt2DPV(pB, iBestB, iBestPV, trueBMass) / GetCt2DPVErr(pB, iBestB, iBestPV, trueBMass);

  (tWriter->BCt3DPV) = GetCt3DPV(pB, iBestB, iBestPV, trueBMass);
  (tWriter->BCt3DPVErr) = GetCt3DPVErr(pB, iBestB, iBestPV, trueBMass);
  (tWriter->BCt3DPVSigmaUnit) = GetCt3DPV(pB, iBestB, iBestPV, trueBMass) / GetCt3DPVErr(pB, iBestB, iBestPV, trueBMass);

  (tWriter->BiSV) = iBestB;
  
  (tWriter->BidGen) = idGenB;
  
  // If no electron has been selected, do not write detailed tag-side information to the ntuple and pass to the next event
  if(iBestEle == -1)
  {
    (tWriter->eleSelected) = 0;
    (tWriter->tagTruth) = -1;
    (tWriter->chargeCorr) = 0;
    tWriter->fill();
    
    return true;
  }
  
  // If we are here, then an electron has been selected. Let's use it.
  (tWriter->eleSelected) = 1;

  // For monitoring purposes
  for(auto iTrk: tracksFromB)
  {
    float dR = deltaR(trkEta->at(iTrk), trkPhi->at(iTrk), eleGsfEta->at(iBestEle), eleGsfPhi->at(iBestEle));
    float dpTOverpT = 2*fabs(trkPt->at(iTrk) - eleGsfPt->at(iBestEle))/(trkPt->at(iTrk) + eleGsfPt->at(iBestEle));
    hEleBTrkDistance->Fill(dR, dpTOverpT);
  }
    
  // Truth information for electron
  int iGenEle = -1;
  int idGenEle = 0;
  int genEleBMot = -1;
  
  // Variables needed to compute the tagging truth (stored into tagTruth)
  int chargeB = 0;
  int chargeCorr = 0;
  int chargeEle = eleCharge->at(iBestEle);
  
  // Use generator information if available
  if(has_gen)
  {
    // Truth information for electron
    iGenEle = GetClosestGen(elePt->at(iBestEle), eleEta->at(iBestEle), elePhi->at(iBestEle));
    if (iGenEle >= 0)
    {
      idGenEle = genId->at(iGenEle);
      genEleBMot = RecursiveLookForMother(iGenEle, allBHadrons);
      if(genEleBMot < 0)
        genEleBMot = -1;
    }
    if(selectionSubStrings[0].compare("BsToJPsiPhi") == 0)
    {
      if(abs(idGenEle) == 11)
        chargeCorr = GetGenLepBsChargeCorrelation(iGenEle, iGenB);
      if(!chargeCorr)
      {
        chargeCorr = GetBsChargeCorrelation(chargeEle, iGenB);
        chargeCorr*=2;
      }
    }
    else if(selectionSubStrings[0].compare("BuToJPsiK") == 0)
    {
      if(abs(idGenEle) == 11)
        chargeCorr = GetGenLepBuChargeCorrelation(iGenEle, iGenB);
      if(!chargeCorr)
      {
        chargeCorr = GetBuChargeCorrelation(chargeEle, iGenB);
        chargeCorr*=2;
      }
    }    
  }  
  else  // If no gen information is available, only the BuToJPsiK channel can be useful
  {
    // Bu is a self-tagging channel - OK!
    if(selectionSubStrings[0].compare("BuToJPsiK") == 0)
    {
      chargeB = 0;
      for(auto iTrack : tracksFromB)
      {
        chargeB += trkCharge->at(iTrack);
      }
      
      chargeCorr = chargeEle * chargeB;
    }    
  }
  
  int tagTruth;
  if(chargeCorr > 0)
    tagTruth = 1;
  else
    tagTruth = 0;
  
  // EleID Variables
  float eleIDNIV2Val = -1;
  float eleIDIV2Val = -1;
  float eleIDHZZV1Val = -1;
  float eleIDNIV2RawVal = -10;
  float eleIDIV2RawVal = -10;
  float eleIDHZZV1RawVal = -10;
  int eleIDNIV2Cat = 0;
  int eleIDIV2Cat = 0;
  int eleIDHZZV1Cat = 0;
  
  for (int iUserInfo = 0; iUserInfo < nUserInfo; iUserInfo++)
  {
//     std::cout << "Looping on userInfo object #" << iUserInfo << std::endl;
    if(useObjType->at(iUserInfo) == PDEnumString::recElectron && useObjIndex->at(iUserInfo) == iBestEle)
    {
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17NoIsoV2Values)
      {
        eleIDNIV2Val = useInfoValue->at(iUserInfo);
      }
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17IsoV2Values)
      {
        eleIDIV2Val = useInfoValue->at(iUserInfo);
      }
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Spring16HZZV1Values)
      {
        eleIDHZZV1Val = useInfoValue->at(iUserInfo);
      }
      
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17NoIsoV2RawValues)
      {
        eleIDNIV2RawVal = useInfoValue->at(iUserInfo);
      }
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17IsoV2RawValues)
      {
        eleIDIV2RawVal = useInfoValue->at(iUserInfo);
      }
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Spring16HZZV1RawValues)
      {
        eleIDHZZV1RawVal = useInfoValue->at(iUserInfo);
      }
      
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17NoIsoV2Categories)
      {
        eleIDNIV2Cat = (int)(useInfoValue->at(iUserInfo));
      }
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17IsoV2Categories)
      {
        eleIDIV2Cat = (int)(useInfoValue->at(iUserInfo));
      }
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Spring16HZZV1Categories)
      {
        eleIDHZZV1Cat = (int)(useInfoValue->at(iUserInfo));
      }
    }
  }

  // CONE variables
  float eleConePtRel = -1;
  float eleConeDR = -1;
  float eleConeEnergyRatio = -1;
  int   eleConeSize = 0;
  float eleConeQ = -1;
  float eleConePt = -1;
  float eleConeNF = 0;
  float eleConeCF = 0;
  int   eleConeNCH = 0;
  std::vector<float> eleConeDx;
  float eleConeAvgDx = 0;
  float eleConeStdDevDx = 0;
  std::vector<float> eleConeDy;
  float eleConeAvgDy = 0;
  float eleConeStdDevDy = 0;
  std::vector<float> eleConeDxy;
  float eleConeAvgDxy = 0;
  float eleConeStdDevDxy = 0;
  std::vector<float> eleConeDz;
  float eleConeAvgDz = 0;
  float eleConeStdDevDz = 0;
  float kappa = 1;
  float drCone = 0.4;
  float axisEta = eleEta->at(iBestEle);
  float axisPhi = elePhi->at(iBestEle);
  TLorentzVector pEle;
  pEle.SetPtEtaPhiM(elePt->at(iBestEle), eleEta->at(iBestEle), elePhi->at(iBestEle), constants::electronMass);
  
  if(verbose)
  {
  std::cout << "Elec: eta = " << eleEta->at(iBestEle) << ", phi = " << elePhi->at(iBestEle) << ", pt = " << elePt->at(iBestEle) << std::endl;
  }
  
  float qCone = 0, ptCone = 0;
  for(int i = 0; i < nConeIterations; i++)
  {
    if(verbose)
    {
      std::cout << "   Iteration " << i << std::endl;
    }
    eleConePtRel = -1;
    eleConeDR = -1;
    eleConeEnergyRatio = -1;
    eleConeSize = 0;
    eleConeQ = -1;
    eleConePt = -1;
    eleConeNF = 0;
    eleConeCF = 0;
    eleConeNCH = 0;
    eleConeDx.clear();
    eleConeAvgDx = 0;
    eleConeStdDevDx = 0;
    eleConeDy.clear();
    eleConeAvgDy = 0;
    eleConeStdDevDy = 0;
    eleConeDxy.clear();
    eleConeAvgDxy = 0;
    eleConeStdDevDxy = 0;
    eleConeDz.clear();
    eleConeAvgDz = 0;
    eleConeStdDevDz = 0;
    
    TLorentzVector pCone(0.,0.,0.,0.);
    for(int iPF=0; iPF < nPF; ++iPF)
    {
      float ptPF = pfcPt->at(iPF);
      float etaPF = pfcEta->at(iPF);
      if(deltaR(etaPF, pfcPhi->at(iPF), axisEta, axisPhi) > drCone)
      {
        continue;
      }
      if(ptPF < 0.2)
      {
        continue;
      }
      if(fabs(etaPF) > 3.0)
      {
        continue;  
      }
      if(std::find(tracksFromB.begin(), tracksFromB.end(), pfcTrk->at(iPF)) != tracksFromB.end())
      {
        continue;
      }
      if(pfcTrk->at(iPF) >= 0)
      {
        int iTrk = pfcTrk->at(iPF);
//         std::cout << "iBestPV = " << iBestPV << ", nPVertices = " << nPVertices << std::endl;
        if(fabs(dZTrk(iTrk, iBestPV)) >= 1.0)
        {
          continue;
        }
        eleConeDx.push_back(dXTrk(iTrk, iBestPV));
        eleConeDy.push_back(dYTrk(iTrk, iBestPV));
        eleConeDxy.push_back(dXYTrk(iTrk, iBestPV));
        eleConeDz.push_back(dZTrk(iTrk, iBestPV));
      }
//       std::cout << "iPF = " << iPF << ", nPF = " << nPF << std::endl;
      TLorentzVector a;
      a.SetPtEtaPhiE(pfcPt->at(iPF), pfcEta->at(iPF), pfcPhi->at(iPF), pfcE->at(iPF));
      pCone += a;
      ++eleConeSize;
      
      qCone += pfcCharge->at(iPF) * pow(ptPF, kappa);
      ptCone += pow(ptPF, kappa);

      if(pfcCharge->at(iPF)==0)
      {
        eleConeNF += pfcE->at(iPF);
      }
      if(abs(pfcCharge->at(iPF))==1)
      {
        eleConeNCH++;
        eleConeCF += pfcE->at(iPF);
      }
    }
    
    if(ptCone != 0)
    {
      qCone /= ptCone;
    }
    else
    {
      qCone = 1;
    }
    qCone *= eleCharge->at(iBestEle);
    if(pCone.E() != 0)
    {
      eleConeCF /= pCone.E();
      eleConeNF /= pCone.E();
    }
    
    eleConeQ = qCone;
    
    eleConePt = pCone.Pt();
    eleConeDR = deltaR(pCone.Eta(), pCone.Phi(), eleEta->at(iBestEle), elePhi->at(iBestEle));
    if(pCone.E() !=0)
    {
      eleConeEnergyRatio = eleE->at(iBestEle) / pCone.E();
    }
    else
    {
      eleConeEnergyRatio = 1;
    }
    pCone -= pEle;
    eleConePtRel = elePt->at(iBestEle) * (pEle.Vect().Unit() * pCone.Vect().Unit());
    
    for(unsigned int i = 0; i < eleConeDx.size(); i++)
    {
      eleConeAvgDx+=eleConeDx[i];
      eleConeAvgDy+=eleConeDy[i];
      eleConeAvgDxy+=eleConeDxy[i];
      eleConeAvgDz+=eleConeDz[i];
    }
    if(eleConeDx.size())
    {
      eleConeAvgDx/=eleConeDx.size();
      eleConeAvgDy/=eleConeDy.size();
      eleConeAvgDxy/=eleConeDxy.size();
      eleConeAvgDz/=eleConeDz.size();
    }
    for(unsigned int i = 0; i < eleConeDx.size(); i++)
    {
      eleConeStdDevDx+=(eleConeDx[i]-eleConeAvgDx)*(eleConeDx[i]-eleConeAvgDx);
      eleConeStdDevDy+=(eleConeDy[i]-eleConeAvgDy)*(eleConeDy[i]-eleConeAvgDy);
      eleConeStdDevDxy+=(eleConeDxy[i]-eleConeAvgDxy)*(eleConeDxy[i]-eleConeAvgDxy);
      eleConeStdDevDz+=(eleConeDz[i]-eleConeAvgDz)*(eleConeDz[i]-eleConeAvgDz);
    }
    if(eleConeDx.size())
    {
      eleConeStdDevDx/=eleConeDx.size();
      eleConeStdDevDx = sqrt(eleConeStdDevDx);
      eleConeStdDevDy/=eleConeDy.size();
      eleConeStdDevDy = sqrt(eleConeStdDevDy);
      eleConeStdDevDxy/=eleConeDxy.size();
      eleConeStdDevDxy = sqrt(eleConeStdDevDxy);
      eleConeStdDevDz/=eleConeDz.size();
      eleConeStdDevDz = sqrt(eleConeStdDevDz);
    }
    
    if(verbose)
    {
      std::cout << "      Cone: eta = " << pCone.Eta() << ", phi = " << pCone.Phi() << ", pt = " << pCone.Pt() <<  std::endl; 
      std::cout << "            eleConePt = " << eleConePt << std::endl;
      std::cout << "            ptRel = " << eleConePtRel << std::endl;
      std::cout << "            eleConeDR = " << eleConeDR << std::endl;
      std::cout << "            eleConeEnergyRatio = " << eleConeEnergyRatio << std::endl;
      std::cout << "            eleConeQ = " << eleConeQ << std::endl;
    }
    
    float distance = deltaR(axisEta, axisPhi, pCone.Eta(), pCone.Phi());
    
    if(verbose)
    {
      std::cout << "            Distance = " << distance << std::endl;
    }
    
    if(distance > coneTolerance)
    {
      axisEta = pCone.Eta();
      axisPhi = pCone.Phi();
    }
    else
    {
      break;
    }
  }
  
  hEleConeDistance->Fill(eleConeDR);
  hEleConePtRel->Fill(eleConePtRel);  

  // CONECLEAN variables
  float eleConeCleanPtRel = -1;
  float eleConeCleanDR = -1;
  float eleConeCleanEnergyRatio = -1;
  int   eleConeCleanSize = 0;
  float eleConeCleanQ = -1;
  float eleConeCleanPt = -1;
  float eleConeCleanNF = 0;
  float eleConeCleanCF = 0;
  int   eleConeCleanNCH = 0;
    std::vector<float> eleConeCleanDx;
  float eleConeCleanAvgDx = 0;
  float eleConeCleanStdDevDx = 0;
  std::vector<float> eleConeCleanDy;
  float eleConeCleanAvgDy = 0;
  float eleConeCleanStdDevDy = 0;
  std::vector<float> eleConeCleanDxy;
  float eleConeCleanAvgDxy = 0;
  float eleConeCleanStdDevDxy = 0;
  std::vector<float> eleConeCleanDz;
  float eleConeCleanAvgDz = 0;
  float eleConeCleanStdDevDz = 0;
  kappa = 1;
  drCone = 0.4;
  axisEta = eleEta->at(iBestEle);
  axisPhi = elePhi->at(iBestEle);
//   TLorentzVector pEle;
//   pEle.SetPtEtaPhiM(elePt->at(iBestEle), eleEta->at(iBestEle), elePhi->at(iBestEle), constants::electronMass);
//   
//   if(verbose)
//   {
//   std::cout << "Elec: eta = " << eleEta->at(iBestEle) << ", phi = " << elePhi->at(iBestEle) << ", pt = " << elePt->at(iBestEle) << std::endl;
//   }
//   
  float qConeClean = 0, ptConeClean = 0;
  TLorentzVector pConeClean(0.,0.,0.,0.);
  for(int i = 0; i < nConeIterations; i++)
  {
    if(verbose)
    {
      std::cout << "   Iteration " << i << std::endl;
    }
    eleConeCleanPtRel = -1;
    eleConeCleanDR = -1;
    eleConeCleanEnergyRatio = -1;
    eleConeCleanSize = 0;
    eleConeCleanQ = -1;
    eleConeCleanPt = -1;
    eleConeCleanNF = 0;
    eleConeCleanCF = 0;
    eleConeCleanNCH = 0;
    eleConeCleanDx.clear();
    eleConeCleanAvgDx = 0;
    eleConeCleanStdDevDx = 0;
    eleConeCleanDy.clear();
    eleConeCleanAvgDy = 0;
    eleConeCleanStdDevDy = 0;
    eleConeCleanDxy.clear();
    eleConeCleanAvgDxy = 0;
    eleConeCleanStdDevDxy = 0;
    eleConeCleanDz.clear();
    eleConeCleanAvgDz = 0;
    eleConeCleanStdDevDz = 0;

    pConeClean.SetPtEtaPhiE(0.,0.,0.,0.);
//     TLorentzVector pConeClean(0.,0.,0.,0.);
    for(int iPF = 0; iPF < nPF; ++iPF)
    {
      float ptPF = pfcPt->at(iPF);
      float etaPF = pfcEta->at(iPF);
      if(deltaR(etaPF, pfcPhi->at(iPF), axisEta, axisPhi) > drCone)
      {
        continue;
      }
      if(ptPF < 0.5)
      {
        continue;
      }
      if(fabs(etaPF) > 3.0)
      {
        continue;  
      }
      if(pfcCharge->at(iPF) == 0)
      {
        continue;
      }
      if(std::find(tracksFromB.begin(), tracksFromB.end(), pfcTrk->at(iPF)) != tracksFromB.end())
      {
        continue;
      }
      int iTrk = pfcTrk->at(iPF);
      if(iTrk < 0)
      {
        continue;
      }
      if(fabs(dZTrk(iTrk, iBestPV)) >= 1.0)
      {
        continue;
      }
      eleConeCleanDx.push_back(dXTrk(iTrk, iBestPV));
      eleConeCleanDy.push_back(dYTrk(iTrk, iBestPV));
      eleConeCleanDxy.push_back(dXYTrk(iTrk, iBestPV));
      eleConeCleanDz.push_back(dZTrk(iTrk, iBestPV));
//       std::cout << "iPF = " << iPF << ", nPF = " << nPF << std::endl;
      TLorentzVector a;
      a.SetPtEtaPhiE(pfcPt->at(iPF), pfcEta->at(iPF), pfcPhi->at(iPF), pfcE->at(iPF));
      pConeClean += a;
      ++eleConeCleanSize;
      
      qConeClean += pfcCharge->at(iPF) * pow(ptPF, kappa);
      ptConeClean += pow(ptPF, kappa);

      if(pfcCharge->at(iPF)==0)
      {
        eleConeCleanNF += pfcE->at(iPF);
      }
      if(abs(pfcCharge->at(iPF))==1)
      {
        eleConeCleanNCH++;
        eleConeCleanCF += pfcE->at(iPF);
      }
    }
    
    if(ptConeClean != 0)
    {
      qConeClean /= ptConeClean;
    }
    else
    {
      qConeClean = 1;
    }
    qConeClean *= eleCharge->at(iBestEle);
    if(pConeClean.E() != 0)
    {
      eleConeCleanCF /= pConeClean.E();
      eleConeCleanNF /= pConeClean.E();
    }
    
    eleConeCleanQ = qConeClean;
    
    eleConeCleanPt = pConeClean.Pt();
    eleConeCleanDR = deltaR(pConeClean.Eta(), pConeClean.Phi(), eleEta->at(iBestEle), elePhi->at(iBestEle));
    if(pConeClean.E() !=0)
    {
      eleConeCleanEnergyRatio = eleE->at(iBestEle) / pConeClean.E();
    }
    else
    {
      eleConeCleanEnergyRatio = 1;
    }
    pConeClean -= pEle;
    eleConeCleanPtRel = elePt->at(iBestEle) * (pEle.Vect().Unit() * pConeClean.Vect().Unit());
    pConeClean += pEle; // for IP sign
    
    for(unsigned int i = 0; i < eleConeCleanDx.size(); i++)
    {
      eleConeCleanAvgDx+=eleConeCleanDx[i];
      eleConeCleanAvgDy+=eleConeCleanDy[i];
      eleConeCleanAvgDxy+=eleConeCleanDxy[i];
      eleConeCleanAvgDz+=eleConeCleanDz[i];
    }
    if(eleConeCleanDx.size())
    {
      eleConeCleanAvgDx/=eleConeCleanDx.size();
      eleConeCleanAvgDy/=eleConeCleanDy.size();
      eleConeCleanAvgDxy/=eleConeCleanDxy.size();
      eleConeCleanAvgDz/=eleConeCleanDz.size();
    }
    for(unsigned int i = 0; i < eleConeCleanDx.size(); i++)
    {
      eleConeCleanStdDevDx+=(eleConeCleanDx[i]-eleConeCleanAvgDx)*(eleConeCleanDx[i]-eleConeCleanAvgDx);
      eleConeCleanStdDevDy+=(eleConeCleanDy[i]-eleConeCleanAvgDy)*(eleConeCleanDy[i]-eleConeCleanAvgDy);
      eleConeCleanStdDevDxy+=(eleConeCleanDxy[i]-eleConeCleanAvgDxy)*(eleConeCleanDxy[i]-eleConeCleanAvgDxy);
      eleConeCleanStdDevDz+=(eleConeCleanDz[i]-eleConeCleanAvgDz)*(eleConeCleanDz[i]-eleConeCleanAvgDz);
    }
    if(eleConeCleanDx.size())
    {
      eleConeCleanStdDevDx/=eleConeCleanDx.size();
      eleConeCleanStdDevDx = sqrt(eleConeCleanStdDevDx);
      eleConeCleanStdDevDy/=eleConeCleanDy.size();
      eleConeCleanStdDevDy = sqrt(eleConeCleanStdDevDy);
      eleConeCleanStdDevDxy/=eleConeCleanDxy.size();
      eleConeCleanStdDevDxy = sqrt(eleConeCleanStdDevDxy);
      eleConeCleanStdDevDz/=eleConeCleanDz.size();
      eleConeCleanStdDevDz = sqrt(eleConeCleanStdDevDz);
    }
    
    if(verbose)
    {
      std::cout << "      ConeClean: eta = " << pConeClean.Eta() << ", phi = " << pConeClean.Phi() << ", pt = " << pConeClean.Pt() <<  std::endl; 
      std::cout << "                 eleConeCleanPt = " << eleConeCleanPt << std::endl;
      std::cout << "                 ptRel = " << eleConeCleanPtRel << std::endl;
      std::cout << "                 eleConeCleanDR = " << eleConeCleanDR << std::endl;
      std::cout << "                 eleConeCleanEnergyRatio = " << eleConeCleanEnergyRatio << std::endl;
      std::cout << "                 eleConeCleanQ = " << eleConeCleanQ << std::endl;
    }
    
    float distance = deltaR(axisEta, axisPhi, pConeClean.Eta(), pConeClean.Phi());
    
    if(verbose)
    {
      std::cout << "            Distance = " << distance << std::endl;
    }
    
    if(distance > coneTolerance)
    {
      axisEta = pConeClean.Eta();
      axisPhi = pConeClean.Phi();
    }
    else
    {
      break;
    }
  }

  hEleConeCleanDistance->Fill(eleConeCleanDR);
  hEleConeCleanPtRel->Fill(eleConeCleanPtRel);

  // Opposite-side variables
  (tWriter->elePt) = elePt->at(iBestEle);
  (tWriter->eleEta) = eleEta->at(iBestEle);
  (tWriter->elePhi) = elePhi->at(iBestEle);
  
  (tWriter->eleCharge) = eleCharge->at(iBestEle);
  
  (tWriter->eleIdGen) = idGenEle;
  (tWriter->eleBMot) = genEleBMot;
  
  (tWriter->eleIDNIV2Val) = eleIDNIV2Val;
  (tWriter->eleIDIV2Val) = eleIDIV2Val;
  (tWriter->eleIDHZZV1Val) = eleIDHZZV1Val;
  
  (tWriter->eleIDNIV2RawVal) = eleIDNIV2RawVal;
  (tWriter->eleIDIV2RawVal) = eleIDIV2RawVal;
  (tWriter->eleIDHZZV1RawVal) = eleIDHZZV1RawVal;
  
  (tWriter->eleIDNIV2Cat) = eleIDNIV2Cat;
  (tWriter->eleIDIV2Cat) = eleIDIV2Cat;
  (tWriter->eleIDHZZV1Cat) = eleIDHZZV1Cat;
  
  (tWriter->eleDx) = dXEle(iBestEle, iBestPV);
  (tWriter->eleDy) = dYEle(iBestEle, iBestPV);
  (tWriter->eleDxy) = dSignEle(iBestEle, pConeClean.Px(), pConeClean.Py())*dXYEle(iBestEle, iBestPV);
  (tWriter->eleDz) = dZEle(iBestEle, iBestPV);
  (tWriter->eleExy) = eleGsfExy->at(iBestEle);
  (tWriter->eleEz) = eleGsfEz->at(iBestEle);
  
  (tWriter->eleDRB) = deltaR(pB.Eta(), pB.Phi(), eleEta->at(iBestEle), elePhi->at(iBestEle));
  (tWriter->elePFIsoScaled) = GetEleRelPFIsoScaled(iBestEle);
  
  (tWriter->eleConePt) = eleConePt;
  (tWriter->eleConePtRel) = eleConePtRel;
  (tWriter->eleConeDR) = eleConeDR;
  (tWriter->eleConeEnergyRatio) = eleConeEnergyRatio;
  (tWriter->eleConeQ) = eleConeQ;
  (tWriter->eleConeSize) = eleConeSize;
  (tWriter->eleConeNF) = eleConeNF;
  (tWriter->eleConeCF) = eleConeCF;
  (tWriter->eleConeNCH) = eleConeNCH;
  
  (tWriter->eleConeAvgDx) = eleConeAvgDx;
  (tWriter->eleConeStdDevDx) = eleConeStdDevDx;
  (tWriter->eleConeAvgDy) = eleConeAvgDy;
  (tWriter->eleConeStdDevDy) = eleConeStdDevDy;
  (tWriter->eleConeAvgDxy) = eleConeAvgDxy;
  (tWriter->eleConeStdDevDxy) = eleConeStdDevDxy;
  (tWriter->eleConeAvgDz) = eleConeAvgDz;
  (tWriter->eleConeStdDevDz) = eleConeStdDevDz;

  (tWriter->eleConeCleanPt) = eleConeCleanPt;
  (tWriter->eleConeCleanPtRel) = eleConeCleanPtRel;
  (tWriter->eleConeCleanDR) = eleConeCleanDR;
  (tWriter->eleConeCleanEnergyRatio) = eleConeCleanEnergyRatio;
  (tWriter->eleConeCleanQ) = eleConeCleanQ;
  (tWriter->eleConeCleanSize) = eleConeCleanSize;
  (tWriter->eleConeCleanNF) = eleConeCleanNF;
  (tWriter->eleConeCleanCF) = eleConeCleanCF;
  (tWriter->eleConeCleanNCH) = eleConeCleanNCH;
  
  (tWriter->eleConeCleanAvgDx) = eleConeCleanAvgDx;
  (tWriter->eleConeCleanStdDevDx) = eleConeCleanStdDevDx;
  (tWriter->eleConeCleanAvgDy) = eleConeCleanAvgDy;
  (tWriter->eleConeCleanStdDevDy) = eleConeCleanStdDevDy;
  (tWriter->eleConeCleanAvgDxy) = eleConeCleanAvgDxy;
  (tWriter->eleConeCleanStdDevDxy) = eleConeCleanStdDevDxy;
  (tWriter->eleConeCleanAvgDz) = eleConeCleanAvgDz;
  (tWriter->eleConeCleanStdDevDz) = eleConeCleanStdDevDz;

  // Tagging truth
  (tWriter->tagTruth) = tagTruth;
  (tWriter->chargeCorr) = chargeCorr;
  
  tWriter->fill();

  return true;
}



void EleMVASecondNtupleProducer::endJob() 
{
  // This runs after the event loop
  tWriter->close();
  
  autoSavedObject = cEleBTrkDistance = CreateCanvas("cEleBTrkDistance", "colz", false, false, false, hEleBTrkDistance);

  autoSavedObject = cWeightedBMass = CreateCanvas("cWeightedBMass", 0, 21, 1, false, false, hWeightedBMass);
  
  autoSavedObject = cEleConeDistance = CreateCanvas("cEleConeDistance", 0, 21, 1, false, true, hEleConeDistance);
  
  autoSavedObject = cEleConePtRel = CreateCanvas("cEleConePtRel", 0, 21, 1, false, true, hEleConePtRel);
    
  autoSavedObject = cEleConeCleanDistance = CreateCanvas("cEleConeCleanDistance", 0, 21, 1, false, true, hEleConeCleanDistance);
  
  autoSavedObject = cEleConeCleanPtRel = CreateCanvas("cEleConeCleanPtRel", 0, 21, 1, false, true, hEleConeCleanPtRel);
  
  return;
}
