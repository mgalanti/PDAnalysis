#define UTIL_USE FULL



#include "PDAnalysis/OSElectronFlavorTagger/interface/OsElectronMvaTagVariablesComputer.h"



OsElectronMvaTagVariablesComputer::OsElectronMvaTagVariablesComputer()
{
  std::cout << "New OsElectronMvaTagVariablesComputer.\n";
  iEle_ = -1;
  iB_ = -1;
  iPV_ = -1;
  
  isInitialized_ = false;
  evtInitialized_ = -1;
}


void OsElectronMvaTagVariablesComputer::computeOsElectronMvaTagVariables()
{
  if(!isInitialized_)
  {
    std::cout << "E R R O R ! OsElectronMvaTagVariablesComputer::computeOsElectronMvaTagVariables(): method called but object indexes are not initialized!\n";
    std::cout << "            Please initialize object indexes before computing Mva variables!\n";
    std::cout << "            Exiting...\n";
    exit(1);
  }
  if(evtInitialized_ != currentEvt)
  {
    std::cout << "E R R O R ! OsElectronMvaTagVariablesComputer::computeOsElectronMvaTagVariables(): object indexes are initialized for a different event!\n";
    std::cout << "            Event inizialized = " << evtInitialized_ << ", current event = " << currentEvt << ".\n";
    std::cout << "            Please fix calling code.\n";
    std::cout << "            Exiting...\n";
    exit(1);
  }
  
  // Compute B features - needed below
  std::vector<int> tracksFromB = tracksFromSV(iB_);
  TLorentzVector pB = GetTLorentzVectorFromJPsiX(iB_);
  
  // Electron kinematics
  TLorentzVector pEle;
  pEle.SetPtEtaPhiM(eleTagVars.elePt, eleTagVars.eleEta, eleTagVars.elePhi, constants::electronMass);

  eleTagVars.elePt = elePt->at(iEle_);
  eleTagVars.eleEta = eleEta->at(iEle_);
  eleTagVars.elePhi = elePhi->at(iEle_);
  eleTagVars.eleExy = eleGsfExy->at(iEle_);
  eleTagVars.eleDz = dZ(iEle_, iPV_);
  eleTagVars.eleEz = eleGsfEz->at(iEle_);
  
  // Get eleID Mva values
  for (int iUserInfo = 0; iUserInfo < nUserInfo; iUserInfo++)
  {
    //     std::cout << "Looping on userInfo object #" << iUserInfo << std::endl;
    if(useObjType->at(iUserInfo) == PDEnumString::recElectron && useObjIndex->at(iUserInfo) == iEle_)
    {
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17NoIsoV2RawValues)
      {
        eleTagVars.eleIDNIV2RawVal = useInfoValue->at(iUserInfo);
      }
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17IsoV2RawValues)
      {
        eleTagVars.eleIDIV2RawVal = useInfoValue->at(iUserInfo);
      }
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Spring16HZZV1RawValues)
      {
        eleTagVars.eleIDHZZV1RawVal = useInfoValue->at(iUserInfo);
      }
      
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17NoIsoV2Values)
      {
        eleTagVars.eleIDNIV2Val = useInfoValue->at(iUserInfo);
      }
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17IsoV2Values)
      {
        eleTagVars.eleIDIV2Val = useInfoValue->at(iUserInfo);
      }
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Spring16HZZV1Values)
      {
        eleTagVars.eleIDHZZV1Val = useInfoValue->at(iUserInfo);
      }
      
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17NoIsoV2Categories)
      {
        eleTagVars.eleIDNIV2Cat = (int)(useInfoValue->at(iUserInfo));
      }
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Fall17IsoV2Categories)
      {
        eleTagVars.eleIDIV2Cat = (int)(useInfoValue->at(iUserInfo));
      }
      if(useInfoType->at(iUserInfo) == PDEnumString::ElectronMVAEstimatorRun2Spring16HZZV1Categories)
      {
        eleTagVars.eleIDHZZV1Cat = (int)(useInfoValue->at(iUserInfo));
      }
    }
  }
  
  eleTagVars.eleDRB = deltaR(pB.Eta(), pB.Phi(), eleTagVars.eleEta, eleTagVars.elePhi);
  eleTagVars.elePFIsoScaled = GetEleRelPFIsoScaled(iEle_);
 
  // CONE variables
  float kappa = 1;
  float drCone = 0.4;

  eleTagVars.eleConePtRel = -1;
  eleTagVars.eleConeDr = -1;
  eleTagVars.eleConeEnergyRatio = -1;
  eleTagVars.eleConeSize = 0;
  eleTagVars.eleConeQ = -1;
  eleTagVars.eleConePt = -1;
  eleTagVars.eleConeNF = 0;
  eleTagVars.eleConeCF = 0;
  eleTagVars.eleConeNCH = 0;
    
  if(verbose)
  {
    std::cout << "I N F O: OsElectronMvaTagVariablesComputer::computeOsElectronMvaTagVariables() - Elec: eta = " << eleTagVars.eleEta << ", phi = " << eleTagVars.elePhi << ", pt = " << eleTagVars.elePt << std::endl;
  }
  
  float qCone = 0, ptCone = 0;
  TLorentzVector pCone(0.,0.,0.,0.);
  
  for(int iPF=0; iPF < nPF; ++iPF)
  {
    float ptPF = pfcPt->at(iPF);
    float etaPF = pfcEta->at(iPF);
    if(deltaR(etaPF, pfcPhi->at(iPF), eleTagVars.eleEta, eleTagVars.elePhi) > drCone)
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
      if(fabs(dZTrk(pfcTrk->at(iPF), iPV_)) >= 1.0)
      {
        continue;
      }
    }
    TLorentzVector a;
    a.SetPtEtaPhiE(pfcPt->at(iPF), pfcEta->at(iPF), pfcPhi->at(iPF), pfcE->at(iPF));
    pCone += a;
    ++eleTagVars.eleConeSize;
    
    qCone += pfcCharge->at(iPF) * pow(ptPF, kappa);
    ptCone += pow(ptPF, kappa);
    
    if(pfcCharge->at(iPF)==0)
    {
      eleTagVars.eleConeNF += pfcE->at(iPF);
    }
    if(abs(pfcCharge->at(iPF))==1)
    {
      eleTagVars.eleConeNCH++;
      eleTagVars.eleConeCF += pfcE->at(iPF);
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
  qCone *= eleCharge->at(iEle_);
  if(pCone.E() != 0)
  {
    eleTagVars.eleConeCF /= pCone.E();
    eleTagVars.eleConeNF /= pCone.E();
  }
  
  eleTagVars.eleConeQ = qCone;
  
  eleTagVars.eleConePt = pCone.Pt();
  eleTagVars.eleConeDr = deltaR(pCone.Eta(), pCone.Phi(), eleEta->at(iEle_), elePhi->at(iEle_));
  if(pCone.E() !=0)
  {
    eleTagVars.eleConeEnergyRatio = eleE->at(iEle_) / pCone.E();
  }
  else
  {
    eleTagVars.eleConeEnergyRatio = 1;
  }
  pCone -= pEle;
  eleTagVars.eleConePtRel = elePt->at(iEle_) * (pEle.Vect().Unit() * pCone.Vect().Unit());
  
  if(verbose)
  {
    std::cout << "I N F O: OsElectronMvaTagVariablesComputer::computeOsElectronMvaTagVariables():\n";
    std::cout << "         Cone: eta = " << pCone.Eta() << ", phi = " << pCone.Phi() << ", pt = " << pCone.Pt() <<  std::endl; 
    std::cout << "               eleConePt = " << eleTagVars.eleConePt << std::endl;
    std::cout << "               ptRel = " << eleTagVars.eleConePtRel << std::endl;
    std::cout << "               eleConeDr = " << eleTagVars.eleConeDr << std::endl;
    std::cout << "               eleConeEnergyRatio = " << eleTagVars.eleConeEnergyRatio << std::endl;
    std::cout << "               eleConeQ = " << eleTagVars.eleConeQ << std::endl;
  }
  
  
  // CONECLEAN variables
  kappa = 1;
  drCone = 0.4;

  eleTagVars.eleConeCleanPtRel = -1;
  eleTagVars.eleConeCleanDr = -1;
  eleTagVars.eleConeCleanEnergyRatio = -1;
  eleTagVars.eleConeCleanSize = 0;
  eleTagVars.eleConeCleanQ = -1;
  eleTagVars.eleConeCleanPt = -1;
  eleTagVars.eleConeCleanNF = 0;
  eleTagVars.eleConeCleanCF = 0;
  eleTagVars.eleConeCleanNCH = 0;
    
  float qConeClean = 0, ptConeClean = 0;
  TLorentzVector pConeClean(0.,0.,0.,0.);
  
  for(int iPF=0; iPF < nPF; ++iPF)
  {
    float ptPF = pfcPt->at(iPF);
    float etaPF = pfcEta->at(iPF);
    if(deltaR(etaPF, pfcPhi->at(iPF), eleTagVars.eleEta, eleTagVars.elePhi) > drCone)
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
    if(pfcTrk->at(iPF) < 0)
    {
      continue;
    }
    if(fabs(dZTrk(pfcTrk->at(iPF), iPV_)) >= 1.0)
    {
      continue;
    }
    TLorentzVector a;
    a.SetPtEtaPhiE(pfcPt->at(iPF), pfcEta->at(iPF), pfcPhi->at(iPF), pfcE->at(iPF));
    pConeClean += a;
    ++eleTagVars.eleConeCleanSize;
    
    qConeClean += pfcCharge->at(iPF) * pow(ptPF, kappa);
    ptConeClean += pow(ptPF, kappa);
    
    if(pfcCharge->at(iPF)==0)
    {
      eleTagVars.eleConeCleanNF += pfcE->at(iPF);
    }
    if(abs(pfcCharge->at(iPF))==1)
    {
      eleTagVars.eleConeCleanNCH++;
      eleTagVars.eleConeCleanCF += pfcE->at(iPF);
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
  qConeClean *= eleCharge->at(iEle_);
  if(pConeClean.E() != 0)
  {
    eleTagVars.eleConeCleanCF /= pConeClean.E();
    eleTagVars.eleConeCleanNF /= pConeClean.E();
  }
  
  eleTagVars.eleConeCleanQ = qConeClean;
  
  eleTagVars.eleConeCleanPt = pConeClean.Pt();
  eleTagVars.eleConeCleanDr = deltaR(pConeClean.Eta(), pConeClean.Phi(), eleEta->at(iEle_), elePhi->at(iEle_));
  if(pConeClean.E() !=0)
  {
    eleTagVars.eleConeCleanEnergyRatio = eleE->at(iEle_) / pConeClean.E();
  }
  else
  {
    eleTagVars.eleConeCleanEnergyRatio = 1;
  }
  pConeClean -= pEle;
  eleTagVars.eleConeCleanPtRel = elePt->at(iEle_) * (pEle.Vect().Unit() * pConeClean.Vect().Unit());
  pConeClean += pEle; // for IP sign
  
  if(verbose)
  {
    std::cout << "I N F O: OsElectronMvaTagVariablesComputer::computeOsElectronMvaTagVariables():\n";
    
    std::cout << "         ConeClean: eta = " << pConeClean.Eta() << ", phi = " << pConeClean.Phi() << ", pt = " << pConeClean.Pt() <<  std::endl; 
    std::cout << "                    eleConeCleanPt = " << eleTagVars.eleConeCleanPt << std::endl;
    std::cout << "                    ptRel = " << eleTagVars.eleConeCleanPtRel << std::endl;
    std::cout << "                    eleConeCleanDr = " << eleTagVars.eleConeCleanDr << std::endl;
    std::cout << "                    eleConeCleanEnergyRatio = " << eleTagVars.eleConeCleanEnergyRatio << std::endl;
    std::cout << "                    eleConeCleanQ = " << eleTagVars.eleConeCleanQ << std::endl;
  }
  
  // Computed at the end because it needs the coneClean
  eleTagVars.eleDxy = dSignEle(iEle_, pConeClean.Px(), pConeClean.Py())*abs(eleGsfDxy->at(iEle_));
}
