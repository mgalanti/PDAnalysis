#define UTIL_USE FULL



#include "PDAnalysis/Ntu/interface/constants.h"
#include "PDAnalysis/Ntu/interface/MGGenTools.h"



MGGenTools::MGGenTools()
{
  setUserParameter( "use_gen"      , "true" );
}



const int MGGenTools::GetClosestGen(const double pt, const double eta, const double phi, double dRMax, double dPtMax)
{
//   std::cout << "MGGenTools::GetClosestGen(...): pt = " << pt << ", eta = " << eta << ", phi = " << phi << std::endl;
//   std::cout << "           Matching parameters: dRMax = " << dRMax << ", dPtMax = " << dPtMax << std::endl;
  int matched = -1;
  double dRMatch = dRMax;
  double dPtMatch = dPtMax;
  
  for(uint iGen = 0; iGen < genId->size(); iGen++)
  {
    if(!IsLongLived(iGen))
    {
      continue;
    }
//     std::cout << "    Gen particle " << iGen << " is long-lived! gen ID = " << genId->at(iGen) << std::endl;
    float dR = deltaR(eta, phi, genEta->at(iGen), genPhi->at(iGen));
    float dPt = fabs(genPt->at(iGen) - pt)/genPt->at(iGen);
    
//     std::cout << "       dR = " << dR << ", dPt = " << dPt << std::endl;
    
    if( dR > dRMatch ) continue;
//     std::cout << "            dR  inside matching window!\n";
    if( dPt > dPtMatch) continue;
//     std::cout << "            dPt inside matching window!\n";
    
    matched = (int)iGen;
//     std::cout << "       matched = " << matched << std::endl;
    dRMatch = dR;
//     std::cout << "       dRMatch = " << dRMatch << std::endl;
  } 
  
  return matched;
}



const int MGGenTools::GetClosestGenNoLL(const double pt, const double eta, const double phi, double dRMax, double dPtMax)
{
//   std::cout << "MGGenTools::GetClosestGenNoLL(...): pt = " << pt << ", eta = " << eta << ", phi = " << phi << std::endl;
//   std::cout << "               Matching parameters: dRMax = " << dRMax << ", dPtMax = " << dPtMax << std::endl;
  int matched = -1;
  double dRMatch = dRMax;
  double dPtMatch = dPtMax;
  
  for(uint iGen = 0; iGen < genId->size(); iGen++)
  {
    float dR = deltaR(eta, phi, genEta->at(iGen), genPhi->at(iGen));
    float dPt = fabs(genPt->at(iGen) - pt)/genPt->at(iGen);
    
//     std::cout << "       dR = " << dR << ", dPt = " << dPt << std::endl;
    
    if( dR > dRMatch ) continue;
//     std::cout << "            dR  inside matching window!\n";
    if( dPt > dPtMatch) continue;
//     std::cout << "            dPt inside matching window!\n";
    
    matched = (int)iGen;
//     std::cout << "       matched = " << matched << std::endl;
    dRMatch = dR;
//     std::cout << "       dRMatch = " << dRMatch << std::endl;
  } 
  
  return matched;
}



const int MGGenTools::GetClosestGenNoLL(const int iSvt, double dRMax, double dPtMax)
{
  if(iSvt < 0 || iSvt >= nSVertices)
  {
    std::cout << "E R R O R ! MGGenTools::GetClosestGenNoLL(const int iSvt, double dRMax, double dPtMax): iSvt outside bounds!\n";
    std::cout << "            Exiting...\n";
    exit(1);
  }
  int type = svtType->at(iSvt);
  TLorentzVector pSvt(0, 0, 0, 0);
  if(type == PDEnumString::svtBsJPsiPhi)
  {
    std::vector<int> tkSsB = tracksFromSV(iSvt);
    int iJPsi = (subVtxFromSV(iSvt)).at(0);
    int iPhi = (subVtxFromSV(iSvt)).at(1);
    std::vector<int> tkJpsi = tracksFromSV(iJPsi);
    std::vector<int> tkPhi = tracksFromSV(iPhi);
    for(uint i = 0; i < tkSsB.size(); i++)
    {
      int j = tkSsB[i];
      float m = constants::kaonMass;
      if(j == tkJpsi[0] || j == tkJpsi[1])
      {
        m = constants::muonMass;
      }
      TLorentzVector a;
      a.SetPtEtaPhiM(trkPt->at(j), trkEta->at(j), trkPhi->at(j), m);
      pSvt += a;
    }
    return GetClosestGenNoLL(pSvt.Pt(), pSvt.Eta(), pSvt.Phi(), dRMax, dPtMax);
  }
  else
  {
    std::cout << "W A R N I N G ! MGGenTools::GetClosestGenNoLL(const int iSvt, double dRMax, double dPtMax): method only works for svt's of type:\n";
    std::cout << "                PDEnumString::svtBsJPsiPhi (i.e. " << PDEnumString::svtBsJPsiPhi <<")\n";
    std::cout << "                svt n. " << iSvt << " is of type " << type << "\n";
    std::cout << "                Returning -1\n";
    return -1;
  }
}



const int MGGenTools::GetClosestGenInList(const double pt, const double eta, const double phi, const std::vector<int>& listGenP, double dRMax, double dPtMax)
{
//   std::cout << "MGGenTools::GetClosestGenInList(...): pt = " << pt << ", eta = " << eta << ", phi = " << phi << std::endl;
//   std::cout << "                 Matching parameters: dRMax = " << dRMax << ", dPtMax = " << dPtMax << std::endl;
  int matched = -1;
  double dRMatch = dRMax;
  double dPtMatch = dPtMax;
  
  for(uint i = 0; i < listGenP.size(); i++)
  {
    uint iGen = listGenP[i];
    float dR = deltaR(eta, phi, genEta->at(iGen), genPhi->at(iGen));
    float dPt = fabs(genPt->at(iGen) - pt)/genPt->at(iGen);
    
//     std::cout << "       dR = " << dR << ", dPt = " << dPt << std::endl;
    
    if( dR > dRMatch ) continue;
//     std::cout << "            dR  inside matching window!\n";
    if( dPt > dPtMatch) continue;
//     std::cout << "            dPt inside matching window!\n";
    
    matched = (int)iGen;
//     std::cout << "       matched = " << matched << std::endl;
    dRMatch = dR;
//     std::cout << "       dRMatch = " << dRMatch << std::endl;
  } 
  
  return matched;
}



const int MGGenTools::GetGenLepBsChargeCorrelation(const int iGenLep, int iGenBs)
{
  if(iGenLep < 0 || iGenLep >= nGenP)
  {
    std::cout << "W A R N I N G ! MGGenTools::GetGenLepBsChargeCorrelation(...): iGenLep outside bounds!\n";
    std::cout << "                Returning 0...\n";
    return 0;
  }
  int idLep = genId->at(iGenLep);
  if(abs(idLep) != 11 && abs(idLep) != 13)
  {
    std::cout << "W A R N I N G ! MGGenTools::GetGenLepBsChargeCorrelation(...): Id of gen particle at index " << iGenLep << " is not 11 or 13!\n";
    std::cout << "                GenpId->at(" << iGenLep << ") = " << idLep << ".\n";
    std::cout << "                Returning 0...\n";
    return 0;
  }
  if(iGenBs < 0 || iGenBs >= nGenP)
  {
    std::cout << "W A R N I N G ! MGGenTools::GetGenLepBsChargeCorrelation(...): iGenBs outside bounds!\n";
    std::cout << "                Returning 0...\n";
    return 0;
  }
  int idBs = genId->at(iGenBs);
  if(abs(idBs) != 531 && abs(idBs) != 533 && abs(idBs) != 10533 && abs(idBs) != 20533)
  {
    std::cout << "W A R N I N G ! MGGenTools::GetGenLepBsChargeCorrelation(...): Id of gen particle at index " << iGenBs << " is not one of: 531, 533, 10533, 20533!\n";
    std::cout << "                GenpId->at(" << iGenBs << ") = " << idBs << ".\n";
    std::cout << "                Returning 0...\n";
    return 0;
  }
  // If we have a Bs, we make sure that it is the one before mixing
  if(abs(idBs) == 531)
  {
    int newIGenBs = RecursiveLookForMotherIds(iGenBs, {-idBs});
    if(newIGenBs >= 0)
    {
      std::cout << "I N F O : MGGenTools::GetGenLepBsChargeCorrelation(...): The Bs meson at index " << iGenBs << " with id " << idBs << " comes from mixing!\n";
      std::cout << "          Found pre-mixing meson at index " << newIGenBs;
      iGenBs = newIGenBs;
      idBs = genId->at(iGenBs);
      std::cout << " with id " << idBs << std::endl;
    }
  }
  // If we have a Bs*, we find the first Bs among the daughters
  if(abs(idBs) > 531)
  {
    int oldIGenBs = iGenBs;
    int oldIdBs = idBs;
    iGenBs = RecursiveLookForDaughterIds(iGenBs, {531}, false);
    if(iGenBs < 0)
    {
      std::cout << "E R R O R ! MGGenTools::GetGenLepBsChargeCorrelation(...): Found a Bs* which does not decay into a Bs!\n";
      std::cout << "            Old iGenBs = " << oldIGenBs << ", old idBs = " << oldIdBs << "\n";
      std::cout << "            New iGenBs = " << iGenBs;
      RecursivePrintDaughters(oldIGenBs, 0, "            ");
      std::cout << "            Exiting...\n";
      exit(1);
    }
    idBs = genId->at(iGenBs);
  }
  std::cout << "iGenBs = " << iGenBs << ", idBs = " << idBs << std::endl;
//   PrintMotherChain(iGenBs);
//   RecursivePrintMothers(iGenBs);
//   int iBQuark = RecursiveLookForMotherIds(iGenBs, {5}, false);
//   std::cout << "iBQuark = " << iBQuark << std::endl;
//   int idBQuark = genId->at(iBQuark);
  if((idLep > 0 && idBs < 0) || (idLep < 0 && idBs > 0))
  {
    return -1;
  }
  else
  {
    return 1;
  }
}



const int MGGenTools::GetGenLepBuChargeCorrelation(const int iGenLep, int iGenBu)
{
  if(iGenLep < 0 || iGenLep >= nGenP)
  {
    std::cout << "W A R N I N G ! MGGenTools::GetGenLepBuChargeCorrelation(...): iGenLep outside bounds!\n";
    std::cout << "                Returning 0...\n";
    return 0;
  }
  int idLep = genId->at(iGenLep);
  if(abs(idLep) != 11 && abs(idLep) != 13)
  {
    std::cout << "W A R N I N G ! MGGenTools::GetGenLepBuChargeCorrelation(...): Id of gen particle at index " << iGenLep << " is not 11 or 13!\n";
    std::cout << "                GenpId->at(" << iGenLep << ") = " << idLep << ".\n";
    std::cout << "                Returning 0...\n";
    return 0;
  }
  if(iGenBu < 0 || iGenBu >= nGenP)
  {
    std::cout << "W A R N I N G ! MGGenTools::GetGenLepBuChargeCorrelation(...): iGenBu outside bounds!\n";
    std::cout << "                Returning 0...\n";
    return 0;
  }
  int idBu = genId->at(iGenBu);
  if(abs(idBu) != 511 && abs(idBu) != 513 && abs(idBu) != 515 && abs(idBu) != 10513 && abs(idBu) != 20513)
  {
    std::cout << "W A R N I N G ! MGGenTools::GetGenLepBuChargeCorrelation(...): Id of gen particle at index " << iGenBu << " is not one of: 511, 513, 515, 10513, 20513!\n";
    std::cout << "                GenpId->at(" << iGenBu << ") = " << idBu << ".\n";
    std::cout << "                Returning 0...\n";
    return 0;
  }
  // If we have a Bu*, we find the first Bu among the daughters
  if(abs(idBu) > 511)
  {
    int oldIGenBu = iGenBu;
    int oldIdBu = idBu;
    iGenBu = RecursiveLookForDaughterIds(iGenBu, {511}, false);
    if(iGenBu < 0)
    {
      std::cout << "E R R O R ! MGGenTools::GetGenLepBuChargeCorrelation(...): Found a Bu* which does not decay into a Bu!\n";
      std::cout << "            Old iGenBu = " << oldIGenBu << ", old idBu = " << oldIdBu << "\n";
      std::cout << "            New iGenBu = " << iGenBu;
      RecursivePrintDaughters(oldIGenBu, 0, "            ");
      std::cout << "            Exiting...\n";
      exit(1);
    }
    idBu = genId->at(iGenBu);
  }
  std::cout << "iGenBu = " << iGenBu << ", idBu = " << idBu << std::endl;
//   PrintMotherChain(iGenBu);
//   RecursivePrintMothers(iGenBu);
//   int iBQuark = RecursiveLookForMotherIds(iGenBu, {5}, false);
//   std::cout << "iBQuark = " << iBQuark << std::endl;
//   int idBQuark = genId->at(iBQuark);
  if((idLep > 0 && idBu < 0) || (idLep < 0 && idBu > 0))
  {
    return -1;
  }
  else
  {
    return 1;
  }
}



const int MGGenTools::GetBsChargeCorrelation(const int charge, int iGenBs)
{
  if(abs(charge) != 1)
  {
    std::cout << "W A R N I N G ! MGGenTools::GetBsChargeCorrelation(...): charge is not +/-1!\n";
    std::cout << "                Returning 0...\n";
    return 0;
  }
  if(iGenBs < 0 || iGenBs >= nGenP)
  {
    std::cout << "W A R N I N G ! MGGenTools::GetBsChargeCorrelation(...): iGenBs outside bounds!\n";
    std::cout << "                Returning 0...\n";
    return 0;
  }
  int idBs = genId->at(iGenBs);
  if(abs(idBs) != 531 && abs(idBs) != 533 && abs(idBs) != 10533 && abs(idBs) != 20533)
  {
    std::cout << "W A R N I N G ! MGGenTools::GetBsChargeCorrelation(...): Id of gen particle at index " << iGenBs << " is not one of: 531, 533, 10533, 20533!\n";
    std::cout << "                GenpId->at(" << iGenBs << ") = " << idBs << ".\n";
    std::cout << "                Returning 0...\n";
    return 0;
  }
  // If we have a Bs, we make sure that it is the one before mixing
  if(abs(idBs) == 531)
  {
    int newIGenBs = RecursiveLookForMotherIds(iGenBs, {-idBs});
    if(newIGenBs >= 0)
    {
      std::cout << "I N F O : MGGenTools::GetBsChargeCorrelation(...): The Bs meson at index " << iGenBs << " with id " << idBs << " comes from mixing!\n";
      std::cout << "          Found pre-mixing meson at index " << newIGenBs;
      iGenBs = newIGenBs;
      idBs = genId->at(iGenBs);
      std::cout << " with id " << idBs << std::endl;
    }
  }
  // If we have a Bs*, we find the first Bs among the daughters
  if(abs(idBs) == 533)
  {
    int oldIGenBs = iGenBs;
    int oldIdBs = idBs;
    iGenBs = RecursiveLookForDaughterIds(iGenBs, {531}, false);
    if(iGenBs < 0)
    {
      std::cout << "E R R O R ! MGGenTools::GetBsChargeCorrelation(...): Found a Bs* which does not decay into a Bs!\n";
      std::cout << "            Old iGenBs = " << oldIGenBs << ", old idBs = " << oldIdBs << "\n";
      std::cout << "            New iGenBs = " << iGenBs;
      RecursivePrintDaughters(oldIGenBs, 0, "            ");
      std::cout << "            Exiting...\n";
      exit(1);
    }
    idBs = genId->at(iGenBs);
  }
  std::cout << "iGenBs = " << iGenBs << ", idBs = " << idBs << std::endl;
//   PrintMotherChain(iGenBs);
//   RecursivePrintMothers(iGenBs);
//   int iBQuark = RecursiveLookForMotherIds(iGenBs, {5}, false);
//   std::cout << "iBQuark = " << iBQuark << std::endl;
//   int idBQuark = genId->at(iBQuark);
  if((charge < 0 && idBs < 0) || (charge > 0 && idBs > 0))
  {
    return -1;
  }
  else
  {
    return 1;
  }
}



const int MGGenTools::GetBuChargeCorrelation(const int charge, int iGenBu)
{
  if(abs(charge) != 1)
  {
    std::cout << "W A R N I N G ! MGGenTools::GetBuChargeCorrelation(...): charge is not +/-1!\n";
    std::cout << "                Returning 0...\n";
    return 0;
  }
  if(iGenBu < 0 || iGenBu >= nGenP)
  {
    std::cout << "W A R N I N G ! MGGenTools::GetBuChargeCorrelation(...): iGenBu outside bounds!\n";
    std::cout << "                Returning 0...\n";
    return 0;
  }
  int idBu = genId->at(iGenBu);
  if(abs(idBu) != 511 && abs(idBu) != 513 && abs(idBu) != 515 && abs(idBu) != 10513 && abs(idBu) != 20513)
  {
    std::cout << "W A R N I N G ! MGGenTools::GetBuChargeCorrelation(...): Id of gen particle at index " << iGenBu << " is not one of: 511, 513, 515, 10513, 20513!\n";
    std::cout << "                GenpId->at(" << iGenBu << ") = " << idBu << ".\n";
    std::cout << "                Returning 0...\n";
    return 0;
  }
  // If we have a Bu*, we find the first Bu among the daughters
  if(abs(idBu) > 511)
  {
    int oldIGenBu = iGenBu;
    int oldIdBu = idBu;
    iGenBu = RecursiveLookForDaughterIds(iGenBu, {511}, false);
    if(iGenBu < 0)
    {
      std::cout << "E R R O R ! MGGenTools::GetBuChargeCorrelation(...): Found a Bu* which does not decay into a Bu!\n";
      std::cout << "            Old iGenBu = " << oldIGenBu << ", old idBu = " << oldIdBu << "\n";
      std::cout << "            New iGenBu = " << iGenBu;
      RecursivePrintDaughters(oldIGenBu, 0, "            ");
      std::cout << "            Exiting...\n";
      exit(1);
    }
    idBu = genId->at(iGenBu);
  }
  std::cout << "iGenBu = " << iGenBu << ", idBu = " << idBu << std::endl;
//   PrintMotherChain(iGenBu);
//   RecursivePrintMothers(iGenBu);
//   int iBQuark = RecursiveLookForMotherIds(iGenBu, {5}, false);
//   std::cout << "iBQuark = " << iBQuark << std::endl;
//   int idBQuark = genId->at(iBQuark);
  if((charge < 0 && idBu < 0) || (charge > 0 && idBu > 0))
  {
    return -1;
  }
  else
  {
    return 1;
  }
}



const int MGGenTools::GetMixStatus(const uint iGen)
{
  int pId = genId->at(iGen);
  
  if(RecursiveLookForDaughterIds(iGen, {-pId}) >= 0)
    return 2;

  if(RecursiveLookForMotherIds(iGen, {-pId}) >= 0)
    return 1;
    
  return 0;  
}



const std::vector<int> MGGenTools::GetAllGenElectrons()
{
  std::vector<int> viGenElectrons;
  for(uint iGen = 0; iGen < genId->size(); iGen++)
  {
    if(abs(genId->at(iGen)) == 11)
    {
      viGenElectrons.push_back(iGen);
    }
  }
  return viGenElectrons;
}



const std::vector<int> MGGenTools::GetAllGenElectronsFromB()
{
  std::vector<int> viGenElectrons;
  for(uint iGen = 0; iGen < genId->size(); iGen++)
  {
    if(abs(genId->at(iGen)) == 11)
    {
      // MG: FIXME: currently excluding bottomonium. 
      if(RecursiveLookForMotherIds(iGen, listBMesonsAndBaryons, false) > -1)
      {
        viGenElectrons.push_back(iGen);
      }
    }
  }
  return viGenElectrons;
}



const std::vector<int> MGGenTools::GetAllLongLivedBHadrons()
{
  std::vector<int> longLivedBHadrons;
  for(uint iGen = 0; iGen < genId->size(); iGen++)
  {
    int pId = abs(genId->at(iGen));
    for(auto longLivedId: listLongLivedBHadrons)
    {
      if(pId == longLivedId)
      {
        longLivedBHadrons.push_back(iGen);
        break;
      }
    }
  }
  return longLivedBHadrons;
}



const void MGGenTools::PrintMotherChain(const int iGen) // Adapted from Alberto Bragagnolo
{
  std::cout << genId->at(iGen) << " << ";
  const std::vector <int>& vMothers = allMothers(iGen);
  uint nMot = vMothers.size();
  if(nMot > 1)
  {
    for(uint i = 0; i < nMot; ++i)
    { 
      int iMot = vMothers[i];
      if(genId->at(iMot) != 21)
      {
        std::cout << "  iMot = " << iMot << ", Id = " << genId->at(iMot) << "\n";
      }
    }
  }
  
  if(nMot == 1)
  {
    PrintMotherChain(vMothers[0]);
  }
  return;
}



const void MGGenTools::RecursivePrintMothers(const unsigned short iGen, int recursionOrder, std::string prepend)
{
  const std::vector <int>& vMothers = allMothers(iGen);
  uint nMot = vMothers.size();
  if(recursionOrder > 99)
  {
    std::cout << "MGGenTools::RecursivePrintMothers(): E R R O R ! Too many recursions!\n";
    std::cout << "                                     recursionOrder = " << recursionOrder << ", genp " << iGen << " has " << nMot << " mothers.\n";
    return;
  }
  for(unsigned int i = 0; i < nMot; i++)
  {
    std::cout << prepend;
//     std::cout << "recursionOrder = " << recursionOrder << " ";
    std::cout << "(" << iGen << ")" << genId->at(iGen) << "<--" << "(" << vMothers[i] << ")" << genId->at(vMothers[i]);
    std::cout << " (charge " << genCharge->at(vMothers[i]) << ")";
    std::cout << " at (" << genVx->at(iGen) << "," << genVy->at(iGen) << "," << genVz->at(iGen) << ")";
    std::cout << " pt,eta,phi,M (" << genPt->at(iGen) << "," << genEta->at(iGen) << "," << genPhi->at(iGen) << "," << genMass->at(iGen) << ")\n";
    // Protection against infinite loops...
    if(vMothers[i] < iGen)
    {
      RecursivePrintMothers(vMothers[i], recursionOrder + 1, prepend + "    ");
    }
  }
}



const void MGGenTools::RecursivePrintDaughters(const unsigned short iGen, int recursionOrder, std::string prepend)
{
  const std::vector <int>& vDaughters = allDaughters(iGen);
  uint nDau = vDaughters.size();
  if(recursionOrder > 99)
  {
    std::cout << "MGGenTools::RecursivePrintDaughters(): E R R O R ! Too many recursions!\n";
    std::cout << "                                       recursionOrder = " << recursionOrder << ", genp " << iGen << " has " << nDau << " daughters.\n";
    return;
  }
  for(unsigned int i = 0; i < nDau; i++)
  {
    std::cout << prepend;
//     std::cout << "recursionOrder = " << recursionOrder << " ";
    std::cout << "(" << iGen << ")" << genId->at(iGen) << "<--" << "(" << vDaughters[i] << ")" << genId->at(vDaughters[i]);
    std::cout << " (charge " << genCharge->at(vDaughters[i]) << ")";
    std::cout << " at (" << genVx->at(iGen) << "," << genVy->at(iGen) << "," << genVz->at(iGen) << ")";
    std::cout << " pt,eta,phi,M (" << genPt->at(iGen) << "," << genEta->at(iGen) << "," << genPhi->at(iGen) << "," << genMass->at(iGen) << ")\n";
    // Protection against infinite loops...
    if(vDaughters[i] > iGen)
    {
      RecursivePrintDaughters(vDaughters[i], recursionOrder + 1, prepend + "    ");
    }
  }
}



const int MGGenTools::RecursiveLookForMotherIds(const unsigned short iGen, const std::vector<int> vIdList, const bool withSign, const int recursionOrder)
{
  const std::vector <int>& vMothers = allMothers(iGen);
  uint nMot = vMothers.size();
  if(recursionOrder == 0 && nMot == 0)
  {
    std::cout << "MGGenTools::RecursiveLookForMotherIds(...): W A R N I N G! Gen particle " << iGen << " has no mothers!\n";
    std::cout << "                                            Returning -65535.\n";
    return -65535;
  }
  if(recursionOrder > 99)
  {
    std::cout << "MGGenTools::RecursiveLookForMotherIds(...): E R R O R ! Too many recursions!\n";
    std::cout << "                                            recursionOrder = " << recursionOrder << ", genp " << iGen << " has " << nMot << " mothers.\n";
    std::cout << "                                            Returning -65535.\n";
    return -65535;
  }
//   bool found = false;
//   int iFound = -65535;
  for(unsigned int i = 0; i < nMot; i++)
  {
    unsigned short iMot =vMothers[i];
    int idMom = genId->at(iMot);
    if(!withSign)
      idMom = abs(idMom);
    for(unsigned int j = 0; j < vIdList.size(); j++)
    {
      int id = vIdList[j];
      if(!withSign)
        id = abs(id);
      if(idMom == id)
      {
        return iMot;
      }
    }
    // Protection against infinite loops...
    if(vMothers[i] < iGen)
    {
      return RecursiveLookForMotherIds(vMothers[i], vIdList, withSign, recursionOrder + 1);
    }
  }
  if(recursionOrder == 0)
  {
    std::cout << "MGGenTools::RecursiveLookForMotherIds(...): W A R N I N G! Reached default return statement! Gen particle = " << iGen << "\n";
    std::cout << "                                            Returning -65535.\n";
  }
  return -65535;
}



const int MGGenTools::RecursiveLookForDaughterIds(const unsigned short iGen, const std::vector<int> vIdList, const bool withSign, const int recursionOrder)
{
  const std::vector <int>& vDaughters = allDaughters(iGen);
  uint nDau = vDaughters.size();
  if(recursionOrder == 0 && nDau == 0)
  {
    std::cout << "MGGenTools::RecursiveLookForDaughterIds(...): W A R N I N G! Gen particle " << iGen << " has no daughters!\n";
    std::cout << "                                              Returning -65535.\n";
    return -65535;
  }
  if(recursionOrder > 99)
  {
    std::cout << "MGGenTools::RecursiveLookForDaughterIds(...): E R R O R ! Too many recursions!\n";
    std::cout << "                                              recursionOrder = " << recursionOrder << ", genp " << iGen << " has " << nDau << " daughters.\n";
    std::cout << "                                              Returning -65535.\n";
    return -65535;
  }
//   bool found = false;
//   int iFound = -65535;
  for(unsigned int i = 0; i < nDau; i++)
  {
    unsigned short iDau =vDaughters[i];
    int idDau = genId->at(iDau);
    if(!withSign)
      idDau = abs(idDau);
    for(unsigned int j = 0; j < vIdList.size(); j++)
    {
      int id = vIdList[j];
      if(!withSign)
        id = abs(id);
      if(idDau == id)
      {
        return iDau;
      }
    }
    // Protection against infinite loops...
    if(vDaughters[i] > iGen)
    {
      return RecursiveLookForDaughterIds(vDaughters[i], vIdList, withSign, recursionOrder + 1);
    }
  }
  if(recursionOrder == 0)
  {
    std::cout << "MGGenTools::RecursiveLookForDaughterIds(...): W A R N I N G! Reached default return statement! Gen particle = " << iGen << "\n";
    std::cout << "                                              Returning -65535.\n";
  }
  return -65535;
}



const bool MGGenTools::IsLongLived(const uint iGen) 
{
  if((int)iGen >= nGenP)
  {
    std::cout << "MGGenTools::IsLongLived(iGen): E R R O R ! Attempt to access a gen particle beyond the end of the collection!\n";
    std::cout << "                               iGen = " << iGen << ", nGenP = " << nGenP << std::endl;
    std::cout << "                               This is a bug. Please check the calling code and fix it.";
    std::cout << "                               Exiting...\n";
    exit(1);
  }
  int genCode = abs(genId->at(iGen));
  for(uint i = 0; i < listLongLived.size(); i++)
  {
    if(genCode == listLongLived[i])
      return true;
  }
  return false;
}
