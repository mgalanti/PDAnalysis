#define UTIL_USE FULL

#include "PDAnalysis/Ntu/interface/constants.h"
#include "PDAnalysis/Ntu/interface/MGSelector.h"

#include <iostream>
#include <utility>



MGSelector::MGSelector()
{
  std::cout << "New MGSelector()...\n";
}



MGSelector::~MGSelector()
{
}



void MGSelector::SetSelectionString(const std::string str)
{
  selectionString = str;
  selectionSubStrings = TokenizeString(selectionString, "_");
  std::cout << "MGSelector::SetSelectionString(): I N F O. Selection strings set to:\n";
  std::cout << "                                  selectionString = " << selectionString << std::endl;
  for(unsigned int i = 0; i < selectionSubStrings.size(); i++)
    std::cout << "                                  selectionSubStrings[" << i << "] = " << selectionSubStrings[i] << std::endl;
}



bool MGSelector::SelectEvent()
{
  selectedObjects.clear();
  return SelectEvent(selectionString, selectedObjects);
}



bool MGSelector::SelectEvent(const std::string selection, std::vector<std::pair<int,int> >& selectedObjectsRef)
{
  // The first part of the selection string defines the kind of signal process to look at
  // The second part of the selection string defines the version of the selection
  std::vector<std::string> selectionSubStrings = TokenizeString(selection, "_");
  bool signalSelected = false;
  
  if(selectionSubStrings[0].compare("BsToJPsiPhi") == 0)
  {
    signalSelected = SelectBsToJPsiPhiEvent(selectionSubStrings[1], selectedObjectsRef);
  }
  else if(selectionSubStrings[0].compare("BuToJPsiK") == 0)
  {
    signalSelected = SelectBuToJPsiKEvent(selection.substr(10, std::string::npos), selectedObjectsRef);
  }
  
//   std::cout << "MGSelector::SelectEvent(): I N F O. Event " << (signalSelected?"passes":"does not pass") << " signal-side selection.\n";
  
  if(!signalSelected)
    return false;

  // If the second part of the selection string starts with "eleTag", then perform also the selection of the OS tag electron
  bool tagSideSelected = false;
  if(selectionSubStrings[1].substr(0,6).compare("eleTag") == 0)
  {
    int iPV = -1;
    int iB = -1;
    for(auto itSelObjects : selectedObjectsRef)
    {
      if(itSelObjects.first == PDEnumString::recPV)
        iPV = itSelObjects.second;
      if(itSelObjects.first == PDEnumString::recSvt)
        iB = itSelObjects.second;
    }
    int iOSElectron = SelectOSElectron(selectionSubStrings[1], iPV, iB);
    if(iOSElectron < 0)
      return false;
    std::pair<int,int> selObject = std::make_pair(PDEnumString::recElectron, iOSElectron);
    selectedObjectsRef.push_back(selObject);

    tagSideSelected = true;
    
//     std::cout << "MGSelector::SelectEvent(): I N F O. Event " << (tagSideSelected?"passes":"does not pass") << " opposite-side selection.\n";
    
    return signalSelected && tagSideSelected;
  }
  
  return signalSelected;
  
  std::cout << "MGSelector::SelectEvent(): W A R N I N G ! Reached default return statement. Event will not be selected...\n";
  std::cout << "                           WE SHOULD NEVER BE HERE! PLEASE CHECK THE CODE FOR BUGS!\n";
  std::cout << "                           The selection string passed to this method is: \"" << selection << "\"\n";
  return false;
}



bool MGSelector::SelectBsToJPsiPhiEvent(const std::string selection, std::vector<std::pair<int,int> >& selectedObjectsRef)
{
  if(selection.compare("eleTagLooseV0") == 0)
  {
    if(SelectHlt("eleTagHLTV0") == false)
      return false;
    int iPV = -1;
    int bestBsToJPsiPhi = SelectBestBsToJPsiPhi("LooseV0", iPV);
    if(bestBsToJPsiPhi < 0) return false;
    if(iPV < 0) return false;
    std::pair<int,int> selObject = std::make_pair(PDEnumString::recSvt, bestBsToJPsiPhi);
    selectedObjectsRef.push_back(selObject);
    selObject = std::make_pair(PDEnumString::recPV, iPV);
     selectedObjectsRef.push_back(selObject);
     return true;
  }
  else if(selection.compare("eleTagTightV1") == 0)
  {
    if(SelectHlt("eleTagHLTV0") == false)
      return false;
    int iPV = -1;
    int bestBsToJPsiPhi = SelectBestBsToJPsiPhi("TightV1", iPV);
    if(bestBsToJPsiPhi < 0) return false;
    if(iPV < 0) return false;
    std::pair<int,int> selObject = std::make_pair(PDEnumString::recSvt, bestBsToJPsiPhi);
    selectedObjectsRef.push_back(selObject);
    selObject = std::make_pair(PDEnumString::recPV, iPV);
    selectedObjectsRef.push_back(selObject);
    return true;
  }
  std::cout << "MGSelector::SelectBsToJPsiPhiEvent(): W A R N I N G ! Reached default return statement. Event will not be selected...\n";
  std::cout << "                                      The selection string passed to this method is: \"" << selection << "\"\n";
  return false;
}



bool MGSelector::SelectBuToJPsiKEvent(const std::string selection, std::vector<std::pair<int, int> >& selectedObjectsRef)
{
  // FIXME: To be written
  std::cout << "MGSelector::SelectBuToJPsiKEvent(): W A R N I N G ! Reached default return statement. Event will not be selected...\n";
  std::cout << "                                    The selection string passed to this method is: \"" << selection << "\"\n";
  return false;
}



bool MGSelector::SelectHlt(const std::string selection)
{
  bool jpsimu = false;
  bool jpsitktk = false;
  bool jpsitk = false;
  
  if(hlt(PDEnumString::HLT_Dimuon0_Jpsi3p5_Muon2_v) || hlt(PDEnumString::HLT_Dimuon0_Jpsi_Muon_v))
  {
    jpsimu = true;
  }
  if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrkTrk_Displaced_v))
  {
    jpsitktk = true;
  }
  if(hlt(PDEnumString::HLT_DoubleMu4_JpsiTrk_Displaced_v))
  {
    jpsitk = true;
  } 
  
  if(selection.compare("eleTagHLTV0") == 0)
  {
    return ((!jpsimu) && (jpsitktk || jpsitk));
  }
  std::cout << "MGSelector::SelectHlt(): W A R N I N G ! Reached default return statement. Event will not be selected...\n";
  return false;
}



bool MGSelector::SelectElectron(const int iEle, const std::string selection)
{
  std::cout << "MGSelector::SelectElectron(): W A R N I N G ! Reached default return statement. Electron will not be selected...\n";
  return false;
}



int MGSelector::SelectOSElectron(const std::string selection, const int iPV, const int iB)
{
  if (nElectrons == 0)
    return -1;
  
  if(selection.compare("eleTagTightV1") == 0 || selection.compare("eleTagLooseV0") == 0)
  {
    std::vector <int> tkSsB = tracksFromSV(iB);
    TLorentzVector tB = GetTLorentzVectorFromJPsiX(iB);
    
    int iOSElectron = -1;
    float ptOSElectron = 2.;
    float dZWrtPV;
    float pt, eta, phi;
    for(int iElectron = 0; iElectron < nElectrons; iElectron++)
    {
      pt = elePt->at(iElectron);
      if(pt < 2.)
        continue;
      eta = eleEta->at(iElectron);
      if(fabs(eta) > 2.4)
        continue;
      phi = elePhi->at(iElectron);
      if(deltaR(tB.Eta(), tB.Phi(), eta, phi) < 0.4) continue;
      dZWrtPV = dZ(iElectron, iPV);
      if(fabs(dZWrtPV) > 1)
        continue;
      // Make sure that the electron is not any of the tracks used to build the B
      // Since the electron track is not a standard track, must check by deltaR, deltaPt
      bool closeTrk = false;
      for(auto iTrk : tkSsB)
      {
        if(deltaR(trkEta->at(iTrk), trkPhi->at(iTrk), eleGsfEta->at(iElectron), eleGsfPhi->at(iElectron)) < 0.05)
        {
          std::cout << "MGSelector::SelectOSElectron(): I N F O. Gsf track of electron " << iElectron << " is close in dR to one of the B tracks.\n";
          std::cout << "                                         trkEta->at(" << iTrk << ") = " << trkEta->at(iTrk) << ", eleGsfEta->at(" << iElectron << ") = " << eleGsfEta->at(iElectron) << std::endl;
          std::cout << "                                         trkPhi->at(" << iTrk << ") = " << trkPhi->at(iTrk) << ", eleGsfPhi->at(" << iElectron << ") = " << eleGsfPhi->at(iElectron) << std::endl;
          std::cout << "                                         Electron will not be considered for OS tagging.\n";
          closeTrk = true;
        }
        if(2*fabs(trkPt->at(iTrk) - eleGsfPt->at(iElectron))/(trkPt->at(iTrk) + eleGsfPt->at(iElectron)) < 0.05)
        {
          std::cout << "MGSelector::SelectOSElectron(): I N F O. Gsf track of electron " << iElectron << " is close in dpT/pT to one of the B tracks.\n"; 
          std::cout << "                                         trkPt->at(" << iTrk << ") = " << trkPt->at(iTrk) << ", eleGsfPt->at(" << iElectron << ") = " << eleGsfPt->at(iElectron) << std::endl;
          std::cout << "                                         Electron will not be considered for OS tagging.\n";
          closeTrk = true;
        }
        if(closeTrk)
          break;
      }
      if(closeTrk)
        continue;
      if(pt > ptOSElectron)
      {
        iOSElectron = iElectron;
        ptOSElectron = pt;
      }
    }
    return iOSElectron;
  } 
  std::cout << "MGSelector::SelectOSElectron(): W A R N I N G ! Reached default return statement. No electrons will not be selected...\n";
  return -1;
}



bool MGSelector::SelectJPsi(const int iJPsi, const std::string selection)
{
  if(selection.compare("TightV0") == 0) // This is the same as IsTightJPsi() in Alberto code.
  {
    if(fabs(svtMass->at(iJPsi) - constants::jPsiMass) > 0.15)
    {
      return false;
    }
    
    std::vector<int> tkJPsi = tracksFromSV(iJPsi);
    TLorentzVector tJPsi(0,0,0,0);

    for(uint i=0; i < tkJPsi.size(); i++)
    {
      int j = tkJPsi[i];
      if(trkPt->at(j) < 3.6) // For JPsiTrktrk
      {
        return false;
      }
      if(fabs(trkEta->at(j)) > 2.4)
      {
        return false;
      }
      TLorentzVector a;
      a.SetPtEtaPhiM(trkPt->at(j), trkEta->at(j), trkPhi->at(j), constants::muonMass);
      tJPsi += a;
    }

    // if(tJPsi.Pt() < 7.0) return false;

    return true;
  }
  std::cout << "MGSelector::SelectJPsi(): W A R N I N G ! Reached default return statement. J/psi will not be selected...\n";
  return false;
}



bool MGSelector::SelectPhi(const int iPhi, const std::string selection)
{
  if(selection.compare("TightV0") == 0) // This is the same as IsTightJPsi() in Alberto code.
  {
    if(fabs(svtMass->at(iPhi) - constants::phiMass) > 0.01)
    {
      return false;
    }
    
    std::vector<int> tkPhi = tracksFromSV(iPhi);
    for(uint i=0; i < tkPhi.size(); i++)
    {
      int j = tkPhi[i];
      if(trkPt->at(j) < 1.0) // For JpsiTrkTrk
      {
        return false;
      }
      if(fabs(trkEta->at(j)) > 2.5)
      {
        return false;
      }

      // FIXME: MG: write a method inside MGRecoTools that does this automatically...
      int kaonHits = trkHitPattern->at(tkPhi[i]);
      kaonHits = (int(kaonHits) / 100) % 10000;
      kaonHits = kaonHits % 100;
      if(kaonHits < 4) return false;
    }
    return true;
  }
  std::cout << "MGSelector::SelectPhi(): W A R N I N G ! Reached default return statement. Phi will not be selected...\n";
  return false;  
}



int MGSelector::SelectBestBsToJPsiPhi(const std::string selection, int& iBestPV)
{
  iBestPV = -1;
  TLorentzVector tB(0,0,0,0);
  
  if(selection.compare("LooseV0") == 0) // This is the same as GetBestBstrange() in Alberto code.
  {
    int index = -1;
    float bestChi2 = 1e9;
    for(int iB = 0; iB < nSVertices; iB++)
    {
       if((svtType->at(iB) != PDEnumString::svtBsJPsiPhi))
       {
         continue;
       }
       if(svtMass->at(iB) < 5.0 || svtMass->at(iB) > 6.0)
       {
         continue;
       }
       if(svtChi2->at(iB) > bestChi2)
       {
         continue;
       }
       
       std::vector<int> tkSsB = tracksFromSV(iB);
       int iJPsi = (subVtxFromSV(iB)).at(0);
       std::vector<int> tkJpsi = tracksFromSV(iJPsi);

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
         tB += a;
       }
       iBestPV = SelectBestPV(iB, tB, "PointingV0");
       index = iB;
       bestChi2 = svtChi2->at(iB);
    }
    return index;
  }
  else if(selection.compare("TightV1") == 0) // This is the same as GetBestBstrangeTight() in Alberto code.
  {
    int index = -1;
    float best = 0.;
    
    for(int iB = 0; iB < nSVertices; iB++)
    {
      if((svtType->at(iB) != PDEnumString::svtBsJPsiPhi))
      {
        continue;
      }
      
      int iJPsi = (subVtxFromSV(iB)).at(0);
      int iPhi = (subVtxFromSV(iB)).at(1);
      
      std::vector<int> tkSsB = tracksFromSV(iB);
      std::vector<int> tkJpsi = tracksFromSV(iJPsi);
      std::vector<int> tkPhi = tracksFromSV(iPhi);
      
      //JPsi
      if(!SelectJPsi(iJPsi, "TightV0"))
      {
        continue;
      }
      //Phi
      if(!SelectPhi(iPhi, "TightV0"))
      {
        continue;
      }
      
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
        tB += a;
      }
      
      //BS
      float bVprob = ChiSquaredProbability(svtChi2->at(iB), svtNDOF->at(iB));
      if(svtMass->at(iB) < 5.2 || svtMass->at(iB) > 5.65)
      {
        continue;
      }
      if(bVprob < 0.001)
      {
        continue;
      }
      if(tB.Pt() < 10.)
      {
        continue;
      }
      
      int iPV = SelectBestPV(iB, tB, "PointingV0");
      
      if(iPV < 0)
      {
        continue;
      }
      
      if(GetCt2D(tB, iB) < 0.02)
      {
        continue;
      }
      
      if(bVprob < best)
      {
        continue;
      }
      iBestPV = iPV;
      index = iB;
      best = bVprob;
    }
    return index;
  }
  
  std::cout << "MGSelector::SelectBestBsToJPsiPhi(): W A R N I N G ! Reached default return statement. Best BsToJpsiPhi will not be selected...\n";
  return -1;
}



int MGSelector::SelectBestPV(const int iSvt, const TLorentzVector& pSvt, const std::string selection)
{
  if(selection.compare("PointingV0") == 0)
  {
    int iBestPV = -1;
    float bestCos = -1;
    
    TVector3 vB(pSvt.Px(),pSvt.Py(),pSvt.Pz());
    TVector3 vSVT(svtX->at(iSvt), svtY->at(iSvt), svtZ->at(iSvt));
    
    for(int iPV = 0; iPV < nPVertices; iPV++)
    {
      if(fabs(svtZ->at(iSvt) - pvtZ->at(iPV)) > 0.5)
      {
        continue;
      }
      
      TVector3 vPV(pvtX->at(iPV), pvtY->at(iPV), pvtZ->at(iPV));
      TVector3 vPointing;
      vPointing = vSVT - vPV;
      float cos = vPointing.Unit() * vB.Unit();
      
      if(cos > bestCos )
      {
        bestCos = cos;
        iBestPV = iPV;
      }      
    }
    
    return iBestPV;
  }

  std::cout << "MGSelector::SelectBestPV(): W A R N I N G ! Reached default return statement. Best PV will not be selected...\n";
  return -1;  
}
