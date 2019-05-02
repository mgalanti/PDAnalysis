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



bool MGSelector::SelectEvent(const std::string selection, std::vector<std::pair<int,int> >& selectedObjects)
{
  if(selection.substr(0, 11).compare("BsToJPsiPhi") == 0)
  {
    return SelectBsToJPsiPhiEvent(selection.substr(12, std::string::npos), selectedObjects);
  }
  else if(selection.substr(0, 9).compare("BuToJPsiK") == 0)
  {
    return SelectBuToJPsiKEvent(selection.substr(10, std::string::npos), selectedObjects);
  }
  std::cout << "MGSelector::SelectEvent(): W A R N I N G ! Reached default return statement. Event will not be selected...\n";
  std::cout << "                           The selection string passed to this method is: \"" << selection << "\"\n";
  return false;
}



bool MGSelector::SelectBsToJPsiPhiEvent(const std::string selection, std::vector<std::pair<int,int> >& selectedObjects)
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
    selectedObjects.push_back(selObject);
    selObject = std::make_pair(PDEnumString::recPV, iPV);
     selectedObjects.push_back(selObject);
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
    selectedObjects.push_back(selObject);
    selObject = std::make_pair(PDEnumString::recPV, iPV);
    selectedObjects.push_back(selObject);
    return true;
  }
  std::cout << "MGSelector::SelectBsToJPsiPhiEvent(): W A R N I N G ! Reached default return statement. Event will not be selected...\n";
  std::cout << "                                      The selection string passed to this method is: \"" << selection << "\"\n";
  return false;
}



bool MGSelector::SelectBuToJPsiKEvent(const std::string selection, std::vector<std::pair<int, int> >& selectedObjects)
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
