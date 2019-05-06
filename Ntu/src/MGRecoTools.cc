#define UTIL_USE FULL

#include <vector>
#include <iostream>

#include "TVector3.h"

#include "PDAnalysis/Ntu/interface/constants.h"
#include "PDAnalysis/Ntu/interface/MGRecoTools.h"



// #define mgDEBUG



short MGRecoTools::FindJetNearToMuon(const unsigned short iMuon, const double deltaRThreshold)
{
  TVector3 pMuon, pJet;
  pMuon.SetPtEtaPhi(muoPt->at(iMuon), muoEta->at(iMuon), muoPhi->at(iMuon));
  short nearestJetIndex = -1;
  double maxDeltaR = 99999;
  double deltaR = 99999;
  for(int iJet = 0; iJet < nJets; iJet++)
  {
    pJet.SetPtEtaPhi(jetPt->at(iJet), jetEta->at(iJet), jetPhi->at(iJet));
    deltaR = pMuon.DeltaR(pJet);
    if(deltaR < maxDeltaR && deltaR < deltaRThreshold)
    {
      nearestJetIndex = iJet;
      maxDeltaR = deltaR;
    }
  }
  return nearestJetIndex;
}



const int MGRecoTools::GetClosestRecoElectron(const double pt, const double eta, const double phi, double dRMax, double dPtMax)
{
//   std::cout << "MGRecoTools::GetClosestRecoElectron(...): pt = " << pt << ", eta = " << eta << ", phi = " << phi << std::endl;
//   std::cout << "                                          Matching parameters: dRMax = " << dRMax << ", dPtMax = " << dPtMax << std::endl;
  int matched = -1;
  double dRMatch = dRMax;
  double dPtMatch = dPtMax;
  
  std::cout << "nElectrons = " << nElectrons << std::endl;
  
  for(int iEle = 0; iEle < nElectrons; iEle++)
  {
    float dR = deltaR(eta, phi, eleEta->at(iEle), elePhi->at(iEle));
    float dPt = fabs(elePt->at(iEle) - pt)/elePt->at(iEle);
    
//     std::cout << "       dR = " << dR << ", dPt = " << dPt << std::endl;
    
    if( dR > dRMatch ) continue;
//     std::cout << "            dR  inside matching window!\n";
    if( dPt > dPtMatch) continue;
//     std::cout << "            dPt inside matching window!\n";
    
    matched = iEle;
//     std::cout << "       matched = " << matched << std::endl;
    dRMatch = dR;
//     std::cout << "       dRMatch = " << dRMatch << std::endl;
  } 
  
  return matched;  
}



const TLorentzVector MGRecoTools::GetTLorentzVectorFromJPsiX(const int iSvt)
{
    int iJPsi = (subVtxFromSV(iSvt)).at(0);
    std::vector<int> tkJpsi = tracksFromSV(iJPsi) ;
    std::vector<int> tkSsB = tracksFromSV(iSvt);
    TLorentzVector t(0,0,0,0);

    for( uint i=0; i<tkSsB.size(); ++i ){
        int j = tkSsB.at(i);
        float m = constants::kaonMass;
        if( (j==tkJpsi.at(0)) || (j==tkJpsi.at(1)) ) m = constants::muonMass;
        TLorentzVector a;
        a.SetPtEtaPhiM( trkPt->at(j), trkEta->at(j), trkPhi->at(j), m );
        t += a;
    }

    return t;
}



float MGRecoTools::GetCt2D(const TLorentzVector& t, const int iSV)
{
    TVector3 vSVT(svtX->at(iSV), svtY->at(iSV), 0.);
    TVector3 vPV(bsX, bsY, 0.);

    TVector3 vPointing = vSVT - vPV;
    TVector3 vBs = t.Vect();

    return constants::BsMass/t.Pt() * (vPointing * vBs.Unit());
}




const double MGRecoTools::dZ(const int iEle, const int iVtx)
{
//   std::cout << "Inside MGRecoTools::dZ(...)\n";
//   std::cout << "       iEle = " << iEle << ", iVtx = " << iVtx << std::endl;
  //  number px = trkVtxPx->at( iTk );
  //  number py = trkVtxPy->at( iTk );
  //  number pz = trkVtxPz->at( iTk );
  number px;
  number py;
  number pz;
 if(eleGsfPx->size())
  {
    px = eleGsfPx->at(iEle);
    py = eleGsfPy->at(iEle);
    pz = eleGsfPz->at(iEle);
  }
  else if(eleGsfPt->size())
  {
    TVector3 pEle;
    pEle.SetPtEtaPhi(eleGsfPt->at(iEle), eleGsfEta->at(iEle), eleGsfPhi->at(iEle));
    px = pEle.Px();
    py = pEle.Py();
    pz = pEle.Pz();
  }
  else  if(eleGsfPxAtVtx->size())
  {
    px = eleGsfPxAtVtx->at(iEle);
    py = eleGsfPyAtVtx->at(iEle);
    pz = eleGsfPzAtVtx->at(iEle);
  }
  else if(eleGsfPtAtVtx->size())
  {
    TVector3 pEle;
    pEle.SetPtEtaPhi(eleGsfPtAtVtx->at(iEle), eleGsfEtaAtVtx->at(iEle), eleGsfPhiAtVtx->at(iEle));
    px = pEle.Px();
    py = pEle.Py();
    pz = pEle.Pz();
  }
  else if(elePx->size())
  {
    px = elePx->at(iEle);
    py = elePy->at(iEle);
    pz = elePz->at(iEle);
  }
  else
  {
    TVector3 pEle;
    pEle.SetPtEtaPhi(elePt->at(iEle), eleEta->at(iEle), elePhi->at(iEle));
    px = pEle.Px();
    py = pEle.Py();
    pz = pEle.Pz();
  }
  
  number pq = modSqua(px, py, 0);
  return (((px * (pvtX->at(iVtx) - bsX))) + (py * (pvtY->at(iVtx) - bsY)) * pz/pq) +
            eleGsfDz->at(iEle) + bsZ - pvtZ->at(iVtx);
}



// double RecoTools::GetTrackValidFraction(const unsigned short iTrack)
// {
//   int valid = GetNumberOfValidTrackerHits(iTrack);
//   int lost  = GetNumberOfLostTrackerHits(iTrack);
//   int lostIn = GetNumberOfLostTrackerInnerHits(iTrack);
//   int lostOut = GetNumberOfLostTrackerOuterHits(iTrack);
//   if ((valid+lost+lostIn+lostOut)==0)
//     return -1;
//   return valid/(1.0*(valid+lost+lostIn+lostOut));
// }



// int RecoTools::GetNumberOfValidTrackerHits(const unsigned short iTrack)
// {
//   int count = 0;
//   
//   unsigned int size = Track_hitPattern->at(iTrack).size();
//   for (unsigned int i = 0; i < size; i++)
//   {
//     unsigned int pattern = trkHitPattern->at(iTrack).at(i);
//     if (pattern == 0) break;
//     if (ValidHitFilter(pattern) && TrackerHitFilter(pattern)) ++count;
//   }
//   return count;
// }



// int RecoTools::GetNumberOfValidTrackerInnerHits(const unsigned short iTrack)
// {
//   int count = 0;
//   unsigned int size = Track_expInnHitPattern->at(iTrack).size();
//   for (unsigned int i = 0; i < size; i++)
//   {
//     unsigned int pattern = Track_expInnHitPattern->at(iTrack).at(i);
//     if (pattern == 0) break;
//     if (ValidHitFilter(pattern) && TrackerHitFilter(pattern)) ++count;
//   }
//   return count;
// }



// int RecoTools::GetNumberOfValidTrackerOuterHits(const unsigned short iTrack)
// {
//   int count = 0;
//   unsigned int size = Track_expOutHitPattern->at(iTrack).size();
//   for (unsigned int i = 0; i < size; i++)
//   {
//     unsigned int pattern = Track_expOutHitPattern->at(iTrack).at(i);
//     if (pattern == 0) break;
//     if (ValidHitFilter(pattern) && TrackerHitFilter(pattern)) ++count;
//   }
//   return count;
// }



// int RecoTools::GetNumberOfLostTrackerHits(const unsigned short iTrack)
// {
//   int count = 0;
//   unsigned int size = Track_hitPattern->at(iTrack).size();
//   for (unsigned int i = 0; i < size; i++)
//   {
//     unsigned int pattern = Track_hitPattern->at(iTrack).at(i);
//     if (pattern == 0) break;
//     if (Type1HitFilter(pattern) && TrackerHitFilter(pattern)) ++count;
//   }
//   return count;
// }



// int RecoTools::GetNumberOfLostTrackerInnerHits(const unsigned short iTrack)
// {
//   int count = 0;
//   unsigned int size = Track_expInnHitPattern->at(iTrack).size();
//   for (unsigned int i = 0; i < size; i++)
//   {
//     unsigned int pattern = Track_expInnHitPattern->at(iTrack).at(i);
//     if (pattern == 0) break;
//     if (Type1HitFilter(pattern) && TrackerHitFilter(pattern)) ++count;
//   }
//   return count;
// }



// int RecoTools::GetNumberOfLostTrackerOuterHits(const unsigned short iTrack)
// {
//   int count = 0;
//   unsigned int size = Track_expOutHitPattern->at(iTrack).size();
//   for (unsigned int i = 0; i < size; i++)
//   {
//     unsigned int pattern = Track_expOutHitPattern->at(iTrack).at(i);
//     if (pattern == 0) break;
//     if (Type1HitFilter(pattern) && TrackerHitFilter(pattern)) ++count;
//   }
//   return count;
// }



// int RecoTools::CountTrackTypedHits(const std::vector<unsigned short>& hitPattern, const hitFilterType typeFilter)
// {
//   int count = 0;
//   unsigned int size = hitPattern.size();
//   for (unsigned int i = 0; i < size; i++)
//   {
//     unsigned int pattern = hitPattern.at(i);
//     if (pattern == 0) break;
//     if (typeFilter(pattern) && TrackerHitFilter(pattern)) ++count;
//   }
//   return count;
// }



// unsigned int RecoTools::GetTrackHitType(const unsigned int hitPattern) const
// {
//   if(hitPattern == 0)
//   {
//     return 999999;
//   } 
//   return ((hitPattern>>HitTypeOffset) & HitTypeMask);
// }



// inline bool RecoTools::ValidHitFilter(const unsigned int hitPattern)
// {
//   if (GetTrackHitType(hitPattern) == 0) return true; 
//   return false;
// }



// inline bool RecoTools::Type1HitFilter(const unsigned int hitPattern)
// {
//   if (GetTrackHitType(hitPattern) == 1) return true; 
//   return false;
// }



// inline bool RecoTools::Type2HitFilter(const unsigned int hitPattern)
// {
//   if (GetTrackHitType(hitPattern) == 2) return true; 
//   return false;
// }



// inline bool RecoTools::Type3HitFilter(const unsigned int hitPattern)
// {
//   if (GetTrackHitType(hitPattern) == 3) return true; 
//   return false;
// }



// inline bool RecoTools::TrackerHitFilter(const unsigned int hitPattern)
// {
//   if(hitPattern == 0)
//   {
//     return false;
//   }  
//   if(((hitPattern>>HitSubDetectorOffset) & HitSubDetectorMask) == 1)
//   {
//     return true;
//   }
//   return false;
// }
