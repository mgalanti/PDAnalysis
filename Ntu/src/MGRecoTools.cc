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



const float MGRecoTools::GetFixedGridRhoFastJetAll()
{
  for(int iUserInfo = 0; iUserInfo < nUserInfo; iUserInfo++)
  {
    if(useObjType->at(iUserInfo) == PDEnumString::eventInfo && useInfoType->at(iUserInfo) == PDEnumString::fixedGridRhoFastjetAll)
    {
      return useInfoValue->at(iUserInfo);
    }
  }
  return 0;
}



const float MGRecoTools::GetFixedGridRhoFastJetAllCalo()
{
  for(int iUserInfo = 0; iUserInfo < nUserInfo; iUserInfo++)
  {
    if(useObjType->at(iUserInfo) == PDEnumString::eventInfo && useInfoType->at(iUserInfo) == PDEnumString::fixedGridRhoFastjetAllCalo)
    {
      return useInfoValue->at(iUserInfo);
    }
  }
  return 0;
}
    


const int MGRecoTools::GetClosestRecoElectron(const double pt, const double eta, const double phi, double dRMax, double dPtMax)
{
//   std::cout << "MGRecoTools::GetClosestRecoElectron(...): pt = " << pt << ", eta = " << eta << ", phi = " << phi << std::endl;
//   std::cout << "                                          Matching parameters: dRMax = " << dRMax << ", dPtMax = " << dPtMax << std::endl;
  int matched = -1;
  double dRMatch = dRMax;
  double dPtMatch = dPtMax;
  
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



float MGRecoTools::GetCt2D(const TLorentzVector& p, const int iSV, const float mass)
{
  TVector3 vSVT(svtX->at(iSV), svtY->at(iSV), 0.);
  TVector3 vPV(bsX, bsY, 0.);
  
  TVector3 vPointing = vSVT - vPV;
  TVector3 vB = p.Vect();
  
  return mass/p.Pt() * (vPointing * vB.Unit());
}



float MGRecoTools::GetCt2DPV(const TLorentzVector& p, const int iSV, const int iPV, const float mass)
{
  TVector3 vSVT(svtX->at(iSV), svtY->at(iSV), 0.);
  TVector3 vPV(pvtX->at(iPV), pvtY->at(iPV), 0.);
  
  TVector3 vPointing = vSVT - vPV;
  TVector3 vBs = p.Vect();
  
  return mass/p.Pt() * (vPointing * vBs.Unit());
}



float MGRecoTools::GetCt2DPVErr(const TLorentzVector& p, const int iSV, const int iPV, const float mass)
{
  TVector3 vSVT(svtX->at(iSV), svtY->at(iSV), 0.);
  TVector3 vPV(pvtX->at(iPV), pvtY->at(iPV), 0.);
  
  TVector3 vPointing = vSVT - vPV;
  TVector3 vB = p.Vect();
  
  TMatrixF covSV(3,3);
  float covSVArray[]={svtSxx->at(iSV),svtSxy->at(iSV),svtSxz->at(iSV),
                      svtSxy->at(iSV),svtSyy->at(iSV),svtSyz->at(iSV), 
                      svtSxz->at(iSV),svtSyz->at(iSV),svtSzz->at(iSV)};
  covSV.SetMatrixArray(covSVArray);
  
  TMatrixF covPV(3,3);
  float covPVArray[]={pvtSxx->at(iPV),pvtSxy->at(iPV),pvtSxz->at(iPV),
                      pvtSxy->at(iPV),pvtSyy->at(iPV),pvtSyz->at(iPV), 
                      pvtSxz->at(iPV),pvtSyz->at(iPV),pvtSzz->at(iPV)};
  covPV.SetMatrixArray(covPVArray);
  
  TMatrixF covTot= covSV+covPV;
  
  float distArray2D[]={float(vPointing.X()),float(vPointing.Y()),0.};
  TVectorF diff2D(3,distArray2D);
  
  if (diff2D.Norm2Sqr()==0) return -1.; //if the secondary vertex is exactly the same as PV 
  
  return mass / p.Pt() * sqrt(covTot.Similarity(diff2D)) / sqrt(diff2D.Norm2Sqr()); 
}



float MGRecoTools::GetCt3DPV(const TLorentzVector& p, const int iSV, const int iPV, const float mass)
{
  TVector3 vSVT(svtX->at(iSV), svtY->at(iSV), svtZ->at(iSV));
  TVector3 vPV(pvtX->at(iPV), pvtY->at(iPV), pvtZ->at(iPV));
  
  TVector3 vPointing = vSVT - vPV;
  TVector3 vB = p.Vect();
  
  return mass/p.P() * vPointing.Dot(vB)/vB.Mag();
}



float MGRecoTools::GetCt3DPVErr(const TLorentzVector& p, const int iSV, const int iPV, const float mass)
{
  TVector3 vSVT(svtX->at(iSV), svtY->at(iSV), svtZ->at(iSV));
  TVector3 vPV(pvtX->at(iPV), pvtY->at(iPV), pvtZ->at(iPV));
  
  TVector3 vPointing = vSVT - vPV;
  TVector3 vB = p.Vect();
  
  TMatrixF covSV(3,3);
  float covSVArray[]={svtSxx->at(iSV),svtSxy->at(iSV),svtSxz->at(iSV),
                      svtSxy->at(iSV),svtSyy->at(iSV),svtSyz->at(iSV), 
                      svtSxz->at(iSV),svtSyz->at(iSV),svtSzz->at(iSV)};
  covSV.SetMatrixArray(covSVArray);
  
  TMatrixF covPV(3,3);
  float covPVArray[]={pvtSxx->at(iPV),pvtSxy->at(iPV),pvtSxz->at(iPV),
                      pvtSxy->at(iPV),pvtSyy->at(iPV),pvtSyz->at(iPV), 
                      pvtSxz->at(iPV),pvtSyz->at(iPV),pvtSzz->at(iPV)};
  covPV.SetMatrixArray(covPVArray);
  
  TMatrixF covTot= covSV+covPV;
  
  float distArray[]={float(vPointing.X()),float(vPointing.Y()),float(vPointing.Z())};
  TVectorF diff(3,distArray);

  if (diff.Norm2Sqr()==0) return -1.; //if the secondary vertex is exactly the same as PV 
  
  return mass/p.P() * sqrt(covTot.Similarity(diff)) / sqrt(diff.Norm2Sqr()); 
}



const double MGRecoTools::dXEle(const int iEle)
{
  float px;
  float py;
  if(eleGsfPx->size())
  {
    px = eleGsfPx->at(iEle);
    py = eleGsfPy->at(iEle);
  }
  else
  {
    TVector3 pEle;
    pEle.SetPtEtaPhi(eleGsfPt->at(iEle), eleGsfEta->at(iEle), eleGsfPhi->at(iEle));
    px = pEle.Px();
    py = pEle.Py();
  }
  float pm = modCart(px, py, 0);
  return -eleGsfDxy->at(iEle) * py / pm;
}



const double MGRecoTools::dYEle(const int iEle)
{
  float px;
  float py;
  if(eleGsfPx->size())
  {
    px = eleGsfPx->at(iEle);
    py = eleGsfPy->at(iEle);
  }
  else
  {
    TVector3 pEle;
    pEle.SetPtEtaPhi(eleGsfPt->at(iEle), eleGsfEta->at(iEle), eleGsfPhi->at(iEle));
    px = pEle.Px();
    py = pEle.Py();
  }
  float pm = modCart(px, py, 0);
  return -eleGsfDxy->at(iEle) * px / pm;
}



const double MGRecoTools::dXEle(const int iEle, const int iVtx)
{
  float px;
  float py;
  if(eleGsfPx->size())
  {
    px = eleGsfPx->at(iEle);
    py = eleGsfPy->at(iEle);
  }
  else
  {
    TVector3 pEle;
    pEle.SetPtEtaPhi(eleGsfPt->at(iEle), eleGsfEta->at(iEle), eleGsfPhi->at(iEle));
    px = pEle.Px();
    py = pEle.Py();
  }
  float pq = modSqua(px, py, 0);
  float pm = sqrt( pq );
  return ( ( ( ( bsX - pvtX->at(iVtx) ) * py ) + ( ( pvtY->at(iVtx) - bsY ) * px ) ) * py / pq ) -
         ( eleGsfDxy->at( iEle ) * py / pm );
}



const double MGRecoTools::dYEle(const int iEle, const int iVtx)
{
  float px;
  float py;
  if(eleGsfPx->size())
  {
    px = eleGsfPx->at(iEle);
    py = eleGsfPy->at(iEle);
  }
  else
  {
    TVector3 pEle;
    pEle.SetPtEtaPhi(eleGsfPt->at(iEle), eleGsfEta->at(iEle), eleGsfPhi->at(iEle));
    px = pEle.Px();
    py = pEle.Py();
  }
  float pq = modSqua( px, py, 0 );
  float pm = sqrt( pq );
  return ( ( ( ( bsY - pvtY->at(iVtx) ) * px ) + ( ( pvtX->at(iVtx) - bsX ) * py ) ) * px / pq ) +
         ( eleGsfDxy->at( iEle ) * px / pm );
}



const double MGRecoTools::dXYEle(const int iEle, const int iVtx)
{
  float vx = pvtX->at(iVtx);
  float vy = pvtY->at(iVtx);
  float px;
  float py;
  if(eleGsfPx->size())
  {
    px = eleGsfPx->at(iEle);
    py = eleGsfPy->at(iEle);
  }
  else
  {
    TVector3 pEle;
    pEle.SetPtEtaPhi(eleGsfPt->at(iEle), eleGsfEta->at(iEle), eleGsfPhi->at(iEle));
    px = pEle.Px();
    py = pEle.Py();
  }
  float pq = modSqua( px, py, 0 );
  float pm = sqrt( pq );
  float dx = ( vx - bsX ) * py / pm;
  float dy = ( vy - bsY ) * px / pm;
  float d1 = eleGsfDxy->at( iEle );
  number xy = sqrt( ( dx * dx ) +     ( dy * dy ) - ( 2 * dx * dy )
                  + ( d1 * d1 ) + ( 2 * d1 * dx ) - ( 2 * d1 * dy ) );
  return ( ( px * dYEle( iEle, iVtx ) ) - ( py * dXEle( iEle, iVtx ) ) > 0 ?
           xy : -xy );
}



const double MGRecoTools::dZEle(const int iEle, const int iVtx)
{
//   std::cout << "Inside MGRecoTools::dZEle(...)\n";
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



const double MGRecoTools::dXTrk(const int iTrk, const int iVtx)
{
    return PDAnalyzerUtil::dX(iTrk, pvtX->at(iVtx), pvtY->at(iVtx));
}



const double MGRecoTools::dYTrk(const int iTrk, const int iVtx)
{
    return PDAnalyzerUtil::dY(iTrk, pvtX->at(iVtx), pvtY->at(iVtx));
}



const double MGRecoTools::dXYTrk(const int iTrk, const int iVtx)
{
    return PDAnalyzerUtil::dXY(iTrk, pvtX->at(iVtx), pvtY->at(iVtx));
}



const double MGRecoTools::dZTrk(const int iTrk, const int iVtx)
{
    return PDAnalyzerUtil::dZ(iTrk, pvtX->at(iVtx), pvtY->at(iVtx), pvtZ->at(iVtx));
}



const int MGRecoTools::dSignEle(const int iEle, const float px, const float py)
{
  return (((dXEle(iEle) * px) + 
           (dYEle(iEle) * py)) > 0 ? 1 : -1);
}



const float MGRecoTools::GetEleEffArea(const float eta, const std::string dataName)
{
  float absEta = fabs(eta);
  if(dataName.compare("effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X") == 0)
  {
    if(absEta < 1.)
      return 0.1440;
    if(absEta < 1.479)
      return 0.1562;
    if(absEta < 2.)
      return 0.1032;
    if(absEta < 2.2)
      return 0.0859;
    if(absEta < 2.3)
      return 0.1116;
    if(absEta < 2.4)
      return 0.1321;
    if(absEta < 2.5)
      return 0.1654;
    else
      return 0;
  }
  else
    return 0;
}



const double MGRecoTools::GetEleRelPFIsoScaled(const int iEle) 
{
  double absEta = eleAbsEta->at(iEle);
  
  // Compute the combined isolation with effective area correction
  const float chad = eleSumCHpt->at(iEle);
  const float nhad = eleSumNHet->at(iEle);
  const float pho  = eleSumPHet->at(iEle);
  const float  eA  = GetEleEffArea(absEta);
  const float rho  = GetFixedGridRhoFastJetAll(); // std::max likes float arguments
  const float iso  = chad + std::max(0.0f, nhad + pho - rho*eA);
  return iso/elePt->at(iEle);
}
