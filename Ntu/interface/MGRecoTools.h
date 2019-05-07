#ifndef MGRecoTools_h
#define MGRecoTools_h

#include "TLorentzVector.h"
#include "TVectorF.h"

#include "PDAnalysis/Ntu/interface/constants.h"
#include "PDAnalysis/Ntu/interface/PDAnalyzerUtil.h"



class MGRecoTools : public virtual PDAnalyzerUtil
{
  public:
//     typedef bool hitFilterType(const unsigned int);
    
    virtual ~MGRecoTools() {};
    
    inline short FindJetNearToMuon(const unsigned short iMuon, const double deltaRThreshold);
    
    
    // GetCt* methods adapted from https://github.com/abragagn/BPHPD-AlbertoUtilities/blob/master/PDAnalysis/Ntu/bin/AlbertoUtil.cc
    float GetCt2D(const TLorentzVector& p, const int iSV, const float mass = constants::BsMass);
    float GetCt2DPV(const TLorentzVector& p, const int iSV, const int iPV, const float mass = constants::BsMass);
    float GetCt2DPVErr(const TLorentzVector& p, const int iSV, const int iPV, const float mass = constants::BsMass);
    float GetCt3DPV(const TLorentzVector& p, const int iSV, const int iPV, const float mass = constants::BsMass);
    float GetCt3DPVErr(const TLorentzVector& p, const int iSV, const int iPV, const float mass = constants::BsMass);
    
    const int GetClosestRecoElectron(const double pt, const double eta, const double phi, double dRMax = 0.12, double dPtMax = 0.3);
    
    // Taken from https://github.com/abragagn/BPHPD-AlbertoUtilities/blob/master/PDAnalysis/Ntu/bin/AlbertoUtil.cc
    const TLorentzVector GetTLorentzVectorFromJPsiX(const int iSvt);
    
    const double dZ(const int iEle, const int iVtx);
   
    // MG: none of the methods below works with the hitPattern implementation in PD ntuples.
//     inline double GetTrackValidFraction(const unsigned short iTrack);
//     inline int GetNumberOfValidTrackerHits(const unsigned short iTrack);
//     inline int GetNumberOfValidTrackerInnerHits(const unsigned short iTrack);
//     inline int GetNumberOfValidTrackerOuterHits(const unsigned short iTrack);
//     inline int GetNumberOfLostTrackerHits(const unsigned short iTrack);
//     inline int GetNumberOfLostTrackerInnerHits(const unsigned short iTrack);
//     inline int GetNumberOfLostTrackerOuterHits(const unsigned short iTrack);
//     inline int CountTrackTypedHits(const std::vector<unsigned short>& hitPattern, const hitFilterType typeFilter);
//     inline unsigned int GetTrackHitType(const unsigned int hitPattern) const;
//     inline bool ValidHitFilter(const unsigned int hitPattern);
//     inline bool Type1HitFilter(const unsigned int hitPattern);
//     inline bool Type2HitFilter(const unsigned int hitPattern);
//     inline bool Type3HitFilter(const unsigned int hitPattern);
//     inline bool TrackerHitFilter(const unsigned int hitPattern);
  private:
//     static const unsigned short HitSize = 11;
//     static const unsigned short HitTypeMask = 0x3;
//     static const unsigned short HitTypeOffset = 0;
//     static const unsigned short HitLayerMask = 0xF;
//     static const unsigned short HitLayerOffset = 3;
//     static const unsigned short HitPatternSize = 25;
//     static const unsigned short HitSideMask = 0x1;
//     static const unsigned short HitSideOffset = 2;
//     static const unsigned short HitSubDetectorMask = 0x1;
//     static const unsigned short HitSubDetectorOffset = 10;
//     static const unsigned short HitSubstrMask = 0x7;
//     static const unsigned short HitSubstrOffset = 7;
};

#endif // MGRecoTools_h
