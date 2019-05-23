#ifndef MGGenTools_h
#define MGGenTools_h

#define ARRAY_SIZE(array) (sizeof((array))/sizeof((array[0])))

#include "TLorentzVector.h"

#include "PDAnalysis/Ntu/interface/MGBaseTools.h"
#include "PDAnalysis/Ntu/interface/PDGenHandler.h"
#include "MGTools/AnalysisTools/interface/MGContainerTools.h"



class MGGenTools: 
    public virtual MGBaseTools,
    public virtual PDGenHandler, 
    public virtual MGContainerTools
{
  public:
    MGGenTools();
    virtual ~MGGenTools() {};
    
    // Matches to the closest long-lived gen particle (see list below)
    const int GetClosestGen(const double pt, const double eta, const double phi, double dRMax = 0.12, double dPtMax = 0.3);
    // Matches to the closest gen particle (no long-lived requirement)
    const int GetClosestGenNoLL(const double pt, const double eta, const double phi, double dRMax = 0.12, double dPtMax = 0.3);
    // To be used for composite particles (belonging to the svt* collection)
    const int GetClosestGenNoLL(const int iSvt, double dRMax = 0.12, double dPtMax = 0.3);
    // Matches to the closest gen particle among those in the vector provided
    const int GetClosestGenInList(const double pt, const double eta, const double phi, const std::vector<int>& listGenP, double dRMax = 0.12, double dPtMax = 0.3);
    // Gets the charge correlation between generated lepton (e or mu) and generated Bs/Bu
    // Returns 1 ("right" correlation) or -1 ("wrong" correlation)
    // If some of the inputs are wrong, also returns 0 (and prints out error)
    const int GetGenLepBsChargeCorrelation(const int iGenLep, int iGenBs);
    const int GetGenLepBuChargeCorrelation(const int iGenLep, int iGenBs);
    
    // Gets the charge correlation between a charge value (can only be +/-1) and generated Bs/Bu
    // Returns 1 ("right" correlation) or -1 ("wrong" correlation)
    // If some of the inputs are wrong, returns 0 (and prints out error)
    const int GetBsChargeCorrelation(const int charge, int iGenBs);
    const int GetBuChargeCorrelation(const int charge, int iGenBu);
    
    const int GetMixStatus(const uint iGen);
    
    const std::vector<int> GetAllGenElectrons();
    const std::vector<int> GetAllGenElectronsFromB();
    const std::vector<int> GetAllLongLivedBHadrons();
    
    const void PrintMotherChain(const int iGen); // From Alberto Bragagnolo
//     const void PrintDaughterTree(const int iGen); // From Alberto Bragagnolo
//     const void PrintDaughterTreePt(const int iGen); // From Alberto Bragagnolo
    const void RecursivePrintMothers(const unsigned short iGen, int recursionOrder = 0, std::string prepend = "");
    const void RecursivePrintDaughters(const unsigned short iGen, int recursionOrder = 0, std::string prepend = "");
    const int RecursiveLookForMotherIds(const unsigned short iGen, const std::vector<int> vIdList, const bool withSign = true, const int recursionOrder = 0);
    const int RecursiveLookForDaughterIds(const unsigned short iGen, const std::vector<int> vIdList, const bool withSign = true, const int recursionOrder = 0);
  
    const bool IsLongLived(const uint iGen);
    
    const std::vector<int> listLongLived =         {11,13,211,321,2212}; //e, mu, pi, K, p  
    const std::vector<int> listLongLivedBHadrons = {511, 521, 531, 541, 5122}; // The ground-state B hadrons
    
    const std::vector<int> listBMesons =      {511, 521, 10511, 10521, 513, 523, 10513, 10523, 20513, 
                                               20523, 515, 525, 531, 10531, 533, 10533, 20533, 535, 541, 
                                               10541, 543, 10543, 20543, 545};
    const std::vector<int> listBBaryons =     {5122, 5112, 5212, 5222, 5114, 5214, 5224, 5132, 5232, 
                                               5312, 5322, 5314, 5324, 5332, 5334, 5142, 5242, 5412, 
                                               5422, 5414, 5424, 5342, 5432, 5434, 5442, 5444, 5512, 
                                               5522, 5514, 5524, 5532, 5534, 5542, 5544, 5554};      
    const std::vector<int> listBottomonium =  {551, 10551, 100551, 110551, 200551, 210551, 553, 10553, 
                                               20553, 30553, 100553, 110553, 120553, 130553, 200553, 210553, 
                                               220553, 300553, 9000553, 9010553, 555, 10555, 20555, 100555, 
                                               110555, 120555, 200555, 557, 100557};
    const std::vector<int> listBMesonsAndBaryons = Join(listBMesons, listBBaryons);
    const std::vector<int> listAllBHadrons =       JoinMany<int>({listBMesons, listBBaryons, listBottomonium});
    
    const std::vector<int> listCMesons =    {411, 421, 10411, 10421, 413, 423, 10413, 10423, 20413, 20423, 
                                             415, 425, 431, 10431, 433, 10433, 20433, 435};
    const std::vector<int> listCBaryons =   {4122, 4222, 4212, 4112, 4224, 4214, 4114, 4232, 4132, 4322, 4312, 
                                             4324, 4314, 4332, 4334, 4412, 4422, 4414, 4424, 4432, 4434, 4444}; 
    const std::vector<int> listCharmonium = {441, 10441, 100441, 443, 10443, 20443, 100443, 30443, 9000443,
                                             9010443, 9020443, 445, 100445};
    const std::vector<int> listCMesonsAndBaryons = Join(listCMesons, listCBaryons);
    const std::vector<int> listAllCHadrons =       JoinMany<int>({listCMesons, listCBaryons, listCharmonium});

  private:
    // Nothing so far...
};



#endif // MGGenTools_h
