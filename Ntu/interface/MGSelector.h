#ifndef MGSelector_H
#define MGSelector_H

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "PDAnalysis/Ntu/interface/MGRecoTools.h"
#include "PDAnalysis/Ntu/interface/PDAnalyzerUtil.h"

class MGSelector: public virtual PDAnalyzerUtil, public virtual MGRecoTools
{  
  public:
    
    MGSelector();
    virtual ~MGSelector();
    
    bool SelectBsToJPsiPhiEvent(const std::string selection, std::vector<std::pair<int,int> >& selectedObjects);
    
    bool SelectHlt(const std::string selection);
    bool SelectElectron(const int iEle, const std::string selection);
    bool SelectJPsi(const int iJPsi, const std::string selection);
    bool SelectPhi(const int iPhi, const std::string selection);
    
    int SelectBestBsToJPsiPhi(const std::string selection, int& iPV);
    int SelectBestPV(const int iSvt, const TLorentzVector& pSvt, const std::string selection);
    
  private:
    // Nothing at the moment
};

#endif
