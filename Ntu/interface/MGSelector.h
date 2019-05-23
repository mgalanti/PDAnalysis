#ifndef MGSelector_H
#define MGSelector_H

#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "MGTools/AnalysisTools/interface/MGStringTools.h"
#include "PDAnalysis/Ntu/interface/MGRecoTools.h"

class MGSelector: 
    public virtual MGRecoTools, 
    public virtual MGStringTools
{  
  public:
    
    MGSelector();
    virtual ~MGSelector();
    
    void SetSelectionString(const std::string str);
    
    bool SelectEvent();  // Default using saved selection string
    bool SelectEvent(const std::string selection, std::vector<std::pair<int, int> >& selectedObjectsRef);
    bool SelectBsToJPsiPhiEvent(const std::string selection, std::vector<std::pair<int,int> >& selectedObjectsRef);
    bool SelectBuToJPsiKEvent(const std::string selection, std::vector<std::pair<int,int> >& selectedObjectsRef);
    
    bool SelectHlt(const std::string selection);
    bool SelectElectron(const int iEle, const std::string selection);
    // Select best electron candidate for OS flavor tagging
    // Adapted from GetOSMuon in:
    // https://github.com/abragagn/BPHPD-AlbertoUtilities/blob/master/PDAnalysis/Ntu/bin/OSMuonMvaTag.cc
    int SelectOSElectron(const std::string selection, const int iPV, const int iB);
    bool SelectJPsi(const int iJPsi, const std::string selection);
    bool SelectPhi(const int iPhi, const std::string selection);
    
    int SelectBestCandidate(const std::string candidateType, const std::string selection, int& iPV);
    int SelectBestBsToJPsiPhi(const std::string selection, int& iPV);
    int SelectBestBuToJPsiK(const std::string selection, int& iPV);
    int SelectBestPV(const int iSvt, const TLorentzVector& pSvt, const std::string selection);
    
    std::vector<std::pair<int, int> > selectedObjects;
    std::string selectionString;
    std::vector<std::string> selectionSubStrings;

  private:
    
};

#endif
