#ifndef PDSecondNtupleData_h
#define PDSecondNtupleData_h

#include "NtuTool/Common/interface/TreeWrapper.h"

class PDSecondNtupleData: public virtual TreeWrapper {

 public:

  PDSecondNtupleData() {
  }
  virtual ~PDSecondNtupleData() {
  }

  void initTree() {
    treeName = "PDsecondTree";
    setBranch( "nSum", &nSum , "nSum/I", &b_nSum );
    setBranch( "nPro", &nPro , "nPro/I", &b_nPro );
  }

  int nSum;
  int nPro;
  TBranch* b_nSum;
  TBranch* b_nPro;

 private:

  PDSecondNtupleData           ( const PDSecondNtupleData& a );
  PDSecondNtupleData& operator=( const PDSecondNtupleData& a );

};

#endif

