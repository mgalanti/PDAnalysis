#define UTIL_USE FULL

#include "PDAnalysis/Ntu/interface/MGBaseTools.h"



void MGBaseTools::beginJob()
{
  getUserParameter( "verbose", verbose );
}

