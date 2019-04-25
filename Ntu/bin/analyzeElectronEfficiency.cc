#define UTIL_USE FULL
#include "ElectronEfficiencyAnalyzer.h"
#include "NtuAnalysis/Read/interface/NtuFlexibleAnalyzer.h"
#include "NtuTool/Read/bin/treeAnalyze.cc"
static NtuFlexibleAnalyzer<ElectronEfficiencyAnalyzer> nfa( "electronEfficiencyAnalysis", "electronEfficiencyAnalyzer" );
 
