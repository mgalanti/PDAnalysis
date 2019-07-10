//#define UTIL_USE FULL
//#include "PDAnalyzer.h"
//#include "NtuAnalysis/Read/interface/NtuFlexibleAnalyzer.h"
//#include "NtuTool/Read/bin/treeAnalyze.cc"
//static NtuFlexibleAnalyzer<PDSecondAnalyzer> nfa( "pdAnalysis", "pdAnalyzer" );
#include "PDAnalysis/EleMVACalibrationProducer/interface/EleMVACalibrationProducer.h"
#include "PDAnalysis/EleMVACalibrationProducer/src/EleMVACalibrationProducer.cc"
// #include "PDAnalysis/Ntu/bin/treeAnalyzeMod.cc"

#include "NtuAnalysis/Read/interface/NtuFlexibleAnalyzer.h"
#include "NtuTool/Read/bin/treeAnalyze.cc"
// static NtuFlexibleAnalyzer<EleMVASecondNtupleReader> snr("eleMVASecondNtupleReading", "eleMVASecondNtupleReader");
static EleMVACalibrationProducer ecp;
