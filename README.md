# PDAnalysis
Analysis code based on the PDAnalysis framework
The repository only contains differences with respect to the standard PDAnalysis package provided by Paolo.
It does not run standalone and needs to be installed on top of full PDAnalysis.
Other dependencies:
- MGTools
- NtuTool
- NtuAnalysis

Instructions to install on top of PDAnalysis 2017V12 on CMSSW_10_2_6:

    Download PDAnalysis according to Paolo installer

In the shell:

$ cd PDAnalysis

$ git init

$ git remote add origin git@github.com:mgalanti/PDAnalysis.git

$ git fetch

$ git reset origin/2017V12-CMSSW_10_2_6

$ git checkout -t origin/2017V12-CMSSW_10_2_6

$ git reset --hard HEAD
