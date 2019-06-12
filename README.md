# PDAnalysis
Analysis code based on the PDAnalysis framework
The repository only contains differences with respect to the standard PDAnalysis package provided by Paolo.
It does not run standalone and needs to be installed on top of full PDAnalysis.
Other dependencies:
- MGTools
- NtuTool
- NtuAnalysis

Instructions to install on top of PDAnalysis 2017V13_02 on CMSSW_10_3_0:

*** W A R N I N G ***  Porting to 2017V13_02 on CMSSW_10_3_0 is still ongoing! 
                       Code does not work yet.

- Download PDAnalysis according to Paolo installer

In the shell:

$ cd PDAnalysis

$ git init

$ git remote add origin git@github.com:mgalanti/PDAnalysis.git

$ git fetch

$ git reset origin/master

$ git reset --hard HEAD
