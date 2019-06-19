# PDAnalysis
Analysis code based on the PDAnalysis framework
The repository only contains differences with respect to the standard PDAnalysis package provided by Paolo.
It does not run standalone and needs to be installed on top of full PDAnalysis.

Instructions to install on top of PDAnalysis 2017V13_02 on CMSSW_10_3_0:

- Use Paolo's installer to setup a write area of git@github.com:mgalanti/NtuTool.git under CMSSW_10_3_0.

Additional steps:

$ cd ${CMSSW_BASE}/src

- Install MGTools:

$ git clone git@github.com:mgalanti/MGTools.git

$ git checkout V13_02-CMSSW_10_3_0-v01-00

- Install NtuTool:

$ git clone git@github.com:mgalanti/NtuTool.git

$ git checkout V13_02-CMSSW_10_3_0-v01-00

- Install NtuAnalysis:

$ git clone git@github.com:mgalanti/NtuAnalysis.git

$ git checkout V13_02-CMSSW_10_3_0-v01-01

- Now, make sure that you downloaded the write version of PDAnalysis according to Paolo installer:
***WARNING: the commands below can be destructive and you may lose any changes you made to the PDAnalysis directory. Try only on a pristine installation of PDAnalysis**

$ cd PDAnalysis

$ git init

$ git remote add origin git@github.com:mgalanti/PDAnalysis.git

$ git fetch

$ git reset origin/master

$ git reset --hard HEAD

$ git checkout V13_02-CMSSW_10_3_0-v01-01

After this, you will end up with a lot of untracked files inside PDAnalysis. This is expected (the git repository only tracks the differences with respect to a standard PDAnalysis installation).
