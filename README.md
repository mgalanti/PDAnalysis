# PDAnalysis
Analysis code based on the PDAnalysis framework
The repository only contains differences with respect to the standard PDAnalysis package provided by Paolo.
It does not run standalone and needs to be installed on top of full PDAnalysis.

Instructions to install on top of PDAnalysis 2017V13_02 on CMSSW_10_3_0:

- Use Paolo's installer to setup a write area of git@github.com:mgalanti/NtuTool.git under CMSSW_10_3_0:

$ ntuAna/NtuProdV13_02/build_write CMSSW_10_3_0 SCRAM_ARCH=slc6_amd64_gcc700 CMSSW_VERSION=CMSSW_10_3_0 FIX_BPHAnalysis=A10 FIX_NtuAnalysis=764 FIX_PDAnalysis=A20

Additional steps:

$ cd ${CMSSW_BASE}/src

- Install MGTools:

$ git clone git@github.com:mgalanti/MGTools.git

$ git checkout V13_02-CMSSW_10_3_0-v01-00

- Install NtuTool:

$ rm -r NtuTool # Or move it outside of src/ if you want to keep it

$ git clone git@github.com:mgalanti/NtuTool.git

$ git checkout V13_02-CMSSW_10_3_0-v01-00

- Install NtuAnalysis:

$ rm -r NtuAnalysis # Or move it outside of src/ if you want to keep it

$ git clone git@github.com:mgalanti/NtuAnalysis.git

$ git checkout V13_02-CMSSW_10_3_0-v01-01

- Now, make sure that you downloaded the write version of PDAnalysis according to Paolo installer:
***WARNING: the commands below can be destructive and you may lose any changes you made to the PDAnalysis directory. Try only on a pristine installation of PDAnalysis**

# Do *NOT* rm -r PDAnalysis !

$ cd PDAnalysis

$ git init

$ git remote add origin git@github.com:mgalanti/PDAnalysis.git

$ git fetch

$ git reset origin/master

$ git reset --hard HEAD

$ git checkout V13_02-CMSSW_10_3_0-v01-01

After this, you will end up with a lot of untracked files inside PDAnalysis. This is expected (the git repository only tracks the differences with respect to a standard PDAnalysis installation).
