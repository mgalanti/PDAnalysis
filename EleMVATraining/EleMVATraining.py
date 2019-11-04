#!/usr/bin/env python

from subprocess import call
from os.path import isfile
import sys
import os
os.environ['TF_CPP_MIN_LOG_LEVEL']='2'

import numpy as np
import keras
import h5py

from keras.models import Sequential
from keras.layers import Dense, Activation, Conv2D, MaxPooling2D, Flatten, Dropout
from keras.regularizers import l2
from keras.optimizers import SGD, Adam
from keras.callbacks import ModelCheckpoint, TensorBoard

from ROOT import TMVA, TFile, TTree, TCut, TChain


##### FUNCTIONS

def getKerasModel(inputDim, modelName, nLayers = 3, layerSize = 200, dropValue = 0.2):
    model = Sequential()
    model.add(Dense(layerSize, activation='relu', kernel_initializer='normal', input_dim=inputDim))
    if dropValue != 0:
        model.add(Dropout(dropValue))

    for i in range(1, nLayers):
        model.add(Dense(layerSize, activation='relu', kernel_initializer='normal'))
        if dropValue != 0:
            model.add(Dropout(dropValue))


    TensorBoard(log_dir='./logs', histogram_freq=0, batch_size=32, write_graph=True, write_grads=False, write_images=False, embeddings_layer_names=None, embeddings_metadata=None)
    # Anything below this point should not be changed in order for the network to work with TMVA, exception: the optimizer
    model.add(Dense(2, activation='softmax'))

    opt = Adam(lr=0.001)
    model.compile(loss='categorical_crossentropy', optimizer=opt, metrics=['accuracy'])
    model.save(modelName)
    model.summary()
    return modelName


# Setup TMVA
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()

# Load data
filenames = [
  '../../../../datasets/EleMVASecondNtuple__BsToJpsiPhiDG0_2018_DCAP__BsToJPsiPhi_eleTagLooseV1__20190808_141613__1.root'
  #, "../../../../datasets/EleMVASecondNtuple__BsToJpsiPhi_2018_DCAP__BsToJPsiPhi_eleTagLooseV1__20190808_140620__1.root"
  #, "../../../../datasets/EleMVASecondNtuple__BsToJpsiPhi_2017_DCAP__BsToJPsiPhi_eleTagLooseV1__20190808_141410__1.root"
  ]

#data = TFile.Open(file)

tree = TChain("EleMVAsecondTree")
for filename in filenames:
  print("Adding file")
  print filename
  #file = TFile.Open(filename)
  tree.Add(filename)
  
nentries = tree.GetEntries()

print "Number of entries: "
print nentries

#tree = data.GetTree()

# Prepare factory
name = 'OsElectronHLTJpsiTrkTrk'

outputName = 'TMVA' + name + '.root'

output = TFile.Open(outputName, 'RECREATE')

factory = TMVA.Factory('TMVAClassification', output,
                       '!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification')

dataloader = TMVA.DataLoader('dataset')

# variable list
varList = [
    ('elePt', 'F')
    ,('eleEta', 'F')
    ,('eleDxy', 'F')
    ,('eleExy', 'F')
    ,('eleDz', 'F')
    ,('eleEz', 'F')
    ,('eleIDNIV2Val', 'F')
    ,('eleIDNIV2Cat', 'I')
    ,('eleDRB', 'F')
    ,('elePFIsoScaled', 'F')
    ,('eleConePt', 'F')
    ,('eleConePtRel', 'F')
    ,('eleConeDR', 'F')
    ,('eleConeEnergyRatio', 'F')
    ,('eleConeQ', 'F')
    ,('eleConeCF', 'F')
    ,('eleConeNCH', 'I')
    # ,('muoJetDFprob', 'F')
    ]

varListClean = [
    ('elePt', 'F')
    ,('eleEta', 'F')
    ,('eleDxy', 'F')
    ,('eleExy', 'F')
    ,('eleDz', 'F')
    ,('eleEz', 'F')
    #,('eleIDNIV2RawVal', 'F')
    ,('eleIDNIV2Val', 'F')
    #,('eleIDNIV2Cat', 'F')
    ,('eleDRB', 'F')
    ,('elePFIsoScaled', 'F')
    ,('eleConeCleanPt', 'F')
    ,('eleConeCleanPtRel', 'F')
    ,('eleConeCleanDR', 'F')
    ,('eleConeCleanEnergyRatio', 'F')
    ,('eleConeCleanQ', 'F')
    #,('eleDxy/eleExy-eleConeCleanAvgDxy/eleConeCleanStdDevDxy', 'F')
    #,('eleDz/eleEz-eleConeCleanAvgDz/eleConeCleanStdDevDz', 'F')
    #,('sqrt(eleConeCleanAvgDxy*eleConeCleanAvgDxy+eleConeCleanAvgDz*eleConeCleanAvgDz)', 'F')
    #,('eleConeCleanAvgDxy', 'F')
    #,('eleConeCleanAvgDz', 'F')
    #,('eleConeCleanStdDevDxy', 'F')
    #,('eleConeCleanStdDevDz', 'F')
    ]

varList = varListClean

# automatic variable counting and adding
nVars = 0
for var in varList:
    dataloader.AddVariable( var[0], var[1] )
    nVars += 1

# prepare dataloader
# Event wise selection
#cut = 'JPsiMuHltBit==0&&eleIDNIV2Val>-0.98'
#cut = 'JPsiMuHltBit==0 && JPsiTrkTrkHltBit==1 &&eleSelected==1&&eleDxy>-999&&eleDxy<999&&eleConeCleanAvgDxy>-999&&eleConeCleanAvgDxy<999&&eleConeCleanStdDevDxy>0&&eleConeCleanStdDevDxy<999&&eleIDNIV2Val>-0.999&&fabs(eleDz)<0.5&&elePt>2.5'

cut = 'JPsiMuHltBit==0 && JPsiTrkTrkHltBit==1 &&eleSelected==1&&eleDxy>-999&&eleDxy<999&&eleConeCleanAvgDxy>-999&&eleConeCleanAvgDxy<999&&eleConeCleanStdDevDxy>0&&eleConeCleanStdDevDxy<999&&eleIDNIV2Val>-0.9999&&fabs(eleDz)<0.2&&fabs(eleDxy)<0.1&&elePt>2.5&&eleDRB>0.4'


cutSgn = cut + '&&tagTruth==1' #Correcly tagged events selection i.e. sign(charge lepton) -> correct flavour 
cutBkg = cut + '&&tagTruth==0' #Uncorrecly tagged events

# same tree, add selection later
dataloader.AddSignalTree(tree)
dataloader.AddBackgroundTree(tree)

# evtWeight variable in the ntuple, address simulation bias
dataloader.SetWeightExpression( 'evtWeight' );

# nBkg = 169226
# nSgn = 393356

#nBkg = '126920'
#nSgn = '295017'
#nBkgTest = '42306'
#nSgnTest = '98339'

## 2018 + 2018DG0 + 2017
#nSgnTot = 530065
#nBkgTot = 366048

# 2018DG0
nSgnTot = 296739
nBkgTot = 202240

nSgn = nSgnTot // 1.25
nSgnTest = nSgnTot - nSgn

nBkg = nBkgTot // 1.25
nBkgTest = nBkgTot - nBkg

nBkg = str(nBkg)
nSgn = str(nSgn)
nBkgTest = str(nBkgTest)
nSgnTest = str(nSgnTest)

dataloaderOpt = 'nTrain_Signal=' + nSgn + ':nTrain_Background=' + nBkg + ':nTest_Signal=' + nSgnTest + ':nTest_Background=' + nBkgTest
dataloaderOpt += ':SplitMode=Random:NormMode=NumEvents:V:'

dataloader.PrepareTrainingAndTestTree(TCut(cutSgn), TCut(cutBkg), dataloaderOpt)

# Create Keras Model
nLayers = 3
layerSize = 200
dropValue = 0.5

modelName = getKerasModel(nVars, 'model' + name + '.h5', nLayers, layerSize, dropValue)
# modelName = 'TrainedModel_DNNOsMuonHLTJpsiMu.h5'
# Book methods
dnnOptions = '!H:!V:NumEpochs=100:TriesEarlyStopping=10:BatchSize=1024:ValidationSize=20%:SaveBestOnly=True'
dnnOptions = dnnOptions + ':Tensorboard=./logs:FilenameModel=' + modelName
#dnnOptions = '!H:!V:NumEpochs=100:TriesEarlyStopping=10:BatchSize=1024:SaveBestOnly=True'
#dnnOptions = dnnOptions + ':FilenameModel=' + modelName
# 

# Preprocessing string creator, loop was for selection of which variable to apply gaussianification
preprocessingOptions = ':VarTransform=N,G,N'

dnnName = 'DNN' + name

factory.BookMethod(dataloader, TMVA.Types.kPyKeras, dnnName, dnnOptions + preprocessingOptions)

# Run training, test and evaluation
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
