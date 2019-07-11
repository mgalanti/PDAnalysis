#!/usr/bin/env python

from ROOT import TMVA, TFile, TTree, TCut
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
file = '../../../../datasets/EleMVASecondNtuple__BsToJpsiPhiDG0_2018_DCAP__BsToJPsiPhi_eleTagLooseV1__20190703_190644__11M_Missing.root'

data = TFile.Open(file)

tree = data.Get('EleMVAsecondTree')

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
    ,('eleIDNIV2Val', 'F')
    #,('eleIDNIV2Cat', 'F')
    ,('eleDRB', 'F')
    ,('elePFIsoScaled', 'F')
    ,('eleConeCleanPt', 'F')
    ,('eleConeCleanPtRel', 'F')
    ,('eleConeCleanDR', 'F')
    ,('eleConeCleanEnergyRatio', 'F')
    ,('eleConeCleanQ', 'F')
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
cut = 'JPsiMuHltBit==0 && JPsiTrkTrkHltBit==1 && eleSelected==1'

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

nBkg = '500000'
nSgn = '600000'
nBkgTest = '46000'
nSgnTest = '82000'


dataloaderOpt = 'nTrain_Signal=' + nSgn + ':nTrain_Background=' + nBkg + ':nTest_Signal=' + nSgnTest + ':nTest_Background=' + nBkgTest
dataloaderOpt += ':SplitMode=Random:NormMode=NumEvents:V:'

dataloader.PrepareTrainingAndTestTree(TCut(cutSgn), TCut(cutBkg), dataloaderOpt)

# Create Keras Model
nLayers = 3
layerSize = 200
dropValue = 0.4

modelName = getKerasModel(nVars, 'model' + name + '.h5', nLayers, layerSize, dropValue)
# modelName = 'TrainedModel_DNNOsMuonHLTJpsiMu.h5'
# Book methods
#dnnOptions = '!H:!V:NumEpochs=50:TriesEarlyStopping=10:BatchSize=1024:ValidationSize=33%:SaveBestOnly=True'
#dnnOptions = dnnOptions + ':Tensorboard=./logs:FilenameModel=' + modelName
dnnOptions = '!H:!V:NumEpochs=50:TriesEarlyStopping=10:BatchSize=1024:SaveBestOnly=True'
dnnOptions = dnnOptions + ':FilenameModel=' + modelName
# 

# Preprocessing string creator, loop was for selection of which variable to apply gaussianification
preprocessingOptions = ':VarTransform=N,G,N'

dnnName = 'DNN' + name

factory.BookMethod(dataloader, TMVA.Types.kPyKeras, dnnName, dnnOptions + preprocessingOptions)

# Run training, test and evaluation
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
