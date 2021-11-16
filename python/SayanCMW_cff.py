import FWCore.ParameterSet.Config as cms

from Analyzers.SayanCMW.SayanCMW_cfi import *

#### Standard analysis for pixel Rereco PbPb 2015 ####

#### centrality 0-5% ####
CPDC05 = defaultCPDC.clone()
CPDC05.centmin = cms.untracked.int32(0)
CPDC05.centmax = cms.untracked.int32(5) ##5
#### centrality 5-10% ####
CPDC510 = defaultCPDC.clone()
CPDC510.centmin = cms.untracked.int32(5)
CPDC510.centmax = cms.untracked.int32(10)
#### centrality 10-20% ####
CPDC1020 = defaultCPDC.clone()
CPDC1020.centmin = cms.untracked.int32(10)
CPDC1020.centmax = cms.untracked.int32(20)
#### centrality 20-30% ####
CPDC2030 = defaultCPDC.clone()
CPDC2030.centmin = cms.untracked.int32(20)
CPDC2030.centmax = cms.untracked.int32(30)
#### centrality 30-40% ####
CPDC3040 = defaultCPDC.clone()
CPDC3040.centmin = cms.untracked.int32(30)
CPDC3040.centmax = cms.untracked.int32(40)
#### centrality 40-50% ####
CPDC4050 = defaultCPDC.clone()
CPDC4050.centmin = cms.untracked.int32(40)
CPDC4050.centmax = cms.untracked.int32(50)
#### centrality 50-70% ####
CPDC5070 = defaultCPDC.clone()
CPDC5070.centmin = cms.untracked.int32(50)
CPDC5070.centmax = cms.untracked.int32(70)
#### centrality 70-90% ####
CPDC70100 = defaultCPDC.clone()
CPDC70100.centmin = cms.untracked.int32(70)
CPDC70100.centmax = cms.untracked.int32(100)
