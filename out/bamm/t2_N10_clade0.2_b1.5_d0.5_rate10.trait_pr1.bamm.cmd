# MCMC parameters
sampleFromPriorOnly = 0                 
runMCMC = 1                             
simulatePriorShifts = 0
loadEventData = 0                       
eventDataInfile = event_data_in.txt
initializeModel = 1                     
overwrite = 0
traitPriorMin = 0
traitPriorMax = 0
betaIsTimeVariablePrior = 1
acceptanceResetFreq = 1000
updateBetaInitScale = 1
updateBetaShiftScale = 1
updateNodeStateScale = 1
updateEventLocationScale = 0.05
updateEventRateScale = 4.0
updateRateEventNumber = 1
updateRateEventPosition = 1
updateRateEventRate = 1
updateRateBeta0 = 1
updateRateBetaShift = 1
updateRateNodeState = 25
updateRateBetaTimeMode = 0
localGlobalMoveRatio = 10.0
betaInit = 0.5
betaShiftInit = 0
initialNumberEvents = 0

# heating might need some adjustment
modeltype = trait
deltaT = 0.01
swapPeriod = 1000
treefile = data/t2_N10_clade0.2_b1.5_d0.5.nwk
traitfile = data/t2_N10_clade0.2_b1.5_d0.5_rate10.trait
mcmcWriteFreq =  1000
printFreq =  1000
eventDataWriteFreq =  10
numberOfGenerations = 10000
numberOfChains = 2
mcmcOutfile = out/bamm/t2_N10_clade0.2_b1.5_d0.5_rate10.trait_pr1_mcmc_out.txt
eventDataOutfile = out/bamm/t2_N10_clade0.2_b1.5_d0.5_rate10.trait_pr1_event_data.txt
runInfoFilename = out/bamm/t2_N10_clade0.2_b1.5_d0.5_rate10.trait_pr1_run_info.txt
chainSwapFileName = out/bamm/t2_N10_clade0.2_b1.5_d0.5_rate10.trait_pr1_chain_swap.txt
###############################################

# Prior block chosen by BAMMtools::setBAMMpriors

expectedNumberOfShifts = 1.0

betaInitPrior = 0.222099409572233

betaShiftPrior = 0.328527612920643

useObservedMinMaxAsTraitPriors = 1

#### End Prior block
######################


