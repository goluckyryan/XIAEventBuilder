# XIAEventBuilder
event builder for XIA PIXIE16 digitizer, from evt to root

## prerequiest
ROOT 6.00+

Codes need to be compiled are in armory/

# armory/DataBlock.h
this is the source file for the class DataBlock, it stored all information from a single data block from pixie16 output. 

# armory/evtReader.h
this is the source file for the class evtReader.
It read the *.evt file (which is same as pixie16 output) and convert each measurement (or data block) from byte into meaningful data and use DataBlock class to store the information.
It can also scan the evt file.

# armory/MergeEVT
this merges all evt files into *_raw.root. 

# armory/EventBuilder
this builds events from *_raw,root to *.root file

# armory/evt2hist
this processes evt file to hstograms. 

# armory/to2root
this build events from *.evt.to files to *.root file (need to check the compactability with EventBuilder)

# armory/pxi-fsu-time-order
this sorting the time from *evt file to *.evt.to.fsu.XXX, where XXX is the time window.
It will search data within XXX time window, if non of the data is from clover, discard.

# armory/xia2root
this is old evt to root for custom pixie DAQ. 

# Analyzer.C/h
this is a TSelector for analysis the *.root file

# PreAnalyzer.C/h
this is a TSelector for consolidate the *.root and save to new *.root file.
It only save added-back gamma and GAGG data.

# PIDAnalyzer.C/h
this is a TSelector for analysing particle-gamma

# GAGGPIDCutCreator.C
this is a PID Cuts Creator
