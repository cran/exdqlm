# exdqlm 0.1.3
- exdqlm 0.1.2 was archived on 2022-10-23 as requires archived package 'dlm'. 'dlm' now back on CRAN.
- Fixes and general improvements
  - changed if/class conditions in exdqlmISVB, exdqlmMCMC to is()
  - added the function is.exdqlm() to utils
  - changed if/class conditions in compPlot, exdqlmChecks, exdqlmForcast, exdqlmPlot to is.exdqlm()
  
# exdqlm 0.1.2
- Fixes and general improvements
  - fixed a bug in the forecast plotting routine
  - dqlm.ind = TRUE automatically sets gam.init = 0 & fix.gamma = TRUE (implemented in check_logics)
  - added a check that gam.init is in the appropriate range in exdqlmISVB and exdqlmMCMC
- Input changes
  - added parameter to specify the percentage of the CrIs in compPlot, exdqlmPlot and exdqlmForecast
- Return changes
  - kt added to transfn_exdqlmISVB return
- Dataset changes
  - Niño 3.4 dataset added
  - BTflow dataset updated
  - monELI dataset removed

# exdqlm 0.1.1
- Documentation updates
  - added author (year) to citation in descriptions
- Argument changes
  - added input verbose to the functions exdqlmISVB & exdqlmMCMC to allow users to suppress progress updates
- Return changes
  - functions exdqlmPlot and compPlot now return a list of the MAP & CrI estimates that are plotted

# exdqlm 0.1.0
- First release.
