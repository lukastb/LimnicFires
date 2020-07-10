# LimnicFires
Effects of wildfire char on dissolved organic matter and  microbial functioning in a fluvial ecosystem

## Summary
Field experiemnt with control/treatment desing in a austrian river carried out in streamside flumes. 5 flumes acted as control and 5 as treatment, i.e. with wildfire char addition. Tiles with biofilm grown locally were added at the beginning of the experiment. During the experiment sampels for DOC and DOM were taken and at the end of the experiemnt (8h) the biofilm on the tiles was scapred off and extracellular enzyme activities (EEA) were measured.

## Technologies
- R 3.5.3
- rstan 2.19.2

## Aim
Build a bayesian model that accounts for (i) correlation over time for variabels like DOM properties (DOM indices; absorbance and fluorescence based) and DOC concentration and (ii) effects of treatment (i.e. flumes with wildfire char vs. flumes without wildfire char)

## Data
- DOC data over time (at the beginning, 2, 4, 6 and 8 hours after the beginning) for control and treatment
- DOM data over time (at the beginning, 2, 4, 6 and 8 hours after the beginning) for control and treatment (only some variables are of interest, i.e. SUVA254)
- EEA data at the end of the experiment for control and treatment
