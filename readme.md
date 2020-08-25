
## Sckrabulis _et al._ - Using metabolic theory to describe temperature and thermal acclimation effects on parasitic infection

### Authors and maintainers

Jason P. Sckrabulis: jason.sckrabulis@gmail.com

Karie A. Altman: karie.altman@gmail.com

Thomas R. Raffel: raffel@oakland.edu

>Conceptualization: JPS, KAA & TRR. Methodology: JPS, KAA, & TRR. Analysis: JPS & TRR. Writing - review & editing: JPS, KAA, & TRR. Funding: TRR. Literature searching: JPS.

### Issues and suggestions

Please email jason.sckrabulis@gmail.com with any issues or suggestions with the data or code supplement, or submit via the GitHub issues tab for this [repository](  https://github.com/jasonsckrabulis/sckrabulis_etal_ribeiroia_mte/issues).
For questions concerning the manuscript, please email the corresponding author at jason.sckrabulis@gmail.com.

### Change log

* July 1, 2020: First full commit

### Citation

Please cite work as:
>TBD
---

### Abstract

Predicting temperature effects on species interactions can be challenging, especially for parasitism where it is difficult to experimentally separate host and parasite thermal performance curves. Prior authors proposed a possible solution based on the metabolic theory of ecology (MTE), using MTE-based equations to describe the thermal mismatch between host and parasite performance curves and account for thermal acclimation responses. Here we use published infection data, supplemented with experiments measuring metabolic responses to temperature in each species, to show that this modeling framework can successfully describe thermal acclimation effects on two different stages of infection in a tadpole-trematode system. All thermal acclimation effects on host performance manifested as changes in one key model parameter (activation energy), with measurements of host respiration generating similar MTE parameter estimates and acclimation effects compared to measurements of the host’s ability to clear encysted parasites. This result suggests that metabolic parameter estimates for whole-body metabolism can sometimes be used to estimate temperature effects on host and parasite performance curves. However, we found different thermal patterns for measurements of host resistance to initial parasite encystment. This result emphasizes potential challenges when applying MTE-based models to complex parasite-host systems with multiple distinct stages of infection.

---

### Repository Contents

* README.md
* data  
   Folder of experimental data as .csv files used in data analysis for manuscript (Please note that metacercaria encystment and clearance data is from Altman et al. 2016 and can be found on [Dryad](https://datadryad.org/stash/dataset/doi:10.5061/dryad.f3k8p) 
   * Uninfected tadpole respiration.csv (dataset for host respiration acclimation experiment)
   * Cerc swimming speed.csv (dataset for cercaria swimming speed experiment)
* code  
   Folder of statistical R code used to analyze data
---

### Variable descriptions

**Uninfected tadpole respiration.csv**

Variable name | Description
--- | ---
Type | Label for respiration jars in the experiment where `tad` = tadpole for measurement & `control` = an empty control jar
ID | Jar label for internal tracking of location
Block | Number representing which block the jar was in
AccInc | Number signifying which acclimation incubator the animal was maintained in during the acclimation period
AccTemp | Temperature at which the tadpole was acclimated measured in Celsius
Stage | Gosner stage of the tadpole (based on **Gosner, K. 1960. A simplified table for staging anuran embryos and larvae with notes on identification. *Herpetologica* 16(3): 183-190**)
Mass | Tadpole mass measured in grams 
PerfRound | Number signifying which of the three sequential performance temperatures the measurement was taken
PerfTemp | Temperature at which respiration measurements were taken
Time | Duration of the respiration period, calculated from start and stop times in minutes
O2 | Change in dissolved oxygen, calculated by the difference between the initial concentration of oxygen in the control jars and the final measurement from each tadpole, in milligrams per liter
O2/Time | Rate of oxygen consumed per minute, calculated by dividing O2 by Time
O2/Time/Mass | Mass-corrected rate of oxygen consumed per minute per tadpole mass, calculated by dividing O2/Time by Mass (units: g/L oxygen per min per g tadpole)
corO2/Time/Mass | Corrected O2/Time/Mass value to obtain mass-corrected rate of oxygen consumed per minute controlling for jar volue, calculated by multiplying O2/Time/Mass by volume of jar 0.57L (units: g oxygen per min per g tadpole)

**Cerc swimming speed.csv**

Variable name | Description
--- | ---
Snail | Number signifying which snail the cercaria originated from
Temperature | Temperature at which the cercariae were recorded
avgSpeed | Average swimming speed for each cercaria, calculated in pixels per second within imageJ by the wrMTrck plugin
avgSpeedMM | Average swimming speed for each cercaria in millimeters per second, calculated using the conversion factor 94.5 pixels per millimeter
