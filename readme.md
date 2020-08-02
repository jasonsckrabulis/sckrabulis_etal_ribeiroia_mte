
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

TBD

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
Stage | Gosner stage of the tadpole
Mass | Tadpole mass measured in grams 
PerfRound | Number signifying which of the three sequential performance temperatures the measurement was taken
PerfTemp | Temperature at which respiration measurements were taken
Time | Duration of the respiration period, calculated from start and stop times in minutes
O2 | Change in dissolved oxygen, calculated by the difference between the initial concentration of oxygen in the control jars and the final measurement from each tadpole, in milligrams per liter
O2/Time | Rate of oxygen consumed per minute, calculated by dividing O2 by Time
O2/Time/Mass | Mass-corrected rate of oxygen consumed per minute per tadpole mass, calculated by dividing O2/Time by Mass (units: g/L oxygen per min per g tadpole)
corO2/Time/Mass | Corrected O2/Time/Mass value to obtain mass-corrected rate of oxygen consumed per minute controlling for jar volue, calculated by multiplying O2/Time/Mass by volume of jar 0.57L (units: g oxygen/min/g tadpole)

**Cerc swimming speed.csv**

Variable name | Description
--- | ---
Snail | Number signifying which snail the cercaria originated from
Temperature | Temperature at which the cercariae were recorded
avgSpeed | Average swimming speed for each cercaria, calculated in pixels per second within imageJ by the wrMTrck plugin
avgSpeedMM | Average swimming speed for each cercaria in millimeters per second, calculated using the conversion factor 94.5 pixels per millimeter
