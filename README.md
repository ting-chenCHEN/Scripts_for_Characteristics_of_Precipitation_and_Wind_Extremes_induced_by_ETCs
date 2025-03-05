# Scripts_for_characteristics_of_Precipitation_and_Wind_Extremes_induced_by_ETCs

This repository contains the scripts used in the paper "Characteristics of Precipitation and Wind Extremes Induced by Extratropical Cyclones in Northeastern North America" by Chen and Di Luca (2025), published in the Journal of Geophysical Research - Atmospheres. For detailed methods of data processing and analysis, please refer to the paper. 

Please also note that the data files that are
(a) created/output by the scripts (**Analysis_Scripts**) and 
(b) used as input for the plotting scripts (**Figure_Scripts**) 
are publicly available on the Federated Research Data Repository (doi: 10.20383/103.01231).

--------------------
GENERAL INFORMATION
--------------------

1-Authors Information: Ting-Chen Chen (a,b) and Alejandro Di Luca (a)

2-Affiliation: 
	a) Centre pour l’Étude et la simulation du climat à l’échelle régionale (ESCER) | Département des sciences de la Terre et de l’atmosphère. Université du Québec à Montréal (UQAM), Montréal, Canada. 
	b) Moody’s, London, United Kingdom.

3- Information about funding sources that supported the collection of the data:

This research has been made as part of the project “Simulation et analyse du climat à haute résolution” thanks to the financial participation of the Government of Québec. This research was enabled in part by support provided by Calcul Québec (https://www.calculquebec.ca) and the Digital Research Alliance of Canada (https://alliancecan.ca). 

A. Di Luca was funded by the Natural Sciences and Engineering Research Council of Canada (NSERC) grant (RGPIN-2020-05631). 

--------------------------------------------------
SHARING/ACCESS INFORMATION
--------------------------------------------------

1-Licenses/Restrictions placed on the data: 
These data are available under a CC BY 4.0 license <https://creativecommons.org/licenses/by/4.0/> 

2-Recommended citation for this repository:
Chen, T.-C., & Di Luca, A. (2025). ETC-induced Precipitation and Wind Extremes analysis (v0.1.0). Zenodo. https://doi.org/10.5281/zenodo.14976915

-----------------------------
FILE OVERVIEW
-----------------------------

Some of the acronyms used in file naming are listed below, but please refer to Chen and Di Luca (2025) for other definitions.

-**WS** refers to 10-m Wind Speeds.

-**TP** refers to total surface precipitation.

-**P98p0** refers to extremes defined using the 98.0 percentile.

-**P99p0** refers to extremes defined using the 99.0 percentile.

-**P99p9** refers to extremes defined using the 99.9 percentile.

-**ERA5** refers to calculations made using ERA5 reanalysis data.

-**IMERG** refers to calculations made using IMERG satellite-based data.

-**NNA** refers to the Northeast North American region (see e.g., Fig. 3 in the paper).

-**Compound** extremes refer to extremes showing simultaneously wind AND precipitation extremes within a 12h window.

-**timeavg** refers to storm extremes calculated using the average Extreme Exceedance metric.

-**timecum** refers to storm extremes calculated using the cumulated Extreme Exceedance metric.

-----------------------------
1- Folder **Analysis_Scripts**:

•	*Find_ExtremesThresholds_TPexample.py*

Finds the local extreme WS or TP threhold values based on the grid-point 98th, 99th, or 99.99th-percentile climatology.
  
•	*Find_local_extreme_records_ETCassocation_ERA5_TP.f90*

For every grid points, records the (a) date time, (b) magnitude, and (c) whether or not it is associated with one or more ETCs, of every hourly location extreme occurrence.

•	*Calcu_averaged_exceedance_ERA5_TP.f90*

For every grid points, calculates the averaged local extreme exceedance on the seasonal or annual scales. 

•	*Calcu_CompoundExtremes_prob_ETCassociation_relaxto12h.f90*

For every grid points, calculates the number of compound extreme events (counted if one hourly precipitation extreme occurs 6hrs prior to 6 hrs after a local wind extreme takes places at the same location). The calculation is separated into ETC-associated or non-ETC associated compounds.  

•	*Obtain_local_ex_timeseries_Montreal_during_HalloweenStorm.f90*

Taking a local grid point in Montreal as a "reference station", records the time series of WS, TP during the passage of the 2019 Halloween Storm, along with the storm intensity and the distance between the station and the storm's center.

•	*TPExtremes_duration_exceedance_maps_HalloweenStorm.f90*

Obtains the TP extreme duration and exceedance (maps and PDF) during the event of 2019 Halloween Storm. 

•	*Calcu_PDF_TPexceedance_duration_ETCs_NNA.f90*

Calculates the PDF of TP exceedance, PDF of TP extreme duration, and the averaged TP exceedance as a function of extreme duration for all ETCs in NNA. 

•	*Calcu_Stormlist_extreme_in_NNA_timecum.f90*

Given a list of all ETC tracks (with varied lifetime), calculates their ETC-associated WS & TP extreme impacts in the NNA region, cumulated over each ETC's lifetime.

•	*Calcu_Stormlist_extreme_in_NNA_timeavg.f90*

Given a list of all ETC tracks (with varied lifetime), calculates their ETC-associated WS & TP extreme impacts in the NNA region, averaged over each ETC's lifetime.

-----------------------------

2- Folder **Figure_Scripts**:

• Fig1:  *Fig1_Fig6.ipynb*

•	Fig2:  *Fig2_Pextreme_ERA5vsISD.py*, *Fig2_Wextreme_ERA5vsISD.py*, *Fig2_Pextreme_IMERGvsISD.py* 
         
•	Fig3:  *Fig3_FirstColumn_Wextreme_seasonal_freq.py*, *Fig3_SecondColumn_Wextreme_seasonal_ETCfreq.py*, *Fig3_ThirdColumn_Wextreme_ETC_ratio.py*

•	Fig4:  *Fig4_FirstColumn_Pextreme_seasonal_freq.py*, *Fig4_SecondColumn_Pextreme_seasonal_ETCfreq.py*, *Fig4_ThirdColumn_Pextreme_ETC_ratio.py*

•	Fig5:  *Fig5.ipynb*

•	Fig6:  *Fig1_Fig6.ipynb*

•	Fig7:  *Fig7_ab.py*, *Fig7_cd.py*, *Fig7_efg.py*

•	Fig8:  *Fig8.ipynb*

•	Fig9:  *Fig9_a.ipynb*, *Fig9_b.ipynb*, *Fig9_c.ipynb*

•	Fig10: *Fig10_ac.ipynb*, *Fig10_bd.ipynb*

•	*funcs_map.py*




