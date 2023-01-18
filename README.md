# hydroVolc v1.1
Repository for matlab glacio/hydrovolcanism 1D coupled model

THIS MODEL HAS NOT BEEN APPROVED FOR PUBLIC DISTRIBUTION BY ALL AUTHORS. PLEASE DO NOT DISTRIBUTE WITHOUT CONSENT.

This repository contains code for the coupled model originally published as
"Rowell, Colin R., A. Mark Jellinek, Sahand Hajimirza, and Thomas J. Aubry. “External Surface Water Influence on Explosive Eruption Dynamics, With Implications for Stratospheric Sulfur Delivery and Volcano-Climate Feedback.” Frontiers in Earth Science 10 (2022). https://www.frontiersin.org/article/10.3389/feart.2022.788294. "

The coupled model consists of three coupled components:
  (1) 1D conduit model for rhyolitic melt (from Hajimirza et al, 2021)
  	Main function: "Conduit_flow_with_nucleation_V6.m"
  (2) 1D model for a subaqueous volcanic jet (Rowell et al., 2022)
  	Main function: "MWIv2.m"
  		- (MWIv1 is deprecated and won't work, v3 is an early version which includes thermal disequilibrium between jet phases, but is not fully functioning)
  (3) 1D model for a subaerial plume (Rowell et al., 2022, and modified from Degruyter and Bonadonna, 2012)
  	Main function: "hmodel.m"

It can be run in full mode with all three components, or with only components 2 and 3, using a "proxy" model to calculate the vent condition instead of running the full conduit model.
See "model_Tutorial_Runs.m" to get started with running the coupled conduit/plume model.

A complete set of model input parameters can be found in the following three scripts:
	getConduitSource.m 			Populates input parameters for the 1D conduit model
	getProxyConduitSourceV2.m 	Populates input parameters for "proxy" conduit (Vent condition with no conduit model run)
	getPlumeSource.m 			Populates input parameters for plume model

atmoFiles/ directory could be properly moved to a separate data directory
	- atmprofile.mat is a generic mid-latitude profile (I think?)
	- the "absWind" versions have absolute (direction removed) wind profiles included. I can include a script to download ERA5 atmospheric profiles if/when it is useful

Overview of main folders:
	Hajimirza_Conduit	Functions of Hajimirza condiut model
	vent				Functions for subaqueous jet model
	water95				IAPWS95 water thermodynamics package (Junglas, 2009)
	utils				Miscellaneous useful functions for quick calculations
	1Dplume_DB2012		Functions of 1D plume model
	atmoFiles			Sample atmospheric profiles
	modelSweeps			Functions to run parameter sweeps and Monte Carlo simulations


Model Sweep Functions:
	conduitPressureSweep.m    	Script used to run "shooting" sweeps of the conduit model to find an appropriate radius for a given mass flux

	conduitPressureSweepQ.m   	Find adjusted Mass Eruption Rate in conduit model for adjusted surface pressures

	getCoupledSweepArrays.m   	Get summary data arrays for parameters sweeps and monte carlo runs

	getHydroVolcSweepSummary.m 	Wrapper function for getCoupledSweepArrays that automatically saves output

	hydroVolcMC.m 		      Function to run monte carlo sweep, optional parallel pool

	runConduitWaterSweep.m    Run a proxy model sweep, varying the mass fraction of conduit water infiltration.

	runCoupledSweep.m 	      Script originally used to run model sweeps over water depth and mass eruption rate



CHANGELOG:
v1.1 December 14, 2022
	- Monte Carlo simulations are now fully functional
	- formalized conduit "proxy" runs as a single model input switch
	- fixed a bug in the plume water vapor pressure model, which occurred for low atmospheric pressures at very high altitude (thanks, Hunga)
	- changed default wind condition to "on"

v1.0 January 25, 2022
	- this is the version submitted along with the published paper 
   "Rowell, Colin R., A. Mark Jellinek, Sahand Hajimirza, and Thomas J. Aubry. “External Surface Water Influence on Explosive Eruption Dynamics, With Implications for Stratospheric Sulfur Delivery and Volcano-Climate Feedback.” Frontiers in Earth Science 10 (2022). https://www.frontiersin.org/article/10.3389/feart.2022.788294. "






