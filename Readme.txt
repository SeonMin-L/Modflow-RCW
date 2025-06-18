This repository provides a set of scripts for evaluating the intake performance of radial collector wells (RCWs) using MODFLOW 6. The model reflects the physical characteristics of RCWs and simulates their hydraulic behavior in groundwater systems.

1. Radial_Well.py
This script implements the physical characteristics of a radial collector well. It is imported and executed within the main workflow to calculate well performance.

2. Gis_FloPy_nom.py
This script generates the required MODFLOW 6 input files (e.g., dis, drn, riv, etc.) based on geospatial input data. All necessary inputs are derived from GIS shapefiles, which define the spatial configuration and boundary conditions of the model.

3. main.py
This is the main execution script that integrates both Gis_FloPy_nom.py and Radial_Well.py. It automates the entire process from model setup to the evaluation of RCW intake performance, making it a convenient entry point for running simulations.

4. model_setting_gis.zip
This folder contains a set of shapefiles that define the spatial and conceptual settings of the model. It serves as the primary input to Gis_FloPy_nom.py

5. pipe.xlsx
This Excel file contains the configuration for radial collector well laterals, including their spatial placement and hydraulic properties. Each sheet within the file corresponds to a single RCW (radial collector well). To simulate multiple RCWs, add additional sheets—one for each well.

Within each sheet, the laterals (screened horizontal pipes) are defined row by row. Each lateral is specified using the following fields:

 layerX, rowX, colX: Grid cell location number where the lateral X is located
 lengthX: Length of the lateral segment in that cell (in meters)

 The sheet also contains global parameters for the corresponding RCW, including:
 head: Target head (m) maintained inside the caisson (central shaft)
 row, col: Grid location of the caisson (all layers at this location are set as inactive)
 well_D: Diameter of the laterals (m)
 Skin_K: Hydraulic conductivity of the filter layer around each lateral (m/day)
 Skin_Thick: Thickness of the filter layer (m)
 pipe_open: Open area ratio of the pipe (dimensionless, 0–1)
 roughness: Relative roughness of the pipe (mm)
 B_coefficient: Head loss coefficient representing convergence and frictional loss near the laterals (mm)

All values in this file are automatically parsed and used by Radial_Well.py to generate the final MODFLOW 6 drn package input.

To run the simulation via main.py, the following requirements must be satisfied:
The model workflow is executed through FloPy, which creates and runs MODFLOW 6 simulations within a designated working directory.
The MODFLOW 6 executable (mf6.exe) must be located in the same working directory where the simulation is executed. FloPy calls this executable during the model run process.
You can download the mf6.exe executable from the official USGS website:
👉 https://water.usgs.gov/water-resources/software/MODFLOW-6/
