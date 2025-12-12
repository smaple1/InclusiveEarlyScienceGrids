# InclusiveEarlyScienceGrids

Analysis code, analysis note, and plotting scripts for the creation of cross section value and uncertainty grids from early science pseudodata.

The scripts required to run the analysis and create the grids are contained within the "Run" directory. The processing of the campaign output files is done by the `Podio_ExtractYield_*.C` scripts, which should be modified such that the paths contained in "filenames" each point to a text file containing a list of files to be processed. An output file containing several ROOT histograms is produced, and should be moved to the "Yields" directory.

The project is organised into the following directories:
- **Run/**: Contains the scripts required to run the analysis and create the grids.
- **Yields/**: Output directory for ROOT histograms generated during the analysis.
- **Tables/**: Contains CSV files with output from analysis notebooks.
- **Plots/**: Output directory for some of the plots generated during the analysis.
- **FinalGrids/**: Contains the final grids after the analysis is complete.
- **PlottingScripts/**: Contains ROOT macros for plotting.

The remaining analysis is done in three jupyter notebooks. You will need PyROOT, LHAPDF, and yadism available.

The notebooks should be ran in the order

1. EvalYields.ipynb
2. EvalXsecPDF.ipynb
3. MakeErrorTables.ipynb

Before running `EvalYields.ipynb`, the settings can be set in the "Settings" cell, where you can set the beam energies, binning scheme (either 5 bins per decade in x and Q^2, or a HERA inspired binning are currently available), and integrated luminosity to scale the events to. Run each cell in the notebook. The output is written to a csv file in the "Tables" directory.

The settings also need to be specified before running `EvalXsecPDF.ipynb` as previously, to read the correct table. A theory card is specified for the cross section predictions - this isn't used unless the compute_predictions() function is ran. By default, `compute_predictions()` is commented out, and the precomputed grids in the PineAPPL_grids directory are used. Run through each cell in the notebook. The output is written to the "Tables" directory, and some plots are written to the "Plots" directory.

The settings need to be specified once more before running `MakeErrorTables.ipynb`, to read the correct tables. The systematic estimates to use are set in the cell below settings. Run each cell (or stop where the `degrade_y_keep_column()` function is defined. The output is written to the FinalGrids directory.

The ROOT macros in the PlottingScripts directory should be modified to point to the desired table in FinalGrids, to plot the intended beam energy configuration, binning scheme, and optimistic/pessimistic systematics.