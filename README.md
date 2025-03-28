# LOCr_toolbox_for_oxygensensors

Title of the dataset:
Oxygen sensor measurements for sodium acetate dilutions and natural wates. 

Version of dataset <include description of update>:
version_1 

Creators (include contact):
Aina McEvoy /email:ainam@chalmers.se
Kathleen Murphy / email:murphyk@chalmers.se 

Contributors:
Tomas McKelvey/ email:tomas.mckelvey@chalmers.se


Related publication:
The dataset is supplement to the publication: New rapid methods for assessing the production and removal of labile organic carbon in water treatment using fluorescence and oxygen measurements . 
https://doi.org/10.1021/acsestwater.5c00153

Description:
Algorithm code and test data for extrakting Labile Oxygen Consumption rate (LOCr) from PreSense sensorvials raw data outputs.
Included are oxygen concnetration data for DWTP-G, DWTP-N and acetate calibration curve. 


Keywords:
labile organic carbon 
drinking water 
PARAFAC
FDOM
oxygen sensors 
biostability


Spatial coverage:
Natural water and drinking water samples collected from different sites in or near the city of Gothenburg,located on the west coast of Sweden.   

Temporal coverage:
Natural water and drinking water collected between September 2022 and August 2023.

This folder contains the following files:

CODE:
AOCoutliers2021.m
AOCsetup2021.m
oxygen_analysis_automated_2021_LOCr_paper.m

Assisiting code and folders with code for data import and figure extraction - passive no adjustment needed
importoxygen.m
"export_fig"

DATA
Data folder "test_dataset" including three folders
"acetate"
"DWTP-G"
"DWTP-N"
info_dataset.docx

All the data files described above contain the following information:

1. one .xlsx file containing the raw data form the PreSens sensor vials including time stamps
2. one autoresults.csv file, created by the algorithm and renewed everytime the algorithm for the corresponding raw data in the foder is run.
3. one output folder, created by the algorithm and renewed everytime the algorithm for the corresponding raw data in the foder is run. This folder containins figures of the different steps included in the algorithm. 

 
Methods, materials and software:
Oxygen concentrations were logged every three minutes for a minimum of three days.
PRESENSE bjknbjkjb


To run the data:
Open 
oxygen_analysis_automated_2021_LOCr_paper.m
AOCsetup2021.m
AOCoutliers2021.m

In oxygen_analysis_automated_2021_LOCr_paper.m un-comment which case to run and comment the remaining cases

Example:
"""

% specify the substrate type for this dataset. It will determine the threshold for detecting peaks

%NOM_substrate = 'test-dataset-DWTP-G';
%NOM_substrate = 'test-dataset-DWTP-N';
NOM_substrate = 'test-dataset-acetate';

"""

In AOCsetup2021.m the parameters such as start and end time, smoothing and names are declared. This has already been done for the three example test data sets.

Example:
"""
    case 'test-dataset-DWTP-G' 
        ExptNum=1;
        cd('test_dataset\DWTP-G\');
        grp={'G-1','G-1','G-1','G-1','G-2','G-2','G-2','G-2','G-3','G-3','G-3','G-3','G-4','G-4','G-4','G-4','X1','X1','X1','X1','X2','X2','X2','X2'};
        StartHour=0.25; %70h [10,95]
        EndHour=70; %89h
        Span1=50; %900 min = about 15 h. Increase this to get more smoothing
        plotlimits=[-0.0005 0.0005]; %Magnification of second plot (1st derivative)
        doscatterplot=true;
        y_firstplot=[100,400];%ylim(y_firstplot(1),y_firstplot(2))
        y_secondplot=[150,400];
        StPt=20 %for plotting derivatives
"""

In AOCoutliers2021.m the outliers are given as inputs, first which sample thehn which replicate. The current inputs are in agreement with the published data 

Example:
"""
   case 'test-dataset-DWTP-G'
        OutlierGrp=cellstr(char('G-3')); %Name the plots that have bad lines
        OutlierCol=[4]; %dk blue=1,red=2, yellow=3,purple=4,green=5,light blue=6

"""


Funding:
Swedish Research Council FORMAS (2021-01411)
Åke och Greta Lissheds Stiftelse (2021-00162)
Gustav Richters Stiftelse (2021-00668) 
Chalmers ICT Area of Advance. 

How to cite this data:
McEvoy, Aina; Paul, Catherine; Modin, Oskar; Mohammadi, Amir Saeid, Murphy, Kathleen (2025). Data for: " New rapid methods for assessing the production and removal of labile organic carbon in water treatment using fluorescence and oxygen measurements https://doi.org/10.1021/acsestwater.5c00153 

This dataset is published under the CC BY (Attribution) license.
This license allows reusers to distribute, remix, adapt, and build upon the material in any medium or format, so long as attribution is given to the creator.
