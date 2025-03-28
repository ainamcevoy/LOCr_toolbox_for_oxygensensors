function [OutlierGrp,OutlierCol] = AOCoutliers2021(NOM_substrate)

%AOCoutliers2021
%   function is used in oxygen_analysis_automated_2021_LOCr_paper.m
%   Names bad sensors in oxygen_analysis_automated_2021_LOCr_paper
%   Data for bad sensors is deleted.
%   First, visually identify outliers from comparison with other sensors in
%   the same treatment.
%   Second, add a case in the list below.
%   OutlierGrp=Name the treatments that have bad sensors
%   OutlierCol=Identifies the corresponding sensor line colors
%           dark_blue=1,red=2, yellow=3,purple=4,green=5,light_blue=6 


switch NOM_substrate
    
    case 'test-dataset-DWTP-G'
        OutlierGrp=cellstr(char('G-3')); %Name the plots that have bad lines
        OutlierCol=[4]; %dk blue=1,red=2, yellow=3,purple=4,green=5,light blue=6
    case 'test-dataset-DWTP-N'
        OutlierGrp=cellstr(char('N-1','N-1', 'N-5', 'N-4', 'N-3', 'N-2', 'N-2')); %Name the plots that have bad lines
        OutlierCol=[3, 1, 2, 3, 1, 2, 4]; %dk blue=1,red=2, yellow=3,purple=4,green=5,light blue=6 
    case 'test-dataset-acetate'
        OutlierGrp=cellstr(char('100','100', '150', '75','75','50', '50', '50', '25','10')); %Name the plots that have bad lines 50 is invalid for these results as only one replicate remains
        OutlierCol=[2, 4, 4, 2, 1, 2, 1, 4, 2,1]; %dk blue=1,red=2, yellow=3,purple=4,green=5,light blue=6 case 'optx009_851'

end

end

