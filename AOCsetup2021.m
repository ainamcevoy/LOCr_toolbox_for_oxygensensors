function SetupObject= AOCsetup2021(NOM_substrate,varargin)
%AOCsetup2021
%   function is used in oxygen_analysis_automated_2021_LOCr_paper.m
%   Specifies default options for plotting, smoothing, and analysis
%   Algorithms affected by:
%        StartHour %data limit (min)
%        EndHour %data limit (max)
%        Span1 %smoothing, bigger number increases smoothing
%   Plotting is affected by:
%        StPt %minimum xaxis value (hours) when plotting derivatives
%        y_firstplot %plot range y axis
%        y_secondplot %plot range y axis
%        doscatterplot %T/F depending on if acetate series is present

%% path to BOX


%defaultsx
ExptNum=[];
StartHour=0.5; %ignore measurements at times<StartHour
EndHour=[];%ignore measurements at times>EndHour
Span1=90;  %increase this to get more smoothing
StPt=5; %minimum xaxis value (hours) when plotting derivatives
y_firstplot=[];
y_secondplot=[];
doscatterplot=false;

if nargin>1
opts=varargin{1};
error('havent implemented options yet')
end

disp(NOM_substrate)

switch NOM_substrate
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
    
    case 'test-dataset-DWTP-N' 
        ExptNum=1;
        cd(['test_dataset\DWTP-N\']);
        grp={'N-1','N-1','N-1','N-1','N-2','N-2','N-2','N-2','N-3','N-3','N-3','N-3','N-4','N-4','N-4','N-4','N-5','N-5','N-5','N-5','N-6','N-6','N-6','N-6'};
        StartHour=0.5; %70h [10,95]
        EndHour=90; %89h
        Span1=90; %900 min = about 15 h. Increase this to get more smoothing
        plotlimits=[-0.0005 0.0005]; %Magnification of second plot (1st derivative)
        doscatterplot=true;
        y_firstplot=[100,400];%ylim(y_firstplot(1),y_firstplot(2))
        y_secondplot=[150,400];
        StPt=20 %for plotting derivatives   

    case 'test-dataset-acetate' 
        ExptNum=2;
        cd(['test_dataset\acetate\']);
        grp={'10','10','10','10','25','25','25','25','50','50','50','50','75','75','75','75','100','100','100','100','150','150','150','150'}
        StartHour=0.5; %70h [10,95]
        EndHour=60; %89h
        Span1=90; %900 min = about 15 h. Increase this to get more smoothing
        plotlimits=[-0.0005 0.0005]; %Magnification of second plot (1st derivative)
        doscatterplot=true;
        y_firstplot=[100,400];%ylim(y_firstplot(1),y_firstplot(2))
        y_secondplot=[150,400];
        StPt=20 %for plotting derivatives   
    
end
SetupObject.ExptNum=ExptNum;
SetupObject.grp=grp;
SetupObject.StartHour=StartHour;
SetupObject.EndHour=EndHour;
SetupObject.Span1=Span1;
SetupObject.plotlimits=plotlimits;
SetupObject.StPt=StPt;
SetupObject.y_firstplot=y_firstplot;
SetupObject.y_secondplot=y_secondplot;
SetupObject.doscatterplot=doscatterplot;
mkdir output/

end



