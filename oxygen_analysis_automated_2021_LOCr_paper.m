%% oxygen_analysis_automated_2021
%current directory for code and export figure
parentDirectory = fileparts(pwd); 
addpath export_fig\


%Import data
clear all
close all
replacefigs=true;
%replacefigs=false;
autowarning=false; %warning if an existing autoresults file will be overwritten
doscatterplot=true; %default will be changed to true if there are calibration data (acetate standards)



% specify the substrate type for this dataset. It will determine the threshold for detecting peaks

%NOM_substrate = 'test-dataset-DWTP-G';
NOM_substrate = 'test-dataset-DWTP-N';
%NOM_substrate = 'test-dataset-acetate';

%% Set up the basic parameters for this dataset
% these are specified in the function AOCsetup2021.m

SetupObject=AOCsetup2021(NOM_substrate);

ExptNum=SetupObject.ExptNum;
grp=SetupObject.grp;
StartHour=SetupObject.StartHour;
EndHour=SetupObject.EndHour;
Span1=SetupObject.Span1;
plotlimits=SetupObject.plotlimits;
StPt=SetupObject.StPt;
y_firstplot=SetupObject.y_firstplot;
y_secondplot=SetupObject.y_secondplot;
doscatterplot=SetupObject.doscatterplot;


%StPt=50 %for plotting derivatives

%% import oxygen dataset
addpath ..\..\
x=importoxygen('*.xlsx');
xbackup=x;
%x=xbackup;



%% save data?
%save Data_expt2021_1.mat x

%% Plot - overview of all data
% Names of the wells on the plate
wellnames=x.Properties.VariableNames(3:end);

% Variables for assigning a color per group
grpu=unique(grp);
col=parula(numel(grpu)); %col=colours

% plot of raw data
figure
hold on % superimpose all the lines
for n=1:numel(wellnames)
    c=col(strcmp(grpu,grp{n}),:); % Assign a color based on the group
    plot(x.TimeMin/60,table2array(x(:,n+2)),'Color',c,'LineWidth',2);
    name=[grp{n} num2str(n)];
    %pause
end
legend1=legend(grp,'Location','eastoutside');
xlabel('time h.')
ylabel('O_2 (100%sat=258uM)')
title(NOM_substrate)
if replacefigs
   fig4export=['output/' NOM_substrate '-raw']
   export_fig(fig4export, '-tif', '-transparent', '-r300');
   savefig(fig4export);
end


%% plot of raw and normalised data
figure

maxvals=[];imaxvals=[];  %maxvals is the highest oxygen value measured by a sensor
nTreats=size(grpu,2);
initialperiod=StartHour*60; %minutes at beginning to ignore
for n=1:numel(grpu)
    idx=find(strcmp(grp,grpu{n}));
    repdata=table2array(x(:,idx+2));
    [maxvals{n} imaxvals{n}]=max(repdata(x.TimeMin<1500,:));
    maxtime=x.TimeMin(end);
    
    subplot(3,nTreats,n)
    plot(x.TimeMin,repdata,'LineWidth',2);
    xlim([initialperiod maxtime])
    title(grpu{n})
    
    subplot(3,nTreats,n+nTreats)
    idx=find(strcmp(grp,grpu{n}));
    plot(x.TimeMin,repdata./maxvals{n},'LineWidth',2);
    title(grpu{n})
    xlim([initialperiod maxtime])
    
end
sgtitle('1. Raw, 2. scaled, 3. scaled + outliers removed')

%% Identify Outliers and Replace with NaN
%specify in next two lines which plot group and which line to remove
%x=xbackup;
[OutlierGrp,OutlierCol] = AOCoutliers2021(NOM_substrate);

Outliertext=' ';
xcor=x;
for i=1:numel(OutlierCol)
    idx4=find(strcmp(grp,OutlierGrp{i}));
    ind_bad=idx4(OutlierCol(i));
    xcor(:,ind_bad+2)=array2table(NaN*ones(size(xcor,1),1));
    Outliertext=[Outliertext OutlierGrp{i} ', '];
end

%Final Row of Plots, with outliers removed
for n=1:numel(grpu)
    idx=find(strcmp(grp,grpu{n}));
    repdatac=table2array(xcor(:,idx+2));
    maxvalsc=max(repdatac(xcor.TimeMin<1500,:));
    
    subplot(3,nTreats,n+2*nTreats)
    idx=find(strcmp(grp,grpu{n}));
    plot(xcor.TimeMin,repdatac./maxvalsc,'LineWidth',2);
    title(grpu{n})
    xlim([initialperiod maxtime])
    %ylim(ylimits1(1),ylimits1(2))
    sgtitle('1. Raw, 2. scaled, 3. scaled + outliers removed')
end
peakindex=cell2mat(imaxvals); %index corresponding to maximum oxygen reading for each sensor
if replacefigs
   %savefig('raw_scaled');
   export_fig('output/raw_scaled', '-tif', '-transparent', '-r300');
end
%% AUTOMATIC DETECTION OF PEAKS

%choose data for the calculations, either 'xcor' or 'all'
%data2plot='all';
data2plot='xcor';

switch data2plot
    case 'xcor'
    xin=xcor; %use outlier-free data
    case 'all'
    xin=x; %use all data
end

%plot derivatives and automatically identify boundaries
%close all
dataout=NaN*ones(numel(wellnames),7);
TailStart=50*60; %time in minutes
[~,EndIndex]=min(abs(x.TimeMin-EndHour*60));
smoothing_time=x.TimeMin(Span1);

%plots %close all
figure
nRep=[];
%for i=1:length(grpu),nRep(i)=sum(strcmp(grpu(i),grp)),end  %number of reps in each unique group, size grpu
for n=1:length(grp),nRep(n)=sum(strcmp(grp(n),grp)),end %number of reps in each unique group matching each vial (24 length)

j=0;
for n=1:numel(wellnames)
    
    j=j+1;  %replicate vials nested within figure
    
    datain=table2array(xin(:,n+2));
    StartIndex=peakindex(n);  %index where the maximum sensor oxygen reading occurred initially
    
    if sum(isnan(datain))<numel(datain) %data are not only NaNs
        %calculate first derivative to see where the curve changes shape
        diferential=gradient(smooth(datain,Span1,'sgolay'));
        diffdat=smooth(diferential,Span1,'sgolay');
        FirstDifferential=diffdat(StartIndex:EndIndex);
        TimeFromStartIndex=xin.TimeMin(StartIndex:EndIndex);
        time2plot=TimeFromStartIndex/60;
        
        %calculate baseline, which is the average derivative of the tail
        %(after 50 h).
        Baseline=mean(FirstDifferential(TimeFromStartIndex>TailStart));
        errordata=(FirstDifferential-Baseline);
        errordata_neg=-errordata;
        
        
        
        %locate the  biggest negative peak in the derivative plots
        [PKSneg,LOCSneg]= findpeaks(errordata_neg(TimeFromStartIndex>initialperiod),find(TimeFromStartIndex>initialperiod,1,'first'):numel(time2plot),'SortStr','descend');%positive peaks in diffdat
        neglocs=LOCSneg(PKSneg>0);
        
        %find the zero value immediately before the big negative peak
        before_data=errordata_neg(1:neglocs);
        ilast_before=find(before_data<0,1,'last');
        if isempty(ilast_before)
            ilast_before =1;
        end
        
        %find the zero value immediately after the big negative peak
        after_data=errordata(neglocs:end);
        ifirst_after=find(after_data>0,1,'first');
        if isempty(ifirst_after)
            ilast_before =1;
        end
        
        %The located points are called Ipeak and are plotted with x o x
        Ipeak=[];
        try Ipeak = [StartIndex+ilast_before StartIndex+neglocs(1) StartIndex+neglocs(1)+ifirst_after]; end
        if  isempty(Ipeak)
            Ipeak=[1,6000]; %give them 10 h
            FirstDifferential=zeros(length(time2plot),1);
        end
try
        dataout(n,:)=[n,Ipeak(1),xin.TimeMin(Ipeak(1))/60,datain(Ipeak(1)),Ipeak(end), xin.TimeMin(Ipeak(3))/60,datain(Ipeak(3))];
catch
end
        
        %plots - initialise colors
        c=col(strcmp(grpu,grp{n}),:); % Assign a color based on the group
        
        %plots - top row
        %start plot at StPt- tells how much data to show in derivative plots,but doesnt affect calculations
        subplot(2,nRep(n),j)
        [ax,h1,h2] = plotyy(xin.TimeMin(StPt:end),datain(StPt:end), xin.TimeMin(StPt:end), diffdat(StPt:end));
        hold on, plot(xin.TimeMin(Ipeak(1)),datain(Ipeak(1)),'kx','MarkerSize',10);
try        hold on, plot(xin.TimeMin(Ipeak(2)),datain(Ipeak(2)),'ko','MarkerSize',3); end
try        hold on, plot(xin.TimeMin(Ipeak(3)),datain(Ipeak(3)),'kx','MarkerSize',10); end
        if j==1
    
            ylabel('oxygen mM')
        end
        title(grp{n})

        subplot(2,nRep(n),j+nRep(n))
        plot(time2plot,errordata);
       
        hold on, plot([min(time2plot),max(time2plot)],[0,0])
        hold on, plot(xin.TimeMin(Ipeak(1))/60,0,'kx','MarkerSize',10);
try        hold on, plot(xin.TimeMin(Ipeak(2))/60,0,'ko','MarkerSize',3); end
try        hold on, plot(xin.TimeMin(Ipeak(3))/60,0,'kx','MarkerSize',10); end
        xlabel('time h.')
        ylim(plotlimits)
  
        title(grp{n})
        if j==1
            ylabel('d/dx oxygen uM')
        end
        
    end


    if mod(j,nRep(n)) == 0 && n<24  %    if mod(n,nRep(j)) == 0 && n<24
        j=0;
       if replacefigs
           fig4export=['output/deriv_' grp{n}];
           export_fig(fig4export, '-tif', '-transparent', '-r300');
          % savefig(fig4export); %%%AINA: LOOK into FIXING THIS BUG XXX
        end
        figure
        set(gcf,'position',[171  126  750  500])
    end
  
end

% % tabulate the resuls in DAT
DAT= array2table(dataout,'VariableNames',{'n' 'init_peak' 't0' 'initial' 'fin_peak' 'tf' 'final'});

%% calculate slope statistics and cumulative oxygen consumption
cumulox=(dataout(:,4)-dataout(:,7));
time2empty=(dataout(:,3)-dataout(:,6));
slope=-cumulox./time2empty;

smean=grpstats(slope',grp','mean'); %mean
sgrp=grpstats(slope',grp','gname'); %groups/treatments
sstd=grpstats(slope',grp','std');   %standard deviation 
sem=grpstats(slope',grp','sem');    %standard error of the mean
snumel=grpstats(slope',grp','numel'); %number of samples

cmean=grpstats(cumulox',grp','mean'); %mean
cgrp=grpstats(cumulox',grp','gname'); %groups/treatments
cstd=grpstats(cumulox',grp','std');   %standard deviation 
cem=grpstats(cumulox',grp','sem');    %standard error of the mean
cnumel=grpstats(cumulox',grp','numel'); %number of samples

%collate data to export to excel
cstype=cellstr(char('sample', 'smean','sstd','sem','cmean','cstd','cem','n'))';
csout=[smean,sstd,sem,cmean,cstd,cem,cnumel];
designlabels=cellstr(char('NOM_substrate','StartHour','EndHour','Span1'));
design=cellstr(char(char(NOM_substrate),char(num2str(StartHour)),char(num2str(EndHour)),char(num2str(Span1))));
DateAndTime=cellstr(datestr(now));

%write to excel


if autowarning==true
    if isfile('autoresults.xls')
         warning('about to overwrite data in the file autoresults.xls! Hit Control-C to cancel...')
         pause
     end
 end

exp_data=[smean sstd sem cmean cstd cem cnumel];
exp_data=num2cell(exp_data);
con_exp=[cgrp  exp_data]
title =cgrp;
full_exp_data=[cstype ; con_exp]

writecell(full_exp_data,'autoresults.xls')


%% plot cumulative oxygen loss (LHS) and rate (RHS)




%choose whether to plot Std Dev or Std Err:
%errormat=[cstd sstd]; %choose this to plot standard deviations
errormat=[cem sem];    %choose this to plot standard errors

%bar plot cumulative oxygen loss (LHS) and rate (RHS)
ff=[90 140 1150 380];
figure('position',ff)
subplot(1,2,1)
bar(1:nTreats,cmean)
hold on
errorbar(1:nTreats,cmean,errormat(:,1),'+','Color','k','LineWidth',2)
set(gca,'XTick',1:nTreats,'XTickLabel',cgrp,'XLim',[-0.5 (nTreats+0.5)])
ylabel({'Cumulative Oxygen Consumption',' (mg L^{-1})'})
y1=ylim;

subplot(1,2,2)
bar(1:nTreats,smean)
hold on
errorbar(1:nTreats,smean,errormat(:,2),'+','Color','k','LineWidth',2)
set(gca,'XTick',1:nTreats,'XTickLabel',sgrp,'XLim',[-0.5 (nTreats+0.5)])
ylabel({'Oxygen consumption rate',' (mg L^{-1} h^{-1})'})
if replacefigs
    fig4export='output/bar_cmean_smean';
    export_fig(fig4export, '-tif', '-transparent', '-r300');
    savefig(fig4export);
end






%% acetate experiments - scatter plots
%only works if dose is numeric, i.e. all non-numeric treatments must be given a number code.
if doscatterplot==true
    include=1:6; %only include the dilution series
    Tstat=table(sgrp,smean,sstd,sem,snumel,cmean,cgrp,cstd,cem,cnumel);
    Tstat(Tstat.snumel==0,:)=[]
    
    Tstat=Tstat(include,:)

    
    grpdose=Tstat.cgrp;
    grpdose=replace(Tstat.cgrp,{'MQ5','MQ4','MQ3','MQ2','MQ1','MQ/in+nut','MQ/in','MQ','EQU1','RÅ','kol','utg','x','gem'},{'5','4','3','2','1','6','7','8','9','10','11','12','13','14'}); %replace non-numeric values with numbers (codes) so plotting can go ahead
    
    %choose how to show error values:
    %errormat=[cstd sstd]; %choose this to plot standard deviations
    errormat=[Tstat.cem Tstat.sem];%choose this to plot standard errors
    
    %produce a scatter plot
    try
        dosenum=(str2num(char(grpdose)))
        fitline='linear'
        
        if ~isempty(dosenum)
            figure
            set(gcf,'position',[171  126  750  500])

            % plot cumlative oxygen vs dose, add linear fit avoiding zero
            subplot(1,2,1)
            plot(dosenum,Tstat.cmean,'o','Color','k','LineWidth',2,'MarkerSize',2)
            hold on
            errorbar(dosenum,Tstat.cmean,errormat(:,1),'.','Color','k','LineWidth',2)
            Poly1 = polyfit(dosenum(2:end),Tstat.cmean(2:end),1);
            yfit1 = Poly1(1)*dosenum+Poly1(2);
            pred1=(Tstat.cmean - Poly1(2))./Poly1(1);
            
            Pred=table(grpdose,pred1);
            disp(Pred)
           
            hold on;
            plot(dosenum(1:end),yfit1(1:end),'r-.');
            %ylabel({'Cumulative Oxygen Consumption',' (µmol L^{-1})'})
            ylabel({'Cumul. O_2 Consumption (µmol L^{-1})'})
            xlabel('Acetate-C (µg L^{-1})')
            
            % plot respiration rate vs dose, add power fit avoiding zero
            subplot(1,2,2)
            plot(dosenum,Tstat.smean,'o','Color','k','LineWidth',2,'MarkerSize',2)
            hold on
            errorbar(dosenum,Tstat.smean,errormat(:,2),'o','Color','k','LineWidth',2)
            p = polyfit(dosenum,Tstat.smean,1); 
            f = polyval(p,dosenum);
            plot(dosenum,Tstat.smean,'o',dosenum,f,'-.', 'Color','b')
            ylabel('O_2 Consumption rate (µmol L^{-1} h^{-1})')
            xlabel('Acetate-C (µg L^{-1})')
            

            if replacefigs
                fig4export='output/calibration';
                export_fig(fig4export, '-tif', '-transparent', '-r300');
                savefig(fig4export);
            end
        else
            fprintf('dosenum is empty')
        end
    catch
        fprintf('Attempt to make a scatter plot of dose vs conc failed')
    end
else
        fprintf('No scatter plots requested')
end


cd ..\..\
    
   