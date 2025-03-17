function [dataout] = importoxygen(fpattern)

% Scan for importable files
Folder = pwd; % Folder will be current working directory
files   = dir(fullfile(Folder, fpattern)); % scan for all files containing ident in name
if size(files,1)==0;error('Nothing to import...');end % error in case no files are found

cnt=0;
for k=1:numel(files)
    if ~contains(files(k).name,{'~','$'})
        cnt=cnt+1;
        %readtable(files(1).name,'FileType','spreadsheet','ReadVariableNames','Range','A13..')
        dat=importdata(files(k).name);
        oxdat=dat.data.OxygenOrginal(:,1:26);
        dattime=datetime(dat.textdata.OxygenOrginal(14:end,1),'InputFormat','dd.MM.yy HH:mm:ss');
        data=table;
        data.time=dattime;
        varnames=dat.textdata.OxygenOrginal(13,2:26);
        varnames=erase(varnames,{'.','/','\','_'});
        for i=1:numel(varnames)
            data.(varnames{i})=oxdat(12:end,i);
        end
        alldata{k}=data;
        filelist{k}=files(k).name;
    end
end

dataout=struct;
dataout.X=alldata;
dataout.filelist=filelist;

if cnt==1
    dataout=dataout.X{1};
end

end
