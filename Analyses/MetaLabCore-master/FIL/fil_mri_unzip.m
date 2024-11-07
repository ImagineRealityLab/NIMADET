function fil_mri_unzip(filename,subject)
% function fil_mri_unzip(filename,subject)
%
% Unzips the data from Charm, where filename is scan ID and subject is the
% subject number in the study (starting at 1).
%
% Note that this function only works on FIL computer as that is the only
% place where I store the zipped data (saving space).
%
% Steve Fleming & Dan Bang, FIL, 07/06/2016

% current directy
cwd = pwd;

% filesep
fs = filesep;

% add import tool
addpath(['D:',fs,'MATLAB',fs,'spm12', fs,'Import_Archive']);

% unzip data in this folder
zipped   = ['D:',fs,'Desktop',fs,'INFABS',fs,'MRI',fs,'data_zipped',fs,filename];

% place unzipped data in this folder
unzipped = ['D:',fs,'Desktop',fs,'INFABS',fs,'MRI', subject];

% error if exist
if exist(unzipped,'dir'); error('--already unzipped'); else; mkdir(unzipped); end;

% error if not at the FIL
[~,name] = system('hostname');
if ~strcmp(name(1:end-1),'motorhead'); error('--only works at the FIL'); end;

% get folder names for unzipping
cd(zipped);
dirinfo = dir;

% loop through
for i = 3:size(dirinfo)
   Import_Archive(dirinfo(i).name,unzipped);
end

cd(cwd);

end