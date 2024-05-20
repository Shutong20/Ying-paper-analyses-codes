% This function is used to automatically process the vacuole properties
% within the same condition that are stored under same large folder name.
% Each vacuole information is stored in the smaller subfolder within
% the large folder, including the individual tiff files of the z stack
% as well as csv file including the center position of that vacuole.

close all
clear

%% Select the large folder of certain condition  and find out the subfolders
disp('Please select the large folder of certain condition containing subfolders of vacuoles.');
folderpath = uigetdir;
[~, Bigfolder_name] = fileparts(folderpath); % Extract the name of big folder
Bigfolder = dir(fullfile(folderpath));
subfolder = Bigfolder([Bigfolder(:).isdir]);  % Only select subfolders within the big folder
subfolder = subfolder(~ismember({subfolder(:).name},{'.','..'}));

addpath(folderpath)
cd(folderpath)

%% Create the excel file for writing the processed data
checkforfile = isfile('Vacuole information.xlsx');
if checkforfile == 0
    header = {'filename','Total number of vacuole','Vacuole volume','Vacuole surface area','X pixel of vacuole center','Y pixel of vacuole center','Z pixel of vacuole center','Shape fitting metric'};
    writecell(header,fullfile(folderpath, 'Vacuole information.xlsx'),'Range','A1');
    N = 0;
else
    N = size(readmatrix('Vacuole information.xlsx'),1);
end

% If there is no subfolder within the big folder, just use the big folder
% as folder path and define the number of subfolder as 1
if ~numel(subfolder)
    N_file = 1;
else
    N_file = numel(subfolder);
end

if numel(subfolder)
    for i = 1:numel(subfolder)
        folder = fullfile(folderpath,subfolder(i).name);
        [n2,V,SA,x,y,z,chi2] = Chan_edge_detection_algorithm_Tong_modified(folder);

        %Save the results to the xlsx file
         AA = strcat('A',num2str(N+2));
         BB = strcat('B',num2str(N+2));
         CC = strcat('C',num2str(N+2));
         writecell({subfolder(i).name},'Vacuole information.xlsx','Range',AA)
         writematrix(n2,'Vacuole information.xlsx','Range',BB)            
         writematrix([V',SA',x,y,z,chi2'],'Vacuole information.xlsx','Range',CC)
         N = size(readmatrix('Vacuole information.xlsx'),1);

    end
    
else 
    folder = fullfile(folderpath);
    [n2,V,SA,x,y,z,chi2] = Chan_edge_detection_algorithm_Tong_modified(folder);
    AA = strcat('A',num2str(N+2));
    BB = strcat('B',num2str(N+2));
    CC = strcat('C',num2str(N+2));
    writecell({Bigfolder_name},'Vacuole information.xlsx','Range',AA)
    writematrix(n2,'Vacuole information.xlsx','Range',BB)            
    writematrix([V',SA',x,y,z,chi2'],'Vacuole information.xlsx','Range',CC)
    N = size(readmatrix('Vacuole information.xlsx'),1);
end



