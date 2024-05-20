% This script contains two subfunctions: Call_directory & Data_extraction
% Call_directory is to navigator through all raw data folders until finding
% the final .xyz file, within which subufunction Data_extraction is called.
% Data_extract is to extract individual GEM trajectories from .xyz file and
% save as a structure variable, which is then used as the input for the
% "Velocity_correlation_MD_simulation.m"

Abs_input = "/Users/shut01/Documents/Levy simulation folder/Glen's simulation results/Glen's simulation results 4 including polysomes and RNA/Raw simulation data";
% Abs_input = '/Volumes/TONG SHU 3/Temp - Ying/Glen simulations/run_ribosome_polysome_v1_npoly6_addrnaEps1/vfr0.30/fip0.25/csRNA5.0/gemd125/s1';
Call_directory(Abs_input)
function Call_directory(input)
    cd(input)
    FolderInfo = dir(input);
    Name = extractfield(FolderInfo,'name');
    Name_type = extractfield(FolderInfo,'isdir');
    Name_type = cell2mat(Name_type); % Transform cell array to logical matrix
    temp = startsWith(Name,'.'); % Get rid of directories with name started with '.'
    Name(temp) = [];
    Name_type(temp) = [];
    
    if ~any(Name_type) % If all file name type is not directory, which is the end level folder
        disp(input);
        Name_xyz_index = ~cellfun(@isempty,regexp(Name,'.xyz$')); % Find file index with name ending with .xyz
        Data_extraction(Name{Name_xyz_index})
        return        
    else
        Name = Name(Name_type); % Only extract the folders
        for i = 1:length(Name)
            Name_i = Name{i};
            update_input = append(input,'/',Name_i);
            Call_directory(update_input)
        end
    end
end


function Data_extraction(filename)
    data = readmatrix(filename,'FileType','text');

    % Extract basic information:
    Nan_column = find(isnan(data(:,1)));
    N_Gem = data(Nan_column(2)-1,1); % Number of GEMs in this simulation
    t_MD = data(Nan_column,3);
    if t_MD(1) == 0
        t = t_MD/t_MD(2); % Time frames run in this simulation
    else
        error('MD simulation time does not starts with 0th frame');
    end

    % Save the position information into a data structure
    GEM_traj = cell(N_Gem,1); % Initiate GEM_traj cell array, with each element being xyz trajectory of one GEM
    for GEM_traj_indi = 1:N_Gem
        GEM_traj{GEM_traj_indi} = zeros(length(t), 3);
    end

    for i = 1:length(Nan_column)
        start_row = Nan_column(i)+1;
        GEM_traj_frame = data(start_row:start_row+N_Gem-1,:);
        for j = 1:N_Gem
            GEM_traj{j}(i,:) = GEM_traj_frame(j,:);
        end
    end

    [VC_norm, VC_sem_norm] = Velocity_correlation_MD_simulation(GEM_traj);
    cd("/Users/shut01/Documents/Levy simulation folder/Glen's simulation results/Glen's simulation results 4 including polysomes and RNA/Analyses_Combined")
%     filename_main = extractBetween(filename,"run_",".xyz");
    filename_main = extractBetween(filename,"gel_",".xyz");
    saveas(gcf,[filename_main{:},'.fig'])
    save([filename_main{:},'.mat'],'VC_norm', 'VC_sem_norm')
end
