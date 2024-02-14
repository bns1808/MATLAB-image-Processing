function [] = image_name(option2,density_range,exp_1,exp_last,Experiment_dir_path);

dir_path = pwd;
main_path = dir_path;

% Calibration  %%%%%%%%%

calib_dir_path = strcat(main_path,'/Calibration');
% Get the list of files in the directory
calibfile_list = dir(fullfile(calib_dir_path, '*.jpg'));

dummy2=size(density_range);
Calimages=dummy2(2);
num_images = Calimages;





if option2==2 || option2==3

file_list{1} = dir(fullfile(Experiment_dir_path{1}, '*.jpg'));
file_list{2} = dir(fullfile(Experiment_dir_path{2}, '*.jpg'));
file_list{3} = dir(fullfile(Experiment_dir_path{3}, '*.jpg'));
file_list{4} = dir(fullfile(Experiment_dir_path{4}, '*.jpg'));

end

    %For Calibration
if option2==1 || option2==3

    for i = 1:length(calibfile_list)

        % Define the old and new file names
        old_file_name = fullfile(calib_dir_path, calibfile_list(i).name);
        new_file_name = fullfile(calib_dir_path, sprintf('IMG%d.jpg', i));
    
        % Rename the file
        movefile(old_file_name, new_file_name);
    end
end

    %For Experiment

    if option2==2 || option2==3

    %Experiment1 to 4
for exp=exp_1:exp_last
    for i = 1:length(file_list{exp})
        old_file_name = fullfile(Experiment_dir_path{exp}, file_list{exp}(i).name);
        new_file_name = fullfile(Experiment_dir_path{exp}, sprintf('IMG%d.jpg', i));

        % Rename the file
        movefile(old_file_name, new_file_name);

    end
end
    end

end
