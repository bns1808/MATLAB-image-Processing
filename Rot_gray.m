function []= Rot_gray(option2,density_range,angle,exp_1,exp_last,Experiment_dir_path);


dir_path = pwd;
main_path = dir_path;

% Calibration  %%%%%%%%%

calib_dir_path = strcat(main_path,'/Calibration');
% Get the list of files in the directory
calibfile_list = dir(fullfile(calib_dir_path, '*.jpg'));

dummy2=size(density_range);
Calimages=dummy2(2);



new_dir_path = strcat(calib_dir_path,'/new_directory');
% Define the new directory for the cropped images
new_dir_path2 = strcat(new_dir_path,'/new_directory');

exp_dir_path{1} =strcat(Experiment_dir_path{1},'/new_directory');
exp_dir_path{2} =strcat(Experiment_dir_path{2},'/new_directory');
exp_dir_path{3} =strcat(Experiment_dir_path{3},'/new_directory');
exp_dir_path{4} =strcat(Experiment_dir_path{4},'/new_directory');


if option2==2 || option2==3

file_list{1} = dir(fullfile(Experiment_dir_path{1}, '*.jpg'));
file_list{2} = dir(fullfile(Experiment_dir_path{2}, '*.jpg'));
file_list{3} = dir(fullfile(Experiment_dir_path{3}, '*.jpg'));
file_list{4} = dir(fullfile(Experiment_dir_path{4}, '*.jpg'));

end


    if ~exist(new_dir_path, 'dir')
        mkdir(new_dir_path);
    end
    for exp=exp_1:exp_last
        if ~exist("exp_dir_path{exp}", 'dir')
            mkdir(exp_dir_path{exp});
        end
    end


 if option2==1 || option2==3
    
    %For calibration

    for i = 1:Calimages

        % Read the original image
        img = imread(fullfile(calib_dir_path, calibfile_list(i).name));

        % Convert the image to grayscale
        gray_img = rgb2gray(img);

        % Apply rotation (if needed)
        gray_img =imrotate(gray_img,angle);

        % Save the grayscale image to the new directory with the same file name
        imwrite(gray_img, fullfile(new_dir_path, calibfile_list(i).name));

    end

  end

    %For experimental images



    if option2==2 || option2==3

    for exp=exp_1:exp_last

    for i = 1:length(file_list{exp})
        % Read the original image
        img = imread(fullfile(Experiment_dir_path{exp}, file_list{exp}(i).name));

        % Convert the image to grayscale
        gray_img = rgb2gray(img);

        % Apply rotation (if needed)
        gray_img =imrotate(gray_img,angle);

        % Convert image to double
        gray_img =im2double(gray_img);


        % Save the grayscale image to the new directory with the same file name
        imwrite(gray_img, fullfile(exp_dir_path{exp}, file_list{exp}(i).name));

    end
    end



    end

end
