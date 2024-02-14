function [idx,idx2,idx3,idx4] = xfdetermine (home_points,option,option2,density_range,row_bottom, row_top, column_left,mm_per_pix,column_right,x_f2,x_f,option3,exp_1,exp_last,Experiment_dir_path)


dir_path = pwd;
main_path = dir_path;

% Calibration  %%%%%%%%%

calib_dir_path = strcat(main_path,'/Calibration');
% Get the list of files in the directory
calibfile_list = dir(fullfile(calib_dir_path, '*.jpg'));

dummy2=size(density_range);
Calimages=dummy2(2);
num_images = Calimages;


% Experiment  %%%%%%%%%

new_dir_path = strcat(calib_dir_path,'/new_directory');
% Define the new directory for the cropped images
new_dir_path2 = strcat(new_dir_path,'/new_directory');

exp_dir_path{1} =strcat(Experiment_dir_path{1},'/new_directory');
exp_dir_path{2} =strcat(Experiment_dir_path{2},'/new_directory');
exp_dir_path{3} =strcat(Experiment_dir_path{3},'/new_directory');
exp_dir_path{4} =strcat(Experiment_dir_path{4},'/new_directory');

exp_dir2_path{1} =strcat(exp_dir_path{1},'/new_directory');
exp_dir2_path{2} =strcat(exp_dir_path{2},'/new_directory');
exp_dir2_path{3} =strcat(exp_dir_path{3},'/new_directory');
exp_dir2_path{4} =strcat(exp_dir_path{4},'/new_directory');

exp_dir3_path{1} =strcat(exp_dir2_path{1},'/new_directory');
exp_dir3_path{2} =strcat(exp_dir2_path{2},'/new_directory');
exp_dir3_path{3} =strcat(exp_dir2_path{3},'/new_directory');
exp_dir3_path{4} =strcat(exp_dir2_path{4},'/new_directory');

exp_dir4_path{1} =strcat(exp_dir2_path{1},'/new_directory2');
exp_dir4_path{2} =strcat(exp_dir2_path{2},'/new_directory2');
exp_dir4_path{3} =strcat(exp_dir2_path{3},'/new_directory2');
exp_dir4_path{4} =strcat(exp_dir2_path{4},'/new_directory2');

if option2==2 || option2==3

file_list{1} = dir(fullfile(Experiment_dir_path{1}, '*.jpg'));
file_list{2} = dir(fullfile(Experiment_dir_path{2}, '*.jpg'));
file_list{3} = dir(fullfile(Experiment_dir_path{3}, '*.jpg'));
file_list{4} = dir(fullfile(Experiment_dir_path{4}, '*.jpg'));

end

if option ~= 5 
    idx=0;
    idx2=0;
    idx3=0;
    idx4=0;
end

for exp=exp_1:exp_last
 %Experiment 1
    cd(exp_dir3_path{exp})
    load("highestpoints.mat")

    saved_myans= (1/3)+saved_myanswer{exp}/x_f;
    figure(exp)
    plot(10:length(saved_myans),saved_myans(10:end,1),'-o')
    hold on
    plot(10:length(saved_myans),saved_myans(10:end,2),'-o')
    hold on 
    three(10:length(saved_myans))=3;
    plot(10:length(saved_myans),three(10:end))

    output_image = ['Nose_positions_exp', num2str(exp),'.jpg'];
    cd(dir_path)
    saveas(gcf,output_image)

    hold off
    clear three
    close all
    cd(exp_dir3_path{exp})
    % Target value
    target = 3;

    % Compute the absolute differences
    diff = abs(saved_myans(:,2) - target);

    % Find the minimum difference
    [~, idx(exp)] = min(diff);

    % Get the closest value
    closest_value = saved_myanswer{exp}(idx(exp),2);


    disp(['For experiment ' num2str(exp) ' use ' num2str(idx(exp)) '.']);

    figure(exp+4)
    cd(exp_dir2_path{exp})
    Variable = ['IMG',num2str(idx(exp)),'_cropped.jpg'];
    imgstodisect=imread(Variable);
    Size=size(imgstodisect);
        rows = Size(1);
        columns =Size(2);


    imshow(imgstodisect)
        hold on
        for col =1:saved_myanswer{exp}(idx(exp),1)
            row=saved_highestpoints_bulk{exp}(idx(exp),col);
            scatter(col, row, 25, 'red', 'filled');
        end
        hold on
        for col =1:saved_myanswer{exp}(idx(exp),2)
            row=saved_highestpoints_dispersed{exp}(idx(exp),col);
            scatter(col, row, 25, 'blue', 'filled');
        end
        y=round(target*x_f-(1/3)*x_f);
        hold on
        scatter(y,1:rows,25,'blue','filled')
        output_image = ['Nose position_Experiment ', num2str(exp),'.jpg'];
        cd(main_path)
        saveas(gcf, output_image);

        
        hold off

clear saved_myanswer

end

cd(main_path)
save('idx.mat', 'idx')

end
