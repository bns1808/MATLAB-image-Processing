function[row_bottom, row_top, column_left,mm_per_pix,column_right,x_f2,x_f]=cropping_imagesubtraction(home_points,option,option2,density_range,exp_1,exp_last,Experiment_dir_path, exp_move, option3,n1, med1)

dir_path = pwd;
main_path = dir_path;

% Calibration  %%%%%%%%%

calib_dir_path = strcat(main_path,'/Calibration');

% Get the list of files in the directory
calibfile_list = dir(fullfile(calib_dir_path, '*.jpg'));

dummy2=size(density_range);
Calimages=dummy2(2);
num_images = Calimages;


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


        % Create the new directory if it doesn't already exist

    if ~exist(new_dir_path2, 'dir')
        mkdir(new_dir_path2);
    end

    for exp=exp_1:exp_last
        if ~exist("exp_dir2_path{exp}", 'dir')
            mkdir(exp_dir2_path{exp});
        end
    end



    images = cell(1, num_images);
    cd(new_dir_path)


    if option==2.5

%Read the first calibration image. This will be used to get the coordinates

    for i = 1:1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Look Over this%%%%%%%%%%
        Variable = ['IMG',num2str(i),'.jpg'];
        image=imread(Variable);
        images{i} = image;
        imtool(image)
    end


    for i = 1:8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Look Over this%%%%%%%%%%

%         Variable = ['IMG',num2str(i),'.jpg'];
%         image=imread(Variable);
%         images{i} = image;
%         imtool(image(home_points(1,2)-100:home_points(1,2)+100,home_points(1,1)-150:home_points(1,1)+300))
%                 cd(exp_dir_path{i})
    end
    end


% The following is a manual way of checking whether the tank has been moved
%     for exp=exp_1:exp_last
%     cd(exp_dir2_path{exp})
%     Variable = 'IMG1.jpg';
%     imtool(image)
%     input("continue?")
%     end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mm_per_pix = 585/(home_points(1,3)-home_points(1,1)); %mm/pix
    x_f = round(150/mm_per_pix);
    x_f2 = 3*x_f;

    row_bottom=home_points(1,2)+10;
    row_top= row_bottom-round(160/mm_per_pix); 
    column_left = home_points(1,1)+round(50/mm_per_pix);
    column_right= home_points(1,1)+round(750/mm_per_pix);



% %%%%%%%Image 2%%%
% %Points for bottom left
%     home_points(2,1:2)=[1 1];
% %Points for top right
%     home_points(2,3:4)=[1 1];
% %Points for top leftgyyyy
%     home_points(2,5:6)=[nan nan];
% 
% %%%%%%%Image 3%%%
% %Points for bottom left
%     home_points(3,1:2)=[1 1];
% %Points for top right
%     home_points(3,3:4)=[1 1];
% %Points for top left
%     home_points(3,5:6)=[nan nan];
% 
% %%%%%%%Image 4%%%
% %Points for bottom left
%     home_points(4,1:2)=[1 1];
% %Points for top right
%     home_points(4,3:4)=[1 1];
% %Points for top left
%     home_points(4,5:6)=[nan nan];
% 
% 
% %%%%%%%Image 5%%%
% %Points for bottom left
%     home_points(5,1:2)=[1 1];
% %Points for top right
%     home_points(5,3:4)=[1 1];
% %Points for top left
%     home_points(5,5:6)=[nan nan];
% 
% %%%%%%%Image 6%%%
% %Points for bottom left
%     home_points(6,1:2)=[1 1];
% %Points for top right
%     home_points(6,3:4)=[1 1];
% %Points for top left
%     home_points(6,5:6)=[nan nan];
%  
% %%%%%%%Image 7%%%
% %Points for bottom left
%     home_points(7,1:2)=[1 1];
% %Points for top right
%     home_points(7,3:4)=[1 1];
% %Points for top left
%     home_points(7,5:6)=[nan nan];
% 
% %%%%%%%Image 8%%%
% %Points for bottom left
%     home_points(8,1:2)=[1 1];
% %Points for top right
%     home_points(8,3:4)=[1 1];
% %Points for top left
%     home_points(8,5:6)=[nan nan];
% 
% %%%%%%%Image 9%%%
% %Points for bottom left
%     home_points(9,1:2)=[1 1];
% %Points for top right
%     home_points(9,3:4)=[1 1];
% %Points for top left
%     home_points(9,5:6)=[nan nan];


% %%%%%%%Image 10 - Experiment 1010 - 0 height%%%
% %Points for bottom left
%     home_points(10,1:2)=[1 1];
% %Points for top right
%     home_points(10,3:4)=[1 1];
% %Points for top left
%     home_points(10,5:6)=[nan nan];
% 
% %%%%%%%Image 11 - Experiment 1020 - 0 height%%%
% %Points for bottom left
%     home_points(11,1:2)=[1 1];
% %Points for top right
%     home_points(11,3:4)=[1 1];
% %Points for top left
%     home_points(11,5:6)=[nan nan];
% 
% %%%%%%%Image 12 - Experiment 1010 - 10 height%%%
% %Points for bottom left
%     home_points(12,1:2)=[1 1];
% %Points for top right
%     home_points(12,3:4)=[1 1];
% %Points for top left
%     home_points(12,5:6)=[nan nan];
% 
% %%%%%%%Image 13 - Experiment 1020 - 10 height%%%
% %Points for bottom left
%     home_points(13,1:2)=[1 1];
% %Points for top right
%     home_points(13,3:4)=[1 1];
% %Points for top left
%     home_points(13,5:6)=[nan nan];


% For cropping

% if move==1
% 
% end


if option==3


    %Read the first image

    Variable = ['IMG',num2str(1),'.jpg'];
    img=imread(Variable);



    %Crop it to the aforementioned dimensions
    cropped_img1 = img(row_top:row_bottom, column_left:column_right, :);

    new_file_name = ['IMG',num2str(1),'_cropped.jpg'] ;

    %Save file 1 (Calibration image 1) - this will not be black

    imwrite(cropped_img1, fullfile(new_dir_path2, new_file_name));

    % Loop through each image for calibration
    cropped_img1_og = cropped_img1;
    cd(new_dir_path)
    for i = 1:Calimages
        % Read the original image
        Variable = ['IMG',num2str(i),'.jpg'];
        img=imread(Variable);
    
        % Crop the image to the selected ROI
        cropped_img = img(row_top:row_bottom, column_left:column_right, :);
        subtracted_img = imsubtract(cropped_img1,cropped_img);
    
        % Save the cropped image to the new directory with a new file name
        new_file_name = ['IMG',num2str(i),'_cropped.jpg'] ;

        if option3==2
                Size=size(subtracted_img);
                x1= round(Size(1)/2);
                y1 = round(Size(2)/2);
                subtracted_img = imresize(subtracted_img, [x1 y1]);
                subtracted_img = imresize(subtracted_img, [Size(1),Size(2)]);
        end

        imwrite(subtracted_img, fullfile(new_dir_path2, new_file_name));
    end

    if option2==2 || option2==3

    for exp=exp_1:exp_last

    cd(exp_dir_path{exp})
    Variable = ['IMG',num2str(1),'.jpg'];
    img=imread(Variable);
    if option3==2
        cropped_img1 = img(row_top+exp_move(exp,2):row_bottom+exp_move(exp,2), column_left+exp_move(exp,1):column_right+exp_move(exp,1), :);
    else
        cropped_img1 = img(row_top:row_bottom, column_left:column_right, :);
    end


    for i = 1:length(file_list{exp})
        % Read the original image
        Variable = ['IMG',num2str(i),'.jpg'];
        img=imread(Variable);

        if option3==2
            cropped_img = img(row_top+exp_move(exp,2):row_bottom+exp_move(exp,2), column_left+exp_move(exp,1):column_right+exp_move(exp,1), :);
            subtracted_img = imsubtract(cropped_img1,cropped_img);
         cd(main_path)
        [subtracted_img] = imageprocessingtool(subtracted_img,n1,med1); 
        cd(exp_dir_path{exp})
        else
        % Crop the image to the selected ROI
        cropped_img = img(row_top:row_bottom, column_left:column_right, :);
        subtracted_img = imsubtract(cropped_img1,cropped_img);
         cd(main_path)
        [subtracted_img] = imageprocessingtool(subtracted_img,n1,med1); 
        cd(exp_dir_path{exp})
        end

        if option3==2
                Size=size(subtracted_img);
                x1= round(Size(1)/2);
                y1 = round(Size(2)/2);
                subtracted_img = imresize(subtracted_img, [x1 y1]);
                subtracted_img = imresize(subtracted_img, [Size(1),Size(2)]);
        end
    
        % Save the cropped image to the new directory with a new file name
        new_file_name = ['IMG',num2str(i),'_cropped.jpg'] ;
        imwrite(subtracted_img, fullfile(exp_dir2_path{exp}, new_file_name));
    end
    subtracted_img = imsubtract(cropped_img1_og,cropped_img1);
    new_file_name = ['IMG',num2str(1),'_cropped.jpg'] ;
    imwrite(subtracted_img, fullfile(main_path, new_file_name));

    end
    
    end

end
