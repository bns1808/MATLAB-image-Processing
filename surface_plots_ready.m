function [] = surface_plots_ready(density_range,max_intensity,min_intensity,min_intensity2,threshold,threshold2,exp_1,exp_last,Experiment_dir_path,x_f,power,power2,mm_per_pix,t_initial)

dir_path = pwd;
main_path = dir_path;

tank_width= 76.2*(10^-3);
porosity = 0.38;

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



file_list{1} = dir(fullfile(Experiment_dir_path{1}, '*.jpg'));
file_list{2} = dir(fullfile(Experiment_dir_path{2}, '*.jpg'));
file_list{3} = dir(fullfile(Experiment_dir_path{3}, '*.jpg'));
file_list{4} = dir(fullfile(Experiment_dir_path{4}, '*.jpg'));


for exp=exp_1:exp_last

    cd(exp_dir3_path{exp})
    load('highestpoints.mat')

     saved_myans= (1/3)+saved_myanswer{exp}/x_f;

    Variable = ['IMG',num2str(1),'_cropped.jpg'];
    cd(exp_dir2_path{exp})
    img=imread(Variable);
    imgstodisect = im2double(img);
    Size=size(img);


    row = Size(1);
    column = Size(2);

    B(1:length(saved_myans))=0;
    D(1:length(saved_myans))=0;

     for i=5:length(saved_myans)

         fakeimage= fakeimages{i};

    % Remove NaN values
    X = 1:saved_myanswer{exp}(i,1);
    Y = row-saved_highestpoints_bulk{exp}(i,1:saved_myanswer{exp}(i,1));
    validData = ~isnan(Y);
    X = X(validData);
    Y = Y(validData);

    W = polyfit(X,Y,power);

    

     for k=1:saved_myanswer{exp}(i,1)
        myarray = fakeimage(1:row,k);
        GGGG=row-round(polyval(W,k));
        if round(row-polyval(W,k))<=0
            continue
        end
        myarray=fillmissing(myarray,'linear');

        for j=1:round(polyval(W,k))
            B(i)=myarray(row+1-j)+B(i);
        end

     end


    clear X Y
    X = 1:saved_myanswer{exp}(i,2);

    Y = row-saved_highestpoints_dispersed{exp}(i,1:min(saved_myanswer{exp}(i,2),column));
    validData = ~isnan(Y);
    X = X(validData);
    Y = Y(validData);
    W = polyfit(X,Y,power);
    clear X

     for k=1:min(saved_myanswer{exp}(i,2),column-1)
        myarray = fakeimage(1:row,k);
        GGGG2=row-round(polyval(W,k));
        
        myarray=fillmissing(myarray,'linear');
        if round(row-polyval(W,k))<=0
            continue
        end
        for j=1:round(polyval(W,k))
            D(i)=myarray(row+1-j)+D(i);
        end
     end
    D(i)=D(i)-B(i);
    Percent(i)=D(i)/(D(i)+B(i));
     end

D_dummy=D;
clear D
D(1:length(D_dummy)-t_initial(exp)+1)=D_dummy(t_initial(exp):end);

t=find(D>0,1,"last");
D=(D*(mm_per_pix^2)*tank_width*porosity)/(t*30);

B_dummy=B;
clear B
B(1:length(B_dummy)-t_initial(exp)+1)=B_dummy(t_initial(exp):end);

B=(B*(mm_per_pix^2)*tank_width*porosity)/(t*30);

     cd(Experiment_dir_path{exp})
     save('Buoyancy.mat','B','D','Percent','saved_myans')
     clear B D Percent saved_myans

end

end

