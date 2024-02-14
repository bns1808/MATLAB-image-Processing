%============Image Processing: Main function================


% Option 1 = Naming Files
% Option 2 = Gray Scale and Rotation
% Option 2.5 = Determine Croping Parameters
% Option 3 = Cropping
% Option 4 = Calibration Curves
% Option 5 = Experiments
% Option 6 = MP4 Video
%%% Option 7 = Trap X - Do Not Use
%%% Option 8 = Trap X Video- Do Not Use
% Option 9 = Determine the 3x_f image
% Option 12 = Window Check
% Option 13 = Buoyancy calc setup
% Option 14 = Area Calc setup 


option=2.5;

if option==1 || option==2 || option ==2.5 || option==4 || option==5 || option == 7 ||option==9
    clc
    clearvars -except option
    close all
end

%Only work on calibration data-> option2=1, Only exp. data -> option2=2,
%Both->option2=3,
option2=3;

%Make it 2 for moving pixels, any other number and it does not move
option3=3;

set(0,'defaultTextInterpreter','latex'); %trying to set the default

l_wid=3;
f_size=25;
scatter_size=100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%User defined%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the range of densities to use for the calibration curve
density_range = [0.0, 0.04*(2/10), 0.04*(4/10),0.04*(5/10), 0.04*(6/10),0.04*(8/10),0.04, 0.05];

%Points for bottom left
    home_points(1,1:2)=[482 2244];
%Points for top right
    home_points(1,3:4)=[3669 444];
%Points for top left
    home_points(1,5:6)=[883 974];

%Rotation defined
angle=-5.75;

%Number of holes for this experiment - Used in Option 14
holes=10;


skip=1;
power=5;
power2=5;

% Define intensity range to be outlined - Bulk Interface
min_intensity = 0.5;
max_intensity = 1;

%Define intensity range to be outlined - Dispersed Interface
min_intensity2 = 0.1;
max_intensity2 = 1;


%Chose threshold for the error tolerance
threshold=40; %For Bulk Interface
threshold2=80; %For Dispersed Interface

exp_1=3;
exp_last=3;

x_f=817;

exp_move(1:4,1:2)=0;

if option3==2

exp_move(1,1)= 0;
exp_move(1,2)= 0;
exp_move(2,1) = 0;
exp_move(2,2) = 0;
exp_move(3,1) = 0;
exp_move(3,2) = 0;
exp_move(4,1) = 0;
exp_move(4,2) = 0;

end

%Window increase
n1 = 15;
med1=round(n1*n1*0.4);
n2= 5;

t_initial(1)=5;
t_initial(2)=2;
t_initial(3)=2;
t_initial(4)=2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

Experiment_dir_path{1} = strcat(main_path,'/Experiment_1010_0elev');
Experiment_dir_path{2} = strcat(main_path,'/Experiment_1020_0elev');
Experiment_dir_path{3} = strcat(main_path,'/Experiment_1010_10elev');
Experiment_dir_path{4} = strcat(main_path,'/Experiment_1020_10elev');

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Break%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Change file names - Changes the files names for all Calibration and
%Experiment (1-4) names to the ones recognized by code

if option == 1

    image_name(option2,density_range,exp_1,exp_last,Experiment_dir_path);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Break%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%To gray scale images. Rotation option possible, change the rotation number
%in user defined section

if option==2

    Rot_gray(option2,density_range,angle,exp_1,exp_last,Experiment_dir_path);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Break%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cropping and Image Subtraction - Use the information here to input into
% user defined

if option ==2.5

[row_bottom, row_top, column_left,mm_per_pix,column_right,x_f2,x_f]=cropping_imagesubtraction(home_points,option,option2,density_range,exp_1,exp_last,Experiment_dir_path, exp_move, option3);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Break%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Crops to the defined windows

if option ==3

[row_bottom, row_top, column_left,mm_per_pix,column_right,x_f2,x_f]=cropping_imagesubtraction(home_points,option,option2,density_range,exp_1,exp_last,Experiment_dir_path, exp_move, option3,n1, med1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Break%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calibration Curves

if option==4

calibration(density_range,power);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Break%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if option==5


experiment_threshholdcalc(option3,density_range,power,power2,max_intensity,min_intensity,min_intensity2,threshold,threshold2,exp_1,exp_last,Experiment_dir_path,x_f)


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Break%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if option==6

video(option3,option2, Experiment_dir_path, exp_1, exp_last,density_range)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Break%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if option==9

    [row_bottom, row_top, column_left,mm_per_pix,column_right,x_f2,x_f]=cropping_imagesubtraction(home_points,option,option2,density_range,exp_1,exp_last,Experiment_dir_path, exp_move, option3);
    cd(main_path)
    [idx] = xfdetermine (home_points,option,option2,density_range,row_bottom, row_top, column_left,mm_per_pix,column_right,x_f2,x_f,option3,exp_1,exp_last,Experiment_dir_path);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Break%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if option==13

[row_bottom, row_top, column_left,mm_per_pix,column_right,x_f2,x_f]=cropping_imagesubtraction(home_points,option,option2,density_range,exp_1,exp_last,Experiment_dir_path, exp_move, option3); 
cd(main_path)
surface_plots_ready(density_range,max_intensity,min_intensity,min_intensity2,threshold,threshold2,exp_1,exp_last,Experiment_dir_path,x_f,power,power2,mm_per_pix,t_initial)

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Break%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if option==14

area_ready(density_range,max_intensity,min_intensity,min_intensity2,threshold,threshold2,exp_1,exp_last,Experiment_dir_path,x_f,power,power2)

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Break%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if option==15
for exp=exp_1:exp_last

cd(exp_dir3_path{exp})
load("highestpoints.mat")
idx=find(saved_myanswer{exp}(:,2)>0.33,1,'last');
X=saved_myanswer{exp}(15:idx,2);
pray=polyfit(15:idx,X,1);
figure(11);plot(15:idx,polyval(pray,15:idx))

for i=12:60
saved_myanswer{exp}(i+1,2)=min(2900,31.62*(i+1)+saved_myanswer{exp}(12,2));
end
hold on;plot(15:idx,saved_myanswer{exp}(15:idx,2))
        cd(exp_dir3_path{exp})
        save('highestpoints.mat', 'saved_equation_bulk', 'saved_equation_dispersed','saved_highestpoints_bulk','saved_highestpoints_dispersed', 'saved_myanswer','fakeimages');
        clear fakeimages

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Break%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Break%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if option == 12
    % Assuming you have an array of images
    images = {0, 0, 0, 0, 0};

    A=18;
    n=3;

    if option3==1
        beg_exp=1;
        end_exp=1;
    elseif option3==2
        beg_exp=2;
        end_exp=2;
    elseif option3==3
        beg_exp=3;
        end_exp=3;
    elseif option3==4
        beg_exp=4;
        end_exp=4;
    elseif option3==5
        beg_exp=1;
        end_exp=5;
    end

    % Experiment 1
    cd(exp_dir2_path1)

    if option2==1

    for i=1:5
        cd(exp_dir2_path1)
         Variable = ['IMG',num2str(A+i*n),'_cropped.jpg'];
         img = imread(Variable);
         img = im2double(img);

         cd(exp_dir3_path1)
         load("highestpoints.mat")

        imshow(img)
         Size=size(img);

        row = Size(1);
        column =Size(2);

        Myanswer = saved_myanswer(A+i*n);
        Myanswer2= saved_myanswer(A+i*n,2);
        highest_points = saved_highestpoints_bulk_exp1(A+i*n,:,:);
        highest_points2= saved_highestpoints_dispersed_exp1(A+i*n,:,:); 

        hold on

        for col =1:Myanswer
            row=highest_points(col);
            scatter(col, row, 25, 'red', 'filled');
        end
        hold on
        for col =1:Myanswer2
            row=highest_points2(col);
            scatter(col, row, 25, 'blue', 'filled');
        end
        hold on

        scatter(Myanswer, 1:row,25, 'red','filled');
        scatter(Myanswer2,1:row, 25, 'blue', 'filled');
        hold off

        % Save the resulting figure as an image file
        output_image = ['Output.IMG',num2str(A+i*n),'_cropped.jpg']; 

        cd(exp_dir3_path1)
        saveas(gcf, output_image);
        close(gcf)
    end
    end
    cd(exp_dir3_path1)
    for i=1:5
         Variable = ['Output.IMG',num2str(A+i*n),'_cropped.jpg']; 
         img = imread(Variable);
         img = im2double(img);
         images{i}=img;
    end

    figure;
    for i = 1:length(images)
        subplot(length(images),1,i); % The first argument is the number of rows, the second is the number of columns, and the third is the index of the current plot
        imshow(images{i});
            pos = get(gca, 'Position');
    pos(2) = pos(2) + 0.0001; % adjust bottom
    pos(4) = pos(4) - 0.0001; % adjust height
    set(gca, 'Position', pos)
    end
end


%Buoyancy Calculations

if option==10

    cd(new_dir_path2)
    load("calibration_data.mat")

    B(1:4)=0;
    B2(1:4)=0;
    D(1:4)=0;
    D2(1:4)=0;

    for exp=exp_1:exp_last
    cd(main_path)
    load("idx.mat")
    cd(exp_dir3_path{exp})
    load("highestpoints.mat")
    %Convert image to double values
    Variable = ['IMG',num2str(idx(exp)),'_cropped.jpg'];
    cd(exp_dir2_path{exp})
    img=imread(Variable);
    imgstodisect = im2double(img);
    Size=size(img);


    row = Size(1);
    column = Size(2);

   
    density = zeros(row,column);
    %Fakeimage used to create a placeholder image, that takes into account the calibration
    fakeimage = zeros(row,column);
    %Display a fake image that stores the values of intensity corrected with
    %calibtration

        Pixel_intensities=nan(1,num_images);
        idx

        for j=1:column
            for k=1:row
                pixel_intensity = imgstodisect(k,j);

                % If pixel intensity = nan, density = nan

                if isnan(Intensities(k,j,1))
                    density(k,j)=nan;
                    continue
                end

                for q=1:length(density_range)
                    Pixel_intensities(q)=Intensities(k,j,q);
                end

                %In case 1 value is NaN, the interp1 function does not
                %work. Hence, we make a "ValidData" variable

                YY = density_range;
                validData = ~isnan(Pixel_intensities);
                YY = YY(validData);
                Pixel_intensities = Pixel_intensities(validData);

                density(k,j) = interp1(Pixel_intensities,YY,pixel_intensity,'linear');

                if density(k,j)<(0.04/20)
                    density(k,j)=0;
                end
                if density(k,j)>0.04 
                density(k,j)=0.04;
                end
                fakeimage(k,j)=density(k,j)*(1/0.04);
            end
        end

        %Fakeimage has NaN values, use neighbhour values to override NaN as an
%average

for j = 1:row
    for k = 1:column
        if isnan(fakeimage(j, k))
            % Define the region around the pixel to average
            minRow = max(1, j - 1);
            maxRow = min(j, j + 1);
            minCol = max(1, k - 1);
            maxCol = min(k, k + 1);

            % Extract the region and remove any NaN values
            region = fakeimage(minRow:maxRow, minCol:maxCol);
            region = region(~isnan(region));

            % Calculate the mean of the region and replace NaN
            if ~isempty(region)  % prevent mean() from NaN output if region is all NaNs
                fakeimage(j, k) = mean(region(:));
            end
        end
    end
end

        target = 3;
        y=round(target*x_f-(1/3)*x_f);
    % Remove NaN values
    X = 1:saved_myanswer{exp}(idx(exp),1);
    Y = row-saved_highestpoints_bulk{exp}(idx(exp),1:saved_myanswer{exp}(idx(exp),1));
    validData = ~isnan(Y);
    X = X(validData);
    Y = Y(validData);

    W = polyfit(X,Y,power);

x=1:saved_myanswer{exp}(idx(exp),1);
polyfun = @(x) polyval(W,x);

% modify it to work only on scalar inputs
scalar_polyfun = @(x) polyfun(x(1));

% define your range
xmin = 1;  % lower limit of x
xmax = saved_myanswer{exp}(idx(exp),1)+25;  % upper limit of x

% use fminbnd to find the x that minimizes your polynomial function
[xval, yval] = fminbnd(polyfun, xmin, xmax);

saved_myanswer{exp}(idx(exp),1)=round(xval);


    clear X
    figure(exp)
imshow(fakeimage)
hold on
    for i=1:saved_myanswer{exp}(idx(exp),1)
        myarray = fakeimage(1:row,i);
        GGGG=row-round(polyval(W,i));
        if round(row-polyval(W,i))<=0
            continue
        end
        myarray=fillmissing(myarray,'linear');

        for j=1:round(polyval(W,i))
            B(exp)=myarray(row+1-j)+B(exp);
        end
    end
    hold on

    saved_myanswer{exp}(idx(exp),2)=y;
    True= min(length(saved_highestpoints_dispersed{exp}),saved_myanswer{exp}(idx(exp),2));
    % Remove NaN values
    clear X Y
    X = 1:saved_myanswer{exp}(idx(exp),2);
    Y = row-saved_highestpoints_dispersed{exp}(idx(exp),1:True);
    validData = ~isnan(Y);
    X = X(validData);
    Y = Y(validData);

    W = polyfit(X,Y,power);
    clear X
    for i=1:True
        myarray = fakeimage(1:row,i);
        GGGG2=row-round(polyval(W,i));
        
        myarray=fillmissing(myarray,'linear');
        if round(row-polyval(W,i))<=0
            continue
        end
        
        for j=1:round(polyval(W,i))
            D(exp)=myarray(row+1-j)+D(exp);
        end
    end
    D(exp)=D(exp)-B(exp);
    hold off

    
    end

    disp(['For experiment 1, the Bulk Buoyancy is ' num2str(B) '.']);
    disp(['For experiment 1, the Dispersed Buoyancy is ' num2str(D) '.']);

    

end

if option==11 
    
    img_list = dir(fullfile(exp_dir2_path2, '*.jpg'));
    YYY=length(file_list2);
    Bulk(1:YYY)=0;
    Storedfitted(1:YYY,1:5)=0;

    cd(exp_dir2_path2)
    Variable = ['IMG',num2str(1),'_cropped.jpg'];
    imgstodisect = imread(Variable);
    imgstodisect = im2double(imgstodisect);
    Size=size(imgstodisect);

    column =Size(2);

    
    for lll=18:YYY

        cd(exp_dir2_path2)
        Variable = ['IMG',num2str(lll),'_cropped.jpg'];
        imgstodisect = imread(Variable);
        imgstodisect = im2double(imgstodisect);
        Size=size(imgstodisect);

        row = Size(1);
        column =Size(2);



        density = zeros(row,column);
        fakeimage = zeros(row,column);
        %Display a fake image that stores the values of intensity corrected with
        %calibtration

        for j=1:column
            for k=1:row
                pixel_intensity = imgstodisect(k,j);
                density(k,j) = polyval(P((((j-1)*row)+k),:), pixel_intensity);
                if density(k,j)<0.0
                    density(k,j)=0;
                end
                if density(k,j)>0.04 
                density(k,j)=0.04;
                end
                fakeimage(k,j)=density(k,j)*(1/0.04);
            end
        end

         for col = 1:column
            for r=1:row
                if r<11
                    m=mean(fakeimage(r:r+10,col));
                end
                if r>10 &&  r<row-10
                    m=mean(fakeimage(r-10:r+10,col));
                end
                if r>row-11
                    m=mean(fakeimage(r-10:r,col));
                end

                if fakeimage(r,col)+0.15<m || fakeimage(r,col)-0.15>m
                    fakeimage(r,col)=m;
                end
            end
            K=polyfit(1:row,fakeimage(1:row,col),5);
            for r=1:row
                Fit(r)=polyval(K,r);
            end
            for r=1:row
                 if fakeimage(r,col)+0.25<Fit(r) || fakeimage(r,col)-0.25>Fit(r)
                     fakeimage(r,col)=Fit(r);
                 end
            end
         end
    end


        % Find highest point in each column within the intensity range
        [height, width] = size(fakeimage);
        highest_points = zeros(1, width);
        highest_points2 = zeros(1, width);
        for col = 1:width
            col_intensity = fakeimage(:, col);
            col_indices = find(col_intensity >= min_intensity & col_intensity <= max_intensity);
            if isempty(col_indices)
                continue;
            end
            R=1;
            A=col_indices(R);
            J=10;
            if J<length(col_indices)
                while A+20<col_indices(J)
                    R=R+1;
                    A=col_indices(R);
                    J=R+9;
                    if J>length(col_indices)
                        A=col_indices(end);
                        break
                    end
                end
            elseif J>=length(col_indices)
                A=col_indices(end);
            end
            
            highest_points(col) = A;
        end
        Myanswer= width;

        for col=1:length(highest_points)-51
            average=mean(highest_points(col:col+50));
            if average>row-15
                Myanswer=col;
                break
            end
        end


        G = polyfit(1:Myanswer, highest_points(1:Myanswer),power);

        Fit=zeros(Myanswer,1);
        highest_pointscal=zeros(Myanswer,1);



        for i=1:Myanswer
            Fit(i)=polyval(G,i);
        end

        for i=1:Myanswer
            highest_pointscal(i)=round(Fit(i));
        end

        residuals=zeros(Myanswer);

        % Remove the outliers
        for i=1:Myanswer
            residuals(i)= abs(highest_points(i) - highest_pointscal(i));
            if residuals(i)>threshold
                highest_points(i) = NaN;
            end
        end 

%For dispersed interface

        for col = 1:width
            col_intensity = fakeimage(:, col);
            col_indices2 = find(col_intensity >= min_intensity2 & col_intensity <= max_intensity);
            if isempty(col_indices2)
                continue;
            end
            R=1;
            A=col_indices(R);
            J=15;
            if J<length(col_indices2)
                while A+30<col_indices2(J)
                    R=R+1;
                    A=col_indices2(R);
                    J=R+9;
                    if J>length(col_indices2)
                        A=col_indices2(end);
                        break
                    end
                end
            elseif J>=length(col_indices2)
                A=col_indices2(end);
            end
            
            highest_points2(col) = A;
        end
        Myanswer2 = width;

        for col=1:length(highest_points2)-51
            average=mean(highest_points2(col:col+50));
            if average>row-40
                Myanswer2=col;
                break
            end
        end


        G2 = polyfit(1:Myanswer2, highest_points2(1:Myanswer2),power2);


        Fit2=zeros(Myanswer2,1);
        highest_pointscal2=zeros(Myanswer2,1);



        for i=1:Myanswer2
            Fit2(i)=polyval(G2,i);
        end

        for i=1:Myanswer2
            highest_pointscal2(i)=round(Fit2(i));
        end

        residuals2=zeros(Myanswer2);

        % Remove the outliers
        for i=1:Myanswer2
            residuals2(i)= abs(highest_points2(i) - highest_pointscal2(i));
            if residuals2(i)>threshold2
                highest_points2(i) = NaN;
            end
        end 



        saved_myanswer(lll,1)=Myanswer;
        saved_myanswer(lll,2)=Myanswer2;


        % Find highest point in each column within the intensity range
        [height, width] = size(fakeimage);

        k=0;
        Z=0;

        for col = 1:column
            for r=1:row
                if r<11
                    m=mean(fakeimage(r:r+10,col));
                end
                if r>10 &&  r<row-10
                    m=mean(fakeimage(r-10:r+10,col));
                end
                if r>row-11
                    m=mean(fakeimage(r-10:r,col));
                end

                if fakeimage(r,col)+0.25<m || fakeimage(r,col)-0.25>m
                    fakeimage(r,col)=m;
                end
            end
            K=polyfit(1:row,fakeimage(1:row,col),5);
            for r=1:row
                Fit(r)=polyval(K,r);
            end
            for r=1:row
                 if fakeimage(r,col)+0.25<Fit(r) || fakeimage(r,col)-0.25>Fit(r)
                     fakeimage(r,col)=Fit(r);
                 end
            end

            K_bulk=polyfit(1:row,fakeimage(1:row,col),5);
            K_disp=K_bulk;

            K_bulk(end)=K_bulk(end)-min_intensity;

            K_disp(end)=K_disp(end)-min_intensity2;

            Roots_bulk=roots(K_bulk);
            Roots_bulk=Roots_bulk(imag(Roots_bulk)== 0 & Roots_bulk >= 1 & Roots_bulk <= row);

            Roots_disp=roots(K_disp);
            Roots_disp=Roots_disp(imag(Roots_disp)== 0 & Roots_disp >= 1 & Roots_disp <= row);

            if k<01
                if isempty(Roots_bulk)
                    k=k+1;
                    mycol=col;
                    continue
                else
                    highest_points(col)=Roots_bulk(1);
                end
            end


            if Z<101
                if isempty(Roots_disp)
                    Z=Z+1;
                    mycol2=col;
                    continue
                else
                    highest_points2(col)=Roots_disp(1);
                end
            end

        end
            
        toprow=highest_points(1:Myanswer);

        % Remove NaN values
        YY = 1:Myanswer;
        validData = ~isnan(toprow(1:Myanswer));
        YY = YY(validData);
        toprow = toprow(validData);


        G = polyfit(YY, toprow,power);

        Fit=zeros(Myanswer,1);
        highest_pointscal=zeros(Myanswer,1);


        for i=1:Myanswer
            Fit(i)=polyval(G,i);
        end

        for i=1:Myanswer
            highest_pointscal(i)=round(Fit(i));
        end

        residuals=zeros(Myanswer);

        % Remove the outliers
        for i=1:Myanswer
            if isnan(highest_points)
                continue
            else
                residuals(i)= abs(highest_points(i) - highest_pointscal(i));
                if residuals(i)>threshold
                    highest_points(i) = NaN;
                end
            end
        end

%For dispersed interface

        toprow2=highest_points2(1:Myanswer2);

        % Remove NaN values
        YY2 = 1:Myanswer2;
        validData2 = ~isnan(toprow2);
        YY2 = YY2(validData2);
        toprow2 = toprow2(validData2);

        G2 = polyfit(YY2, toprow2,power2);


        Fit2=zeros(Myanswer2,1);
        highest_pointscal2=zeros(Myanswer2,1);



        for i=1:Myanswer2
            Fit2(i)=polyval(G2,i);
        end

        for i=1:Myanswer2
            highest_pointscal2(i)=round(Fit2(i));
        end

        residuals2=zeros(Myanswer2);

        % Remove the outliers
        for i=1:Myanswer2
            if isnan(highest_points)
                continue
            else
                residuals2(i)= abs(highest_points2(i) - highest_pointscal2(i));
                if residuals2(i)>threshold2
                    highest_points2(i) = NaN;
                end
            end
        end 
%         imshow(imgstodisect)
%         hold on
%         for col =1:Myanswer
%             row=highest_points(col);
%             scatter(col, row, 25, 'red', 'filled');
%         end
%         hold on
%         for col =1:Myanswer2
%             row=highest_points2(col);
%             scatter(col, row, 25, 'blue', 'filled');
%         end
%         hold off
% 
%         % Save the resulting figure as an image file
%         output_image = ['Output.IMG',num2str(lll),'_cropped.jpg'];


        cd(exp_dir3_path2)
%         saveas(gcf, output_image);
%         close(gcf)

        %saved power
        saved_equation_bulk_exp2(lll,1:power+1)= G;
        saved_equation_dispersed_exp2(lll,1:power+1)= G2;
        saved_highestpoints_bulk_exp2(lll,:,:)= highest_points;
        saved_highestpoints_dispersed_exp2(lll,:,:)= highest_points2;
        saved_myanswer(lll,3)=mycol;
        saved_myanswer(lll,4)=mycol2;
        

end

if option == 7

    cd(new_dir_path2)
    load("calibration_data.mat")

    if ~exist("exp_dir4_path1", 'dir')
        mkdir(exp_dir4_path1);
    end

    if ~exist("exp_dir4_path2", 'dir')
        mkdir(exp_dir4_path2);
    end

    if ~exist("exp_dir4_path3", 'dir')
        mkdir(exp_dir4_path3);
    end

    if ~exist("exp_dir4_path4", 'dir')
        mkdir(exp_dir4_path4);
    end

    %     %Experiment 1

    if option3==1 || option3==5

    img_list = dir(fullfile(exp_dir2_path1, '*.jpg'));
    YYY=length(file_list1);
    Bulk(1:YYY)=0;
    Storedfitted(1:YYY,1:5)=0;

    cd(exp_dir2_path1)
    Variable = ['IMG',num2str(1),'_cropped.jpg'];
    imgstodisect = imread(Variable);
    imgstodisect = im2double(imgstodisect);
    Size=size(imgstodisect);

    column =Size(2);
    row = Size(1);

    saved_equation_bulk_exp1= zeros(YYY,power+1);
    saved_equation_dispersed_exp1= zeros(YYY,power+1);
    saved_highestpoints_bulk_exp1= zeros(YYY,column);
    saved_highestpoints_dispersed_exp1= zeros(YYY,column);
    savenoseposition=zeros(YYY);
    savefittedequation_bulk=zeros(YYY,power+1);
    savefittedequation_dispersed=zeros(YYY,power+1);


    for lll=1:YYY

        cd(exp_dir2_path1)
        Variable = ['IMG',num2str(lll),'_cropped.jpg'];
        imgstodisect = imread(Variable);
        imgstodisect = im2double(imgstodisect);
        Size=size(imgstodisect);

        row = Size(1);
        column =Size(2);
        density = zeros(row,column);
        fakeimage = zeros(row,column);
        %Display a fake image that stores the values of intensity corrected with
        %calibtration

        Fit = zeros(row);

        for j=1:column
            for k=1:row
                pixel_intensity = imgstodisect(k,j);
                density(k,j) = polyval(P((((j-1)*row)+k),:), pixel_intensity);
                if density(k,j)<0.0
                    density(k,j)=0;
                end
                if density(k,j)>0.04 
                density(k,j)=0.04;
                end
                fakeimage(k,j)=density(k,j)*(1/0.04);
            end
        end


        % Find highest point in each column within the intensity range
        [height, width] = size(fakeimage);
        myheight = zeros(1, width);
        myheight2 = zeros(1, width);
        dummycheck=0;


        %Trapezoidal Integral
        for col = 1:column
            for r=1:row
                if r<11
                    m=mean(fakeimage(r:r+10,col));
                end
                if r>10 &&  r<row-10
                    m=mean(fakeimage(r-10:r+10,col));
                end
                if r>row-11
                    m=mean(fakeimage(r-10:r,col));
                end

                if fakeimage(r,col)+0.25<m || fakeimage(r,col)-0.25>m
                    fakeimage(r,col)=m;
                end
            end
            K=polyfit(1:row,fakeimage(1:row,col),5);
            for r=1:row
                Fit(r)=polyval(K,r);
            end
            for r=1:row
                 if fakeimage(r,col)+0.4<Fit(r) || fakeimage(r,col)-0.4>Fit(r)
                     fakeimage(r,col)=Fit(r);
                 end
            end

            myheight(col)=trapz(1:row,fakeimage(1:row,col));

            myheight2(col)=myheight(col)+trapz(1:row-round(myheight(col)),fakeimage(1:row-round(myheight(col)),col));

%             if col<column-200 && col>200
%                 if mean(myheight2(col-10:col)) > mean(myheight2(col-200:col))
%                     mycol=col-200;
%                     dummycheck=dummycheck+1;
%                     if dummycheck==10
%                         break
%                     end
%                 end
%             else
%                 mycol=column;
%             end
        end

        for col=51:column
            if mean(myheight2(col-50:col)) < 1.1*mean(myheight2(col:column))
                mycol=col-50;
                break
            end
        end

        K2=polyfit(1:mycol,myheight(1:mycol),power);
        K2_2=polyfit(1:mycol,myheight2(1:mycol),power);



        fig = figure('visible', 'off'); % Create an invisible figure
        imshow(imgstodisect)
        hold on
        for col =1:mycol
            rowing=row - round(myheight(col));
            scatter(col, rowing, 25, 'red', 'filled');
        end
        for col =1:mycol
            rowing=row - round(myheight2(col));
            scatter(col, rowing, 25, 'blue', 'filled');
        end
        hold off

        % Save the resulting figure as an image file
        output_image = ['Output.IMG',num2str(lll),'_cropped.jpg'];

        cd(exp_dir4_path1)
        exportgraphics(gca, output_image);

        %saved power
        saved_highestpoints_bulk_exp1(lll,:,:)= myheight;
        saved_highestpoints_dispersed_exp1(lll,:,:)= myheight2;
        savenoseposition(lll)=mycol;
        savefittedequation_bulk(lll,:)=K2;
        savefittedequation_dispersed(lll,:)=K2_2;

    end
    save('highestpoints.mat','saved_highestpoints_bulk_exp1','saved_highestpoints_dispersed_exp1');
    end

        %     %Experiment 2

    if option3==2 || option3==5

    img_list = dir(fullfile(exp_dir2_path2, '*.jpg'));
    YYY=length(file_list2);
    Bulk(1:YYY)=0;
    Storedfitted(1:YYY,1:5)=0;

    cd(exp_dir2_path2)
    Variable = ['IMG',num2str(1),'_cropped.jpg'];
    imgstodisect = imread(Variable);
    imgstodisect = im2double(imgstodisect);
    Size=size(imgstodisect);

    column =Size(2);
    row = Size(1);

    saved_equation_bulk_exp1= zeros(YYY,power+1);
    saved_equation_dispersed_exp1= zeros(YYY,power+1);
    saved_highestpoints_bulk_exp1= zeros(YYY,column);
    saved_highestpoints_dispersed_exp1= zeros(YYY,column);
    savenoseposition=zeros(YYY);
    savefittedequation_bulk=zeros(YYY,power+1);
    savefittedequation_dispersed=zeros(YYY,power+1);

    for lll=10:YYY

        cd(exp_dir2_path2)
        Variable = ['IMG',num2str(lll),'_cropped.jpg'];
        imgstodisect = imread(Variable);
        imgstodisect = im2double(imgstodisect);
        Size=size(imgstodisect);

        row = Size(1);
        column =Size(2);
        density = zeros(row,column);
        fakeimage = zeros(row,column);
        %Display a fake image that stores the values of intensity corrected with
        %calibtration

        Fit = zeros(row);

        for j=1:column
            for k=1:row
                pixel_intensity = imgstodisect(k,j);
                density(k,j) = polyval(P((((j-1)*row)+k),:), pixel_intensity);
                if density(k,j)<0.0
                    density(k,j)=0;
                end
                if density(k,j)>0.04 
                density(k,j)=0.04;
                end
                fakeimage(k,j)=density(k,j)*(1/0.04);
            end
        end


        % Find highest point in each column within the intensity range
        [height, width] = size(fakeimage);
        myheight = zeros(1, width);
        myheight2 = zeros(1, width);
        dummycheck=0;


        %Trapezoidal Integral
        for col = 1:column
            for r=1:row
                if r<11
                    m=mean(fakeimage(r:r+10,col));
                end
                if r>10 &&  r<row-10
                    m=mean(fakeimage(r-10:r+10,col));
                end
                if r>row-11
                    m=mean(fakeimage(r-10:r,col));
                end

                if fakeimage(r,col)+0.25<m || fakeimage(r,col)-0.25>m
                    fakeimage(r,col)=m;
                end
            end
            K=polyfit(1:row,fakeimage(1:row,col),5);
            for r=1:row
                Fit(r)=polyval(K,r);
            end
            for r=1:row
                 if fakeimage(r,col)+0.4<Fit(r) || fakeimage(r,col)-0.4>Fit(r)
                     fakeimage(r,col)=Fit(r);
                 end
            end

            myheight(col)=trapz(1:row,fakeimage(1:row,col));

            myheight2(col)=myheight(col)+trapz(1:row-round(myheight(col)),fakeimage(1:row-round(myheight(col)),col));

        end


        for col=51:column
            if mean(myheight2(col-50:col)) < 1.1*mean(myheight2(col:column))
                mycol=col-50;
                break
            end
        end

        K2=polyfit(1:mycol,myheight(1:mycol),power);
        K2_2=polyfit(1:mycol,myheight2(1:mycol),power);


        fig = figure('visible', 'off'); % Create an invisible figure
        imshow(imgstodisect)
        hold on
        for col =1:mycol
            rowing=row - round(myheight(col));
            scatter(col, rowing, 25, 'red', 'filled');
        end
        for col =1:mycol
            rowing=row - round(myheight2(col));
            scatter(col, rowing, 25, 'blue', 'filled');
        end
        hold off

        % Save the resulting figure as an image file
        output_image = ['Output.IMG',num2str(lll),'_cropped.jpg'];

        cd(exp_dir4_path2)
        exportgraphics(gca, output_image);

        %saved power
        saved_highestpoints_bulk_exp1(lll,:,:)= myheight;
        saved_highestpoints_dispersed_exp1(lll,:,:)= myheight2;
        savenoseposition(lll)=mycol;
        savefittedequation_bulk(lll,:)=K2;
        savefittedequation_dispersed(lll,:)=K2_2;

    end
    save('highestpoints.mat','saved_highestpoints_bulk_exp1','saved_highestpoints_dispersed_exp1','savenoseposition');
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Testing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     %Experiment 4

    if option3==4 || option3==5
    
    img_list = dir(fullfile(exp_dir2_path4, '*.jpg'));
    YYY=length(file_list4);
    Bulk(1:YYY)=0;
    Storedfitted(1:YYY,1:5)=0;

    cd(exp_dir2_path4)
    Variable = ['IMG',num2str(1),'_cropped.jpg'];
    imgstodisect = imread(Variable);
    imgstodisect = im2double(imgstodisect);
    Size=size(imgstodisect);

    row = Size(1);
    column =Size(2);

    saved_equation_bulk_exp4= zeros(YYY,power+1);
    saved_equation_dispersed_exp4= zeros(YYY,power+1);
    saved_highestpoints_bulk_exp4= zeros(YYY,column);
    saved_highestpoints_dispersed_exp4= zeros(YYY,column);


    for lll=10:YYY

        cd(exp_dir2_path4)
        Variable = ['IMG',num2str(lll),'_cropped.jpg'];
        imgstodisect = imread(Variable);
        imgstodisect = im2double(imgstodisect);
        Size=size(imgstodisect);

        row = Size(1);
        column =Size(2);
        density = zeros(row,column);
        fakeimage = zeros(row,column);
        %Display a fake image that stores the values of intensity corrected with
        %calibtration

        for j=1:column
            for k=1:row
                pixel_intensity = imgstodisect(k,j);
                density(k,j) = polyval(P((((j-1)*row)+k),:), pixel_intensity);
                if density(k,j)<0.0
                    density(k,j)=0;
                end
                if density(k,j)>0.04 
                density(k,j)=0.04;
                end
                fakeimage(k,j)=density(k,j)*(1/0.04);
            end
        end


        % Find highest point in each column within the intensity range
        [height, width] = size(fakeimage);
        highest_points = zeros(1, width);
        highest_points2 = zeros(1, width);

        k=0;
        Z=0;

        for col = 1:column
            for r=1:row
                if r<11
                    m=mean(fakeimage(r:r+10,col));
                end
                if r>10 &&  r<row-10
                    m=mean(fakeimage(r-10:r+10,col));
                end
                if r>row-11
                    m=mean(fakeimage(r-10:r,col));
                end

                if fakeimage(r,col)+0.25<m || fakeimage(r,col)-0.25>m
                    fakeimage(r,col)=m;
                end
            end
            K=polyfit(1:row,fakeimage(1:row,col),5);
            for r=1:row
                Fit(r)=polyval(K,r);
            end
            for r=1:row
                 if fakeimage(r,col)+0.25<Fit(r) || fakeimage(r,col)-0.25>Fit(r)
                     fakeimage(r,col)=Fit(r);
                 end
            end
            K_bulk=polyfit(1:row,fakeimage(1:row,col),5);
            K_disp=K_bulk;

            K_bulk(end)=K_bulk(end)-min_intensity;

            K_disp(end)=K_disp(end)-min_intensity2;

            Roots_bulk=roots(K_bulk);
            Roots_bulk=Roots_bulk(imag(Roots_bulk)== 0 & Roots_bulk >= 1 & Roots_bulk <= row);

            Roots_disp=roots(K_disp);
            Roots_disp=Roots_disp(imag(Roots_disp)== 0 & Roots_disp >= 1 & Roots_disp <= row);

            if k<01
                if isempty(Roots_bulk)
                    k=k+1;
                    mycol=col;
                    continue
                else
                    highest_points(col)=Roots_bulk(1);
                end
            end


            if Z<101
                if isempty(Roots_disp)
                    Z=Z+1;
                    mycol2=col;
                    continue
                else
                    highest_points2(col)=Roots_disp(1);
                end
            end
        end


        G = polyfit(1:col, highest_points(1:col),power);

        Fit=zeros(mycol,1);
        highest_pointscal=zeros(mycol,1);


        for i=1:mycol
            Fit(i)=polyval(G,i);
        end

        for i=1:mycol
            highest_pointscal(i)=round(Fit(i));
        end

        residuals=zeros(mycol);

        % Remove the outliers
        for i=1:mycol
            residuals(i)= abs(highest_points(i) - highest_pointscal(i));
            if residuals(i)>threshold
                highest_points(i) = NaN;
            end
        end 

%For dispersed interface


        G2 = polyfit(1:mycol2, highest_points2(1:mycol2),power2);


        Fit2=zeros(mycol2,1);
        highest_pointscal2=zeros(mycol2,1);



        for i=1:mycol2
            Fit2(i)=polyval(G2,i);
        end

        for i=1:mycol2
            highest_pointscal2(i)=round(Fit2(i));
        end

        residuals2=zeros(mycol2);

        % Remove the outliers
        for i=1:mycol2
            residuals2(i)= abs(highest_points2(i) - highest_pointscal2(i));
            if residuals2(i)>threshold2
                highest_points2(i) = NaN;
            end
        end 
        imshow(imgstodisect)
        hold on
        for col =1:mycol
            row=highest_points(col);
            scatter(col, row, 25, 'red', 'filled');
        end
        hold on
        for col =1:mycol2
            row=highest_points2(col);
            scatter(col, row, 25, 'blue', 'filled');
        end
        hold off

        % Save the resulting figure as an image file
        output_image = ['Output.IMG',num2str(lll),'_cropped.jpg'];

        cd(exp_dir3_path4)
        saveas(gcf, output_image);
        close(gcf)

        %saved power

        saved_equation_bulk_exp4(lll,1:power+1)= G;
        saved_equation_dispersed_exp4(lll,1:power+1)= G2;
        saved_highestpoints_bulk_exp4(lll,:,:)= highest_points;
        saved_highestpoints_dispersed_exp4(lll,:,:)= highest_points2;

    end

    save('highestpoints.mat', 'saved_equation_bulk_exp4', 'saved_equation_dispersed_exp4','saved_highestpoints_bulk_exp4','saved_highestpoints_dispersed_exp4');
    
    end

end

if option==8
    %For Experiment 1
    cd(exp_dir4_path1)


    file_list1_vid= dir(fullfile(exp_dir4_path1, '*.jpg'));
    num_images = length(file_list1_vid);

    output_video='output_video.mp4';
    video_writer = VideoWriter(output_video,'MPEG-4');
    video_writer.FrameRate = 5; % Set the desired frame rate (frames per second)
    open(video_writer);

    for i = 1:num_images
        image = imread(fullfile(exp_dir4_path1, file_list1_vid(i).name));
        writeVideo(video_writer, image);
    end
    close(video_writer);

end

