

function [] = experiment_threshholdcalc(option3,density_range,power,power2,max_intensity,min_intensity,min_intensity2,threshold,threshold2,exp_1,exp_last,Experiment_dir_path,x_f)

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



file_list{1} = dir(fullfile(Experiment_dir_path{1}, '*.jpg'));
file_list{2} = dir(fullfile(Experiment_dir_path{2}, '*.jpg'));
file_list{3} = dir(fullfile(Experiment_dir_path{3}, '*.jpg'));
file_list{4} = dir(fullfile(Experiment_dir_path{4}, '*.jpg'));

   cd(new_dir_path2)
   load("calibration_data.mat")

    for exp=exp_1:exp_last
        if ~exist("exp_dir3_path{exp}", 'dir')
            mkdir(exp_dir3_path{exp});
        end
    end


    for exp=exp_1:exp_last

    YYY=length(file_list{exp});
    Bulk(1:YYY)=0;
    Storedfitted(1:YYY,1:5)=0;

    cd(exp_dir2_path{exp})
    Variable = ['IMG',num2str(1),'_cropped.jpg'];
    imgstodisect = imread(Variable);
    imgstodisect = im2double(imgstodisect);
    Size=size(imgstodisect);

    column =Size(2);

    saved_equation_bulk{exp}= zeros(YYY,power+1);
    saved_equation_dispersed{exp}= zeros(YYY,power+1);
    saved_highestpoints_bulk{exp}= zeros(YYY,column);
    saved_highestpoints_dispersed{exp}= zeros(YYY,column);
    saved_myanswer{exp}= zeros(YYY,4);
    fakeimages{YYY} = zeros(Size(1),Size(2));

    number=0;

%     if exist("highestpoints.mat", 'file')
%             load('highestpoints.mat')
%     end

    for lll=5:YYY



        cd(exp_dir2_path{exp})
  

        %Convert image to double values
        Variable = ['IMG',num2str(lll),'_cropped.jpg'];
        imgstodisect = imread(Variable);
        imgstodisect = im2double(imgstodisect);
        Size=size(imgstodisect);

        row = Size(1);
        column =Size(2);
        density = zeros(row,column);

        %Fakeimage used to create a placeholder image, that takes into
        %account the calibration
        fakeimage = nan(row,column);
        %Display a fake image that stores the values of intensity corrected with
        %calibtration

        Pixel_intensities=nan(1,num_images);

%         if lll>5  && lll<28
%              continue
%         end

        lll
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

                if density(k,j)<(0.04/15)
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



        % Find highest point in each column within the intensity range
        [height, width] = size(fakeimage);
        highest_points(1,1:width) = height; %The bulk interface
        highest_points2(1,1:width) = height; %The dispersed Interface
        for col = 1:width
            col_intensity = fakeimage(:, col);
            col_indices = find(col_intensity >= min_intensity & col_intensity <= max_intensity);
            if isempty(col_indices) %If none meet criteria, leave at 0
                continue;
            end
            R=1;

            A=col_indices(R);

            if length(col_indices)>R+2
            A2=col_indices(R+2);

            if A2>A+15
                A=A2;
                R=R+2;
            end
            end

%             J=15;

            J=min(5,round((height-A)*0.1));


            if J<length(col_indices) && J>0
                while A+(3*J)<col_indices(J)
                    R=R+1;
                    A=col_indices(R);
                    J=J+1;
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
%Myanswer is the nose position for bulk interface

        Myanswer= width;

%Looks at average of the next 50 rows, to see if they are
%increasing/decreasing

%The number in "row-number" is chosen based on judgement

        for col=1:length(highest_points)-101
            average=mean(highest_points(col:col+100));
            if average>row-10
                Myanswer=col+100;
                break
            end
        end

        Myanswer_forcurve = min(Myanswer+100,width);

        G = polyfit(1:Myanswer_forcurve, highest_points(1:Myanswer_forcurve),power);

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

        if number==0
            number=mean(sum(fakeimage(:,width-300:width)))+1.5;
        end

        for col = 1:width
            col_intensity = fakeimage(:, col);
            col_indices2 = find(col_intensity >= min_intensity2 & col_intensity <= max_intensity);
            if isempty(col_indices2)
                continue;
            end
            R=1;
            A=col_indices2(R);

            if length(col_indices2)>=R+2
            A2=col_indices2(R+2);

            if A2>A+15
                A=A2;
                R=R+2;
            end
            end

%             J=15;

            J=min(10,round((height-A)*0.2));

            if J<10
                xrayxray=0;
            end

            if J<length(col_indices2) && J>0
                 while A+(2*J)<col_indices2(J)
                    R=R+1;
                    A=col_indices2(R);
                    J=J+1;
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


        for col=1:length(highest_points2)-101
            average=mean(sum(fakeimage(:,col:col+100)));
            if average<number
                Myanswer2=col+100;
                break
            end
        end
        

        Myanswer_forcurve2 = min(Myanswer2+100,width);

        G2 =  polyfit(1:Myanswer_forcurve2, highest_points2(1:Myanswer_forcurve2),power2);


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
       




        



%         % Find highest point in each column within the intensity range
%         [height, width] = size(fakeimage);
% 
%         k=0;
%         Z=0;
% 
%         for col = 1:column
%             for r=1:row
%                 if r<11
%                     m=mean(fakeimage(r:r+10,col));
%                 end
%                 if r>10 &&  r<row-10
%                     m=mean(fakeimage(r-10:r+10,col));
%                 end
%                 if r>row-11
%                     m=mean(fakeimage(r-10:r,col));
%                 end
% 
%                 if fakeimage(r,col)+0.25<m || fakeimage(r,col)-0.25>m
%                     fakeimage(r,col)=m;
%                 end
%             end
%             K=polyfit(1:row,fakeimage(1:row,col),5);
%             for r=1:row
%                 Fit(r)=polyval(K,r);
%             end
%             for r=1:row
%                  if fakeimage(r,col)+0.25<Fit(r) || fakeimage(r,col)-0.25>Fit(r)
%                      fakeimage(r,col)=Fit(r);
%                  end
%             end

%             topbulkrow(:)=fakeimage(1:row, col);
% 
%             % Remove NaN values
%             YY_bulk = 1:row;
%             validData = ~isnan(topbulkrow(1:row));
%             YY_bulk = YY_bulk(validData);
%             topbulkrow = topbulkrow(validData);
% 
% 
% 
%             K_bulk=polyfit(YY_bulk,topbulkrow,power);
%             clear topbulkrow validData YY_bulk 
%             K_disp=K_bulk;
% 
%             K_bulk(end)=K_bulk(end)-min_intensity;
% 
%             K_disp(end)=K_disp(end)-min_intensity2;
% 
%             Roots_bulk=roots(K_bulk);
%             Roots_bulk=Roots_bulk(imag(Roots_bulk)== 0 & Roots_bulk >= 1 & Roots_bulk <= row);
% 
%             Roots_disp=roots(K_disp);
%             Roots_disp=Roots_disp(imag(Roots_disp)== 0 & Roots_disp >= 1 & Roots_disp <= row);

%             if k<01
%                 if isempty(Roots_bulk)
%                     k=k+1;
%                     mycol=col;
%                     continue
%                 else
%                     highest_points(col)=Roots_bulk(1);
%                 end
%             end
% 
% 
%             if Z<101
%                 if isempty(Roots_disp)
%                     Z=Z+1;
%                     mycol2=col;
%                     continue
%                 else
%                     highest_points2(col)=Roots_disp(1);
%                 end
%             end
% 
%         end
%             
        toprow=highest_points(1:Myanswer);

        % Remove NaN values
        YY = 1:Myanswer;
        validData = ~isnan(toprow(1:Myanswer));
        YY = YY(validData);
        toprow_test = height-toprow(validData);
        toprow = toprow(validData);


        G = polyfit(YY, toprow_test,power);

        x=1:Myanswer;
        polyfun = @(x) polyval(G,x);

        % modify it to work only on scalar inputs
        scalar_polyfun = @(x) polyfun(x(1));

        % define your range
        xmin = 1;  % lower limit of x
        xmax = Myanswer; % upper limit of x

        [xval, yval] = fminbnd(polyfun, xmin, xmax);


        while yval > 1 && xmax > Myanswer+50

        xmax = xmax+1; % upper limit of x
        [xval, yval] = fminbnd(polyfun, xmin, xmax);
        end

        if yval <= 1
            xmax = xmax-1; % upper limit of x
            [xval, yval] = fminbnd(polyfun, xmin, xmax);
        end




        if xval ~= Myanswer && yval>0
            Myanswer = round(xval);
            disp("Myanswer changed")
        end

        clear G

        G = polyfit(YY, toprow,power);

        Fit=zeros(Myanswer,1);
        highest_pointscal=zeros(Myanswer,1);


        if Myanswer+50>width
            display(["Breaking at " num2str(lll)])
            break
        end


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
        toprow2_test = height-toprow2(validData2);
        toprow2 = toprow2(validData2);


        G2 = polyfit(YY2, toprow2_test,power2);

        x2=1:Myanswer2;
        polyfun = @(x2) polyval(G2,x2);

        % modify it to work only on scalar inputs
        scalar_polyfun = @(x2) polyfun(x2(1));

        % define your range
        xmin = 1;  % lower limit of x
        xmax = Myanswer2; % upper limit of x

        % use fminbnd to find the x that minimizes your polynomial function
        [xval2, yval2] = fminbnd(polyfun, xmin, xmax);


        while yval2 > 1 && xmax > Myanswer2+50

        xmax = xmax+1; % upper limit of x
        [xval2, yval2] = fminbnd(polyfun, xmin, xmax);

        end

        if yval2 <= 1
            xmax = xmax-1; % upper limit of x
            [xval2, yval2] = fminbnd(polyfun, xmin, xmax);
        end




        if xval2 ~= Myanswer2 && yval2>0
            Myanswer2 = round(xval2);
            disp("Myanswer2 changed")
        end

       
        clear G2

        G2 = polyfit(YY2, toprow2,power2);
        Fit2=zeros(Myanswer2,1);
        highest_pointscal2=zeros(Myanswer2,1);

        if Myanswer2+50>width
            display(["Breaking at " num2str(lll)])
            break
        end


        for i=1:Myanswer2
            Fit2(i)=polyval(G2,i);
        end

        for i=1:Myanswer2
            highest_pointscal2(i)=round(Fit2(i));
        end

        residuals2=zeros(Myanswer2);

        % Remove the outliers
        for i=1:Myanswer2
            if isnan(highest_points2)
                continue
            else
                residuals2(i)= abs(highest_points2(i) - highest_pointscal2(i));
                if residuals2(i)>threshold2
                    highest_points2(i) = NaN;
                end
            end
        end 


        for col=1:column
            if isnan(fakeimage(:,col))
            fakeimage(:,col)=fillmissing(fakeimage(:,col));
            end
        end

        highest_points2=fillmissing(highest_points2,"linear");


        if Myanswer2 < Myanswer
            Myanswer2 = Myanswer +10;
        end
      


mmm=fillmissing(highest_points,'linear');
highest_points=mmm;
newImage = fakeimage;

mmm=round((row+mmm)/2);

rice=0.9;
for c = 1:Myanswer
    m=round(mmm(c));
    n=find(fakeimage(m:row,c)>=rice,10,'first');
    n=n+m-1;
    if length(n) >= 10
    for r = n(5):row
        if fakeimage(r,c) < rice
            % Get index range of 5 pixels above and below, avoiding indices out of range
            range_below = max(1,r-5):r-1;
            range_above = r+1:min(row,r+5);
            
            % Select only pixels with value greater than 0.5
            pixels_below = fakeimage(range_below,c);
            pixels_above = fakeimage(range_above,c);
            selected_pixels_below = pixels_below(pixels_below > rice);
            selected_pixels_above = pixels_above(pixels_above > rice);
            
            % Combine selected pixels
            selected_pixels = [selected_pixels_below; selected_pixels_above];
            
            % Calculate mean of selected pixels, if there are any
            if ~isempty(selected_pixels)
                newImage(r,c) = mean(selected_pixels);
            end
        end
    end
    end

end
    fakeimage=newImage;


    
mmm= highest_points;

rice=0.5;
for c = 1:Myanswer
    m=round(mmm(c));
    n=find(fakeimage(m:row,c)>=rice,10,'first');
    n=n+m-1;
    if length(n) >= 10
    for r = n(5):row
        if fakeimage(r,c) < rice
            % Get index range of 5 pixels above and below, avoiding indices out of range
            range_below = max(1,r-5):r-1;
            range_above = r+1:min(row,r+5);
            
            % Select only pixels with value greater than 0.5
            pixels_below = fakeimage(range_below,c);
            pixels_above = fakeimage(range_above,c);
            selected_pixels_below = pixels_below(pixels_below > rice);
            selected_pixels_above = pixels_above(pixels_above > rice);
            
            % Combine selected pixels
            selected_pixels = [selected_pixels_below; selected_pixels_above];
            
            % Calculate mean of selected pixels, if there are any
            if ~isempty(selected_pixels)
                newImage(r,c) = mean(selected_pixels);
            end
        end
    end
    end

end
    fakeimage=newImage;


mmm=fillmissing(highest_points2,'linear');
highest_points2=mmm;


rice=0.1;



for c = 1:Myanswer2
    m=round(mmm(c));
    n=find(fakeimage(m:row,c)>=rice,10,'first');
    n=n+m-1;
    if length(n) >= 10
    for r = n(10):row
        if fakeimage(r,c) < rice
            % Get index range of 5 pixels above and below, avoiding indices out of range
            range_below = max(1,r-5):r-1;
            range_above = r+1:min(row,r+5);
            
            % Select only pixels with value greater than 0.5
            pixels_below = fakeimage(range_below,c);
            pixels_above = fakeimage(range_above,c);
            selected_pixels_below = pixels_below(pixels_below > rice);
            selected_pixels_above = pixels_above(pixels_above > rice);
            
            % Combine selected pixels
            selected_pixels = [selected_pixels_below; selected_pixels_above];
            
            % Calculate mean of selected pixels, if there are any
            if ~isempty(selected_pixels)
                newImage(r,c) = mean(selected_pixels);
            end
        end
    end
    end

end
    fakeimage=newImage;


        cd(exp_dir3_path{exp})
%         saveas(gcf, output_image);
%         close(gcf)

        saved_myanswer{exp}(lll,1)=Myanswer;
        saved_myanswer{exp}(lll,2)=Myanswer2;


        saved_highestpoints_bulk{exp}(lll,:,:)= highest_points(1:column);
        saved_highestpoints_dispersed{exp}(lll,:,:)= highest_points2(1:column);
        fakeimages{lll}=fakeimage;
%         saved_myanswer(lll,3)=mycol;
%         saved_myanswer(lll,4)=mycol2;


    end
        cd(exp_dir3_path{exp})
        save('highestpoints.mat', 'saved_equation_bulk', 'saved_equation_dispersed','saved_highestpoints_bulk','saved_highestpoints_dispersed', 'saved_myanswer','fakeimages','-v7.3');
        clear fakeimages
    end
end

