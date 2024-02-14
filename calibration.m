function [] = calibration(density_range,power)

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


 cd(new_dir_path2)    
    imgs = cell(num_images, 1);

 %Each of the calibration images are read through and saved

    for i = 1:num_images
        Variable = ['IMG',num2str(i),'_cropped.jpg'];
        image = imread(Variable);
        image = im2double(image);
        imgs{i} = image;
    end

% Initialize arrays to store the pixel intensities and corresponding

    cellArray=imgs;

% Get dimensions of data in cells
[row_dim, col_dim] = size(cellArray{1});

% Initialize a matrix to keep track of whether the points at each location are increasing monotonically

non_monotonic_points = false(row_dim, col_dim);


counter=0;

% Check each location for monotonicity
for i = 1:row_dim
    for j = 1:col_dim

        clear point_values 
        point_values = cellfun(@(x) x(i, j), cellArray);

        %Get the difference between each subsequent calibration images
        non_monotonic_diff = diff(point_values);
        non_increasing_index = find(non_monotonic_diff <= 0); % '+1' because 'diff' reduces the index by 1
        if any(non_monotonic_diff <= 0)
            num_non_increasing_points = sum(non_monotonic_diff <= 0);
            if num_non_increasing_points > 1 %More than 1 
                non_monotonic_points(i, j) = true;
                counter=counter+1;
            else
                % Find the non-increasing point and assign NaN only to it
                idx1 = non_increasing_index;
                idx2 = non_increasing_index + 1;
            
                pixel_vals1 = point_values;
                pixel_vals1(idx1) = [];
                pixel_vals2 = point_values;
                pixel_vals2(idx2) = [];
            
                is_mono1 = all(diff(pixel_vals1) > 0);
            
                is_mono2 = all(diff(pixel_vals2) > 0);

                if is_mono1 && is_mono2
                    % If removing either point results in a monotonically increasing sequence, remove the point that is less similar to its neighbors
                    val1 = point_values(idx1);
                    val2 = point_values(idx2);

                    if idx2+1 <= num_images && idx1-1 ~= 0

                        slope = abs(point_values(idx2+1)-point_values(idx1-1))/(density_range(idx2+1)-density_range(idx1-1));

                        diff1 = abs(point_values(idx1)-(slope*density_range(idx1)));
                        diff2 = abs(point_values(idx2)-(slope*density_range(idx2)));

                    if diff1 < diff2
                        cellArray{idx2}(i,j) = NaN;
                    else
                        cellArray{idx1}(i,j) = NaN;
                    end
               
                    elseif idx1-1 == 0
                        diff1 = abs(val1 - point_values(idx2+1));
                        diff2 = abs(val2 - point_values(idx2+1));
                        if diff1 <= diff2
                            cellArray{idx2}(i,j) = NaN;
                        else
                            cellArray{idx1}(i,j) = NaN;
                        end
                    else
                        diff1 = abs(val1 - point_values(idx1-1)); 
                        diff2 = abs(val2 - point_values(idx1-1));
                        
                        if diff1 < diff2
                            cellArray{idx2}(i,j) = NaN;
                        else
                            cellArray{idx1}(i,j) = NaN;
                        end
                    
                    end

                elseif is_mono1
                    cellArray{idx1}(i,j) = NaN;
                elseif is_mono2
                    cellArray{idx2}(i,j) = NaN;
                else
                    % If removing neither point results in a monotonically increasing sequence, set all values to NaN
                    for k = 1:num_images
                        cellArray{k}(i,j) = NaN;
                        counter=counter+1;
                    end
                end
            end
        end
    end
end

% Replace non-monotonic points with NaN in each cell
for k = 1:numel(cellArray)
    cellArray{k}(non_monotonic_points) = NaN;
end

% % Check each location for monotonicity. This is done a second time 
% for i = 1:row_dim
%     for j = 1:col_dim
%         point_values = cellfun(@(x) x(i, j), cellArray);
%         if sum(isnan(point_values))>1
%             continue
%         end
%         n = numel(point_values);
%         non_monotonic_diff = nan(n-1, 1);  % Initialize output
%  
%         for k = 2:n
%             if isnan(point_values(k))
%                 continue;  % Skip if either of the current or previous element is NaN
%             elseif isnan(point_values(k-1))
%                 non_monotonic_diff(k-1)=point_values(k)- point_values(k-2);
%             else
%                 non_monotonic_diff(k-1) = point_values(k) - point_values(k-1);  % Compute difference
%             end
%         end
%         if any(non_monotonic_diff <= 0)
%             non_monotonic_points(i, j) = true;
%         end
%     end
% end
% 
% for k = 1:numel(cellArray)
%     cellArray{k}(non_monotonic_points) = NaN;
% end

%For user check. 
% Check percentage of total points that are NaN
totalPoints = row_dim * col_dim * numel(cellArray);
numNaN = sum(cellfun(@(x) sum(isnan(x(:))), cellArray));
percentNaN = (numNaN / totalPoints) * 100;  % Convert to percentage
percentlost = (counter/(row_dim*col_dim))*100;

disp(['For calibration, the total NaN values are ' num2str(percentNaN) '%.']);
disp(['For calibration, the total lost curves are ' num2str(percentlost) '%.']);


% for k = 1:numel(cellArray)
%     
%     nan_mask = ~isnan(cellArray{k});
%     % Display the mask as an image
%     figure;
%     imshow(nan_mask);
%     title(['NaN Locations for Image ', num2str(k)]);
% end

imgs= cellArray;

[row_dim, col_dim] = size(imgs{1});

% Flatten images
flattened_imgs = cellfun(@(x) x(:), imgs, 'UniformOutput', false);

    intensities = nan(numel(imgs{1}),num_images);

    % Loop over all pixels
for i = 1:row_dim*col_dim
    % Get the pixel values from all images
    pixel_vals = cellfun(@(x) x(i), flattened_imgs);
    
    % Check for NaN values
    not_nan_inds = ~isnan(pixel_vals);
    
    % Only proceed if there are more than 'power' number of non-NaN values
    if sum(not_nan_inds) > length(density_range)-2
        % Fit a polynomial curve to the pixel values vs. density
        intensities(i,:)=pixel_vals;
    end

end

% Reshape P to have the same dimensions as the original images
Intensities = reshape(intensities, [row_dim, col_dim, length(density_range)]);

    % Save the intensity/density arrays to a file
    save('calibration_data.mat', 'intensities', 'Intensities');
end
