for image=1:1

% Suppose this is your grayscale image data matrix
exp=;
x_f = 687;

% Set thresholds for different intensity levels
low_threshold = 0.125;
mid_threshold = 0.3;
high_threshold = 0.5;

I=fakeimages{exp};
I2 = I;
% 
% Assign each pixel an intensity level based on these thresholds
I(I2<=low_threshold) = 1; % Black
I(I2>low_threshold & I<=mid_threshold) = 2; % Blue
I(I2>mid_threshold & I<=high_threshold) = 3; % Yellow
I(I2>high_threshold) = 4; % Red

% Create a new colormap for these intensity levels
cmap = [0 0 0; 0 0 1; 1 1 0; 1 0 0]; % RGB values for black, blue, yellow, red
colormap(cmap);

% Display the image
imagesc(I);

% Specify conversion factor from pixels to mm
mm_per_pixel = 1/x_f;

% Generate the corresponding new X-axis scale (in mm)
new_X_scale = 0.5:0.5:5;

% Generate the corresponding new X-axis scale (in mm)
new_Y_scale = fliplr(0:0.5:1);


% Apply new scale
set(gca, 'XTick', round(0.5*x_f-(1/3)*x_f):round(0.5*x_f):round(5*x_f-(1/3)*x_f));
set(gca, 'XTickLabel', new_X_scale);

% Apply new scale
set(gca, 'YTick', 0:round(0.5*x_f):round(x_f));
set(gca, 'YTickLabel', new_Y_scale);

% Add labels and title
xlabel('$X^{*}$', 'Interpreter', 'latex');
ylabel('$Z^{*}$', 'Interpreter', 'latex');

% Crop the window to desired size
xlim([0 round(4.4*x_f)]);

ylim([0 round(x_f)])

pbaspect([5 1 1]);  % Make the width 5 times the height

         % Save the resulting figure as an image file
         output_image = ['Output_',num2str(exp),'.png'];
         saveas(gcf, output_image);
         close(gcf)



end
