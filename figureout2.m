clc
clear
hold on;

option=2;
option2=1;

main_path = pwd;

cd(main_path)

time= [10 21 29];
distance = [2 3 4];

density = 1010; %or 1020

x_f_chose=4;

if density == 1020
    time= [8 15 22];
end

%Non_dimensionalize

xf= 0.15; %m
phi=0.38;
mu= 10^(-3);
D = 0.003;
density_not=998;

k=((D^2)*(phi^3))/(180*((1-phi)^2));
rho= mu/density_not;

gs=(9.81*(density(1)-density_not(1)))/density_not;
Qs= 10^(-6);
t= 0.0762;
qs=Qs/t;

T=((((xf)^3)*phi*rho)/(k*gs*qs))^0.5;

H = ((xf*qs*phi*rho)/(k*gs))^(0.5);

% Row titles
experimentNames = {'0', '2', '4', '6', '8', '10'};


% Set thresholds for different intensity levels
low_threshold = 0.1;
mid_threshold = 0.25;
mid_threshold2 = 0.5;
high_threshold = 0.75;
% high_threshold2 =0.6;
% high_threshold3 =0.7;
% high_threshold4 =0.8;
% high_threshold5 =0.9;

% subplot(length(Exp_dir),(length(something)),n); % 2 rows, 4 columns, 1st plot



n=0;

%User Input - Add or remove the folder for experiments needed
Exp_dir{1}= strcat(main_path, '/Experiment_0holes_5deg_new/Experiment_1010_0elev');
Exp_dir{2}= strcat(main_path, '/Experiment_2holes_5deg_new/Experiment_1010_0elev');
Exp_dir{3}= strcat(main_path, '/Experiment_4holes_5deg/Experiment_1010_0elev');
Exp_dir{4}= strcat(main_path, '/Experiment_6holes_5deg_new/Experiment_1010_0elev');
Exp_dir{5}= strcat(main_path, '/Experiment_8holes_5deg/Experiment_1010_0elev');
Exp_dir{6}= strcat(main_path, '/Experiment_10holes_5deg/Experiment_1010_0elev');


x_f(1)=825;
x_f(2)=654;
x_f(3)=687;
x_f(4)=687;
x_f(5)=697;
x_f(6)=687;


corr(1)=1;
corr(2)=1;
corr(3)=2;
corr(4)=4;
corr(5)=1;
corr(6)=2;

exp_dir_path{1} =strcat(Exp_dir{1},'/new_directory');
exp_dir_path{2} =strcat(Exp_dir{2},'/new_directory');
exp_dir_path{3} =strcat(Exp_dir{3},'/new_directory');
exp_dir_path{4} =strcat(Exp_dir{4},'/new_directory');
exp_dir_path{5} =strcat(Exp_dir{5},'/new_directory');
exp_dir_path{6} =strcat(Exp_dir{6},'/new_directory');

exp_dir2_path{1} =strcat(exp_dir_path{1},'/new_directory');
exp_dir2_path{2} =strcat(exp_dir_path{2},'/new_directory');
exp_dir2_path{3} =strcat(exp_dir_path{3},'/new_directory');
exp_dir2_path{4} =strcat(exp_dir_path{4},'/new_directory');
exp_dir2_path{5} =strcat(exp_dir_path{5},'/new_directory');
exp_dir2_path{6} =strcat(exp_dir_path{6},'/new_directory');


exp_dir3_path{1} =strcat(exp_dir2_path{1},'/new_directory');
exp_dir3_path{2} =strcat(exp_dir2_path{2},'/new_directory');
exp_dir3_path{3} =strcat(exp_dir2_path{3},'/new_directory');
exp_dir3_path{4} =strcat(exp_dir2_path{4},'/new_directory');
exp_dir3_path{5} =strcat(exp_dir2_path{5},'/new_directory');
exp_dir3_path{6} =strcat(exp_dir2_path{6},'/new_directory');

if option==1


for exp=1:length(Exp_dir)
    cd(exp_dir3_path{exp})
    load("highestpoints.mat")
for image=1:length(time)

n=n+1;

% Column titles (time stamps)
timeStamps = {num2str(round((time(1)*300)/(T))/10), num2str(round((time(2)*300)/(T))/10), num2str(round((time(3)*300)/(T))/10)};

hSub=subplot(length(Exp_dir),(length(time)),n);
I3=fakeimages{time(image)+corr(exp)};
Size=size(I3);
crop(n,1)= 4*(x_f(exp));
crop(n,2)= round(0.8*x_f(exp));
% I = flipud(I);

if isempty(I3)
    continue
end
I4= I3(end-crop(n,2)+1:end,1:crop(n,1));

I=imresize(I4, [(600/1.3545)*(0.15/H) 3000]);
I2=I;
% Assign each pixel an intensity level based on these thresholds
% I(I2<=low_threshold) = 1; % Black
% I(I2>low_threshold & I<=mid_threshold) = 2; % Blue
% I(I2>mid_threshold & I<=mid_threshold2) = 3;
% I(I2>mid_threshold2 & I<=high_threshold) = 4; % Yellow
% I(I2>high_threshold & I<=high_threshold2) = 5; %Red
% I(I2>high_threshold2) = 6; % White


I(I2<=low_threshold) = 1; % Black
I(I2>low_threshold & I<=mid_threshold) = 2; % Blue
I(I2>mid_threshold & I<=mid_threshold2) = 3;
I(I2>mid_threshold2 & I<=high_threshold) = 4; % Yellow
I(I2>high_threshold) = 5; % White

% Create a new colormap for these intensity levels
% cmap = [0 0 0; 0 0 1;0 1 0; 1 1 0; 1 0 0; 1 1 1]; % RGB values for black, blue,green, yellow, red, white

% cmap = [0 0 0; 0 0 0.8; 0 0.8 0.8; 0 0.8 0;0.8 0.8 0; 0.8 0 0];
cmap = [0 0 0; 0 0 0.8;0.8 0.8 0;0.8 0 0; 1 1 1];
colormap(cmap);

% Display the image
imagesc(I);



% Generate the corresponding new X-axis scale (in mm)
new_X_scale = 0.5:0.5:4.5;

% Generate the corresponding new X-axis scale (in mm)
new_Y_scale = fliplr(0:round((0.2/1.3545)*(0.15/H)):round((0.8/1.3545)*(0.15/H)));


% Apply new scale
set(gca, 'XTick', round(0.5*750-(1/3)*750):round(0.5*750):round(5*750-(1/3)*750));


if exp<length(Exp_dir)
    set(gca, 'XTickLabel', '','FontSize',1);
else
    set(gca, 'XTickLabel', new_X_scale,'FontSize',18);
end



% Apply new scale
set(gca, 'YTick', 0:round((0.2/1.3545)*750*(0.15/H)):round(750*(0.8/1.3545)*(0.15/H)));

if image>1
   set(gca, 'YTickLabel','','FontSize',1);
else
   set(gca, 'YTickLabel', new_Y_scale);
end



if exp>length(Exp_dir)-1
% Add labels and title
xlabel('$x^{*}$', 'Interpreter', 'latex');
end

if image == 1
ylabel('$z^{*}$', 'Interpreter', 'latex','Rotation',0,'HorizontalAlignment','center','VerticalAlignment','middle');


hYLabel = ylabel('$z^{*}$', 'Interpreter', 'latex', 'Rotation', 0, ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
             
% Get the current position of the ylabel
currentPos = get(hYLabel, 'Position');

if exp==6
    newPos = currentPos + [170 0 0];
else
% Adjust the x-coordinate of the position (change -1 to any desired offset value)
newPos = currentPos + [-250 0 0];
end

% Set the new position
set(hYLabel, 'Position', newPos);

clear newPos

end

% Crop the window to desired size
xlim([0 round(4*750)]);

ylim([0 round((0.8/1.3545)*750*(0.15/H))])

    % Add column titles for the topmost plots
    if n <= 3
        title(timeStamps{n});
    end

    if mod(n,3) == 1
        axPos = get(hSub, 'position'); 
        text('Units', 'normalized', 'Position', [-0.2, 0.5], ...
             'String', experimentNames{ceil(n/3)}, ...
             'Rotation', 0, ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', ...
             'Parent', hSub,'FontSize',18);
    end

box on
set(gca,'fontsize',18)


end
end

% Adding a Y label for the entire figure using text
annotation('textbox',[0.05, 0.5, 0, 0], 'String', 'Number of Holes', ...
            'HorizontalAlignment', 'center', ...
           'VerticalAlignment', 'middle', 'EdgeColor', 'none','FontSize',18);

% Adding a X label for the entire figure using text
annotation('textbox',[0.5, 0.98, 0, 0], 'String', '$t^{*}$', 'Interpreter', 'latex', ...
            'HorizontalAlignment', 'center', ...
           'VerticalAlignment', 'middle', 'EdgeColor', 'none','FontSize',18);


end

if option==2

for exp=1:length(Exp_dir)
    cd(exp_dir3_path{exp})
    load("highestpoints.mat")

    cd(Exp_dir{exp})
    load("Buoyancy.mat")

experimentNames = {'0 Holes', '2 Holes', '4 Holes', '6 Holes', '8 Holes', '10 Holes'};

timeStamps = {num2str(x_f_chose)};


    idx=find(saved_myans(:,2) >= 4,1,'first');
    Z(exp)= idx;
hSub=subplot(length(Exp_dir),1,exp);

%Correction

if exp==1
    Z(exp)=Z(exp)+1;
elseif exp==2
    Z(exp)=Z(exp)+0;
elseif exp==3
    Z(exp)=Z(exp)+1;
elseif exp==4
    Z(exp)=Z(exp)+1;
elseif exp==5
    Z(exp)=Z(exp)+1;
elseif exp==6
    Z(exp)=Z(exp)+3;
end

I3=fakeimages{Z(exp)};
Size=size(I3);
crop(exp,1)= 4*(x_f(exp));
crop(exp,2)= round(0.8*x_f(exp));
% I = flipud(I);

if isempty(I3)
    continue
end
I4= I3(end-crop(exp,2)+1:end,1:crop(exp,1));

I=imresize(I4, [600*(0.15/H) 3000]);
I2=I;
% Assign each pixel an intensity level based on these thresholds
I(I2<=low_threshold) = 1; % Black
I(I2>low_threshold & I<=mid_threshold) = 2; % Blue
I(I2>mid_threshold & I<=mid_threshold2) = 3;
I(I2>mid_threshold2 & I<=high_threshold) = 4; % Yellow
I(I2>high_threshold) = 5; % White

% Create a new colormap for these intensity levels
% cmap = [0 0 0; 0 0 1;0 1 0; 1 1 0; 1 0 0; 1 1 1]; % RGB values for black, blue,green, yellow, red, white

% cmap = [0 0 0; 0 0 0.8; 0 0.8 0.8; 0 0.8 0;0.8 0.8 0; 0.8 0 0];
cmap = [0 0 0; 0 0 0.8;0.8 0.8 0;0.8 0 0; 1 1 1];
% Create a new colormap for these intensity levels
% cmap = [0 0 0; 0 0 1; 1 1 0; 1 0 0]; % RGB values for black, blue, yellow, red
colormap(cmap);

% Display the image
imagesc(I);



% Generate the corresponding new X-axis scale (in mm)
new_X_scale = 0.5:0.5:4.5;

% Generate the corresponding new X-axis scale (in mm)
new_Y_scale = fliplr(0:round(0.2*(0.15/H)):round(0.8*(0.15/H)));


% Apply new scale
set(gca, 'XTick', round(0.5*750-(1/3)*750):round(0.5*750):round(5*750-(1/3)*750));


if exp<length(Exp_dir)
    set(gca, 'XTickLabel', '','FontSize',1);
else
    set(gca, 'XTickLabel', new_X_scale,'FontSize',18);
end


% Apply new scale
set(gca, 'YTick', 0:round(0.2*750*(0.15/H)):round(750*0.6*(0.15/H)));
set(gca, 'YTickLabel', new_Y_scale);


% Add labels and title
if exp==length(Exp_dir)
xlabel('$x^{*}$', 'Interpreter', 'latex');
end
ylabel('$z^{*}$', 'Interpreter', 'latex','Rotation',0,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',18);


hYLabel = ylabel([experimentNames{exp}, '             $z^{*}$'], 'Interpreter', 'latex', 'Rotation', 0, ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle','FontSize',18);
             
% Get the current position of the ylabel
currentPos = get(hYLabel, 'Position');

if exp<length(Exp_dir)
   newPos = currentPos + [-175 0 0];
else
   newPos = currentPos + [-75 0 0];
end
    

box on

% Set the new position
set(hYLabel, 'Position', newPos);
% Crop the window to desired size
xlim([0 round(4*750)]);

ylim([0 round(0.8*750*(0.15/H))])

%     % Add column titles for the topmost plots
%     if exp <= 1
%         title(timeStamps{1});
%     end

%     if mod(n,1) == 1
%         axPos = get(hSub, 'position'); 
%         text('Units', 'normalized', 'Position', [-0.2, 0.5], ...
%              'String', experimentNames{ceil(exp)}, ...
%              'Rotation', 0, ...
%              'HorizontalAlignment', 'center', ...
%              'VerticalAlignment', 'middle', ...
%              'Parent', hSub,'FontSize',18);
%     end

box on
set(gca,'fontsize',18)


end



% % Adding a Y label for the entire figure using text
% annotation('textbox',[0.05, 0.5, 0, 0], 'String', 'Number of Holes', ...
%             'HorizontalAlignment', 'center', ...
%            'VerticalAlignment', 'middle', 'EdgeColor', 'none','FontSize',18);

% % Adding a X label for the entire figure using text
% annotation('textbox',[0.5, 0.98, 0, 0], 'String', '$t^{*}$', 'Interpreter', 'latex', ...
%             'HorizontalAlignment', 'center', ...
%            'VerticalAlignment', 'middle', 'EdgeColor', 'none','FontSize',18);
% 

end

