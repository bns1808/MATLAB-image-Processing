clc
clear
close all
hold on;

option = 1;

ess=1;

main_path = pwd;

cd(main_path)


%User Input - Specify density/nozzle height for plotting

density(1)= 1010;%% Use the two densities 1010 or 1020
density(2)=1020; %% Use the two densities 1010 or 1020
elev(1) = 0; %%Use the two heights 0 or 10 for elevation (in cm)
elev(2) = 0;

if option == 1

%User Input - Add or remove the folder for experiments needed

Exp_dir{1}= strcat(main_path, '/Experiment_0holes_5deg_new');
Exp_dir{2}= strcat(main_path, '/Experiment_2holes_5deg_new');
Exp_dir{3}= strcat(main_path, '/Experiment_4holes_5deg');
Exp_dir{4}= strcat(main_path, '/Experiment_6holes_5deg_new');
Exp_dir{5}= strcat(main_path, '/Experiment_8holes_5deg');
Exp_dir{6}= strcat(main_path, '/Experiment_10holes_5deg');



x_f(1)=775;
x_f(2)=682;
x_f(3)=687;
x_f(4)=698;
x_f(5)=697;
x_f(6)=688;




%User input for number of Holes per experiment

H(1) = 0;
H(2) = 2;
H(3) = 4;
if density(1)==1020
    H(4)=5;
else
H(4) = 6;
end
H(5) = 8;
H(6) = 10;




elseif option==2

%User Input - Add or remove the folder for experiments needed

Exp_dir{1}= strcat(main_path, '/Experiment_2holes_-5deg_new');
Exp_dir{2}= strcat(main_path, '/Experiment_2holes_-2.5deg_new');
Exp_dir{3}= strcat(main_path, '/Experiment_2holes_0deg');
Exp_dir{3}= strcat(main_path, '/Experiment_2holes_2.5deg');
Exp_dir{4}= strcat(main_path, '/Experiment_2holes_5deg_new');


%User Input - Specify density/nozzle height for plotting



x_f(1)=813;
x_f(2)=787;
% x_f(3)=835;
x_f(3)=658;
x_f(4)=682;



%User input for number of Holes per experiment

H(1) = -5;
H(2) = -2.5;
H(3) = 0;
H(3) = 2.5;
H(4) = 5;
end






%User input - threshold - Choose how many data points to have

threshold = [2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4.0];

if option==2
    threshold= [2.4 2.6 2.8 3 3.2 3.4 3.6 3.8 4.0 4.2];
end

starting = 5;
ending = 25;
max_ending = 45;

imagenumber=0;

Y = threshold;
X=H;

all1=4;

for exp=1:length(Exp_dir)
    x_f_mm(exp)=150/x_f(exp);
end

for exp=1:length(Exp_dir)
    Variable=['/Experiment_', num2str(density(1)),'_',num2str(elev(1)),'elev'];
    Exp2_dir{exp}= strcat(Exp_dir{exp}, Variable);
end

if ess==2
    for exp=1:length(Exp_dir)
    Variable=['/Experiment_', num2str(density(2)),'_',num2str(elev(2)),'elev'];
    Exp3_dir{exp}= strcat(Exp_dir{exp}, Variable);
    end
end



for exp=1:length(Exp_dir)
    cd(Exp2_dir{exp})
    load('Buoyancy.mat')
    y{exp} = saved_myans(:,2);
    y{exp} = round((y{exp}/(threshold(2)-threshold(1))))*(threshold(2)-threshold(1));
    dispersion{exp}=D(starting:ending);
    notdispersions{exp}=B;
    Side{exp}= Percent(starting:ending);
    Sides{exp}=Percent;
    dispersions{exp}=D;
    buoyancies{exp}=D+B;
    trying_new{exp}= diff(D)./diff(buoyancies{exp});
    P_Holder=abs((Percent(ending)-Percent(starting)))/28;
    speed{exp}=saved_myans(:,1);
    speed_disp{exp}=saved_myans(:,2);

    for i=2:ending-starting+1
        prev_non_nan_index = find(~isnan(Side{exp}(1:i-1)), 1, 'last'); % finds the last non-NaN index before the specified index
        prev_non_nan_value = Side{exp}(prev_non_nan_index);
        if Side{exp}(i)+0.1 <  prev_non_nan_value
            Side{exp}(i)=NaN;
            dispersion{exp}(i)=NaN;
        end
    end
Side{exp} = fillmissing(Side{exp},'linear');
dispersion{exp}=fillmissing(dispersion{exp},'linear');

end
if ess==2
for exp=1:length(Exp_dir)
    cd(Exp3_dir{exp})
    load('Buoyancy.mat')
    y2{exp} = saved_myans(:,2);
    y2{exp} = round((y2{exp}/(threshold(2)-threshold(1))))*(threshold(2)-threshold(1));
    dispersion2{exp}=D(starting:ending);
    Side2{exp}= Percent(starting:ending);
    Sides2{exp}=Percent;
    dispersions2{exp}=D;
    buoyancies2{exp}=D+B;
    P_Holder=abs((Percent(ending)-Percent(starting)))/28;

    speed2{exp}=saved_myans(:,1);
    speed_disp2{exp}=saved_myans(:,2);

    for i=2:ending-starting+1
        prev_non_nan_index = find(~isnan(Side2{exp}(1:i-1)), 1, 'last'); % finds the last non-NaN index before the specified index
        prev_non_nan_value = Side2{exp}(prev_non_nan_index);
        if Side2{exp}(i)+0.1 <  prev_non_nan_value
            Side2{exp}(i)=NaN;
            dispersion2{exp}(i)=NaN;
        end
    end
Side2{exp} = fillmissing(Side2{exp},'linear');
dispersion2{exp}=fillmissing(dispersion2{exp},'linear');

end
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

gs2=(9.81*(density(2)-density_not))/density_not;
Qs= 10^(-6);
t= 0.0762;
qs=Qs/t;

T=((((xf)^3)*phi*rho)/(k*gs*qs))^0.5;
if ess==2
T2=((((xf)^3)*phi*rho)/(k*gs2*qs))^0.5;
end

for exp=1:length(Exp_dir)
    T_fix(exp) = find(dispersions{exp}>0,1,'last');
    if ess==2
        T_fix2(exp) = find(dispersions2{exp}>0,1,'last');
    end
end

T_fix= T_fix*30;
T_fix= T_fix/T;

At(1)=round(((density(1)-density_not)/density_not)*100)/100;
At(2)=round(((density(2)-density_not)/density_not)*100)/100;

if ess==2
T_fix2 = T_fix2*30;
T_fix2= T_fix2/T2;
end

   figure(imagenumber+1)
box on
set(gca,'fontsize',18)
for exp=1:length(Exp_dir)
 
    endings=max(ending,length(Sides{exp})-2);
    endings = min(endings,max_ending);
    hold on;plot(((starting:endings)*30)/T,Sides{exp}(starting:endings),'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3)
end

subplot(2,1,1)

imagenumber=imagenumber+1;
xlabel('$t^{*}$', 'Interpreter', 'latex');
ylabel('$\overline{B}^{*}_{disp}$', 'Interpreter', 'latex');

if option==1
    legend('0', '2', '4', '6', '8', '10','location','northwest');
elseif option==2
    legend(['-5', char(176)], ['-2.5', char(176)], ['2.5', char(176)], ['5', char(176)],'location','northwest');
end

xlim([(60*(5/2))/T ((max_ending/2)*60)/T])
maxValue = max(cellfun(@(x) max(x(:)), Sides));
ylim([0.1 maxValue*1.1])
box on
set(gca,'fontsize',18)
hold off


for exp=1:length(Exp_dir)
    figure(imagenumber+1)
    endings=max(ending,length(dispersions{exp})-2);
    ending2 = find(dispersions{exp}>0,1,'last');
    endings = min(endings,max_ending);
    endings = min(endings,ending2);
    hold on;plot(((starting:endings)*30)/T,T_fix(exp)*dispersions{exp}(starting:endings),'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3)
end
imagenumber=imagenumber+1;
box on
set(gca,'fontsize',18)

xlabel('$t^{*}$', 'Interpreter', 'latex');
ylabel('$B^{*}_{disp}$', 'Interpreter', 'latex');

for jk=1:length(Exp_dir)
    jkk=find(~isnan(trying_new{exp}),1,"first");
    for jk=5:jkk
    if trying_new{exp}(jkk)==0 || isempty(trying_new{exp}(jkk))
        trying_new= max(trying_new{exp}(jk-1),trying_new{exp}(jk+1));
    end
    end
end

if option==1
    legend('0', '2', '4', '6', '8', '10','location','northwest');
elseif option==2
    legend(['-5', char(176)], ['-2.5', char(176)], ['2.5', char(176)], ['5', char(176)],'location','northwest');
end

xlim([(60*(5/2))/T ((max_ending/2)*60)/T])
% maxValue = max(cellfun(@(x) max(x(:)), dispersions));
ylim([0.03 0.3])

box on
set(gca,'fontsize',18)
hold off
    figure(imagenumber+1)
    box on
set(gca,'fontsize',18)

for exp=1:length(Exp_dir)

    hold on;plot((starting:ending)/2,dispersion{exp})
end
imagenumber=imagenumber+1;
box on
set(gca,'fontsize',18)
xlabel('$t^{*}$', 'Interpreter', 'latex');
ylabel('$B^{*}$ Dispersed', 'Interpreter', 'latex');

if option==1
    legend('0', '2', '4', '6', '8', '10','location','northwest');
elseif option==2
    legend(['-5', char(176)], ['-2.5', char(176)], ['2.5', char(176)], ['5', char(176)],'location','northwest');
end

h = (elev(1)/100)/(((xf*qs*phi*rho)/(k*gs))^(0.5));
h2= (elev(2)/100)/(((xf*qs*phi*rho)/(k*gs2))^(0.5));

xlim([5/2 (max_ending/2)+1])
maxValue = max(cellfun(@(x) max(x(:)), dispersion));
ylim([0 maxValue*1.1])

hold off


for exp=1:length(Exp_dir)
for i=1:length(threshold)
    idx=find(y{exp} >= threshold(i),1,'first');
    Z(i,exp)= Sides{exp}(idx);
    S(i,exp) = buoyancies{exp}(idx);
    if ess==2
    idx=find(y2{exp} >= threshold(i),1,'first');    
    Z2_2(i,exp)=Sides2{exp}(idx);
    S2(i,exp) = buoyancies2{exp}(idx);
    end
    TT(i,exp)=trying_new{exp}(idx);
end
S(:,exp)=S(:,exp)*T_fix(exp);
if ess==2
S2(:,exp)=S2(:,exp)*T_fix2(exp);
end
end
surf(X,Y,Z)

box on
set(gca,'fontsize',18)

if option==1
xlabel('$n_{holes}$', 'Interpreter', 'latex','FontSize',15);
elseif option==2
    xlabel('$\theta$', 'Interpreter', 'latex');
end
ylabel('$x^{*}$','Interpreter','latex','FontSize',15);
zlabel('$\overline{B}^{*}_{disp}$', 'Interpreter', 'latex','FontSize',15);


box on
set(gca,'fontsize',18)

figure(imagenumber+1)
box on
set(gca,'fontsize',18)
for exp=1:length(Exp_dir)
    endings=max(ending,length(Sides{exp})-2);
    endings = min(endings,max_ending);
    ending2 = find(speed{exp}>0.5,1,'last');
    endings = min(endings,max_ending);
    endings = min(endings,ending2);
    hold on;plot(speed{exp}(starting:endings),((starting:endings)*30)/T,'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3)
end

if option==1
    legend('0', '2', '4', '6', '8', '10','location','northwest');
elseif option==2
    legend(['-5', char(176)], ['-2.5', char(176)], ['2.5', char(176)], ['5', char(176)],'location','northwest');
end

ylim([(60*(5/2))/T 4.5])
% maxValue = max(cellfun(@(x) max(x(:)), speed_disp));
% xlim([0 maxValue*1.1])
xlim([1 5])


    xlabel('$x^{*}_{nose, Bulk}$', 'Interpreter', 'latex');

ylabel('$t^{*}$', 'Interpreter', 'latex');


imagenumber=imagenumber+1;
if option==1
    legend('0', '2', '4', '6', '8', '10','location','northwest');
elseif option==2
    legend(['-5', char(176)], ['-2.5', char(176)], ['2.5', char(176)], ['5', char(176)],'location','northwest');
end


figure (imagenumber+1)
box on
set(gca,'fontsize',18)
imagenumber=imagenumber+1;

Y2=starting:ending;
Y2=(Y2*30)/T;

for exp=1:length(Exp_dir)
for i=1:(ending-starting)+1
    Z2(i,exp)= Side{exp}(i);
end
end



surf(X,Y2,Z2)
if option==1
xlabel('$n_{holes}$', 'Interpreter', 'latex','FontSize',15);
elseif option==2
        xlabel('$\theta$', 'Interpreter', 'latex');
end
ylabel('$t^{*}$', 'Interpreter', 'latex','FontSize',15);
zlabel('$\overline{B}^{*}_{disp}$', 'Interpreter', 'latex','FontSize',15);
box on
set(gca,'fontsize',18)

figure (imagenumber+1)
box on
set(gca,'fontsize',18)
imagenumber=imagenumber+1;


S=(((S).^2))*300;

if ess==2
    S2=(((S2).^2))*300;
end

f=9;

plot(H,Z(f,:),'LineWidth', 2,'Color','red')
if ess==2
    H2=H;
    if density(2)==1020
        idx=find(H==6,1,"first");
        H2(idx)=5;
    end
hold on;plot(H,Z2_2(f,:),'LineWidth',2,'Color','blue')
end


hold on
box on
set(gca,'fontsize',18)
scatter(H(1:end),Z(f,1:end),S(f,1:end),'filled','MarkerfaceColor','red')

if ess==2
    scatter(H(1:end),Z2_2(f,1:end),S2(f,1:end),'filled','MarkerFaceColor','blue')
end


if option==2
xlim([-5 max(H)+3]);
ylim([0 0.7]);
elseif option==1
xlim([0 max(H)+3]);

if ess==2
ylim([0 0.7]); 
else
    ylim([0 0.7])
end
end
if option==1 
xlabel('$n_{holes}$', 'Interpreter', 'latex');
elseif option==2
    xlabel('$\theta$', 'Interpreter', 'latex');
end

ylabel('$\overline{B}^{*}_{disp}$', 'Interpreter', 'latex');
% legend('x_{f} = 4');

if ess==2
    if max(elev)>min(elev)
        legend(['elev = ',num2str(round(h*10)/10)], ['elev = ',num2str(round(h2*10)/10)], 'Interpreter', 'latex');
    else
    legend(['At = ', num2str(At(1)) ],['At = ', num2str(At(2)) ], 'Interpreter', 'latex','location','northeast');
    end
end

figure(imagenumber+1)
box on
set(gca,'fontsize',18)

for exp=1:length(Exp_dir)
    figure(imagenumber+1)
    endings=max(ending,length(Sides{exp})-2);
    endings = min(endings,max_ending);
    ending2 = find(speed_disp{exp}>0.5,1,'last');
    endings = min(endings,max_ending);
    endings = min(endings,ending2);
    hold on;plot(speed_disp{exp}(starting:endings),((starting:endings)*30)/T,'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3)
end

if option==1
    legend('0', '2', '4', '6', '8', '10','location','northwest');
elseif option==2
    legend(['-5', char(176)], ['-2.5', char(176)], ['2.5', char(176)], ['5', char(176)],'location','northwest');
end

ylim([(60*(5/2))/T 4.5])
% maxValue = max(cellfun(@(x) max(x(:)), speed_disp));
% xlim([0 maxValue*1.1])
xlim([1 5])

xlabel('$x^{*}_{nose, dispersed}$', 'Interpreter', 'latex');
ylabel('$t^{*}$', 'Interpreter', 'latex');

box on
set(gca,'fontsize',18)


imagenumber=imagenumber+1;





for exp=1:length(Exp_dir)
    cd(Exp2_dir{exp})
    load('Area.mat')

    
    A_dispersion{exp}=A_D(starting:ending);
    A_notdispersions{exp}=A_B;
    A_Side{exp}= A_Percent(starting:ending);
    A_Sides{exp}=A_Percent;
    A_dispersions{exp}=A_D;
    A_buoyancies{exp}=A_D+A_B;
    P_Holder=abs((A_Percent(ending)-A_Percent(starting)))/28;

    speed{exp}=saved_myans(:,1);
    speed_disp{exp}=saved_myans(:,2);

    for i=2:ending-starting+1
        prev_non_nan_index = find(~isnan(A_Side{exp}(1:i-1)), 1, 'last'); % finds the last non-NaN index before the specified index
        prev_non_nan_value = Side{exp}(prev_non_nan_index);
        if A_Side{exp}(i)+0.1 <  prev_non_nan_value
            A_Side{exp}(i)=NaN;
            A_dispersion{exp}(i)=NaN;
        end
    end
A_Side{exp} = fillmissing(A_Side{exp},'linear');
A_dispersion{exp}=fillmissing(A_dispersion{exp},'linear');

end

%%%%%%%%%%%%%%%%%%%%%%%%\
%%%%
%%


    figure(imagenumber+1)


for exp=1:length(Exp_dir)

    endings=max(ending,length(A_Sides{exp})-2);
    endings = min(endings,max_ending);
    hold on;plot(((starting:endings)*30)/T,A_Sides{exp}(starting:endings),'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3)
end
imagenumber=imagenumber+1;
xlabel('$t^{*}$', 'Interpreter', 'latex');
ylabel('$\overline{A}^{*}_{disp}$', 'Interpreter', 'latex');

if option==1
    legend('0', '2', '4', '6', '8', '10','location','northwest');
elseif option==2
    legend(['-5', char(176)], ['0', char(176)], ['2.5', char(176)], ['5', char(176)],'location','northwest');
end

xlim([(60*(5/2))/T ((max_ending/2)*60)/T])
maxValue = max(cellfun(@(x) max(x(:)), A_Sides));
ylim([0.1 0.9])
box on
set(gca,'fontsize',18)

hold off

figure(imagenumber+1)

% Position the plot at position 1,1 in a 2x2 subplot
subplot(2, 2, 1);
box on
set(gca, 'fontsize', 18)

for exp = 1:length(Exp_dir)
    endings = max(ending, length(A_dispersions{exp}) - 2);
    ending2 = find(A_dispersions{exp} > 0, 1, 'last');
    endings = min(endings, max_ending);
    endings = min(endings, ending2);
    hold on; plot(((starting:endings) * 30) / T, T_fix(exp) * A_dispersions{exp}(starting:endings), 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3)
end

imagenumber=imagenumber+1;
xlabel('$t^{*}$', 'Interpreter', 'latex');
ylabel('$A^{*}_{disp}$', 'Interpreter', 'latex');

if option==1
    legend('0', '2', '4', '6', '8', '10','location','northwest');
elseif option==2
    legend(['-5', char(176)], ['0', char(176)], ['2.5', char(176)], ['5', char(176)],'location','northwest');
end


%  maxValue = max(cellfun(@(x) max(x(:)), dispersions));
ylim([0 1])
xlim([(60*(5/2))/T ((max_ending/2)*60)/T])
box on
set(gca,'fontsize',18)

subplot(2,2,3)


for exp=1:length(Exp_dir)
    endings=max(ending,length(dispersions{exp})-2);
    ending2 = find(dispersions{exp}>0,1,'last');
    endings = min(endings,max_ending);
    endings = min(endings,ending2);
    hold on;plot(((starting:endings)*30)/T,T_fix(exp)*dispersions{exp}(starting:endings),'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3)
end
imagenumber=imagenumber+1;
box on
set(gca,'fontsize',18)

xlabel('$t^{*}$', 'Interpreter', 'latex');
ylabel('$B^{*}_{disp}$', 'Interpreter', 'latex');

if option==1
    legend('0', '2', '4', '6', '8', '10','location','northwest');
elseif option==2
    legend(['-5', char(176)], ['0', char(176)], ['2.5', char(176)], ['5', char(176)],'location','northwest');
end

xlim([(60*(5/2))/T ((max_ending/2)*60)/T])
% maxValue = max(cellfun(@(x) max(x(:)), dispersions));
ylim([0 1])

box on
set(gca,'fontsize',18)
hold off


% Position the plot at position 1,1 in a 2x2 subplot
subplot(2, 2, 2);
box on
set(gca, 'fontsize', 18)

for exp = 1:length(Exp_dir)
    endings = max(ending, length(A_notdispersions{exp}) - 2);
    ending2 = find(A_notdispersions{exp} > 0, 1, 'last');
    endings = min(endings, max_ending);
    endings = min(endings, ending2);
    hold on; plot(((starting:endings) * 30) / T, T_fix(exp) * A_notdispersions{exp}(starting:endings), 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3)
end

imagenumber=imagenumber+1;
xlabel('$t^{*}$', 'Interpreter', 'latex');
ylabel('$A^{*}_{bulk}$', 'Interpreter', 'latex');

if option==1
    legend('0', '2', '4', '6', '8', '10','location','northwest');
elseif option==2
    legend(['-5', char(176)], ['0', char(176)], ['2.5', char(176)], ['5', char(176)],'location','northwest');
end


%  maxValue = max(cellfun(@(x) max(x(:)), dispersions));
ylim([0 1])
xlim([(60*(5/2))/T ((max_ending/2)*60)/T])
box on
set(gca,'fontsize',18)

subplot(2,2,4)

box on
set(gca, 'fontsize', 18)

for exp = 1:length(Exp_dir)
    endings = max(ending, length(notdispersions{exp}) - 2);
    ending2 = find(notdispersions{exp} > 0, 1, 'last');
    endings = min(endings, max_ending);
    endings = min(endings, ending2);
    hold on; plot(((starting:endings) * 30) / T, T_fix(exp) * notdispersions{exp}(starting:endings), 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3)
end

imagenumber=imagenumber+1;
xlabel('$t^{*}$', 'Interpreter', 'latex');
ylabel('$B^{*}_{bulk}$', 'Interpreter', 'latex');

if option==1
    legend('0', '2', '4', '6', '8', '10','location','northwest');
elseif option==2
    legend(['-5', char(176)], ['0', char(176)], ['2.5', char(176)], ['5', char(176)],'location','northwest');
end


%  maxValue = max(cellfun(@(x) max(x(:)), dispersions));
ylim([0 1])
xlim([(60*(5/2))/T ((max_ending/2)*60)/T])
box on
set(gca,'fontsize',18)

imagenumber=imagenumber+1;

figure(imagenumber+1);

for exp = 1:length(Exp_dir)
subplot (3,2,exp)
endings = max(ending, length(notdispersions{exp}) - 2);
    ending2 = find(notdispersions{exp} > 0, 1, 'last');
    endings = min(endings, max_ending);
    endings = min(endings, ending2);
    hold on; plot(((starting:endings) * 30) / T, Sides{exp}(starting:endings), 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3)
    hold on; plot(((starting:endings) * 30) / T, A_Sides{exp}(starting:endings), 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 3)
xlabel('$t^{*}$', 'Interpreter', 'latex');
ylabel('$\overline{B}^{*}_{disp}$, $\overline{A}^{*}_{disp}$ ', 'Interpreter', 'latex');

    legend('Buoyancy', 'Area','location','northwest');

title([num2str(H(exp)),' Holes'])

%  maxValue = max(cellfun(@(x) max(x(:)), dispersions));
ylim([0 1])
xlim([0.6 ((max_ending/2)*60)/T])
box on
set(gca,'fontsize',18)
end


