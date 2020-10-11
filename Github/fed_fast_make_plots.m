clear all
% by default, we assume a full stomach (fed) and youth (full NAD)
simulation = 1;  % fed (insulin=1, AA=1), no perturbation
aging = [1 1];  % young
% aging = [0.25 0.4];  % NAD 25% of what it used to be; SIRT scaled by 1/3

% no drug
% [amplitude, timing, duration, rise time]
stac_params = [0 6 2 1]; 
nad_params = [0 6 2 1];

% simulate a well-fed young mouse
clock_metab
ind2 = [116:163]-1; x_fed=x;
ind = [120:167]-1; x_fed=x;
bmal_fed = x(ind2,73);
prot_bmal_fed = x(ind,78);  
prot_cry_fed = x(ind,75);  
prot_per_fed = x(ind,74);  
ror_fed = x(ind2,72); 
rev_fed = x(ind2,71); 
cry_fed = x(ind2,70);
per_fed = x(ind2,69);
cb_fed = x(ind,80);
mtorc1_fed = x(ind,17);
mtorc1_fed2 = x(ind2,17);
ampk_fed = x(ind,8);
sirt_fed = x(ind,56);
irs_fed = x(ind,6); % IRS_pS636
akt_fed = x(ind,12);  % Akt_pT308_pS473
tsc_fed = x(ind,14);  % TSC1_TSC2_pT1462
time_fed = t(ind)-t(ind(1));

% now don't feed it
simulation = 2;  % fasting (insulin=0.1, AA=0.5), no perturbation
clock_metab;  x_fast=x;
bmal_fast = x(ind,73);
prot_bmal_fast = x(ind,78);
prot_cry_fast = x(ind,75);  
prot_per_fast = x(ind,74);  
ror_fast = x(ind,72);
rev_fast = x(ind,71);
cry_fast = x(ind,70);
per_fast = x(ind,69);
cb_fast = x(ind,80);
mtorc1_fast = x(ind,17);
ampk_fast = x(ind,8);
sirt_fast = x(ind,56);
irs_fast = x(ind,6); % IRS_pS636
akt_fast = x(ind,12);  % Akt_pT308_pS473
tsc_fast = x(ind,14);  % TSC1_TSC2_pT1462
time_fast = t(ind)-t(ind(1));

simulation = 3;  % eat half the day
clock_metab;   x_half1=x;
bmal_half1 = x(ind,73); 
prot_bmal_half1 = x(ind,78);
prot_cry_half1 = x(ind,75);  
prot_per_half1 = x(ind,74);  
ror_half1 = x(ind,72);
rev_half1 = x(ind,71);
cry_half1 = x(ind,70);
per_half1 = x(ind,69);
cb_half1 = x(ind,80);
mtorc1_half1 = x(ind,17);
ampk_half1 = x(ind,8);
sirt_half1 = x(ind,56);
irs_half1 = x(ind,6); % IRS_pS636
akt_half1 = x(ind,12);  % Akt_pT308_pS473
tsc_half1 = x(ind,14);  % TSC1_TSC2_pT1462
time_half1 = t(ind)-t(ind(1));

simulation = 4;  % eat the other half the day
clock_metab;  x_half2=x;
bmal_half2 = x(ind,73);
prot_bmal_half2 = x(ind,78);
prot_cry_half2 = x(ind,75);  
prot_per_half2 = x(ind,74); 
ror_half2 = x(ind+24,72);  
rev_half2 = x(ind,71);  
cry_half2 = x(ind,70);
per_half2 = x(ind,69);
cb_half2 = x(ind,80);
mtorc1_half2 = x(ind,17);   
ampk_half2 = x(ind,8);
sirt_half2 = x(ind,56);
irs_half2 = x(ind,6); % IRS_pS636
akt_half2 = x(ind,12);  % Akt_pT308_pS473
tsc_half2 = x(ind,14);  % TSC1_TSC2_pT1462
time_half2 = t(ind)-t(ind(1));

% close all
% plot baseline mRNA results to compare with data
figure(1), clf
% h=subplot(3,2,[1:2]);
% set(h,'position',[.3 .7 .43 .25]);
subplot(3,2,1)
l1 = plot(time_fed,bmal_fed,'-b');
ylabel('Bmal1')
circadian_plots_appearance
axis([0 48 0 2.5])
text(1,2.3,'A','Fontsize',18)
% now add experimental data
bmal1=readtable('bmal1.txt');
bmal1.Properties.VariableNames{'Var1'} = 'time';
bmal1.Properties.VariableNames{'Var2'} = 'value';
bmal1=table2array(bmal1);
bmalc = spline(bmal1(:,1),bmal1(:,2),[1:48])';
hold on
plot([1:48],bmalc,'ok','LineWidth',1,'MarkerFaceColor','k')

subplot(3,2,2)
l1 = plot(time_fed,mtorc1_fed2,'-b');
ylabel('mTORC1')
circadian_plots_appearance
axis([0 48 2.5 5])
text(1,4.8,'B','Fontsize',18)

subplot(3,2,3)
l1 = plot(time_fed,per_fed,'-b');
ylabel('Per2')
circadian_plots_appearance
axis([0 48 0 2.5])
text(1,2.3,'C','Fontsize',18)
% now add experimental data
per_data=readtable('per.txt');
per_data.Properties.VariableNames{'Var1'} = 'time';
per_data.Properties.VariableNames{'Var2'} = 'value';
per_data=table2array(per_data);
per_data2 = spline(per_data(:,1),per_data(:,2),[1:48])';
hold on
plot(per_data2,'ok','LineWidth',1,'MarkerFaceColor','k')

subplot(3,2,4)
l1 = plot(time_fed,cry_fed,'-b');
ylabel('Cry1')
circadian_plots_appearance
axis([0 48 0 2.5])
text(1,2.3,'D','Fontsize',18)
% now add experimental data
cry_data=readtable('cry.txt');
cry_data.Properties.VariableNames{'Var1'} = 'time';
cry_data.Properties.VariableNames{'Var2'} = 'value';
cry_data=table2array(cry_data);
cry_data2 = spline(cry_data(:,1),cry_data(:,2),[1:48])';
hold on
plot([1:48],cry_data2,'ok','LineWidth',1,'MarkerFaceColor','k')

subplot(3,2,5)
l1 = plot(time_fed,rev_fed,'-b');
ylabel('Rev-Erb')
circadian_plots_appearance
axis([0 48 0 5])
text(1,4.5,'E','Fontsize',18)
% now add experimental data
rev_data=readtable('REV.txt');
rev_data.Properties.VariableNames{'Var1'} = 'time';
rev_data.Properties.VariableNames{'Var2'} = 'value';
rev_data=table2array(rev_data);
rev_data2 = spline(rev_data(:,1),rev_data(:,2),[1:48])';
hold on
plot(rev_data2,'ok','LineWidth',1,'MarkerFaceColor','k')

subplot(3,2,6)
l1 = plot(time_fed,ror_fed,'-b');
ylabel('Ror')
circadian_plots_appearance
axis([0 48 0 3])
text(1,2.7,'F','Fontsize',18)
% now add experimental data
ror_data=readtable('ROR.txt');
ror_data.Properties.VariableNames{'Var1'} = 'time';
ror_data.Properties.VariableNames{'Var2'} = 'value';
ror_data=table2array(ror_data);
ror_data2 = spline(ror_data(:,1),ror_data(:,2),[1:48])';
hold on
plot(ror_data2,'ok','LineWidth',1,'MarkerFaceColor','k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% now make plots to compare 4 feeding cases
figure(2), clf
subplot(3,2,1)
l1 = plot(time_fed,ones(1,length(ind)));
hold on
l2 = plot(time_fast,0.1*ones(1,length(ind)));
l3 = plot(time_half1,0.1+0.9*get_P(time_half1,6,12,0.1,24));
l4 = plot(time_half2,0.1+0.9*get_P(time_half2,18,12,0.1,24));
ylabel('Insulin')
legend({'Fed','Fast','Day-feeding','Night-feeding'},'Location','best')
circadian_plots_appearance
axis([0 48 0 1.1])
text(1,1,'A','Fontsize',18)

subplot(3,2,2)
l1 = plot(time_fed,mtorc1_fed,'-b');
hold on
l2 = plot(time_fast,mtorc1_fast,'-r');
l3 = plot(time_half1,mtorc1_half1,'-k');
l4 = plot(time_half2,mtorc1_half2,'-g');
ylabel('mTORC1')
circadian_plots_appearance
axis([0 48 0 5])
text(1,4.6,'B','Fontsize',18)

subplot(3,2,3)
refval = mean(prot_bmal_fed);
l1 = plot(time_fed,prot_bmal_fed/refval,'-b');
hold on
l2 = plot(time_fast,prot_bmal_fast/refval,'-r');
l3 = plot(time_half1,prot_bmal_half1/refval,'-k');
l4 = plot(time_half2,prot_bmal_half2/refval,'-g');
ylabel('BMAL1')
circadian_plots_appearance
axis([0 48 0 2])
text(1,1.8,'C','Fontsize',18)

subplot(3,2,4)
refval = mean(cb_fed);
l1 = plot(time_fed,cb_fed/refval,'-b');
hold on
l2 = plot(time_fast,cb_fast/refval,'-r');
l3 = plot(time_half1,cb_half1/refval,'-k');
l4 = plot(time_half2,cb_half2/refval,'-g');
ylabel('CLOCK-BMAL1')
circadian_plots_appearance
axis([0 48 0 2])
text(1,1.8,'D','Fontsize',18)

subplot(3,2,5)
refval = mean(prot_per_fed);
l1 = plot(time_fed,prot_per_fed/refval,'-b');
hold on
l2 = plot(time_fast,prot_per_fast/refval,'-r');
l3 = plot(time_half1,prot_per_half1/refval,'-k');
l4 = plot(time_half2,prot_per_half2/refval,'-g');
ylabel('PER2')
circadian_plots_appearance
axis([0 48 0 5.5])
text(1,5,'E','Fontsize',18)

subplot(3,2,6)
refval = mean(prot_cry_fed);
l1 = plot(time_fed,prot_cry_fed/refval,'-b');
hold on
l2 = plot(time_fast,prot_cry_fast/refval,'-r');
l3 = plot(time_half1,prot_cry_half1/refval,'-k');
l4 = plot(time_half2,prot_cry_half2/refval,'-g');
ylabel('CRY1')
circadian_plots_appearance
axis([0 48 0 2.5])
text(1,2.3,'F','Fontsize',18)
