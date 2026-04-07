clear all; clc; format long;

data = readtable('PhotoelectricEffectData.csv', 'NumHeaderLines', 11);

% Frequency of (two shades of) violet, blue, green, and yellow light respectively
f = table2array(data(3:7, 5));


% Values for stopping potential imported from CSV file. 
V_0_im = cell(3, 1);

V_0_im{1} = table2array(data(3:7, 6));
V_0_im{2} = table2array(data(3:7, 7));
V_0_im{3} = table2array(data(3:7, 8));

% Converted into a 5x3 numeric array to neatly display data across each trial: each column is for one trial where the frequency decreases as you move down the rows
V_0_tr = cell2mat(V_0_im');

% Values averaged across each row to give single value for stopping potential for each frequency
V_0 = mean(V_0_tr, 2);



% Digital Multimeter Resolution, 0.5%
sigma_inst = 0.005 .* V_0;

% Statistical Error
sigma_stat = std(V_0_tr, 0, 2) ./ sqrt(3);

% Overall Uncertainty
sigma_V0 = sqrt((sigma_inst).^2 + (sigma_stat).^2);



% Stopping Potential vs Frequency Plot
figure(1)
hold on
plot(f, V_0, 'ko', MarkerFaceColor= 'k')

p = polyfit(f, V_0, 1);
V0_fit = polyval(p, f);

plot(f, V0_fit, 'b')


residuals = V_0 - V0_fit;                
chi2 = sum((residuals ./ sigma_V0).^2);  % weighted chi-squared
num_params = 2;                          % slope + intercept
nu = length(V_0) - num_params;           % degrees of freedom
chi2_reduced = chi2 / nu;   

scale_factor = sqrt(chi2_reduced);
sigma_V0_scaled = sigma_V0 * scale_factor;

chi2_scaled = sum((residuals ./ sigma_V0_scaled) .^2);
chi2_reduced_scaled = chi2_scaled / nu;

sigma_sys = 0.055;                                          % Initially set to 0.01 arbitrarily and fine-tuned to get chi2_reduced approximately equal to 1. 
sigma_tuned = sqrt((sigma_V0).^2 + (sigma_sys).^2);

chi2_tuned = sum((residuals ./ sigma_tuned).^2);
chi2_reduced_tuned = chi2_tuned / nu;


sigma_total = sqrt((sigma_V0).^2 + (sigma_sys).^2);

errorbar(f, V_0, sigma_total, 'ko', MarkerFaceColor= 'k')


xlim([min(f), max(f)])

xlabel('Frequency (Hz)', fontsize = 18)
ylabel('Stopping Potential, V_0 (V)', fontsize = 18)

title('Stopping Potential vs Frequency')
set(gca, fontsize = 16)

grid on
exportgraphics(gcf, 'V0vsF.png', 'ContentType', 'vector' )



figure(2)
plot(V0_fit, residuals, 'ko', MarkerFaceColor = 'k')
yline(0, 'k--')

xlabel('Stopping Voltage (V)', fontsize = 18)
ylabel('Residuals', fontsize = 18)
title('Residuals vs Stopping Voltage')


set(gca, fontsize = 16)
grid on
exportgraphics(gcf, 'ResidualsVsV0.png', 'ContentType', 'vector')



% Weights
w = 1 ./ sigma_total.^2;                                                        % Inverse variance weighting formula. Assigns a weight to each uncertainty such that the larger uncertainties are penalised more than the smaller ones. 

% Weighted means
f_bar = sum(w .* f) / sum(w);
V0_bar = sum(w .* V_0) / sum(w);

% Weighted slope
m = sum(w .* (f - f_bar) .* (V_0 - V0_bar)) / sum(w .* (f - f_bar).^2);

% Weighted intercept
c = V0_bar - m * f_bar;

% Uncertainties
sigma_m = sqrt(1 / sum(w .* (f - f_bar).^2));
sigma_c = sqrt(sum(w .* f.^2) / (sum(w) * sum(w .* (f - f_bar).^2)));



% Charge of an electron
e = 1.60217663 *10^(-19);

% Plank's Constant (J): gradient of graph multiplied by e.
h = p(1) * e;

% Work Function (eV): y-intercept of graph.
phi = - p(2);

% Propagate to h and phi
h = m * e;
sigma_h = sigma_m * e;

phi = - c;
sigma_phi = sigma_c;

fprintf('h = %.4e ± %.4e J*s\n', h, sigma_h);
fprintf('phi = %.4e ± %.4e eV\n', phi, sigma_phi);






%{

% Linear Regression Model: Ideal for weighted fit and calculating uncertainty in slope and intercept
model = fitlm(f, V_0, 'Weights', 1 ./ sigma_total.^2);                                      % V_0 = mf + c

m = model.Coefficients.Estimate(2);
sigma_m = model.Coefficients.SE(2);

c = model.Coefficients.Estimate(1);
sigma_c = model.Coefficients.SE(1);



% Plot using model
figure(3)
hold on
errorbar(f, V_0, sigma_total, 'ko', MarkerFaceColor= 'k')

V_0_fit = predict(model, f);

plot(f, V_0_fit, 'r-')

xlim([0, max(f)])

xlabel('Frequency (Hz)', fontsize = 18)
ylabel('Stopping Potential, V_0 (V)', fontsize = 18)

title('Stopping Potential vs Frequency')
set(gca, fontsize = 16)

grid on

%}