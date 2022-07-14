% plot radial dose fall off

% load data
load("result.mat")
load('ct.mat', 'ct');
load('tplan_full_orig.mat', 'full_tplan_orig');
load('basedata/brachy_LDR.mat')

% matRad dose fall off
matrad_dose = resultGUI.physicalDose;
mdose = matrad_dose(:,:,2);
% seed position
p = full_tplan_orig{4,2}(1);
% get max dose value
[c,i] = max(mdose(:));
[i1,i2] = ind2sub(size(mdose), i);
line = mdose(i1,:);
line(line>1000) = 1000;

% z slice
mdose_z = squeeze(matrad_dose(i2,:,:));
line_z = mdose_z(i2,:);
x_val_line_z = [-5, 0, 5, 10, 15];

% imported dose fall off
ind = [3,4,6,8,9,11,13];
x_val = p+machine.data.RadialDoseDistance(ind)*10;
dose_rate_val = machine.data.SourceStrengthImplanted*machine.data.lambda...
    *machine.data.RadialDoseValue(ind) .* machine.data.RadialDoseDistance(ind).^(-2) .* ... 
    machine.data.AnisotropyFactorValue;
dose_val = dose_rate_val * 0.01 * machine.data.SourceIsotopeHalfLife*24/0.693147;

% plot
figure
imagesc(mdose)
colorbar

figure
hold on 
title('dose fall-off')
xlabel('distance in mm')
ylabel('dose in Gy')
plot(x_val, dose_val, '*g')
plot(ct.x,line,'b'); %dose fall off matrad
plot(p+x_val_line_z, line_z, 'or')
plot([p-1e-7,p], ylim,'-r');

legend('specification','matRad x dir', 'matRad z dir', 'seed position')


