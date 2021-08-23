clear
clc
close all

out = csvread("Output2.csv");
%load("n100.mat");
%out = csvread("OutputS1L1.csv");
opt = csvread("Optimal.csv");
opt2 = csvread("Optimal2.csv");
vinf = csvread("Vinf.csv");
vinf2 = csvread("Vinf2.csv");

Vinf1 = vinf(:, 1:3);
Vinf2 = vinf(:, 4:6);
Vinf3 = vinf(:, 7:9);
Vinf4 = vinf(:, 10:12);
Vinf5 = vinf(:, 13:15);
Vinf6 = vinf(:, 16:18);
Vinf7 = vinf(:, 19:21);
Vinf8 = vinf(:, 22:24);

Vinf12 = vinf2(:, 1:3);
Vinf22 = vinf2(:, 4:6);
Vinf32 = vinf2(:, 7:9);
Vinf42 = vinf2(:, 10:12);
Vinf52 = vinf2(:, 13:15);
Vinf62 = vinf2(:, 16:18);
Vinf72 = vinf2(:, 19:21);
Vinf82 = vinf2(:, 22:24);


for i = 1:length(vinf(:, 1))
    VinfE(i) = norm(Vinf1(i, :));
    VinfE2(i) = norm(Vinf12(i, :));
    VinfM(i) = norm(Vinf2(i, :));
    VinfM2(i) = norm(Vinf22(i, :));
end

dim1 = out(1, 1)
dim2 = out(1, 2)
dim3 = out(1, 3)
dim4 = out(1, 4)

nums = out(2:end, 1);
cords = out(2:end, 2:5)';

%dV = zeros(dim1, dim2, dim3, dim4);

d1 = linspace(0, 2*pi-2*pi/dim1, dim1);
d2 = linspace(90, 300-210/dim2, dim2)*24*3600;
d3 = linspace(23, 30-7/dim3, dim3)*30*24*3600;
d4 = linspace(90, 300-210/dim4, dim4)*24*3600;

dVopt = opt(:, 1);
phiOpt = opt(:, 2);
dT1Opt = opt(:, 3);
dT2Opt = opt(:, 4);
dT3Opt = opt(:, 5);

dVopt2 = opt2(:, 1);
phiOpt2 = opt2(:, 2);
dT1Opt2 = opt2(:, 3);
dT2Opt2 = opt2(:, 4);
dT3Opt2 = opt2(:, 5);

count = 1;
countDV = length(nums);
dVBal = nums;

[phiBal, ind1] = sort(cords(1, :));
[dT1Bal, ind2] = sort(cords(2, :));
[dT2Bal, ind3] = sort(cords(3, :));
[dT3Bal, ind4] = sort(cords(4, :));

dVPhi = dVBal(ind1);
dVdT1 = dVBal(ind2);
dVdT2 = dVBal(ind3);
dVdT3 = dVBal(ind4);

for i = 1:countDV
    phi(i) = d1(phiBal(i));
    dT1(i) = d2(dT1Bal(i));
    dT2(i) = d3(dT2Bal(i));
    dT3(i) = d4(dT3Bal(i));
end

% dV output
figure;
plot(1:countDV, sort(dVBal), 'kx');
title("Initial Solver Results \DeltaV Sorted")
ylabel("\DeltaV (km/s)")
xlabel("Solution")

% dV Combined
figure;
plot(1:countDV, sort(dVBal), 'kx');
hold on
plot(1:length(dVopt), sort(dVopt), 'rx');
title("Inital & Optimized Results \DeltaV Sorted")
ylabel("\DeltaV (km/s)")
xlabel("Solution")
legend("Inital", "Optimal")

% dV Initial Opt & Opt2
figure;
plot(1:countDV, sort(dVBal), 'kx');
hold on
plot(1:length(dVopt), sort(dVopt), 'rx');
hold on
plot(1:length(dVopt2), sort(dVopt2), 'bx');
title("Inital, \DeltaV Optimized, & Taxi Optimized Results \DeltaV Sorted")
ylabel("\DeltaV (km/s)")
xlabel("Solution")
legend("Inital", "\DeltaV Optimal", "Taxi Optimal")

% dV Optimal
figure;
plot(1:length(dVopt), sort(dVopt), 'kx');
title("\Delta V Optimal")

% Regular Run
figure;
title("Independent Variable Domain Analysis")
subplot(2, 2, 1)
plot(1:countDV, phi, 'kx')
title("Phi")
ylabel("Rad")

subplot(2, 2, 2)
plot(1:countDV, dT1/3600/24, 'kx')
title("dT1")
ylabel("Days")

subplot(2, 2, 3)
plot(1:countDV, dT2/3600/24/30, 'kx')
title("dT2")
ylabel("Months")

subplot(2, 2, 4)
plot(1:countDV, dT3/3600/24, 'kx')
title("dT3")
ylabel("Days")

% Combined Regular & Optimized Run
figure;
subplot(2, 2, 1)
plot(1:length(dVopt), sort(phiOpt), 'kx')
hold on 
plot(1:countDV, phi, 'rx')
title("Phi Combined")
ylabel("Rad")

subplot(2, 2, 2)
plot(1:length(dVopt), sort(dT1Opt)/3600/24, 'kx')
hold on
plot(1:countDV, dT1/3600/24, 'rx')
title("dT1 Combined")
ylabel("Days")

subplot(2, 2, 3)
plot(1:length(dVopt), sort(dT2Opt)/3600/24/30, 'kx')
hold on
plot(1:countDV, dT2/3600/24/30, 'rx')
title("dT2 Combined")
ylabel("Months")

subplot(2, 2, 4)
plot(1:length(dVopt), sort(dT3Opt)/3600/24, 'kx')
hold on
plot(1:countDV, dT3/3600/24, 'rx')
title("dT3 Combined")
ylabel("Days")
legend("Inital", "Optimal")

% Optimized Run
[dVOptSort, ind] = sort(dVopt);
figure;
subplot(2, 2, 1)
plot(phiOpt(ind), dVOptSort, 'kx')
ylim([0 0.1])
title("Phi Optimal")
xlabel("Rad")
ylabel("\DeltaV (km/s)")

subplot(2, 2, 2)
plot(dT1Opt(ind)/3600/24, dVOptSort, 'kx')
ylim([0 0.1])
title("dT1 Optimal")
xlabel("Days")
ylabel("\DeltaV (km/s)")

subplot(2, 2, 3)
plot(dT2Opt(ind)/3600/24/30, dVOptSort, 'kx')
ylim([0 0.1])
title("dT2 Optimal")
xlabel("Months")
ylabel("\DeltaV (km/s)")

subplot(2, 2, 4)
plot(dT3Opt(ind)/3600/24, dVOptSort, 'kx')
ylim([0 0.1])
title("dT3 Optimal")
xlabel("Days")
ylabel("\DeltaV (km/s)")

figure; 
scatter(dT1Opt/3600/24, dT3Opt/3600/24, 50, dVopt, 'filled')
title("Travel Time Analysis")
xlabel("dT1 (Days)")
ylabel("dT3 (Days)")
color = colorbar;
color.Label.String = 'Delta V (km/s)';

figure; 
scatter(dT1Opt/3600/24, dT3Opt/3600/24, 50, dVopt, 'filled')
title("Travel Time Analysis")
xlabel("dT1 (Days)")
ylabel("dT3 (Days)")
color = colorbar;
color.Label.String = 'Delta V (km/s)';

figure; 
scatter(dT1Opt/3600/24, dT3Opt/3600/24, 50, VinfE, 'filled')
title("V_\infty Earth Analysis")
xlabel("dT1 (Days)")
ylabel("dT3 (Days)")
color = colorbar;
color.Label.String = 'V_\infty Earth (km/s)';

figure; 
scatter(dT1Opt/3600/24, dT3Opt/3600/24, 50, VinfM, 'filled')
title("V_\infty Mars Analysis")
xlabel("dT1 (Days)")
ylabel("dT3 (Days)")
color = colorbar;
color.Label.String = 'V_\infty Mars (km/s)';

%Earth
mu = 3.986e5;
alt = 400;
Vc = sqrt(mu/(6378 + alt));

V = sqrt(VinfE.^2 + (sqrt(2)*Vc)^2);
V2 = sqrt(VinfE2.^2 + (sqrt(2)*Vc)^2);
dVTaxi = V - Vc;
dVTaxi2 = V2 - Vc;

% Mars
mu = 4.2828e4;
alt = 200;
Vc = sqrt(mu/(3396 + alt));

V = sqrt(VinfM.^2 + (sqrt(2)*Vc)^2);
V2 = sqrt(VinfM2.^2 + (sqrt(2)*Vc)^2);
dVTaxiM = V - Vc;
dVTaxiM2 = V2 - Vc;

figure; 
plot(1:length(VinfE), sort(dVTaxi2) - sort(dVTaxi), 'kx')
xlabel("Solution")
ylabel("Taxi \DeltaV Difference (km/s)")
title("Earth Taxi Comparison Optimizer2 vs Optimizer 1")

figure; 
plot(1:length(VinfE), sort(dVTaxiM2) - sort(dVTaxiM), 'kx')
xlabel("Solution")
ylabel("Taxi \DeltaV Difference (km/s)")
title("Mars Taxi Comparison Optimizer2 vs Optimizer 1")

figure; 
plot(1:length(VinfE), sort(dVopt2) - sort(dVopt), 'kx')
xlabel("Solution")
ylabel("Cycler \DeltaV Difference (km/s)")
title("Cycler Comparison Optimizer2 vs Optimizer 1")