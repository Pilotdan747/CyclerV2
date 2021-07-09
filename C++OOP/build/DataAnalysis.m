clear
clc
close all

out = csvread("Output.csv");
%load("n100.mat");
%out = csvread("OutputS1L1.csv");
opt = csvread("Optimal.csv");

dim1 = out(1)
dim2 = out(2)
dim3 = out(3)
dim4 = out(4)

nums = out(2:end);

dV = zeros(dim1, dim2, dim3, dim4);

d1 = linspace(0, 2*pi-2*pi/dim1, dim1);
d2 = linspace(90, 300-210/dim2, dim2)*24*3600;
d3 = linspace(23, 30-7/dim3, dim3)*30*24*3600;
d4 = linspace(90, 300-210/dim4, dim4)*24*3600;

dVopt = opt(:, 1);
phiOpt = opt(:, 2);
dT1Opt = opt(:, 3);
dT2Opt = opt(:, 4);
dT3Opt = opt(:, 5);

count = 1;
countDV = 0;
for i=1:dim1
    for j=1:dim2
        for k=1:dim3
            for l=1:dim4
                num = nums(count);
                
                dV(i, j, k, l) = num;
                count = count+1;
                
                if (num < 1)
                    countDV = countDV + 1;
                    cords(:, countDV) = [i; j; k; l];
                    cordsLength(countDV) = norm(cords(:, countDV));
                    dVBal(countDV) = num;
                    
                end
            end
        end
    end
end

mindV = min(min(min(min(dV))));


% for i = 1:countDV
%     for j = i:countDV
%         dist(i, j) = norm(cords(:, i) - cords(:, j));
%     end
% end
% 
% lngth = zeros(countDV, 1);
% for i = 1:countDV
%     set = [i];
%     for j = 1:countDV
%         if dist(i, j) == 1
%             set = [set, j];
%             lngth(i) = lngth(i) + 1;
%         else
%             set = [set, nan];
%         end
%     end
%     set = sort(set);
%     neighboors(:, i) = set;
% end
% 
% dist2 = dist.^2;

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
title("\Delta V")

% dV Combined
figure;
plot(1:countDV, sort(dVBal), 'kx');
hold on
plot(1:length(dVopt), sort(dVopt), 'rx');
title("\Delta V Combined")

% dV Optimal
figure;
plot(1:length(dVopt), sort(dVopt), 'kx');
title("\Delta V Optimal")

% Regular Run
figure;
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
xlabel("dT1 (Days)")
ylabel("dT3 (Days)")
color = colorbar;
color.Label.String = 'Delta V (km/s)';
