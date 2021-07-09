clear
clc
close all

R = csvread("SingleCycleR.csv");
V = csvread("SingleCycleV.csv");

mu = 1.32712440e11;
re = 1.495979e8;
rm = 2.279483e8;

theta = linspace(0, 2*pi, 1000);

figure;
plot(re*cos(theta)/re, re*sin(theta)/re, 'b');
hold on
plot(rm*cos(theta)/re, rm*sin(theta)/re, 'r');
hold on

for i = 1:4
    OE = SV2OE(mu, R(i, :), V(i, :));
    
    r = OE(1)^2/mu./(1 + OE(2)*cos(theta));
    x = r.*cos(theta + OE(5));
    y = r.*sin(theta + OE(5));
    
    plot(x/re, y/re);
    hold on
end

legend("Earth", "Mars", "Orbit1", "Orbit2", "Orbit3", "Orbit4")
axis equal
xlabel("X (AU)");
ylabel("Y (AU)");

figure;
for i = 1:4
    OE = SV2OE(mu, R(i, :), V(i, :));
    
    r = OE(1)^2/mu./(1 + OE(2)*cos(theta));
    x = r.*cos(theta + OE(5));
    y = r.*sin(theta + OE(5));
    
    subplot(2, 2, i)
    plot(re*cos(theta)/re, re*sin(theta)/re, 'b');
    hold on
    plot(rm*cos(theta)/re, rm*sin(theta)/re, 'r');
    hold on
    plot(x/re, y/re);
    axis equal
    xlabel("X (AU)");
    ylabel("Y (AU)");
end

