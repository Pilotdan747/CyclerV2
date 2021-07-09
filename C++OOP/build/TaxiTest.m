clear
clc
close all


% Leave Earth
mu = 3.986e5;
alt = 400;
Vc = sqrt(mu/(6378 + alt));


Vinf = 2.8;

%V^2 = Vinf^2 + Vesc^2
V = sqrt(Vinf^2 + (sqrt(2)*Vc)^2);
dV = V - Vc;

% Leave Mars