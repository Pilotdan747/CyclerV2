clear
close all
clc


out = csvread("Vinf.csv");
out2 = csvread("Optimal.csv");

size = length(out(:, 1));

Vinf1 = out(:, 1:3);
Vinf2 = out(:, 4:6);
Vinf3 = out(:, 7:9);
Vinf4 = out(:, 10:12);
Vinf5 = out(:, 13:15);
Vinf6 = out(:, 16:18);
Vinf7 = out(:, 19:21);
Vinf8 = out(:, 22:24);

dV = out2(:, 1);
[dVSort, ind] = sort(dV);

muEarth = 3.986e5;
muMars = 4.2828e4;

Re = 6378;
Rm = 3396;

% Turning Angles
for i = 1:size
    delt1(i) = abs(atan2(Vinf3(i, 2), Vinf3(i, 1)) - atan2(Vinf2(i, 2), Vinf2(i, 1)));
    delt1Max(i) = 2*asin(1/(1 + Rm*norm(Vinf2(i, :))/muMars));
    
    delt2(i) = abs(atan2(Vinf5(i, 2), Vinf5(i, 1)) - atan2(Vinf4(i, 2), Vinf4(i, 1)));
    delt2Max(i) = 2*asin(1/(1 + Rm*norm(Vinf4(i, :))/muMars));
    
    delt3(i) = abs(atan2(Vinf7(i, 2), Vinf7(i, 1)) - atan2(Vinf6(i, 2), Vinf6(i, 1)));
    delt3Max(i) = 2*asin(1/(1 + Rm*norm(Vinf6(i, :))/muEarth));
    
    delt4(i) = abs(atan2(Vinf1(i, 2), Vinf1(i, 1)) - atan2(Vinf8(i, 2), Vinf8(i, 1)));
    delt4Max(i) = 2*asin(1/(1 + Rm*norm(Vinf8(i, :))/muEarth));
end

figure; plot(1:size, delt1(ind), 1:size, delt1Max(ind))
title("Mars 1")
figure; plot(1:size, delt2(ind), 1:size, delt2Max(ind))
title("Mars 2")
figure; plot(1:size, delt3(ind), 1:size, delt3Max(ind))
title("Earth 1")
figure; plot(1:size, delt4(ind), 1:size, delt4Max(ind))
title("Earth 2")

figure; plot(1:size, dVSort);

% Vinf Delta V

for i = 1:size
    dV1(i) = abs(norm(Vinf3(i, :)) - norm(Vinf2(i, :)));
    dV2(i) = abs(norm(Vinf5(i, :)) - norm(Vinf4(i, :)));
    dV3(i) = abs(norm(Vinf7(i, :)) - norm(Vinf6(i, :)));
    dV4(i) = abs(norm(Vinf1(i, :)) - norm(Vinf8(i, :)));
end

figure; plot(1:size, [dV1(ind); dV2(ind); dV3(ind); dV4(ind)]);
legend("Mars1", "Mars2", "Earth1", "Earth2")
