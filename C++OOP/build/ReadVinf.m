clear
clc
close all

test = csvread("Vinf.csv");

VinfE1 = test(:, 1);
VinfM2 = test(:, 2);
VinfM31 = test(:, 3);
VinfM32 = test(:, 4); % Always 0
VinfM41 = test(:, 5);
VinfM42 = test(:, 6); % Always 0
VinfM5 = test(:, 7);
VinfE6 = test(:, 8);
VinfE71 = test(:, 9); %Sometimes 0
VinfE72 = test(:, 10);
VinfE81 = test(:, 11); %Sometimes 0
VinfE82 = test(:, 12);

