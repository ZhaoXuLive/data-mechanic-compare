close all
clc
clear

% load data

% double change
load('1853.mat');
startpoint1 = 8000;
endpoint1 = 11000;
data_len1 = size(q_est(startpoint1:endpoint1, :), 1);
time1 = linspace(0,(data_len1-1)*0.01,data_len1)';
velocity1 = hypot(qd_est(startpoint1:endpoint1,1), qd_est(startpoint1:endpoint1,2)) * 3.6;
figure,
plot(time1, velocity1);
figure,
plot(q_est(startpoint1:endpoint1, 1), q_est(startpoint1:endpoint1, 2));

% quarter
load('2108.mat');
startpoint1 = 7000;
endpoint1 = 11000;
data_len1 = size(q_est(startpoint1:endpoint1, :), 1);
time1 = linspace(0,(data_len1-1)*0.01,data_len1)';
velocity1 = hypot(qd_est(startpoint1:endpoint1,1), qd_est(startpoint1:endpoint1,2)) * 3.6;
figure,
plot(time1, velocity1);
figure,
plot(q_est(startpoint1:endpoint1, 1), q_est(startpoint1:endpoint1, 2));

% quarter
load('2213.mat');
startpoint1 = 4000;
endpoint1 = 11000;
data_len1 = size(q_est(startpoint1:endpoint1, :), 1);
time1 = linspace(0,(data_len1-1)*0.01,data_len1)';
velocity1 = hypot(qd_est(startpoint1:endpoint1,1), qd_est(startpoint1:endpoint1,2)) * 3.6;
figure,
plot(time1, velocity1);
figure,
plot(q_est(startpoint1:endpoint1, 1), q_est(startpoint1:endpoint1, 2));

% quarter
load('2241.mat');
startpoint1 = 3000;
endpoint1 = 7000;
data_len1 = size(q_est(startpoint1:endpoint1, :), 1);
time1 = linspace(0,(data_len1-1)*0.01,data_len1)';
velocity1 = hypot(qd_est(startpoint1:endpoint1,1), qd_est(startpoint1:endpoint1,2)) * 3.6;
figure,
plot(time1, velocity1);
figure,
plot(q_est(startpoint1:endpoint1, 1), q_est(startpoint1:endpoint1, 2));

% 3/4 circle
load('2259.mat');
startpoint1 = 5000;
endpoint1 = 20000;
data_len1 = size(q_est(startpoint1:endpoint1, :), 1);
time1 = linspace(0,(data_len1-1)*0.01,data_len1)';
velocity1 = hypot(qd_est(startpoint1:endpoint1,1), qd_est(startpoint1:endpoint1,2)) * 3.6;
figure,
plot(time1, velocity1);
figure,
plot(q_est(startpoint1:endpoint1, 1), q_est(startpoint1:endpoint1, 2));

% circle + line
load('2322.mat');
startpoint1 = 1000;
endpoint1 = 7000;
data_len1 = size(q_est(startpoint1:endpoint1, :), 1);
time1 = linspace(0,(data_len1-1)*0.01,data_len1)';
velocity1 = hypot(qd_est(startpoint1:endpoint1,1), qd_est(startpoint1:endpoint1,2)) * 3.6;
figure,
plot(time1, velocity1);
figure,
plot(q_est(startpoint1:endpoint1, 1), q_est(startpoint1:endpoint1, 2));

% straight
load('2325.mat');
startpoint1 = 6000;
endpoint1 = 8000;
data_len1 = size(q_est(startpoint1:endpoint1, :), 1);
time1 = linspace(0,(data_len1-1)*0.01,data_len1)';
velocity1 = hypot(qd_est(startpoint1:endpoint1,1), qd_est(startpoint1:endpoint1,2)) * 3.6;
figure,
plot(time1, velocity1);
figure,
plot(q_est(startpoint1:endpoint1, 1), q_est(startpoint1:endpoint1, 2));