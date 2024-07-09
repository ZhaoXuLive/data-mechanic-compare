%% Sparse Identification of Nonlinear Dynamics for discrepancy modeling & learning missing physics

% ���������ģ����������
% % ���������ʵģ�ͣ����� + EPS* ���죩����
% % ������߲��������� + EPS* ���� + ����������
% �����ݷ�Ϊѵ�����Ͳ��Լ�
% % ������������ź�
% % ��ȥ�ź��Ը������ 
% ��ѵ�������źŽ��лع�
% ��ȷ���Ĳ���ع�ģ����ǿ����ģ��
% Ԥ����ǿ������ģ�Ͳ��������
% �����µĲ������ݽ�������ͬ��

addpath('./dynamics'); % system dynamics
addpath('./mat'); % data mat
addpath('./filter'); % system filter
addpath('./pool'); % pool dynamics
addpath('./regression'); % regression method
addpath('./getdata'); % get data

%% Generage true dynamics (z) and measurement dynamics (y)

close all
clear all
system = 'Maaav';
polyorder = 3;
lambda = 5;
cutoff = 4;
dt = 0.01;

fitMode = 1;
dataGenerateMode = 1;
poolMode = 1;
sparsifyMode = 3;

% articulated vehicle parameter
avp.m = [11142; 7651; 11142;];
avp.iz = [143570; 37800; 143570;];
avp.lj = [11.4; 9.4; 11.4;];
avp.la = [2.5; 4; 9.9; 1.5; 7.9; 1.5; 7.4; 8.9];
avp.lg = [5.862; 4.643; 5.538;];
avp.h = [0.9625; 0.9625; 0.98; 0.98; 0.98; 0.98; 0.9625; 0.9625;];
avp.d = [10000; 10000];
avp.ca = 3065;

getTrainData;

if fitMode == 1
    eval(['discrepancyDynamics_',system])
end
if fitMode == 2
    eval(['completeDynamics_',system])
end

% SINDy_discrep_IDphysics;
% GPR_discrep_IDphysics;
% NN_discrep_IDphysics;
