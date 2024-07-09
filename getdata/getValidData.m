%% load data

% Validation1 
load('2213.mat')
startpoint_v1 = 5000;
endpoint_v1 = 9000;
qqd_est_out = [q_est qd_est];
data_len_v1 = size(qqd_est_out(startpoint_v1:endpoint_v1,:),1);
time_v1 = linspace(0,(data_len_v1-1)*0.01,data_len_v1)';
s1.qqd_v1 = qqd_est_out(startpoint_v1:endpoint_v1, :);
s1.u_v1 = steer_m(startpoint_v1:endpoint_v1, :) * deg;
u_v1 = zeros(data_len_v1, 8);
u_v1(:, 1) = s1.u_v1(:,1);
u_v1(:, 2) = 0.7 * s1.u_v1(:,1);
u_v1(:, 3:6) = s1.u_v1(:,2:5);
u_v1(:, 7) = 0.7 * s1.u_v1(:,6);
u_v1(:, 8) = s1.u_v1(:,6);
qqd_v1 = s1.qqd_v1;
% low filter
qqd_tmp_v1 = [];
for i = 1:10
    qqd_tmp_v1(:,i) = butterfilt(s1.qqd_v1(:,i), cutoff, 4, 'Low', time_v1, 'off');
end
u_tmp_v1 = [];
for i = 1:8
    u_tmp_v1(:,i) = butterfilt(u_v1(:,i), cutoff, 4, 'Low', time_v1, 'off');
end
qqd_v1 = qqd_tmp_v1;
u_v1 = u_tmp_v1;

% Validation2
load('2259.mat')
startpoint_v2 = 5000;
endpoint_v2 = 14000;
qqd_est_out = [q_est qd_est];
data_len_v2 = size(qqd_est_out(startpoint_v2:endpoint_v2,:),1);
time_v2 = linspace(0,(data_len_v2-1)*0.01,data_len_v2)';
s1.qqd_v2 = qqd_est_out(startpoint_v2:endpoint_v2, :);
s1.u_v2 = steer_m(startpoint_v2:endpoint_v2, :) * deg;
u_v2 = zeros(data_len_v2, 8);
u_v2(:, 1) = s1.u_v2(:,1);
u_v2(:, 2) = 0.7 * s1.u_v2(:,1);
u_v2(:, 3:6) = s1.u_v2(:,2:5);
u_v2(:, 7) = 0.7 * s1.u_v2(:,6);
u_v2(:, 8) = s1.u_v2(:,6);
qqd_v2 = s1.qqd_v2;
% low filter
qqd_tmp_v2 = [];
for i = 1:10
    qqd_tmp_v2(:,i) = butterfilt(s1.qqd_v2(:,i), cutoff, 4, 'Low', time_v2, 'off');
end
u_tmp_v2 = [];
for i = 1:8
    u_tmp_v2(:,i) = butterfilt(u_v2(:,i), cutoff, 4, 'Low', time_v2, 'off');
end
qqd_v2 = qqd_tmp_v2;
u_v2 = u_tmp_v2;

% % Validation3
% load('2237.mat')
% startpoint_v3 = 7000;
% endpoint_v3 = 11000;
% qqd_est_out = [q_est qd_est];
% data_len_v3 = size(qqd_est_out(startpoint_v3:endpoint_v3,:),1);
% time_v3 = linspace(0,(data_len_v3-1)*0.01,data_len_v3)';
% s1.qqd_v3 = qqd_est_out(startpoint_v3:endpoint_v3, :);
% s1.u_v3 = steer_m(startpoint_v3:endpoint_v3, :) * deg;
% u_v3 = zeros(data_len_v3, 8);
% u_v3(:, 1) = s1.u_v3(:,1);
% u_v3(:, 2) = 0.7 * s1.u_v3(:,1);
% u_v3(:, 3:6) = s1.u_v3(:,2:5);
% u_v3(:, 7) = 0.7 * s1.u_v3(:,6);
% u_v3(:, 8) = s1.u_v3(:,6);
% qqd_v3 = s1.qqd_v3;
% % low filter
% qqd_tmp_v3 = [];
% for i = 1:10
%     qqd_tmp_v3(:,i) = butterfilt(s1.qqd_v3(:,i), cutoff, 4, 'Low', time_v3, 'off');
% end
% u_tmp_v3 = [];
% for i = 1:8
%     u_tmp_v3(:,i) = butterfilt(u_v3(:,i), cutoff, 4, 'Low', time_v3, 'off');
% end
% qqd_v3 = qqd_tmp_v3;
% u_v3 = u_tmp_v3;

ValidY = [qqd_v1; qqd_v2;];% qqd_v3];
ValidU = [u_v1; u_v2;];% u_v3];
dt = time_v1(2) - time_v1(1);
ValidT = [time_v1; time_v2+time_v1(end)+dt;];% time_v3+time_v1(end)+time_v2(end)+dt+dt;];
sepV = [1 endpoint_v1-startpoint_v1+1; ...
    endpoint_v1-startpoint_v1+2 endpoint_v1-startpoint_v1+2+endpoint_v2-startpoint_v2;];% ...
%     endpoint_v1-startpoint_v1+3+endpoint_v2-startpoint_v2 ...
%     endpoint_v1-startpoint_v1+3+endpoint_v2-startpoint_v2+endpoint_v3-startpoint_v3;];

if dataGenerateMode == 1
    Validx1 = Maaav_single(time_v1, qqd_v1, avp, u_v1, n);
    Validx2 = Maaav_single(time_v2, qqd_v2, avp, u_v2, n);
%     Validx3 = Maaav_single(time_v3, qqd_v3, avp, u_v3, n);
end
if dataGenerateMode == 2
    Validx1 = Maaav_iter(time_v1, qqd_v1, avp, u_v1, n);
    Validx2 = Maaav_iter(time_v2, qqd_v2, avp, u_v2, n);
%     Validx3 = Maaav_iter(time_v3, qqd_v3, avp, u_v3, n);
end
ValidX = [Validx1; Validx2;];% Validx3];
