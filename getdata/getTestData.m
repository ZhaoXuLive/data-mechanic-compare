%% load data

% Test1 
load('2108.mat')
startpoint_t1 = 5000;
endpoint_t1 = 11000;
qqd_est_out = [q_est qd_est];
data_len_t1 = size(qqd_est_out(startpoint_t1:endpoint_t1,:),1);
time_t1 = linspace(0,(data_len_t1-1)*0.01,data_len_t1)';
s1.qqd_t1 = qqd_est_out(startpoint_t1:endpoint_t1, :);
s1.u_t1 = steer_m(startpoint_t1:endpoint_t1, :) * deg;
u_t1 = zeros(data_len_t1, 8);
u_t1(:, 1) = s1.u_t1(:,1);
u_t1(:, 2) = 0.7 * s1.u_t1(:,1);
u_t1(:, 3:6) = s1.u_t1(:,2:5);
u_t1(:, 7) = 0.7 * s1.u_t1(:,6);
u_t1(:, 8) = s1.u_t1(:,6);
qqd_t1 = s1.qqd_t1;
% low filter
qqd_tmp_t1 = [];
for i = 1:10
    qqd_tmp_t1(:,i) = butterfilt(s1.qqd_t1(:,i), cutoff, 4, 'Low', time_t1, 'off');
end
u_tmp_t1 = [];
for i = 1:8
    u_tmp_t1(:,i) = butterfilt(u_t1(:,i), cutoff, 4, 'Low', time_t1, 'off');
end
qqd_t1 = qqd_tmp_t1;
u_t1 = u_tmp_t1;

% Test2
load('2237.mat')
startpoint_t2 = 8300;
endpoint_t2 = 12000;
qqd_est_out = [q_est qd_est];
data_len_t2 = size(q_est(startpoint_t2:endpoint_t2,:),1);
time_t2 = linspace(0,(data_len_t2-1)*0.01,data_len_t2)';
qqd_est_out = [q_est qd_est];
s1.qqd_t2 = qqd_est_out(startpoint_t2:endpoint_t2, :);
s1.u_t2 = steer_m(startpoint_t2:endpoint_t2, :) * deg;
u_t2 = zeros(data_len_t2, 8);
u_t2(:, 1) = s1.u_t2(:,1);
u_t2(:, 2) = 0.7 * s1.u_t2(:,1);
u_t2(:, 3:6) = s1.u_t2(:,2:5);
u_t2(:, 7) = 0.7 * s1.u_t2(:,6);
u_t2(:, 8) = s1.u_t2(:,6);
qqd_t2 = s1.qqd_t2;
% low filter
qqd_tmp_t2 = [];
for i = 1:10
    qqd_tmp_t2(:,i) = butterfilt(s1.qqd_t2(:,i), cutoff, 4, 'Low', time_t2, 'off');
end
u_tmp_t2 = [];
for i = 1:8
    u_tmp_t2(:,i) = butterfilt(u_t2(:,i), cutoff, 4, 'Low', time_t2, 'off');
end
qqd_t2 = qqd_tmp_t2;
u_t2 = u_tmp_t2;

% % Test3 quarter
% load('2259.mat')
% startpoint_t3 = 6000;
% endpoint_t3 = 12000;
% qqd_est_out = [q_est qd_est];
% data_len_t3 = size(q_est(startpoint_t3:endpoint_t3,:),1);
% time_t3 = linspace(0,(data_len_t3-1)*0.01,data_len_t3)';
% qqd_est_out = [q_est qd_est];
% s1.qqd_t3 = qqd_est_out(startpoint_t3:endpoint_t3, :);
% s1.u_t3 = steer_m(startpoint_t3:endpoint_t3, :) * deg;
% u_t3 = zeros(data_len_t3, 8);
% u_t3(:, 1) = s1.u_t3(:,1);
% u_t3(:, 2) = 0.7 * s1.u_t3(:,1);
% u_t3(:, 3:6) = s1.u_t3(:,2:5);
% u_t3(:, 7) = 0.7 * s1.u_t3(:,6);
% u_t3(:, 8) = s1.u_t3(:,6);
% qqd_t3 = s1.qqd_t3;
% % low filter
% qqd_tmp_t3 = [];
% for i = 1:10
%     qqd_tmp_t3(:,i) = butterfilt(s1.qqd_t3(:,i), cutoff, 4, 'Low', time_t3, 'off');
% end
% u_tmp_t3 = [];
% for i = 1:8
%     u_tmp_t3(:,i) = butterfilt(u_t3(:,i), cutoff, 4, 'Low', time_t3, 'off');
% end
% qqd_t3 = qqd_tmp_t3;
% u_t3 = u_tmp_t3;

TestY = [qqd_t1; qqd_t2;];% qqd_t3];
TestU = [u_t1; u_t2;];% u_t3];
dt = time_t1(2) - time_t1(1);
TestT = [time_t1; time_t2+time_t1(end)+dt;];% time_t3+time_t1(end)+time_t2(end)+dt+dt;];
sepT = [1 endpoint_t1-startpoint_t1+1; ...
    endpoint_t1-startpoint_t1+2 endpoint_t1-startpoint_t1+2+endpoint_t2-startpoint_t2;];% ...
%     endpoint_t1-startpoint_t1+3+endpoint_t2-startpoint_t2 ...
%     endpoint_t1-startpoint_t1+3+endpoint_t2-startpoint_t2+endpoint_t3-startpoint_t3;];

if dataGenerateMode == 1
    Testx1 = Maaav_single(time_t1, qqd_t1, avp, u_t1, n);
    Testx2 = Maaav_single(time_t2, qqd_t2, avp, u_t2, n);
%     Testx3 = Maaav_single(time_t3, qqd_t3, avp, u_t3, n);
end
if dataGenerateMode == 2
    Testx1 = Maaav_iter(time_t1, qqd_t1, avp, u_t1, n);
    Testx2 = Maaav_iter(time_t2, qqd_t2, avp, u_t2, n);
%     Testx3 = Maaav_iter(time_t3, qqd_t3, avp, u_t3, n);
end
TestX = [Testx1; Testx2;];% Testx3];
