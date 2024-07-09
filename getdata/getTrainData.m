%% load data

deg = pi / 180;

% 25m÷±Ω«Õ‰≤‚ ‘1  
load('2122.mat')
startpoint1 = 3800;
endpoint1 = 7800;
data_len1 = size(qqd_est_out(startpoint1:endpoint1,:),1);
time1 = linspace(0,(data_len1-1)*0.01,data_len1)';
s1.qqd1 = qqd_est_out(startpoint1:endpoint1, :);
s1.u1 = steer_m(startpoint1:endpoint1, :) * deg;
u1 = zeros(data_len1, 8);
u1(:, 1) = s1.u1(:,1);
u1(:, 2) = 0.7 * s1.u1(:,1);
u1(:, 3:6) = s1.u1(:,2:5);
u1(:, 7) = 0.7 * s1.u1(:,6);
u1(:, 8) = s1.u1(:,6);
qqd1 = s1.qqd1;
% low filter
qqd_tmp1 = [];
for i = 1:10
    qqd_tmp1(:,i) = butterfilt(s1.qqd1(:,i), cutoff, 4, 'Low', time1, 'off');
end
u_tmp1 = [];
for i = 1:8
    u_tmp1(:,i) = butterfilt(u1(:,i), cutoff, 4, 'Low', time1, 'off');
end
qqd1 = qqd_tmp1;
u1 = u_tmp1;

% 25m÷±Ω«Õ‰≤‚ ‘2
load('2146.mat')
startpoint2 = 2000;
endpoint2 = 5000;
data_len2 = size(qqd_est_out(startpoint2:endpoint2,:),1);
time2 = linspace(0,(data_len2-1)*0.01,data_len2)';
s1.qqd2 = qqd_est_out(startpoint2:endpoint2, :);
s1.u2 = steer_m(startpoint2:endpoint2, :) * deg;
u2 = zeros(data_len2, 8);
u2(:, 1) = s1.u2(:,1);
u2(:, 2) = 0.7 * s1.u2(:,1);
u2(:, 3:6) = s1.u2(:,2:5);
u2(:, 7) = 0.7 * s1.u2(:,6);
u2(:, 8) = s1.u2(:,6);
qqd2 = s1.qqd2;
% low filter
qqd_tmp2 = [];
for i = 1:10
    qqd_tmp2(:,i) = butterfilt(s1.qqd2(:,i), cutoff, 4, 'Low', time2, 'off');
end
u_tmp2 = [];
for i = 1:8
    u_tmp2(:,i) = butterfilt(u2(:,i), cutoff, 4, 'Low', time2, 'off');
end
qqd2 = qqd_tmp2;
u2 = u_tmp2;

% 25m÷±Ω«Õ‰≤‚ ‘3
load('2148.mat')
startpoint3 = 2000;
endpoint3 = 5000;
data_len3 = size(qqd_est_out(startpoint3:endpoint3,:),1);
time3 = linspace(0,(data_len3-1)*0.01,data_len3)';
s1.qqd3 = qqd_est_out(startpoint3:endpoint3, :);
s1.u3 = steer_m(startpoint3:endpoint3, :) * deg;
u3 = zeros(data_len3, 8);
u3(:, 1) = s1.u3(:,1);
u3(:, 2) = 0.7 * s1.u3(:,1);
u3(:, 3:6) = s1.u3(:,2:5);
u3(:, 7) = 0.7 * s1.u3(:,6);
u3(:, 8) = s1.u3(:,6);
qqd3 = s1.qqd3;
% low filter
qqd_tmp3 = [];
for i = 1:10
    qqd_tmp3(:,i) = butterfilt(s1.qqd3(:,i), cutoff, 4, 'Low', time3, 'off');
end
u_tmp3 = [];
for i = 1:8
    u_tmp3(:,i) = butterfilt(u3(:,i), cutoff, 4, 'Low', time3, 'off');
end
qqd3 = qqd_tmp3;
u3 = u_tmp3;

% 15m÷±Ω«Õ‰≤‚ ‘1  
load('2202.mat')
startpoint4 = 3000;
endpoint4 = 7000;
data_len4 = size(qqd_est_out(startpoint4:endpoint4,:),1);
time4 = linspace(0,(data_len4-1)*0.01,data_len4)';
s1.qqd4 = qqd_est_out(startpoint4:endpoint4, :);
s1.u4 = steer_m(startpoint4:endpoint4, :) * deg;
u4 = zeros(data_len4, 8);
u4(:, 1) = s1.u4(:,1);
u4(:, 2) = 0.7 * s1.u4(:,1);
u4(:, 3:6) = s1.u4(:,2:5);
u4(:, 7) = 0.7 * s1.u4(:,6);
u4(:, 8) = s1.u4(:,6);
qqd4 = s1.qqd4;
% low filter
qqd_tmp4 = [];
for i = 1:10
    qqd_tmp4(:,i) = butterfilt(s1.qqd4(:,i), cutoff, 4, 'Low', time4, 'off');
end
u_tmp4 = [];
for i = 1:8
    u_tmp4(:,i) = butterfilt(u4(:,i), cutoff, 4, 'Low', time4, 'off');
end
qqd4 = qqd_tmp4;
u4 = u_tmp4;

% 15m÷±Ω«Õ‰≤‚ ‘2
load('2241.mat')
startpoint5 = 3000;
endpoint5 = 7000;
qqd_est_out = [q_est qd_est];
data_len5 = size(qqd_est_out(startpoint5:endpoint5,:),1);
time5 = linspace(0,(data_len5-1)*0.01,data_len5)';
s1.qqd5 = qqd_est_out(startpoint5:endpoint5, :);
s1.u5 = steer_m(startpoint5:endpoint5, :) * deg;
u5 = zeros(data_len5, 8);
u5(:, 1) = s1.u5(:,1);
u5(:, 2) = 0.7 * s1.u5(:,1);
u5(:, 3:6) = s1.u5(:,2:5);
u5(:, 7) = 0.7 * s1.u5(:,6);
u5(:, 8) = s1.u5(:,6);
qqd5 = s1.qqd5;
% low filter
qqd_tmp5 = [];
for i = 1:10
    qqd_tmp5(:,i) = butterfilt(s1.qqd5(:,i), cutoff, 4, 'Low', time5, 'off');
end
u_tmp5 = [];
for i = 1:8
    u_tmp5(:,i) = butterfilt(u5(:,i), cutoff, 4, 'Low', time5, 'off');
end
qqd5 = qqd_tmp5;
u5 = u_tmp5;

% 15m‘≤«˙œﬂ≤‚ ‘
load('2314.mat')
startpoint6 = 21000;
endpoint6 = 30000;
data_len6 = size(qqd_est_out(startpoint6:endpoint6,:),1);
time6 = linspace(0,(data_len6-1)*0.01,data_len6)';
s1.qqd6 = qqd_est_out(startpoint6:endpoint6, :);
s1.u6 = steer_m(startpoint6:endpoint6, :) * deg;
u6 = zeros(data_len6, 8);
u6(:, 1) = s1.u6(:,1);
u6(:, 2) = 0.7 * s1.u6(:,1);
u6(:, 3:6) = s1.u6(:,2:5);
u6(:, 7) = 0.7 * s1.u6(:,6);
u6(:, 8) = s1.u6(:,6);
qqd6 = s1.qqd6;
% low filter
qqd_tmp6 = [];
for i = 1:10
    qqd_tmp6(:,i) = butterfilt(s1.qqd6(:,i), cutoff, 4, 'Low', time6, 'off');
end
u_tmp6 = [];
for i = 1:8
    u_tmp6(:,i) = butterfilt(u6(:,i), cutoff, 4, 'Low', time6, 'off');
end
qqd6 = qqd_tmp6;
u6 = u_tmp6;

% º«¬º∆ ºµ„
sep = [1 endpoint1-startpoint1+1; ...
    endpoint1-startpoint1+2 endpoint1-startpoint1+2+endpoint2-startpoint2; ...
    endpoint1-startpoint1+3+endpoint2-startpoint2 endpoint1-startpoint1+3+endpoint2-startpoint2+endpoint3-startpoint3; ...
    endpoint1-startpoint1+4+endpoint2-startpoint2+endpoint3-startpoint3 ...
    endpoint1-startpoint1+4+endpoint2-startpoint2+endpoint3-startpoint3+endpoint4-startpoint4; ...
    endpoint1-startpoint1+5+endpoint2-startpoint2+endpoint3-startpoint3+endpoint4-startpoint4 ...
    endpoint1-startpoint1+5+endpoint2-startpoint2+endpoint3-startpoint3+endpoint4-startpoint4+endpoint5-startpoint5; ...
    endpoint1-startpoint1+6+endpoint2-startpoint2+endpoint3-startpoint3+endpoint4-startpoint4+endpoint5-startpoint5 ...
    endpoint1-startpoint1+6+endpoint2-startpoint2+endpoint3-startpoint3+endpoint4-startpoint4+endpoint5-startpoint5+endpoint6-startpoint6;];
    