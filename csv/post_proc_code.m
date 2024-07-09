clc;close all;clear;
deg = 180/pi;

% zz_gps = [36.151446 119.998397];
% zz_gps = [36.330000000000000 1.202723000000000e+02];
zz_gps = [36.330000 120.272488];


% data_record = readtable('20210407161203_file1_v2.csv');
data_record = readtable('11282237.csv');
filename = '11282237';
% data_record = readtable('file1(1).csv');
% filename = 'file1(1)';
startpoint = 20;

autoMode = table2array(data_record(startpoint:end,1));
selfGuideMode = table2array(data_record(startpoint:end,2));
manualMode = table2array(data_record(startpoint:end,3));
MC1IsOK = table2array(data_record(startpoint:end,4));
MC2IsOK = table2array(data_record(startpoint:end,5));
yawRate1_m = table2array(data_record(startpoint:end,6));
spdEast1_m = table2array(data_record(startpoint:end,7));
spdNorth1_m = table2array(data_record(startpoint:end,8));
spd1_m = table2array(data_record(startpoint:end,9));
yaw1_m = table2array(data_record(startpoint:end,10));
lon1_m = table2array(data_record(startpoint:end,11));
lat1_m = table2array(data_record(startpoint:end,12));
yawRate3_m = table2array(data_record(startpoint:end,13));
spdEast3_m = table2array(data_record(startpoint:end,14));
spdNorth3_m = table2array(data_record(startpoint:end,15));
% vehicleVelocity3 = table2array(data_record(startpoint:end,16));
yaw3_m = table2array(data_record(startpoint:end,17));
lon3_m = table2array(data_record(startpoint:end,18));
lat3_m = table2array(data_record(startpoint:end,19));
art1_m = table2array(data_record(startpoint:end,20));
art2_m = table2array(data_record(startpoint:end,21));
steer_m = table2array(data_record(startpoint:end,22:27));
steer_cmd = table2array(data_record(startpoint:end,28:33));
q_debug = table2array(data_record(startpoint:end,34:38));
qd_debug = table2array(data_record(startpoint:end,39:43));
ff_error_debug = table2array(data_record(startpoint:end,44:49));
fb_str = table2array(data_record(startpoint:end,50:55));
fb_error_debug = table2array(data_record(startpoint:end,56:61));

% intsteercmd = table2array(data_record(startpoint:end,62:67));
ANGref = table2array(data_record(startpoint:end,68:73));

q_est = q_debug;
qd_est = qd_debug;
u1 = [cos(q_debug(:,3)), sin(q_debug(:,3))];
q_est(:,1:2) = q_debug(:,1:2) - 11.4 * u1;
% q_est(:,1) = q_debug(:,1) - 11.4 * cos(q_debug(:,3));
% q_est(:,2) = q_debug(:,2) - 11.4 * sin(q_debug(:,3));
u2 = [-qd_debug(:,3).*sin(q_debug(:,3)),qd_debug(:,3).*cos(q_debug(:,3))];
qd_est(:,1:2) = qd_debug(:,1:2) - 11.4 * u2;

data_len = size(data_record(startpoint:end,:),1);
time = linspace(0,(data_len-1)*0.01,data_len)';

[x_m,y_m] = LonLat2XY(lat1_m,lon1_m,zz_gps(1),zz_gps(2));
[x_m3,y_m3] = LonLat2XY(lat3_m,lon3_m,zz_gps(1),zz_gps(2));

if MC2IsOK(end-1000) == 1
    steer_m = flip(steer_m,2);
%     yawRate1_m = yawRate1_m + 0.8467;
%     yawRate3_m = yawRate3_m + 0.5832;
%     steer_cmd = flip(steer_cmd,2);
% else
%     yawRate1_m = yawRate1_m + 0.5832;
%     yawRate3_m = yawRate3_m + 0.8467;

end

rwa = max(abs(yawRate3_m))/max(abs(yawRate1_m));
figure('Name','Yaw Rate RWA');
hold on;
plot(time,yawRate1_m,"-","LineWidth",2);
% plot(time,qd_est(:,4)*deg,"--","LineWidth",2);
plot(time,yawRate3_m,"-.","LineWidth",2);% 
% plot(yaw1_m,"b--");
% plot(q_est(:,3),"g.-");
legend("Yaw Rate 1","Yaw Rate 3");
grid on;
title(['Yaw Rate, RWA = ',num2str(rwa)]);
xlabel('time');
ylabel('Yaw Rate [deg/s]');
set(gca,'FontName','Times New Rome','FontSize',15);
savefig(['result/' filename '_YawRate_RWA']);
%% Axis Position Cal
l = [8.9 7.4 1.5 7.9 1.5 7.4 1.5 0];
ld = diag(l,0);
l_2 = 9.4;
l_3 = 8.9;
x_J1 = q_est(:,1);
y_J1 = q_est(:,2);
x_J2 = x_J1 - l_2.*cos(q_est(:,4));
y_J2 = y_J1 - l_2.*sin(q_est(:,4));
x_A8 = x_J2 - l_3.*cos(q_est(:,5));
y_A8 = y_J2 - l_3.*sin(q_est(:,5));
x0 = [x_J1 x_J1 x_J1 x_J2 x_J2 x_A8 x_A8 x_A8];
y0 = [y_J1 y_J1 y_J1 y_J2 y_J2 y_A8 y_A8 y_A8];
yaw = [q_est(:,3) q_est(:,3) q_est(:,3) q_est(:,4) q_est(:,4) q_est(:,5) q_est(:,5) q_est(:,5) ];
x = x0 + cos(yaw)*ld;
y = y0 + sin(yaw)*ld;

% %% Offtrack Calculate
% % 
% x_ref = x(:,1);
% y_ref = y(:,1);
% traj = [x_ref y_ref];
% % % point = [x_point y_point];
% section = diff(traj); % 各路径段向量
% sectionLength = hypot(section(:,1),section(:,2));% 各路径段长度
% s_traj = zeros(size(sectionLength));
% for i = 1:size(section,1)
%     s_traj(i) = sum(sectionLength(1:i));
% end
% s_traj = [0;s_traj];
% 
% s = zeros(size(x));
% offtrack = zeros(size(y));
% s_J1 = zeros(size(x_J1));
% s_J2 = zeros(size(x_J2));
% offtrack_J1 = zeros(size(y_J1));
% offtrack_J2 = zeros(size(y_J2));
% 
% temp_start = 1;
% temp_end = data_len;
% [s_J1(temp_start:temp_end),offtrack_J1(temp_start:temp_end)] = segment_project(x_ref(temp_start:temp_end),y_ref(temp_start:temp_end),s_traj(temp_start:temp_end),x_J1(temp_start:temp_end),y_J1(temp_start:temp_end));
% [s_J2(temp_start:temp_end),offtrack_J2(temp_start:temp_end)] = segment_project(x_ref(temp_start:temp_end),y_ref(temp_start:temp_end),s_traj(temp_start:temp_end),x_J2(temp_start:temp_end),y_J2(temp_start:temp_end));
% for j = 1:8
%     [s(temp_start:temp_end,j),offtrack(temp_start:temp_end,j)] = segment_project(x_ref(temp_start:temp_end),y_ref(temp_start:temp_end),s_traj(temp_start:temp_end),x(temp_start:temp_end,j),y(temp_start:temp_end,j));
% end  
% %%
% % linetype = ["-" "--" ":" "-." "--." ":." ];%"-." "m-x"];
% linetype = ["-" "--" ":" "-." "--." ":." ];
% figure('Name','Offtrack');
% hold on;
% % axis equal;
% % temp_start = 5800;
% % temp_end = data_len;
% % temp_start = 4500;
% % temp_start = 2500;% 1919
% % temp_start = 5000;% 1853
% % temp_start = 2700;% 1854
% % temp_start = 5700;% 2108
% % temp_start = 4500;% 2202
% temp_start = 2500;% 2315
% % temp_start = 3500;% 1928
% 
% 
% 
% temp_end = data_len;
% 
% 
% temp = [1 3 4 5 6 8];
% for i = 1:6
%     plot(s(temp_start:temp_end,temp(i)),offtrack(temp_start:temp_end,temp(i)),linetype(i),'MarkerIndices',1:100:data_len,'Linewidth', 1);
% end
% 
% plot(s_J1(temp_start:temp_end),offtrack_J1(temp_start:temp_end),'b--','Linewidth', 2);
% plot(s_J2(temp_start:temp_end),offtrack_J2(temp_start:temp_end),'g:','Linewidth', 2);
% % plot(s(temp_start:temp_end,1),offtrack(temp_start:temp_end,1),'Linewidth', 2);
% % plot(s(temp_start:temp_end,8),offtrack(temp_start:temp_end,8),'Linewidth', 2);
% % legend("Joint 1","Joint 2","Axle 1","Axle 8");
% 
% % legend("Axle 1","Axle 2","Axle 3","Axle 4","Axle 5","Axle 6","Axle 7","Axle 8","Joint 1","Joint 2");
% legend("Axle 1","Axle 3","Axle 4","Axle 5","Axle 6","Axle 8","Joint 1","Joint 2");
% 
% grid on;
% title('Offtrack');
% xlabel('s [m]');
% ylabel('Offtrack [m]');
% set(gcf,'Position',[100,100,800,420]);
% set(gca,'FontName','Times New Rome','FontSize',15);
% 
% savefig(['result/' filename '_Offtrack']);

%%
% 

figure('Name','traj');
hold on;
axis equal;
traj_a_1 = plot(x(:,1),y(:,1),"k","LineWidth",2);
traj_a_8 = plot(x(:,8),y(:,8),"g--","LineWidth",2);
% plot(x_m,y_m,"b.");
% plot(x_m3,y_m3,"r.");
traj_j_1 = plot(x_J1,y_J1,"b:","LineWidth",2);
traj_j_2 = plot(x_J2,y_J2,"r-.","LineWidth",2);
% 
% dtRows_j_1 = [dataTipTextRow("Dstance",s_J1'),dataTipTextRow("Offtrack",offtrack_J1')];
% traj_j_1.DataTipTemplate.DataTipRows(end+1:end+2) = dtRows_j_1;
% dtRows_j_2 = [dataTipTextRow("Dstance",s_J2'),dataTipTextRow("Offtrack",offtrack_J2')];
% traj_j_2.DataTipTemplate.DataTipRows(end+1:end+2) = dtRows_j_2;
% dtRows_a_1 = [dataTipTextRow("Dstance",s(:,1)'),dataTipTextRow("Offtrack",offtrack(:,1)')];
% traj_a_1.DataTipTemplate.DataTipRows(end+1:end+2) = dtRows_a_1;
% dtRows_a_8 = [dataTipTextRow("Dstance",s(:,8)'),dataTipTextRow("Offtrack",offtrack(:,8)')];
% traj_a_8.DataTipTemplate.DataTipRows(end+1:end+2) = dtRows_a_8;
% 
% plot(q_est(:,1),q_est(:,2),"g.");
% legend("measurement J1","measurement J2","KF");
legend("Axle 1","Axle 8","J1","J2");

grid on;
title('Trajectory of Joints');
xlabel('X [m]');
ylabel('Y [m]');
grid minor
set(gca,'FontName','Times New Rome','FontSize',15);
savefig(['result/' filename '_Traj']);

%%
angle_off_set = 0;
temp_offset = 0;
figure('Name','yaw');
subplot(2,1,1);
hold on;
plot(time,wrap_angle(yaw1_m/deg + angle_off_set)*deg - angle_off_set*deg,"b--");
plot(time,wrap_angle(q_est(:,3) + angle_off_set)*deg - angle_off_set*deg - temp_offset,"g.-");
% 
% plot(yaw1_m,"b--");
% plot(q_est(:,3),"g.-");
legend("measurement","KF");
grid on;
title('Yaw1');
xlabel('time');
ylabel('Yaw 1 [deg]');
subplot(2,1,2);
hold on;
plot(time,wrap_angle(yaw3_m/deg + angle_off_set)*deg - angle_off_set*deg,"b--");
plot(time,wrap_angle(q_est(:,5) + angle_off_set)*deg - angle_off_set*deg +temp_offset,"g.-");
legend("measurement","KF");
grid on;
title('Yaw3');
xlabel('time');
ylabel('Yaw 3 [deg]');
savefig(['result/' filename '_Yaw']);


figure('Name','art');
subplot(2,1,1);
hold on;
% angle_off_set = 0;
plot(time,wrap_angle(art1_m/deg)*deg,"b--");
plot(time,wrap_angle(q_est(:,3) - q_est(:,4) )*deg,"g.-");
legend("measurement","KF");
grid on;
title('Art1');
xlabel('time');
ylabel('art 1 [deg]');
subplot(2,1,2);
hold on;
% angle_off_set = 0;
plot(time,wrap_angle(art2_m/deg)*deg,"b--");
plot(time,wrap_angle(q_est(:,4) - q_est(:,5) )*deg,"g.-");
legend("measurement","KF");
grid on;
title('Art2');
xlabel('time');
ylabel('art 2 [deg]');
savefig(['result/' filename '_Art']);



%%

gps_vn1 = diff(y_J1)./diff(time);
gps_vn3 = diff(y_J2)./diff(time);
gps_ve1 = diff(x_J1)./diff(time);
gps_ve3 = diff(x_J2)./diff(time);
gps_vn1 = [gps_vn1(1) ; gps_vn1];
gps_ve1 = [gps_ve1(1) ; gps_ve1];


figure('Name','speed');
subplot(3,1,1);
hold on;
plot(time,spdEast1_m,"b--");
plot(time,qd_est(:,1),"g.-");
% plot(time,gps_ve1,"r--");

legend("measurement","KF");
grid on;
title('spdEast1');
xlabel('time');
ylabel('spdEast1 [m/s]');
subplot(3,1,2);
hold on;
plot(time,spdNorth1_m,"b--");
plot(time,qd_est(:,2),"g.-");
% plot(time,gps_vn1,"r--");

legend("measurement","KF");
grid on;
title('spdNorth1');
xlabel('time');
ylabel('spdNorth1 [m/s]');
subplot(3,1,3);
hold on;
plot(time,spd1_m,"b--");
legend("measurement");
grid on;
title('spd1');
xlabel('time');
ylabel('spd1 [m/s]');

savefig(['result/' filename '_Speed']);


figure('Name','yawrate');
subplot(2,1,1);
hold on;
plot(time,yawRate1_m,"b--");
plot(time,qd_est(:,3)*deg,"g.-");
legend("measurement","KF");
grid on;
title('yawRate1');
xlabel('time');
ylabel('yawRate1 [deg/s]');
subplot(2,1,2);
hold on;
% plot(time,yawRate3_m/deg,"b--");
plot(time,yawRate3_m,"b--");

plot(time,qd_est(:,5)*deg,"g.-");
legend("measurement","KF");
grid on;
title('yawRate3');
xlabel('time');
ylabel('yawRate3 [deg/s]');
savefig(['result/' filename '_Yawrate']);

%%

figure('Name','steer and steer cmd');
subplot(3,2,1);
hold on;
plot(time,steer_m(:,1),"b--");
% plot(time,steer_cmd(:,1),"g.-");
% plot(time,ff_error_debug(:,1)*deg,"r.-");
% plot(time,fb_str(:,1)*deg,"r.-");
title('Steer 1');
ylabel('Steer [deg]');
legend('real steer');
xlabel('time');
subplot(3,2,2);
hold on;
plot(time,steer_m(:,2),"b--");
plot(time,steer_cmd(:,2),"g.-");
plot(time,wrap_angle(ff_error_debug(:,2))*deg,"r.-");
% plot(time,fb_str(:,2)*deg,"k.-");
title('Steer 2');
ylabel('Steer [deg]');
legend('real steer','steer command','feedforward steer','feedback steer');
xlabel('time');
subplot(3,2,3);
hold on;
plot(time,steer_m(:,3),"b--");
plot(time,steer_cmd(:,3),"g.-");
plot(time,wrap_angle(ff_error_debug(:,3))*deg,"r.-");
% plot(time,fb_str(:,3)*deg,"k.-");
title('Steer 3');
ylabel('Steer [deg]');
legend('real steer','steer command','feedforward steer','feedback steer');
xlabel('time');
subplot(3,2,4);
hold on;
plot(time,steer_m(:,4),"b--");
plot(time,steer_cmd(:,4),"g.-");
plot(time,wrap_angle(ff_error_debug(:,4))*deg,"r.-");
% plot(time,fb_str(:,4)*deg,"k.-");
title('Steer 4');
ylabel('Steer [deg]');
legend('real steer','steer command','feedforward steer','feedback steer');
xlabel('time');
subplot(3,2,5);
hold on;
plot(time,steer_m(:,5),"b--");
plot(time,steer_cmd(:,5),"g.-");
plot(time,wrap_angle(ff_error_debug(:,5))*deg,"r.-");
% plot(time,fb_str(:,5)*deg,"k.-");
title('Steer 5');
ylabel('Steer [deg]');
legend('real steer','steer command','feedforward steer','feedback steer');
xlabel('time');
subplot(3,2,6);
hold on;
plot(time,steer_m(:,6),"b--");
plot(time,steer_cmd(:,6),"g.-");
plot(time,wrap_angle(ff_error_debug(:,6))*deg,"r.-");
% plot(time,fb_str(:,6)*deg,"k.-");
title('Steer 6');
ylabel('Steer [deg]');
legend('real steer','steer command','feedforward steer','feedback steer');
xlabel('time');

savefig(['result/' filename '_SteerCmd']);

% figure('Name','steer and steer cmd');
% subplot(3,2,1);
% hold on;
% plot(time,steer_m(:,1),"b--");
% % plot(time,steer_cmd(:,1),"g.-");
% % plot(time,ff_error_debug(:,1)*deg,"r.-");
% % plot(time,fb_str(:,1)*deg,"r.-");
% title('Steer 1');
% ylabel('Steer [deg]');
% legend('real steer');
% xlabel('time');
% subplot(3,2,2);
% hold on;
% plot(time,steer_m(:,2),"b--");
% plot(time,steer_cmd(:,2),"g.-");
% plot(time,ff_error_debug(:,2)*deg,"r.-");
% plot(time,fb_error_debug(:,2)*deg,"m.-");
% plot(time,fb_str(:,2)*deg,"k.-");
% title('Steer 2');
% ylabel('Steer [deg]');
% legend('real steer','steer command','feedforward steer','feedback error','feedback steer');
% xlabel('time');
% subplot(3,2,3);
% hold on;
% plot(time,steer_m(:,3),"b--");
% plot(time,steer_cmd(:,3),"g.-");
% plot(time,ff_error_debug(:,3)*deg,"r.-");
% plot(time,fb_error_debug(:,3)*deg,"m.-");
% plot(time,fb_str(:,3)*deg,"k.-");
% title('Steer 3');
% ylabel('Steer [deg]');
% legend('real steer','steer command','feedforward steer','feedback error','feedback steer');
% xlabel('time');
% subplot(3,2,4);
% hold on;
% plot(time,steer_m(:,4),"b--");
% plot(time,steer_cmd(:,4),"g.-");
% plot(time,ff_error_debug(:,4)*deg,"r.-");
% plot(time,fb_error_debug(:,4)*deg,"m.-");
% plot(time,fb_str(:,4)*deg,"k.-");
% title('Steer 4');
% ylabel('Steer [deg]');
% legend('real steer','steer command','feedforward steer','feedback error','feedback steer');
% xlabel('time');
% subplot(3,2,5);
% hold on;
% plot(time,steer_m(:,5),"b--");
% plot(time,steer_cmd(:,5),"g.-");
% plot(time,ff_error_debug(:,5)*deg,"r.-");
% plot(time,fb_error_debug(:,5)*deg,"m.-");
% plot(time,fb_str(:,5)*deg,"k.-");
% title('Steer 5');
% ylabel('Steer [deg]');
% legend('real steer','steer command','feedforward steer','feedback error','feedback steer');
% xlabel('time');
% subplot(3,2,6);
% hold on;
% plot(time,steer_m(:,6),"b--");
% plot(time,steer_cmd(:,6),"g.-");
% plot(time,ff_error_debug(:,6)*deg,"r.-");
% plot(time,fb_error_debug(:,6)*deg,"m.-");
% plot(time,fb_str(:,6)*deg,"k.-");
% title('Steer 6');
% ylabel('Steer [deg]');
% legend('real steer','steer command','feedforward steer','feedback error','feedback steer');
% xlabel('time');
% savefig(['result/' filename '_SteerCmd']);
% figure();
% subplot(3,2,1);
% plot(time,intsteercmd(:,1),"k.-");
% 
fb_start= 3000;
fb_end = data_len - 1000;
figure('name','feedback error')
subplot(3,2,2);
% plot(time,intsteercmd(:,2),"k.-");
plot(time(fb_start:fb_end),fb_error_debug(fb_start:fb_end,2)*deg,"m.-");
title('axle 2');
legend('fb error');
ylabel('feedback error [deg]');
xlabel('time');

subplot(3,2,3);
% plot(time,intsteercmd(:,3),"k.-");
plot(time(fb_start:fb_end),fb_error_debug(fb_start:fb_end,3)*deg,"m.-");
title('axle 3');
legend('fb error');
ylabel('feedback error [deg]');
xlabel('time')

subplot(3,2,4);
% plot(time,intsteercmd(:,4),"k.-");
plot(time(fb_start:fb_end),fb_error_debug(fb_start:fb_end,4)*deg,"m.-");
title('axle 4');
legend('fb error');
ylabel('feedback error [deg]');
xlabel('time')

subplot(3,2,5);
% plot(time,intsteercmd(:,5),"k.-");
plot(time(fb_start:fb_end),fb_error_debug(fb_start:fb_end,5)*deg,"m.-");
title('axle 5');
legend('fb error');
ylabel('feedback error [deg]');
xlabel('time')

subplot(3,2,6);
% plot(time,intsteercmd(:,6),"k.-");
plot(time(fb_start:fb_end),fb_error_debug(fb_start:fb_end,6)*deg,"m.-");
title('axle 6');
legend('fb error');
ylabel('feedback error [deg]');
xlabel('time')

savefig(['result/' filename '_feedbackerr']);


offset = pi;

figure('name','Angref');
subplot(3,2,2);
hold on;
plot(time,wrap_angle(ANGref(:,2) - fb_error_debug(:,2)+offset)*deg-offset*deg,"k.-");
plot(time,wrap_angle(ANGref(:,2)+offset)*deg-offset*deg,"m.-");
legend('speed angle','ANGref');
ylabel('angle [deg]');
xlabel('time')

subplot(3,2,3);
hold on;
plot(time,wrap_angle(ANGref(:,3) - fb_error_debug(:,3)+offset)*deg-offset*deg,"k.-");
plot(time,wrap_angle(ANGref(:,3)+offset)*deg-offset*deg,"m.-");
legend('speed angle','ANGref');
ylabel('angle [deg]');
xlabel('time')

subplot(3,2,4);
hold on;
plot(time,wrap_angle(ANGref(:,4) - fb_error_debug(:,4)+offset)*deg-offset*deg,"k.-");
plot(time,wrap_angle(ANGref(:,4)+offset)*deg-offset*deg,"m.-");
legend('speed angle','ANGref');
ylabel('angle [deg]');
xlabel('time')

subplot(3,2,5);
hold on;
plot(time,wrap_angle(ANGref(:,5) - fb_error_debug(:,5)+offset)*deg-offset*deg,"k.-");
plot(time,wrap_angle(ANGref(:,5)+offset)*deg-offset*deg,"m.-");
legend('speed angle','ANGref');
ylabel('angle [deg]');
xlabel('time')

subplot(3,2,6);
hold on;
plot(time,wrap_angle(ANGref(:,6) - fb_error_debug(:,6)+offset)*deg-offset*deg,"k.-");
plot(time,wrap_angle(ANGref(:,6)+offset)*deg-offset*deg,"m.-");
legend('speed angle','ANGref');
ylabel('angle [deg]');
xlabel('time')

%% 投影函数 整段
function [s,offtrack] = segment_project(x_ref,y_ref,s_traj,x,y)
    s = zeros(size(x));
    offtrack = zeros(size(y));
    for i = 1:size(x,1)
        index_start = 1;
        index_end = size(x,1);
        [s(i),offtrack(i)] = cart2frenet(x_ref(index_start:index_end),y_ref(index_start:index_end),s_traj(index_start:index_end),x(i),y(i));
    end
end

%% 投影函数 单点

function [s, d, err_all] = cart2frenet(x_ref, y_ref, s_traj, x_point, y_point)
traj = [x_ref y_ref];
section = diff(traj);
vector = [x_point - x_ref(1:end-1), y_point - y_ref(1:end-1)];
sectionLength = sum(section.^2, 2);
t = dot(vector, section, 2) ./ sectionLength;
t = max(0.0, min(t, 1.0));
point_proj = traj(1:end-1, :) + t .* section;
err_all = hypot(point_proj(:, 1) - x_point, point_proj(:, 2) - y_point);
[err, index] = min(err_all);
point_proj_out = point_proj(index, :);
s = s_traj(index) + dot(vector(index, :), section(index, :)) / sqrt(sectionLength(index));
s = max(0, s);
d = hypot(x_point - point_proj_out(1), y_point - point_proj_out(2));
dir = (x_point - point_proj_out(1)) * section(index, 2) - (y_point - point_proj_out(2)) * section(index, 1);
d = d * sign(dir);
end
% 
% %% 投影函数 单点
% 
% function [s,d,err_all] = cart2frenet(x_ref,y_ref,s_traj,x_point,y_point)
% % function [s,d] = cart2frenet(traj,section,sectionLength,x_point,y_point)
% % x_ref = x_ref(index_start:index_end);
% % y_ref = y_ref(index_start:index_end);
% traj = [x_ref y_ref];
% % point = [x_point y_point];
% section = diff(traj); % 各路径段向量
% vector = [x_point - x_ref y_point - y_ref];% 被投影点到各路径点的距离
% vector = vector(1:end-1,:);
% sectionLength = sum(section.*section,2);% 各路径段长度平方
% t = (vector(:,1).*section(:,1)+vector(:,2).*section(:,2)) ./ sectionLength;% 被投影点在各路径的投影与相应路径段长度的比值
% point_proj = zeros(size(section));
% 
% % t_temp = t(end);
% % t(t < 0.0) = 0;
% % t(t > 1.0) = ;
% 
% for i = 1:size(t,1)
%     if ((t(i) > 0.0 && t(i) < 1.0) || (t(i) >= 1.0 && i == size(t,1)) || (t(i) <= 0.0 && i <= 20))
%         point_proj(i,:) = traj(i,:) + t(i)*section(i,:);
%     elseif(t(i) <= 0.0)
%         point_proj(i,:) = traj(i,:);
%     else%(t(i) >= 1.0)
%         point_proj(i,:) = traj(i+1,:);
%     end
% end
% err_all = hypot(point_proj(:,1) - x_point, point_proj(:,2) - y_point);
% [err,index] = min(err_all);
% point_proj_out = [point_proj(index,1),point_proj(index,2)];
% % s = s_traj(index_start+index) + vector(index,1).*section(index,1)+vector(index,2).*section(index,2);
% s = s_traj(index) + vector(index,1).*section(index,1)+vector(index,2).*section(index,2)/sqrt(sectionLength(index));
% if( s <0)
%     s =0;
%     d = 0;
% else
%     d = hypot(x_point - point_proj_out(1), y_point - point_proj_out(2));
%     dir = (x_point - point_proj_out(1))*section(index,2) - (y_point - point_proj_out(2))*section(index,1);
%     d = d*sign(dir);
% end
% % vector = point - section;
% 
% end


%%
function [x_out, y_out] = LonLat2XY(lat,lon,origin_lat,origin_lon)
origin_y = 0;
origin_x = 0;
a = 6378137; 
b = 6356752.3142; 
f = (a-b)/a;
lat_rad = deg2rad(lat);
lon_rad = deg2rad(lon);
origin_lat_rad = deg2rad(origin_lat);
origin_lon_rad = deg2rad(origin_lon);
%  https://www.mathworks.com/help/aerotbx/ug/lla2flat.html
RN = a/(sqrt(1-(2*f-f^2)*(sin(origin_lat_rad))^2));
RM = (RN*(1-(2*f-f^2)))/((1-(2*f-f^2)*(sin(origin_lat_rad))^2));
dlat = lat_rad - origin_lat_rad;
dlon = lon_rad - origin_lon_rad;
dy = dlat/atan(1/RM);
dx = dlon/atan(1/(RN*cos(origin_lat_rad)));
x_out = dx + origin_x;
y_out = dy + origin_y;
end