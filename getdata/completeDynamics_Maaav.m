% Copyright 2024, All Rights Reserved
% Code by Chao Liang
% For Paper, "xxxxxx"
% by Chao Liang

n = 10;

%% True dynamics (simulated) and measured states

ty1 = time1;
y1 = qqd1;
ty2 = time2;
y2 = qqd2;
ty3 = time3;
y3 = qqd3;
ty4 = time4;
y4 = qqd4;
ty5 = time5;
y5 = qqd5;
ty6 = time6;
y6 = qqd6;

%% Multi acticulated all axle vehicle system, lagrange

if dataGenerateMode == 1
    tx1 = time1;
    x1 = Maaav_single(time1, y1, avp, u1, n);
    tx2 = time2;
    x2 = Maaav_single(time2, y2, avp, u2, n);
    tx3 = time3;
    x3 = Maaav_single(time3, y3, avp, u3, n);
    tx4 = time4;
    x4 = Maaav_single(time4, y4, avp, u4, n);
    tx5 = time5;
    x5 = Maaav_single(time5, y5, avp, u5, n);
    tx6 = time6;
    x6 = Maaav_single(time6, y6, avp, u6, n);
end
if dataGenerateMode == 2
    tx1 = time1;
    x1 = Maaav_iter(time1, y1, avp, u1, n);
    tx2 = time2;
    x2 = Maaav_iter(time2, y2, avp, u2, n);
    tx3 = time3;
    x3 = Maaav_iter(time3, y3, avp, u3, n);
    tx4 = time4;
    x4 = Maaav_iter(time4, y4, avp, u4, n);
    tx5 = time5;
    x5 = Maaav_iter(time5, y5, avp, u5, n);
    tx6 = time6;
    x6 = Maaav_iter(time6, y6, avp, u6, n);
end

% figure,
% for i = 1:n
%     subplot(2,n/2,i)
%     plot(tx1,x1(:,i),'b--','LineWidth',1), hold on
%     plot(tx1,y1(:,i),'k','LineWidth',1)
%     axis equal
%     legend('ideal', 'real')
% end
% sgtitle('state data 1')
% 
% figure, 
% plot(x1(:,1), x1(:,2), 'b--','LineWidth',1), hold on
% plot(y1(:,1), y1(:,2), 'k','LineWidth',1), hold on
% axis equal
% legend('ideal', 'real')
% sgtitle('trajectory of ideal and real system 1')
% 
% figure,
% for i = 1:n
%     subplot(2,n/2,i)
%     plot(tx2,x2(:,i),'b--','LineWidth',1), hold on
%     plot(tx2,y2(:,i),'k','LineWidth',1)
%     axis equal
%     legend('ideal', 'real')
% end
% sgtitle('state data 2')
% 
% figure, 
% plot(x2(:,1), x2(:,2), 'b--','LineWidth',1), hold on
% plot(y2(:,1), y2(:,2), 'k','LineWidth',1), hold on
% axis equal
% legend('ideal', 'real')
% sgtitle('trajectory of ideal and real system 2')
% 
% figure,
% for i = 1:n
%     subplot(2,n/2,i)
%     plot(tx3,x3(:,i),'b--','LineWidth',1), hold on
%     plot(tx3,y3(:,i),'k','LineWidth',1)
%     axis equal
%     legend('ideal', 'real')
% end
% sgtitle('state data 3')
% 
% figure, 
% plot(x3(:,1), x3(:,2), 'b--','LineWidth',1), hold on
% plot(y3(:,1), y3(:,2), 'k','LineWidth',1), hold on
% axis equal
% legend('ideal', 'real')
% sgtitle('trajectory of ideal and real system 3')
%
% figure,
% for i = 1:n
%     subplot(2,n/2,i)
%     plot(tx4,x4(:,i),'b--','LineWidth',1), hold on
%     plot(tx4,y4(:,i),'k','LineWidth',1)
%     axis equal
%     legend('ideal', 'real')
% end
% sgtitle('state data 4')
% 
% figure, 
% plot(x4(:,1), x4(:,2), 'b--','LineWidth',1), hold on
% plot(y4(:,1), y4(:,2), 'k','LineWidth',1), hold on
% axis equal
% legend('ideal', 'real')
% sgtitle('trajectory of ideal and real system 4')
%
% figure,
% for i = 1:n
%     subplot(2,n/2,i)
%     plot(tx5,x5(:,i),'b--','LineWidth',1), hold on
%     plot(tx5,y5(:,i),'k','LineWidth',1)
%     axis equal
%     legend('ideal', 'real')
% end
% sgtitle('state data 5')
% 
% figure, 
% plot(x5(:,1), x5(:,2), 'b--','LineWidth',1), hold on
% plot(y5(:,1), y5(:,2), 'k','LineWidth',1), hold on
% axis equal
% legend('ideal', 'real')
% sgtitle('trajectory of ideal and real system 5')
%
% figure,
% for i = 1:n
%     subplot(2,n/2,i)
%     plot(tx6,x6(:,i),'b--','LineWidth',1), hold on
%     plot(tx6,y6(:,i),'k','LineWidth',1)
%     axis equal
%     legend('ideal', 'real')
% end
% sgtitle('state data 6')
% 
% figure, 
% plot(x6(:,1), x6(:,2), 'b--','LineWidth',1), hold on
% plot(y6(:,1), y6(:,2), 'k','LineWidth',1), hold on
% axis equal
% legend('ideal', 'real')
% sgtitle('trajectory of ideal and real system 6')

%% compute Derivative

for i = 1:n
    dy1(:,i) = gradient(y1(:,i), dt);
    dy2(:,i) = gradient(y2(:,i), dt);
    dy3(:,i) = gradient(y3(:,i), dt);
    dy4(:,i) = gradient(y4(:,i), dt);
    dy5(:,i) = gradient(y5(:,i), dt);
    dy6(:,i) = gradient(y6(:,i), dt);
end

% figure,
% for i = 1:n
%     subplot(2,n/2,i)
%     plot(ty1, dy1(:,i),'k','Linewidth',[2]);
%     hold on;
%     plot(ty1, y1(:,i),'r--','Linewidth',[2])
%     legend('dy','y')
%     title('Detrivative of real Data 1')
% end
% 
% figure,
% for i = 1:n
%     subplot(2,n/2,i)
%     plot(ty2, dy2(:,i),'k','Linewidth',[2]);
%     hold on;
%     plot(ty2, y2(:,i),'r--','Linewidth',[2])
%     legend('dy','y')
%     title('Detrivative of real Data 2')
% end
% 
% figure,
% for i = 1:n
%     subplot(2,n/2,i)
%     plot(ty3, dy3(:,i),'k','Linewidth',[2]);
%     hold on;
%     plot(ty3, y3(:,i),'r--','Linewidth',[2])
%     legend('dy','y')
%     title('Detrivative of real Data 3')
% end
% 
% figure,
% for i = 1:n
%     subplot(2,n/2,i)
%     plot(ty4, dy4(:,i),'k','Linewidth',[2]);
%     hold on;
%     plot(ty4, y4(:,i),'r--','Linewidth',[2])
%     legend('dy','y')
%     title('Detrivative of real Data 4')
% end
% 
% figure,
% for i = 1:n
%     subplot(2,n/2,i)
%     plot(ty5, dy5(:,i),'k','Linewidth',[2]);
%     hold on;
%     plot(ty5, y5(:,i),'r--','Linewidth',[2])
%     legend('dy','y')
%     title('Detrivative of real Data 5')
% end
% 
% figure,
% for i = 1:n
%     subplot(2,n/2,i)
%     plot(ty6, dy6(:,i),'k','Linewidth',[2]);
%     hold on;
%     plot(ty6, y6(:,i),'r--','Linewidth',[2])
%     legend('dy','y')
%     title('Detrivative of real Data 6')
% end

X = [x1; x2; x3; x4; x5; x6;];
Y = [y1; y2; y3; y4; y5; y6;];
dY = [dy1; dy2; dy3; dy4; dy5; dy6;];
U = [u1; u2; u3; u4; u5; u6;]; 
T = [time1; time2+time1(end)+dt; time3+time1(end)+time2(end)+dt+dt;
     time4+time1(end)+time2(end)+time3(end)+dt+dt+dt;
     time5+time1(end)+time2(end)+time3(end)+time4(end)+dt+dt+dt+dt;
     time6+time1(end)+time2(end)+time3(end)+time4(end)+time5(end)+dt+dt+dt+dt+dt;];
