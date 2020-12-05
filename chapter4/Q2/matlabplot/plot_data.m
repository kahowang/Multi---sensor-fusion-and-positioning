
delimiterIn = ',';
headerlinesIn = 0;
gt = importdata('/../gnss-ins-sim/imu_data/data7/gt.txt', delimiterIn, headerlinesIn);
pose = importdata('/../gnss-ins-sim/imu_data/data7/pose.txt', delimiterIn, headerlinesIn);
for i = 1:size(pose,1)
    temp = pose(i,8:10).*pose(i,8:10);
    temp = sum(temp);
    temp = temp^0.5;
    phi = pose(i,8:10)/temp*(-1.0)*(2*pi-temp);
    a = pose(i,8:10) - gt(i,8:10);
    a = sum(a.*a);
    b = phi - gt(i,8:10);
    b = sum(b.*b);
    if b < a
        pose(i,8:10) = phi;
    end
end
sv = importdata('/../gnss-ins-sim/imu_data/data7/sv.txt', delimiterIn, headerlinesIn);
sv = sv';

figure(1)
subplot(3,2,1)
plot(gt(:,1),gt(:,2),'r--',gt(:,1),gt(:,3),'g--',gt(:,1),gt(:,4),'b--',...
    pose(:,1),pose(:,2),'r',pose(:,1),pose(:,3),'g',pose(:,1),pose(:,4),'b','LineWidth',2);
xlabel('time (sec)','FontSize',12,'FontWeight','bold')
ylabel('position (m)','FontSize',12,'FontWeight','bold')
legend('x','y','z','Location','Best');
%legend('x_g_t','y_g_t','z_g_t','x','y','z','Location','SouthWest','FontSize',15);
set(gca,'FontSize',12);

subplot(3,2,2)
plot(gt(:,1),gt(:,5),'r--',gt(:,1),gt(:,6),'g--',gt(:,1),gt(:,7),'b--',...
    pose(:,1),pose(:,5),'r',pose(:,1),pose(:,6),'g',pose(:,1),pose(:,7),'b','LineWidth',2);
xlabel('time (sec)','FontSize',12,'FontWeight','bold')
ylabel('velocity (m/s)','FontSize',12,'FontWeight','bold')
legend('vx','vy','vz','Location','Best');
set(gca,'FontSize',12);

subplot(3,2,3)
plot(gt(:,1),gt(:,8),'r--',gt(:,1),gt(:,9),'g--',gt(:,1),gt(:,10),'b--',...
    pose(:,1),pose(:,8),'r',pose(:,1),pose(:,9),'g',pose(:,1),pose(:,10),'b','LineWidth',2);
xlabel('time (sec)','FontSize',12,'FontWeight','bold')
ylabel('q (rad)','FontSize',12,'FontWeight','bold')
legend('qx','qy','qz','Location','Best');
set(gca,'FontSize',12);

subplot(3,2,4)
plot(gt(:,1),gt(:,11),'r--',gt(:,1),gt(:,12),'g--',gt(:,1),gt(:,13),'b--',...
    pose(:,1),pose(:,11),'r',pose(:,1),pose(:,12),'g',pose(:,1),pose(:,13),'b','LineWidth',2);
xlabel('time (sec)','FontSize',12,'FontWeight','bold')
ylabel('bg (rad/s)','FontSize',12,'FontWeight','bold')
legend('bgx','bgy','bgz','Location','Best');
set(gca,'FontSize',12);

subplot(3,2,5)
plot(gt(:,1),gt(:,14),'r--',gt(:,1),gt(:,15),'g--',gt(:,1),gt(:,16),'b--',...
    pose(:,1),pose(:,14),'r',pose(:,1),pose(:,15),'g',pose(:,1),pose(:,16),'b','LineWidth',2);
xlabel('time (sec)','FontSize',12,'FontWeight','bold')
ylabel('ba (m/s^2)','FontSize',12,'FontWeight','bold')
legend('bax','bay','baz','Location','Best');
set(gca,'FontSize',12);

