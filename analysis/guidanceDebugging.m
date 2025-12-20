clear; clc; close all;

load("guidanceLogs.mat");
aim = guidanceLog_aim;
sim = guidanceLog_sim;
figure
plot(sim(:,1),sim(:,4),'bo:',LineWidth=1)
hold on
plot(aim(:,1),aim(:,4),'rx:',LineWidth=1)
xlabel('Sim Time (s)')
ylabel('Closing Velocity (m/s)')
legend('Simple Model','Full Model','location','best')
title('Closing Velocity Solution')
grid on

figure
subplot(3,1,1)
plot(sim(:,1),sim(:,6),'bo:',LineWidth=1)
hold on
plot(aim(:,1),aim(:,6),'rx:',LineWidth=1)
grid on

title('LOS Vector Solution')
ylabel('X')
subplot(3,1,2)
plot(sim(:,1),sim(:,7),'bo:',LineWidth=1)
hold on
plot(aim(:,1),aim(:,7),'rx:',LineWidth=1)
grid on

ylabel('Y')
subplot(3,1,3)
plot(sim(:,1),sim(:,8),'bo:',LineWidth=1)
hold on
plot(aim(:,1),aim(:,8),'rx:',LineWidth=1)
grid on

ylabel('Z')
xlabel('Sim Time(s)')

figure
subplot(3,1,1)
plot(sim(:,1),sim(:,9),'bo:',LineWidth=1)
hold on
plot(aim(:,1),aim(:,9),'rx:',LineWidth=1)
grid on

title('Inertial LOS Rates')
ylabel('X')
subplot(3,1,2)
plot(sim(:,1),sim(:,10),'bo:',LineWidth=1)
grid on

hold on
plot(aim(:,1),aim(:,10),'rx:',LineWidth=1)
ylabel('Y')
subplot(3,1,3)
plot(sim(:,1),sim(:,11),'bo:',LineWidth=1)
grid on

hold on
plot(aim(:,1),aim(:,11),'rx:',LineWidth=1)
ylabel('Z')
xlabel('Sim Time(s)')
