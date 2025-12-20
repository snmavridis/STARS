%close all
global hpbw

N = 12;
angmax = 45;
%-------------------------------
theta_deg = -angmax:0.2:angmax;
theta = pi/180*theta_deg';

Null2Null_bw = 2/N * 180/pi;
hpbw = 0.886*Null2Null_bw;
hpbw_onesided = hpbw/2;

Nend = (N-1)/2;
phi = pi*sin(theta)*(-Nend:Nend);

s = exp(1j*phi');

Sum = sum(s);
Dif = 1j*(sum(s(1:N/2,:))-sum(s(N/2+1:N,:)));

rho = real(Dif.*conj(Sum))./abs(Sum).^2;
rho2 = real(Dif)./abs(Sum);

peak = max(Sum);
Sum = abs(Sum/peak);
Dif = real(Dif/peak);
Sum2D = Sum.' * Sum;

figure(1)
plot(theta_deg,Sum,'LineWidth',2)
grid on
xlabel('Angle (deg)')
ylabel('Magnitude')
title('Sum Pattern')

figure(2)
plot(theta_deg,20*log10(Sum),'LineWidth',2)
grid on
xlabel('Angle (deg)')
ylabel('dB')
axis([-angmax angmax -30 0])
title('Normalized Sum Power Pattern')

figure(3)
plot(theta_deg,Dif,'LineWidth',2)
grid on
xlabel('Angle (deg)')
ylabel('Voltage')
title('Difference Pattern')

figure(4)
plot(theta_deg,Sum,theta_deg,Dif,'r','LineWidth',2)
grid on
xlabel('Angle (deg)')
ylabel('Voltage')
title('Sum & Difference Patterns')

figure(5)
plot(theta_deg,rho,theta_deg,rho2,'r','LineWidth',2)
hold on
plot(theta_deg,tan(15.7*theta_deg*pi/180),'k*')
grid on
xlabel('Angle (deg)')
ylabel('Ratio')
title('Monopulse Ratio')
angmax2 = floor(2*Null2Null_bw)/2;
axis([-angmax2 angmax2 -2 2])

figure(6)
mesh(theta_deg,theta_deg,abs(Sum2D))
colormap('jet')
xlabel('El')
ylabel('Az')
title('Sum Pattern')

%---------
azdata = theta_deg;
eldata = theta_deg;
sumdata = Sum2D;
difazdata = Dif;
difeldata = Dif;

