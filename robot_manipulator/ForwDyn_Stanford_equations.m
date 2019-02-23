% Test program to get inverse dynamic resuls for Stanford
clc
clear all


% alpha=[0;-90;90];
% alpha=deg2rad(alpha);
a=zeros(3,1);

TypeOfJoints='rrp';
% m=[4;2;2];
I=zeros(3,3);
Pc=[0;0;0];


syms t theta1 theta2 dtheta1 dtheta2 ddtheta1 ddtheta2 f df ddf m1 m2 m3 h r alpha2 alpha3
m=[m1,m2,m3];
alpha=[0;alpha2;alpha3];
% theta1='0.785398163397448*(1+6*exp(-t/0.6)-8*exp(-t/0.8))';
% theta2='0.785398163397448*(1+6*exp(-t/0.6)-8*exp(-t/0.8))';
% f='0.4*(1+6*exp(-t/0.6)-6*exp(-t/0.8))';
% 
% dtheta1='0.785398163397448*(-6*exp(-t/0.6)/0.6+8*exp(-t/0.8)/0.8)';
% dtheta2='0.785398163397448*(-6*exp(-t/0.6)/0.6+8*exp(-t/0.8)/0.8)';
% df='0.4*(-6*exp(-t/0.6)/0.6+6*exp(-t/0.8)/0.8)';
% 
% ddtheta1='0.785398163397448*(6*exp(-t/0.6)/0.36-8*exp(-t/0.8)/0.64)';
% ddtheta2='0.785398163397448*(6*exp(-t/0.6)/0.36-8*exp(-t/0.8)/0.64)';
% ddf='0.4*(6*exp(-t/0.6)/0.36-6*exp(-t/0.8)/0.64)';

% time=0:0.01:3;
% for i=1:length(time)

%     theta=[subs(theta1,t,time(i));subs(theta2,t,time(i));0];
%     d=[0.4;0.1;subs(f,t,time(i))];
    theta=[theta1;theta2;0];
    d=[h;r;f];

    JointVelocities=[dtheta1;dtheta2;df];

    JointAcceleration=[ddtheta1;ddtheta2;ddf];

    tou=InvDyn_radianversion(JointVelocities,JointAcceleration,alpha,a,d,theta,TypeOfJoints,m,I,Pc);
    tou=subs(tou,alpha2,-1.570796326794897);
    tou=subs(tou,alpha3,1.570796326794897);
    tou=simple(tou)
%     T(i,:)=[time(i) tou(i,:)];
% end

% plot(time,tou(:,1),time,tou(:,2),time,tou(:,3))
M=[f^2*m3+m2*r^2-cos(theta2)^2*f^2*m3+m3*r^2,-r*cos(theta2)*m3*f,-r*sin(theta2)*m3
    -f*m3*cos(theta2)*r, f^2*m3,0
    m3*(-sin(theta2)*r),0,m3];
V=[2*sin(theta2)*f^2*m3*cos(theta2)*dtheta1*dtheta2+r*sin(theta2)*m3*dtheta2^2*f-2*r*cos(theta2)*m3*dtheta2*df+2*f*m3*dtheta1*df-2*cos(theta2)^2*f*m3*dtheta1*df
    -1/5*f*m3*(5*cos(theta2)*dtheta1^2*sin(theta2)*f-10*dtheta2*df)
1/5*m3*(-5*dtheta1^2*f-5*dtheta2^2*f+5*dtheta1^2*f*cos(theta2)^2)];

G=[0
   -1/5*f*m3*(49*sin(theta2)) 
1/5*m3*(49*cos(theta2))];

dt=0.01;
time=0:dt:3;
for i=1:length(time)
 






% tou1=2*sin(theta2)*f^2*m3*cos(theta2)*dtheta1*dtheta2+r*sin(theta2)*m3*dtheta2^2*f-r*sin(theta2)*m3*ddf-2*r*cos(theta2)*m3*dtheta2*df+f^2*m3*ddtheta1+2*f*m3*dtheta1*df-r*cos(theta2)*m3*ddtheta2*f+m2*ddtheta1*r^2-cos(theta2)^2*f^2*m3*ddtheta1-2*cos(theta2)^2*f*m3*dtheta1*df+m3*ddtheta1*r^2;
% tou2=-1/5*f*m3*(5*cos(theta2)*ddtheta1*r+49*sin(theta2)-5*ddtheta2*f+5*cos(theta2)*dtheta1^2*sin(theta2)*f-10*dtheta2*df);
% tou3=1/5*m3*(-5*sin(theta2)*ddtheta1*r-5*dtheta1^2*f+49*cos(theta2)-5*dtheta2^2*f+5*ddf+5*dtheta1^2*f*cos(theta2)^2);
