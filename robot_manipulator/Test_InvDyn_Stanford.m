% This program cab be used to get inverse dynamic and forward dynamic 
% resuls for the Stanford maipulator
%
%
%       Peyman Yadmellat
%       April 4, 2010
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


clc
clear all
close all
% Initializing D-H parameters for the Stanford manipulator
m=[4;2;2];
alpha=[0;-90;90];
a=zeros(3,1);
% determining joint type of the Stanford manipulator
TypeOfJoints='rrp';
% setting the Stanford manipulator configuration including mass, inertia
% tensor, and position of center of mass
m=[4;2;2];
I1=zeros(3,3);
I2=zeros(3,3);
I3=zeros(3,3);
I=[I1 I2 I3];

Pc1=[0;0;0];
Pc2=[0;0;0];
Pc3=[0;0;0];
Pc=[Pc1 Pc2 Pc3];

% defining joint positions, velocities and accelerations
syms t theta1 theta2 dtheta1 dtheta2 ddtheta1 ddtheta2 f df ddf
theta1='45*(1+6*exp(-t/0.6)-8*exp(-t/0.8))';
theta2='45*(1+6*exp(-t/0.6)-8*exp(-t/0.8))';
f='0.4*(1+6*exp(-t/0.6)-6*exp(-t/0.8))';

dtheta1='45*(-6*exp(-t/0.6)/0.6+8*exp(-t/0.8)/0.8)';
dtheta2='45*(-6*exp(-t/0.6)/0.6+8*exp(-t/0.8)/0.8)';
df='0.4*(-6*exp(-t/0.6)/0.6+6*exp(-t/0.8)/0.8)';

ddtheta1='45*(6*exp(-t/0.6)/0.36-8*exp(-t/0.8)/0.64)';
ddtheta2='45*(6*exp(-t/0.6)/0.36-8*exp(-t/0.8)/0.64)';
ddf='0.4*(6*exp(-t/0.6)/0.36-6*exp(-t/0.8)/0.64)';


dt=0.001;    %step time
time=0:dt:3; 

%initializing numerical joint positions and velocities for using in Euler
%integration method
NumDtheta(1,:)=zeros(1,3);
Numtheta(1,:)=[deg2rad(subs(theta1,t,time(1))),...
    deg2rad(subs(theta2,t,time(1))),subs(f,t,time(1))];

% "for" loop to compute desired joint positions, joint velocities, 
% joint accelerations, and torque along time using inverse dynamic. Also
% this loop is simultaneously used to compute joint position, velocities
% and acceleratin resulting from forward dynamic (simulation) by applying
% the torque profiles abtained form inverse dynamic
for i=1:length(time)
%inverse dynamic part:

    %Computing time-varying D-H parameters against time
    theta=[subs(theta1,t,time(i));subs(theta2,t,time(i));0];
    d=[0.4;0.1;subs(f,t,time(i))];
    
    %Computing desired joint velocities by substituting time
    JointVelocities=[subs(dtheta1,t,time(i));subs(dtheta2,t,time(i))...
        ;subs(df,t,time(i))];
    
    %Computing desired joint accelerations by substituting time
    JointAcceleration=[subs(ddtheta1,t,time(i));subs(ddtheta2,t,time(i))...
        ;subs(ddf,t,time(i))];
    
    %Computing torque by using inverse dynamic
    Tau(i,:)=InvDyn(JointVelocities,JointAcceleration,alpha,a,d,theta,...
        TypeOfJoints,m,I,Pc);
    
%Forward dynamic part:
    
    %seting variables needed to compute mass matrix, velocity vector 
    %and gravity vector for the Stanford manipulator 
    NumTheta1=Numtheta(i,1);
    NumTheta2=Numtheta(i,2);
    NumF=Numtheta(i,3);
    NumDtheta1=NumDtheta(i,1);
    NumDtheta2=NumDtheta(i,2);
    NumDf=NumDtheta(i,3);
    
    %computing mass matrix for the Stanford manipulator
    M=StanfordMass(NumTheta1,NumTheta2,NumF);
    detM(i)=det(M); %computing determinan of mass matrix
    
    %computing velocity vector for the Stanford manipulator
    V(i,:)=StanfordVelocity(NumTheta1,NumTheta2,NumF,NumDtheta1,...
        NumDtheta2,NumDf);
    %computing gravity vector for the Stanford manipulator
    G(i,:)=StanfordGravity(NumTheta1,NumTheta2,NumF);
    
    %using forward dynamic to compute joint acceleration by applying
    % the torque abtained from inverse dynamic stage
    NumDDtheta(i,:)=M\(Tau(i,:)-V(i,:)-G(i,:)).';
    
    %using Euler integration method to obtain joint velocities
    NumDtheta(i+1,:)=NumDtheta(i,:)+NumDDtheta(i,:)*dt;
    
    %using Euler integration method to obtain joint positions
    Numtheta(i+1,:)=Numtheta(i,:)+NumDtheta(i,:)*dt+...
        0.5*NumDDtheta(i,:)*dt^2;    
end

figure(1)
subplot(311)
plot(time, subs(theta1,t,time))
grid on
% Create ylabel
ylabel('\theta_1(t) (\circ)');

% Create title
title('Joint positions of the Stanford manipulator','FontWeight','bold');

subplot(312)
plot(time, subs(theta2,t,time))
grid on;

% Create ylabel
ylabel('\theta_2(t) (\circ)');

subplot(313)
plot(time, subs(f,t,time))
grid on
% Create xlabel
xlabel('Time (sec)');

% Create ylabel
ylabel('f(t) (m)');

figure(2)
subplot(311)
plot(time, subs(dtheta1,t,time))
grid on
% Create ylabel
ylabel('\theta\prime_1(t) (\circ/sec)');

% Create title
title('Joint velocities of the Stanford manipulator','FontWeight','bold');

subplot(312)
plot(time, subs(dtheta2,t,time))
grid on;

% Create ylabel
ylabel('\theta\prime_2(t) (\circ/sec)');

subplot(313)
plot(time, subs(df,t,time))
grid on
% Create xlabel
xlabel('Time (sec)');

% Create ylabel
ylabel('f\prime(t) (m/sec)');

figure(3)
subplot(311)
plot(time, subs(ddtheta1,t,time))
grid on
% Create ylabel
ylabel('\theta\prime\prime_1(t) (\circ/sec^2)');

% Create title
title('Joint accelerations of the Stanford manipulator',...
    'FontWeight','bold');

subplot(312)
plot(time, subs(ddtheta2,t,time))
grid on;

% Create ylabel
ylabel('\theta\prime\prime_2(t) (\circ/sec^2)');

subplot(313)
plot(time, subs(ddf,t,time))
grid on
% Create xlabel
xlabel('Time (sec)');

% Create ylabel
ylabel('f\prime\prime(t) (m/sec^2)');
figure(4)
plot(time,Tau(:,1),time,Tau(:,2),time,Tau(:,3))

% Create xlabel
xlabel('Time (sec)');

% Create ylabel
ylabel('\tau');

% Create title
title('Torque applied of each link of the Stanford manipulator',...
    'FontWeight','bold');
% Create legend
legend('1st link','2nd link','3rd link');

figure(5)
subplot(311)
plot(time, subs(theta1,time),'--')
hold on
plot(time,Numtheta(1:length(time),1)*180/pi)
% Create ylabel
ylabel('\theta_1(t) (\circ)');

% Create title
title('Simulation: Joint positions',...
    'FontWeight','bold');
% Create legend
legend('Desired positions','Simulation results');
subplot(312)
plot(time, subs(theta2,time),'--')
hold on
plot(time,Numtheta(1:length(time),2)*180/pi)
% Create ylabel
ylabel('\theta_2(t) (\circ)');
subplot(313)
plot(time, subs(f,time),'--')
hold on
plot(time,Numtheta(1:length(time),3))
% Create ylabel
ylabel('f(t) (m)');
% Create xlabel
xlabel('Time (sec)');

figure(6)
subplot(311)
plot(time, subs(dtheta1,time),'--')
hold on
plot(time,NumDtheta(1:length(time),1)*180/pi)
% Create ylabel
ylabel('\theta\prime_1(t) (\circ/sec)');
% Create title
title('Simulation: Joint velocities',...
    'FontWeight','bold');
% Create legend
legend('Desired velocities','Simulation results');

subplot(312)
plot(time, subs(dtheta2,time),'--')
hold on
plot(time,NumDtheta(1:length(time),2)*180/pi)
% Create ylabel
ylabel('\theta\prime_2(t) (\circ/sec)');
subplot(313)
plot(time, subs(df,time),'--')
hold on
plot(time,NumDtheta(1:length(time),1))
% Create ylabel
ylabel('f\prime(t) (m/sec)');
% Create xlabel
xlabel('Time (sec)');


figure(7)
subplot(311)
plot(time, subs(ddtheta1,time),'--')
hold on
plot(time,NumDDtheta(1:length(time),1)*180/pi)
% Create ylabel
ylabel('\theta\prime\prime_1(t) (\circ/sec^2)');
% Create title
title('Simulation: Joint accelerations',...
    'FontWeight','bold');
% Create legend
legend('Desired accelerations','Simulation results');
subplot(312)
plot(time, subs(ddtheta2,time),'--')
hold on
plot(time,NumDDtheta(1:length(time),2)*180/pi)
% Create ylabel
ylabel('\theta\prime\prime_2(t) (\circ/sec^2)');

subplot(313)
plot(time, subs(ddf,time),'--')
hold on
plot(time,NumDDtheta(1:length(time),3))
% Create ylabel
ylabel('f\prime\prime(t) (m/sec^2)');
% Create xlabel
xlabel('Time (sec)');

figure(8)
plot(time,detM(1:length(time)),'r')

