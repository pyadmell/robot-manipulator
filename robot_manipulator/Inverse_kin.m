function [f,theta1,theta2]=Inverse_kin(position)

%  [f,theta1,theta2]=Inverse_kin(position) returns the possible solutions
%  for f, theta1, and theta2 corresponding to the numerical values of 
%  position relevant to end-efector frame for Stanford manipulator
%           
%       position=[Px;Py;Pz]
%
%       Peyman Yadmellat
%       February 20, 2010
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

sym theta1 theta2 f 
sym Px Py Pz

Pxn=position(1); % numerical value of Px
Pyn=position(2); % numerical value of Py
Pzn=position(3); % numerical value of Pz

e1=subs('f*cos(theta1)*sin(theta2)-0.1*sin(theta1)-Px','Px',Pxn);
e2=subs('f*sin(theta1)*sin(theta2)+0.1*cos(theta1)-Py','Py',Pyn);
e3=subs('f*cos(theta2)+0.4-Pz','Pz',Pzn);

temp=solve(e1,e2,e3);
f=temp.f;
theta1=temp.theta1;
theta2=temp.theta2;
