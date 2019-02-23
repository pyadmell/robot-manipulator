function V=StanfordVelocity(NumTheta1,NumTheta2,NumF,NumDtheta1,...
    NumDtheta2,NumDf)
% This function computes velocity vector of the Stanford manipulator given
% joint positions and velocities
%
%
%       Peyman Yadmellat
%       April 4, 2010
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
syms t theta1 theta2 dtheta1 dtheta2 ddtheta1 ddtheta2 f df ddf

r=0.1;
m=[4;2;2];
m1=m(1);
m2=m(2);
m3=m(3);
V=subs([2*sin(theta2)*f^2*m3*cos(theta2)*dtheta1*dtheta2+...
    r*sin(theta2)*m3*dtheta2^2*f-2*r*cos(theta2)*m3*dtheta2*df+...
    2*f*m3*dtheta1*df-2*cos(theta2)^2*f*m3*dtheta1*df
    -1/5*f*m3*(5*cos(theta2)*dtheta1^2*sin(theta2)*f-10*dtheta2*df)
1/5*m3*(-5*dtheta1^2*f-5*dtheta2^2*f+5*dtheta1^2*f*cos(theta2)^2)],...
{theta1 theta2 f dtheta1 dtheta2 df},...
{NumTheta1 NumTheta2 NumF NumDtheta1 NumDtheta2 NumDf});
