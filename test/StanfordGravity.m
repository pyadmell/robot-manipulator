function G=StanfordGravity(NumTheta1,NumTheta2,NumF)
% This function computes gravity vector of the Stanford manipulator given
% joint positions
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
G=subs([0
   -1/5*f*m3*(49*sin(theta2)) 
1/5*m3*(49*cos(theta2))],{theta1 theta2 f},{NumTheta1,NumTheta2,NumF});
