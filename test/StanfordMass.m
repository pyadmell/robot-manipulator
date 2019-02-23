function M=StanfordMass(NumTheta1,NumTheta2,NumF)
% This function computes mass matrix of the Stanford manipulator given
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
M=subs([f^2*m(3)+m(2)*r^2-cos(theta2)^2*f^2*m(3)+m(3)*r^2,...
    -r*cos(theta2)*m(3)*f,-r*sin(theta2)*m(3)
    -f*m(3)*cos(theta2)*r, f^2*m(3),0
    m(3)*(-sin(theta2)*r),0,m(3)],{theta1 theta2 f},...
    {NumTheta1 NumTheta2 NumF});
