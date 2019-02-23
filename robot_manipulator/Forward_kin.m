
%  T=Forward_kin(alpha,a,d,theta,frame1,frame2) gives the transformation T 
%  which transforms vector defined in {frame2} to their description in
%  {frame1} corresponding to D-H parameters {alpha,a,d & theta}
%           
%       alpha is the angle between Z_{i} to Z_{i+1} measured about X_{i} in
%               degree 
%       
%       a     is the distance from Z_{i} to Z_{i+1} measured along X_{i}
%       
%       d     is the distance from X_{i-1} to X_{i} measured along Z_{i} 
%       
%       theta is the angle between X_{i-1} to X_{i} measured about Z_{i} in
%             degree
%
%
%           frame1      frame1    frame1+1          frame2-1
%                 T =         T x         T x ... x         T
%           frame2    frame1+1    frame1+2            frame2
%
%
%
%       Peyman Yadmellat
%       February 20, 2010
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function T=Forward_kin(alpha,a,d,theta,frame1,frame2)

T=eye(4);
for i=frame2:-1:frame1+1
    T=[cos(deg2rad(theta(i))) -sin(deg2rad(theta(i))) 0 a(i)
        sin(deg2rad(theta(i)))*cos(deg2rad(alpha(i))) cos(deg2rad(theta(i)))*cos(deg2rad(alpha(i))) ...
        -sin(deg2rad(alpha(i))) -sin(deg2rad(alpha(i)))*d(i)
        sin(deg2rad(theta(i)))*sin(deg2rad(alpha(i))) cos(deg2rad(theta(i)))*sin(deg2rad(alpha(i))) ...
        cos(deg2rad(alpha(i))) cos(deg2rad(alpha(i)))*d(i)
        0 0 0 1]*T;
end
