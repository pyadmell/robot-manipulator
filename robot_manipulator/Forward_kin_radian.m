
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
function T=Forward_kin_radian(alpha,a,d,theta,frame1,frame2)

T=eye(4);
for i=frame2:-1:frame1+1
    T=[cos((theta(i))) -sin((theta(i))) 0 a(i)
        sin((theta(i)))*cos((alpha(i))) cos((theta(i)))*cos((alpha(i))) ...
        -sin((alpha(i))) -sin((alpha(i)))*d(i)
        sin((theta(i)))*sin((alpha(i))) cos((theta(i)))*sin((alpha(i))) ...
        cos((alpha(i))) cos((alpha(i)))*d(i)
        0 0 0 1]*T;
end
