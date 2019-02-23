function Tou=InvDyn_radianversion(JointVelocities,JointAcceleration,alpha,a,d,theta,TypeOfJoints,m,I,Pc)
% Robotic manipulators - Project part 4
%
% This function computes the Cartesian velocity of the end-effector, given
% a set of velocities and manipulator's configuration
%
% According to mechanical type of each joint, Two sets of equations 
% can be used to compute the Cartesian velocities of end-effector: 
%   
% In the case of Revolute Joints, link to link velocities can be calculated
% as given below:
% 
%        i+1      i+1  i      .     i+1^
%           w    =   R  w + theta      Z                        (R-1)
%            i+1    i    i   i+1    i+1 
%----------------------------------------------------
%      
%        i+1      i+1   i    i     i
%           v    =   R ( v +  w  x  P   )                       (R-2)
%            i+1    i     i    i     i+1
%
% And in case of Prismatic Joints, link to link velocities can be 
% calculated as given below: 
%
%        i+1      i+1  i  
%           w    =   R  w                                        (P-1)
%            i+1    i    i 
%----------------------------------------------------
%      
%        i+1      i+1   i    i     i       .   i+1^ 
%           v    =   R ( v +  w  x  P   )+ d      Z              (P-2)
%            i+1    i     i    i     i+1      i+1    i+1     
%       
%                                           N       N
% Also to rotate resulting velocities, i.e.  v  and  w, into base
%                                           N       N
% coordinates, following equation should be used:
% 
%           0     0  N
%            v  =  R  v                                          
%             N   N    N
%     ------------------------
%           0     0  N
%            w  =  R  w                                          
%             N   N    N
%
%------------------------------------------------------------------------- 
% Input arguments: 
%
%   JointVelocities is an either column or row vector of velocity. 
%   Default value is set as JointVelocities=[35;20;0.1].
%   
%   alpha,a,d,and theta are column vectors corresponding D-H parameters
%   Default values of D-H parameters are set as 
%   alpha=[0;-90;90], a=[0;0;0], d=[0.4;0.1;0.3536], and 
%   theta=[9.22;115.1;0].
%   
%   TypeOfJoints is a string containing the mechanical type of each joint.
%   Prismatic joint is represented by 'P', and 'R' stands for revolute
%   joint. Mechanical type of joints should be written in order. For 
%   example, 'RPRR' (equivalently ['R','P','R','R']) means that first
%   joint is revolute, second one is prismatic, and the other two are
%   revolute.   
%
% Note that the default values will be considered only in a case of missing
% argument(s). In other word, the default values will be considered
% automatically if those arguments are not entered.
%------------------------------------------------------------------------- 
% Output arguments:
% v	returns the column vector of linear velocities of end-effector 
%   or the last frame of manipulator in meter/second
% w	returns the column vector of angular velocity of end-effector
%   or the last frame of manipulator in degree/second
%-------------------------------------------------------------------------
%   
% Peyman Yadmellat
% March 29, 2010
% -------------------------------
% -------------------------------

% % default type of joints is set for Stanford manipulator
% defaultTypeOfJoints = 'RRP';
% 
% %default DH parameters considered as DH parameters of Stanford manipulator
% default_alpha=[0;-90;90];
% default_a=zeros(3,1);
% default_d=[0.4;0.1;0.3536];
% default_theta=[9.22;115.1;0];
% 
% %default joint velocities are set as given in project definition
% defaultJointVelocities=[35;20;0.1];
% 
% % Handle missing arguments and set each missed argument as it should be
% % for Standford manipulator
% if nargin<6,TypeOfJoints = defaultTypeOfJoints;
%     if nargin<5,theta=default_theta;
%         if nargin<4,d=default_d;
%             if nargin<3,a=default_a;
%                 if nargin<2,alpha=default_alpha;
%                     if nargin<1,JointVelocities=defaultJointVelocities;
%                     end,end,end,end,end,end
% 
TypeOfJoints=upper(TypeOfJoints); %Set to uppercase


numJoints=length(TypeOfJoints); %number of joints

G=9.8;

v_new=[0;0;0];
w_new=zeros(3,1);

dv_new=[0;0;G];
dw_new=zeros(3,1);



% Computing linear and angular velocities of end effector
% respect to end effector frame
for i=1:numJoints
    T=Forward_kin_radian(alpha,a,d,theta,i-1,i);
    
    % Extracting the rotation matrix of joint i respect to its description
    % in frame {i-1}
    R=T(1:3,1:3).';
   
    % Extracting the translation matrix of joint i-1 respect to its
    % description in frame {i}
    P=T(1:3,4);

    % this set of formulas computes linear and angular velocities of each
    % joint regarding the equation (3) in report in which these formulas
    % switch between two possible cases according to the mechanical type
    % of each joint. The switching parameter is directly assigned by
    % comparing the mechanical type of each joint to possible Prismatic and
    % Revolute cases by using "isequal" command
    v=v_new;
    w=w_new;
    v_new=R*(v+cross(w,P))+...
        isequal(TypeOfJoints(i),'P')*[0;0;JointVelocities(i)];
    w_new=R*w+isequal(TypeOfJoints(i),'R')*[0;0;JointVelocities(i)];
    
    dv=dv_new;
    dw=dw_new;
    dw_new=R*dw+isequal(TypeOfJoints(i),'R')*(cross(R*w,[0;0;JointVelocities(i)])+[0;0;JointAcceleration(i)]);
    dv_new=R*(dv+cross(dw,P)+cross(w,cross(w,P)))+...
        isequal(TypeOfJoints(i),'P')*(cross(2*w_new,[0;0;JointVelocities(i)])+[0;0;JointAcceleration(i)]);
    
      
    dvc=cross(w_new,Pc)+cross(w_new,w_new+Pc)+dv_new;
    
    F(:,i)=m(i)*dvc;
    N(:,i)=I*dw_new+cross(w_new,I*w_new);
end

f=zeros(3,1);
n=zeros(3,1);
R=eye(3);
P=zeros(3,1);
n=N(:,numJoints);
f=F(:,numJoints);
Tou(:,numJoints)=isequal(TypeOfJoints(i),'P')*f(3)+isequal(TypeOfJoints(i),'R')*n(3);
    
for i=numJoints-1:-1:1
    T=Forward_kin_radian(alpha,a,d,theta,i,i+1);
    R=T(1:3,1:3);
    P=T(1:3,4);
    n=N(:,i)+R*n+cross(Pc,F(:,i))+cross(P,R*f);
    f=R*f+F(:,i);
    Tou(:,i)=isequal(TypeOfJoints(i),'P')*f(3)+isequal(TypeOfJoints(i),'R')*n(3);
    
end


% T=Forward_kin(alpha,a,d,theta,0,n);
% 
% %Computing the rotation matrix to rotate velocities of the end-effector
% %respect to base coordinates
% R=T(1:3,1:3);
% 
% % Returns linear velocity respect to base frame in m/s
% v=R*v;
% 
% % Returns angular velocity respect to base frame in degree/sec
% w=180/pi*R*w;
