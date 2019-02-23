function Tau=InvDyn(JointVelocities,JointAcceleration,alpha,a,d,theta,...
    TypeOfJoints,m,I,Pc)

%This function computes the torque applied by each link to its neighbor.
%------------------------------------------------------------------------- 
% Input arguments: 
%
%   JointVelocities is a column vector of velocity, 
%                   it should be in rad/sec. 
%
%   JointAcceleration is a column vector of joint acceleration
%                     it should be in rad/sec^2. 
%   
%   alpha,a,d,and theta are column vectors corresponding D-H parameters
%   
%   TypeOfJoints is a string containing the mechanical type of each joint.
%   Prismatic joint is represented by 'P', and 'R' stands for revolute
%   joint. Mechanical type of joints should be written in order. For 
%   example, 'RPRR' (equivalently ['R','P','R','R']) means that first
%   joint is revolute, second one is prismatic, and the other two are
%   revolute.   
%
%   m is a column vector of mass of each link
%
%   I is a matrix containing all inertia tensor matrices which they should
%    place in one matrix in order, for example I=[I1,,I2,...In]
%   
%   Pc should contain all position of center of mass vectors in order, for
%   example Pc=[Pc1,Pc2,...,Pcn]
%
%------------------------------------------------------------------------- 
% Output arguments:
% Tau	returns the torque applied by each link to its neighbor.
%-------------------------------------------------------------------------
%   
% Peyman Yadmellat
% April 04, 2010
% -------------------------------
% -------------------------------


TypeOfJoints=upper(TypeOfJoints); %Set to uppercase

numJoints=length(TypeOfJoints); %number of joints

G=[0;0;9.8];

%Initializing angular and linear velocities
v_new=[0;0;0];
w_new=zeros(3,1);

%Initializing angular and linear accelerations
dv_new=G;
dw_new=zeros(3,1);



% Computing the force and torque acting on each link, i.e. F adn N
% Outward recursions
for i=1:numJoints
    T=Forward_kin(alpha,a,d,theta,i-1,i);
    
    % Extracting the rotation matrix of joint i respect to its description
    % in frame {i-1}
    R=T(1:3,1:3).';
   
    % Extracting the translation matrix of joint i-1 respect to its
    % description in frame {i}
    P=T(1:3,4);

    % this set of formulas computes linear and angular velocities of each
    % joint
    v=v_new;
    w=w_new;
    v_new=R*(v+cross(w,P))+...
        isequal(TypeOfJoints(i),'P')*[0;0;JointVelocities(i)];
    w_new=R*w+isequal(TypeOfJoints(i),'R')*[0;0;JointVelocities(i)*pi/180];
    
    % this set of formulas computes linear and angular accelerations of 
    % each joint
    dv=dv_new;
    dw=dw_new;
    dw_new=R*dw+isequal(TypeOfJoints(i),'R')*(cross(R*w,...
        [0;0;JointVelocities(i)*pi/180])+...
        [0;0;JointAcceleration(i)*pi/180]);
    dv_new=R*(dv+cross(dw,P)+cross(w,cross(w,P)))+...
        isequal(TypeOfJoints(i),'P')*(cross(2*w_new,...
        [0;0;JointVelocities(i)])+[0;0;JointAcceleration(i)]);
    
    %Computing linear acceleration of the center of mass of each link   
    dvc=cross(w_new,Pc(1:numJoints,i))+...
        cross(w_new,w_new+Pc(1:numJoints,i))+dv_new;
    
    %Computing Force and torque acting on each link
    F(:,i)=m(i)*dvc;
    N(:,i)=I(1:numJoints,(i-1)*numJoints+1:i*numJoints)*dw_new+...
        cross(w_new,I(1:numJoints,(i-1)*numJoints+1:i*numJoints)*w_new);
end

% Inward recursions:

%initializing force and torque exerted on each link
f=zeros(3,1);
n=zeros(3,1);
R=eye(3);
P=zeros(3,1);
n=N(:,numJoints);
f=F(:,numJoints);
Tau(:,numJoints)=isequal(TypeOfJoints(i),'P')*f(3)+...
    isequal(TypeOfJoints(i),'R')*n(3);

% Computing the torque applied  by each link on its neighbor    
for i=numJoints-1:-1:1
    T=Forward_kin(alpha,a,d,theta,i,i+1);
    R=T(1:3,1:3);
    P=T(1:3,4);
    n=N(:,i)+R*n+cross(Pc(1:numJoints,i),F(:,i))+cross(P,R*f);
    f=R*f+F(:,i);
    Tau(:,i)=isequal(TypeOfJoints(i),'P')*f(3)+isequal(TypeOfJoints(i),'R')*n(3);
    
end

