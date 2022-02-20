function T_EE = DHcalc(a,alpha,d,theta,offset)
% calculates end effector position and orientation using the 
% Denavit-Hartenberg Matrix method

% INPUTS
% a - translational offsets of the robot joint axes, 7x1 vector
% alphar - angular offsets of the robot joint axes, 7x1 vector
% d - lengths of the robot arms
% theta - input angles to servos, 7x1 vector
% o - offset of the robot joints, used for calibration only, rads

% OUTPUTS
% T_EE - DH configuration of the end effector transformed with respect to
% the base

for i = 1:length(a)
    T{i}(1,1) = cos(theta(i) + offset(i));
    T{i}(1,2) = -cos(alpha(i))*sin(theta(i) + offset(i));
    T{i}(1,3) = sin(alpha(i))*sin(theta(i) + offset(i));
    T{i}(1,4) = a(i)*cos(theta(i) + offset(i));
    
    T{i}(2,1) = sin(theta(i) + offset(i));
    T{i}(2,2) = cos(alpha(i))*cos(theta(i) + offset(i));
    T{i}(2,3) = -sin(alpha(i))*cos(theta(i) + offset(i));
    T{i}(2,4) = a(i)*sin(theta(i) + offset(i));

    T{i}(3,1) = 0;
    T{i}(3,2) = sin(alpha(i));
    T{i}(3,3) = cos(alpha(i));
    T{i}(3,4) = d(i);
    
    T{i}(4,1) = 0;
    T{i}(4,2) = 0;
    T{i}(4,3) = 0;
    T{i}(4,4) = 1;
end
T_EE = T{1}*T{2}*T{3}*T{4}*T{5}*T{6}*T{7};
end