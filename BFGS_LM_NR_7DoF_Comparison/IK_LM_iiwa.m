function [q,k,qcheck,time] = IK_LM_iiwa(a,alpha,d,q0,o,PRdEE,phi_max,lambda_cur,v,WEE)
% calculates servo input angles to acheive a desired end effector position
% and orientation using the Levenberg-Marquardt (damped least squares)
% method

% INPUTS
% a = axis offset of robot joints
% alpha = angle offset of robot joints
% d = length of robot arms
% q0 = servo input angles initial guess
% o = angular offset of servos (typically zero)
% PRdEE = end effector desired cartesian position and orientation which has
% the form
% [nxdEE; nydEE; nzdEE; sxdEE; sydEE; szdEE; axdEE; aydEE; azdEE; PxdEE;
% PydEE; PzdEE]
% phi_max = maximum acceptable value for phi (error function)
% lambda = intial guess for lambda, damping factor
% v = scaling factor for lambda calculation, determined experimentally
% WEE = weighted matrix for end effector cartesian position and
% orientation, input as a 1x12 vector

% OUTPUTS
% q = 1x7 vector of servo angles to achieve desired end effector position
% and orientation, PRdEE

q_cur = q0';
qcheck(:,1) = q_cur;
k = 2;

%% END EFFECTOR
tic
phi_prevEE = phicalcEE(PRdEE,WEE,d,q_cur,a,alpha,o);
while phi_prevEE > phi_max
    % update lambda
    lambda_cur = lambdacalcEE(q_cur,o,d,a,alpha,PRdEE,lambda_cur,phi_prevEE,WEE,v);
    
    % update delta
    delta_cur = deltacalcEE(q_cur,o,d,a,alpha,PRdEE,lambda_cur);
    
    % update q
    q_cur = q_cur + delta_cur;
    
    % update phi
    phi_prevEE = phicalcEE(PRdEE,WEE,d,q_cur,a,alpha,o);
    qcheck(:,k) = wrapToPi(q_cur');
    k = k+1;
end
k = 1:k-1;
q = wrapToPi(q_cur);
time = toc;

% calculate phi for the end effector
function phiEE = phicalcEE(d_out,W_EE_vec,d,q,a,alpha,o)

T_temp = DH(a,alpha,d,q,o);
T_cur = reshape(T_temp(1:3,:),[12,1]);

% create residuals vector
e = d_out' - T_cur;
% create weighted matrix
W_EE = diag(W_EE_vec);
% calculate phi
phiEE = 0.5*e'*W_EE*e;
end

function J = Jcalc(q,d)
% calculate the Jacobian
d7 = d(7);
J(1,1) = cos(q(7))*(sin(q(6))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2))) + cos(q(6))*(cos(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))))) - sin(q(7))*(sin(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))));
J(2,1) = cos(q(7))*(sin(q(6))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2))) + cos(q(6))*(cos(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3))))) - sin(q(7))*(sin(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3))));
J(3,1) = 0;
J(4,1) = - sin(q(7))*(sin(q(6))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2))) + cos(q(6))*(cos(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))))) - cos(q(7))*(sin(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))));
J(5,1) = - sin(q(7))*(sin(q(6))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2))) + cos(q(6))*(cos(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3))))) - cos(q(7))*(sin(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3))));
J(6,1) = 0;
J(7,1) = cos(q(6))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2))) - sin(q(6))*(cos(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))));
J(8,1) = cos(q(6))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2))) - sin(q(6))*(cos(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3))));
J(9,1) = 0;
J(10,1) = (1199*cos(q(6))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2))))/10000 - (21*sin(q(1))*sin(q(2)))/50 - (1199*sin(q(6))*(cos(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3)))))/10000 + (2*sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))))/5 - (2*cos(q(4))*sin(q(1))*sin(q(2)))/5;
J(11,1) = (21*cos(q(1))*sin(q(2)))/50 + (1199*cos(q(6))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2))))/10000 - (1199*sin(q(6))*(cos(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)))))/10000 + (2*sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))))/5 + (2*cos(q(1))*cos(q(4))*sin(q(2)))/5;
J(12,1) = 0;

J(1,2) = cos(q(7))*(sin(q(6))*(cos(q(1))*cos(q(2))*cos(q(4)) + cos(q(1))*cos(q(3))*sin(q(2))*sin(q(4))) - cos(q(6))*(cos(q(5))*(cos(q(1))*cos(q(2))*sin(q(4)) - cos(q(1))*cos(q(3))*cos(q(4))*sin(q(2))) + cos(q(1))*sin(q(2))*sin(q(3))*sin(q(5)))) + sin(q(7))*(sin(q(5))*(cos(q(1))*cos(q(2))*sin(q(4)) - cos(q(1))*cos(q(3))*cos(q(4))*sin(q(2))) - cos(q(1))*cos(q(5))*sin(q(2))*sin(q(3)));
J(2,2) = cos(q(7))*(sin(q(6))*(cos(q(2))*cos(q(4))*sin(q(1)) + cos(q(3))*sin(q(1))*sin(q(2))*sin(q(4))) - cos(q(6))*(cos(q(5))*(cos(q(2))*sin(q(1))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(1))*sin(q(2))) + sin(q(1))*sin(q(2))*sin(q(3))*sin(q(5)))) + sin(q(7))*(sin(q(5))*(cos(q(2))*sin(q(1))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(1))*sin(q(2))) - cos(q(5))*sin(q(1))*sin(q(2))*sin(q(3)));
J(3,2) = cos(q(7))*(cos(q(6))*(cos(q(5))*(sin(q(2))*sin(q(4)) + cos(q(2))*cos(q(3))*cos(q(4))) - cos(q(2))*sin(q(3))*sin(q(5))) - sin(q(6))*(cos(q(4))*sin(q(2)) - cos(q(2))*cos(q(3))*sin(q(4)))) - sin(q(7))*(sin(q(5))*(sin(q(2))*sin(q(4)) + cos(q(2))*cos(q(3))*cos(q(4))) + cos(q(2))*cos(q(5))*sin(q(3)));
J(4,2) = cos(q(7))*(sin(q(5))*(cos(q(1))*cos(q(2))*sin(q(4)) - cos(q(1))*cos(q(3))*cos(q(4))*sin(q(2))) - cos(q(1))*cos(q(5))*sin(q(2))*sin(q(3))) - sin(q(7))*(sin(q(6))*(cos(q(1))*cos(q(2))*cos(q(4)) + cos(q(1))*cos(q(3))*sin(q(2))*sin(q(4))) - cos(q(6))*(cos(q(5))*(cos(q(1))*cos(q(2))*sin(q(4)) - cos(q(1))*cos(q(3))*cos(q(4))*sin(q(2))) + cos(q(1))*sin(q(2))*sin(q(3))*sin(q(5))));
J(5,2) = cos(q(7))*(sin(q(5))*(cos(q(2))*sin(q(1))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(1))*sin(q(2))) - cos(q(5))*sin(q(1))*sin(q(2))*sin(q(3))) - sin(q(7))*(sin(q(6))*(cos(q(2))*cos(q(4))*sin(q(1)) + cos(q(3))*sin(q(1))*sin(q(2))*sin(q(4))) - cos(q(6))*(cos(q(5))*(cos(q(2))*sin(q(1))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(1))*sin(q(2))) + sin(q(1))*sin(q(2))*sin(q(3))*sin(q(5))));
J(6,2) = - cos(q(7))*(sin(q(5))*(sin(q(2))*sin(q(4)) + cos(q(2))*cos(q(3))*cos(q(4))) + cos(q(2))*cos(q(5))*sin(q(3))) - sin(q(7))*(cos(q(6))*(cos(q(5))*(sin(q(2))*sin(q(4)) + cos(q(2))*cos(q(3))*cos(q(4))) - cos(q(2))*sin(q(3))*sin(q(5))) - sin(q(6))*(cos(q(4))*sin(q(2)) - cos(q(2))*cos(q(3))*sin(q(4))));
J(7,2) = cos(q(6))*(cos(q(1))*cos(q(2))*cos(q(4)) + cos(q(1))*cos(q(3))*sin(q(2))*sin(q(4))) + sin(q(6))*(cos(q(5))*(cos(q(1))*cos(q(2))*sin(q(4)) - cos(q(1))*cos(q(3))*cos(q(4))*sin(q(2))) + cos(q(1))*sin(q(2))*sin(q(3))*sin(q(5)));
J(8,2) = cos(q(6))*(cos(q(2))*cos(q(4))*sin(q(1)) + cos(q(3))*sin(q(1))*sin(q(2))*sin(q(4))) + sin(q(6))*(cos(q(5))*(cos(q(2))*sin(q(1))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(1))*sin(q(2))) + sin(q(1))*sin(q(2))*sin(q(3))*sin(q(5)));
J(9,2) = - sin(q(6))*(cos(q(5))*(sin(q(2))*sin(q(4)) + cos(q(2))*cos(q(3))*cos(q(4))) - cos(q(2))*sin(q(3))*sin(q(5))) - cos(q(6))*(cos(q(4))*sin(q(2)) - cos(q(2))*cos(q(3))*sin(q(4)));
J(10,2) = (21*cos(q(1))*cos(q(2)))/50 + (1199*cos(q(6))*(cos(q(1))*cos(q(2))*cos(q(4)) + cos(q(1))*cos(q(3))*sin(q(2))*sin(q(4))))/10000 + (1199*sin(q(6))*(cos(q(5))*(cos(q(1))*cos(q(2))*sin(q(4)) - cos(q(1))*cos(q(3))*cos(q(4))*sin(q(2))) + cos(q(1))*sin(q(2))*sin(q(3))*sin(q(5))))/10000 + (2*cos(q(1))*cos(q(2))*cos(q(4)))/5 + (2*cos(q(1))*cos(q(3))*sin(q(2))*sin(q(4)))/5;
J(11,2) = (21*cos(q(2))*sin(q(1)))/50 + (1199*cos(q(6))*(cos(q(2))*cos(q(4))*sin(q(1)) + cos(q(3))*sin(q(1))*sin(q(2))*sin(q(4))))/10000 + (1199*sin(q(6))*(cos(q(5))*(cos(q(2))*sin(q(1))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(1))*sin(q(2))) + sin(q(1))*sin(q(2))*sin(q(3))*sin(q(5))))/10000 + (2*cos(q(2))*cos(q(4))*sin(q(1)))/5 + (2*cos(q(3))*sin(q(1))*sin(q(2))*sin(q(4)))/5;
J(12,2) = (2*cos(q(2))*cos(q(3))*sin(q(4)))/5 - (2*cos(q(4))*sin(q(2)))/5 - (1199*sin(q(6))*(cos(q(5))*(sin(q(2))*sin(q(4)) + cos(q(2))*cos(q(3))*cos(q(4))) - cos(q(2))*sin(q(3))*sin(q(5))))/10000 - (1199*cos(q(6))*(cos(q(4))*sin(q(2)) - cos(q(2))*cos(q(3))*sin(q(4))))/10000 - (21*sin(q(2)))/50;

J(1,3) = - sin(q(7))*(cos(q(5))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(4))*sin(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)))) - cos(q(7))*(cos(q(6))*(sin(q(5))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(4))*cos(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)))) - sin(q(4))*sin(q(6))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3))));
J(2,3) = sin(q(7))*(cos(q(5))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + cos(q(4))*sin(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3)))) + cos(q(7))*(cos(q(6))*(sin(q(5))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*cos(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3)))) - sin(q(4))*sin(q(6))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))));
J(3,3) = - sin(q(7))*(cos(q(3))*cos(q(5))*sin(q(2)) - cos(q(4))*sin(q(2))*sin(q(3))*sin(q(5))) - cos(q(7))*(cos(q(6))*(cos(q(3))*sin(q(2))*sin(q(5)) + cos(q(4))*cos(q(5))*sin(q(2))*sin(q(3))) + sin(q(2))*sin(q(3))*sin(q(4))*sin(q(6)));
J(4,3) = sin(q(7))*(cos(q(6))*(sin(q(5))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(4))*cos(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)))) - sin(q(4))*sin(q(6))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)))) - cos(q(7))*(cos(q(5))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(4))*sin(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3))));
J(5,3) = cos(q(7))*(cos(q(5))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + cos(q(4))*sin(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3)))) - sin(q(7))*(cos(q(6))*(sin(q(5))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*cos(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3)))) - sin(q(4))*sin(q(6))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))));
J(6,3) = sin(q(7))*(cos(q(6))*(cos(q(3))*sin(q(2))*sin(q(5)) + cos(q(4))*cos(q(5))*sin(q(2))*sin(q(3))) + sin(q(2))*sin(q(3))*sin(q(4))*sin(q(6))) - cos(q(7))*(cos(q(3))*cos(q(5))*sin(q(2)) - cos(q(4))*sin(q(2))*sin(q(3))*sin(q(5)));
J(7,3) = sin(q(6))*(sin(q(5))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(4))*cos(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)))) + cos(q(6))*sin(q(4))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)));
J(8,3) = - sin(q(6))*(sin(q(5))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*cos(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3)))) - cos(q(6))*sin(q(4))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3)));
J(9,3) = sin(q(6))*(cos(q(3))*sin(q(2))*sin(q(5)) + cos(q(4))*cos(q(5))*sin(q(2))*sin(q(3))) - cos(q(6))*sin(q(2))*sin(q(3))*sin(q(4));
J(10,3) = (2*sin(q(4))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3))))/5 + (1199*sin(q(6))*(sin(q(5))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(4))*cos(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)))))/10000 + (1199*cos(q(6))*sin(q(4))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3))))/10000;
J(11,3) = - (2*sin(q(4))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))))/5 - (1199*sin(q(6))*(sin(q(5))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*cos(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3)))))/10000 - (1199*cos(q(6))*sin(q(4))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))))/10000;
J(12,3) = (1199*sin(q(6))*(cos(q(3))*sin(q(2))*sin(q(5)) + cos(q(4))*cos(q(5))*sin(q(2))*sin(q(3))))/10000 - (2*sin(q(2))*sin(q(3))*sin(q(4)))/5 - (1199*cos(q(6))*sin(q(2))*sin(q(3))*sin(q(4)))/10000;

J(1,4) = cos(q(7))*(sin(q(6))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*cos(q(6))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2)))) + sin(q(5))*sin(q(7))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2)));
J(2,4) = - cos(q(7))*(sin(q(6))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*cos(q(6))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2)))) - sin(q(5))*sin(q(7))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2)));
J(3,4) = sin(q(5))*sin(q(7))*(cos(q(2))*cos(q(4)) + cos(q(3))*sin(q(2))*sin(q(4))) - cos(q(7))*(sin(q(6))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) + cos(q(5))*cos(q(6))*(cos(q(2))*cos(q(4)) + cos(q(3))*sin(q(2))*sin(q(4))));
J(4,4) = cos(q(7))*sin(q(5))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2))) - sin(q(7))*(sin(q(6))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*cos(q(6))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2))));
J(5,4) = sin(q(7))*(sin(q(6))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*cos(q(6))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2)))) - cos(q(7))*sin(q(5))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2)));
J(6,4) = sin(q(7))*(sin(q(6))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) + cos(q(5))*cos(q(6))*(cos(q(2))*cos(q(4)) + cos(q(3))*sin(q(2))*sin(q(4)))) + cos(q(7))*sin(q(5))*(cos(q(2))*cos(q(4)) + cos(q(3))*sin(q(2))*sin(q(4)));
J(7,4) = cos(q(6))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) + cos(q(5))*sin(q(6))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2)));
J(8,4) = - cos(q(6))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*sin(q(6))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2)));
J(9,4) = cos(q(5))*sin(q(6))*(cos(q(2))*cos(q(4)) + cos(q(3))*sin(q(2))*sin(q(4))) - cos(q(6))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2)));
J(10,4) = (1199*cos(q(6))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))))/10000 + (2*cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))))/5 + (1199*cos(q(5))*sin(q(6))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2))))/10000 - (2*cos(q(1))*sin(q(2))*sin(q(4)))/5;
J(11,4) = - (1199*cos(q(6))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))))/10000 - (2*cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))))/5 - (2*sin(q(1))*sin(q(2))*sin(q(4)))/5 - (1199*cos(q(5))*sin(q(6))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2))))/10000;
J(12,4) = (1199*cos(q(5))*sin(q(6))*(cos(q(2))*cos(q(4)) + cos(q(3))*sin(q(2))*sin(q(4))))/10000 - (1199*cos(q(6))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))))/10000 - (2*cos(q(2))*sin(q(4)))/5 + (2*cos(q(3))*cos(q(4))*sin(q(2)))/5;

J(1,5) = - sin(q(7))*(cos(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)))) - cos(q(6))*cos(q(7))*(sin(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3))));
J(2,5) = sin(q(7))*(cos(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3)))) + cos(q(6))*cos(q(7))*(sin(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))));
J(3,5) = sin(q(7))*(cos(q(5))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) + sin(q(2))*sin(q(3))*sin(q(5))) + cos(q(6))*cos(q(7))*(sin(q(5))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) - cos(q(5))*sin(q(2))*sin(q(3)));
J(4,5) = cos(q(6))*sin(q(7))*(sin(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)))) - cos(q(7))*(cos(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3))));
J(5,5) = cos(q(7))*(cos(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3)))) - cos(q(6))*sin(q(7))*(sin(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))));
J(6,5) = cos(q(7))*(cos(q(5))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) + sin(q(2))*sin(q(3))*sin(q(5))) - cos(q(6))*sin(q(7))*(sin(q(5))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) - cos(q(5))*sin(q(2))*sin(q(3)));
J(7,5) = sin(q(6))*(sin(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3))));
J(8,5) = -sin(q(6))*(sin(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))));
J(9,5) = -sin(q(6))*(sin(q(5))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) - cos(q(5))*sin(q(2))*sin(q(3)));
J(10,5) = (1199*sin(q(6))*(sin(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)))))/10000;
J(11,5) = -(1199*sin(q(6))*(sin(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3)))))/10000;
J(12,5) = -(1199*sin(q(6))*(sin(q(5))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) - cos(q(5))*sin(q(2))*sin(q(3))))/10000;

J(1,6) = cos(q(7))*(cos(q(6))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2))) - sin(q(6))*(cos(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)))));
J(2,6) = -cos(q(7))*(cos(q(6))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2))) - sin(q(6))*(cos(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3)))));
J(3,6) = cos(q(7))*(sin(q(6))*(cos(q(5))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) + sin(q(2))*sin(q(3))*sin(q(5))) + cos(q(6))*(cos(q(2))*cos(q(4)) + cos(q(3))*sin(q(2))*sin(q(4))));
J(4,6) = -sin(q(7))*(cos(q(6))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2))) - sin(q(6))*(cos(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)))));
J(5,6) = sin(q(7))*(cos(q(6))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2))) - sin(q(6))*(cos(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3)))));
J(6,6) = -sin(q(7))*(sin(q(6))*(cos(q(5))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) + sin(q(2))*sin(q(3))*sin(q(5))) + cos(q(6))*(cos(q(2))*cos(q(4)) + cos(q(3))*sin(q(2))*sin(q(4))));
J(7,6) = - sin(q(6))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2))) - cos(q(6))*(cos(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3))));
J(8,6) = sin(q(6))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2))) + cos(q(6))*(cos(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))));
J(9,6) = cos(q(6))*(cos(q(5))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) + sin(q(2))*sin(q(3))*sin(q(5))) - sin(q(6))*(cos(q(2))*cos(q(4)) + cos(q(3))*sin(q(2))*sin(q(4)));
J(10,6) = - (1199*sin(q(6))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2))))/10000 - (1199*cos(q(6))*(cos(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)))))/10000;
J(11,6) = (1199*sin(q(6))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2))))/10000 + (1199*cos(q(6))*(cos(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3)))))/10000;
J(12,6) =(1199*cos(q(6))*(cos(q(5))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) + sin(q(2))*sin(q(3))*sin(q(5))))/10000 - (1199*sin(q(6))*(cos(q(2))*cos(q(4)) + cos(q(3))*sin(q(2))*sin(q(4))))/10000;

J(1,7) = - sin(q(7))*(sin(q(6))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2))) + cos(q(6))*(cos(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3))))) - cos(q(7))*(sin(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3))));
J(2,7) = sin(q(7))*(sin(q(6))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2))) + cos(q(6))*(cos(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))))) + cos(q(7))*(sin(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))));
J(3,7) = cos(q(7))*(sin(q(5))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) - cos(q(5))*sin(q(2))*sin(q(3))) + sin(q(7))*(cos(q(6))*(cos(q(5))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) + sin(q(2))*sin(q(3))*sin(q(5))) - sin(q(6))*(cos(q(2))*cos(q(4)) + cos(q(3))*sin(q(2))*sin(q(4))));
J(4,7) = sin(q(7))*(sin(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)))) - cos(q(7))*(sin(q(6))*(sin(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) + cos(q(1))*cos(q(4))*sin(q(2))) + cos(q(6))*(cos(q(5))*(cos(q(4))*(sin(q(1))*sin(q(3)) - cos(q(1))*cos(q(2))*cos(q(3))) - cos(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(3))*sin(q(1)) + cos(q(1))*cos(q(2))*sin(q(3)))));
J(5,7) = cos(q(7))*(sin(q(6))*(sin(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) - cos(q(4))*sin(q(1))*sin(q(2))) + cos(q(6))*(cos(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) + sin(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))))) - sin(q(7))*(sin(q(5))*(cos(q(4))*(cos(q(1))*sin(q(3)) + cos(q(2))*cos(q(3))*sin(q(1))) + sin(q(1))*sin(q(2))*sin(q(4))) - cos(q(5))*(cos(q(1))*cos(q(3)) - cos(q(2))*sin(q(1))*sin(q(3))));
J(6,7) = cos(q(7))*(cos(q(6))*(cos(q(5))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) + sin(q(2))*sin(q(3))*sin(q(5))) - sin(q(6))*(cos(q(2))*cos(q(4)) + cos(q(3))*sin(q(2))*sin(q(4)))) - sin(q(7))*(sin(q(5))*(cos(q(2))*sin(q(4)) - cos(q(3))*cos(q(4))*sin(q(2))) - cos(q(5))*sin(q(2))*sin(q(3)));
J(7,7) = 0;
J(8,7) = 0;
J(9,7) = 0;
J(10,7) = 0;
J(11,7) = 0;
J(12,7) = 0;
end

% calculate delta for the end effector
function deltaEE = deltacalcEE(q,o,d,a,alpha,d_out,lambda)
J = Jcalc(q,d);
T_temp = DH(a,alpha,d,q,o);
T_cur = reshape(T_temp(1:3,:),[12,1]);

% create residuals vector
e = d_out' - T_cur;

% calculate delta
lhs = transpose(J)*J + lambda^2*eye(7);
rhs = transpose(J)*e;
deltaEE = lhs\rhs;
end

% calculate lambda
function lambdaEE = lambdacalcEE(q_cur,o,d,a,alpha,d_out,lambda_r,phi_prev,W_EE_vec,v)
% compute new lambda given prior phi value and current q value

% try case 1 lambda = lambda_{r-1} / v
lambda_temp = lambda_r/v;

% update delta
delta_temp = deltacalcEE(q_cur,o,d,a,alpha,d_out,lambda_temp);

% update q
q_temp = q_cur + delta_temp;

% update phi
phi_temp = phicalcEE(d_out,W_EE_vec,d,q_temp,a,alpha,o);

if phi_temp < phi_prev
    lambdaEE = lambda_temp;
    %disp('case 1 used')
else
    % if case 1 fails, try case 2
    lambda_temp = lambda_r;    
    delta_temp = deltacalcEE(q_cur,o,d,a,alpha,d_out,lambda_temp);
    q_temp = q_cur + delta_temp;
    phi_temp = phicalcEE(d_out,W_EE_vec,d,q_temp,a,alpha,o);
    if phi_temp < phi_prev
        lambdaEE = lambda_temp;
        %disp('case 2 used')
    else
        % if case 2 fails, try case 3
        while phi_temp <= phi_prev
            % keep multiplying lambda_prev by v until it satisfies error criteria
            % this is the equivalent of solving for w
            lambda_temp = lambda_temp*v;
            delta_temp = deltacalcEE(q_cur,o,d,a,alpha,d_out,lambda_temp);
            q_temp = q_cur + delta_temp;
            phi_temp = phicalcEE(d_out,W_EE_vec,d,q_temp,a,alpha,o);
            %disp('case 3 used')
        end
        lambdaEE = lambda_temp;
    end
end
end



% calculates end effector position and orientation using the 
% Denavit-Hartenberg Matrix method
function T_EE = DH(a,alpha,d,q,o)
for i = 1:length(a)
    T{i}(1,1) = cos(q(i) + o(i));
    T{i}(1,2) = -cos(alpha(i))*sin(q(i) + o(i));
    T{i}(1,3) = sin(alpha(i))*sin(q(i) + o(i));
    T{i}(1,4) = a(i)*cos(q(i) + o(i));
    
    T{i}(2,1) = sin(q(i) + o(i));
    T{i}(2,2) = cos(alpha(i))*cos(q(i) + o(i));
    T{i}(2,3) = -sin(alpha(i))*cos(q(i) + o(i));
    T{i}(2,4) = a(i)*sin(q(i) + o(i));

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

end
