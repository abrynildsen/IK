function [q,k,qcheck,time] = IK_NR_iiwa(q0,epsilon,PRd,d,a,alphar,o)
% calculates the input servo angles required to achieve a desired DH
% configuration for the 7 axis Kuka iiwa using the NR method

% INPUTS
% q0 - initial estimate for servo input angles, rads
% epsilon - acceptable error 
% PRd - desired DH configuration for the Kuka iiwa end effector
% d - lengths of the robot arms
% a - translational offsets of the robot joint axes, 7x1 vector
% alphar - angular offsets of the robot joint axes, 7x1 vector
% o - offset of the robot joints, used for calibration only, rads

% OUTPUTS
% q - servo angles required to achieve desired end effector position and
% orientation
% k - number of iterations
% qcheck - q values for each iteration

%----------------------------------------------------------------------

% ITERATION 1
q_cur = q0';
qcheck(:,1) = q_cur;
k = 2;
tic
% calculate Jacobian pseudoinverse
J = Jcalc(q_cur,d);
Jp = pinv(J'*J)*J';

% calculate current step DH configuration
T_cur = DH(a,alphar,d,q_cur,o);

% convert to vector
PR_cur = reshape(T_cur(1:3,:),[12,1]);


% calculate residual
e = PRd' - PR_cur;

% calculate delta q
dq = Jp*e;

% update q for step k+1
q_new = q_cur + dq;

% calculate updated DH configuration
T_new = DH(a,alphar,d,q_new,o);

% convert to vector
PR_new = reshape(T_new(1:3,:),[12,1]);

% calculate current step DH configuration
T_cur = DH(a,alphar,d,q_cur,o);
etest = epsCalc(T_cur,T_new);


% update q
q_cur = q_new;


% ITERATIONS 2-k
while etest > epsilon
    % calculate Jacobian pseudoinverse
    J = Jcalc(q_cur,d);
    Jp = pinv(J'*J)*J';% J'*J\J';
    
    % calculate current step DH configuration
    T_cur = DH(a,alphar,d,q_cur,o);
    
    % convert to vector
    PR_cur = reshape(T_cur(1:3,:),[12,1]);
    
    % calculate residual
    e = PRd' - PR_cur;
    
    % calculate delta q
    dq = Jp*e;
    
    % update q for step k+1
    q_new = q_cur + dq;
    
    % calculate updated DH configuration
    T_new = DH(a,alphar,d,q_new,o);
    
    % convert to vector
    PR_new = reshape(T_new(1:3,:),[12,1]);
    
    % calculate epsilon
    etest = epsCalc(T_cur,T_new);
    
    % update q
    q_cur = q_new;
    
    % save current q value for plot
    qcheck(:,k) = wrapToPi(q_cur);
    
    % update step
    k = k + 1;
end
k = 1:k-1;
q = wrapToPi(q_cur);

time = toc;

function eps = epsCalc(T_cur,T_new)
eps = norm(T_cur - T_new)/norm(T_cur);
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

function T_EE = DH(a,alpha,d,theta,offset)
% calculates end effector position and orientation using the 
% Denavit-Hartenberg Matrix method
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


end