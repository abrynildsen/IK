clc;clear;close all

% parameters
L_EE = 0;            % end effector length
d1 = 0.36;           % length of robot arm 1
d2 = 0;              % length of robot arm 2
d3 = 0.42;           % length of robot arm 3
d4 = 0;              % length of robot arm 4
d5 = 0.4;            % length of robot arm 5
d6 = 0;              % length of robot arm 6
d7 = 0.1199+L_EE;    % length of robot arm 7
d = [d1 d2 d3 d4 d5 d6 d7];

% axis translation offset values, m
a1 = 0;   % axis offset of arm 1
a2 = 0;   % axis offset of arm 1
a3 = 0;   % axis offset of arm 1
a4 = 0;   % axis offset of arm 1
a5 = 0;   % axis offset of arm 1
a6 = 0;   % axis offset of arm 1
a7 = 0;   % axis offset of arm 1
a = [a1 a2 a3 a4 a5 a6 a7];

% alpha offset values, rad
alpha1 = pi/2;   % angle offset of axis 1
alpha2 = -pi/2;  % angle offset of axis 2
alpha3 = -pi/2;  % angle offset of axis 3
alpha4 = pi/2;   % angle offset of axis 4
alpha5 = pi/2;   % angle offset of axis 5
alpha6 = -pi/2;  % angle offset of axis 6
alpha7 = 0;      % angle offset of axis 7
alphar = [alpha1 alpha2 alpha3 alpha4 alpha5 alpha6 alpha7];

% coordinate frame axis calibration offset values, rad
o1 = pi; % axis 1 offset
o2 = 0;  % axis 2 offset
o3 = 0;  % axis 3 offset
o4 = 0;  % axis 4 offset
o5 = 0;  % axis 5 offset
o6 = 0;  % axis 6 offset
o7 = 0;  % axis 7 offset
o = [o1 o2 o3 o4 o5 o6 o7];

% joint limits for the LBR iiwa
qlims = [-170 170;
         -120 120;
         -170 170;
         -120 120;
         -170 170;
         -120 120;
         -175 175];

% end effector weights
WnxEE = 1;  % weight for end effector orientation nx component 
WnyEE = 1;  % weight for end effector orientation ny component 
WnzEE = 1;  % weight for end effector orientation nz component 
WsxEE = 1;  % weight for end effector orientation sx component 
WsyEE = 1;  % weight for end effector orientation sy component 
WszEE = 1;  % weight for end effector orientation sz component 
WaxEE = 1;  % weight for end effector orientation ax component 
WayEE = 1;  % weight for end effector orientation ay component 
WazEE = 1;  % weight for end effector orientation az component 
WPxEE = 1;  % weight for end effector position x component 
WPyEE = 1;  % weight for end effector position y component 
WPzEE = 1;  % weight for end effector position z component 
W = [WnxEE WnyEE WnzEE WsxEE WsyEE WszEE WaxEE WayEE WazEE WPxEE WPyEE WPzEE];

% set max and min limits for error
errorlims = [-45 45];

% set number of iterations
num_its = 10;

% if MATLAB produces warnings for ill-conditioned matrices, try changing 
% lambda0 and v 
phi_max = 1e-10;     % define acceptable error
v = 1.1;             % scaling factor, determined experimentally
lambda0 = 5e-5;      % initial estimate for lambda


epsilon = 1e-6;   % define acceptable error
alpha0 = 1;        % intial guess for step size
beta = 0.5;        % step size scaling factor
% values for c1, c2
% note: 0 < c1 < c2 < 1
c1 = 1e-4;

epsilon_NR = 1e-6; % acceptable error for NR method

% initial estimate of hessian
H0 = eye(7);


% allocate empty matrices
BFGS_time = zeros(1,num_its);
LM_time = zeros(1,num_its);
NR_time = zeros(1,num_its);

BFGS_its = zeros(1,num_its);
LM_its = zeros(1,num_its);
NR_its = zeros(1,num_its);

RSSRf_BFGS = zeros(1,num_its);
RSSPf_BFGS = zeros(1,num_its);
RSSRf_LM = zeros(1,num_its);
RSSPf_LM = zeros(1,num_its);
RSSRf_NR = zeros(1,num_its);
RSSPf_NR = zeros(1,num_its);


tic
for n = 1:num_its
    % put in arbitrary angle targets, radians
    q1_target = randi(qlims(1,:))*pi/180;
    q2_target = randi(qlims(2,:))*pi/180;
    q3_target = randi(qlims(3,:))*pi/180;
    q4_target = randi(qlims(4,:))*pi/180;
    q5_target = randi(qlims(5,:))*pi/180;
    q6_target = randi(qlims(6,:))*pi/180;
    q7_target = randi(qlims(7,:))*pi/180;
    q_target = [q1_target q2_target q3_target q4_target q5_target q6_target q7_target];
    qtarget(:,n) = q_target';
    
    % give the angle an initial guess, add some error
    q1_cur = q1_target + randi(errorlims)*pi/180;
    q2_cur = q2_target + randi(errorlims)*pi/180;
    q3_cur = q3_target + randi(errorlims)*pi/180;
    q4_cur = q4_target + randi(errorlims)*pi/180;
    q5_cur = q5_target + randi(errorlims)*pi/180;
    q6_cur = q6_target + randi(errorlims)*pi/180;
    q7_cur = q7_target + randi(errorlims)*pi/180;
    q0 = sat([q1_cur q2_cur q3_cur q4_cur q5_cur q6_cur q7_cur],qlims);
    
    % set up desired orientation and position for end effector
    Td = DHcalc(a,alphar,d,q_target,o);
    PRd = reshape(Td(1:3,:),[1,12]);
    
    % calculate q using the BFGS method
    [q_BFGS(:,n), k_BFGS, qcheck_BFGS,BFGS_time(n)] = IK_BFGS_iiwa(q0',H0,epsilon,alpha0,beta,c1,W,PRd,d,a,alphar,o);
    qd_BFGS = rad2deg(q_BFGS);
    qdcheck_BFGS = rad2deg(qcheck_BFGS);
    
    % calculate q using the LM method
    [q_LM(:,n), k_LM, qcheck_LM, LM_time(n)] = IK_LM_iiwa(a,alphar,d,q0',o,PRd,phi_max,lambda0,v,W);
    qd_LM = rad2deg(q_LM);
    qdcheck_LM = rad2deg(qcheck_LM);
    
    % calculate q using the NR method
    [q_NR(:,n), k_NR, qcheck_NR, NR_time(n)] = IK_NR_iiwa(q0',epsilon_NR,PRd,d,a,alphar,o);
    qd_NR = rad2deg(q_NR);
    qdcheck_NR = rad2deg(qcheck_NR);
    
    % determine residual sum of squares for each algorithm
    [RSSR_BFGS,RSSP_BFGS,n_BFGS] = RSScalc(PRd',qcheck_BFGS,a,alphar,d,o);
    [RSSR_LM,RSSP_LM,n_LM] = RSScalc(PRd',qcheck_LM,a,alphar,d,o);
    [RSSR_NR,RSSP_NR,n_NR] = RSScalc(PRd',qcheck_NR,a,alphar,d,o);
    
    % save number of iterations for each algorithm
    BFGS_its(n) = length(k_BFGS);
    LM_its(n) = length(k_LM);
    NR_its(n) = length(k_NR);
    
    % save RSS at final iteration for each algorithm
    RSSRf_BFGS(n) = RSSR_BFGS(end);
    RSSPf_BFGS(n) = RSSP_BFGS(end);
    RSSRf_LM(n) = RSSR_LM(end);
    RSSPf_LM(n) = RSSP_LM(end);
    RSSRf_NR(n) = RSSR_NR(end);
    RSSPf_NR(n) = RSSP_NR(end);
end


clc
x_tickmarks = 1:num_its;
figure(1)
set(gcf,'position',[200,100,800,600]);
subplot(3,1,1)
plot(x_tickmarks, BFGS_its,'*')
ax = gca;
ax.FontSize = 12; 
title('BFGS','FontSize',22)
subplot(3,1,2)
plot(x_tickmarks, LM_its,'*')
ax = gca;
ax.FontSize = 12; 
title('LM','FontSize',22)
ylabel('Number of Iterations','FontSize',22)
subplot(3,1,3)
plot(x_tickmarks, NR_its,'*')
ax = gca;
ax.FontSize = 12; 
title('NR','FontSize',22)
xlabel('Trial Number','FontSize',22)

figure(2)
set(gcf,'position',[200,100,800,600]);
subplot(3,1,1)
plot(x_tickmarks, BFGS_time,'*')
ax = gca;
ax.FontSize = 12; 
title('BFGS','FontSize',22)
subplot(3,1,2)
plot(x_tickmarks, LM_time,'*')
ax = gca;
ax.FontSize = 12; 
title('LM','FontSize',22)
ylabel('Convergence Time (s)','FontSize',22)
subplot(3,1,3)
plot(x_tickmarks, NR_time,'*')
ax = gca;
ax.FontSize = 12; 
title('NR','FontSize',22)
xlabel('Trial Number','FontSize',22)

figure(3)
set(gcf,'position',[200,100,800,600]);
subplot(3,1,1)
plot(x_tickmarks, RSSRf_BFGS,'*')
ax = gca;
ax.FontSize = 12; 
title('BFGS','FontSize',22)
subplot(3,1,2)
plot(x_tickmarks, RSSRf_LM,'*')
ax = gca;
ax.FontSize = 12; 
title('LM','FontSize',22)
ylabel('Orientation RSS','FontSize',22)
subplot(3,1,3)
plot(x_tickmarks, RSSRf_NR,'*')
ax = gca;
ax.FontSize = 12; 
title('NR','FontSize',22)
xlabel('Trial Number','FontSize',22)

figure(4)
set(gcf,'position',[200,100,800,600]);
subplot(3,1,1)
plot(x_tickmarks, RSSPf_BFGS,'*')
ax = gca;
ax.FontSize = 12; 
title('BFGS','FontSize',22)
subplot(3,1,2)
plot(x_tickmarks, RSSPf_LM,'*')
ax = gca;
ax.FontSize = 12; 
title('LM','FontSize',22)
ylabel('Position RSS','FontSize',22)
subplot(3,1,3)
plot(x_tickmarks, RSSPf_NR,'*')
ax = gca;
ax.FontSize = 12; 
title('NR','FontSize',22)
xlabel('Trial Number','FontSize',22)

% statistical analysis
BFGS_summary = ['BFGS: ',num2str(1000*mean(BFGS_time)),' milliseconds'];
BFGS_summary2 = ['BFGS: ',num2str(1000*std(BFGS_time)),' milliseconds'];
LM_summary = ['LM: ',num2str(1000*mean(LM_time)),' milliseconds'];
LM_summary2 = ['LM: ',num2str(1000*std(LM_time)),' milliseconds'];
NR_summary = ['NR: ',num2str(1000*mean(NR_time)),' milliseconds'];
NR_summary2 = ['NR: ',num2str(1000*std(NR_time)),' milliseconds'];

BFGS_summary3 = ['BFGS: ',num2str(mean(RSSPf_BFGS))];
BFGS_summary4 = ['BFGS: ',num2str(std(RSSPf_BFGS))];
LM_summary3 = ['LM: ',num2str(mean(RSSPf_LM))];
LM_summary4 = ['LM: ',num2str(std(RSSPf_LM))];
NR_summary3 = ['NR: ',num2str(mean(RSSPf_NR))];
NR_summary4 = ['NR: ',num2str(std(RSSPf_NR))];

BFGS_summary5 = ['BFGS: ',num2str(mean(RSSRf_BFGS))];
BFGS_summary6 = ['BFGS: ',num2str(std(RSSRf_BFGS))];
LM_summary5 = ['LM: ',num2str(mean(RSSRf_LM))];
LM_summary6 = ['LM: ',num2str(std(RSSRf_LM))];
NR_summary5 = ['NR: ',num2str(mean(RSSRf_NR))];
NR_summary6 = ['NR: ',num2str(std(RSSRf_NR))];

BFGS_summary7 = ['BFGS: ',num2str(mean(BFGS_its))];
BFGS_summary8 = ['BFGS: ',num2str(std(BFGS_its))];
LM_summary7 = ['LM: ',num2str(mean(LM_its))];
LM_summary8 = ['LM: ',num2str(std(LM_its))];
NR_summary7 = ['NR: ',num2str(mean(NR_its))];
NR_summary8 = ['NR: ',num2str(std(NR_its))];

fprintf('AVERAGE COMPUTATION TIME')
fprintf('\n')
fprintf('------------------------')
fprintf('\n')
disp(BFGS_summary)
disp(LM_summary)
disp(NR_summary)
fprintf('\n')
fprintf('COMPUTATION TIME STANDARD DEVIATION')
fprintf('\n')
fprintf('------------------------')
fprintf('\n')
disp(BFGS_summary2)
disp(LM_summary2)
disp(NR_summary2)
fprintf('\n')

fprintf('POSITION: AVERAGE RESIDUAL SUM OF SQUARES')
fprintf('\n')
fprintf('------------------------')
fprintf('\n')
disp(BFGS_summary3)
disp(LM_summary3)
disp(NR_summary3)
fprintf('\n')
fprintf('POSITION: RESIDUAL SUM OF SQUARES STANDARD DEIVIATION')
fprintf('\n')
fprintf('------------------------')
fprintf('\n')
disp(BFGS_summary4)
disp(LM_summary4)
disp(NR_summary4)
fprintf('\n')

fprintf('ORIENTATION: AVERAGE RESIDUAL SUM OF SQUARES')
fprintf('\n')
fprintf('------------------------')
fprintf('\n')
disp(BFGS_summary5)
disp(LM_summary5)
disp(NR_summary5)
fprintf('\n')
fprintf('ORIENTATION: RESIDUAL SUM OF SQUARES STANDARD DEIVIATION')
fprintf('\n')
fprintf('------------------------')
fprintf('\n')
disp(BFGS_summary6)
disp(LM_summary6)
disp(NR_summary6)
fprintf('\n')

fprintf('AVERAGE NUMBER OF ITERATIONS')
fprintf('\n')
fprintf('------------------------')
fprintf('\n')
disp(BFGS_summary7)
disp(LM_summary7)
disp(NR_summary7)
fprintf('\n')
fprintf('STANDARD DEVIATION: NUMBER OF ITERATIONS')
fprintf('\n')
fprintf('------------------------')
fprintf('\n')
disp(BFGS_summary8)
disp(LM_summary8)
disp(NR_summary8)
fprintf('\n')

% plot all configurations overlaid just for laughs
% verify results using the MATLAB robotics toolbox
% https://www.mathworks.com/help/robotics/ug/build-a-robot-step-by-step.html

L_EE = 0; % end effector length
% d values for plot
d = [0.36 0 0.42 0 0.4 0 0.1199+L_EE];
% axis translation offset values, m
a = [0 0 0 0 0 0 0];

% alpha offset values, rad
al = [pi/2 -pi/2 -pi/2 pi/2 pi/2 -pi/2 0];


% plot both solutions, with elbow + EE specification vs. only EE
% specification
figure(5)
set(gcf,'position',[200,100,1800,1200]);
for n = 1:num_its
    robotExact = buildRobot(qtarget(:,n),a,al,d);
    show(robotExact);
    ax = gca;
    ax.FontSize = 12; 
    hold on
end

% set limits for robot plot
zlim([-0.5,1.5])
xlim([-1 1])
ylim([-1 1])


