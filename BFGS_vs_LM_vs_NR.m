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

% input joint limits
qlims = [-170 170;
         -120 120;
         -170 170;
         -120 120;
         -170 170;
         -120 120;
         -175 175].*pi/180;

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

% put in arbitrary angle targets, radians
q1_target = 0*pi/180;
q2_target = -7*pi/180;
q3_target = 0*pi/180;
q4_target = -70*pi/180;
q5_target = 0*pi/180;
q6_target = 120*pi/180;
q7_target = 0*pi/180;
q_target = [q1_target q2_target q3_target q4_target q5_target q6_target q7_target];

% set max and min limits for error
errorlims = [-45 45];
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
nxdEE = Td(1,1);
nydEE = Td(2,1);
nzdEE = Td(3,1);
sxdEE = Td(1,2);
sydEE = Td(2,2);
szdEE = Td(3,2);
axdEE = Td(1,3);
aydEE = Td(2,3);
azdEE = Td(3,3);
PxdEE = Td(1,4);
PydEE = Td(2,4);
PzdEE = Td(3,4);
PRd = [nxdEE nydEE nzdEE sxdEE sydEE szdEE axdEE aydEE azdEE PxdEE PydEE PzdEE];

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

% if MATLAB produces warnings for ill-conditioned matrices, try changing 
% lambda0 and v 
phi_max = 1e-20;     % define acceptable error
v = 1.1;             % scaling factor, determined experimentally
lambda0 = 5e-5;      % initial estimate for lambda


epsilon = 1e-12;   % define acceptable error
alpha0 = 1;        % intial guess for step size
beta = 0.5;        % step size scaling factor
% values for c1, c2
% note: 0 < c1 < c2 < 1
c1 = 1e-4;

epsilon_NR = 1e-12; % acceptable error for NR method

H = eye(7);

% calculate q using the BFGS method
[q_BFGS, k_BFGS, qcheck_BFGS,BFGS_time] = IK_BFGS_iiwa(q0',H,epsilon,alpha0,beta,c1,W,PRd,d,a,alphar,o);
qd_BFGS = rad2deg(q_BFGS);
qdcheck_BFGS = rad2deg(qcheck_BFGS);

% calculate q using the LM method
[q_LM, k_LM, qcheck_LM, LM_time] = IK_LM_iiwa(a,alphar,d,q0',o,PRd,phi_max,lambda0,v,W);
qd_LM = rad2deg(q_LM);
qdcheck_LM = rad2deg(qcheck_LM);

% calculate q using the NR method
[q_NR, k_NR, qcheck_NR, NR_time] = IK_NR_iiwa(q0',epsilon_NR,PRd,d,a,alphar,o);
qd_NR = rad2deg(q_NR);
qdcheck_NR = rad2deg(qcheck_NR);

tickfontsize = 16;
figure(1)
subplot(1,3,1)
set(gcf,'position',[200,100,2200,1200]);
plot(k_BFGS,qdcheck_BFGS(1,:),k_BFGS,qdcheck_BFGS(2,:),k_BFGS,qdcheck_BFGS(3,:),k_BFGS,qdcheck_BFGS(4,:),k_BFGS,qdcheck_BFGS(5,:),k_BFGS,qdcheck_BFGS(6,:),k_BFGS,qdcheck_BFGS(7,:),'LineWidth',2)
ax = gca;
ax.FontSize = tickfontsize; 
title('BFGS method','FontSize',22)
%legend('q_1','q_2','q_3','q_4','q_5','q_6','q_7','FontSize',16,'Location','northeastoutside')
%xlabel('iteration','FontSize',22)
ylabel('angle, degrees','FontSize',22)
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));

subplot(1,3,2)
plot(k_LM,qdcheck_LM(1,:),k_LM,qdcheck_LM(2,:),k_LM,qdcheck_LM(3,:),k_LM,qdcheck_LM(4,:),k_LM,qdcheck_LM(5,:),k_LM,qdcheck_LM(6,:),k_LM,qdcheck_LM(7,:),'LineWidth',2)
ax = gca;
ax.FontSize = tickfontsize; 
title('LM method','FontSize',22)
xlabel('iteration','FontSize',22)
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));

subplot(1,3,3)
plot(k_NR,qdcheck_NR(1,:),k_NR,qdcheck_NR(2,:),k_NR,qdcheck_NR(3,:),k_NR,qdcheck_NR(4,:),k_NR,qdcheck_NR(5,:),k_NR,qdcheck_NR(6,:),k_NR,qdcheck_NR(7,:),'LineWidth',2)
ax = gca;
ax.FontSize = tickfontsize; 
title('NR method','FontSize',22)
legend('q_1','q_2','q_3','q_4','q_5','q_6','q_7','FontSize',16,'Location','northeastoutside')
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));

% display initial estimate error
initial_error_summary = ['Initial estimate error: '];
initial_error = num2str(rad2deg(q_target - q0')');

% calculate total error, display computation time for both methods
BFGS_error = sum(sum(abs(DHcalc(a,alphar,d,q_target,o) - DHcalc(a,alphar,d,q_BFGS,o))));
BFGS_summary = ['BFGS method computation time: ',num2str(1000*BFGS_time),' milliseconds'];
BFGS_summary2 = ['BFGS method total error: ',num2str(BFGS_error)];

LM_error = sum(sum(abs(DHcalc(a,alphar,d,q_target,o) - DHcalc(a,alphar,d,q_LM,o))));
LM_summary = ['LM method computation time: ',num2str(1000*LM_time),' milliseconds'];
LM_summary2 = ['LM method total error: ',num2str(LM_error)];

NR_error = sum(sum(abs(DHcalc(a,alphar,d,q_target,o) - DHcalc(a,alphar,d,q_NR,o))));
NR_summary = ['NR method computation time: ',num2str(1000*NR_time),' milliseconds'];
NR_summary2 = ['NR method total error: ',num2str(NR_error)];

disp(initial_error_summary)
disp(initial_error)
fprintf('\n')
disp(BFGS_summary)
disp(LM_summary)
disp(NR_summary)
fprintf('\n')
disp(BFGS_summary2)
disp(LM_summary2)
disp(NR_summary2)
fprintf('\n')

% compare computation times, output summary
if BFGS_time > LM_time
    ratio = BFGS_time/LM_time;
    summary = ['BFGS method took ',num2str(ratio),' times longer than LM method'];
    disp(summary)
else
    ratio = LM_time/BFGS_time;
    summary = ['LM method took ',num2str(ratio),' times longer than BFGS method'];
    disp(summary)
end

if BFGS_time > NR_time
    ratio = BFGS_time/NR_time;
    summary = ['BFGS method took ',num2str(ratio),' times longer than NR method'];
    disp(summary)
else
    ratio = NR_time/BFGS_time;
    summary = ['NR method took ',num2str(ratio),' times longer than BFGS method'];
    disp(summary)
end

if LM_time > NR_time
    ratio = LM_time/NR_time;
    summary = ['LM method took ',num2str(ratio),' times longer than NR method'];
    disp(summary)
else
    ratio = NR_time/LM_time;
    summary = ['NR method took ',num2str(ratio),' times longer than LM method'];
    disp(summary)
end

% convert final angles to degrees
qd_target = rad2deg(q_target)';

results = table(qd_target,qd_BFGS,qd_LM,qd_NR);
fprintf('\n')
disp('Servo angles at final iteration:')
fprintf('\n')
disp(results)

[RSSR_BFGS,RSSP_BFGS,n_BFGS] = RSScalc(PRd',qcheck_BFGS,a,alphar,d,o);
[RSSR_LM,RSSP_LM,n_LM] = RSScalc(PRd',qcheck_LM,a,alphar,d,o);
[RSSR_NR,RSSP_NR,n_NR] = RSScalc(PRd',qcheck_NR,a,alphar,d,o);

figure(2)
set(gcf,'position',[200,100,2200,1200]);
subplot(3,1,1)
plot(n_BFGS,RSSR_BFGS,'x--',n_BFGS,RSSP_BFGS,'*--','LineWidth',1)
ax = gca;
ax.FontSize = tickfontsize; 
yline(0)
title('BFGS Method: Residual Sum of Squares','FontSize',22)
ylabel('Squared residual','FontSize',16)
xlabel('Iteration','FontSize',16)
legend('Rotation ','Position ','FontSize',16,'Location','northeast')
BFGS_lims = ceil(max(max(RSSR_BFGS),max(RSSP_BFGS)));
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
ylim([-0.5,BFGS_lims])


subplot(3,1,2)
plot(n_LM,RSSR_LM,'x--',n_LM,RSSP_LM,'*--','LineWidth',1)
ax = gca;
ax.FontSize = tickfontsize; 
yline(0)
title('LM Method: Residual Sum of Squares','FontSize',22)
ylabel('Squared residual','FontSize',16)
xlabel('Iteration','FontSize',16)
legend('Rotation ','Position ','FontSize',16,'Location','northeast')
LM_lims = ceil(max(max(RSSR_LM),max(RSSP_LM)));
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
ylim([-0.5,LM_lims])

subplot(3,1,3)
plot(n_NR,RSSR_NR,'x--',n_NR,RSSP_NR,'*--','LineWidth',1)
ax = gca;
ax.FontSize = tickfontsize; 
yline(0)
title('NR Method: Residual Sum of Squares','FontSize',22)
ylabel('Squared residual','FontSize',16)
xlabel('Iteration','FontSize',16)
legend('Rotation ','Position ','FontSize',16,'Location','northeast')
NR_lims = ceil(max(max(RSSR_NR),max(RSSP_NR)));
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
ylim([-0.5,NR_lims])

% Verification
% verify results using the MATLAB robotics toolbox
% https://www.mathworks.com/help/robotics/ug/build-a-robot-step-by-step.html

L_EE = 0; % end effector length
% d values for plot
d = [0.36 0 0.42 0 0.4 0 0.1199+L_EE];
% axis translation offset values, m
a = [0 0 0 0 0 0 0];

% alpha offset values, rad
al = [pi/2 -pi/2 -pi/2 pi/2 pi/2 -pi/2 0];

robotExact = buildRobot(q_target,a,al,d);
BFGS_robot = buildRobot(q_BFGS,a,al,d);
LM_robot = buildRobot(q_LM,a,al,d);
NR_robot = buildRobot(q_NR,a,al,d);

% plot both solutions, with elbow + EE specification vs. only EE
% specification
figure(3)
set(gcf,'position',[200,100,1800,1200]);
show(BFGS_robot);
hold on
show(LM_robot);
hold on
show(NR_robot);
hold on
show(robotExact);
ax = gca;
ax.FontSize = tickfontsize; 

% set limits for robot plot
zlim([0,1])
xlim([-0.5 0.5])
ylim([-0.5 0.5])

%% Alternative plots
% use these to put in papers, etc.
figure(4)
set(gcf,'position',[200,100,2200,1200]);
plot(k_BFGS,qdcheck_BFGS(1,:),k_BFGS,qdcheck_BFGS(2,:),k_BFGS,qdcheck_BFGS(3,:),k_BFGS,qdcheck_BFGS(4,:),k_BFGS,qdcheck_BFGS(5,:),k_BFGS,qdcheck_BFGS(6,:),k_BFGS,qdcheck_BFGS(7,:),'LineWidth',2)
title('BFGS method','FontSize',22)
legend('q_1','q_2','q_3','q_4','q_5','q_6','q_7','FontSize',16,'Location','northeastoutside')
xlabel('iteration','FontSize',22)
ylabel('angle, rads','FontSize',22)
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));

figure(5)
set(gcf,'position',[200,100,2200,1200]);
plot(k_LM,qdcheck_LM(1,:),k_LM,qdcheck_LM(2,:),k_LM,qdcheck_LM(3,:),k_LM,qdcheck_LM(4,:),k_LM,qdcheck_LM(5,:),k_LM,qdcheck_LM(6,:),k_LM,qdcheck_LM(7,:),'LineWidth',2)
title('LM method','FontSize',22)
legend('q_1','q_2','q_3','q_4','q_5','q_6','q_7','FontSize',16,'Location','northeastoutside')
xlabel('iteration','FontSize',22)
ylabel('angle, rads','FontSize',22)
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));

figure(6)
set(gcf,'position',[200,100,2200,1200]);
plot(k_NR,qdcheck_NR(1,:),k_NR,qdcheck_NR(2,:),k_NR,qdcheck_NR(3,:),k_NR,qdcheck_NR(4,:),k_NR,qdcheck_NR(5,:),k_NR,qdcheck_NR(6,:),k_NR,qdcheck_NR(7,:),'LineWidth',2)
title('NR method','FontSize',22)
legend('q_1','q_2','q_3','q_4','q_5','q_6','q_7','FontSize',16,'Location','northeastoutside')
xlabel('iteration','FontSize',22)
ylabel('angle, rads','FontSize',22)
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
