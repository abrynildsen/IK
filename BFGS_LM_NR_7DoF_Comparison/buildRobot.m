function robot = buildRobot(q,a,al,d)
% uses the MATLAB robotics toolbox to create a visual of the Kuka iiwa
% configuration

% create the robot rigidBodyTree
robot = rigidBodyTree;

% create the base
arm0 = rigidBody('arm0');
jnt0 = rigidBodyJoint('jnt0','fixed'); 
setFixedTransform(jnt0,[0 0 0 pi],"dh");
arm0.Joint = jnt0;
addBody(robot,arm0,'base'); % add arm 0 to base

% create the 1st arm
arm1 = rigidBody('arm1');
jnt1 = rigidBodyJoint('jnt1','revolute');
jnt1.HomePosition = q(1);  
setFixedTransform(jnt1,[a(1) al(1) d(1) 0],"dh");
arm1.Joint = jnt1;
addBody(robot,arm1,'arm0'); % add arm 1 to arm 0

% create the 2nd arm
arm2 = rigidBody('arm2');
jnt2 = rigidBodyJoint('jnt2','revolute');
jnt2.HomePosition = q(2);  
setFixedTransform(jnt2,[a(2) al(2) d(2) 0],"dh");
arm2.Joint = jnt2;
addBody(robot,arm2,'arm1'); % add arm 2 to arm 1

% create the 3rd arm
arm3 = rigidBody('arm3');
jnt3 = rigidBodyJoint('jnt3','revolute');
jnt3.HomePosition = q(3);          % rotation about z-axis
setFixedTransform(jnt3,[a(3) al(3) d(3) 0],"dh");
arm3.Joint = jnt3;
addBody(robot,arm3,'arm2'); % add arm 3 to arm 2

% create the 4th arm
arm4 = rigidBody('arm4');
jnt4 = rigidBodyJoint('jnt4','revolute');
jnt4.HomePosition = q(4);          % rotation about z-axis
setFixedTransform(jnt4,[a(4) al(4) d(4) 0],"dh");
arm4.Joint = jnt4;
addBody(robot,arm4,'arm3'); % add arm 4 to arm 3

% create the 5th arm
arm5 = rigidBody('arm5');
jnt5 = rigidBodyJoint('jnt5','revolute');
jnt5.HomePosition = q(5);          % rotation about z-axis
setFixedTransform(jnt5,[a(5) al(5) d(5) 0],"dh");
arm5.Joint = jnt5;
addBody(robot,arm5,'arm4'); % add arm 5 to arm 4

% create the 6th arm
arm6 = rigidBody('arm6');
jnt6 = rigidBodyJoint('jnt6','revolute');
jnt6.HomePosition = q(6);          % rotation about z-axis
setFixedTransform(jnt6,[a(6) al(6) d(6) 0],"dh");
arm6.Joint = jnt6;
addBody(robot,arm6,'arm5'); % add arm 6 to arm 5

% create the 7th arm
arm7 = rigidBody('arm7');
jnt7 = rigidBodyJoint('jnt7','revolute');
jnt7.HomePosition = q(7);          % rotation about z-axis
setFixedTransform(jnt7,[a(7) al(7) d(7) 0],"dh");
arm7.Joint = jnt7;
addBody(robot,arm7,'arm6'); % add arm 7 to arm 6

% create the end effector
EE = rigidBody('EE');
jnt8 = rigidBodyJoint('jnt8','revolute');
setFixedTransform(jnt8,[0 0 0 0],"dh");
EE.Joint = jnt8;
addBody(robot,EE,'arm7'); % add end effector to arm 7
