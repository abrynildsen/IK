function [RSSR, RSSP, n] = RSScalc(Td,q,a,alphar,d,offset)
% calculates residual sum of squares for position and orientation of the
% end effector over iterations

% INPUTS
% Td - desired end effector DH configuration, input as a 12x1 vector
% q - servo input angle vectors input as a 7xi matrix where i is the
% iteration
% d - lengths of the robot arms
% a - translational offsets of the robot joint axes, 7x1 vector
% alphar - angular offsets of the robot joint axes, 7x1 vector
% o - offset of the robot joints, used for calibration only, rads

% OUTPUTS
% RSSP - residual sum of squares for end effector position
% RSSR - residual sum of squares for end effector orientation
% n - number of iterations

[~,c] = size(q);
for i = 1:c
    % calculate DH configuration for step 
    Tc_mat = DHcalc(a,alphar,d,q(:,i),offset);
    % convert to a vector
    Tc_vec = reshape(Tc_mat(1:3,:),[12,1]);
    % calculate residual for step
    e = Td - Tc_vec;
    % square the residuals
    Rsqr = e.^2;
    % sum square of residuals for orientation
    RSSR(i) = sum(Rsqr(1:9));
    % sum square of residuals for position
    RSSP(i) = sum(Rsqr(10:12));
end
n = 1:c;

end