% This function describes the activation dynamics of the ideal arm torque
% actuators.
%
% Author: Antoine Falisse
% Date: 12/19/2018
% 
function dadt = ArmActivationDynamics(e,a)
tau = 0.035;
dadt = (e-a)./tau;

end
