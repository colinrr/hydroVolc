function [value,isterminal,direction] = stopHybrid(y,numGS,u_0)
% CR 2021 - I/O changes for GSD

m_d       = y(3);
m_v       = y(4);
m_l       = y(5);
m_s       = sum(y(6:6+numGS-1));
psi       = y(6+numGS);
angle     = y(7+numGS);
m         = m_d + m_v + m_l + m_s;
u         = psi/m;

% Locate the height when velocity passes through zero
value = u*sin(angle)-0.1;     % Detect velocity = 0
    %the 0.001 helps for stability
isterminal = 1;   % Stop the integration
direction = 0;   % detect all zeros (default)