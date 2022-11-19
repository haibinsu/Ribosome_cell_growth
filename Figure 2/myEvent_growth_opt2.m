function [value, isterminal, direction] = myEvent_growth_opt2(t, y, Mth, par)

% value = y(1) - 2*R0 ;
NR = par(4);
% MR = 2.454545454545455e+04;
NP = par(5);

value = NR*y(1) + NP*y(2) + y(3) - 2*Mth;
% value = MR*y(1) + NP*y(2) + y(3) - 2*Mth;

isterminal = 1;   % Stop the integration
direction  = 0;

end