function dx = tmFunc(t, x, P1, P2,P3, heli, par)

tht0 = interp1(P1{1}, P1{2}, t);

thtP = interp1(P2{1}, P2{2}, t);
% Q=max(0,-1.4916e+03*t+10441)-9.3127e2;

% Q = max(0, -4250*t + 8500);% - 9.3127e2;
% par.p=Q;

% omg=interp1(P4{1},P4{2},t);
% x(5)=omg;
% x = [x(1:4); omg; x(5)];

% Repurposing P3 for cyclic computation
par.thtc = interp1(P3{1}, P3{2}, t);
par.thts = interp1(P3{1}, P3{3}, t);

dx = LDM_state_eqns(x, [tht0; thtP], heli, par, 1);
% dx = dx([1:4,6]);
end