function [v, t, d] = myE(t, x, P2)%, lnsP, tmsP)
    thtP = interp1(P2{1}, P2{2}, t);
    v = x(2) + (74.956*thtP*180/pi + 1502.6)*1e-3;
    t = 1;
    d = 1;
end