function r = EOC_aprx(U,XU,YU)
% Approximate convergence rate r
% Inputs:
% U = {U_4h; U_2h; U_h} where h is the finest mesh, 4h is most coarse
% XU = {XU_4h; XU_2h; XU_h} where each XU is from meshgrid(xu,yu)
% YU = {YU_4h; YU_2h; YU_h}
n = length(U) - 2;
r = zeros(n,1);
for k = 1:n
    U4h = U{k};
    U2h = interp2(XU{k+1},YU{k+1},U{k+1}, XU{k}, YU{k}, 'spline');
    Uh = interp2(XU{k+2},YU{k+2},U{k+2}, XU{k}, YU{k}, 'spline');
    Cu = log( abs( U4h - U2h )./ abs(U2h - Uh) ) /log(2);
    r(k) = mean(mean(Cu,'omitnan'), 'omitnan');
end
