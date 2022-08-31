function [ru,rv,rp] = apprx_conv(S,grd,N)
% Approximate convergence rates ru, rv, rp, for U, V, and P.
% Inputs:
% S   = {U2nh, V2nh, P2nh;
%          . ,  .  ,  .  ;
%        U2h , V2h , P2h ;
%         Uh ,  Vh ,  Ph  } where h is the finest mesh, 2nh is most coarse
% grd = {XU2nh, YU2nh, XV2nh, YV2nh, XP2nh, YP2nh;
%          .  ,   .  ,   .  ,   .  ,   .  ,   .  ;
%        XU2h , YU2h , XV2h , YV2h , XP2h , YP2h ;
%         XUh ,  YUh , XVh  , YVh  , XPh  ,  YPh  }
% N   = [NxL_2nh, NxM_2nh, Nx_2nh, NyL_2nh, NyM_2nh, Ny_2nh;
%           .   ,    .   ,    .  ,    .   ,    .   ,   .   ;
%        NxL_2h , NxM_2h , Nx_2h , NyL_2h , NyM_2h , Ny_2h ;
%         NxL_h ,  NxM_h ,  Nx_h ,  NyL_h ,  NyM_h ,  Ny_h  ]


U = S(:,1); V = S(:,2); P = S(:,3);
XU = grd(:,1); YU = grd(:,2); 
XV = grd(:,3); YV = grd(:,4);
XP = grd(:,5); YP = grd(:,6);


n = length(U) - 2;
ru = zeros(n,1);
rv = zeros(n,1);
rp = zeros(n,1);
for k = 1:n
    NC = N(k,:);
    
    U4h = U{k};
    U2h = interp2(XU{k+1},YU{k+1},U{k+1}, XU{k}, YU{k}, 'spline');
    Uh = interp2(XU{k+2},YU{k+2},U{k+2}, XU{k}, YU{k}, 'spline');
    
    V4h = V{k};
    V2h = interp2(XV{k+1},YV{k+1},V{k+1}, XV{k}, YV{k}, 'spline');
    Vh = interp2(XV{k+2},YV{k+2},V{k+2}, XV{k}, YV{k}, 'spline');
    
    P4h = P{k};
    P2h = interp2(XP{k+1},YP{k+1},P{k+1}, XP{k}, YP{k}, 'spline');
    Ph = interp2(XP{k+2},YP{k+2},P{k+2}, XP{k}, YP{k}, 'spline');
    
    [U4h,V4h,P4h] = reshapeFinal(NC,U4h,V4h,P4h);
    [U2h,V2h,P2h] = reshapeFinal(NC,U2h,V2h,P2h);
    [Uh,Vh,Ph] = reshapeFinal(NC,Uh,Vh,Ph);

    Cu = log( abs( U4h - U2h )./ abs(U2h - Uh) ) /log(2);
    ru(k) = mean(mean(Cu,'omitnan'), 'omitnan');
    
    Cv = log( abs( V4h - V2h )./ abs(V2h - Vh) ) /log(2);
    rv(k) = mean(mean(Cv,'omitnan'), 'omitnan');
    
    Cp = log( abs( P4h - P2h )./ abs(P2h - Ph) ) /log(2);
    rp(k) = mean(mean(Cp,'omitnan'), 'omitnan');
end

ru = [NaN; NaN; ru];
rv = [NaN; NaN; rv];
rp = [NaN; NaN; rp];
