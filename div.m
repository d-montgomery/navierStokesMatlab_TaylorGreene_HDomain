function D = div(N,U,V,XU,YV)
% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

D = zeros(Ny,Nx);

dx = XU(2:Ny+1 , 2:Nx+1) - XU(2:Ny+1 , 1:Nx);
dy = YV(2:Ny+1 , 2:Nx+1) - YV(1:Ny , 2:Nx+1);

ue = U(2:Ny+1 , 2:Nx+1); uw = U(2:Ny+1 , 1:Nx);
vn = V(2:Ny+1 , 2:Nx+1); vs = V(1:Ny , 2:Nx+1);

% Wash Channel
D(1:Ny,1:NxL) = (ue(1:Ny,1:NxL) - uw(1:Ny,1:NxL))./dx(1:Ny,1:NxL)...
    + (vn(1:Ny,1:NxL) - vs(1:Ny,1:NxL))./dy(1:Ny,1:NxL);

% Injury Channel
D(NyL+1:NyM,NxL+1:NxM) = (ue(NyL+1:NyM,NxL+1:NxM) - uw(NyL+1:NyM,NxL+1:NxM))./dx(NyL+1:NyM,NxL+1:NxM)...
    + (vn(NyL+1:NyM,NxL+1:NxM) - vs(NyL+1:NyM,NxL+1:NxM))./dy(NyL+1:NyM,NxL+1:NxM);

% Blood Channel
D(1:Ny,NxM+1:Nx) = (ue(1:Ny,NxM+1:Nx) - uw(1:Ny,NxM+1:Nx))./dx(1:Ny,NxM+1:Nx)...
    + (vn(1:Ny,NxM+1:Nx) - vs(1:Ny,NxM+1:Nx))./dy(1:Ny,NxM+1:Nx);