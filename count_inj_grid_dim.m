% Domain H within [0, Wleft +Wmid +Wright] x [0, Hlow +Hmid +Hup]
% WL = 2*pi/3; Wmid = 2*pi/3; WR = 2*pi/3; 
% Hlow  = 2*pi/3 ; Hmid = 2*pi/3; Hup    = 2*pi/3;

WL = 1e-4; Wmid = 1.5e-4; WR = 1e-4; 
Hlow  = .75e-4; Hmid = .2e-4; Hup   = .75e-4;

Nx = [300, 350, 400, 450, 500]'; 
Ny = Nx;
NU_y = 0*Nx;
NU_x = NU_y;
NV_y = NU_y;
NV_x = NU_y;

GRID = 1;
PLOT = 0;

% Get grid for FVM 
for i = 1:length(Nx)
    [XU,YU,XV,YV,xp,yp,N] = NSGridHnoGhst(Nx(i),Ny(i),WL,Wmid,WR,Hlow,Hmid,Hup,GRID,PLOT);
    
    % Unpack Various N values for H-domain
    NxL = N(1); NxM = N(2); 
    NyL = N(4); NyM = N(5); 
    
    % COUNT DOES NOT INCLUDE ANY NODES (GHOST) OUTSIDE INJ. CHANNEL
    NU_y(i) = length(XU(NyL+2:NyM+1,NxL+1));
    NU_x(i) = length(XU(NyL+2,NxL+1:NxM+1));

    NV_y(i) = length(XV(NyL+1:NyM+1,NxL+2));
    NV_x(i) = length(XV(NyL+2,NxL+2:NxM+1));
end

total_U = NU_x.*NU_y;
total_V = NV_x.*NV_y;
table(Nx,Ny,NU_x,NU_y,total_U,NV_x,NV_y,total_V)