function[U,V,P] = Correction_Ghst(U,V,P0,PHI,N,dt,xp,yp)

% Unpack Various N values for H-domain
NxL = N(1); NxM = N(2); Nx = N(3);
NyL = N(4); NyM = N(5); Ny = N(6);

% update Pressure
P = P0 + PHI;
% Zero out interior of H
P(1:NyL,NxL+1:NxM) = 0;
P(NyM+1:Ny,NxL+1:NxM) = 0;

% update U - Wash Channel
i = 2 : Ny+1;
j = 2 : NxL;
U(i,j) = U(i,j) - dt * ( PHI(i-1,j) - PHI(i-1,j-1) )./( xp(j+1)-xp(j) );

% update U - Injury Channel
i = NyL+2 : NyM;
j = NxL+1 : NxM+1;
U(i,j) = U(i,j) - dt * ( PHI(i-1,j) - PHI(i-1,j-1) )./( xp(j+1)-xp(j) );

% update U - Blood Channel
i = 2 : Ny+1;
j = NxM+2 : Nx;
U(i,j) = U(i,j) - dt * ( PHI(i-1,j) - PHI(i-1,j-1) )./( xp(j+1)-xp(j) );

% update V - Wash Channel
i = 2 : Ny;
j = 2 : NxL+1;
V(i,j) = V(i,j) - dt * ( PHI(i,j-1) - PHI(i-1,j-1) )./( yp(i+1)-yp(i) )';

% update V - Injury Channel
i = NyL+2 : NyM;
j = NxL+2 : NxM+1;
V(i,j) = V(i,j) - dt * ( PHI(i,j-1) - PHI(i-1,j-1) )./( yp(i+1)-yp(i) )';


% update V - Blood Channel
i = 2 : Ny;
j = NxM+2 : Nx+1;
V(i,j) = V(i,j) - dt * ( PHI(i,j-1) - PHI(i-1,j-1) )./( yp(i+1)-yp(i) )';