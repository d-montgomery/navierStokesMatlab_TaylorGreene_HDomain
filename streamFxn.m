function PSI = streamFxn(U,V,xp,yp,Nx,Ny)
% Inputs:
% U and V are the velocites with ghost cells 
% xp and yp are the MAC grid locations of pressure
% Nx and Ny are the number of spatial cells

PSI = zeros(Ny, Nx);

% Compute first row
for I = 2:Ny
    % Compensate for Ghost Points
    i = I + 1;
    
    dy = yp(i) - yp(i-1);
    
    uw =  U(i-1,1);
    ue = U(i-1,2);
    unw = U(i,1);
    une = U(i,2);
    
    un = 0.5 * ( unw + une);
    up = 0.5 * ( uw + ue );
    
    PSI(I,1) = PSI(I-1,1) +  0.5 * dy *( un + up );
end


% Compute all other columns
I = 1:Ny;
for J = 2:Nx
    j = J+1;
    dx = xp(j) - xp(j-1);
    
    vn = V(I+1 , j-1);
    vs = V(I , j-1);
    vne = V(I+1 , j);
    vse = V(I  , j);
    
    ve = 0.5 * ( vne + vse);
    vp = 0.5 * ( vn + vs );
    
    PSI(I,J) = PSI(I,J-1) - 0.5 * dx .* ( ve + vp );
end

    
    
    
    
    
    
