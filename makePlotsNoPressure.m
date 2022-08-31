function [] = makePlotsNoPressure(U,V,Ue,Ve,XU,YU,XV,YV,t)
% Compute Error in x-eqt

% Up = reshape(U, Nx+1,Ny+2)';
xAbsErr = abs(Ue - U);

% Compute Error in y-eqt
% Vp = reshape(V, Nx+2,Ny+1)';
yAbsErr = abs(Ve - V);


figure(1)
% Plot the U solution
subplot(2,3,1)
surf(XU,YU,Ue, 'EdgeColor', 'none')
xlabel('x')
ylabel('y')
zlabel('U')
title('Exact Solution')
axis([XU(1,1), XU(end, end), YU(1,1), YU(end,end), -1, 1])

subplot(2,3,2)
surf(XU,YU,U, 'EdgeColor', 'none')
xlabel('x')
ylabel('y')
zlabel('U')
title('Numerical Solution')
axis([XU(1,1), XU(end, end), YU(1,1), YU(end,end), -1, 1])

subplot(2,3,3)
surf(XU,YU,xAbsErr, 'EdgeColor', 'none')
xlabel('x')
ylabel('y')
zlabel('Error')
title('Error')


% Plot the V Solution
subplot(2,3,4)
surf(XV,YV,Ve, 'EdgeColor', 'none')
xlabel('x')
ylabel('y')
zlabel('V')
title('Exact Solution')
axis([XV(1,1), XV(end, end), YV(1,1), YV(end,end), -1, 1])

subplot(2,3,5)
surf(XV,YV,V, 'EdgeColor', 'none')
xlabel('x')
ylabel('y')
zlabel('V')
title('Numerical Solution')
axis([XV(1,1), XV(end, end), YV(1,1), YV(end,end), -1, 1])

subplot(2,3,6)
surf(XV,YV,yAbsErr, 'EdgeColor', 'none')
xlabel('x')
ylabel('y')
zlabel('Error')
title('Error')
        
sgtitle(['$\vec{u}(x,y,$', num2str(t), '$)$'], 'FontSize', 18, 'Interpreter', 'latex')
