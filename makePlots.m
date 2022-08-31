function [] = makePlots(U,V,P,Ue,Ve,Pe,XU,YU,XV,YV,XP,YP,t)
% Compute Error in u-eqt
xAbsErr = abs(Ue - U);

% Compute Error in v-eqt
yAbsErr = abs(Ve - V);

% Compute Error in P
pAbsErr = abs(Pe - P);

% Plot the U solution
subplot(3,3,1)
surf(XU,YU,Ue, 'EdgeColor', 'none')
xlabel('x')
ylabel('y')
zlabel('U')
title('Exact Solution')
% axis([XU(1,1), XU(end, end), YU(1,1), YU(end,end), -1, 1])

subplot(3,3,2)
surf(XU,YU,U, 'EdgeColor', 'none')
xlabel('x')
ylabel('y')
zlabel('U')
title('Numerical Solution')
% axis([XU(1,1), XU(end, end), YU(1,1), YU(end,end), -1, 1])

subplot(3,3,3)
surf(XU,YU,xAbsErr, 'EdgeColor', 'none')
xlabel('x')
ylabel('y')
zlabel('Error')
title('Error')


% Plot the V Solution
subplot(3,3,4)
surf(XV,YV,Ve, 'EdgeColor', 'none')
xlabel('x')
ylabel('y')
zlabel('V')
title('Exact Solution')
% axis([XV(1,1), XV(end, end), YV(1,1), YV(end,end), -1, 1])

subplot(3,3,5)
surf(XV,YV,V, 'EdgeColor', 'none')
xlabel('x')
ylabel('y')
zlabel('V')
title('Numerical Solution')
% axis([XV(1,1), XV(end, end), YV(1,1), YV(end,end), -1, 1])

subplot(3,3,6)
surf(XV,YV,yAbsErr, 'EdgeColor', 'none')
xlabel('x')
ylabel('y')
zlabel('Error')
title('Error')
        
% Plot the P Solution
subplot(3,3,7)
% pcolor(XP,YP,Pe)
surf(XP,YP,Pe)
shading flat
xlabel('x')
ylabel('y')
zlabel('P')
title('Exact Solution')
% axis([XP(1,1), XP(end, end), YP(1,1), YP(end,end), -1, 1])

subplot(3,3,8)
% pcolor(XP,YP,P)
surf(XP,YP,P)
shading flat
xlabel('x')
ylabel('y')
zlabel('P')
title('Numerical Solution')
% axis([XP(1,1), XP(end, end), YP(1,1), YP(end,end), -1, 1])

subplot(3,3,9)
% pcolor(XP,YP,pAbsErr)
surf(XP,YP,pAbsErr)
shading flat
xlabel('x')
ylabel('y')
zlabel('Error')
title('Error')
sgtitle(['$\vec{u}(x,y,$', num2str(t), '$)\;\; \& \;\; p(x,y,$', num2str(t),...
    '$)$'], 'FontSize', 18, 'Interpreter', 'latex')
