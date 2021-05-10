syms theta_1(t) L_1 m_1 g

x_1 = L_1*sin(theta_1);
y_1 = -L_1*cos(theta_1);

vx_1 = diff(x_1);
vy_1 = diff(y_1);

ax_1 = diff(vx_1);
ay_1 = diff(vy_1);

% define equation of motion
syms T_1 

eqx_1 = m_1*ax_1(t) == -T_1*sin(theta_1(t))

eqy_1 = m_1*ay_1(t) == T_1*cos(theta_1(t)) - m_1*g

% eqx_2 = m_2*ax_2(t) == -T_2*sin(theta_2(t))
% 
% eqy_2 = m_2*ay_2(t) == T_2*cos(theta_2(t)) - m_2*g

Tension = solve(eqy_1, T_1);

eqRed_1 = subs(eqx_1, T_1, Tension);
eqRed_2 = subs(eqy_1, T_1, Tension);

%% solve the system

L_1 = 1;
%L_2 = 1.5;
m_1 = 2;
%m_2 = 1;
g = 9.8;
eqn_1 = subs(eqRed_1)

eqn_2 = subs(eqRed_2)

[V,S] = odeToVectorField(eqn_1,eqn_2);

M = matlabFunction(V,'vars',{'t','Y'});

initCond = [pi/4 0];
sols = ode45(M,[0 10],initCond);

plot(sols.x,sols.y)
legend('\theta_1','d\theta_1/dt')
title('Solutions of State Variables')
xlabel('Time (s)')
ylabel('Solutions (rad or rad/s)')

%% animation of single double pendulum

x_1 = @(t) L_1*sin(deval(sols,t,1));
y_1 = @(t) -L_1*cos(deval(sols,t,1));
%x_2 = @(t) L_1*sin(deval(sols,t,3))+L_2*sin(deval(sols,t,1));
%y_2 = @(t) -L_1*cos(deval(sols,t,3))-L_2*cos(deval(sols,t,1));

fanimator(@(t) plot(x_1(t),y_1(t),'ro','MarkerSize',m_1*5,'MarkerFaceColor','r'));
axis equal;

hold on;
fanimator(@(t) plot([0 x_1(t)],[0 y_1(t)],'r-'));
%fanimator(@(t) plot(x_2(t),y_2(t),'go','MarkerSize',m_2*10,'MarkerFaceColor','g'));
%fanimator(@(t) plot([x_1(t) x_2(t)],[y_1(t) y_2(t)],'g-'));

fanimator(@(t) text(-0.3,0.3,"Timer: "+num2str(t,2)));
hold off;


playAnimation