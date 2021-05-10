syms theta_1(t) theta_2(t) theta_3(t) L_1 L_2 L_3 m_1 m_2 m_3 g

x_1 = L_1*sin(theta_1);
y_1 = -L_1*cos(theta_1);
x_2 = x_1 + L_2*sin(theta_2);
y_2 = y_1 - L_2*cos(theta_2);
x_3 = x_2 + L_3*sin(theta_3);
y_3 = y_2 - L_3*cos(theta_3);

vx_1 = diff(x_1);
vy_1 = diff(y_1);
vx_2 = diff(x_2);
vy_2 = diff(y_2);
vx_3 = diff(x_3);
vy_3 = diff(y_3);

ax_1 = diff(vx_1);
ay_1 = diff(vy_1);
ax_2 = diff(vx_2);
ay_2 = diff(vy_2);
ax_3 = diff(vx_3);
ay_3 = diff(vy_3);

% define equation of motion
syms T_1 T_2 T_3

eqx_1 = m_1*ax_1(t) == -T_1*sin(theta_1(t)) + T_2*sin(theta_2(t))

eqy_1 = m_1*ay_1(t) == T_1*cos(theta_1(t)) - T_2*cos(theta_2(t)) - m_1*g

eqx_2 = m_2*ax_2(t) == -T_2*sin(theta_2(t)) + T_3*cos(theta_3(t))

eqy_2 = m_2*ay_2(t) == T_2*cos(theta_2(t)) - T_3*cos(theta_3(t))- m_2*g

eqx_3 = m_3*ax_3(t) == -T_3*sin(theta_3(t))

eqy_3 = m_3*ay_3(t) == T_3*cos(theta_3(t)) - m_3*g

Tension_T1_T2 = solve([eqx_1 eqy_1],[T_1 T_2]);

eqx_2_sub = subs(eqx_2, T_2, Tension_T1_T2.T_2);
eqy_2_sub = subs(eqy_2, T_2, Tension_T1_T2.T_2);

Tension_T3 = solve(eqx_2_sub, T_3);

eqRed_1 = subs(eqx_3, [T_3 T_2 T_1], [Tension_T3 Tension_T1_T2.T_2 Tension_T1_T2.T_1]);
eqRed_2 = subs(eqy_3, [T_3 T_2 T_1], [Tension_T3 Tension_T1_T2.T_2 Tension_T1_T2.T_1]);
eqRed_3 = subs(eqx_2, [T_3 T_2 T_1], [Tension_T3 Tension_T1_T2.T_2 Tension_T1_T2.T_1]);
eqRed_4 = subs(eqy_2, [T_3 T_2 T_1], [Tension_T3 Tension_T1_T2.T_2 Tension_T1_T2.T_1]);

%% solve the system

L_1 = 1;
L_2 = 1.5;
L_3 = 2;
m_1 = 2;
m_2 = 1;
m_3 = 0.5;
g = 9.8;

eqn_1 = subs(eqRed_1)

eqn_2 = subs(eqRed_2)

eqn_3 = subs(eqRed_3)

eqn_4 = subs(eqRed_4)

[V,S] = odeToVectorField(eqn_1,eqn_2, eqn_3, eqn_4);

M = matlabFunction(V,'vars',{'t','Y'});

initCond = [pi/4 0 pi/6 0 pi/8 0];
%opts = odeset('RelTol',1e-2,'AbsTol',1e-4);
sols = ode45(M,[0 10],initCond);   %opts

plot(sols.x, sols.y)
legend('\theta_3','d\theta_3/dt', '\theta_2','d\theta_2/dt','\theta_1','d\theta_1/dt')
title('Solutions of State Variables')
xlabel('Time (s)')
ylabel('Solutions (rad or rad/s)')

%% animation of single double pendulum

x_1 = @(t) L_1*sin(deval(sols,t,5));
y_1 = @(t) -L_1*cos(deval(sols,t,5));
x_2 = @(t) L_1*sin(deval(sols,t,5))+L_2*sin(deval(sols,t,3));
y_2 = @(t) -L_1*cos(deval(sols,t,5))-L_2*cos(deval(sols,t,3));
x_3 = @(t) L_1*sin(deval(sols,t,5))+L_2*sin(deval(sols,t,3))+L_3*sin(deval(sols,t,1));
y_3 = @(t) -L_1*cos(deval(sols,t,5))-L_2*cos(deval(sols,t,3))-L_3*cos(deval(sols,t,1));

fanimator(@(t) plot(x_1(t),y_1(t),'ro','MarkerSize',m_1*10,'MarkerFaceColor','r'));
axis equal;

hold on;
fanimator(@(t) plot([0 x_1(t)],[0 y_1(t)],'r-'));
fanimator(@(t) plot(x_2(t),y_2(t),'go','MarkerSize',m_2*10,'MarkerFaceColor','g'));
fanimator(@(t) plot([x_1(t) x_2(t)],[y_1(t) y_2(t)],'g-'));
fanimator(@(t) plot(x_3(t),y_3(t),'bo','MarkerSize',m_3*10,'MarkerFaceColor','g'));
fanimator(@(t) plot([x_2(t) x_3(t)],[y_2(t) y_3(t)],'b-'));

fanimator(@(t) text(-0.3,0.3,"Timer: "+num2str(t,2)));
hold off;

playAnimation