syms m0 m1 m2 m3 m4 m5 m6 t
syms x(t) gemma(t) I0xx I0yy I0zz I1xx I1yy I1zz I2xx I2yy I2zz I3xx I3yy I3zz
syms q11(t) q12(t) q13(t) q21(t) q22(t) q23(t)%[rad]
% syms q11 q12 q13  q21 q22 q23 
syms dx(t)  dgemma(t) dq11(t) dq12(t) dq13(t)  dq21(t) dq22(t) dq23(t) 
syms a1 a2 a3 a4  a5 a6 a34 bd a33
syms I b
   
x=x(t);
gemma=gemma(t);
q11 = q11(t); 
q12 = q12(t); q13 = q13(t); q21 = q21(t); q22 = q22(t); q23 = q23(t);
dx=dx(t);dgemma=dgemma(t);dq11 = dq11(t); dq12 = dq12(t); dq13 = dq13(t); dq21 = dq21(t); dq22 = dq22(t); dq23 = dq23(t);
assume([x, gemma,q11,q12,q13,q21,q22,q23,dx,dgemma,dq11,dq12,dq13,dq21,dq22,dq23], 'real')
% x=0;gemma=0;
% q11 = pi/4;q12 = pi/4;q13 = -pi/4;
% q21 = pi/4;q22 = pi/4;q23 = -pi/4;%初始位形
% 鍩虹鏁版嵁 
a1=a4;
a2=a5;
a3=a6;
m0=10.18;m1=0.5;m2=0.6;m3=0.6;
m4=m1;m5=m2;m6=m3;
a1=0.6;a2=0.5;a3=0.4;
a4=a1;a5=a2;a6=a3;
b=0.258;
% x1c=0;       y1c=-a1/2;        z1c=0;
% x2c=a2/2;       y2c=0;        z2c=0;
% x3c=0;       y3c=0;        z3c=a3/2;
% I0xx=10;I0yy=10;I0zz=10;
I0xx=10.4;I0yy=10.4;I0zz=10.4;
I1xx=2*10^(-4);I1yy=2*10^(-4);I1zz=2*10^(-4);
I2xx=3.5*10^(-4);I2yy=3.5*10^(-4);I2zz=2*10^(-4);
I3xx=3.5*10^(-4);I3yy=3.7*10^(-4);I3zz=2*10^(-4);
% % I1xx=0.0048;I1yy=0.0036;I1zz=0.08;
% I2xx=1.2552;I2yy=1.254;I2zz=0.85;
% I3xx=0.0048;I3yy=0.0036;I3zz=0.03;      

I0=[I0xx 0 0;
    0 I0yy 0;
    0 0 I0zz];
I1=[I1xx 0 0;
    0 I1yy 0;
    0 0 I1zz];
I2=[I2xx 0 0;
    0 I2yy 0;
    0 0 I2zz];
I3=[I3xx 0 0;
    0 I3yy 0;
    0 0 I3zz];
I4=I1;
I5=I2;
I6=I3; 
It=[2.08*10^(-4),0,0;
    0,2.08*10^(-4),0;
    0,0,3.33*10^(-4)];
 I=[I1 I2 I3 I4 I5 I6 It ];
 m=[m1 m2 m3 m4 m5 m6 ];
 mt=0.2;

q=[x;0;0;0;0;gemma;q11;q12;q13;q21;q22;q23;];
dq=[dx;0;0;0;0;dgemma;dq11;dq12;dq13;dq21;dq22;dq23;];
thata=[q11;q12;q13;q21;q22;q23;];
dtheta=[dq11;dq12;dq13;dq21;dq22;dq23;];
m_tot=(m0+m1+m2+m3+m4+m5+m6);
m_tott=(m0+m1+m2+m3+m4+m5+m6+mt);




%% 运动学
%惯性坐标系-基座 i-b
Rb_I=[1,0,0; 
      0,cos(gemma),-sin(gemma);
      0,sin(gemma),cos(gemma);];

Tb_I=[1,0,0,x; 
      0,cos(gemma),-sin(gemma),2;
       0,sin(gemma),cos(gemma),2;
       0,0,0,1;];
%左臂
%基座-左边连接点 b-b1
Rb1_b=[1,0,0; 
      0,cos(gemma),-sin(gemma);
      0,sin(gemma),cos(gemma);];
Tb1_b=[1,0,0,x; 
      0,cos(gemma),-sin(gemma),1;
      0,sin(gemma),cos(gemma),2;
      0,0,0,1];
% Tb1_I=Tb_I*Tb1_b;


%左臂 基座--1
A1_b1=[cos(q11) -sin(q11) 0 a1*cos(q11);
       sin(q11) cos(q11)  0 a1*sin(q11);
       0        0       1 0;
       0        0       0 1];


% A1_0=simplify(A1_0);%b1-1
%1-2
A2_1=[cos(q12) -sin(q12) 0 a2*cos(q12);
       sin(q12) cos(q12) 0 a2*sin(q12);
       0        0       1 0;
       0        0       0 1];





% 2-3 
A3_2=[cos(q13) -sin(q13) 0 a3*cos(q13);
       sin(q13) cos(q13)  0 a3*sin(q13);
       0        0       1 0;
       0        0       0 1];

% 3-e1 关节3-end effctor
Ae1_3=[1 0 0 0;
      0 1 0 0;
      0 0 1 b/2;
      0 0 0 1];



% 节点齐次变换
A1_0=Tb_I*Tb1_b*A1_b1;
A2_0=A1_0*A2_1;
A3_0=A2_0*A3_2;
Ae1_0=A3_0*Ae1_3;


% 连杆末端坐标
rb=[x;2;2;]
p0=[0 0 0]';
p1=A1_0(1:3,4);
p2=A2_0(1:3,4);
p3=A3_0(1:3,4);
pe1=Ae1_0(1:3,4);

% p=[p1 p2 p3 p4 p5 p6];
% 质心齐次变换
A1c_1=[1 0 0 a1/2;
       0 1 0 0;
       0 0 1 0;   
       0 0 0 1];
A2c_2=[1 0 0 a2/2;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];
A3c_3=[1 0 0 a3/2;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];



A1c_0=A1_0*A1c_1;
A2c_0=A2_0*A2c_2;
A3c_0=A3_0*A3c_3;



% 质心变换矩阵
R1=A1c_0(1:3,1:3);
R2=A2c_0(1:3,1:3);
R3=A3c_0(1:3,1:3);


% 节点旋转矩阵
R1_0=A1_0(1:3,1:3);
R2_0=A2_0(1:3,1:3);
R3_0=A3_0(1:3,1:3);

%右臂
%基座-右臂连接点b2
Rb2_b=[0,1,0;
       -1,0,0;
       0,0,-1;];
Tb2_b=[1,0,0,x; 
      0,cos(gemma),-sin(gemma),3;
      0,sin(gemma),cos(gemma),2;
      0,0,0,1];
% Tb2_I=Tb_I*Tb2_b;

%b2-4
A4_b2=[cos(q21) -sin(q21) 0 a4*cos(q21);
       sin(q21) cos(q21)  0 a4*sin(q21);
       0        0       1 0;
       0        0       0 1];

%4-5
A5_4=[cos(q22) -sin(q22) 0 a5*cos(q22);
       sin(q22) cos(q22)  0 a5*cos(q22);
       0        0       1 0;
       0        0       0 1];



% 5-6 
A6_5=[cos(q23) -sin(q23) 0 a6*cos(q23);
       sin(q23) cos(q23)  0 a6*cos(q23);
       0        0       1 0;
       0        0       0 1];

% 6-e2 关节6-end effctor
Ae2_6=[1 0 0 0;
      0 1 0 0;
      0 0 1 b/2;
      0 0 0 1];
% 节点齐次变换
A4_0=Tb_I*Tb2_b*A4_b2;
A5_0=A4_0*A5_4;
A6_0=A5_0*A6_5;
Ae2_0=A6_0*Ae2_6;


% 连杆末端坐标
rb=[x;2;2;]
p0=[0 0 0]';
r0=Tb_I(1:3,4);%基座
p4=A4_0(1:3,4);
p5=A5_0(1:3,4);
p6=A6_0(1:3,4);
pe2=Ae2_0(1:3,4);

% p=[p1 p2 p3 p4 p5 p6];
% 质心齐次变换
A4c_4=[1 0 0 a1/2;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];
A5c_5=[1 0 0 a2/2;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];
A6c_6=[1 0 0 a3/2;
       0 1 0 0;
       0 0 1 0;
       0 0 0 1];

A4c_0=A4_0*A4c_4;
A5c_0=A5_0*A5c_5;
A6c_0=A6_0*A6c_6;



% 质心旋转矩阵
R4=A4c_0(1:3,1:3);
R5=A5c_0(1:3,1:3);
R6=A6c_0(1:3,1:3);

% 节点旋转矩阵
R4_0=A4_0(1:3,1:3);
R5_0=A5_0(1:3,1:3);
R6_0=A6_0(1:3,1:3);



p1c=A1c_0(1:3,4);
p2c=A2c_0(1:3,4);
p3c=A3c_0(1:3,4);
p4c=A4c_0(1:3,4);
p5c=A5c_0(1:3,4);
p6c=A6c_0(1:3,4);
pc=[p1c p2c p3c p4c p5c p6c];%连杆质心
% R=[R1 R2 R3 R4 R5 R6 ];%质心旋转矩阵

%模块
Ae1_b=[1 0 0 0;
      0 1 0 0;
      0 0 1 b/2;
      0 0 0 1];
%模块质心在惯性坐标系下
Ab_0=A3_0*Ae1_3*Ae1_b;
pb=Ab_0(1:3,4);
R7=Ab_0(1:3,1:3);
R=[R1 R2 R3 R4 R5 R6 R7];%质心旋转矩阵

%% 动力学?   
z0=[0 0 1]';
z1=R1_0*z0;
z2=R2_0*z0;
z3=R3_0*z0;
z4=R4_0*z0;
z5=R5_0*z0;
z6=R6_0*z0;
z7=R7*z0;

O=[0 0 0]';
JT1=[myCrossProduct(z1,p1c-Tb1_b(1:3,4)),zeros(3,2)];
JT2=[myCrossProduct(z1,p1c-Tb1_b(1:3,4)),myCrossProduct(z2,p2c-p2),zeros(3,1)];
JT3=[myCrossProduct(z1,p1c-Tb1_b(1:3,4)),myCrossProduct(z2,p2c-p2),myCrossProduct(z3,p3c-p3)];
JT4=[myCrossProduct(z4,p4c-Tb2_b(1:3,4)),zeros(3,2)];
JT5=[myCrossProduct(z4,p4c-Tb2_b(1:3,4)),myCrossProduct(z5,p5c-p5),zeros(3,1)];
JT6=[myCrossProduct(z4,p4c-Tb2_b(1:3,4)),myCrossProduct(z5,p5c-p5),myCrossProduct(z6,p6c-p6)];
JT7=[myCrossProduct(z5,pb-Tb2_b(1:3,4)),myCrossProduct(z6,pb-p5),myCrossProduct(z7,pb-pe1)];
%模块算在1臂上
JT=[JT1 JT2 JT3 JT4 JT5 JT6 JT7];

JR1=[z1,zeros(3,2)];
JR2=[z1,z2,zeros(3,1)];
JR3=[z1,z2,z3];
JR4=[z4,zeros(3,2)];
JR5=[z4,z5,zeros(3,1)];
JR6=[z4,z5,z6];
JR7=[z5,z6,z7];
JR=[JR1 JR2 JR3 JR4 JR5 JR6 JR7 ];

%路径总时间 tf=10s
tf=10;%(s) 
st=t/tf;dst=1/tf;ddst=0;
%%  模块 关于s的贝塞尔曲线线速度插值
% 定义三维贝塞尔曲线的控制点
%第一条
P0 = [0, 0, 0];%初始
P1 = [0.1, 0.2, 0.3];
P2 = [0.2, -0.2, 0.7];
P3 = [0.5, 0, 1];
% 经过(0.15,0,0.5)

% 定义参数 s 的范围
st = linspace(0, 1, 100);

% 计算贝塞尔曲线上的点
[x_m, y_m, z_m, dx, dy, dz, ddx, ddy, ddz] = bezier_volcitycurve(P0, P1, P2, P3);
dm_t=[dx;dy;dz]*dst;%模块x、y、z方向速度 关于t
ddm_t=[ddx;ddy;ddz]*dst^2;%模块x、y、z方向加速度 关于t
% 绘制三维贝塞尔曲线

figure;
plot3(x_m, y_m, z_m, 'b-', 'LineWidth', 2);
hold on;
plot3([P0(1), P1(1), P2(1), P3(1)], [P0(2), P1(2), P2(2), P3(2)], [P0(3), P1(3), P2(3), P3(3)], 'ro--');
title('三维贝塞尔曲线');
xlabel('X');
ylabel('Y');
zlabel('Z');
legend('贝塞尔曲线', '控制点');
grid on;

% 绘制速度和加速度
figure;
subplot(2,1,1);
plot(10*st, dm_t(1,:), 'r-', 10*st, dm_t(2,:), 'g-', 10*st, dm_t(3,:), 'b-');
title('速度');
legend('dx/ds', 'dy/ds', 'dz/ds');
grid on;

subplot(2,1,2);
plot(10*st, ddm_t(1,:), 'r-', 10*st, ddm_t(2,:), 'g-', 10*st, ddm_t(3,:), 'b-');  
title('加速度');
legend('d^2x/ds^2', 'd^2y/ds^2', 'd^2z/ds^2');
grid on;

%%  模块 关于s的贝塞尔曲线角速度插值
% 定义时间点和对应的三维角速度数据（以弧度为单位）
% t = [0, 5, 7, 10]; % 时间点
omega_rad = [0, 0, 0;0.2,0.12,0.06;0.4,0.16,0.08; 0.5, 0.2, 0.1; ]; % 三维角速度数据（弧度）

% 定义贝塞尔曲线的控制点
P4 = omega_rad(1, :);
P5 = omega_rad(2, :);
P6 = omega_rad(3, :);
P7 = omega_rad(4, :);
% 定义参数 s 的范围
% st = linspace(0, 1, 100);
% 计算贝塞尔曲线上的点
[omega_x, omega_y, omega_z] = bezier_angularcurve(P4, P5, P6, P7);
% 绘制三维角速度插值结果
figure;
plot3(omega_x, omega_y, omega_z, 'b-', 'LineWidth', 2);
hold on;
plot3(omega_rad(:,1), omega_rad(:,2), omega_rad(:,3), 'ro');
title('三维角速度（弧度）贝塞尔曲线插值');
xlabel('Omega X (rad/s)');
ylabel('Omega Y (rad/s)');
zlabel('Omega Z (rad/s)');
legend('插值曲线', '原始数据点');
grid on;

%% 由模块位置得到end-effector1位置
[v_ee1, w_ee1,v_ee2, w_ee2,x_ee1,y_ee1,z_ee1,x_ee2,y_ee2,z_ee2] = cal_end_effector_velocity(x_m, y_m, z_m, dm_t, omega_x, omega_y, omega_z)
T_ee1_m=[1 0 0 0;
         0 1 0 0;
         0 0 1 0.258/2;
         0 0 0 1]; 
v_ee1x=v_ee1(1,:);
v_ee1y=v_ee1(2,:);
v_ee1z=v_ee1(3,:);
w_ee1x=w_ee1(1,:);
w_ee1y=w_ee1(2,:);
w_ee1z=w_ee1(3,:);

%% 由模块位置得到end-effector2位置
T_ee2_m=[1 0 0 0;
         0 1 0 0;
         0 0 1 -0.258/2;
         0 0 0 1];
v_ee2x=v_ee2(1,:);
v_ee2y=v_ee2(2,:);
v_ee2z=v_ee2(3,:);
w_ee2x=w_ee2(1,:);
w_ee2y=w_ee2(2,:);
w_ee2z=w_ee2(3,:);
%%整合末端点速度
dot_vee=[v_ee1,;w_ee1;v_ee2;w_ee2];


%% 逆运动学
%计算初始值
rb(:,1)=[0;2;2];%基座初始位置
r_e1(:,1)=[x_ee1(1);y_ee1(1);z_ee1(1)];%1臂末端
r_e2(:,1)=[x_ee2(1);y_ee2(1);z_ee2(1)];%2臂末端
J_b1=[eye(3) -skew(r_e1(:,1)-rb(:,1));
         zeros(3,3) eye(3)]; 
J_b2=[eye(3) -skew(r_e2(:,1)-rb(:,1));
         zeros(3,3) eye(3)];
r0g=(m1*(p1c-Tb_I(1:3,4))+m2*(p2c-Tb_I(1:3,4))+m3*(p3c-Tb_I(1:3,4))+m4*(p4c-Tb_I(1:3,4))+m5*(p5c-Tb_I(1:3,4))+m6*(p6c-Tb_I(1:3,4))+mt*(pb-Tb_I(1:3,4)))/m_tott;
r0g=vpa(r0g,5);
r0g(:,1)=[0.142;0.791;0.523];
Hs(:,:,1)=[49.447, -3.9842, -2.9335;
         -3.9838,22.126,-16.185;
          -2.9335,-16.185,39.563];
H_q1(:,:,1)=[-1.8393, -0.097963, -0.12505;
            -7.0233,0.15261,-0.12505;
            6.7213, -0.012537,  0.19256];
JT_w1(:,:,1)=[-4.7233, 0.043431, -0.084853;
              1.2233, 0.056569,  0.084853;
              0,0,0];
H_q2(:,:,1)=[-1.5944, 0, -0.12505;
             -6.6052, -0.44213,-0.12505;
             15.047,1.0414, 0.31393];
JT_w2(:,:,1)=[-4.4819,-0.3,-0.084853;
              1.0819, 0, 0.084853;
              0,0,0];
z10=[0;0;1];z20=[0;0;1];z30=[0;0;1];
z40=[0;0;1];z50=[0;0;1];z60=[0;0;1];
bmv11=-skew(r0g(:,1))/(Hs(:,:,1))*H_q1(:,:,1)-JT_w1(:,:,1);
bmv12=-skew(r0g(:,1))/(Hs(:,:,1))*H_q2(:,:,1)-JT_w2(:,:,1);
bmw11=-inv(Hs(:,:,1))*H_q1(:,:,1);
bmw12=-inv(Hs(:,:,1))*H_q2(:,:,1);
J_bmv=[bmv11  bmv12];%此处应用数值计算
J_bmv=vpa(J_bmv,5);
J_bmw=[bmw11  bmw12];%此处应用数值计算
J_bmw=vpa(J_bmv,5);
Jg(:,:,1)=([J_b1;J_b2]*[J_bmv;J_bmw]+[myCrossProduct(z10,[-0.7828;0.2828;0.1290]),myCrossProduct(z20,[0.1414;-3.1414;-3.8710]),myCrossProduct(z30,[0.42;-3.42;-3.87]),zeros(3,3);
                                       z10,z20,z30,zeros(3,3);
                                       zeros(3,3),myCrossProduct(z40,[0.78;0.28;-0.129]),myCrossProduct(z50,[0.28;0.28;-0.129]),myCrossProduct(z60,[0;0;-0.13]);
                                       zeros(3,3),z40,z50,z60]);

Jg(:,:,1)=vpa(Jg(:,:,1),5);
weiinvJg(:,:,1)=Jg(:,:,1)'/(Jg(:,:,1)*Jg(:,:,1)');
theta(:,1)=zeros(6,1);
x_base(1)=0;y_base(1)=2;z_base(1)=2;
v_basex(1)=0;v_basey(1)=0;v_basez(1)=0;
w_basex(1)=0;w_basey(1)=0;w_basez(1)=0;
vw_basex(1)=0;vw_basey(1)=0;vw_basez(1)=0;
pe1(:,1)=[0.0026;0.0061;-0.12];
pe2(:,1)=[0.0035;0.0056;0.14];

dt=1;

for i=1:dt:100
    xe_current = [x_ee1(i);y_ee1(i);z_ee1(i);x_ee2(i);y_ee2(i);z_ee2(i)];
    p1 =[(3*cos(theta(1,i)))/5 + 2*x_base(i);
        cos(w_basex(i)) - 2*sin(w_basex(i)) + (3*sin(theta(1,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2))/5 + 2;
        2*cos(w_basex(i)) + sin(w_basex(i)) + (6*cos(w_basex(i))*sin(w_basex(i))*sin(theta(1,i)))/5 + 2];
    p2 =[(3*cos(theta(1,i)))/5 + 2*x_base(i) - (sin(theta(1,i))*sin(theta(2,i)))/2 + (cos(theta(1,i))*cos(theta(1,i)))/2;
         cos(w_basex(i)) - 2*sin(w_basex(i)) + (3*sin(theta(1,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2))/5 + (cos(theta(1,i))*sin(theta(2,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2))/2 + (cos(theta(2,i))*sin(theta(1,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2))/2 + 2;
         2*cos(w_basex(i)) + sin(w_basex(i)) + (6*cos(w_basex(i))*sin(w_basex(i))*sin(theta(1,i)))/5 + cos(w_basex(i))*cos(theta(1,i))*sin(w_basex(i))*sin(theta(2,i)) + cos(w_basex(i))*cos(theta(2,i))*sin(w_basex(i))*sin(theta(1,i)) + 2];
   p3 =[(3*cos(theta(1,i)))/5 + 2*x_base(i) - (2*cos(theta(3,i))*(sin(theta(1,i))*sin(theta(2,i)) - cos(theta(1,i))*cos(theta(2,i))))/5 - (sin(theta(1,i))*sin(theta(2,i)))/2 - (2*sin(theta(3,i))*(cos(theta(1,i))*sin(theta(2,i)) + cos(theta(2,i))*sin(theta(1,i))))/5 + (cos(theta(1,i))*cos(theta(2,i)))/2;
        cos(w_basex(i)) - 2*sin(w_basex(i)) + (2*cos(theta(3,i))*(cos(theta(1,i))*sin(theta(2,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2) + cos(theta(2,i))*sin(theta(1,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2)))/5 - (2*sin(theta(3,i))*(sin(theta(1,i))*sin(theta(2,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2) - cos(theta(1,i))*cos(theta(2,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2)))/5 + (3*sin(theta(1,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2))/5 + (cos(theta(1,i))*sin(theta(2,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2))/2 + (cos(theta(2,i))*sin(theta(1,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2))/2 + 2;
        2*cos(w_basex(i)) + sin(w_basex(i)) + (2*cos(theta(3,i))*(2*cos(w_basex(i))*cos(theta(1,i))*sin(w_basex(i))*sin(theta(2,i)) + 2*cos(w_basex(i))*cos(theta(2,i))*sin(w_basex(i))*sin(theta(1,i))))/5 + (2*sin(theta(3,i))*(2*cos(w_basex(i))*cos(theta(1,i))*cos(theta(2,i))*sin(w_basex(i)) - 2*cos(w_basex(i))*sin(w_basex(i))*sin(theta(1,i))*sin(theta(2,i))))/5 + (6*cos(w_basex(i))*sin(w_basex(i))*sin(theta(1,i)))/5 + cos(w_basex(i))*cos(theta(1,i))*sin(w_basex(i))*sin(theta(2,i)) + cos(w_basex(i))*cos(theta(2,i))*sin(w_basex(i))*sin(theta(1,i)) + 2];
    p4 =[(3*cos(theta(4,i)))/5 + 2*x_base(i);
        3*cos(w_basex(i)) - 2*sin(w_basex(i)) + (3*sin(theta(4,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2))/5 + 2;
        2*cos(w_basex(i)) + 3*sin(w_basex(i)) + (6*cos(w_basex(i))*sin(w_basex(i))*sin(theta(4,i)))/5 + 2];
 

 
p5 =[(3*cos(theta(4,i)))/5 + 2*x_base(i) - (cos(theta(5,i))*sin(theta(4,i)))/2 + (cos(theta(4,i))*cos(theta(5,i)))/2;
      3*cos(w_basex(i)) - 2*sin(w_basex(i)) + (3*sin(theta(4,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2))/5 + (cos(theta(5,i))*sin(theta(4,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2))/2 + (cos(theta(4,i))*cos(theta(5,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2))/2 + 2;
      2*cos(w_basex(i)) + 3*sin(w_basex(i)) + (6*cos(w_basex(i))*sin(w_basex(i))*sin(theta(4,i)))/5 + cos(w_basex(i))*cos(theta(4,i))*cos(theta(5,i))*sin(w_basex(i)) + cos(w_basex(i))*cos(theta(5,i))*sin(w_basex(i))*sin(theta(4,i)) + 2];
 

p6 =[(3*cos(theta(4,i)))/5 + 2*x_base(i) - (cos(theta(5,i))*sin(theta(4,i)))/2 - (2*cos(theta(6,i))*(cos(theta(4,i))*sin(theta(5,i)) + cos(theta(5,i))*sin(theta(4,i))))/5 - (2*cos(theta(6,i))*(sin(theta(4,i))*sin(theta(5,i)) - cos(theta(4,i))*cos(theta(5,i))))/5 + (cos(theta(4,i))*cos(theta(5,i)))/2;
      3*cos(w_basex(i)) - 2*sin(w_basex(i)) + (2*cos(theta(6,i))*(cos(theta(4,i))*sin(theta(5,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2) + cos(theta(5,i))*sin(theta(4,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2)))/5 - (2*cos(theta(6,i))*(sin(theta(4,i))*sin(theta(5,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2) - cos(theta(4,i))*cos(theta(5,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2)))/5 + (3*sin(theta(4,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2))/5 + (cos(theta(5,i))*sin(theta(4,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2))/2 + (cos(theta(4,i))*cos(theta(5,i))*(cos(w_basex(i))^2 - sin(w_basex(i))^2))/2 + 2;
      2*cos(w_basex(i)) + 3*sin(w_basex(i)) + (2*cos(theta(6,i))*(2*cos(w_basex(i))*cos(theta(4,i))*cos(theta(5,i))*sin(w_basex(i)) - 2*cos(w_basex(i))*sin(w_basex(i))*sin(theta(4,i))*sin(theta(5,i))))/5 + (2*cos(theta(6,i))*(2*cos(w_basex(i))*cos(theta(4,i))*sin(w_basex(i))*sin(theta(5,i)) + 2*cos(w_basex(i))*cos(theta(5,i))*sin(w_basex(i))*sin(theta(4,i))))/5 + (6*cos(w_basex(i))*sin(w_basex(i))*sin(theta(4,i)))/5 + cos(w_basex(i))*cos(theta(4,i))*cos(theta(5,i))*sin(w_basex(i)) + cos(w_basex(i))*cos(theta(5,i))*sin(w_basex(i))*sin(theta(4,i)) + 2];
 
 
    g_pseudo = computeGeneralizedJacobian(x_base,xe_current,p1,p2,p3,p4,p5,p6);
   z1(:,i)=[0;0;1];z2(:,i)=[0;0;1];z3(:,i)=[0;0;1];
   z4(:,i)=[0;0;-1];z5(:,i)=[0;0;-1];z6(:,i)=[0;0;-1];

    dtheta=Jg_pseudo*[v_ee1x(i);v_ee1y(i);v_ee1z(i);w_ee1x(i);w_ee1y(i);w_ee1z(i);v_ee2x(i);v_ee2y(i);v_ee2z(i);w_ee2x(i);w_ee2y(i);w_ee2z(i)];
    [J_bm_v, J_bm_w] = computeBaseJacobian(R,JR,JT)
    v_base = [J_bm_v; J_bm_w] * dtheta; % 式15
    v0 = base_vel(1:3);
    vw_basex = base_vel(4);
    %积分更新位置
    x_base_new = x_base + v0 * dt;
    w_basex=w_basex+vw_basex*dt;
    % 更新关节角度和基座位姿
     theta = theta + dtheta * dt; % 更新机械臂关节
    x_base = x_base_new;
    
    
end

