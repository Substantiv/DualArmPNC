function [J_bm_v, J_bm_w] = computeBaseJacobian(R,JR,JT)
% 定义参数
I0xx=10.4;I0yy=10.4;I0zz=10.4;
I1xx=2*10^(-4);I1yy=2*10^(-4);I1zz=2*10^(-4);
I2xx=3.5*10^(-4);I2yy=3.5*10^(-4);I2zz=2*10^(-4);
I3xx=3.5*10^(-4);I3yy=3.7*10^(-4);I3zz=2*10^(-4);
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
m0=10.18;m1=0.5;m2=0.6;m3=0.6;
m4=m1;m5=m2;m6=m3;
a1=0.6;a2=0.5;a3=0.4;
a4=a1;a5=a2;a6=a3;
b=0.258;
m=[m1 m2 m3 m4 m5 m6 ];
mt=0.2;

% 输入:
%   q_arm: 当前关节角度向量 params: 包含动力学参数的结构体，字段包括:
%       m0, mT, I0, IT - 基座和末端的质量及惯性张量 m_j{i}: 第i段的各关节质量 (cell数组) I_j{i}:
%       第i段的各关节惯性张量 (cell数组) r0: 基座位置 (全局坐标系) r0_j{i}: 第i段各关节相对于基座的位置
%       (cell数组) r0T: 末端相对于基座的位置 J_Tj{i}, J_Rj{i}: 第i段的平动/转动雅可比矩阵 (cell数组)
%       J_TT, J_RT: 末端的平动/转动雅可比矩阵
% 输出:
%   J_bm_v: 式(15)中的基座平动雅可比矩阵 (3x(n1+n2)) J_bm_w: 式(15)中的基座转动雅可比矩阵
%   (3x(n1+n2))

% --- 步骤1: 计算总质量 M ---
m_tott=(m0+m1+m2+m3+m4+m5+m6+mt);

% --- 步骤2: 计算质心位置 r0g --- 假设rg是全局质心，需通过各部分质心计算 此处简化为各部分质心加权平均
r0g=(m1*(p1c-Tb_I(1:3,4))+m2*(p2c-Tb_I(1:3,4))+m3*(p3c-Tb_I(1:3,4))+m4*(p4c-Tb_I(1:3,4))+m5*(p5c-Tb_I(1:3,4))+m6*(p6c-Tb_I(1:3,4))+mt*(pb-Tb_I(1:3,4)))/m_tott;

% --- 步骤3: 计算平动惯量 J_Tomega^i (式中的J_{Tω}^i) ---
JT_w0=zeros(3,3);
for j=1:3
    JT_w0=JT_w0+m(j)*JT(:,3*j-2:3*j);
end
JT_w1=JT_w0+mt*JT(:,19:21);
JT_w1=vpa(JT_w1,5);

JT_w02=zeros(3,3);
for j=4:6
    JT_w02=JT_w02+m(j)*JT(:,3*j-2:3*j);
end
JT_w2=JT_w02;
JT_w2=vpa(JT_w2,5);
% --- 步骤4: 计算转动惯量 H_omega (式中的Hω) ---
Hw=zeros(3,3);
for j=1:6
    Hw=Hw+(R(:,3*j-2:3*j)*I(:,3*j-2:3*j)*R(:,3*j-2:3*j).'-m(j)*skew(pc(:,j)-Tb_I(1:3,4))*skew(pc(:,j)-Tb_I(1:3,4)));
end
HW=Hw+I0+R(:,19:21)*I(:,19:21)*R(:,19:21)-mt*skew(pb-Tb_I(1:3,4))*skew(pb-Tb_I(1:3,4));
Hs=m_tott*skew(r0g)*skew(r0g)+HW;
Hs=vpa(Hs,5);

% --- 步骤5: 计算耦合惯量 H_omega_phi^i (式中的Hωφ^i) ---
F3=zeros(3,3);
for j=1:3
    F3=F3+R(:,3*j-2:3*j)*I(:,3*j-2:3*j)*R(:,3*j-2:3*j).'*JR(:,3*j-2:3*j)+m(j)*skew(pc(:,j)-Tb_I(1:3,4))*JT(:,3*j-2:3*j);
end
H_wpsi1=F3+R(:,19:21)*I(:,19:21)*R(:,19:21)*JR(:,19:21)+mt*skew(pb-Tb_I(1:3,4))*JT(:,19:21);
H_wpsi1=vpa(H_wpsi1,5);
F4=zeros(3,3);
for j=4:6
    F4=F4+R(:,3*j-2:3*j)*I(:,3*j-2:3*j)*R(:,3*j-2:3*j).'*JR(:,3*j-2:3*j)+m(j)*skew(pc(:,j)-Tb_I(1:3,4))*JT(:,3*j-2:3*j);
end
H_wpsi2=F4;
H_wpsi2=vpa(H_wpsi2,5);

% --- 步骤6: 计算 H_q^i (式中的H_q^i = H_omega_phi^i - r̃0g * J_Tomega^i) ---
H_q1=H_wpsi1-skew(r0g)*JT_w1;
H_q1=vpa(H_q1,5);
H_q2=H_wpsi2-skew(r0g)*JT_w2;
H_q2=vpa(H_q2,5);

% --- 步骤7: 计算基座雅可比 J_bm_v 和 J_bm_w (式15中的J_{bm_v}^i, J_{bm_w}^i) ---
% 合并两个机械臂段的雅可比
n1 = size(3, 2); % 第一段的关节数
n2 = size(3, 2); % 第二段的关节数

J_bm_v = zeros(3, n1+n2);
J_bm_w = zeros(3, n1+n2);

bmv11=-skew(r0g)/Hs*H_q1-JT_w1;
bmv12=-skew(r0g)/Hs*H_q2-JT_w2;
bmw11=-inv(Hs)*H_q1;
bmw12=-inv(Hs)*H_q2;
J_bmv=[bmv11  bmv12];%此处应用数值计算
J_bmv=vpa(J_bmv,5);
J_bmw=[bmw11  bmw12];%此处应用数值计算
J_bmw=vpa(J_bmv,5);

end