function Jg_pseudo = computeGeneralizedJacobian(x_base,xe_current,p1,p2,p3,p4,p5,p6)
    Jg_pseudo = computeGeneralizedJacobian(x_base,xe_current,p1,p2,p3,p4,p5,p6,J_bmv,J_bmw);

% ��������ſɱȾ���ʽ16�� p1 p2 p3 p4 p5 p6 % ĩ���ڻ�е������ϵ�е�λ��
z1=[0;0;1];z2=[0;0;1];z3=[0;0;1];
z4=[0;0;-1];z5=[0;0;-1];z6=[0;0;-1];

% ���� r_{0e} = r_e - r_0 ���˴�r_eΪȫ�����꣬r_0Ϊ����λ�ã�
r_0e1 = xe_current(1:3) - x_base;
r_0e2 = xe_current(4:6) - x_base;

% ���� J_b^i ����ʽ16��
J_b1 = [eye(3), -skew(r_0e1);
        zeros(3), eye(3)];
J_b2= [eye(3), -skew(r_0e2);
        zeros(3), eye(3)];

Jb=[J_b1; J_b2];

% �����е���ſɱ� J_Te, J_Re �������ģ��ʵ�֣�
[J_Te] = [myCrossProduct(z1,xe_current(1:3)-p1),myCrossProduct(z2,xe_current(1:3)-p2),myCrossProduct(z3,xe_current(1:3)-p3),zeros(3,3);
          z1,z2,z3,zeros(3,3)];
[J_Re]= [zeros(3,3),myCrossProduct(z4,xe_current(4:6)-p4),myCrossProduct(z5,xe_current(4:6)-p5),myCrossProduct(z6,xe_current(4:6)-p6);
        zeros(3,3),z4,z5,z6];
    
[J_bm_v, J_bm_w] = computeBaseJacobian(R,JR,JT);

% �����ſɱ� J_g = J_b * [J_bm_v; J_bm_w] + [J_Te; J_Re]
Jg = J_b * [J_bm_v;J_bm_w] + [J_Te; J_Re];
Jg_pseudo = Jg'/(Jg*Jg')
end