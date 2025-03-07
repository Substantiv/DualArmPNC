function Jg_pseudo = computeGeneralizedJacobian(x_base,xe_current,p1,p2,p3,p4,p5,p6)
    % 计算广义雅可比矩阵（式16）
    % p1 p2 p3 p4 p5 p6 % 末端在机械臂坐标系中的位置
     z1=[0;0;1];z2=[0;0;1];z3=[0;0;1];
   z4=[0;0;-1];z5=[0;0;-1];z6=[0;0;-1];
    % 计算 r_{0e} = r_e - r_0 （此处r_e为全局坐标，r_0为基座位置）
    r_0e1 = xe_current(1:3) - x_base; 
    r_0e2 = xe_current(4:6) - x_base; 
    % 构建 J_b^i 矩阵（式16）
    J_b1 = [eye(3), -skew(r_0e1); 
           zeros(3), eye(3)];
    J_b2= [eye(3), -skew(r_0e2); 
           zeros(3), eye(3)];
    Jb=[J_b1; J_b2];
    
    % 计算机械臂雅可比 J_Te, J_Re （需具体模型实现）
    [J_Te] = [myCrossProduct(z1,xe_current(1:3)-p1),myCrossProduct(z2,xe_current(1:3)-p2),myCrossProduct(z3,xe_current(1:3)-p3),zeros(3,3);
             z1,z2,z3,zeros(3,3)];

    [J_Re]= [zeros(3,3),myCrossProduct(z4,xe_current(4:6)-p4),myCrossProduct(z5,xe_current(4:6)-p5),myCrossProduct(z6,xe_current(4:6)-p6);
               zeros(3,3),z4,z5,z6];
    [J_bm_v, J_bm_w] = computeBaseJacobian(R,JR,JT);
    % 广义雅可比 J_g = J_b * [J_bm_v; J_bm_w] + [J_Te; J_Re]
    Jg = J_b * [J_bm_v;J_bm_w] + [J_Te; J_Re]; 
    Jg_pseudo=Jg'/(Jg*Jg')
end