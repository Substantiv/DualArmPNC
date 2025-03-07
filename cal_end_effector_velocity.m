function [v_ee1, w_ee1,v_ee2, w_ee2,x_ee1,y_ee1,z_ee1,x_ee2,y_ee2,z_ee2] = cal_end_effector_velocity(x_m, y_m, z_m, dm_t, omega_x, omega_y, omega_z)
    % 定义从模块到末端执行器1的变换矩阵
    T_ee1_m = [1 0 0 0;
               0 1 0 0;
               0 0 1 0.258/2;
               0 0 0 1];
    % 预分配变量
    num_points = length(x_m);
    x_ee1 = zeros(1, num_points);
    y_ee1 = zeros(1, num_points);
    z_ee1 = zeros(1, num_points);
    p_m_ee1 = zeros(3, num_points);
    J_ee1_m = cell(1, num_points);
    vw_ee1 = zeros(6, num_points);

    % 循环计算每个点的相关信息
    for i = 1:num_points
        % 计算从全局坐标系到模块的变换矩阵
        T_g_m{i} = [cos(omega_x(i))*cos(omega_y(i))  cos(omega_x(i))*sin(omega_y(i))*sin(omega_z(i))-sin(omega_x(i))*cos(omega_z(i))    cos(omega_x(i))*sin(omega_y(i))*cos(omega_z(i))+sin(omega_x(i))*sin(omega_z(i))   x_m(i);
                    sin(omega_x(i))*cos(omega_y(i))  -sin(omega_x(i))*sin(omega_y(i))*sin(omega_z(i))+cos(omega_x(i))*cos(omega_z(i))  -sin(omega_x(i))*sin(omega_y(i))*cos(omega_z(i))-cos(omega_x(i))*sin(omega_z(i))   y_m(i);
                    -sin(omega_y(i))                  cos(omega_y(i))*sin(omega_z(i)) cos(omega_y(i))*cos(omega_z(i)) z_m(i);
                    0 0 0 1];
        % 计算从全局坐标系到末端执行器1的变换矩阵
        T_g_ee1{i} = T_g_m{i} * inv(T_ee1_m);
        % 提取末端执行器1的位置信息
        x_ee1(i) = T_g_ee1{i}(1,4);
        y_ee1(i) = T_g_ee1{i}(2,4);
        z_ee1(i) = T_g_ee1{i}(3,4);
        % 计算模块到末端执行器1的位置向量
        p_m_ee1(:,i) = [x_ee1(i)-x_m(i);
                        y_ee1(i)-y_m(i);
                        z_ee1(i)-z_m(i)];
        % 计算雅可比矩阵
        J_ee1_m{i} = [eye(3,3) -skew(p_m_ee1(:,i));
                      zeros(3,3) eye(3,3)]; 
        % 计算末端执行器1的速度和角速度
        vw_ee1(:,i) = inv(J_ee1_m{i}) * [dm_t(1,i);dm_t(2,i);dm_t(3,i);omega_x(i);omega_y(i);omega_z(i)];
    end

    % 提取末端执行器1的线速度和角速度
    v_ee1 = vw_ee1(1:3,:);
    w_ee1 = vw_ee1(4:6,:);
   
T_ee2_m=[1 0 0 0;
         0 1 0 0;
         0 0 1 -0.258/2;
         0 0 0 1];
for i=1:100
     T_g_m{i}=[cos(omega_x(i))*cos(omega_y(i))  cos(omega_x(i))*sin(omega_y(i))*sin(omega_z(i))-sin(omega_x(i))*cos(omega_z(i))    cos(omega_x(i))*sin(omega_y(i))*cos(omega_z(i))+sin(omega_x(i))*sin(omega_z(i))   x_m(i);
                  sin(omega_x(i))*cos(omega_y(i))  -sin(omega_x(i))*sin(omega_y(i))*sin(omega_z(i))+cos(omega_x(i))*cos(omega_z(i))  -sin(omega_x(i))*sin(omega_y(i))*cos(omega_z(i))-cos(omega_x(i))*sin(omega_z(i))   y_m(i);
                  -sin(omega_y(i))                  cos(omega_y(i))*sin(omega_z(i)) cos(omega_y(i))*cos(omega_z(i)) z_m(i);
                  0 0 0 1];
    T_g_ee2{i}=T_g_m{i}*inv(T_ee2_m);
    x_ee2(i)=T_g_ee2{i}(1,4);
    y_ee2(i)=T_g_ee2{i}(2,4);
    z_ee2(i)=T_g_ee2{i}(3,4);
    p_m_ee2(:,i)=[x_ee2(i)-x_m(i);
                  y_ee2(i)-y_m(i);
                  z_ee2(i)-z_m(i)];
  J_ee2_m{i}=[eye(3,3) -skew(p_m_ee2(:,i));
              zeros(3,3) eye(3,3)] ; 
  vw_ee2(:,i)=inv(J_ee2_m{i})*[dm_t(1,i);dm_t(2,i);dm_t(3,i);omega_x(i);omega_y(i);omega_z(i)]
end
% end-effector2速度
v_ee2=vw_ee2(1:3,:);
w_ee2=vw_ee2(4:6,:);


end