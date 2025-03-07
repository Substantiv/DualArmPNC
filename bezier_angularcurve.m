function [omega_x, omega_y, omega_z] = bezier_angularcurve(P4, P5, P6, P7)
    % 定义参数 s 的范围
    st = linspace(0, 1, 100);

    % 计算贝塞尔曲线上的点
   omega_x = (1-st).^3 * P4(1) + 3*(1-st).^2 .* st * P5(1) + 3*(1-st) .* st.^2 * P6(1) + st.^3 * P7(1);
omega_y = (1-st).^3 * P4(2) + 3*(1-st).^2 .* st * P5(2) + 3*(1-st) .* st.^2 * P6(2) + st.^3 * P7(2);
omega_z = (1-st).^3 * P4(3) + 3*(1-st).^2 .* st * P5(3) + 3*(1-st) .* st.^2 * P6(3) + st.^3 * P7(3);
end