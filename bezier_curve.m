function [x_m, y_m, z_m, dx, dy, dz, ddx, ddy, ddz] = bezier_curve(P0, P1, P2, P3)
    % 定义参数 s 的范围
    st = linspace(0, 1, 100);

    % 计算贝塞尔曲线上的点
    x_m = (1 - st).^3 * P0(1) + 3 * (1 - st).^2 .* st * P1(1) + 3 * (1 - st) .* st.^2 * P2(1) + st.^3 * P3(1);
    y_m = (1 - st).^3 * P0(2) + 3 * (1 - st).^2 .* st * P1(2) + 3 * (1 - st) .* st.^2 * P2(2) + st.^3 * P3(2);
    z_m = (1 - st).^3 * P0(3) + 3 * (1 - st).^2 .* st * P1(3) + 3 * (1 - st) .* st.^2 * P2(3) + st.^3 * P3(3);

    % 计算速度（一阶导数）
    dx = 3 * (1 - st).^2 * (P1(1) - P0(1)) + 6 * (1 - st) .* st * (P2(1) - P1(1)) + 3 * st.^2 * (P3(1) - P2(1));
    dy = 3 * (1 - st).^2 * (P1(2) - P0(2)) + 6 * (1 - st) .* st * (P2(2) - P1(2)) + 3 * st.^2 * (P3(2) - P2(2));
    dz = 3 * (1 - st).^2 * (P1(3) - P0(3)) + 6 * (1 - st) .* st * (P2(3) - P1(3)) + 3 * st.^2 * (P3(3) - P2(3));

    % 计算加速度（二阶导数）
    ddx = 6 * (1 - st) * (P2(1) - 2 * P1(1) + P0(1)) + 6 * st * (P3(1) - 2 * P2(1) + P1(1));
    ddy = 6 * (1 - st) * (P2(2) - 2 * P1(2) + P0(2)) + 6 * st * (P3(2) - 2 * P2(2) + P1(2));
    ddz = 6 * (1 - st) * (P2(3) - 2 * P1(3) + P0(3)) + 6 * st * (P3(3) - 2 * P2(3) + P1(3));
end