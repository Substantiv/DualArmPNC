function C = myCrossProduct(A, B)  
    % 检查输入向量是否为三维  
    % if length(A) ~= 3 || length(B) ~= 3  
    %     error('输入向量必须是三维的。');  
    % end  
    % 
    % 计算叉乘  
    C = [A(2)*B(3) - A(3)*B(2);  
         A(3)*B(1) - A(1)*B(3);  
         A(1)*B(2) - A(2)*B(1)];  
end