clear all
close all
clc

%% 系统参数初始化
% 状态变量初始化
x1 = 0.06;
x2 = 0.08;

% 自适应参数初始化
theta1j = 0.01 * ones(9, 1);
theta2j = 0.01 * ones(9, 1);

% 控制器参数
a1c = 0;
taud = 0.05;
r = 0.8;
c1 = 0.6;
c2 = 0.6;

b1 = 0.6;
b2 = 0.6;
tau1 = 45;
tau2 = 43;
rr1 =133;%113无噪声 133bucuo

% 时间参数
t = 0;
h = 0.005;


% 自适应律参数
beta1 = 0.14;
beta2 = 0.12;
beta3 = 0.16;
mu1 = 5.2;
mu2 = 5.3;
mu3 = 5.6;

% 饱和函数参数
alpha_adaptive = 0.88;
epsilon_star = 0.0006;  % 饱和边界
gamma_rho = 0.0035;      % 衰减率

%% 迭代学习控制参数
max_iterations =20;    % 最大迭代次数
ILC_gamma = 0.05;       % ILC学习增益
ILC_alpha = 0.08;       % ILC遗忘因子

% 自适应PI参数
chi1 = 10; 
chi2 = 14;            % 自适应增益
delta1 = 0.90; 
delta2 = 0.90;         % 自适应参数


% PI控制器参数
mu11 = 0.65; 
mu12 = 0.26;
mu21 = 0.68; 
mu22 = 0.36;

% 系统增益边界
g1min = 1; 
g1max = 1;
g1 = 1;
g2min =-0.8; 
g2max = 0.8;
varsigma_2=0.06;
iota1=0.01;

% 分数阶参数
varpho1 = 0.05; 
varpho2 = 0.08; 
alpha_val = 0.9;       % 分数阶阶次

noise_enabled =true;          
noise_intensity = 0.0005;     

%% 初始化ILC存储变量
ILC_a11_history = cell(max_iterations, 1);
ILC_a12_history = cell(max_iterations, 1);
ILC_a1_history = cell(max_iterations, 1);

ILC_u1_history = cell(max_iterations, 1);
ILC_u2_history = cell(max_iterations, 1);
ILC_u_history = cell(max_iterations, 1);

ILC_e1_history = cell(max_iterations, 1);
ILC_e2_history = cell(max_iterations, 1);
ILC_x1_history = cell(max_iterations, 1);
ILC_x2_history = cell(max_iterations, 1);
ILC_yd_history = cell(max_iterations, 1);
ILC_learning_update = 0;

ILC_theta1j_history = cell(max_iterations, 1);
ILC_theta2j_history = cell(max_iterations, 1);


% 性能指标存储
RMSE_e1 = zeros(1, max_iterations);
RMSE_e2 = zeros(1, max_iterations);
RMSE_z1 = zeros(1, max_iterations);
RMSE_z2 = zeros(1, max_iterations);
Max_e1 = zeros(1, max_iterations);
Max_e2 = zeros(1, max_iterations);
Max_z1 = zeros(1, max_iterations);
Max_z2 = zeros(1, max_iterations);


% 能量消耗存储
J = zeros(1, max_iterations);  % 目标函数（累计误差）
R = zeros(1, max_iterations);  % 总目标函数（累计误差）
energy_consumption = zeros(1, max_iterations);

% PI增益存储
kp1_history = cell(max_iterations, 1);
ki1_history = cell(max_iterations, 1);
kp2_history = cell(max_iterations, 1);
ki2_history = cell(max_iterations, 1);
k01_history = cell(max_iterations, 1);
k02_history = cell(max_iterations, 1);

% 自适应参数存储
deltap1_history = cell(max_iterations, 1);
deltap2_history = cell(max_iterations, 1);
deltaki1_history = cell(max_iterations, 1);
deltaki2_history = cell(max_iterations, 1);
deltak01_history = cell(max_iterations, 1);
deltak02_history = cell(max_iterations, 1);

% 噪声存储 - 在循环外初始化
if noise_enabled
    noise1_history = zeros(max_iterations, 4000);
    noise2_history = zeros(max_iterations, 4000);
end


for iteration = 1:max_iterations
    fprintf('正在进行第 %d 次迭代...\n', iteration);
    
    % 重置状态变量（保持初始条件一致）
    if iteration > 1
        x1 = 0.06;
        x2 = 0.08;
        a1c = 0;
        t = 0;
    end
      
    % 预分配存储空间
    n_steps = 4000;
    
    a11_ilc = zeros(1, n_steps);
    a12_ilc = zeros(1, n_steps);
    a1_ilc = zeros(1, n_steps);
    
    u1_ilc = zeros(1, n_steps);
    u2_ilc = zeros(1, n_steps);
    u_ilc = zeros(1, n_steps);
    
    e1_store = zeros(1, n_steps);
    e2_store = zeros(1, n_steps);
    z1_store = zeros(1, n_steps);
    z2_store = zeros(1, n_steps);
    u_store = zeros(1, n_steps);
    a1_store = zeros(1, n_steps);
    yd_store = zeros(1, n_steps);
    x1_store = zeros(1, n_steps);
    x2_store = zeros(1, n_steps);
    
    
    % PI增益存储
    kp1_iter = zeros(1, n_steps);
    ki1_iter = zeros(1, n_steps);
    kp2_iter = zeros(1, n_steps);
    ki2_iter = zeros(1, n_steps);
    k01_iter = zeros(1, n_steps);
    k02_iter = zeros(1, n_steps);
    
    % 自适应参数存储
    deltap1_iter = zeros(1, n_steps);
    deltap2_iter = zeros(1, n_steps);
    deltaki1_iter = zeros(1, n_steps);
    deltaki2_iter = zeros(1, n_steps);
    deltak01_iter = zeros(1, n_steps);
    deltak02_iter = zeros(1, n_steps);
    
    % 预分配dottheta1j和dottheta2j
    dottheta1j = zeros(9, n_steps);
    dottheta2j = zeros(9, n_steps);
    
    % 初始化当前迭代的theta
    theta1j_iter = zeros(9, n_steps+1);
    theta2j_iter = zeros(9, n_steps+1);
    theta1j_iter(:,1) = theta1j;
    theta2j_iter(:,1) = theta2j;
    
    for n = 1:n_steps
        t(n+1) = t(n) + h;
        
        % 模糊基函数计算
        kk = 1;
        for ii1 = -2:2:2
            for ii2 = -2:2:2 
                kp1(kk) = gaussmf(x1(n), [0.9, ii1]) * gaussmf(x2(n), [0.9, ii2]);
                kk = kk + 1;
            end
        end
       
        pp1 = sum(kp1);
        psi1(:, n) = kp1' / pp1;
    
        % 参考信号生成
        yd(n) = 0.6* sin(t(n));
        dotyd(n) = 0.6* cos(t(n));
        yd_store(n) = yd(n);
        
        % 跟踪误差计算
        e1(n) = x1(n) - yd(n);
        e2(n) = x2(n) - a1c(n);
        
        % 存储误差用于ILC
        e1_store(n) = e1(n);
        e2_store(n) = e2(n);
        x1_store(n) = x1(n);
        x2_store(n) = x2(n);
        
        
        % 时滞误差积分
        if t(n) > taud
            e11(n) = h * sum(e1(1:n));
        else
            e11(n) = 0;
        end
    
        if t(n) > taud
            e21(n) = h * sum(e2(1:n));
        else
            e21(n) = 0;
        end
    
        % 饱和函数计算
        rho_tk(n) = epsilon_star * exp(-gamma_rho*t(n));
        sat_e1(n) = enhanced_saturate(e1(n), rho_tk(n));
        sat_e2(n) = enhanced_saturate(e2(n), rho_tk(n));

        % 屏障Lyapunov函数变量
        z1(n) = mu11 * e1(n) + mu12 * e11(n) - rho_tk(n) * sat_e1(n);
        z2(n) = mu21 * e2(n) + mu22 * e21(n) - rho_tk(n) * sat_e2(n);
        
        % 屏障函数
        barrier1(n) = b1^2 - z1(n)^2;
        barrier2(n) = b2^2 - z2(n)^2;
        
        % 自适应参数更新
        deltap1(n) = -(mu11 * delta1 * g1max) / (g1min * barrier1(n)^2) - ...
                     (mu11 * g1max^2 * delta1 * norm(theta1j_iter(:, n)) * norm(psi1(:, n))^2) / ...
                     (2 * chi1^2 * barrier1(n)^2) + 1.5 / g1min * z1(n)^2;
        
        deltap2(n) = -(mu21 * delta2 * g2max) / (g2min * barrier2(n)^2) + ...
                     1.5 / g2min * z2(n)^2 - ...
                     (mu21 * g2max^2 * delta2 * norm(theta2j_iter(:, n)) * norm(psi1(:, n))^2) / ...
                     (2 * chi2^2 * barrier2(n)^2);
                 
         deltaki1(n) = mu12/mu11 * deltap1(n);
         deltaki2(n) = mu22/mu21 * deltap2(n);

       
        % 存储自适应参数
        deltap1_iter(n) = deltap1(n);
        deltap2_iter(n) = deltap2(n);
        deltaki1_iter(n) = deltaki1(n);
        deltaki2_iter(n) = deltaki2(n);
       
        
        % 在主循环内部，在计算z1(n)和z2(n)后添加：
        z1_store(n) = z1(n);
        z2_store(n) = z2(n);

       
        kp1_val = -1/mu11 - 1/g1min * b1^2; 
        kp2_val = -1/mu21 - 1/g1min * b2^2;%+1/(mu21*g2min)*(varsigma_2*iota1^2)/(2*iteration); 
        
        % 存储PI增益
        kp1_iter(n) = kp1_val + deltap1(n);
        ki1_iter(n) = mu12/mu11 * deltap1(n);
        kp2_iter(n) = kp2_val + deltap2(n);
        ki2_iter(n) = mu22/mu21 * deltap2(n);
        
        deltaki1(n) = mu12/mu11 * deltap1(n);%辅助增益
        deltaki2(n) = mu22/mu11 * deltap2(n);
         
        k01 = 1/mu11 * kp1_val;
        k02 = 1/mu21 * kp2_val;
        
        % 存储k01和k02
        k01_iter(n) = k01;
        k02_iter(n) = k02;
        
        deltak01(n) = 1/mu11 * deltap1(n);
        deltak02(n) = 1/mu21 * deltap2(n);
        
        deltak01_iter(n) = deltak01(n);
        deltak02_iter(n) = deltak02(n);
 
        % 虚拟控制律
        a11(n) = - (k01 + deltak01(n)) * rho_tk(n) * sat_e1(n);
        a12(n) = kp1_val * z1(n) + deltap1(n) * z1(n);
        a1(n) = a11(n) + a12(n);
        
        % 存储a1用于分析
        a1_store(n) = a1(n);
        a11_store(n) = a11(n);
        a12_store(n) = a12(n);
        
        
        % 存储 u1 和 u2
        a11_ilc(n) = a11(n);
        a12_ilc(n) = a12(n);
        a1_ilc(n) = ILC_learning_update;
        
        % 分数阶滤波器
        dota1c(n) = -rr1 * (a1c(n) - a1(n));
        
        % 分数阶积分
        jj1 = (t(n+1) - t(1:n)).^(alpha_val-1) .* dota1c(1:n) / gamma(alpha_val);
        JJ1 = h * sum(jj1);
        a1c(n+1) = a1c(n) + JJ1;
      
        % 自适应律
        adaptive1_term = (mu11 * beta1 * g1max^2 * delta1 * z1(n) * psi1(:, n)) / ...
                         (2 * chi1^2 * barrier1(n)^2 * (1 - varpho1));
        adaptive2_term = (mu21 * beta2 * g2max^2 * delta2 * z2(n) * psi1(:, n)) / ...
                         (2 * chi2^2 * barrier2(n)^2 * (1 - varpho2));
        
        dottheta1j(:, n) = adaptive1_term - (varpho1/(1-varpho1)) * theta1j_iter(:, n);
        dottheta2j(:, n) = adaptive2_term - (varpho2/(1-varpho2)) * theta2j_iter(:, n);
        
        % 分数阶积分更新自适应参数
        if n > 1
            time_diff = t(n+1) - t(1:n);
            aaa1 = time_diff.^(alpha_val-1) / gamma(alpha_val);
            
            integral_term = zeros(9, 1);
            for i = 1:n
                if i <= size(dottheta1j, 2)
                    integral_term = integral_term + aaa1(i) * dottheta1j(:, i);
                end
            end
            
            theta1j_iter(:, n+1) = theta1j_iter(:, 1) + h * integral_term;
            
            aaa2 = time_diff.^(alpha_val-1) / gamma(alpha_val);
            
            integral_term2 = zeros(9, 1);
            for i = 1:n
                if i <= size(dottheta2j, 2)
                    integral_term2 = integral_term2 + aaa2(i) * dottheta2j(:, i);
                end
            end
            
            theta2j_iter(:, n+1) = theta2j_iter(:, 1) + h * integral_term2;
        else
            theta1j_iter(:, n+1) = theta1j_iter(:, n);
            theta2j_iter(:, n+1) = theta2j_iter(:, n);
        end
        
        % 控制律
        u1(n) =  - (k02 + deltak02(n)) * rho_tk(n) * sat_e2(n);
        u2(n) = kp2_val * z2(n) + deltap2(n) * z2(n);
        u(n) = u1(n) + u2(n);
        
        % 应用ILC更新
        if iteration > 1
            prev_e1 = ILC_e1_history{iteration-1};
            prev_e2 = ILC_e2_history{iteration-1};
            
            if n > 1
                ILC_learning_update = ILC_gamma * prev_e1(n) + ILC_alpha * ILC_learning_update;
            end
            
            u(n) = u(n) + ILC_learning_update;
        end
        
        % 存储控制输入
        u_store(n) = u(n);
        u1_store(n) = u1(n);
        u2_store(n) = u2(n);
        % 存储 u1 和 u2
        u1_ilc(n) = u1(n);
        u2_ilc(n) = u2(n);
        u_ilc(n) = ILC_learning_update;
        
        % 噪声
        if noise_enabled
            noise1_val = noise_intensity * randn;
            noise2_val = noise_intensity * randn;
            % 存储噪声值
            noise1_history(iteration, n) = noise1_val;
            noise2_history(iteration, n) = noise2_val;
        else
            noise1_val = 0; 
            noise2_val = 0;
        end
        
        % 系统动力学
        f1(n) = 0;
        f2(n) = -x1(n) * x2(n)^2;
        g2(n) = 0.8 * sin(x1(n)^2* x2(n));
        pho1(n) =0.001*sin(t(n));
        
        dotx1(n) = g1 * x2(n) + f1(n) + noise1_val;  
        dotx2(n) =u(n) + f2(n) +pho1(n) + noise2_val;
        
        % 分数阶积分更新状态
        if n > 1
            time_diff = t(n+1) - t(1:n);
            j1_coeff = time_diff.^(alpha_val-1) / gamma(alpha_val);
            j2_coeff = time_diff.^(alpha_val-1) / gamma(alpha_val);
            
            J1 = h * sum(j1_coeff .* dotx1(1:n));
            J2 = h * sum(j2_coeff .* dotx2(1:n));
            
            x1(n+1) = x1(1) + J1;
            x2(n+1) = x2(1) + J2;
        else
            x1(n+1) = x1(n);
            x2(n+1) = x2(n);
        end
    end
    
    % 存储本次迭代的数据
    ILC_u_history{iteration} = u_store;
    ILC_e1_history{iteration} = e1_store;
    ILC_e2_history{iteration} = e2_store;
    ILC_x1_history{iteration} = x1_store;
    ILC_x2_history{iteration} = x2_store;
    ILC_yd_history{iteration} = yd_store;
    ILC_z1_history{iteration} = z1_store;
    ILC_z2_history{iteration} = z2_store;
    
    %存储本次迭代的数据
    ILC_u1_history{iteration} = u1_ilc;
    ILC_u2_history{iteration} = u2_ilc;
    ILC_a11_history{iteration} = a11_ilc;
    ILC_a12_history{iteration} = a12_ilc;
    ILC_a1_history{iteration} = a1_store;
    
    % 存储自适应参数
    ILC_theta1j_history{iteration} = theta1j_iter;
    ILC_theta2j_history{iteration} = theta2j_iter;

    
    % 存储PI增益
    kp1_history{iteration} = kp1_iter;
    ki1_history{iteration} = ki1_iter;
    kp2_history{iteration} = kp2_iter;
    ki2_history{iteration} = ki2_iter;
    k01_history{iteration} = k01_iter;
    k02_history{iteration} = k02_iter;
    
    % 存储自适应参数
    deltap1_history{iteration} = deltap1_iter;
    deltap2_history{iteration} = deltap2_iter;
    deltaki1_history{iteration} = deltaki1_iter;
    deltaki2_history{iteration} = deltaki2_iter;
    deltak01_history{iteration} = deltak01_iter;
    deltak02_history{iteration} = deltak02_iter;
    
    m1=0.001;m2=0.05;m3=0.06;
    
    % 计算性能指标
    RMSE_e1(iteration) = sqrt(mean(e1_store.^2));
    RMSE_e2(iteration) = sqrt(mean(e2_store.^2));
    RMSE_z1(iteration) = sqrt(mean(z1_store.^2));
    RMSE_z2(iteration) = sqrt(mean(z2_store.^2));
    Max_z1(iteration) = max(abs(z1_store));
    Max_z2(iteration) = max(abs(z2_store));
    J(iteration) = m1*sum(z1_store.^2) + m1*sum(z2_store.^2)+m3*Max_z1(iteration)+m3*Max_z2(iteration);  % 目标函数（累计误差平方和）Tracking Accuracy
    energy_consumption(iteration) = m2*sum(u_store.^2) * h+m2*sum(a1_store.^2) * h;  % 能量消耗（控制输入的平方积分）
    R(iteration)=J(iteration) +energy_consumption(iteration)+m3*Max_z1(iteration)+m3*Max_z2(iteration);
    
    fprintf('迭代 %d 完成: e1_RMSE=%.4f, e2_RMSE=%.4f\n', iteration, RMSE_e1(iteration), RMSE_e2(iteration));
end


% 添加性能指标文本
fprintf('\nBarrier Lyapunov Variable Performance:\n');
fprintf('RMSE of z1 improved from %.4f to %.4f (Improvement Rate: %.2f%%)\n', ...
    RMSE_z1(1), RMSE_z1(end), (RMSE_z1(1)-RMSE_z1(end))/RMSE_z1(1)*100);
fprintf('RMSE of z2 improved from %.4f to %.4f (Improvement Rate: %.2f%%)\n', ...
    RMSE_z2(1), RMSE_z2(end), (RMSE_z2(1)-RMSE_z2(end))/RMSE_z2(1)*100);

fprintf('\nIterative Learning Control Completed!\n');
fprintf('RMSE of e1 improved from %.4f to %.4f (Improvement Rate: %.2f%%)\n', ...
    RMSE_e1(1), RMSE_e1(end), (RMSE_e1(1)-RMSE_e1(end))/RMSE_e1(1)*100);
fprintf('RMSE of e2 improved from %.4f to %.4f (Improvement Rate: %.2f%%)\n', ...
    RMSE_e2(1), RMSE_e2(end), (RMSE_e2(1)-RMSE_e2(end))/RMSE_e2(1)*100);
fprintf('Objective Function J improved from %.4f to %.4f (Improvement Rate: %.2f%%)\n', ...
    J(1), J(end), (J(1)-J(end))/J(1)*100);

%% 辅助函数
function sat_val = enhanced_saturate(x, rho)
    if abs(x) > rho
        sat_val = sign(x);
    else
        sat_val = x/rho;
    end
end

