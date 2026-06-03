%% Force-and-moment-based Model Predictive Control
% ==========================================================================
% 基于力-力矩的模型预测控制 (Locomotion MPC) — 论文核心算法实现
% ==========================================================================
%
% 本文件实现了 Junheng Li & Quan Nguyen (2021 CDC) 论文的核心算法:
% "Force-and-moment-based model predictive control for achieving highly
%  dynamic locomotion on bipedal robots"
%
% 核心思想:
%   1. 将机器人简化为单刚体动力学 (SRBD: Single Rigid Body Dynamics)
%      - 身体 = 单刚体, 腿 = 无质量连杆
%      - 13维状态空间: [欧拉角, 位置, 角速度, 线速度, 重力]
%   2. 同时优化地面反力(GRF)和力矩(GRM)
%      - 12维控制量: [F_R, F_L, M_R, M_L]  (每足6维)
%   3. 将MPC转化为稠密QP (Quadratic Program)
%      - 消去状态变量 → 仅剩控制量
%      - 使用 CasADi + qpOASES 求解
%   4. 约束设计:
%      - 摩擦锥: |F_xy| ≤ μF_z
%      - 线足约束 (Line-Foot): 约束CoP在支撑域内
%      - 力/力矩/关节力矩限幅
%      - 步态调度: 摆动腿控制量清零
%   5. 输出GRFM → 经JacobianMapping转换为关节力矩τ
%
% 调用方式: Simulink模型中每个MPC采样步调用一次 (dt=0.04s)
% 输入 uin: [xdes(12), x(12), q(10), foot(12)] = 46维
% 输出 u:  GRFM(12) = [F_R, F_L, M_R, M_L]

% CITE AS:
% @inproceedings{li2021force,
%  title={Force-and-moment-based model predictive control for achieving highly dynamic locomotion on bipedal robots},
% author={Li, Junheng and Nguyen, Quan},
% booktitle={2021 60th IEEE Conference on Decision and Control (CDC)},
% pages={1024--1030},
% year={2021},
% organization={IEEE}
%}

function u = LocomotionMPC_qpoases(uin)
tic  % 开始计时 (函数总耗时)

%% ===== 第1步: 加载MPC全局参数 =====
global i_MPC_var dt_MPC_vec gait x_traj_IC Contact_Jacobian Rotm_foot addArm
% i_MPC_var:    当前MPC阶段索引 (由AssignMPCStage.m每步更新)
% dt_MPC_vec:   每个MPC阶段的步长向量 (当前固定为 [0.04,...,0.04])
% gait:         步态模式 (0=站立, 1=行走)
% x_traj_IC:    初始参考轨迹 (k=1时使用)
% Contact_Jacobian: 接触雅可比函数句柄 Jc(q_R, q_L, R_body)
% Rotm_foot:     足端姿态函数句柄
% addArm:       是否包含手臂 (影响质量和惯量)

k = i_MPC_var;   % 当前MPC时域索引 (horizon index)
h = 10;          % 预测时域长度 (horizon length): 向前预测10步
g = 9.81;        % 重力加速度 (m/s²)

%% ===== 第2步: QP问题规模预定义 =====
import casadi.*
numVar = 12*h;    % 优化变量数 = 12(每步控制量) × 10(时域) = 120
numCons = 34*h;   % 约束数量   = 34(每步约束)   × 10(时域) = 340
% 每步34个约束的构成:
%   摩擦锥: 8个  + 力限幅: 2个 + 线足: 12个 + 关节力矩: 10个 + Mx=0: 2个 = 34个

%% ===== 第3步: 解析输入信号 =====
% Simulink输入 uin 的结构 (共46维):
%   uin(1:12)  = xdes: 期望状态 [eul(3), p(3), omega(3), v(3)]
%   uin(13:24) = x:    当前状态 [eul(3), p(3), omega(3), v(3)]
%   uin(25:34) = q:    关节角度 [q_right(5), q_left(5)]
%   uin(35:46) = foot: 足端位置和速度 [p_L(3), v_L(3), p_R(3), v_R(3)]

xdes = uin(1:12);   % 期望状态 (desired SRBD states)
x    = uin(13:24);  % 当前实际状态 (actual SRBD states)
q    = uin(25:34);  % 10个关节角度 [hip1_R, hip2_R, thigh_R, calf_R, toe_R,
                    %               hip1_L, hip2_L, thigh_L, calf_L, toe_L]
foot = uin(35:46);  % 足端位置/速度 [p_L(3), v_L(3), p_R(3), v_R(3)]

%% ===== 第4步: 姿态处理 =====
% 提取欧拉角 [roll; pitch; yaw]
eul     = x(1:3);     % 当前姿态 (ZYX欧拉角)
eul_des = xdes(1:3);  % 期望姿态

% 计算身体旋转矩阵 R ∈ SO(3)
% eul2rotm使用ZYX约定: R = Rz(yaw) * Ry(pitch) * Rx(roll)
% flip(eul')将 [roll,pitch,yaw] 翻转为 [yaw,pitch,roll] 以匹配ZYX顺序
R = eul2rotm(flip(eul'));

%% ===== 第5步: 偏航角环绕处理 (Yaw Wrap-Around Correction) =====
% 问题: 当机器人转向超过π时, 欧拉角会从π跳变到-π (或反之)
% 这会导致控制器误以为需要瞬间回转, 产生错误的大幅控制信号
% 解决: 将当前偏航角调整到与期望偏航角的"最短路径"上
%
% 示例: eul_des(3) = 170°, eul(3) = -170°
%       实际只需要左转20°, 但数值上看是 -340° 的差距
%       通过+360°修正: eul(3) = 190°, 差距 = -20° (正确)
yaw_correction = 0;
while yaw_correction == 0
    if eul_des(3,1) - eul(3,1) > pi
        eul(3,1) = eul(3,1) + 2*pi;   % 当前偏航角落后太多, +360°
    elseif eul_des(3,1) - eul(3,1) < -pi
        eul(3,1) = eul(3,1) - 2*pi;   % 当前偏航角超前太多, -360°
    else
        yaw_correction = 1;            % 修正完成
    end
end

%% ===== 第6步: 状态增广 (加入重力) =====
% 将13个状态合并: [eul(3), p(3), omega(3), v(3), g(1)]
% g作为常值状态, 简化重力在动力学方程中的表达
xdes = [xdes; g];   % 期望状态 (13维), g期望值始终=9.81
x    = [x; g];      % 当前状态 (13维)

%% ===== 第7步: 生成参考轨迹 =====
% 首次调用时, 使用初始参考轨迹 x_traj_IC 初始化
if k == 1
    x_traj = x_traj_IC;
end

% 生成未来h步的参考状态轨迹 x_traj (13×h)
% 原理: 基于当前状态 + 期望速度, 线性外推未来h步的期望位置和姿态
%       x_traj(j, i+1) = x(j) + xdes(6+j) * Σ(dt)  (位置/姿态)
%       x_traj(6+j,i+1) = xdes(6+j)                  (速度/角速度不变)
x_traj = Calc_x_traj(xdes, x, h, k); % desired trajectory

% 足端参考位置 (简化处理: 假设足端在髋关节正下方)
% 格式: [p_R(3); p_L(3)] = 6维
foot_traj = [foot(1:3); foot(7:9)]; % assume foot under hip
% 备选方案: 使用 Calc_foot_traj_3Dwalking 获得更精确的足端位置估计
% 但该函数在3D行走场景下可能存在bug, 因此默认使用简化版本

%% ===== 第8步: 机器人物理参数 =====
% 摩擦系数: 决定摩擦锥的陡峭程度
mu = 0.5; % friction coefficient

% 根据是否包含手臂选择不同的质量/惯量参数
if addArm
    % 完整人形机器人 (包含手臂)
    m = 16.4;  % 总质量 (kg): 躯干 + 手臂 + 髋部 + 大腿
    Ib = diag([0.932, 0.9420, 0.0711]); % SRBD转动惯量 (kg·m²)
    % Ib对角线: [I_xx, I_yy, I_zz]
    % I_xx = 0.932:  绕X轴 (前后翻滚)
    % I_yy = 0.942:  绕Y轴 (左右侧翻)
    % I_zz = 0.0711: 绕Z轴 (偏航转动) — 最小, 因为身体扁平
else
    % 仅腿部机器人 (无手臂)
    % m = 躯干(5.75) + 2×(髋1(0.835) + 髋2(0.764) + 大腿(1.613) + 小腿(0.12) + 脚(0.08))
    m = 5.75 + 2*(0.835 + 0.764 + 1.613 + 0.12 + 0.08); % ≈ 12.574 kg
    Ib = diag([0.5413, 0.5200, 0.0691]); % SRBD转动惯量 (不含手臂)
end

% 足底法向力限制 (N)
Fmax = 500;  % 上限: 防止过大的接触力导致关节力矩饱和
Fmin = 0;    % 下限: 法向力不能为负 (地面只能推, 不能拉)

%% ===== 第9步: 接触雅可比 & 关节力矩映射 =====
% 将3×3旋转矩阵R展平为9×1向量, 传给Contact_Jacobian函数
RR = reshape(R, [9, 1]);

% 调用预编译的接触雅可比函数, 得到当前构型下的 Jc (12×10)
% Jc = [Jr_R(3×5),  zeros(3×5);    ← 右足位置对关节的偏导
%       zeros(3×5), Jr_L(3×5);     ← 左足位置对关节的偏导
%       Jo_R(3×5),  zeros(3×5);    ← 右足姿态对关节的偏导
%       zeros(3×5), Jo_L(3×5)]     ← 左足姿态对关节的偏导
% 注意: 每只脚的位姿只取决于同侧腿的5个关节
Jc = Contact_Jacobian(q(1), q(2), q(3), q(4), q(5), ...
                       q(6), q(7), q(8), q(9), q(10), ...
                       RR(1), RR(2), RR(3), RR(4), RR(5), ...
                       RR(6), RR(7), RR(8), RR(9));

% 构建力-力矩到关节扭矩的映射矩阵 (10×12)
% contact_mapping = Jc' = 雅可比矩阵的转置
% τ = contact_mapping * [F_R; F_L; M_R; M_L]
%   = [Jr_R'*F_R + Jo_R'*M_R;        ← 右腿5个关节力矩
%      Jr_L'*F_L + Jo_L'*M_L]        ← 左腿5个关节力矩
contact_mapping = [blkdiag(Jc(1:3, :)', Jc(4:6, :)'), ...
                   blkdiag(Jc(7:9, :)', Jc(10:12, :)')];

%% ===== 第10步: 构建状态空间动力学模型 (SRBD) =====
% 动力学在 世界坐标系 (World Frame) 中描述
%
% 连续时间状态方程:  dx/dt = Ac·x + Bc·u
%
% 状态 x (13维) = [φ, θ, ψ,      ← 欧拉角 roll, pitch, yaw
%                   px, py, pz,   ← CoM位置 (世界系)
%                   ωx, ωy, ωz,   ← 角速度 (世界系)
%                   vx, vy, vz,   ← 线速度 (世界系)
%                   g]            ← 重力加速度 (常数)
%
% 控制 u (12维) = [F_Rx, F_Ry, F_Rz,    ← 右脚地面反力
%                   F_Lx, F_Ly, F_Lz,    ← 左脚地面反力
%                   M_Rx, M_Ry, M_Rz,    ← 右脚地面反力矩
%                   M_Lx, M_Ly, M_Lz]    ← 左脚地面反力矩

% 获取足端在世界坐标系中的旋转矩阵
% R_foot = [R_foot_R(3×3); R_foot_L(3×3)], 共6×3
R_foot = Rotm_foot(q(1), q(2), q(3), q(4), q(5), ...
                    q(6), q(7), q(8), q(9), q(10), ...
                    RR(1), RR(2), RR(3), RR(4), RR(5), ...
                    RR(6), RR(7), RR(8), RR(9));
R_foot_R = R_foot(1:3, :);  % 右脚在世界系中的旋转矩阵
R_foot_L = R_foot(4:6, :);  % 左脚在世界系中的旋转矩阵

% 转动惯量的世界坐标系变换
% I_world = R · I_body · Rᵀ  (相似变换)
I = R * Ib * R';  % 世界系中的转动惯量张量 (3×3)

% 初始化离散化矩阵的cell数组
B     = repmat({zeros(13, 12)}, h, 1);  % h个离散B矩阵
A_hat = repmat({zeros(13, 13)}, h, 1);  % h个离散A矩阵

%% ===== 第10.1步: 构建连续/离散系统矩阵 =====
for i = 1:h
    %% --- 连续A矩阵 Ac (13×13) ---
    % 欧拉角运动学: dΘ/dt = E^{-1}(Θ) · ω
    % 对于 ZYX 欧拉角约定, 角速度 ω 与欧拉角速率 Θ̇ 的关系为:
    %   ω = E(Θ) · Θ̇
    %   E = [ cos(ψ)cos(θ),  -sin(ψ),  0 ;
    %         sin(ψ)cos(θ),   cos(ψ),  0 ;
    %        -sin(θ),         0,       1 ]
    % 因此: Θ̇ = E^{-1} · ω  (E\eye(3) = inv(E))
    %
    % Ac 结构:
    % ┌─────────────┬─────────────┬─────────────┬─────────────┬──────┐
    % │  0_{3×3}    │  0_{3×3}    │  E^{-1}     │  0_{3×3}    │ 0    │ ← dΘ/dt = E⁻¹·ω
    % │  0_{3×3}    │  0_{3×3}    │  0_{3×3}    │  I_{3×3}    │ 0    │ ← dP/dt = v
    % │  0_{3×3}    │  0_{3×3}    │  0_{3×3}    │  0_{3×3}    │ 0    │ ← dω/dt = ... (见Bc)
    % │  0_{3×3}    │  0_{3×3}    │  0_{3×3}    │  0_{3×3}    │[0;0;-1]│ ← dv/dt = -g
    % │  0_{1×3}    │  0_{1×3}    │  0_{1×3}    │  0_{1×3}    │ 0    │ ← dg/dt = 0
    % └─────────────┴─────────────┴─────────────┴─────────────┴──────┘
    %
    % 关键元素说明:
    % - Ac(1:3, 7:9) = E^{-1}:  角速度 → 欧拉角速率 (运动学)
    % - Ac(4:6, 10:12) = eye(3): 线速度 → 位置变化 (运动学)
    % - Ac(10:12, 13) = [0;0;-1]: 重力贡献 dv_z/dt = -g
    % - 角加速度(dω/dt)和线加速度(dv/dt)的动力学部分由 Bc·u 提供

    Ac = [zeros(3,3), zeros(3,3), ...
          [cos(x_traj(3,i))*cos(x_traj(2,i)), -sin(x_traj(3,i)), 0;
           sin(x_traj(3,i))*cos(x_traj(2,i)),  cos(x_traj(3,i)), 0;
          -sin(x_traj(2,i)),                   0,                 1] \ eye(3), ...
          zeros(3,3), zeros(3,1);
          zeros(3,3), zeros(3,3), zeros(3,3), eye(3), zeros(3,1);
          zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,1);
          zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3), [0; 0; -1];
          zeros(1,13)];

    %% --- 连续B矩阵 Bc (13×12) ---
    % Bc 结构 (控制量 u 对状态导数的影响):
    % ┌───────────────────────────────────────────────────────────┐
    % │ 0_{3×3}      0_{3×3}       0_{3×3}       0_{3×3}       │ ← Θ̇不受力直接影响
    % │ 0_{3×3}      0_{3×3}       0_{3×3}       0_{3×3}       │ ← Ṗ不受力直接影响
    % │ I⁻¹·skew(rR) I⁻¹·skew(rL)  I⁻¹·I_{3×3}   I⁻¹·I_{3×3}  │ ← ω̇ = I⁻¹(r×F+M)
    % │ I_{3×3}/m    I_{3×3}/m     0_{3×3}       0_{3×3}       │ ← V̇ = ΣF/m
    % │ 0_{1×3}      0_{1×3}       0_{1×3}       0_{1×3}       │ ← ġ = 0
    % └───────────────────────────────────────────────────────────┘
    %
    % 牛顿-欧拉动力学:
    %
    % 角动量方程 (绕CoM):
    %   I·ω̇ = r_R × F_R + r_L × F_L + M_R + M_L
    %   ω̇   = I⁻¹·[skew(r_R)·F_R + skew(r_L)·F_L + M_R + M_L]
    %
    % 线动量方程:
    %   m·v̇ = F_R + F_L - m·g·e_z  (重力部分在Ac中)
    %   v̇   = (F_R + F_L)/m
    %
    % 各分量说明:
    % - Bc(7:9, 1:3):   I\skew(r_R), 右脚力产生的角加速度
    % - Bc(7:9, 4:6):   I\skew(r_L), 左脚力产生的角加速度
    % - Bc(7:9, 7:9):   I\eye(3) = I⁻¹, 右脚力矩产生的角加速度
    % - Bc(7:9, 10:12): I\eye(3) = I⁻¹, 左脚力矩产生的角加速度
    % - Bc(10:12, 1:3):  eye(3)/m, 右脚力产生的线加速度
    % - Bc(10:12, 4:6):  eye(3)/m, 左脚力产生的线加速度
    %
    % r_R = foot_traj(1:3) - x_traj(4:6) = 右脚位置 - CoM位置 (世界系)
    % r_L = foot_traj(4:6) - x_traj(4:6) = 左脚位置 - CoM位置 (世界系)

    Bc = [zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3);
          zeros(3,3), zeros(3,3), zeros(3,3), zeros(3,3);
          I\skew(-x_traj(4:6,i) + foot_traj(1:3,1)), ...  % I⁻¹·skew(r_R)
          I\skew(-x_traj(4:6,i) + foot_traj(4:6,1)), ...  % I⁻¹·skew(r_L)
          I\eye(3), I\eye(3);                             % I⁻¹ (力矩直接贡献)
          eye(3)/m, eye(3)/m, zeros(3), zeros(3);        % (1/m)·I
          zeros(1,12)];

    %% --- 离散化 (前向欧拉法: Forward Euler) ---
    % x_{k+1} = x_k + dt · dx/dt
    %         = x_k + dt · (Ac·x_k + Bc·u_k)
    %         = (I + Ac·dt)·x_k + (Bc·dt)·u_k
    %         = A_hat_k · x_k + B_k · u_k
    %
    % 其中 dt 从 dt_MPC_vec 中读取, 索引为 i+k-1
    % 当前全部固定为 0.04s, 但框架支持可变步长
    B{i}     = Bc * dt_MPC_vec(i + k - 1);         % 离散B矩阵 (13×12)
    A_hat{i} = eye(13) + Ac * dt_MPC_vec(i + k - 1); % 离散A矩阵 (13×13)
end

%% ===== 第11步: 稠密QP公式转换 (Condensed QP Formulation) =====
% 参考文献: Jerez et al., "A condensed and sparse QP formulation for
%            predictive control", CDC-ECC 2011.
%
% 基本原理: 消去状态变量, 将MPC优化问题转化为仅关于控制量的QP
%
% 预测模型 (h步):
%   X = Aqp · x_0 + Bqp · U
%
% 其中:
%   X = [x₁; x₂; ...; x_h]    (13h × 1)  预测状态序列
%   U = [u₀; u₁; ...; u_{h-1}] (12h × 1)  控制输入序列
%
%   Aqp = [A₁;                          (13h × 13) 状态转移矩阵
%          A₂·A₁;
%          A₃·A₂·A₁;
%          ...
%          A_h·...·A₁]
%
%   Bqp = [B₁        0      ...  0  ;   (13h × 12h) 控制→状态映射
%          A₂·B₁    B₂     ...  0  ;
%          A₃·A₂·B₁ A₃·B₂ ...  0  ;
%          ...       ...   ...  ...;
%          (A_h·...·A₂)·B₁ ...  B_h]

%% --- 第11.1步: 构建参考状态向量 y ---
% 将所有h步的参考状态展平为一个长向量
y = reshape(x_traj, [13*h, 1]);  % (130 × 1)

%% --- 第11.2步: 构建 Aqp (状态自由响应矩阵) ---
% Aqp{i} = A₁·A₂·...·A_i (从初始状态 → 第i步状态的传递矩阵)
Aqp = repmat({zeros(13, 13)}, h, 1);
Aqp{1} = A_hat{1};                    % A₁
for i = 2:h
    Aqp{i} = Aqp{i-1} * A_hat{i};     % A_{i-1}·...·A₁ · A_i
end
Aqp = cell2mat(Aqp);  % 将cell数组拼接为 (130 × 13) 的大矩阵

%% --- 第11.3步: 构建 Bqp (控制→状态强迫响应矩阵) ---
% Bqp是 分块下三角矩阵 (block lower triangular):
% Bqp(i,j) = 第j步控制量对第i步状态的影响
%
% 注意: 代码中使用 A_hat{i}^(i-j) 近似替代正确的 A_i·A_{i-1}·...·A_{j+1}
% 这是一个 时不变近似 (time-invariant approximation):
%   - 假设所有中间A矩阵相等 (当dt固定且姿态变化缓慢时近似成立)
%   - 在MPC的滚动时域框架下, 每步重新计算会修正近似误差
Bqp = repmat({zeros(13, 12)}, h, h);

for i = 1:h
    Bqp{i, i} = B{i};  % 对角块: 第i步控制对第i步状态的直接影响
    for j = 1:h-1
        % 下三角块: 第j步控制量经过(i-j)步传播对第i步状态的影响
        % A_hat{i}^(i-j): A矩阵的(i-j)次幂 (时不变近似)
        Bqp{i, j} = A_hat{i}^(i-j) * B{j};
    end
end

% 上三角块置零: 未来的控制不影响过去的状态 (因果性)
for i = 1:h-1
    for j = i+1:h
        Bqp{i, j} = zeros(13, 12);
    end
end
Bqp = cell2mat(Bqp);  % 将cell数组拼接为 (130 × 120) 的大矩阵

%% ===== 第12步: QP目标函数权重 =====
% MPC代价函数:
%   J = Σ_{i=0}^{h-1} [ ||x_i - x_{des,i}||²_Q + ||u_i||²_R ]
%
% 其中:
%   Q = diag(L1): 状态跟踪权重 (惩罚偏离期望状态)
%   R = diag(alpha): 控制量惩罚权重 (惩罚过大的控制输入)

% --- 状态跟踪权重 L1 (13维) ---
% L1 = [roll, pitch, yaw,  px, py, pz,  ωx, ωy, ωz,  vx, vy, vz,  g]
L1 = [850, 600, 250, ...   % 欧拉角权重 (roll最大→防侧翻, yaw最小→允许转向)
      300, 300, 350, ...   % CoM位置权重 (高度z最重要→维持站立高度)
      1, 1, 1, ...         % 角速度权重 (平滑旋转)
      1, 1, 1, ...         % 线速度权重 (速度跟踪)
      0];                  % 重力权重 (重力为常数, 不跟踪)

% --- 控制量惩罚权重 alpha (12维) ---
% alpha = [F_R(3), F_L(3), M_R(3), M_L(3)] × 1e-6
% 非常小的惩罚: 避免奇异, 但不显著影响跟踪性能
alpha = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1] * 1e-6;  % humanoid

% 将权重沿时域复制h次, 构建块对角权重矩阵
L10 = repmat(L1, 1, h);
L   = diag(L10);           % Q (130×130) — 状态权重矩阵
alpha10 = repmat(alpha, 1, h);
K   = diag(alpha10);       % R (120×120) — 控制权重矩阵

%% --- 第12.1步: QP目标函数数学推导 ---
% 代价函数展开:
%   J = (X - Y)ᵀ·Q·(X - Y) + Uᵀ·R·U
%     = (Aqp·x₀ + Bqp·U - Y)ᵀ·Q·(Aqp·x₀ + Bqp·U - Y) + Uᵀ·R·U
%
% 按U整理:
%   J = Uᵀ·(Bqpᵀ·Q·Bqp + R)·U + 2·(Aqp·x₀ - Y)ᵀ·Q·Bqp·U + const
%
% 其中:
%   H = 2·(Bqpᵀ·Q·Bqp + R)   ← Hessian矩阵 (二次项)
%   f = 2·Bqpᵀ·Q·(Aqp·x₀ - Y) ← 梯度向量 (一次项)
%   const = (Aqp·x₀ - Y)ᵀ·Q·(Aqp·x₀ - Y)  ← 常数项 (不影响优化)

Hd = 2 * (Bqp' * L * Bqp + K);  % H (120×120) — QP的Hessian矩阵
fd = 2 * Bqp' * L * (Aqp * x - y);  % f (120×1) — QP的梯度向量

%% ===== 第13步: MPC约束构建 =====
% qpOASES约束格式: lbA ≤ A·u ≤ ubA
% 共340个约束 (每步34个 × h=10步)
%
% 约束层次:
%   行 1-80:   摩擦锥约束 (8×h) — 防止足底打滑
%   行 81-100:  力幅值约束 (2×h) — 法向力不超限
%   行 101-220: 线足约束  (12×h) — CoP在支撑域内
%   行 221-320: 关节力矩  (10×h) — 力矩不超电机限幅
%   行 321-340: Mx=0约束  (2×h) — 踝关节无力矩

A     = DM(numCons, numVar);  % CasADi稀疏矩阵 (340×120)
bigNum  = 1e6;                 % 大数 (近似∞, 用于无上界约束)
smallNum = -bigNum;            % 小数 (近似-∞, 用于无下界约束)
lba   = zeros(numCons, 1);
uba   = lba;

%% --- 约束1: 摩擦锥约束 (Friction Cone) ---
% 物理原理: 为防止足底打滑, 水平力不得超过法向力乘以摩擦系数
%   |F_x| ≤ μ·F_z,  |F_y| ≤ μ·F_z
%
% 每个足的4个线性化约束 (金字塔近似):
%   F_x - μ·F_z ≤ 0        (向前)
%  -F_x - μ·F_z ≤ 0        (向后)
%   F_y - μ·F_z ≤ 0        (向左)
%  -F_y - μ·F_z ≤ 0        (向右)
%
% A_mu 结构 (8×12, 单步):
%   每行对应一个约束, F_R的4个 + F_L的4个  = 8个约束/步
A_mu = [1, 0, -mu, zeros(1,9);        ...  % F_Rx - μ·F_Rz ≤ 0
        0, 1, -mu, zeros(1,9);        ...  % F_Ry - μ·F_Rz ≤ 0
        1, 0,  mu, zeros(1,9);        ...  % -F_Rx - μ·F_Rz ≤ 0 → F_Rx + μ·F_Rz ≥ 0
        0, 1,  mu, zeros(1,9);        ...  % -F_Ry - μ·F_Rz ≤ 0
        zeros(1,3), 1, 0, -mu, zeros(1,6); ...  % F_Lx - μ·F_Lz ≤ 0
        zeros(1,3), 0, 1, -mu, zeros(1,6); ...  % F_Ly - μ·F_Lz ≤ 0
        zeros(1,3), 1, 0,  mu, zeros(1,6); ...  % -F_Lx - μ·F_Lz ≤ 0
        zeros(1,3), 0, 1,  mu, zeros(1,6)];     % -F_Ly - μ·F_Lz ≤ 0

% 沿时域复制: 每个预测步独立约束
A_mu_h = blkdiag(A_mu, A_mu, A_mu, A_mu, A_mu, ...
                   A_mu, A_mu, A_mu, A_mu, A_mu);  % (80×120)

% 摩擦锥约束的上下界:
%   前4个约束(每足):   -∞ ≤ F_xy - μF_z ≤ 0    (水平力不超过摩擦)
%   后4个约束(每足):   0  ≤ F_xy + μF_z ≤ +∞   (负向水平力约束)
lba_mu = repmat([smallNum; smallNum; 0; 0; ...     % 右脚不等式前半
                  smallNum; smallNum; 0; 0], ...    % 左脚不等式前半
                  h, 1);
uba_mu = repmat([0; 0; bigNum; bigNum; ...          % 右脚不等式前半上界
                  0; 0; bigNum; bigNum], ...         % 左脚不等式前半上界
                  h, 1);

%% --- 约束2: 法向力幅值限制 (Force Limit) ---
% 地面反力法向分量约束: F_min ≤ F_z ≤ F_max
% F_min = 0:  地面不能提供拉力
% F_max = 500: 防止过大的力导致关节力矩超出电机限制
A_f   = [0, 0, 1, zeros(1,9); ...           % 右脚Fz
         zeros(1,3), 0, 0, 1, zeros(1,6)];   % 左脚Fz
A_f_h = blkdiag(A_f, A_f, A_f, A_f, A_f, ...
                 A_f, A_f, A_f, A_f, A_f);   % (20×120)
lba_force = repmat([Fmin; Fmin], h, 1);      % 下界: 0N (不能拉)
uba_force = repmat([Fmax; Fmax], h, 1);      % 上界: 500N

%% --- 约束3: 线足约束 (Line-Foot Constraints) ---
% 参考文献: Ding et al., "Orientation-Aware Model Predictive Control
%            with Footstep Adaptation for Dynamic Humanoid Walking",
%            IEEE-RAS Humanoids 2022.
%
% 核心思想: 将足底矩形支撑区域简化为两条正交的接触线
%   线1 (沿Y轴, 脚宽方向, 长度 2*lh): 约束CoP在X方向不超出支撑域
%   线2 (沿X轴, 脚长方向, 长度 2*lt): 约束CoP在Y方向不超出支撑域
%
% 压力中心 (CoP) 约束:
%   -lh ≤ M_y/F_z ≤ lh   →  约束CoP的前后位置在支撑域内
%   -lt ≤ M_x/F_z ≤ lt   →  约束CoP的左右位置在支撑域内
%
% 线足参数:
lt = 0.09 - 0.01;  % 脚长半轴 = 0.08m (实际脚长约0.16m, 留0.02m余量)
lh = 0.06 - 0.02;  % 脚宽半轴 = 0.04m (实际脚宽约0.08m, 留0.04m余量)

% A_LF1: 线足线性化约束 (4×12)
% 约束 CoP 在足底支撑多边形内 (线性化后)
A_LF1 = [-lh*[0,0,1]*R_foot_R', zeros(1,3),  [0,1,0]*R_foot_R', zeros(1,3);   % 右脚 Y方向
         -lt*[0,0,1]*R_foot_R', zeros(1,3), -[0,1,0]*R_foot_R', zeros(1,3);   % 右脚 -Y方向
          zeros(1,3), -lh*[0,0,1]*R_foot_L', zeros(1,3),  [0,1,0]*R_foot_L';  % 左脚 Y方向
          zeros(1,3), -lt*[0,0,1]*R_foot_L', zeros(1,3), -[0,1,0]*R_foot_L']; % 左脚 -Y方向

% A_LF2: 耦合摩擦的线足约束 (8×12)
A_LF2 = [ [0, lt, -mu*lt]*R_foot_R', zeros(1,3),  [0, -mu, -1]*R_foot_R', zeros(1,3);
           zeros(1,3), [0, lt, -mu*lt]*R_foot_L', zeros(1,3),  [0, -mu, -1]*R_foot_L';
          [0, -lt, -mu*lt]*R_foot_R', zeros(1,3), [0, -mu, -1]*R_foot_R', zeros(1,3);
           zeros(1,3), [0, -lt, -mu*lt]*R_foot_L', zeros(1,3), [0, -mu, -1]*R_foot_L';
          [0, lh, -mu*lh]*R_foot_R', zeros(1,3),  [0, mu, 1]*R_foot_R', zeros(1,3);
           zeros(1,3), [0, lh, -mu*lh]*R_foot_L', zeros(1,3),  [0, mu, 1]*R_foot_L';
          [0, -lh, -mu*lh]*R_foot_R', zeros(1,3), [0, mu, -1]*R_foot_R', zeros(1,3);
           zeros(1,3), [0, -lh, -mu*lh]*R_foot_L', zeros(1,3), [0, mu, -1]*R_foot_L'];

% 合并线足约束: 12个/步
uba_LF = repmat(zeros(12,1), h, 1);          % 上界: 全零
lba_LF = repmat(ones(12,1)*smallNum, h, 1);  % 下界: -∞
A_LF   = [A_LF1; A_LF2];                     % 合并为12×12
A_LF_h = blkdiag(A_LF, A_LF, A_LF, A_LF, A_LF, ...
                   A_LF, A_LF, A_LF, A_LF, A_LF);  % (120×120)

%% --- 约束4: 关节力矩约束 (Joint Torque Limit) ---
% τ = contact_mapping · u
% 约束: τ_min ≤ contact_mapping · u ≤ τ_max
%
% 关节力矩上限 (Nm), 对应10个关节:
% [hip1_R, hip2_R, thigh_R, knee_R, toe_R, hip1_L, hip2_L, thigh_L, knee_L, toe_L]
%  33.5    33.5    50       50      33.5    33.5    33.5    50       50      33.5
A_tau   = contact_mapping;  % (10×12) — 力/力矩 → 关节力矩映射
A_tau_h = blkdiag(A_tau, A_tau, A_tau, A_tau, A_tau, ...
                    A_tau, A_tau, A_tau, A_tau, A_tau);  % (100×120)

%% --- 约束4.1: 步态调度 (Gait Schedule) ---
% 关键: 摆动腿的控制输入必须为零 (摆动腿不接触地面)
% gaitSchedule 生成一个 0/1 掩码向量:
%   1 → 支撑腿 (Control ON: 可以产生力/力矩)
%   0 → 摆动腿 (Control OFF: 力/力矩清零)
%
% 对于 walking (gait=1): R-L-R-L-R / L-R-L-R-L 交替支撑
%   第1-5步: 右脚支撑(Rm=[0;0;0;0;0;1;1;1;1;1])
%   第6-10步: 左脚支撑(Lm=[1;1;1;1;1;0;0;0;0;0])
if gait == 1
    gaitm = gaitSchedule(k, 1);  % 行走模式: 生成步态掩码 (100×1)
elseif gait == 0
    gaitm = ones(10*h, 1);       % 站立模式: 双腿始终支撑
end

% 约束上下界通过掩码实现:
%   支撑腿: -τ_max ≤ τ ≤ τ_max  (正常限幅)
%   摆动腿: 0 ≤ τ ≤ 0          (强制清零)
uba_tau = repmat([33.5; 33.5; 50; 50; 33.5; ...  % 右脚关节力矩上限
                   33.5; 33.5; 50; 50; 33.5], ...  % 左脚关节力矩上限
                   h, 1) .* gaitm;   % 摆动腿 → 上界=0
lba_tau = -uba_tau .* gaitm;         % 摆动腿 → 下界=0

%% --- 约束5: 足部Mx=0约束 (Foot Moment about X-axis) ---
% 踝关节在绕X轴 (roll方向) 不产生力矩
% 这是因为双足机器人的踝关节通常是2自由度 (pitch + yaw), 不包含roll
% M_Rx = [1,0,0] * R_foot_R' * M_R = 0
Moment_selection = [1, 0, 0];  % 选择Mx分量
A_M = [zeros(1,3), zeros(1,3), Moment_selection*R_foot_R', zeros(1,3);  % 右脚Mx
       zeros(1,3), zeros(1,3), zeros(1,3), Moment_selection*R_foot_L']; % 左脚Mx
A_M_h = blkdiag(A_M, A_M, A_M, A_M, A_M, ...
                 A_M, A_M, A_M, A_M, A_M);  % (20×120)
uba_M = zeros(2*h, 1);   % 上界 = 0
lba_M = uba_M;           % 下界 = 0  (等式约束: Mx = 0)

%% --- 第13.1步: 约束汇总 (Constraint Aggregation) ---
% 将所有约束拼接到 A 矩阵和 lba/uba 向量中
% qpOASES格式: lba ≤ A·u ≤ uba

A(1:80,   :)   = A_mu_h;    % 行 1-80:   摩擦锥 (8×10)
A(81:100,  :)   = A_f_h;     % 行 81-100:  力限幅 (2×10)
A(101:220, :)   = A_LF_h;    % 行 101-220: 线足   (12×10)
A(221:320, :)   = A_tau_h;   % 行 221-320: 关节力矩(10×10)
A(321:340, :)   = A_M_h;     % 行 321-340: Mx=0   (2×10)

lba = [lba_mu; lba_force; lba_LF; lba_tau; lba_M];  % (340×1)
uba = [uba_mu; uba_force; uba_LF; uba_tau; uba_M];  % (340×1)

%% ===== 第14步: CasADi + qpOASES 配置与求解 =====
% 将MATLAB矩阵转换为CasADi DM (Dense Matrix) 对象
Hsize = size(Hd);
fsize = size(fd);
H = DM(Hsize);
H = DM(Hd);   % Hessian矩阵 (120×120)
g_qp = DM(fsize);
g_qp = DM(fd);   % 梯度向量 (120×1)

% 构建qpOASES求解器结构
qp = struct;
qp.h = H.sparsity();   % Hessian稀疏模式
qp.a = A.sparsity();   % 约束矩阵稀疏模式

% qpOASES求解器选项:
%   enableEqualities:  启用等式约束处理
%   printLevel: 'low'  仅输出关键信息
%   enableFullLITests: 禁用完整线性无关测试 (加速求解)
opts = struct('enableEqualities', 1, 'printLevel', 'low', ...
             'enableFullLITests', 0);

% 通过CasADi的conic接口创建qpOASES求解器实例
S = conic('S', 'qpoases', qp, opts);

% 调用求解器: 求解稠密QP
%   min  ½·uᵀ·H·u + g_qpᵀ·u
%   s.t. lba ≤ A·u ≤ uba
r = S('h', H, 'g', g_qp, 'a', A, 'lba', lba, 'uba', uba);

%% ===== 第15步: 输出结果 =====
tic
x_opt = r.x;  % QP最优解 (120×1): 未来h步的控制序列

% 显示求解信息
disp(' ');
disp(['MPC Time Step: ', num2str(k)]);
disp('QPOASES-MPC Solve Time:');
toc

% 提取当前控制步的GRFM (仅使用第一步)
% 滚动时域: 只实施第一个控制量 u₀, 下一步重新求解
GRFM = full(x_opt);
GRFM = GRFM(1:12);  % 取前12个元素: [F_R(3), F_L(3), M_R(3), M_L(3)]

% 输出GRFM, 后续由 JacobianMapping.m 转换为关节力矩 τ
% 注释掉的备选方案 (直接在MPC内部完成力矩映射):
%   u = -contact_mapping * GRFM(1:12);
u = [GRFM(1:12)];

disp('QPOASES-MPC Function Total Time:');
toc
end


%% ========================================================================
%%                        内部辅助函数 (Nested Functions)
%% ========================================================================

%% --- Calc_foot_traj_3Dwalking: 3D行走的足端轨迹预测 ---
% 注意: 此函数目前未被默认调用 (因为可能存在3D行走bug)
% 当前默认使用简化的足端假设 (foot under hip)
%
% 功能: 基于当前支撑足位置和期望速度, 预测未来h步的足端位置
%
% 输入:
%   xdes(12) - 期望SRBD状态
%   x_traj(13×h) - 参考状态轨迹
%   foot(6)  - 当前足端位置 [p_L(3); p_R(3)]
%   h, k     - 时域长度和当前索引
%   R        - 身体旋转矩阵
%
% 输出:
%   foot_traj(6×h) - 预测的足端位置 [p_R(1:3,:); p_L(1:3,:)]
%
% 步态支撑序列 (walking):
%   第1-5步 (i_MPC_gait ∈ [1,5]):  右腿支撑, 左腿摆动 → R-L-R
%   第6-10步 (i_MPC_gait ∈ [6,10]): 左腿支撑, 右腿摆动 → L-R-L
function foot_traj = Calc_foot_traj_3Dwalking(xdes, x_traj, foot, h, k, R)
global dt_MPC_vec
i_MPC_gait = rem(k, h);  % 当前在步态周期中的位置 (1~10)

foot_traj = zeros(6, h);

% 情况1: 支撑序列 R → L → R (右腿当前支撑)
if 1 <= i_MPC_gait && i_MPC_gait <= 5   % stance sequence: R-L-R
    current_R = foot(1:3);  % 当前右脚 (锚定足)

    % 预测后续落脚点:
    % 下一步左脚位置 = 当前右脚 + 期望速度 × 半步时长
    next_L = current_R + xdes(10) * dt_MPC_vec(k) * h / 2;
    % 再下一步右脚位置 = 左脚 + 期望速度 × 半步时长
    next_R = next_L + xdes(10) * dt_MPC_vec(k+5) * h / 2;

    % 阶段1: 右脚锚定 (swing_idx = 6-i_MPC_gait 步)
    for i = 1:6 - i_MPC_gait
        foot_traj(4:6, i) = R * current_R;  % 右脚在世界系中的位置
    end
    % 阶段2: 左脚锚定 (固定5步)
    for i = 7 - i_MPC_gait : 11 - i_MPC_gait
        foot_traj(1:3, i) = R * next_L;
    end
    % 阶段3: 下一右脚锚定 (剩余步数)
    if i_MPC_gait > 1
        for i = 12 - i_MPC_gait : h
            foot_traj(4:6, i) = R * next_R;
        end
    end

% 情况2: 支撑序列 L → R → L (左腿当前支撑)
else  % stance sequence: L-R-L
    if i_MPC_gait == 0
        i_MPC_gait = 10;
    end
    i_MPC_gait = i_MPC_gait - h/2;  % 映射到 [1,5] 区间
    current_L = foot(1:3);  % 当前左脚 (锚定足)

    % 注意: n 变量在此函数中未定义 (这是一个已知bug)
    % 正确应为: next_R = current_L + xdes(10) * dt_MPC_vec(k) * h/2;
    next_R = current_L + n * xdes(10) * dt_MPC_vec(k) * h / 2;
    next_L = next_R + n * xdes(10) * dt_MPC_vec(k+5) * h / 2;

    % 阶段1: 左脚锚定
    for i = 1:6 - i_MPC_gait
        foot_traj(1:3, i) = R * current_L;
    end
    % 阶段2: 右脚锚定 (5步)
    for i = 7 - i_MPC_gait : 11 - i_MPC_gait
        foot_traj(4:6, i) = R * next_R;
    end
    % 阶段3: 下一左脚锚定
    if i_MPC_gait > 1
        for i = 12 - i_MPC_gait : h
            foot_traj(1:3, i) = R * next_L;
        end
    end
end
end


%% --- Calc_x_traj: 计算参考状态轨迹 ---
% 功能: 基于当前状态 + 期望速度, 线性外推未来h步的期望状态
%
% 输入:
%   xdes(13) - 期望状态 [eul(3), p(3), ω(3), v(3), g]
%   x(13)    - 当前实际状态
%   h        - 预测时域长度
%   k        - 当前MPC阶段索引
%
% 输出:
%   x_traj(13×h) - 参考状态轨迹
%     行1-3: 期望欧拉角
%     行4-6: 期望CoM位置
%     行7-9: 期望角速度 (恒定)
%     行10-12: 期望线速度 (恒定)
%     行13: 重力 (常数 9.81)
%
% 轨迹生成逻辑:
%   如果期望速度为0 (xdes(6+j)=0): 使用期望位置作为目标 (位置控制)
%   如果期望速度非0:            从当前位置外推 (速度控制)
function x_traj = Calc_x_traj(xdes, x, h, k)
global dt_MPC_vec
for i = 0:h-1
    for j = 1:6  % 遍历 [eul(3), p(3)] 共6个位置/姿态状态
        if xdes(6+j) == 0
            % 速度为零 → 期望位置恒定不变 (如: 站立不动)
            x_traj(j, i+1) = xdes(j) + xdes(6+j) * sum(dt_MPC_vec(k:k+i));
        else
            % 速度非零 → 从当前位置开始外推 (如: 匀速行走)
            x_traj(j, i+1) = x(j) + xdes(6+j) * sum(dt_MPC_vec(k:k+i));
        end
        % 速度/角速度分量: 保持期望速度不变
        x_traj(6+j, i+1) = xdes(6+j);
    end
    % 重力分量: 始终保持为常数 g
    x_traj(13, i+1) = xdes(13);
end
end


%% --- gaitSchedule: 步态调度函数 ---
% 功能: 生成h步的步态掩码, 标记每一步中哪些腿是支撑腿
%
% 输入:
%   i - 当前MPC阶段索引 (i_MPC_var)
%   gait - 步态模式 (当前仅支持 walking)
%
% 输出:
%   gaitm (10*h × 1) - 步态掩码向量
%     每个元素 ∈ {0, 1}
%     1 = 支撑腿, 对应的控制量生效
%     0 = 摆动腿, 对应的控制量清零
%
% 行走步态模式 (10步周期):
%   Rm = [0;0;0;0;0; 1;1;1;1;1]  ← 右脚支撑 (前5步右腿撑)
%   Lm = [1;1;1;1;1; 0;0;0;0;0]  ← 左脚支撑 (后5步左腿撑)
%
%   索引: [R关节1-5, L关节1-5] × h步
%
%   步态周期示意 (10步):
%   ┌──────────────────────────────────────────────┐
%   │ 步1-5 (右腿支撑): R关节有控制, L关节清零    │
%   │ 步6-10(左腿支撑): L关节有控制, R关节清零    │
%   └──────────────────────────────────────────────┘
function gaitm = gaitSchedule(i, gait)
h = 10;
k = rem(i, h);       % 当前步在10步周期中的位置 (1~10)
if k == 0
    k = 10;
end

% 定义支撑掩码模板 (每个关节维度)
% Rm: 右脚支撑时, 哪些关节可以产生力矩
%     [R_hip1, R_hip2, R_thigh, R_knee, R_toe, L_hip1, ...]
%      0        0        0         0        0       1   ...  ← 左腿关节可作用
Rm = [0; 0; 0; 0; 0;  1; 1; 1; 1; 1];  % 右脚支撑 → 左脚关节激活
Lm = [1; 1; 1; 1; 1;  0; 0; 0; 0; 0];  % 左脚支撑 → 右脚关节激活

% 循环移位索引: 根据当前步k旋转步态模板
% j = [2-k, 3-k, 4-k, 5-k, 6-k, 7-k, 8-k, 9-k, 10-k, 11-k]
% 将负索引映射回1~10范围 (MATLAB不支持负索引)
j = [2-k, 3-k, 4-k, 5-k, 6-k, ...
      7-k, 8-k, 9-k, 10-k, 11-k];
j(j <= 0) = j(j <= 0) + h;

% 分配步态: 前5步右脚支撑(Rm), 后5步左脚支撑(Lm)
IOI{j(1)}  = Rm;
IOI{j(2)}  = Rm;
IOI{j(3)}  = Rm;
IOI{j(4)}  = Rm;
IOI{j(5)}  = Rm;
IOI{j(6)}  = Lm;
IOI{j(7)}  = Lm;
IOI{j(8)}  = Lm;
IOI{j(9)}  = Lm;
IOI{j(10)} = Lm;

% 垂直拼接: 将10个10×1向量拼接为 100×1 的完整步态掩码
gaitm = vertcat(IOI{1}, IOI{2}, IOI{3}, IOI{4}, IOI{5}, ...
                 IOI{6}, IOI{7}, IOI{8}, IOI{9}, IOI{10});
end


%% --- skew: 反对称矩阵 (Skew-Symmetric Matrix) ---
% 功能: 将3D向量转换为对应的反对称矩阵, 用于叉积计算
%
% 输入: v = [vx; vy; vz]
% 输出: skew(v) = [ 0   -vz   vy;
%                    vz   0   -vx;
%                   -vy  vx    0  ]
%
% 数学性质: skew(v) · w = v × w (叉积等价于反对称矩阵乘法)
%
% 在MPC中的用途:
%   力矩 τ = r × F → τ = skew(r) · F
%   所以: I⁻¹·skew(r)·F = I⁻¹·(r×F) = 由力F产生的角加速度
function A = skew(v)
A = [0,    -v(3),  v(2);
     v(3),  0,    -v(1);
    -v(2),  v(1),  0];
end
