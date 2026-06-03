clear;
clc;
%% HECTOR OPEN SOURCE SIMULATION SOFTWARE IN MATLAB/SIMULINK%%
% PLEASE READ LICENSE AGREEMENT BEFORE PROCEEDING

%% Run this script before simulation to load parameters
% 在运行 Simulink 仿真之前必须执行此脚本，以加载所有全局参数

%% Pre-formulate functions (预构建符号函数)
% addpath(genpath('casadi-windows-matlabR2016a-v3.5.1'))
import casadi.*
% generate global symbolic functions (生成全局符号函数):
global Contact_Jacobian Rotm_foot
% 调用 Formulate_Contact_Jacobian 函数，预构建足端接触雅可比矩阵和旋转矩阵
% Contact_Jacobian: 描述关节角速度到足端线速度的映射关系
% Rotm_foot: 足端在世界坐标系下的旋转矩阵
[Contact_Jacobian,Rotm_foot]=Formulate_Contact_Jacobian;

%% General (sim world physics) — 仿真物理环境参数
% 世界坐标系下的平移阻尼，模拟空气阻力等外部阻尼
world_damping = 1e-3;
% 世界坐标系下的旋转阻尼
world_rot_damping = 1e-3;
% 关节刚度 (0 表示关节无弹性恢复力)
joint_stiffness = 0;
% 关节阻尼 (模拟关节摩擦力)
joint_damping = 0.001;
% elastic stiff ground (弹性硬地面接触模型):
% 地面法向接触刚度 (N/m)，高刚度模拟硬地面
contact_stiffness = 1e5;
% 地面法向接触阻尼 (Ns/m)，吸收碰撞能量
contact_damping = 1e3;
% 接触点可视化球体半径 (m)
contact_point_radius = 0.001; % contact cloud visual
% 踝关节被动刚度/阻尼 (此处设为0，踝关节完全由主动控制驱动)
ankle_stiffness = 0.0;
ankle_damping = 0.0;
% elastic joint hard stop (关节限位弹性硬止挡):
% 关节限位刚度 (N/rad)，当关节接近限位时产生弹性恢复力
limit_stiffness = 1e3;
% 关节限位阻尼 (Ns/rad)，吸收限位碰撞能量
limit_damping = 1e2;
% ground friction (static/kinematic) — 地面摩擦参数
% 静摩擦系数
mu_s = 1.0;
% 动摩擦系数
mu_k = 1.0;
% 静/动摩擦转换的临界速度 (m/s)，低于此速度视为静摩擦
mu_vth = 0.001; % critical velocity for transition (m/s)

%% Controller global params — MPC控制器全局参数
% 声明全局变量，供其他脚本/函数共享
global p_foot_w p_foot_pre t_start gaitcycle dt_MPC N gait dt_MPC_vec ...
    acc_t current_height addArm

% N: MPC 总预测步长 (占位值，后续可能被修改)
N = 5000; % place holder for total MPC prediction length
% p_foot_w: 初始足端位姿 [x, y, z, roll, pitch, yaw]^T (世界坐标系)
p_foot_w = zeros(6,1); % initial foot position
% t_start: 仿真起始时间 (s)
t_start = 0; % simulation start time

% MPC setup (MPC 控制参数设置):
% gaitcycle: 步态周期，机器人完成一步的时间 (s)
gaitcycle = 0.4; % gait period
% h: MPC 预测时域长度 (horizon)，即未来预测的步数
h = 10; % horizon length
% dt_MPC: MPC 采样时间 (s)，即每两次优化求解之间的时间间隔
%         gaitcycle/10 = 0.04s，一个步态周期内进行 10 次优化
dt_MPC = gaitcycle/10; % MPC sample time
% p_foot_pre: 足端轨迹预存储矩阵 (6行 × N+10列)
%             存储未来 N 步的足端参考轨迹
p_foot_pre = zeros(6,N+10); % foot matrix initialization

% Gait setup (步态模式选择):
% gait = 0: 站立模式 (Standing)，机器人保持原地不动
% gait = 1: 行走模式 (Walking)，机器人按步态周期行走
gait = 1; % Standing = 0; walking = 1;

%% grid surface definition (if using uneven terrain) — 地形网格定义
% 注意: 使用非平坦地形时，需要在 Hector Simulink 模型中将接触面切换为
%        "Mesh Grid Surface" 模式
% make sure to switch to "Mesh Grid Surface" in Hector model
% 以下为正弦波地形示例 (sine wave example)
% X 方向网格点 (范围 -1m 到 5m，分辨率 0.01m)
x_grid = -1:0.01:5;
% Y 方向网格点 (范围 -1m 到 1m，分辨率 0.01m)
y_grid = -1:0.01:1;
% 生成 XY 网格矩阵
[X,Y] = meshgrid(x_grid,y_grid);
% 地面高度按余弦函数变化: z = 0.05*cos(2πX) - 0.05
% 振幅 ±0.05m (5cm)，形成沿 X 方向的波纹地形
z_heights = [0.05*cos(2*pi*X)-0.05]';

%% contact cloud (接触点云定义):
% 足底接触面用离散点云近似，而非连续曲面
% foot contact approximated by point clouds:
% npt: 每条线上点的数量
npt = 5; % number of contact points per line
% 接触点云由3条线组成，每条线5个点，共15个点:
%   第1行: 足底中心线 (y=0)，X方向从 -0.06m 到 0.09m
%   第2行: 足底左侧线 (y=0.01m)，即偏左 1cm
%   第3行: 足底右侧线 (y=-0.01m)，即偏右 1cm
% X 范围 -0.06~0.09m 表示脚长约 15cm
contact_cloud = [ [[linspace(-0.06,0.09,npt)]',ones(npt,1)*0, ones(npt,1)*0];
                    [[linspace(-0.06,0.09,npt)]',ones(npt,1)*0.01, ones(npt,1)*0];
                    [[linspace(-0.06,0.09,npt)]',ones(npt,1)*-0.01, ones(npt,1)*0]];
% 地面标记网格 (用于可视化地面参考点):
mark_distance=1;            % 标记点间距 (m)
mark_x=-5:mark_distance:5;  % X 方向标记坐标
mark_y=-5:mark_distance:5;  % Y 方向标记坐标
% 构建标记云网格点的 X 坐标矩阵 (重复排列)
mark_cloud_x=repmat(mark_x,length(mark_y),1);
mark_cloud_x=reshape(mark_cloud_x,[length(mark_x)*length(mark_y),1]);
% 构建标记云网格点的 Y 坐标矩阵 (转置重复)
mark_cloud_y=repmat(mark_y',length(mark_x),1);
% 合并为标记云矩阵 [X, Y, Z=0] (所有标记点在地面 z=0 平面上)
mark_cloud=[mark_cloud_x,mark_cloud_y,zeros(length(mark_x)*length(mark_y),1)];

%% vairable MPC dt (可变MPC时间步长 — 功能即将实现)
% 当前版本使用固定步长，可变步长功能将在后续版本中完善
% feature will be implemented soon
% 参考函数 AssignMPCStage(t)
% See function AssignMPCStage(t) %
global fixed_MPC_dt
% fixed_MPC_dt: 1 = 使用固定 MPC 时间步长; 0 = 可变步长
fixed_MPC_dt = 1; % fixed dt = 1;
% dt_MPC_vec: 每个MPC阶段的时长向量，调用 define_dt_MPC 函数生成
dt_MPC_vec = define_dt_MPC(dt_MPC); % vector defining MPC dts
% acc_t: 累积时间数组，记录每个 MPC 阶段对应的仿真时间
%        acc_t(k) = 前 k-1 个阶段的总时长
acc_t(:,1) = zeros(length(dt_MPC_vec)+1,1);
for ith = 2:length(dt_MPC_vec)+1
    % 计算累积时间: 第 ith 个阶段起点的时刻
    acc_t(ith,1) = sum(dt_MPC_vec(1:ith-1));
end

%% Initial Condition/Parameter of Robot State/Joints — 机器人初始状态与关节参数
% --- 身体 (Body) 初始状态 ---
% 身体初始 X/Y 水平位置 (m)
body_x0 = 0;
body_y0 = 0;
% 身体初始高度 (m)，0.55m 为机器人站立时的基准高度
body_z0 = 0.55; % m
% current_height: 当前身体高度 (全局变量，MPC 中动态更新)
current_height = body_z0;
% 身体初始线速度 [vx, vy, vz] (m/s)，静止启动
body_v0 = [0 0 0];
% 身体初始姿态角 [roll, pitch, yaw] (度)
body_R = [0 0 0]; % deg.

% --- 腿部关节初始状态 ---
% hip1: 髋关节 Rz (绕Z轴旋转 — 腿部内外旋转)
%       控制腿部的内外八字角度
hip1_q0 = 0;        % 初始关节角度 (度)
hip1_dq0 = 0;       % 初始关节角速度 (度/秒)
hip_max = 30;        % 关节上限 (度)
hip_min = -30;       % 关节下限 (度)

% hip2: 髋关节 Rx (绕X轴旋转 — 腿部侧向摆动)
%       控制腿部的左右展开/收拢
%       左右腿分别设置 (L=左腿, R=右腿)
hip2_q0_L = 0;       % 左髋侧摆初始角度
hip2_dq0_L = 0;      % 左髋侧摆初始角速度
hip2_q0_R = 0;       % 右髋侧摆初始角度
hip2_dq0_R = 0;      % 右髋侧摆初始角速度
hip2_max = 18;        % 髋侧摆上限 (度)
hip2_min = -18;       % 髋侧摆下限 (度)

% thigh: 大腿关节 Ry (绕Y轴旋转 — 腿部前后摆动)
%        控制大腿的前后运动 (俯仰方向)
thigh_q0 = 45*1;     % 初始角度: 前倾 45° (度)
thigh_dq0 = 0;       % 初始角速度
thigh_max = 120;      % 大腿前摆上限 (度)
thigh_min = -120;     % 大腿后摆下限 (度)

% calf: 小腿/膝关节 Ry (绕Y轴旋转 — 膝关节屈伸)
%       控制小腿的弯曲与伸展
calf_q0 = -90*1;     % 初始角度: 后收 90° (度)，形成弯曲站立姿态
calf_dq0 = 0;        % 初始角速度
knee_max = -15;       % 膝关节伸展上限 (度)，-15°表示微微弯曲
knee_min = -160;      % 膝关节弯曲下限 (度)

% toe/ankle: 脚踝/脚趾关节 Ry (绕Y轴旋转 — 踝关节)
%            控制足端相对于小腿的角度
toe_q0 = 45*1;       % 初始角度: 前倾 45° (度)
toe_dq0 = 0;         % 初始角速度
ankle_max = 75;       % 踝关节上限 (度)
ankle_min = -75;      % 踝关节下限 (度)

% --- 手臂关节参数 ---
% arms (手臂模块):
% addArm = 1: 启用手臂模块; addArm = 0: 禁用手臂模块
addArm = 1;

% shoulder Roll: 肩部 Rx 关节 (绕X轴 — 手臂侧向抬举)
shoulderx_q0 = 0;    % 初始角度
shoulderx_min = -15;  % 下限 (手臂下垂内收)
shoulderx_max = 90;   % 上限 (手臂侧举)

% shoulder Pitch: 肩部 Ry 关节 (绕Y轴 — 手臂前后摆动)
shouldery_q0 = 0;    % 初始角度
shouldery_min = -90;  % 下限 (手臂后摆)
shouldery_max = 90;   % 上限 (手臂前举)

% Elbow: 肘关节 Ry (绕Y轴 — 肘部屈伸)
elbow_q0 = -90;      % 初始角度: 弯曲 90°
elbow_min = -150;     % 下限 (最大弯曲)
elbow_max = -10;      % 上限 (接近伸直)

%% Fncs — 自定义辅助函数
% define_dt_MPC: 生成 MPC 每个阶段的时间步长向量
%   输入 vec: 基础步长向量 (当前仅包含 dt_MPC)
%   输出 out: 每个 MPC 阶段对应的步长列向量
%   逻辑说明:
%     1. 将每个基础步长重复 5 次 (形成 5 个相同步长的子阶段)
%     2. 如果 fixed_MPC_dt 为真，用固定 dt_MPC 填充剩余步长至约 5000 步
function out = define_dt_MPC(vec)
global dt_MPC fixed_MPC_dt
    out = []; % 初始化输出向量为空
    % 第一步: 将输入的基础步长每个重复5次
    for i=1:length(vec)
        out = [out;vec(i)*ones(5,1)];
    end
    % 第二步: 若使用固定步长模式，用 dt_MPC 填充剩余步数
    %         5*(1000-i) ≈ 5000 个总时间步
    if boolean(fixed_MPC_dt)
        out = [out;dt_MPC*ones(5*(1000-i),1)]; % 剩余步数统一使用固定 dt_MPC
    end
end
