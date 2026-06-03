# HECTOR-MPC
- 这是我的HECTOR机器人MPC的matlab代码学习
- 包括了开源的：
    - 📂Hector_Simulation：论文对应的库
    - 📂casadi代码库：CAS (符号代数)+AD (自动微分)，底层 C++ 编写、开源 (LGPL 协议)、面向非线性优化 + 最优控制的符号 + 自动微分科学计算库



## 项目解读

### 一、Hector_Simulation/ — 核心仿真工程

这是整个项目的心脏，也是一个独立的 Git 仓库。包含以下子结构：

#### 1.1 Hector_Simulation.slx — Simulink 主模型文件（631 KB）

这是整个机器人的 Simscape Multibody 物理模型，包含：
- 机器人的所有刚体连杆（躯干、大腿、小腿、脚、手臂等）
- 关节定义（髋部、膝盖、脚踝、肩部、肘部）
- 与地面的接触模型（弹性接触力/摩擦力）
- MPC 控制器模块、摆动腿 PD 控制模块、手臂控制模块
- 可以通过注释/取消注释来切换是否包含手臂组件

#### 1.2 Initiate_Simulation.m — 仿真初始化脚本

运行 Simulink 模型前必须先执行的脚本，负责：
- 导入 CasADi 库
- 预计算全局符号函数 Contact_Jacobian 和 Rotm_foot（足部接触雅可比矩阵和旋转矩阵）
- 设置物理世界参数（阻尼、刚度、摩擦系数等）
- 定义 MPC 参数（步态周期 gaitcycle = 0.4s，预测时域 h = 10，MPC 采样时间）
- 设置机器人初始状态（躯干高度 0.55m，各关节初始角度/限位）
- 可选定义不平坦地形（如正弦波网格）
- 生成足部接触点云

#### 1.3 Functions/ — 控制器函数库

这是项目的算法核心，包含以下关键文件：

| 文件 | 功能 |
|------|------|
| `LocomotionMPC_qpoases.m` | 纯运动的力-矩 MPC 控制器。基于 2021 CDC 论文，实现了 13 维状态空间（欧拉角、位置、角速度、线速度、重力）的模型预测控制。将 MPC 转化为稠密 QP 问题，通过 CasADi 调用 qpOASES 求解。约束包含：摩擦锥、力限幅、线足约束（Line Foot）、关节力矩限幅、足部力矩 Mx=0 |
| `LocoManipulationMPC_qpoases.m` | 运动+操作（Loco-manipulation）MPC 控制器。在 LocomotionMPC 基础上增加了：变质量负载处理（Variable_Mass）、负载接触时间调度（payloadContact）、负载导致的姿态补偿。对应 2023 年 arXiv 论文 |
| `Formulate_Contact_Jacobian.m` | 使用符号计算推导足部接触雅可比矩阵 Jc = [dr; do]/dq。一次性预计算，编译为 matlabFunction 供后续调用。输出足部在世界系下的位置/姿态对关节角的导数 |
| `JacobianMapping.m` | 将 MPC 输出的 GRFM（地面反力+力矩）映射为关节扭矩：τ = Jcᵀ · (-GRFM) |
| `Swing_PD.m` | 摆动腿 PD 控制器。根据步态调度决定哪条腿是摆动腿，使用启发式落足策略计算期望落脚点（基于 Raibert 方法），然后通过逆运动学 + PD 控制跟踪期望轨迹 |
| `Arm_Controller.m` | 手臂 PD 控制器。有两种模式：1) 启发式摆臂（与对侧腿同步摆动）；2) 固定姿态。通过 PD 控制器跟踪期望关节角度 |
| `AssignMPCStage.m` | 根据仿真时间 t 确定当前 MPC 阶段索引 i_MPC_var、当前 MPC 步长 dt_MPC、当前步态索引 i_gait（支撑腿/摆动腿判断） |
| `Variable_Mass.m` | 负载质量随时间变化的时间表函数（0.5kg → 4kg → 0.5kg，模拟搬运不同重量物体） |

**Functions/misc/ — 辅助工具函数**

| 文件 | 功能 |
|------|------|
| `Rx.m` / `Ry.m` / `Rz.m` | 绕 X / Y / Z 轴的基本旋转矩阵生成 |
| `foot_to_joint.m` | 足部位置 → 关节角度的逆运动学求解。使用几何方法计算髋部 roll/pitch、膝盖角度，脚踝角度通过补偿保持足底水平 |

#### 1.4 STL files/ — 机器人 3D 模型文件

包含 HECTOR 机器人的所有机械结构 STL/3MF/SLDPRT 文件，用于 Simscape 可视化：

- 腿部：L_calf.STL, R_calf.STL（小腿）、L_thigh.STL, R_thigh.STL（大腿）、L_foot.STL, R_foot.STL（脚部）、L_hip1/2.STL, R_hip1/2.STL（髋部）
- 手臂：upperarm.STL、forearm.STL、shoulder1.STL 及其镜像
- 躯干：body.STL
- 电机：A1motor.STL（Unitree A1 电机）
- 连接件：各类 linkage、plate、bracket 零件
- 杂物：Dorito.stl（可能是测试用物体）、Hector_picture.jpg（机器人照片）

## 控制架构总结

整个系统的控制流程如下：

```
期望状态 x_des ──→ [MPC 控制器] ──→ GRFM (地面反力+力矩)
                        │                    │
                        │                    │
                   qpOASES QP 求解      JacobianMapping
                   (CasADi conic)       (τ = Jcᵀ · u)
                                             │
                                        关节扭矩 τ
                                             │
                              ┌──────────────┼──────────────┐
                              │              │              │
                        [支撑腿]      [摆动腿 PD]      [手臂 PD]
                         MPC扭矩      逆运动学+PD      摆臂/固定
                              │              │              │
                              └──────────────┼──────────────┘
                                             │
                                    Simscape 物理模型
```

**核心技术特征：**

1. **力-矩混合 MPC**：同时优化地面反力（GRF）和力矩（GRM），不单独依赖 ZMP
2. **单刚体动力学（SRBD）**：将机器人简化为单个刚体+无质量腿，大大降低计算量
3. **线足约束**：将足部简化为两条接触线，约束 CoP 必须在支撑区域内
4. **步态调度**：10 步预测时域内 R-L-R-L 交替支撑/摆动