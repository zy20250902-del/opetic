# MATLAB 空心高斯光束（HGB）传播与 87Rb 梯度力

`hgb_rb87_simulation.m` 是一个可直接运行的 MATLAB 脚本，完成以下工作：

1. 生成给定环半径的空心高斯光束，在 `z=0.05 m` 放置薄透镜，并用傅里叶域菲涅尔传递函数计算沿 `z` 轴的传播；
2. 输出 `z = 0.2 m` 的 `xy` 光强图及 `y = 0` 的 `xz` 光强图；
3. 计算 `z = 0.2 m` 平面内 87Rb 原子的径向光势阱和径向梯度力，并绘制光势阱；
4. 计算 `r = 0` 轴线上 87Rb 原子的轴向光势阱和轴向梯度力，并绘制光势阱。

## 使用方法

在 MATLAB 中切换到仓库目录，然后运行：

```matlab
hgb_rb87_simulation
```

脚本会创建 `output/` 目录并生成：

| 文件 | 内容 |
|---|---|
| `HGB_xy_z0p2m.png` | `z=0.2 m` 的 xy 面光强 |
| `HGB_xz.png` | `y=0` 的 xz 面光强 |
| `Rb87_radial_potential_z0p2m.png` | `z=0.2 m` 时随径向位置变化的光势阱 |
| `Rb87_axial_potential_r0.png` | `r=0` 时随 z 变化的光势阱 |
| `Rb87_radial_force_z0p2m.csv` | `z=0.2 m` 时的光强、光势和径向力数据 |
| `Rb87_axial_force_r0.csv` | `r=0` 时的光强、光势和轴向力数据 |

## 模型与参数

初始复振幅采用标准空心高斯形式

```text
E(rho,0) ∝ (rho/w0)^m exp[-rho²/w0²]
```

其中 `m=hollowOrder`。输入参数 `ringRadius` 表示**初始光强最大值所在的环半径**，脚本利用
`w0 = ringRadius*sqrt(2/m)` 自动换算束腰参数。默认功率为 50 mW、环半径为 0.50 mm、
阶数为 2。这些参数均可在脚本开头修改。

传播采用傍轴菲涅尔传递函数

```text
H(fx,fy;z) = exp[-i pi lambda z (fx²+fy²)]
```

薄透镜位于 `zLens=0.05 m`，其相位传递函数为

```text
Tlens(x,y) = exp[-i pi (x²+y²)/(lambda f)]
```

脚本利用初始高斯包络的瑞利长度和 ABCD 复光束参数自动求解焦距 `f`，使透镜后的复光束
参数在 `zTarget=0.2 m` 处满足束腰条件。二次方程一般存在强、弱两个聚焦解，脚本选取短焦距
解以获得明显的聚焦效果。运行结束时会在命令行打印实际使用的焦距；`HGB_xz.png` 中的白色
虚线和点线分别标记透镜位置与目标束腰位置。

梯度力采用 87Rb D2 `5S1/2(F=2) -> 5P3/2(F'=3)` 780.241 nm 循环跃迁的远失谐、
低饱和二能级近似：

```text
Udip = [hbar Gamma²/(8 Delta)] I/Isat
Fgrad = -grad(Udip)
```

默认 `Gamma=2*pi*6.065 MHz`、`Isat=1.669 mW/cm²`、红失谐 `Delta=-2*pi*1 GHz`。
若实验所用波长或跃迁不同，必须相应修改 `lambda`、`Gamma`、`Isat` 和 `detuning`；对于远离
共振的光阱，更严谨的做法是以该波长处的动态极化率替代二能级公式。

两个 CSV 同时以焦耳 (`dipole_potential_J`) 和等效温度微开尔文
(`dipole_potential_uK = Udip/kB`) 给出光势。两幅光势阱曲线使用微开尔文为纵坐标；默认红失谐
使高光强区域的光势为负，零光强处为势能零点。

## 数值计算提示

- `windowSize` 应明显大于整个传播区间内的光斑，以避免 FFT 周期边界造成卷绕误差。
- 增大 `N` 可改善横向采样但会增加运行时间；增大 `Nz` 可改善轴向力的数值微分。
- CSV 使用 SI 单位，并在列名中注明单位。
- 在圆对称轴 `r=0` 处，径向力因对称性为零；轴线上表格报告的是沿传播方向的轴向力。
