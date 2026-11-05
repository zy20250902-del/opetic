# 两束相向 HGB 蓝失谐光晶格仿真

本仓库提供一个不依赖 Optimization Toolbox 的 MATLAB 程序
`simulate_hgb_blue_detuned_lattice.m`。这里的 **HGB** 指空心高斯光束
（Hollow Gaussian Beam），其腰面场振幅定义为

\[
E_m(r,0)=A_m\left(\frac{r^2}{w_0^2}\right)^m
\exp\left(-\frac{r^2}{w_0^2}\right),\qquad m\ge 1.
\]

两束光的束腰分别位于 `z=-5 cm` 和 `z=+5 cm`，传播方向分别为 `+z`
和 `-z`。程序计算：

- 从 `-5 cm` 到 `+5 cm` 的 `xz` 光束传播包络；
- `z=0` 附近可分辨驻波周期的干涉光晶格；
- `z=0` 的 `xy` 光强与蓝失谐光偶极势；
- 在固定总功率和有限参数范围内，对 HGB 阶数、束腰及两束功率分配进行网格优化，使
  `z=0` 暗中心到周围亮环之间的势垒（光势阱深度）最大。

## 运行

```matlab
result = simulate_hgb_blue_detuned_lattice();
```

默认结果写入 `results_hgb/`：

- `hgb_propagation.png`：宏观传播包络和中心附近光晶格；
- `xy_intensity_and_potential.png`：`z=0` 的光强与势阱；
- `optimization.png`：所有候选参数的势垒；
- `optimization_history.csv`：优化历史；
- `optimized_parameters.csv`：最优 HGB 参数；
- `hgb_lattice_result.mat`：完整数值数据。

## 物理假设与参数

默认使用 `760 nm` 激光和 Rb-87 D2 线的二能级远失谐近似。因为激光频率高于
原子共振频率，偶极势为正，原子被低光强区域吸引。势能采用

\[
U=\frac{3\pi c^2\Gamma}{2\omega_0^3\Delta}I.
\]

请在主函数开头修改 `cfg`，尤其是：

- `laserWavelength`、`transitionWavelength`、`linewidthHz`：激光与原子参数；
- `totalPower`：两束光总功率（优化过程中保持不变，避免无界的“增大功率”解）；
- `orderCandidates`、`waistCandidates`、`powerFractionCandidates`：有限优化范围；
- `sourceZ`：两个束腰位置，默认 `[-0.05, +0.05] m`；
- `relativePhase`：两束光在对称中心的相对相位，默认 `0` 对应中心亮环最强。

“最高势阱”必须在固定总功率和给定束腰下限等工程约束中定义；否则只要无限增加功率、
无限缩小束腰或让激光无限接近共振，理论势垒就没有有限最大值。本程序优化的是中心暗区到
最近亮环峰值的横向势垒，而不是原子温度、散射率或三维俘获寿命。HGB 中心轴始终为暗轴，
所以若需要封闭的三维陷阱，还需加入端帽光束或其他轴向约束。靠近共振时二能级远失谐公式
也不再可靠，应改用包含精细结构、偏振和散射的原子模型。
