# 390 nm 相向传播 HGB 光束：光强与梯度力

`simulate_hgb_390nm.m` 模拟两束波长为 **390 nm** 的空心高斯光束
（Hollow Gaussian Beam, HGB）相向传播：第一束的束腰位于 `z1` 并沿 `+z`
传播，第二束的束腰位于 `z2` 并沿 `-z` 传播。`z1`、`z2` 是函数输入参数，
只需满足 `z1 < z2`。

## 运行示例

```matlab
% 单位为 m：两束光分别从 -5 cm 和 +5 cm 相向传播
result = simulate_hgb_390nm(-0.05, 0.05);

% 修改为 -3 cm 和 +4 cm
result = simulate_hgb_390nm(-0.03, 0.04);
```

不传参数时默认使用 `z1=-0.05 m`、`z2=0.05 m`。`xy` 观察面自动设为
`(z1+z2)/2`。波束阶数、束腰、功率、相对相位和采样范围可在函数开头的
`cfg` 结构中调整。

## 计算内容

程序输出：

1. `z1` 到 `z2` 的宏观 `xz` 光强包络；
2. 观察面附近、能分辨 `lambda/2 = 195 nm` 驻波周期的 `xz` 光强；
3. 观察位置的 `xy` 光强；
4. `xz` 平面中的轴向梯度力 `Fz`；
5. `xy` 平面中的横向梯度力 `Fx`、`Fy` 和力矢量图。

结果保存在 `results_hgb_390nm/`：

- `intensity_maps.png`；
- `gradient_force_maps.png`；
- `axial_profiles.png`；
- `simulation_summary.csv`；
- `hgb_390nm_result.mat`。

## 锶-88 原子的梯度力模型

程序默认计算基态 **Sr-88** 原子所受的光学梯度力。采用
`5s^2 1S0 -> 5s5p 1P1` 强跃迁的二能级近似，默认跃迁波长为
`460.862 nm`、自然线宽为 `Gamma/(2*pi)=32.0 MHz`，并同时保留旋转项和
反旋转项来计算 390 nm 处的动态极化率：

\[
\alpha(\omega)=\frac{3\pi\epsilon_0c^3\Gamma}{\omega_0^3}
\left(\frac{1}{\omega_0-\omega}+\frac{1}{\omega_0+\omega}\right).
\]

随后由动态极化率将光强转换成偶极势和梯度力：

\[
U=-\frac{\alpha}{2\epsilon_0c}I,
\qquad \mathbf F_{grad}=-\nabla U.
\]

390 nm 相对 461 nm 锶强跃迁为蓝失谐，因此该二能级模型给出 `alpha<0`，原子被
推向低光强区域。程序使用 Sr-88 质量进一步给出 `accelerationXZ_Z`、
`accelerationXY_X` 和 `accelerationXY_Y`。最大光学力与最大加速度同时写入
`simulation_summary.csv`。

这里的二能级结果适合用于模型和量级估算。精密实验需要用包含锶原子所有相关跃迁、
偏振、Zeeman 子能级和实际光频率的动态极化率替换该近似；光强计算本身不受原子模型影响。

理想 HGB 的几何轴上光强严格为零，因此代码既保存真正的轴上结果
`onAxisIntensity`/`onAxisForceZ`，也给出穿过亮环位置的轴向曲线
`brightRingIntensityZ`/`brightRingForceZ`，避免把物理上的零轴线误认为数值错误。

宏观厘米尺度图不能同时解析 195 nm 的驻波周期，因此宏观图显示相位平均包络，
观察面附近则使用独立的细网格显示真实相干干涉和 `Fz`。
