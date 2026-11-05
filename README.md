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

## 梯度力模型

使用粒子的实动态极化率 `alpha` 将光强转换成偶极势和梯度力：

\[
U=-\frac{\alpha}{2\epsilon_0c}I,
\qquad \mathbf F_{grad}=-\nabla U.
\]

蓝失谐情况下 `alpha<0`，粒子被推向低光强区域。代码中的
`cfg.dynamicPolarizabilitySI=-1.0e-39 C*m^2/V` 只是用于演示量级；计算具体
原子、分子或纳米粒子时，必须替换为该对象在 **390 nm**、对应内态和偏振下的
动态极化率。光强分布不受这个参数影响，但梯度力与它成正比。

理想 HGB 的几何轴上光强严格为零，因此代码既保存真正的轴上结果
`onAxisIntensity`/`onAxisForceZ`，也给出穿过亮环位置的轴向曲线
`brightRingIntensityZ`/`brightRingForceZ`，避免把物理上的零轴线误认为数值错误。

宏观厘米尺度图不能同时解析 195 nm 的驻波周期，因此宏观图显示相位平均包络，
观察面附近则使用独立的细网格显示真实相干干涉和 `Fz`。
