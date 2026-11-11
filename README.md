# 695 nm 高斯光束中心遮挡仿真（MATLAB）

本项目使用 MATLAB 和傅里叶域菲涅耳传递函数，模拟沿 **+z** 方向传播的 695 nm TEM00
高斯光束。默认在 `z = 15 cm` 放置半径 `0.4 mm` 的不透明圆盘，并观察至 `z = 25 cm`
的衍射与自重构。

## 运行方法

在 MATLAB 中将当前文件夹切换到项目目录，然后执行：

```matlab
gaussian_beam_obstruction
```

程序不依赖额外工具箱，结果保存在 `output/`：

- `xy_intensity_planes.png`：z=0、15 cm（遮挡后）、25 cm 的 xy 光强；
- `longitudinal_xz_intensity.png`：y=0 的 xz 传播光强；
- `xyz_intensity_3d.png`：xyz 三维光强点云；
- `intensity_data.mat`：x、y、z 坐标、完整三维光强以及仿真参数，可用于后处理。

## 修改参数

函数接受名称-值参数。默认束腰半径为 1 mm、计算窗口宽 8 mm、横向采样数 512。例如：

```matlab
gaussian_beam_obstruction( ...
    'Waist', 0.8e-3, ...
    'ObstacleRadius', 0.3e-3, ...
    'Samples', 768, ...
    'OutputDirectory', 'my_result');
```

所有距离参数均使用米，完整参数包括 `Wavelength`、`Waist`、`ObstacleZ`、
`ObstacleRadius`、`FinalZ`、`Window`、`Samples`、`ZSamples` 和 `OutputDirectory`。
