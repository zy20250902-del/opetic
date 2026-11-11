%% 空心高斯光束的菲涅尔传播及 87Rb 原子梯度力
% 直接运行本脚本。所有可调参数集中在下一节，结果保存在 output/ 中。
clear; close all; clc;

%% 用户参数
lambda = 780.241e-9;         % 光波长 (m)，87Rb D2 跃迁
power = 50e-3;               % 光束总功率 (W)
ringRadius = 0.50e-3;        % z=0 光强环的半径 r (m)
hollowOrder = 2;             % 空心高斯阶数 m（正整数）
zTarget = 0.20;              % 径向力的计算平面 (m)
zMax = 0.40;                 % x-z 图及轴向力的最大传播距离 (m)
zLens = 0.05;                % 薄透镜位置 (m)

N = 512;                     % 每个横向维度的采样点数
windowSize = 4e-3;           % 计算窗口边长 (m)，兼顾聚焦束腰的横向采样
Nz = 401;                    % z 方向采样点数

% 87Rb D2 5S1/2(F=2)->5P3/2(F'=3) 循环跃迁的二能级近似参数
Gamma = 2*pi*6.065e6;        % 自然线宽（角频率，rad/s）
Isat = 1.669e-3/1e-4;        % 饱和光强：1.669 mW/cm^2 -> W/m^2
detuning = -2*pi*1.0e9;      % 激光失谐 Delta=omega_L-omega_0 (rad/s)

assert(hollowOrder >= 1 && hollowOrder == round(hollowOrder), ...
    'hollowOrder 必须为正整数。');
assert(zTarget >= 0 && zTarget <= zMax, 'zTarget 必须在 [0,zMax] 内。');
assert(zLens > 0 && zLens < zTarget, 'zLens 必须在 0 与 zTarget 之间。');

%% 初始空心高斯光场
% 定义 E(rho,0) ~ (rho/w0)^m exp[-rho^2/w0^2]。
% 对应光强最大值位于 rho=w0*sqrt(m/2)，故按给定环半径反求 w0。
dx = windowSize/N;
x = ((-N/2):(N/2-1))*dx;
[X,Y] = meshgrid(x,x);
rho = hypot(X,Y);
w0 = ringRadius*sqrt(2/hollowOrder);

amplitudeShape = (rho/w0).^hollowOrder .* exp(-(rho/w0).^2);
intensityShape = amplitudeShape.^2;
I0 = power*intensityShape/(sum(intensityShape,'all')*dx^2);
U0 = sqrt(I0);               % 采用 sqrt(W/m^2) 作为归一化复振幅单位

%% 菲涅尔传播，并在 z=0.05 m 放置薄透镜
% 傅里叶域菲涅尔公式：H(fx,fy;z)=exp[-i*pi*lambda*z*(fx^2+fy^2)]。
f = ((-N/2):(N/2-1))/windowSize;
[FX,FY] = meshgrid(f,f);
spatialFrequencySquared = FX.^2 + FY.^2;
spectrum0 = fftshift(fft2(ifftshift(U0)));
propagateBeforeLens = @(z) fftshift(ifft2(ifftshift(spectrum0 .* ...
    exp(-1i*pi*lambda*z*spatialFrequencySquared))));

% 初始高斯包络的 q 参数为 i*zR。根据 ABCD 定律选择透镜焦距，使透镜后的
% q 参数传播至 zTarget 时实部为零，即让新的束腰位于 zTarget。
rayleighRange = pi*w0^2/lambda;
distanceAfterLens = zTarget-zLens;
inverseQAtLens = 1/(zLens+1i*rayleighRange);
focusDiscriminant = 1-4*distanceAfterLens^2*imag(inverseQAtLens)^2;
assert(focusDiscriminant >= 0, ...
    '当前参数下单个薄透镜无法在 zTarget 形成实束腰，请调整束腰或透镜位置。');
% 二次方程有强、弱两个解；选取短焦距解以形成清晰的聚焦束腰。
inverseQAfterLensReal = (-1-sqrt(focusDiscriminant))/(2*distanceAfterLens);
focalLength = 1/(real(inverseQAtLens)-inverseQAfterLensReal);

fieldBeforeLens = propagateBeforeLens(zLens);
lensPhase = exp(-1i*pi*rho.^2/(lambda*focalLength));
fieldAfterLens = fieldBeforeLens.*lensPhase;
spectrumAfterLens = fftshift(fft2(ifftshift(fieldAfterLens)));
propagateAfterLens = @(distance) fftshift(ifft2(ifftshift(spectrumAfterLens .* ...
    exp(-1i*pi*lambda*distance*spatialFrequencySquared))));

z = linspace(0,zMax,Nz).';
xzIntensity = zeros(Nz,N);
centerRow = N/2 + 1;
for iz = 1:Nz
    if z(iz) <= zLens
        Uz = propagateBeforeLens(z(iz));
    else
        Uz = propagateAfterLens(z(iz)-zLens);
    end
    xzIntensity(iz,:) = abs(Uz(centerRow,:)).^2;
end
Utarget = propagateAfterLens(zTarget-zLens);
Itarget = abs(Utarget).^2;

%% 绘制 xy 面与 xz 面光强
outputDir = fullfile(pwd,'output');
if ~exist(outputDir,'dir'), mkdir(outputDir); end

figure('Color','w','Name','HGB xy intensity');
imagesc(x*1e3,x*1e3,Itarget/1e4);
axis image; axis xy; colorbar;
xlabel('x (mm)'); ylabel('y (mm)');
title(sprintf('HGB intensity at z = %.3f m (W/cm^2)',zTarget));
colormap turbo;
exportgraphics(gcf,fullfile(outputDir,'HGB_xy_z0p2m.png'),'Resolution',200);

figure('Color','w','Name','HGB xz intensity');
imagesc(x*1e3,z,xzIntensity/1e4);
axis xy; colorbar;
xlabel('x (mm)'); ylabel('z (m)');
title('HGB Fresnel propagation, y = 0 (W/cm^2)');
colormap turbo;
hold on;
plot(x([1 end])*1e3,[zLens zLens],'w--','LineWidth',1.2);
plot(x([1 end])*1e3,[zTarget zTarget],'w:','LineWidth',1.2);
text(x(8)*1e3,zLens,' lens','Color','w','VerticalAlignment','bottom');
text(x(8)*1e3,zTarget,' waist','Color','w','VerticalAlignment','bottom');
hold off;
exportgraphics(gcf,fullfile(outputDir,'HGB_xz.png'),'Resolution',200);

%% z=0.2 m 时不同径向位置的光势阱与梯度力
% 远失谐、低饱和二能级近似：Udip=(hbar*Gamma^2/(8*Delta))*I/Isat。
hbar = 1.054571817e-34;
kB = 1.380649e-23;           % 玻尔兹曼常数 (J/K)，用于将光势换算为温度
dipoleCoefficient = hbar*Gamma^2/(8*detuning*Isat); % U = coefficient*I

positiveSide = centerRow:N;
rSamples = x(positiveSide).';
radialIntensity = Itarget(centerRow,positiveSide).';
radialPotential = dipoleCoefficient*radialIntensity;
radialPotentialMicroK = radialPotential/kB*1e6;
radialForce = -gradient(radialPotential,rSamples);

radialTable = table(rSamples,radialIntensity,radialPotential, ...
    radialPotentialMicroK,radialForce, ...
    'VariableNames',{'r_m','intensity_W_m2','dipole_potential_J', ...
    'dipole_potential_uK','radial_force_N'});
writetable(radialTable,fullfile(outputDir,'Rb87_radial_force_z0p2m.csv'));

figure('Color','w','Name','87Rb radial optical potential');
plot(rSamples*1e3,radialPotentialMicroK,'LineWidth',1.6);
grid on;
xlabel('r (mm)'); ylabel('U_{dip}/k_B (\muK)');
title(sprintf('^{87}Rb radial optical potential at z = %.3f m',zTarget));
exportgraphics(gcf,fullfile(outputDir, ...
    'Rb87_radial_potential_z0p2m.png'),'Resolution',200);

%% r=0 时不同 z 位置的光势阱与轴向梯度力
axisIntensity = xzIntensity(:,centerRow);
axisPotential = dipoleCoefficient*axisIntensity;
axisPotentialMicroK = axisPotential/kB*1e6;
axialForce = -gradient(axisPotential,z);

axialTable = table(z,axisIntensity,axisPotential,axisPotentialMicroK,axialForce, ...
    'VariableNames',{'z_m','intensity_W_m2','dipole_potential_J', ...
    'dipole_potential_uK','axial_force_N'});
writetable(axialTable,fullfile(outputDir,'Rb87_axial_force_r0.csv'));

figure('Color','w','Name','87Rb axial optical potential');
plot(z,axisPotentialMicroK,'LineWidth',1.6);
grid on;
xlabel('z (m)'); ylabel('U_{dip}/k_B (\muK)');
title('^{87}Rb axial optical potential at r = 0');
exportgraphics(gcf,fullfile(outputDir, ...
    'Rb87_axial_potential_r0.png'),'Resolution',200);

fprintf('计算完成。\n');
fprintf('  输入功率积分: %.6g W\n',sum(I0,'all')*dx^2);
fprintf('  薄透镜: z = %.4g m, f = %.6g m; 目标束腰: z = %.4g m\n', ...
    zLens,focalLength,zTarget);
fprintf('  光强图、光势阱图及两个数据表已写入: %s\n',outputDir);
