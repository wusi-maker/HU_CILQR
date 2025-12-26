(4) 详细的公式推导过程

### 2.1 离散动力学与代价目标（Discrete Dynamics and Cost objectives）
- 状态扰动动力学近似（一阶泰勒展开）
  - 将一个非线性系统化通过泰勒展开作对应的线性系统 
  $$x_{k+1}=f(x_k,u_k)$$
  $$x_{k+1}+\delta x_{k+1}=f\left(x_{k}+\delta x_{k}, u_{k}+\delta u_{k}\right)\approx f\left(x_{k}, u_{k}\right)+\left.\frac{\partial f}{\partial x}\right|_{x_{k}, u_{k}}\left(x-x_{k}\right)+\left.\frac{\partial f}{\partial u}\right|_{x_{k}, u_{k}}\left(u-u_{k}\right)$$
$$\delta x_{k+1} =A(x_k,u_k)\delta x_k+B(x_k,u_k)\delta u_k$$
  - 其中iLQR\CiLQR\DDP的最终目标都是得到一个最优的输入使得对应的代价函数最小，即存在一个对应的 $U=\{u_0,u_1,\dots,u_{m-1}\}$能够使得  $J(x_0,U)=\mathcal{L}(x_N)+\sum_{k=0}^{N-1}\mathcal{x_k,u_k}$最小
### 2.2 贝尔曼最优性原理（The Bellman Optimality Principle）
- 公式3：终端代价函数
$$V_{N}\left(x_{N}\right)=\mathcal{L}\left(x_{N}\right)$$
- 公式4：递推代价函数（最优性方程）
$$
V_{k}\left(x_{k}\right)=\min _{u}\left\{\mathcal{L}\left(x_{k}, u_{k}\right)+V_{k+1}\left(f\left(x_{k}, u_{k}\right)\right)\right\}=\min _{u}\left\{\mathcal{Q}(x_k,u_k)\right\}$$       
令 $\mathcal{Q}(x_k,u_k)=\mathcal{L}\left(x_{k}, u_{k}\right)+V_{k+1}\left(f\left(x_{k}, u_{k}\right)\right)$
故$\delta V=\min _{\delta u}\{\delta Q(x, u)\}$
- 公式5：代价函数扰动形式
$$V_{k}+\delta V_{k}=V_{k}(x_{k}+\delta x_{k})\approx V(x_{k})+\left.\frac{\partial V}{\partial x}\right|_{x_{k}}\left( x-x_{k}\right) +\left.\frac{1}{2}\left( x-x_{k}\right) ^{T}\frac{\partial ^{2}V}{\partial x^{2}}\right| _{x_{k}}\left( x-x_{k}\right)$$

### 2.3 反向传播（Backward pass）
- 公式6：Q函数扰动展开（全量形式）
$$\begin{aligned} Q_{k}+\delta Q_{k}= & Q\left(x_{k}+\delta x, u_{k}+\delta u\right) \\ \approx & Q\left(x_{k}, u_{k}\right)+\left.\frac{\partial Q}{\partial x}\right|_{x_{k}, u_{k}}\left(x-x_{k}\right)+\left.\frac{\partial Q}{\partial u}\right|_{x_{k}, u_{k}}\left(u-u_{k}\right) \\ & +\left.\frac{1}{2}\left(x-x_{k}\right)^{T} \frac{\partial^{2} Q}{\partial x^{2}}\right|_{x_{k}, u_{k}}\left(x-x_{k}\right)+\left.\frac{1}{2}\left(u-u_{k}\right)^{T} \frac{\partial^{2} Q}{\partial u^{2}}\right|_{x_{k}, u_{k}}\left(u-u_{k}\right) \\ & +\left.\frac{1}{2}\left(u-u_{k}\right)^{T} \frac{\partial^{2} Q}{\partial u \partial x}\right|_{x_{k}, u_{k}}\left(x-x_{k}\right)+\left.\frac{1}{2}\left(x-x_{k}\right)^{T} \frac{\partial^{2} Q}{\partial x \partial u}\right|_{x_{k}, u_{k}}\left(u-u_{k}\right) \end{aligned}$$
- 公式7：将上述Q函数扰动转化为简化形式，其对应的表达式如下所示
$$\delta Q_{k}\left(x_{k}, u_{k}\right)=Q_{x} \delta x+Q_{u} \delta u+\frac{1}{2} \delta x^{T} Q_{x x} \delta x+\frac{1}{2} \delta u^{T} Q_{u u} \delta u+\frac{1}{2} \delta x^{T} Q_{x u} \delta u+\frac{1}{2} \delta u^{T} Q_{u x} \delta x$$
-  公式8：Q函数扰动矩阵形式1
$$\delta Q_{k}\left(x_{k}, u_{k}\right)=\frac{1}{2}\left[\begin{array}{l} \delta x_{k} \\ \delta u_{k} \end{array}\right]^{T}\left[\begin{array}{ll}{Q_{x x}} & Q_{x u} \\ Q_{u x} & Q_{u u} \end{array}\right]\left[\begin{array}{l} \delta x_{k} \\ \delta u_{k} \end{array}\right]+\left[\begin{array}{l} Q_{x} \\ Q_{u} \end{array}\right]^{T}\left[\begin{array}{l} \delta x_{k} \\ \delta u_{k} \end{array}\right]$$
- 公式9：Q函数扰动矩阵形式2（含常数项维度）
$$\delta Q\left(x_{k}, u_{k}\right) \approx \frac{1}{2}\left[\begin{array}{c}1 \\ \delta x_{k} \\ \delta u_{k}\end{array}\right]^{T}\left[\begin{array}{ccc}0 & Q_{x}^{T} & Q_{u}^{T} \\ Q_{x} & Q_{x x} & Q_{x u} \\ Q_{u} & Q_{x u}^{T} & Q_{u u}\end{array}\right]\left[\begin{array}{c}1 \\ \delta x_{k} \\ \delta u_{k}\end{array}\right]$$
#### 根据上述会中得到的展开式即可得到对应的Q函数偏导数性质
- $$Q_{u x}=Q_{x u}^{T}$$
- 公式11：Q函数一阶偏导数（状态方向）
$$Q_{x}=\frac{\partial Q}{\partial x}=\frac{\partial \mathcal{L}}{\partial x}+\frac{\partial f^{T}}{\partial x} \frac{\partial V_{k+1}}{\partial x}$$

- 公式12：Q函数一阶偏导数（控制方向）
$$Q_{u}=\frac{\partial Q}{\partial u}=\frac{\partial \mathcal{L}}{\partial u}+\frac{\partial f^{T}}{\partial u} \frac{\partial V_{k+1}}{\partial x}$$

- 式13：Q函数二阶偏导数（状态-状态）
$$Q_{x x}=\frac{\partial^{2} Q}{\partial x^{2}}=\frac{\partial^{2} \mathcal{L}}{\partial x^{2}}+\frac{\partial f^{T}}{\partial x} \frac{\partial^{2} V_{k+1}}{\partial x^{2}} \frac{\partial f}{\partial x}+\frac{\partial V_{k+1}}{\partial x} \frac{\partial^{2} f}{\partial x^{2}}$$
- 公式14：Q函数二阶偏导数（控制-控制）
$$Q_{u u}=\frac{\partial^{2} Q}{\partial u^{2}}=\frac{\partial^{2} \mathcal{L}}{\partial u^{2}}+\frac{\partial f^{T}}{\partial u} \frac{\partial^{2} V_{k+1}}{\partial x^{2}} \frac{\partial f}{\partial u}+\frac{\partial V_{k+1}}{\partial x} \frac{\partial^{2} f}{\partial u^{2}}$$

- 公式15：Q函数二阶偏导数（状态-控制）
$$Q_{x u}=\frac{\partial^{2} Q}{\partial x \partial u}=\frac{\partial^{2} \mathcal{L}}{\partial x \partial u}+\frac{\partial f^{T}}{\partial x} \frac{\partial^{2} V_{k+1}}{\partial x^{2}} \frac{\partial f}{\partial u}+\frac{\partial V_{k+1}}{\partial x} \frac{\partial^{2} f}{\partial x \partial u}$$

#### 根据上述代价函数的定义以及V函数和Q函数的定义，其中最优控制扰动求解如下（一阶条件）
$$\frac{\partial \delta Q}{\partial \delta u}=Q_{u}+\frac{1}{2}Q_{u x}\delta x+\frac{1}{2}Q_{x u}^{T}\delta x+Q_{uu}\delta u=0$$

- 最优控制扰动表达式由Q函数进行求导得到对应的极值，对应表达式如下
$$\begin{aligned} \delta u^{*} & =-Q_{u u}^{-1}\left(Q_{u x} \delta x_{k}+Q_{u}\right) \\ & =K \delta x+d \end{aligned}$$
$$d_{i}=-Q_{u u}^{-1} Q_{u}\\
K=-Q_{u u}^{-1} Q_{u x}$$


#### 公式20：价值函数扰动展开（代入最优控制）
$$\begin{aligned} \delta V & =\delta Q\left(\delta x, \delta u^{*}\right)=\frac{1}{2}\left[\begin{array}{c} 1 \\ \delta x_{k} \\ \left(K \delta x_{k}+d\right) \end{array}\right]^{T}\left[\begin{array}{ccc} 0 & Q_{x}^{T} & Q_{u}^{T} \\ Q_{x} & Q_{x x} & Q_{u x} \\ Q_{u} & Q_{x u}^{T} & Q_{u u} \end{array}\right]\left[\begin{array}{c} 1 \\ \delta x_{k} \\ \left(K \delta x_{k}+d\right) \end{array}\right] \\ & =\left(Q_{x}+K^{T} Q_{u}+Q_{u x}^{T} d\right)^{T} \delta x_{k}+\frac{1}{2} \delta x_{k}^{T}\left(Q_{x x}+K^{T} Q_{u u} K+K^{T} Q_{u x}+Q_{u x}^{T} K\right) \delta x_{k}+\frac{1}{2} d^{T} Q_{u u} d+d^{T} Q_{u} \end{aligned}$$
根据上述得到的最终表达式以及一开始定义的V函数和Q函数的表达式即可得到对应的一阶偏导和二阶偏导的表达式$$\frac{\partial V}{\partial x}=Q_{x}+K^{T} Q_{u u} d+K^{T} Q_{u}+Q_{u x}^{T} d$$
$$\frac{\partial^{2} V}{\partial x^{2}}=Q_{x x}+K^{T} Q_{u u} K+K^{T} Q_{u x}+Q_{u x}^{T} K$$
$$\Delta V=\frac{1}{2} d^{T} Q_{u u} d+d^{T} Q_{u}$$

### 2.4 正向传播（Forward pass）
#### 公式24：状态偏差定义
$$\delta x_{k}=\overline{x}_{k}-x_{k}$$

#### 公式25：控制扰动更新（含线搜索步长α）
$$\delta u_{k}=K_{k} \delta x_{k}+\alpha d_{k}$$

#### 公式26：控制量更新
$$\overline{u}_{k}=u_{k}+\delta u_{k}$$

#### 公式27：状态量更新（动力学模型）
$$\overline{x}_{k+1}=f\left(\overline{x}_{k}, \overline{u}_{k}\right)$$

#### 公式28：线搜索评价指标z定义
$$z=\frac{J(X, U)-J(\overline{X}, \overline{U})}{-\Delta V(\alpha)}$$

#### 公式29：价值增量ΔV（修正形式）
$$\Delta V=\frac{1}{2} d^{T} Q_{u u} d-d^{T} Q_{u}$$

#### 公式30：带步长的价值增量ΔV(α)
$$\Delta V(\alpha )=\sum _{k=0}^{N-1} \alpha d_{k}^{T} Q_{u}+\alpha ^{2} \frac {1}{2}d_{k}^{T} Q_{u u} d_{k}$$

#### 公式31：总代价函数J定义
$$J(X, U) = \mathcal{L}(x_N) + \sum_{k=0}^{N-1} \mathcal{L}(x_k, u_k)$$

## 3. 补充说明（Side note: If we use RL?）
无直接公式，核心思路：

