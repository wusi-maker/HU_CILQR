# CILQR / iLQR 推导公式（对应 `CILQR.m`）

本文给出与你这份 `CILQR.m`（iLQR/CILQR）实现一一对应的推导公式（LaTeX），按“模型 → 代价 → 二阶展开 → 反向递推 → 前向更新”的顺序组织。

---

## 1. 系统动力学（与 `vehicle_dynamics` 对应）

状态与控制：
$
x_k =
\begin{bmatrix}
p_x\\ p_y\\ v\\ \theta
\end{bmatrix},\quad
u_k =
\begin{bmatrix}
a\\ \omega
\end{bmatrix}
$

离散时间非线性动力学：
$
x_{k+1} = f(x_k,u_k) =
x_k + 
\begin{bmatrix}
v_k\cos\theta_k\\
v_k\sin\theta_k\\
a_k\\
\omega_k
\end{bmatrix}\Delta t
$

线性化（与 `get_linearized_dynamics` 对应）：
$
\delta x_{k+1} \approx A_k\,\delta x_k + B_k\,\delta u_k
$
$
A_k=\frac{\partial f}{\partial x}\Big|_{(x_k,u_k)},\qquad
B_k=\frac{\partial f}{\partial u}\Big|_{(x_k,u_k)}
$

对该模型显式写出：
$
A_k=
\begin{bmatrix}
1 & 0 & \Delta t\cos\theta_k & -\Delta t\,v_k\sin\theta_k\\
0 & 1 & \Delta t\sin\theta_k & \ \Delta t\,v_k\cos\theta_k\\
0 & 0 & 1 & 0\\
0 & 0 & 0 & 1
\end{bmatrix},
\quad
B_k=
\begin{bmatrix}
0 & 0\\
0 & 0\\
\Delta t & 0\\
0 & \Delta t
\end{bmatrix}
$

---

## 2. 优化目标（软约束形式）

总目标：
$
J = \sum_{k=0}^{N-1}\ell(x_k,u_k) + \ell_f(x_N)
$   

### 2.1 运行代价（与 `running_cost` 主体对应）

车道/速度/航向的二次型（用偏差 \(\Delta x_k\)）：
$
\Delta x_k =
\begin{bmatrix}
0\\
p_{y,k}-y_{\text{ref}}\\
v_k-v_{\text{ref}}\\
\theta_k-\theta_{\text{ref}}
\end{bmatrix}
$
$
\ell_{\text{track}}(x_k,u_k)
=\frac12\,\Delta x_k^\top Q\,\Delta x_k+\frac12\,u_k^\top R\,u_k
$

“向前推进”的线性奖励（代码里体现在 $l_{x(1)}$ 上）可写成：
$
\ell_{\text{prog}}(x_k)= - w_{\text{prog}}\,p_{x,k}
$

“禁止倒车”（软约束）常用写法之一（对 \(v<0\) 的二次惩罚）：
$
\ell_{v^-}(x_k)=
\begin{cases}
w_{v^-}\,v_k^2, & v_k<0\\
0, & v_k\ge 0
\end{cases}
$

避障软约束（指数障碍）：
$
d_k=\left\|p_k-p_{\text{obs}}\right\|,\quad p_k=\begin{bmatrix}p_{x,k}\\p_{y,k}\end{bmatrix}
$
一种等价的“势场”写法：
$
\ell_{\text{obs}}(x_k)= w_{\text{obs}}\exp\!\bigl(-(d_k-d_{\text{safe}})\bigr)
$

于是运行代价合在一起：
$
\ell(x_k,u_k)=
\ell_{\text{track}}(x_k,u_k)+\ell_{\text{prog}}(x_k)+\ell_{v^-}(x_k)+\ell_{\text{obs}}(x_k)
$

### 2.2 终端代价（与 `terminal_cost` 对应）

到达前方终点 + 保持车道/速度/航向：
$
\ell_f(x_N)=\frac12\,(x_N-x_{\text{goal}})^\top Q_f (x_N-x_{\text{goal}})
$
其中：
$
x_{\text{goal}}=
\begin{bmatrix}
p_{x,\text{goal}}\\ 0\\ v_{\text{ref}}\\ 0
\end{bmatrix}
$

---

## 3. iLQR 的二阶近似与 Q-函数展开（反向传播核心）

定义名义轨迹 $\{\bar x_k,\bar u_k\}$，扰动为：
$
\delta x_k = x_k-\bar x_k,\qquad \delta u_k=u_k-\bar u_k
$

对运行代价二阶展开（在 $(\bar x_k,\bar u_k)$ 处）：
$
\ell(x_k,u_k)\approx
\ell_k
+\ell_{x,k}^\top\delta x_k
+\ell_{u,k}^\top\delta u_k
+\frac12\delta x_k^\top \ell_{xx,k}\delta x_k
+\frac12\delta u_k^\top \ell_{uu,k}\delta u_k
+\delta u_k^\top \ell_{ux,k}\delta x_k
$

价值函数二阶近似：
$
V_{k+1}(x_{k+1})\approx
V_{k+1}
+V_{x,k+1}^\top\delta x_{k+1}
+\frac12\delta x_{k+1}^\top V_{xx,k+1}\delta x_{k+1}
$

代入线性化动态 $\delta x_{k+1}=A_k\delta x_k+B_k\delta u_k$，得到 Q-函数二阶形式，并定义：
$
Q_x = \ell_x + A^\top V_x,\qquad
Q_u = \ell_u + B^\top V_x
$
$
Q_{xx} = \ell_{xx} + A^\top V_{xx}A,\qquad
Q_{uu} = \ell_{uu} + B^\top V_{xx}B
$
$
Q_{ux} = \ell_{ux} + B^\top V_{xx}A
$

---

## 4. 最优控制增量、反馈增益（与 `k_ff` 和 `K_gains` 对应）

令二次 Q-函数对 \(\delta u\) 最小，解一阶最优条件：
$
\frac{\partial}{\partial \delta u}
\left(
Q_u^\top\delta u + \frac12\delta u^\top Q_{uu}\delta u + \delta u^\top Q_{ux}\delta x
\right)=0
$

得到控制律（前馈 + 反馈）：
$
\delta u^\star = k + K\,\delta x
$
其中：
$$   
k = -Q_{uu}^{-1}Q_u,\qquad
K = -Q_{uu}^{-1}Q_{ux}
$$

---

## 5. 价值函数反向递推更新（与 `V_x`,`V_xx` 更新对应）

把 \(\delta u^\star\) 代回，得到：
$
V_x = Q_x + K^\top Q_{uu}k + K^\top Q_u + Q_{ux}^\top k
$
$
V_{xx} = Q_{xx} + K^\top Q_{uu}K + K^\top Q_{ux} + Q_{ux}^\top K
$

终端初始化：
$
V_{x,N}=\ell_{f,x}(x_N),\qquad
V_{xx,N}=\ell_{f,xx}(x_N)
$   

---

## 6. 前向更新（与 “Update Nominal Trajectory” 对应）

前向 rollout（步长 \(\alpha\)）：
$
u_k \leftarrow \bar u_k + \alpha\left(k_k + K_k\,\delta x_k\right)
$
$
x_{k+1}\leftarrow f(x_k,u_k),\qquad
\delta x_{k+1}\leftarrow x_{k+1}-\bar x_{k+1}
$

---

## 7. 避障项的梯度（与 `lx(1:2)` 中“排斥力”对应）

若采用势函数：
$
\ell_{\text{obs}}(p)=w_{\text{obs}}e^{-(d-d_{\text{safe}})},\quad d=\|p-p_{\text{obs}}\|
$

则梯度：
$
\nabla_p d=\frac{p-p_{\text{obs}}}{\|p-p_{\text{obs}}\|}
$
$
\nabla_p \ell_{\text{obs}}=-w_{\text{obs}}e^{-(d-d_{\text{safe}})}\nabla_p d
$

