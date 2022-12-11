### Parallel Computing

添加了并行计算，加快运行速度；调整了代码参数并记录

###  Algorithm Modification on Chang18

在第2版中利用chang18年的算法，修改代码

### Second EQ Version

估计方程：$\mathbf{g}_{i j}(\boldsymbol{\theta})=z_{i j}\left(y_{i}-\mu_{i}-\boldsymbol{x}_{i}^{\mathrm{T}} \boldsymbol{\beta}-y_{j}+\mu_{j}+\boldsymbol{x}_{j}^{\mathrm{T}} \boldsymbol{\beta}\right)$

目标函数：$\widehat{\boldsymbol{\theta}}_{m}=\arg \min _{\boldsymbol{\theta} \in \boldsymbol{\Theta}} \max _{\lambda \in \widehat{\Lambda}_{m}(\boldsymbol{\theta})}\left\{\frac{1}{m} \sum_{i<j} \log _{\star}\left[1+\boldsymbol{\lambda}^{\mathrm{T}} \mathbf{g}_{i j}(\boldsymbol{\theta})\right]-\sum_{i<j}^{1 \leq i, j \leq n} P_{2, \eta_{2}}\left(\left|\lambda_{i j}\right|\right)\right.\left.+\sum_{s<t}^{1 \leq s, t \leq n} P_{1, \eta_{1}}\left(\left|\xi_{s t}\right|\right)+\boldsymbol{\delta}^{\mathrm{T}} \mathrm{D} \boldsymbol{\xi}\right\}$

其中：$D\xi=0$使得对所有$i、j、k$成立$\xi_{j k}=\mu_{j}-\mu_{k}=\left(\mu_{i}-\mu_{k}\right)-\left(\mu_{i}-\mu_{j}\right)=\xi_{i k}-\xi_{i j}, \quad \text { for } \quad i<j<k$

### Init

估计方程：$\mathbf{g}_{\mathbf{i}}(\boldsymbol{\theta})=z_{i}\left(y_{i}-z_{i}^{T} \boldsymbol{\theta}\right)$

目标函数：$\widehat{\boldsymbol{\theta}}_{n}=\arg \min _{\boldsymbol{\theta} \in \boldsymbol{\Theta}} \max _{\lambda \in \widehat{\Lambda}_{n}(\boldsymbol{\theta})}\{\sum_{i=1}^{n} \log _{\star}\left[1+\boldsymbol{\lambda}^{\mathrm{T}} \mathbf{g}_{\mathbf{i}}(\boldsymbol{\theta})\right]+n \sum_{i<j}^{1 \leq i, j \leq n} P_{1, \eta_{1}}\left(\left|\mu_{i}-\mu_{j}\right|\right)\left.-n \sum_{s=1}^{r} P_{2, \eta_{2}}\left(\left|\lambda_{s}\right|\right)\right\}$
