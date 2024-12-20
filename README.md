# ScalarDiffractionAnaylsis
Requires: gsl

The GNU `gsl` Library could be download from Github directly. 
## Folder Symmetric
圆环区域上的光圈。i.e. `u=u(rho,z)` and ,`u(rho,0)=1` for `rho` in `[a,b]` and `0` otherwise.

***REMARK***: 该分支的 symmetric case 代码存在问题 `exp(mikz)`. `RingInput` 分支是正确的。

## Folder Angular
增加了 Theta 方向，输入的图像是圆环且有部分角度为`0`. (类似望远镜的副镜支架)
