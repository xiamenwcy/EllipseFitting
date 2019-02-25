# EllipseFitting

注意：本代码参考了 https://github.com/seisgo/EllipseFit  在此表示感谢！

基于代数距离的椭圆拟合，主要参考论文《[Andrew W. Fitzgibbon, Maurizio Pilu, and Robert B. Fisher
Direct least-squares fitting of ellipses,
IEEE Transactions on Pattern Analysis and Machine Intelligence, 21(5), 476--480, May 1999](http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/FITZGIBBON/ELLIPSE/)》

本程序使用C++完成。其中我们使用了**Eigen**与**clapack**库。其中Eigen易于表达矩阵，和matlab用法类似，是个强大的C\++线性代数库。而CLAPACK是线性代数包Lapack面向C/c++的接口。里面包含了很丰富的线性代数算法，包括**广义特征值求解接口**，而且速度很快。我们希望将二者结合起来使用。

## Eigen的安装
Eigen直接以源代码的方式提供给用户，因此我们从[官网](http://eigen.tuxfamily.org/index.php?title=Main_Page)上下载下后，直接在工程中包含其头文件路径即可。具体可参考：http://blog.csdn.net/abcjennifer/article/details/7781936

## clapack的安装
请查看[官网](http://icl.cs.utk.edu/lapack-for-windows/clapack/index.html#build)，里面包含了详细的使用与安装步骤。

也可以使用我们已经编译了的vc2010和vc2013的库，可以[点击下载](http://7xs15g.com1.z0.glb.clouddn.com/Clapack/clapack.rar)。

尽管clapack面向c语言，因此需要我们在包含头文件的时候，记得加上extern "C".但是最新的版本（比如CLAPACK 3.2.1）已经为我们在头文件中加上了这些限制符，因此最新的版本可以兼容c和c\++，所以直接在项目包含头文件即可。

比如像下面一样：

```
//Eigen
#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>

//clapack,必须放在Eigen后面

#include <f2c.h>
#include <clapack.h>

```
而且应该注意**Eigen与CLAPACK混合使用的时候，CLAPACK的头文件要加在Eigen的后面。否则会出错**。
## 实例
![image](https://img-blog.csdnimg.cn/20190225093925349.png?x-oss-process=image/watermark,type_ZmFuZ3poZW5naGVpdGk,shadow_10,text_aHR0cHM6Ly9ibG9nLmNzZG4ubmV0L3hpYW1lbnRpbmd0YW8=,size_16,color_FFFFFF,t_70)
