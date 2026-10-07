// ============================================================================
// 模块总览：ePlace 全局布局的 FFT 加速模块（头文件）
// ----------------------------------------------------------------------------
// 【职责】
//   实现 ePlace（静电场类比布局算法）中「由单元密度求电场」这一核心步骤：
//   把每个 bin 的单元面积密度当成电荷密度 ρ，解一次泊松方程得到电势 φ 与电场 E，
//   单元再沿 E 方向被推开，从而实现全局均匀化。
//   唯一调用方是 EPlace/eplace.cpp 的 EPlacer_2D::densityGradientUpdate()。
//
// 【为什么用 FFT】
//   静电学里库仑斥力是 O(N²) 的两两求和（每对 bin 之间按 1/r 相互作用），
//   N 为 bin 数，规模一大就不可接受。等价地，电势满足泊松方程
//       ∇²φ = −ρ ，    E = −∇φ
//   而 cos / sin 基函数是拉普拉斯算子 ∇² 的特征函数，变换到频域后 ∇² 退化成
//   「逐点乘一个数」，于是求解只需三步：正变换 → 逐点除法 → 反变换，
//   复杂度从 O(N²) 降到 O(N log N)。这就是本模块存在的全部理由。
//
// 【频域里做的那个除法是什么】
//   设谱系数下标为 (i, j)，对应波数（见 FFT_2D::init()）
//       w_x = π·i / binCntX ，  w_y = π·j / binCntY（bin 非正方形时再乘尺寸比）
//   基函数 cos(w_x x)·cos(w_y y) 是 ∇² 的特征函数，特征值为 −(w_x² + w_y²)，因此
//       ∇²φ = −ρ   ⇒   −(w_x² + w_y²)·φ̂ = −ρ̂   ⇒   φ̂ = ρ̂ / (w_x² + w_y²)
//   即：频域里除以的正是拉普拉斯算子的特征值（的相反数抵消后的结果）。
//   (i, j) = (0, 0) 是直流（DC）分量，特征值为 0 不能除，代码里直接置 0，
//   等价于把电势的整体均值（规范自由度）定为 0，对电场 E = −∇φ 没有任何影响。
//
// 【电场是怎么从电势来的】
//   注意这里并不是「先算出 φ 的空域网格、再做数值差分」，而是在频域直接做解析求导：
//   对 cos 求导得到 sin，所以「在哪个方向求导，那个方向的基就变成 sin 基」，
//       Ê_x = +w_x · φ̂    →  反变换用 ddsct2d（X 方向 sin 基、Y 方向 cos 基）
//       Ê_y = +w_y · φ̂    →  反变换用 ddcst2d（X 方向 cos 基、Y 方向 sin 基）
//       φ   =      φ̂      →  反变换用 ddct2d （两个方向都是 cos 基）
//   这样得到的是谱精度的导数，不需要任何差分模板，也不损失阶数。
//
// 【主要接口】class FFT_2D
//   FFT_2D(binCntX, binCntY, binSizeX, binSizeY) 构造：固定网格规模 + 预建三角函数表
//   updateDensity(x, y, density)   写入 (x, y) 处 bin 的电荷密度 ρ
//   doFFT()                        正变换 → 频域除法 → 三次反变换（φ、E_x、E_y）
//   getElectroPhi(x, y)            取该 bin 的电势 φ
//   getElectroForce(x, y)          取该 bin 的电场 (E_x, E_y)
//
// 【典型调用方式】（摘自 EPlace/eplace.cpp::densityGradientUpdate）
//   FFT_2D fft(binDimension.x, binDimension.y, binStep.x, binStep.y);
//   for (每个 bin) fft.updateDensity(i, j, eDensity * invertedBinArea); // ρ = 面积/面积
//   fft.doFFT();
//   for (每个 bin) { bins[i][j]->E   = fft.getElectroForce(i, j);
//                    bins[i][j]->phi = fft.getElectroPhi(i, j); }
//   之后调用方把每个单元覆盖到的 bin 上的 E 按「重叠面积」加权累加，得到密度梯度 ∇D。
//
// 【底层来源】
//   fftsg.cpp   ：Ooura FFT 包的 1D 部分（cdft / rdft / ddct / ddst），第三方代码
//   fftsg2d.cpp ：2D 部分（ddct2d / ddsct2d / ddcst2d），用行列分离法复用 1D
//   fftsg3d.cpp ：3D 部分（ePlace-MS 遗留），当前未参与构建
// ============================================================================

#ifndef __REPLACE_FFT__
#define __REPLACE_FFT__

#include <vector>
#include <algorithm>

namespace replace {

// FFT_2D：2D 泊松方程 ∇²φ = −ρ 的谱方法求解器（DCT/DST 版）
//
// 数据布局：binDensity_ 等二维数组的第一维是 X（长度 binCntX_）、第二维是 Y
// （长度 binCntY_），即 a[x][y]；与底层 ddct2d(n1 = binCntX_, n2 = binCntY_, ...)
// 的约定一致，也和调用方 bins[x][y] 的顺序一致。
//
// 生命周期：构造（分配网格 + 建表）→ updateDensity() 若干次 → doFFT() 一次
//           → getElectroForce() / getElectroPhi() 若干次 → 析构。
// 网格尺寸在构造时固定，之后不可更改；doFFT() 可以被反复调用（每次迭代重算一遍）。
class FFT_2D {
  public:
    // 默认构造：只把 bin 数与 bin 尺寸清零（不分配数组，见下方 疑似问题 标注）
    FFT_2D();
    // 功能构造：给定 X/Y 方向的 bin 个数与 bin 的物理尺寸，立即 init() 分配内存并建表
    FFT_2D(int binCntX, int binCntY, float binSizeX, float binSizeY);
    // 析构：释放四张二维网格，并清空三角函数表 / 工作表
    ~FFT_2D();

    // input func
    // 输入接口：把 (x, y) 处 bin 的电荷密度 ρ 写进 binDensity_（未做边界检查）
    void updateDensity(int x, int y, float density);

    // do FFT
    // 核心：一次调用完成「正向 DCT → 频域除法求 φ̂、Ê → 三次反变换」，
    // 结束后 electroPhi_ / electroForceX_ / electroForceY_ 即可被读取
    void doFFT();

    // returning func
    // 取 (x, y) 处 bin 的电场 (E_x, E_y)，即 E = −∇φ 的两个分量
    std::pair<float, float> getElectroForce(int x, int y) const;
    // 取 (x, y) 处 bin 的电势 φ
    float getElectroPhi(int x, int y) const;

  private:
    // 2D array; width: binCntX_, height: binCntY_;
    // No hope to use Vector at this moment...
    float** binDensity_;
    float** electroPhi_;
    float** electroForceX_;
    float** electroForceY_;

    // 上面四个二维数组即四张网格：ρ（bin 密度）、φ（电势）、E_x、E_y（电场分量），
    // 尺寸均为 X×Y；doFFT() 里它们被原地（in-place）变换，正变换时是谱系数、反变换后是空域值

    // cos/sin table (prev: w_2d)
    // length:  max(binCntX, binCntY) * 3 / 2
    // 三角函数表：Ooura 库的 cos/sin 预计算表（内部变量名为 w），
    // 前 nw 项是复数 FFT 用的旋转因子，其后 nc 项是 ddct/ddst 用的余弦表；
    // 只做 DCT/DST 时长度取 max(n1, n2) * 3 / 2 即可
    std::vector<float> csTable_;

    // wx. length:  binCntX_
    // X 方向的波数 w_x = π·i / binCntX_（已按 X 方向归一化），及其平方
    std::vector<float> wx_;
    std::vector<float> wxSquare_;

    // wy. length:  binCntY_
    // Y 方向的波数 w_y = π·j / binCntY_ × (binSizeY_ / binSizeX_)，及其平方
    // （乘尺寸比是为了把 Y 方向的频率折合到与 X 相同的物理尺度上；
    //   但这个比值方向存疑，见 fft.cpp 中 init() 处的 疑似问题 标注）
    std::vector<float> wy_;
    std::vector<float> wySquare_;

    // work area for bit reversal (prev: ip)
    // length: round(sqrt( max(binCntX_, binCntY_) )) + 2
    // 位反转用的工作区（Ooura 库内部叫 ip），ip[0] / ip[1] 还兼作
    //「已生成的旋转因子表 / 余弦表的规模」的缓存：若当前规模不够大就在变换前重建
    std::vector<int> workArea_;

    // X/Y 方向的 bin 个数（即频域采样点数；Ooura 库要求为 2 的幂）
    int binCntX_;
    int binCntY_;
    // 单个 bin 的物理宽 / 高（布局数据库中的 binStep）
    float binSizeX_;
    float binSizeY_;

    // 分配四张二维网格、建三角函数表、算好 w_x / w_y 及其平方
    void init();
};

// ----------------------------------------------------------------------------
// 以下是底层 Ooura FFT 包的内部例程声明（定义在 fftsg.cpp），第三方代码，
// 不建议修改也不逐行注释：三角函数表生成（makewt / makect / makeipt）、
// 位反转（bitrv2 系列）、以及 split-radix 复数 FFT 的各种蝶形核（cft* 系列）。
// ----------------------------------------------------------------------------

void makewt(int nw, int *ip, float *w);
void cftfsub(int n, float *a, int *ip, int nw, float *w);
void cftbsub(int n, float *a, int *ip, int nw, float *w);
void makect(int nc, int *ip, float *c);
void rftfsub(int n, float *a, int nc, float *c);
void rftbsub(int n, float *a, int nc, float *c);
void dctsub(int n, float *a, int nc, float *c);
void dstsub(int n, float *a, int nc, float *c);
void makeipt(int nw, int *ip);
void bitrv2(int n, int *ip, float *a);
void bitrv216(float * a);
void bitrv208(float * a);
void cftf1st(int n, float *a, float *w);
void cftrec4(int n, float *a, int nw, float *w);
void cftleaf(int n, int isplt, float *a, int nw, float *w);
void cftfx41(int n, float *a, int nw, float *w);
void cftf161(float * a, float * w);
void cftf081(float * a, float * w);
void cftf040(float * a);
void cftx020(float * a);
void bitrv2conj(int n, int *ip, float *a);
void bitrv216neg(float * a);
void bitrv208neg(float * a);
void cftb1st(int n, float *a, float *w);
void cftb040(float * a);
int cfttree(int n, int j, int k, float *a, int nw, float *w);
void cftmdl1(int n, float *a, float *w);
void cftmdl2(int n, float *a, float *w);
void cftf162(float * a, float * w);
void cftf082(float * a, float * w);
void cdft(int n, int isgn, float *a, int *ip, float *w);
void cdft2d_sub(int n1, int n2, int isgn, float **a, float *t, int *ip,
	float *w);
void rdft(int n, int isgn, float *a, int *ip, float *w);
void rdft2d_sub(int n1, int isgn, float **a);
void ddxt2d_sub(int n1, int n2, int ics, int isgn, float **a, float *t, int *ip,
	float *w);

/// 1D FFT ////////////////////////////////////////////////////////////////
// 1D 变换（定义见 fftsg.cpp）：
//   cdft：复数 DFT；ddct：实数余弦变换 DCT；ddst：实数正弦变换 DST
//   isgn = -1 表示正变换，isgn = +1 表示反变换；a 为原地输入/输出数组
// 这三者是 2D 封装的全部基石：2D 变换只是把它们按行、按列各跑一遍
void cdft(int n, int isgn, float *a, int *ip, float *w);
void ddct(int n, int isgn, float *a, int *ip, float *w);
void ddst(int n, int isgn, float *a, int *ip, float *w);

/// 2D FFT ////////////////////////////////////////////////////////////////
// 2D 变换（定义见 fftsg2d.cpp），命名规则为 dd[X][Y]2d，其中字母 c = cos 基（DCT）、
// s = sin 基（DST），第一个字母对应第一维 n1（本项目即 X 方向），第二个字母对应 n2（Y）：
//   ddct2d  = X 用 cos、Y 用 cos  → 电势 φ 的反变换（以及密度 ρ 的正变换）
//   ddsct2d = X 用 sin、Y 用 cos  → 电场 E_x 的反变换（对 X 方向求过导）
//   ddcst2d = X 用 cos、Y 用 sin  → 电场 E_y 的反变换（对 Y 方向求过导）
// 参数：n1/n2 为两维长度，isgn 同 1D，a 为 a[n1][n2] 二维数组（原地），
//       t 为临时缓冲（传 NULL 时由库内部分配释放），ip/w 为工作表与三角函数表
// 实现上是行列分离（row-column）：先对每一行做 1D 变换，再对每一列做 1D 变换
void cdft2d(int, int, int, float **, float *, int *, float *);
void rdft2d(int, int, int, float **, float *, int *, float *);
void ddct2d(int, int, int, float **, float *, int *, float *);
void ddst2d(int, int, int, float **, float *, int *, float *);
void ddsct2d(int n1, int n2, int isgn, float **a, float *t, int *ip, float *w);
void ddcst2d(int n1, int n2, int isgn, float **a, float *t, int *ip, float *w);
//! 疑似问题：cdft2d / rdft2d / ddst2d（以及上面的 cdft2d_sub / rdft2d_sub）
//! 在 fftsg.cpp 与 fftsg2d.cpp 中都没有实现，一旦调用会链接失败；
//! 目前 ePlace 只用到 ddct2d / ddsct2d / ddcst2d 三个，所以没有暴露出来。

/// 3D FFT ////////////////////////////////////////////////////////////////
// 3D 变换：ePlace-MS（多层级/三维版）遗留接口，定义见 fftsg3d.cpp。
// 注意 fftsg3d.cpp 并未被 EPlace/CMakeLists.txt 纳入构建，且文件中使用的
// prec 类型在本仓库里没有定义，因此这一组接口目前实际不可用（详见该文件头部说明）。
void cdft3d(int, int, int, int, float ***, float *, int *, float *);
void rdft3d(int, int, int, int, float ***, float *, int *, float *);
void ddct3d(int, int, int, int, float ***, float *, int *, float *);
void ddst3d(int, int, int, int, float ***, float *, int *, float *);
void ddscct3d(int, int, int, int isgn, float ***, float *, int *, float *);
void ddcsct3d(int, int, int, int isgn, float ***, float *, int *, float *);
void ddccst3d(int, int, int, int isgn, float ***, float *, int *, float *);

}

#endif
