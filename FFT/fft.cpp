// ============================================================================
// 模块总览：FFT_2D 的实现 —— 用 DCT/DST 谱方法解 2D 泊松方程 ∇²φ = −ρ
// ----------------------------------------------------------------------------
// 【这个文件在干什么】
//   输入：一张 binCntX_ × binCntY_ 的 bin 密度网格 ρ（由调用方用 updateDensity 填入，
//         物理含义是「该 bin 内单元面积 / bin 面积」，相当于面电荷密度）。
//   输出：同一张网格上的电势 φ 与电场 E = (E_x, E_y) = −∇φ。
//   全局布局器拿 E 去推单元（沿 E 方向移动即互相散开），见 EPlace/eplace.cpp。
//
// 【数学背景】
//   直接用库仑定律两两求和是 O(N²)；改用泊松方程的谱解法：
//       ∇²φ = −ρ ，  E = −∇φ
//   DCT 的基函数 cos(w_x x)·cos(w_y y) 是 ∇² 的特征函数，特征值 −(w_x² + w_y²)，
//   于是频域里 φ̂ = ρ̂ / (w_x² + w_y²)（一步除法就完成了整个泊松求解），
//   再乘波数即得电场谱 Ê_x = w_x·φ̂、Ê_y = w_y·φ̂，反变换回空域即可。
//   总复杂度 O(N log N)，N = binCntX_ × binCntY_。
//
// 【doFFT() 的四个步骤】（详见该函数处的行内注释）
//   ① ddct2d(isgn = −1)：把 ρ 正变换成谱系数 ρ̂（原地覆盖 binDensity_）
//   ② 归一化与边界项修正：乘 4/(binCntX·binCntY)，并把 k = 0 那一维乘 0.5
//   ③ 频域求解：φ̂ = ρ̂/(w_x² + w_y²)，Ê_x = w_x·φ̂，Ê_y = w_y·φ̂；DC 分量置 0
//   ④ 三次反变换：φ 用 ddct2d(+1)、E_x 用 ddsct2d(+1)、E_y 用 ddcst2d(+1)
//
// 【为什么反变换还分三种】
//   对 cos 基求导会变成 sin 基：在哪个方向求了导，那个方向就必须用 sin 基反变换，
//   所以 φ 用全 cos（ddct2d），E_x 用「X 是 sin、Y 是 cos」（ddsct2d），
//   E_y 用「X 是 cos、Y 是 sin」（ddcst2d）。这是本文件最容易被看漏的一处细节。
//
// 【调用方式】
//   FFT_2D fft(binDimension.x, binDimension.y, binStep.x, binStep.y);
//   …fft.updateDensity(i, j, eDensity);…  fft.doFFT();  …fft.getElectroForce(i, j);…
//   全局布局每迭代一次就完整跑一遍上面的流程。
// ============================================================================

#include <cstdlib>
#include <cmath>
#include <cfloat>

#include <iostream>

#include "fft.h"

#define REPLACE_FFT_PI 3.141592653589793238462L 

namespace replace {


// 默认构造：只清零 bin 数与 bin 尺寸
//! 疑似问题：这里没有把 binDensity_ / electroPhi_ / electroForceX_ / electroForceY_
//! 初始化为 nullptr，它们此时是野指针；一旦对这个默认构造出来的对象调用析构函数，
//! 就会 delete 野指针（未定义行为）。此外本类没有定义拷贝构造与拷贝赋值（违反
//! Rule of Three），按值拷贝/返回会导致同一批指针被重复释放。均属于历史遗留，未修改。
FFT_2D::FFT_2D()
  : binCntX_(0), binCntY_(0), binSizeX_(0.0), binSizeY_(0.0) {}

// 功能构造：记录网格规模与 bin 物理尺寸，然后 init() 分配内存、建立三角函数表与波数表
FFT_2D::FFT_2D(int binCntX, int binCntY, float binSizeX, float binSizeY)
  : binCntX_(binCntX), binCntY_(binCntY), 
  binSizeX_(binSizeX), binSizeY_(binSizeY) {
  init();   
}

// 析构：逐行释放四张二维网格，再释放指针数组本身，最后清空各 std::vector 缓冲
//! 疑似问题：数组是用 new float[binCntY_] 分配的，这里却用 delete（标量形式）释放，
//! 正确写法应为 delete[]；对数组用标量 delete 属未定义行为（实践中某些编译器/平台
//! 能侥幸跑通）。同理末尾的 delete(binDensity_) 也应为 delete[]。仅标注，未修改。
FFT_2D::~FFT_2D() {
  using std::vector;
  for(int i=0; i<binCntX_; i++) {
    delete(binDensity_[i]);
    delete(electroPhi_[i]);
    delete(electroForceX_[i]);
    delete(electroForceY_[i]);
  }
  delete(binDensity_);
  delete(electroPhi_);
  delete(electroForceX_);
  delete(electroForceY_);


  csTable_.clear();
  wx_.clear();
  wxSquare_.clear();
  wy_.clear();
  wySquare_.clear();
  
  workArea_.clear();
}


// init()：一次性准备工作
//   1) 分配四张 X×Y 的二维网格（bin 密度 ρ、电势 φ、电场 E_x / E_y）并清零；
//      注意第一维是 X（binCntX_ 个指针），第二维是 Y（每个指针指向 binCntY_ 个 float）
//   2) 分配 Ooura 库需要的三角函数表 csTable_ 与位反转工作表 workArea_（均初始化为 0，
//      其中 workArea_[0] / [1] 为 0 表示「表还没生成」，首次变换时会自动 makewt / makect）
//   3) 预计算两个方向的波数及其平方，供 doFFT() 的频域除法直接使用
void
FFT_2D::init() {
  binDensity_ = new float*[binCntX_];
  electroPhi_ = new float*[binCntX_];
  electroForceX_ = new float*[binCntX_];
  electroForceY_ = new float*[binCntX_];

  for(int i=0; i<binCntX_; i++) {
    binDensity_[i] = new float[binCntY_];
    electroPhi_[i] = new float[binCntY_];
    electroForceX_[i] = new float[binCntY_];
    electroForceY_[i] = new float[binCntY_];

    for(int j=0; j<binCntY_; j++) {
      binDensity_[i][j] 
        = electroPhi_[i][j] 
        = electroForceX_[i][j] 
        = electroForceY_[i][j] 
        = 0.0f;  
    }
  }

  csTable_.resize( std::max(binCntX_, binCntY_) * 3 / 2, 0 );

  wx_.resize( binCntX_, 0 );
  wxSquare_.resize( binCntX_, 0);
  wy_.resize( binCntY_, 0 );
  wySquare_.resize( binCntY_, 0 );

  workArea_.resize( round(sqrt(std::max(binCntX_, binCntY_))) + 2, 0 );
 
  // X 方向的波数：w_x = π·i / binCntX_（i 为频域下标，已按 X 方向的 bin 数归一化）
  for(int i=0; i<binCntX_; i++) {
    wx_[i] = REPLACE_FFT_PI * static_cast<float>(i) 
      / static_cast<float>(binCntX_);
    wxSquare_[i] = wx_[i] * wx_[i]; 
  }

  // Y 方向的波数：w_y = π·j / binCntY_，再乘 binSizeY_/binSizeX_ 做尺度折合
  //! 疑似问题：若按物理波数推导，X 方向 w_x ∝ π·i/(binCntX·binSizeX)、
  //! Y 方向 w_y ∝ π·j/(binCntY·binSizeY)，统一到 X 的物理尺度后 Y 应乘
  //! (binSizeX_/binSizeY_)，而这里乘的是其倒数 (binSizeY_/binSizeX_)。
  //! 当 bin 为正方形（binSizeX_ == binSizeY_，ePlace 里的常见情形）时两者相同、无影响，
  //! 但 bin 非正方形时 Y 方向的电势与电场会被反向缩放。仅标注，未修改。
  for(int i=0; i<binCntY_; i++) {
    wy_[i] = REPLACE_FFT_PI * static_cast<float>(i)
      / static_cast<float>(binCntY_) 
      * static_cast<float>(binSizeY_) 
      / static_cast<float>(binSizeX_);
    wySquare_[i] = wy_[i] * wy_[i];
  }
}

// 写入 (x, y) 处 bin 的电荷密度 ρ。调用方传入的是「bin 内单元面积 / bin 面积」
// （eplace.cpp 中还包含 filler 面积，filler 也要参与产生斥力）
//! 疑似问题：未做 (x, y) 的边界检查，越界会直接越界写；下同（两个 getter 亦如此）
void
FFT_2D::updateDensity(int x, int y, float density) {
  binDensity_[x][y] = density;
}

// 取 (x, y) 处 bin 的电场 (E_x, E_y) = −∇φ。调用方按「单元与 bin 的重叠面积」
// 把它加权累加，得到该单元的密度梯度
std::pair<float, float> 
FFT_2D::getElectroForce(int x, int y) const {
  return std::make_pair(
      electroForceX_[x][y],
      electroForceY_[x][y]);
}

// 取 (x, y) 处 bin 的电势 φ（DC 分量已被置 0，故 φ 的整体均值为 0）
float
FFT_2D::getElectroPhi(int x, int y) const {
  return electroPhi_[x][y]; 
}

// 注意：原作者在此把整个 std 命名空间引入全局作用域，
// 因此后面的 round() / sqrt() 等都没有加 std:: 限定（保持原样，未修改）
using namespace std;

// ----------------------------------------------------------------------------
// doFFT()：本模块的核心，一次调用完成整个泊松求解
//
//   ① 正变换   ρ ──ddct2d(−1)──▶ ρ̂            （原地覆盖 binDensity_）
//   ② 归一化   乘 4/(Nx·Ny)，并把 k = 0 的那一维乘 0.5
//   ③ 频域解   φ̂ = ρ̂ / (w_x² + w_y²)          ← 除的是拉普拉斯算子的特征值
//              Ê_x = w_x·φ̂ ， Ê_y = w_y·φ̂     ← 解析求导（cos → sin），不是数值差分
//   ④ 反变换   φ̂ ──ddct2d(+1)──▶ φ            （两维都是 cos 基）
//              Ê_x ──ddsct2d(+1)──▶ E_x       （X 维 sin 基 —— 在 X 方向求的导）
//              Ê_y ──ddcst2d(+1)──▶ E_y       （Y 维 sin 基 —— 在 Y 方向求的导）
//
// 完成后 electroPhi_ / electroForceX_ / electroForceY_ 即为可直接读取的空域结果。
// ----------------------------------------------------------------------------
void
FFT_2D::doFFT() {
  // ① 正向 2D DCT：把空域的 bin 密度 ρ 变成谱系数 ρ̂
  //    t 传 NULL 表示让库内部分配临时缓冲（每轮迭代都会 malloc/free 一次，
  //    若追求性能可预先分配一块复用 —— 这里保持原样未改）
  ddct2d(binCntX_, binCntY_, -1, binDensity_, 
      NULL, (int*) &workArea_[0], (float*)&csTable_[0]);

  // ② Ooura 的 DCT 约定：k = 0 的那一维要额外乘 0.5 才与标准 DCT-II 一致，
  //    因此 j = 0（Y 方向第 0 个频率）的所有系数乘 0.5
  for(int i = 0; i < binCntX_; i++) {
    binDensity_[i][0] *= 0.5;
  }

  //    同理 i = 0（X 方向第 0 个频率）的所有系数乘 0.5；
  //    两个循环叠加后 [0][0] 共乘 0.25，正是「两个维度都取了 k=0」的正确结果
  for(int i = 0; i < binCntY_; i++) {
    binDensity_[0][i] *= 0.5;
  }

  //    整体归一化系数 4/(binCntX·binCntY)：配合 Ooura 正反变换的 2/n 归一化约定，
  //    使「正变换 → 反变换」正好还原原函数（两次各带 2/n，这里预先补掉）
  for(int i = 0; i < binCntX_; i++) {
    for(int j = 0; j < binCntY_; j++) {
      binDensity_[i][j] *= 4.0 / binCntX_ / binCntY_;
    }
  }

  // ③ 在频域逐点求解：对每个频率 (i, j) 做除法与求导
  for(int i = 0; i < binCntX_; i++) {
    float wx = wx_[i];      // X 方向波数 w_x
    float wx2 = wxSquare_[i]; // w_x²

    for(int j = 0; j < binCntY_; j++) {
      float wy = wy_[j];        // Y 方向波数 w_y
      float wy2 = wySquare_[j]; // w_y²

      float density = binDensity_[i][j];
      float phi = 0;
      float electroX = 0, electroY = 0;

      // (i, j) = (0, 0) 是 DC（直流）分量：拉普拉斯算子在该频率上的特征值为 0，
      // 无法做除法，直接置 0 —— 物理上等价于把电势的整体均值（规范自由度）定为 0，
      // 由于 E = −∇φ，加一个常数电势不会改变电场，所以这个选择不影响结果
      if(i == 0 && j == 0) {
        phi = electroX = electroY = 0.0f;
      }
      else {
        // 频域求解的核心两行：
        //   ∇²φ = −ρ，而 ∇² 在 cos(w_x x)·cos(w_y y) 上的特征值是 −(w_x² + w_y²)，
        //   故 φ̂ = ρ̂ / (w_x² + w_y²)   ← 这个除法就是「解泊松方程」本身
        //   E = −∇φ：对 cos 求导变 sin，故 Ê_x = +w_x·φ̂、Ê_y = +w_y·φ̂
        //   （是频域里的解析求导，不是回到空域再做数值差分）
        //! 疑似问题：下方被注释掉的 lutong 推导版与现行实现差常数因子 ——
        //! 注释版 denom = (w_x² + w_y²)/4（相当于 φ̂ 再乘 4），且电场还要再乘 1/2；
        //! 现行实现分别是它的 1/4 与 1/2，两者的比例并不自洽。可能是不同 DCT 归一化
        //! 约定下的等价写法，也可能是历史遗留的不一致。仅标注，未修改。
        //////////// lutong
        //  denom =
        //  wx2 / 4.0 +
        //  wy2 / 4.0 ;
        // a_phi = a_den / denom ;
        ////b_phi = 0 ; // -1.0 * b / denom ;
        ////a_ex = 0 ; // b_phi * wx ;
        // a_ex = a_phi * wx / 2.0 ;
        ////a_ey = 0 ; // b_phi * wy ;
        // a_ey = a_phi * wy / 2.0 ;
        ///////////
        phi = density / (wx2 + wy2);
        electroX = phi * wx;
        electroY = phi * wy;
      }
      // 存下该频率上的谱系数（下一步会被反变换原地覆盖成空域值）
      electroPhi_[i][j] = phi;
      electroForceX_[i][j] = electroX;
      electroForceY_[i][j] = electroY;
    }
  }

  // ④ 反变换（isgn = +1）：把三个谱场变回空域
  // Inverse DCT
  //   φ̂ 两个方向都是 cos 基 → 用 ddct2d（c + c）
  ddct2d(binCntX_, binCntY_, 1, 
      electroPhi_, NULL, 
      (int*) &workArea_[0], (float*) &csTable_[0]);
  //   Ê_x 是对 X 求过导的，X 方向应为 sin 基、Y 方向仍是 cos 基 → ddsct2d（s + c）
  ddsct2d(binCntX_, binCntY_, 1, 
      electroForceX_, NULL, 
      (int*) &workArea_[0], (float*) &csTable_[0]);
  //   Ê_y 是对 Y 求过导的，X 方向仍是 cos 基、Y 方向应为 sin 基 → ddcst2d（c + s）
  ddcst2d(binCntX_, binCntY_, 1, 
      electroForceY_, NULL, 
      (int*) &workArea_[0], (float*) &csTable_[0]);

}


}
