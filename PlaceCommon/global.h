#ifndef GLOBAL_H
#define GLOBAL_H

// ============================================================================
// global.h —— 全工程公共头文件：基础类型 + 数值工具 + 计时/调试工具
//
// 内容分四块：
//   1. 常量与宏：浮点比较精度 EPS、取整宏 INT_CONVERT / INT_DOWN、shell 彩色输出
//   2. 几何类型：POS_2D / POS_3D（点）、VECTOR_2D / VECTOR_3D（向量）、CRect（矩形）
//      约定：打印时点用 () 、向量用 []，方便在日志里区分
//   3. 数值工具：带容差的浮点比较、浮点乘除、矩形重叠面积
//   4. 工程工具：计时（time_start/time_end）、调试输出、fastExp、快速选择 getKth
//
// 注意：这个头文件被几乎所有模块 include，改动它会影响全工程编译。
// ============================================================================

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <assert.h>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <map>
#include <list>
#include <set>
#include <limits.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <iomanip>
#include "string.h"
#include "arghandler.h"
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace std;
const string padding(30, '=');          //! 日志分隔线（30 个等号），打印小节标题用
#define EPS 1.0E-15 // for float number comparison //! 浮点比较容差；注意比 float 的有效精度(~1e-7)还小，对 float 而言基本等价于精确比较
#define DOUBLE_MAX __DBL_MAX__
#define DOUBLE_MIN __DBL_MIN__ //!!!!! double min > 0!!! //! 注意：这是最小的**正**规格化 double，不是最负数！想表示"负无穷大"要用 -DOUBLE_MAX
#define FLOAT_MAX __FLT_MAX__
#define FLOAT_MIN __FLT_MIN__ //!!!!! double min > 0!!! //! 同上，是最小的正 float
#define INT_CONVERT(a) (int)(1.0 * (a) + 0.5f) //! 四舍五入取整
#define INT_DOWN(a) (int)(a)                   //! 直接截断取整（向零取整，对负数不是 floor）
#define NEGATIVE_MAX_EXP -300                  //! e^-300 已小到可视为 0，用于指数下溢保护
//! The followings are colors used for log(in shell)
#define RESET "\033[0m"
#define BLACK "\033[30m"              /* Black */
#define RED "\033[31m"                /* Red */
#define GREEN "\033[32m"              /* Green */
#define YELLOW "\033[33m"             /* Yellow */
#define BLUE "\033[34m"               /* Blue */
#define MAGENTA "\033[35m"            /* Magenta */
#define CYAN "\033[36m"               /* Cyan */
#define WHITE "\033[37m"              /* White */
#define BOLDBLACK "\033[1m\033[30m"   /* Bold Black */
#define BOLDRED "\033[1m\033[31m"     /* Bold Red */
#define BOLDGREEN "\033[1m\033[32m"   /* Bold Green */
#define BOLDYELLOW "\033[1m\033[33m"  /* Bold Yellow */
#define BOLDBLUE "\033[1m\033[34m"    /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m" /* Bold Magenta */
#define BOLDCYAN "\033[1m\033[36m"    /* Bold Cyan */
#define BOLDWHITE "\033[1m\033[37m"   /* Bold White */

//! 单元/布局行的朝向：N/R90(W)/S/R180(E) 及其翻转（F 前缀 = flipped）
enum ORIENT
{
    OR_N,
    OR_W,
    OR_S,
    OR_E,
    OR_FN,
    OR_FW,
    OR_FS,
    OR_FE
    //         0    1    2    3    4     5     6     7
    //         |    -    |    -    |     -     |     -
};

//! 把 Bookshelf 里的朝向字符串（"N"/"FS" 等）转成 ORIENT 枚举值；无法识别时退回 N
inline int orientInt(char *c) // transform orientation string to int (enum ORIENT)
{
    assert(strlen(c) <= 2);
    if (c[0] == 'F')
    {
        switch (c[1])
        {
        case 'N':
            return 4;
        case 'W':
            return 5;
        case 'S':
            return 6;
        case 'E':
            return 7;
        }
    }
    else
    {
        switch (c[0])
        {
        case 'N':
            return 0;
        case 'W':
            return 1;
        case 'S':
            return 2;
        case 'E':
            return 3;
        }
    }
    return 0; // Use "N" if any problem occurs.
}
//! 2D 点 / 偏移量。既用来存坐标，也用来存 pin 相对单元的 offset
struct POS_2D // POS means postition, POS_2D can be used to store coordinates, offsets
{
    float x;
    float y;
    POS_2D() { SetZero(); };
    POS_2D(float _x, float _y)
    {
        x = _x;
        y = _y;
    }
    inline void SetZero()
    {
        x = y = 0.0; //!! 0.0!!!!
    }
    friend inline std::ostream &operator<<(std::ostream &os, const POS_2D &pos)
    {
        os << "(" << pos.x << "," << pos.y << ")";
        return os;
    }
};

//! 2D 向量（梯度、电场等）。语义上和 POS_2D 不同，所以单独定义一份类型
struct VECTOR_2D
{
    float x;
    float y;
    VECTOR_2D()
    {
        SetZero();
    }
    inline void SetZero()
    {
        x = y = 0.0; //!! 0.0!!!!
    }
    friend inline std::ostream &operator<<(std::ostream &os, const VECTOR_2D &vec)
    {
        os << "[" << vec.x << "," << vec.y << "]"; // [] for vectors and () for pos
        return os;
    }
};

//! 2D 整数向量，主要用于存 bin 网格下标（如 binStartIdx / binDimension）
struct VECTOR_2D_INT
{
    int x;
    int y;
    VECTOR_2D_INT()
    {
        SetZero();
    }
    inline void SetZero()
    {
        x = y = 0; //!! 0.0!!!!
    }
    friend inline std::ostream &operator<<(std::ostream &os, const VECTOR_2D_INT &vec)
    {
        os << "[" << vec.x << "," << vec.y << "]"; // [] for vectors and () for pos
        return os;
    }
};

//! 3D 向量，带加减、数乘、点积。布局里的梯度、电场、坐标都用它
struct VECTOR_3D
{
    float x;
    float y;
    float z;
    VECTOR_3D()
    {
        SetZero();
    }
    inline void SetZero()
    {
        x = y = z = 0.0; //!! 0.0!!!!
    }
    friend inline std::ostream &operator<<(std::ostream &os, const VECTOR_3D &vec)
    {
        os << "[" << vec.x << "," << vec.y << "," << vec.z << "]"; // [] for vectors and () for pos
        return os;
    }
    inline VECTOR_3D operator+(const VECTOR_3D &rhs)
    {
        VECTOR_3D v;
        v.x = this->x + rhs.x;
        v.y = this->y + rhs.y;
        v.z = this->z + rhs.z;
        return v;
    }
    inline VECTOR_3D operator-(const VECTOR_3D &rhs) const
    {
        VECTOR_3D v;
        v.x = this->x - rhs.x;
        v.y = this->y - rhs.y;
        v.z = this->z - rhs.z;
        return v;
    }
    inline VECTOR_3D operator*(float c)
    {
        VECTOR_3D v;
        v.x = this->x * c;
        v.y = this->y * c;
        v.z = this->z * c;
        return v;
    }
    inline float operator*(const VECTOR_3D &rhs) const
    {
        return (this->x * rhs.x + this->y * rhs.y + this->z * rhs.z);
    }
};

//! 3D 点，继承自 VECTOR_3D（复用运算），只是打印格式改成 ()
struct POS_3D : public VECTOR_3D
{
    // float x;
    // float y;
    // float z;
    POS_3D() { SetZero(); };
    POS_3D(float _x, float _y, float _z)
    {
        x = _x;
        y = _y;
        z = _z;
    }
    inline void SetZero()
    {
        x = y = z = 0.0; //!! 0.0!!!!
    }
    friend inline std::ostream &operator<<(std::ostream &os, const POS_3D &pos)
    {
        os << "(" << pos.x << "," << pos.y << "," << pos.z << ")";
        return os;
    }
};

//! 三分量 bool，用来标记 X/Y/Z 上某个指数项是否已下溢为 0（WA/LSE 模型用）
struct VECTOR_3D_BOOL
{
    bool x;
    bool y;
    bool z;
    VECTOR_3D_BOOL()
    {
        SetZero();
    }
    inline void SetZero()
    {
        x = y = z = false;
    }
    friend inline std::ostream &operator<<(std::ostream &os, const VECTOR_3D_BOOL &vec)
    {
        os << "b " << vec.x << "," << vec.y << "," << vec.z << " b"; // [] for vectors and () for pos
        return os;
    }
};

//! 轴对齐矩形（axis-aligned rectangle），用左下 ll + 右上 ur 表示。
//! coreRegion、bin、单元包围盒、布局行全都用它
class CRect
{
public:
    CRect()
    {
        Init();
    }
    void Print()
    {
        cout << "lower left: " << ll << " to upper right: " << ur << "\n";
    }
    void Init()
    {
        ll.SetZero();
        ur.SetZero();
    }
    POS_2D ll; // ll: lower left coor
    POS_2D ur; // ur: upper right coor
    float getWidth()
    {
        float width = ur.x - ll.x;
        assert(width > 0.0);
        return width;
    }
    float getHeight()
    {
        float height = ur.y - ll.y;
        assert(height > 0.0);
        return height;
    }
    POS_2D getCenter()
    {
        POS_2D center = ll;
        center.x += 0.5 * this->getWidth();
        center.y += 0.5 * this->getHeight();
        return center;
    }
    float getArea()
    {
        return getHeight() * getWidth();
    }
    bool inside(POS_2D &point)
    {
        return (point.x >= ll.x) && (point.x <= ur.x) && (point.y >= ll.y) && (point.y <= ur.y);
    }
    friend inline std::ostream &operator<<(std::ostream &os, const CRect &rect)
    {
        os << "CRect Size: " << rect.ur.x - rect.ll.x << "," << rect.ur.y - rect.ll.y << endl;
        return os;
    }
};

    //! 浮点乘/除/平方的简单包装。目前就是直接运算，历史上是为统一加溢出检查留的口子
    inline float float_mul(float a, float b) // a*b
{
    float c = a * b;
    return c;
}

inline float float_div(float a, float b) // a/b
{
    float c = a / b;
    return c;
}

inline float float_square(float a) // a^2
{
    return a * a;
}

//! 带容差的浮点比较族：布局里大量「是否比 bin 小」「面积是否为 0」这类判断，
//! 直接用 == 会被浮点误差坑，统一走这里（容差 EPS = 1e-15）
inline bool float_greater(float a, float b) // return true if a > b
{
    return a - b > 1.0 * EPS;
}

inline bool float_less(float a, float b) // return true if a < b
{
    return a - b < -1.0 * EPS;
}
inline bool float_equal(float a, float b)
{
    return std::fabs(a - b) < EPS;
}

inline bool float_lessorequal(float a, float b)
{
    return float_less(a, b) || float_equal(a, b);
}

inline bool float_greaterorequal(float a, float b)
{
    return float_greater(a, b) || float_equal(a, b);
}

inline POS_2D POS_2D_add(POS_2D a, POS_2D b)
{
    a.x += b.x;
    a.y += b.y;
    return a;
}

inline POS_2D POS_2D_scale(POS_2D a, float scaleFactor)
{
    a.x *= scaleFactor;
    a.y *= scaleFactor;
    return a;
}

inline double seconds()
{
    //! 返回进程的 CPU 时间（用户态 + 内核态），注意和下面 time_start/time_end 的墙钟时间口径不同
#ifdef WIN32 // Windows
    struct __timeb64 tstruct;
    _ftime64(&tstruct);
    return (double)tstruct.time + 0.001 * tstruct.millitm;
#else // Linux
    rusage time;
    getrusage(RUSAGE_SELF, &time);
    // return (double)(1.0*time.ru_utime.tv_sec+0.000001*time.ru_utime.tv_usec);	// user time

    return (double)(1.0 * time.ru_utime.tv_sec + 0.000001 * time.ru_utime.tv_usec + // user time +
                    1.0 * time.ru_stime.tv_sec + 0.000001 * time.ru_stime.tv_usec); // system time
#endif

    // clock() loop is about 72min. (or 4320 sec)
    // return double(clock())/CLOCKS_PER_SEC;
}

//! 记录当前墙钟时间到 *time_cost（gettimeofday，精度微秒）
inline void time_start(double *time_cost)
{
    struct timeval time_val;
    time_t time_secs;
    suseconds_t time_micro;
    gettimeofday(&time_val, NULL);
    time_micro = time_val.tv_usec;
    time_secs = time_val.tv_sec;
    *time_cost = (double)time_micro / 1000000 + time_secs;
    return;
}

//! 用当前墙钟时间减去 *time_cost，把「耗时」写回同一变量（原地复用，成对使用）
inline void time_end(double *time_cost)
{
    struct timeval time_val;
    time_t time_secs;
    suseconds_t time_micro;
    gettimeofday(&time_val, NULL);
    time_micro = time_val.tv_usec;
    time_secs = time_val.tv_sec;
    *time_cost = (double)time_micro / 1000000 + time_secs - *time_cost;
    return;
}

//! 段错误定位用：加 -segDebug 参数时，在每个阶段入口打印 checkpoint 名字，
//! 程序崩了就能从最后打印的 checkpoint 判断挂在哪个阶段
inline void segmentFaultCP(string checkpointname)
{
    if (gArg.CheckExist("segDebug"))
    {
        cout << endl
             << padding << checkpointname << padding << endl;
    }
}

//! 加 -debug 参数时才打印的标量调试输出
inline void debugOutput(string caption, double value)
{
    if (gArg.CheckExist("debug"))
    {
        cout << caption << ": " << value << endl;
    }
}

// ----------------------------------------------------------------------------
// 快速指数近似：利用 lim(n→∞)(1 + a/n)^n = e^a，取 n = 1024 = 2^10，
// 于是「自乘 10 次」就能得到 (1 + a/1024)^1024 ≈ e^a，只用了 1 次除法 + 10 次乘法。
//
// ePlace 每轮要给每个 pin 算若干次 e^[]，是热点中的热点，所以宁可牺牲精度换速度。
// 代价：只在 |a| 较小时准，a 很负时误差大——因此在指数小于 NEGATIVE_MAX_EXP 时
// 直接判定为 0（见 expZeroFlg*），既避免下溢也避开这个近似的失真区间。
// ----------------------------------------------------------------------------
inline float fastExp(float a)
{
    a = 1.0f + a / 1024.0f;
    a *= a;
    a *= a;
    a *= a;
    a *= a;
    a *= a;
    a *= a;
    a *= a;
    a *= a;
    a *= a;
    a *= a;
    return a;
}

//! 求两条线段 [x1,x2] 与 [x3,x4] 的重叠长度。无重叠时返回 0
//! 疑似问题：有重叠时返回的是 (overlapStart − overlapEnd)，符号是**负**的。
//!   单独用这个函数拿到的长度是负数；但 getOverlapArea 里两个负长度相乘正好得正面积，
//!   所以面积计算侥幸是对的，直接调用方则会被坑
inline double getOverlap(double x1, double x2, double x3, double x4) // two lines: x1->x2 and x3->x4
{
    assert(x1 <= x2);
    assert(x3 <= x4);

    // overlapStart: start point of the overlap line
    double overlapStart = max(x1, x3);
    double overlapEnd = min(x2, x4);

    if (overlapStart >= overlapEnd)
    {
        return 0;
    }
    else
    {
        return (overlapStart - overlapEnd);
    }
}

//! 两个矩形的重叠面积：X 方向重叠长度 × Y 方向重叠长度。
//! 任一方向无重叠就直接返回 0（提前退出，省一次计算）
//! 这是密度计算（单元摊到 bin）里调用最频繁的函数之一
inline double getOverlapArea(double left1, double bottom1, double right1, double top1,
                             double left2, double bottom2, double right2, double top2)
{
    assert(left1 <= right1);
    assert(bottom1 <= top1);
    assert(left2 <= right2);
    assert(bottom2 <= top2);

    double rangeH;
    rangeH = getOverlap(left1, right1, left2, right2);
    if (rangeH == 0)
        return 0;

    double rangeV;
    rangeV = getOverlap(bottom1, top1, bottom2, top2);
    if (rangeV == 0)
        return 0;

    return (rangeH * rangeV);
}

inline double getOverlapArea_2D(CRect rect1, CRect rect2)
{
    double left1 = rect1.ll.x;
    double bottom1 = rect1.ll.y;
    double right1 = rect1.ur.x;
    double top1 = rect1.ur.y;
    double left2 = rect2.ll.x;
    double bottom2 = rect2.ll.y;
    double right2 = rect2.ur.x;
    double top2 = rect2.ur.y;
    return getOverlapArea(left1, bottom1, right1, top1, left2, bottom2, right2, top2);
}

inline double getOverlapArea_2D(POS_2D ll1, POS_2D ur1, POS_2D ll2, POS_2D ur2)
{
    CRect rect1;
    CRect rect2;

    rect1.ll = ll1;
    rect1.ur = ur1;

    rect2.ll = ll2;
    rect2.ur = ur2;

    return getOverlapArea_2D(rect1, rect2);
}

//! 三维向量的 L2 范数
inline float L2NORM(VECTOR_3D a)
{
    return sqrt(float_square(a.x) + float_square(a.y) + float_square(a.z));
}

//! 快速排序的划分（partition）步骤，把 array[begin] 作为枢轴放到正确位置并返回其下标。
//! 原函数名 partion 是拼写错误（应为 partition）。
//! 疑似问题：`int key = array[begin]` 把 float 截成了 int，比较用的是整数化的枢轴，
//!   结果不正确（除非数组元素恰好是整数）
inline int partion(vector<float> &array, int begin, int end)
{
    int start = begin;
    int key = array[begin];
    while (begin < end)
    {
        while (begin < end && array[end] >= key)
        {
            end--;
        }
        while (begin < end && array[begin] <= key)
        {
            begin++;
        }
        swap(array[begin], array[end]);
    }
    swap(array[start], array[begin]);
    return begin;
}

//! 快速选择（quickselect）：平均 O(n) 求出第 K 小的元素（K 从 0 起），不改排序整体。
//! 注意：partion 会**原地修改** array 的元素顺序
//! 疑似问题：若 K 越界（>= size）会死循环
inline float getKth(vector<float> &array, int K) // return the Kth smallest number in an array (start from 0)
{
    int left = 0;
    int right = array.size() - 1;
    int index = -1;
    while (index != K)
    {
        index = partion(array, left, right);
        if (index > K)
        {
            right = index - 1;
        }
        else if (index < K)
        {
            left = index + 1;
        }
        else
        {
            break;
        }
    }
    return array[index];
}
#endif