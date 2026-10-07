#ifndef QPLACE_H
#define QPLACE_H
//! qplace: quadratic placement, used to create an initial placement for ePlace
//! 二次布局：解 A·x = b 得到单元坐标。只做线长最小化，不含密度/扩散约束，
//! 因此结果是重叠的初始解，后续交给 ePlace 的全局布局（含密度力）展开。
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers> //! BiCGSTAB 在这个 unsupported 模块里
#include "global.h"
#include "placedb.h"
#include "plot.h"
#include "arghandler.h"

using Eigen::VectorXf;
typedef Eigen::SparseMatrix<float, Eigen::RowMajor> SMatrix; //! 行主序稀疏矩阵，float 精度换取更快的求解
typedef Eigen::Triplet<float> T;                             //! (row, col, value) 三元组，用于拼装稀疏矩阵

//! 权重计算时 pin 间距的下限（防止两 pin 重合导致权重除零/爆炸）。
//! 同时起到平滑作用：该值越大，二次线长在近距离处越接近线性，早期迭代更稳定。
#define MIN_DISTANCE 25.0 /* 10.0 */ /* 5.0 */ /* 1.0 */

class QPPlacer;

class QPPlacer
{
public:
    QPPlacer()
    {
        db = NULL;
    }
    QPPlacer(PlaceDB *_db)
    {
        db = _db;
    }
    PlaceDB *db; //! 非拥有指针，生命周期由调用方管理

    //! 主入口：迭代「建矩阵 → 解 A·x=b → 更新坐标」，直到残差收敛或达到 IPiteCount 次
    void quadraticPlacement();

    //! 构造 X / Y 两个方向的 A、x、b，并顺带把当前坐标填入 x 作为求解初值
    void createSparseMatrix(SMatrix &X_A, SMatrix &Y_A, VectorXf &X_x, VectorXf &Y_x, VectorXf &X_b, VectorXf &Y_b); //! solve Ax=b for X and Y coordinates. Create A,x and b, see kraftwerk2 for more details
    //! 把解向量写回 module 中心坐标
    void updateModuleLocation(VectorXf &X_x, VectorXf &Y_x);
};
#endif
