#ifndef NESTEROV_HPP
#define NESTEROV_HPP

// ============================================================================
// nesterov.hpp —— Nesterov 加速梯度优化器（ePlace 的全局布局求解引擎）
//
// Nesterov 加速梯度（NAG）维护两组解：
//   · main_solution      x_k ：真正的解
//   · reference_solution y_k ：外推点（在 x_k 基础上沿「上一步的移动方向」再外推一段）
// 每步：
//   1. 在 y_k 处求梯度（而不是在 x_k 处，这就是 Nesterov 与普通动量法的区别）
//      y_k = x_k + (a_k − 1)/a_{k+1} · (x_k − x_{k−1})
//   2. x_{k+1} = y_k − step · ∇f(y_k)
//   3. 动量系数递推 a_{k+1} = (1 + √(1 + 4·a_k²)) / 2
//
// 步长策略三选一（由命令行参数切换）：
//   · vanilla（默认） 1/L，L 用最近两点的梯度差/位置差在线估计
//   · bb    （-bb）   Barzilai-Borwein 步长，取长短两种 BB 步长里更稳的那个
//   · bktrk （-bktrk）回溯线搜索：先试一步，若 Lipschitz 估计变差就缩小步长重试
//
// 重要约定：placer->getGradient() 返回的**不是梯度本身**，而是已经取好负号、
// 乘好 preconditioner 的「本步应当移动的方向与距离」，所以更新式写作
// `新位置 = 旧位置 + gradient × step_size`（加号，不是减号）。
// ============================================================================

#include <vector>
#include <cstdint>
#include <cstddef>

#include "global.h"
#include "eplace.h"
#include "plot.h"
#include "optutils.hpp"
#include "Eigen/Dense"

#define MAX_ITERATION 1000  //! 迭代轮数硬上限，防止不收敛时死循环
#define BKTRK_EPS 0.95      //! 回溯线搜索的接受阈值：新步长不小于旧步长的 95% 就收

//! 一次迭代所需的三份数据：主解、外推解、梯度
template <typename T>
class NSIter
{
public:
    void resize(size_t length);
    std::vector<T> main_solution;
    std::vector<T> reference_solution;
    std::vector<T> gradient;
};

template <typename T>
void NSIter<T>::resize(size_t length)
{
    main_solution.resize(length);
    reference_solution.resize(length);
    gradient.resize(length);
}

//! ePlace 专用的 Nesterov 优化器，通过 placer 指针与 EPlacer_2D 交互
template <typename T>
class EplaceNesterovOpt : public FirstOrderOptimizer<T>
{
public:
    EplaceNesterovOpt(EPlacer_2D *placer) : placer(placer) {};
    // void opt();
private:
    bool stop_condition();
    void opt_step();
    void init();
    void wrap_up();
    void opt_step_bb();
    void opt_step_bktrk();
    void opt_step_vanilla();
    EPlacer_2D *placer;
    NSIter<T> cur_iter, last_iter;
    size_t iter_count;
    float NS_opt_param; // ak //! Nesterov 动量系数 a_k
};

// ----------------------------------------------------------------------------
// 收敛判据：不同阶段的阈值不同
//   mGP        ：溢出率 τ 降到 targetOverflow（默认 0.10）以下，或超过 MAX_ITERATION
//                退出时把迭代轮数记到 placer->mGPIterationCount，供 switch2cGP 回调 λ 用
//   FILLERONLY ：固定跑 20 轮（只是重撒 filler，不需要严格收敛）
//   cGP        ：τ 降到 0.07 以下（比 mGP 更严，毕竟只剩单元了）
// ----------------------------------------------------------------------------
template <typename T>
bool EplaceNesterovOpt<T>::stop_condition()
{
    float targetOverflow;
    if (gArg.CheckExist("targetOverflow"))
    {
        gArg.GetFloat("targetOverflow", &targetOverflow);
    }
    else
    {
        targetOverflow = 0.10f;
    }

    float cGPtargetOverflow = 0.07f;

    bool judge;
    switch (placer->placementStage)
    {
    case mGP:
        judge = (placer->globalDensityOverflow < targetOverflow) || (iter_count > MAX_ITERATION);
        if (judge)
        {
            placer->mGPIterationCount = iter_count;
        }
        // return (placer->globalDensityOverflow < targetOverflow) || (iter_count > MAX_ITERATION);
        return judge;
        break;
    case FILLERONLY:
        return iter_count >= 20;
        break;
    case cGP:
        return (placer->globalDensityOverflow < cGPtargetOverflow) || (iter_count > MAX_ITERATION);
        break;
    default:
        cerr << "INCORRECT PLACEMENT STAGE!\n";
        exit(0);
    }

    // if (placer->placementStage == mGP)
    // {
    //     bool judge = (placer->globalDensityOverflow < targetOverflow) || (iter_count > MAX_ITERATION);
    //     if (judge)
    //     {
    //         placer->mGPIterationCount = iter_count;
    //     }
    //     // return (placer->globalDensityOverflow < targetOverflow) || (iter_count > MAX_ITERATION);
    //     return judge;
    // }
    // else if (placer->placementStage == FILLERONLY)
    // {
    //     return iter_count >= 20;
    // }
    // else if (placer->placementStage == cGP)
    // {
    //     return (placer->globalDensityOverflow < cGPtargetOverflow) || (iter_count > MAX_ITERATION);
    // }
    // else
    // {
    //     cerr << "INCORRECT PLACEMENT STAGE!\n";
    //     exit(0);
    // }

    // return (placer->globalDensityOverflow < targetOverflow) || (iter_count > MAX_ITERATION);
}

template <typename T>
void EplaceNesterovOpt<T>::opt_step()
{
    printf("Iter %d\n", iter_count);

    //! 按命令行参数选择步长策略：bb / bktrk / 默认 vanilla
    if (gArg.CheckExist("bb"))
    {
        opt_step_bb();
    }
    else if (gArg.CheckExist("bktrk"))
    {
        opt_step_bktrk();
    }
    else
    {
        opt_step_vanilla();
    }
    //! 每步之后都要：根据线长变化调整惩罚因子 λ，并打印本轮指标
    placer->updatePenaltyFactor();
    placer->showInfo();

    //! 加 -fullPlot 时按阶段周期性出图（FILLERONLY 阶段每轮都出，方便观察 filler 扩散）
    switch (placer->placementStage)
    {
    case mGP:
        if (iter_count % 10 == 0 && gArg.CheckExist("fullPlot"))
        {
            PLOTTING::plotEPlace_2D("mGP Iter-" + to_string(iter_count), placer);
        }
        break;
    case FILLERONLY:
        if (gArg.CheckExist("fullPlot"))
        {
            PLOTTING::plotEPlace_2D("FILLERONLY Iter-" + to_string(iter_count), placer);
        }
        break;
    case cGP:
        if (iter_count % 10 == 0 && gArg.CheckExist("fullPlot"))
        {
            PLOTTING::plotEPlace_2D("cGP Iter-" + to_string(iter_count), placer);
        }
        break;
    default:
        cerr << "INCORRECT PLACEMENT STAGE!\n";
        exit(0);
    }

    iter_count++;
}

// template <typename T>
// void EplaceNesterovOpt<T>::reset(){
//     iter_count = 0;
// }

template <typename T>
void EplaceNesterovOpt<T>::init()
{
    iter_count = 0;
    NS_opt_param = 1; //! a_0 = 1
    cur_iter.main_solution = placer->getPosition(); //! 以当前布局作为优化起点
    printf("main length %d\n", cur_iter.main_solution.size());
}

template <typename T>
void EplaceNesterovOpt<T>::wrap_up()
{
    //! 把最终解写回数据库。
    //! 疑似问题：循环里每步 setPosition 写的是 reference_solution，而这里写的是 main_solution，
    //!   两者在最后一步并不相等（相差一个动量外推项），最后一次写回会引入一个额外位移
    placer->setPosition(cur_iter.main_solution);
}

// ----------------------------------------------------------------------------
// Barzilai-Borwein 步长版本
//   BB 法把「拟牛顿」思想简化成标量步长：用最近一步的 s = Δx、y = Δg 估计 Hessian 的标量近似
//     长步长：sᵀs / sᵀy        短步长：sᵀy / yᵀy
//   短步长通常更稳定，优先取它；若非正（说明曲率方向异常）则退回 min(‖s‖/‖y‖, 1/L) 兜底。
// ----------------------------------------------------------------------------
template <>
void EplaceNesterovOpt<VECTOR_3D>::opt_step_bb()
{
    cur_iter.reference_solution = placer->getPosition();
    placer->totalGradientUpdate(); //?
    cur_iter.gradient = placer->getGradient();
    float step_size;
    NSIter<VECTOR_3D> new_iter;
    size_t length = cur_iter.main_solution.size();
    new_iter.resize(length);
    if (iter_count == 0)
    {
        // step_size = 0.01;
        step_size = 1; //! 首轮没有历史信息，先迈一大步
    }
    else
    {
        float lipschitz_constant = calc_lipschitz_constant(cur_iter.reference_solution, last_iter.reference_solution, cur_iter.gradient, last_iter.gradient);
        float lip_step = 1 / lipschitz_constant;
        Eigen::VectorXf s(3 * length);
        Eigen::VectorXf y(3 * length);
        // move them into eigen vector
        for (size_t i = 0; i < length; i++)
        {
            s[3 * i] = cur_iter.reference_solution[i].x - last_iter.reference_solution[i].x;
            s[3 * i + 1] = cur_iter.reference_solution[i].y - last_iter.reference_solution[i].y;
            s[3 * i + 2] = cur_iter.reference_solution[i].z - last_iter.reference_solution[i].z;
            y[3 * i] = cur_iter.gradient[i].x - last_iter.gradient[i].x;
            y[3 * i + 1] = cur_iter.gradient[i].y - last_iter.gradient[i].y;
            y[3 * i + 2] = cur_iter.gradient[i].z - last_iter.gradient[i].z;
        }
        float s_norm = s.norm();
        float y_norm = y.norm();
        float s_dot_y = s.dot(y);
        float bb_long_step = pow(s_norm, 2) / s_dot_y;
        float bb_short_step = s_dot_y / pow(y_norm, 2);
        if (bb_short_step > 0)
        {
            step_size = bb_short_step;
        }
        else
        {
            step_size = min(s_norm / y_norm, lip_step);
        }
        printf("lip=%f,bb_short=%f,bb_long=%f\n", lip_step, bb_short_step, bb_long_step);
    }
    printf("stepSize : %f\n", step_size);

    //! Nesterov 动量系数递推 a_{k+1} = (1 + √(1 + 4a_k²)) / 2
    float new_NS_opt_param = (1 + sqrt(4 * float_square(NS_opt_param) + 1)) / 2; // ak+1

    // perform one step
    //! 核心两步（三个 opt_step_* 版本共用这一段）：
    //!   x_{k+1} = y_k + step · dir        （dir 是已取负号的移动方向，所以这里是加号）
    //!   y_{k+1} = x_{k+1} + (a_k−1)/a_{k+1} · (x_{k+1} − x_k)   ← 动量外推
    for (size_t idx = 0; idx < length; idx++)
    {
        VECTOR_3D &new_position = new_iter.main_solution[idx];
        VECTOR_3D &new_reference_position = new_iter.reference_solution[idx];

        VECTOR_3D gradient = cur_iter.gradient[idx];
        VECTOR_3D cur_position = cur_iter.main_solution[idx];
        VECTOR_3D cur_reference_position = cur_iter.reference_solution[idx];
        new_position = cur_reference_position + gradient * step_size;
        new_reference_position = new_position + (new_position - cur_position) * ((NS_opt_param - 1) / new_NS_opt_param);
    }
    //! 注意：写回的是外推点 reference_solution —— 下一步就在外推点求梯度，这才是 Nesterov 的精髓
    placer->setPosition(new_iter.reference_solution);
    last_iter = cur_iter;
    cur_iter = new_iter;
    NS_opt_param = new_NS_opt_param;
}

// ----------------------------------------------------------------------------
// 回溯线搜索版本：先按 1/L 试一步，用新位置重新算一次 L；
// 若新步长相比旧步长没有明显退化（≥ BKTRK_EPS × 旧步长）就接受，否则缩小再试。
// 代价：每轮要多算一次梯度（totalGradientUpdate），比 vanilla 慢，但更稳。
// 疑似问题：循环里没有重试次数上限，若条件一直不满足会死循环
// ----------------------------------------------------------------------------
template <>
void EplaceNesterovOpt<VECTOR_3D>::opt_step_bktrk()
{
    cur_iter.reference_solution = placer->getPosition();
    placer->totalGradientUpdate(); //?
    cur_iter.gradient = placer->getGradient();
    float step_size;
    NSIter<VECTOR_3D> new_iter;
    size_t length = cur_iter.main_solution.size();
    new_iter.resize(length);

    // initial step size
    if (iter_count == 0)
    {
        // step_size = 0.01;
        step_size = 1;
    }
    else
    {
        float lipschitz_constant = calc_lipschitz_constant(cur_iter.reference_solution, last_iter.reference_solution, cur_iter.gradient, last_iter.gradient);
        step_size = 1 / lipschitz_constant;
    }

    printf("stepSize : %f\n", step_size);

    float new_NS_opt_param = (1 + sqrt(4 * float_square(NS_opt_param) + 1)) / 2; // ak+1

    while (true)
    {
        // first move cells
        for (size_t idx = 0; idx < length; idx++)
        {
            VECTOR_3D &new_position = new_iter.main_solution[idx];
            VECTOR_3D &new_reference_position = new_iter.reference_solution[idx];

            VECTOR_3D gradient = cur_iter.gradient[idx];
            VECTOR_3D cur_position = cur_iter.main_solution[idx];
            VECTOR_3D cur_reference_position = cur_iter.reference_solution[idx];
            new_position = cur_reference_position + gradient * step_size;
            new_reference_position = new_position + (new_position - cur_position) * ((NS_opt_param - 1) / new_NS_opt_param);
        }
        //! 真走一步，再用新位置的梯度重估 L
        placer->setPosition(new_iter.reference_solution);
        placer->totalGradientUpdate();
        new_iter.gradient = placer->getGradient();
        float new_lipschitz_constant = calc_lipschitz_constant(new_iter.reference_solution, cur_iter.reference_solution, new_iter.gradient, cur_iter.gradient);
        float new_step_size = 1 / new_lipschitz_constant;
        if (BKTRK_EPS * step_size <= new_step_size)
        {
            step_size = new_step_size;
            break; //! 步长没有明显退化 → 接受这一步
        }
        printf("bktrk!\n");
        step_size = new_step_size; //! 缩小步长，下一轮循环重走
    }
    // update iter variable
    last_iter = cur_iter;
    cur_iter = new_iter;
    NS_opt_param = new_NS_opt_param;
}

// ----------------------------------------------------------------------------
// 最朴素的 Nesterov：步长固定取 1/L（L 由最近两点的梯度差/位置差在线估计）。
// 这是默认路径（不加 -bb / -bktrk 时走这里）。
// ----------------------------------------------------------------------------
template <>
void EplaceNesterovOpt<VECTOR_3D>::opt_step_vanilla()
{
    cur_iter.reference_solution = placer->getPosition();
    placer->totalGradientUpdate(); //?
    cur_iter.gradient = placer->getGradient();
    float step_size;
    NSIter<VECTOR_3D> new_iter;
    size_t length = cur_iter.main_solution.size();
    new_iter.resize(length);
    if (iter_count == 0)
    {
        // step_size = 0.01;
        step_size = 1;
    }
    else
    {
        float lipschitz_constant = calc_lipschitz_constant(cur_iter.reference_solution, last_iter.reference_solution, cur_iter.gradient, last_iter.gradient);
        step_size = 1 / lipschitz_constant;
    }
    printf("stepSize : %f\n", step_size);

    float new_NS_opt_param = (1 + sqrt(4 * float_square(NS_opt_param) + 1)) / 2; // ak+1

    // perform one step
    for (size_t idx = 0; idx < length; idx++)
    {
        VECTOR_3D &new_position = new_iter.main_solution[idx];
        VECTOR_3D &new_reference_position = new_iter.reference_solution[idx];

        VECTOR_3D gradient = cur_iter.gradient[idx];
        VECTOR_3D cur_position = cur_iter.main_solution[idx];
        VECTOR_3D cur_reference_position = cur_iter.reference_solution[idx];
        new_position = cur_reference_position + gradient * step_size;
        new_reference_position = new_position + (new_position - cur_position) * ((NS_opt_param - 1) / new_NS_opt_param);
    }
    placer->setPosition(new_iter.reference_solution);
    last_iter = cur_iter;
    cur_iter = new_iter;
    NS_opt_param = new_NS_opt_param;
}

#endif