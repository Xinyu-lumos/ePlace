#ifndef OPT_UTILS_HPP
#define OPT_UTILS_HPP

#include <vector>
#include <cstddef>
#include <assert.h>
#include "global.h"
/*
    T should override + - * / operators
*/

// ----------------------------------------------------------------------------
// 估计 Lipschitz 常数 L：‖∇f(x_t) − ∇f(x_{t−1})‖ / ‖x_t − x_{t−1}‖
//
// 用途：一阶优化里步长上界通常取 1/L。ePlace 的目标函数（尤其是带 γ 的平滑线长）
// 的 L 没有解析表达式，就用最近两点的「梯度差 / 位置差」来在线估计。
// 这也是 RePlAce / ePlace 论文里动态调整步长的标准做法。
// ----------------------------------------------------------------------------
template <typename T>
float calc_lipschitz_constant(const std::vector<T> &xt, const std::vector<T> &xt_1, const std::vector<T> &dxt, const std::vector<T> &dxt_1)
{
    assert(xt.size() == xt_1.size());
    assert(dxt.size() == dxt_1.size());
    float squared_sum_denominator = 0;
    float squared_sum_numerator = 0;
    for (size_t idx = 0; idx < xt.size(); idx++)
    {
        //! 分子：梯度变化的平方和；分母：位置变化的平方和
        squared_sum_numerator += (dxt[idx] - dxt_1[idx]) * (dxt[idx] - dxt_1[idx]);
        squared_sum_denominator += (xt[idx] - xt_1[idx]) * (xt[idx] - xt_1[idx]);
    
    }
    //! 疑似问题：分母为 0（两步位置完全没变）时会得到 inf/nan
    return sqrt(squared_sum_numerator) / sqrt(squared_sum_denominator);
}

#endif