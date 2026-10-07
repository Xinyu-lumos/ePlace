#ifndef OPT_HPP
#define OPT_HPP

// ============================================================================
// opt.hpp —— 一阶优化器的抽象基类（模板方法模式）
//
// 基类把「优化主循环」固定下来：init → while(!stop_condition) opt_step → wrap_up，
// 四个步骤是纯虚函数，由子类实现。这样换优化算法（Nesterov / Adam / Momentum）
// 只需替换子类，主循环逻辑不用重写。
//
// 模板参数 T 是「位置/梯度」的元素类型，本工程实际用的是 VECTOR_3D。
// ============================================================================

#include <vector>
#include <cstdint>
#include <cstddef>

template<typename T>
class FirstOrderOptimizer{
public:
    void opt();
private:
    virtual void init()=0;            //! 初始化：分配解向量、记录初始位置等
    virtual void wrap_up()=0;         //! 收尾：把最终解写回数据库
    virtual void opt_step()=0;        //! 单步更新：算梯度、走一步
    virtual bool stop_condition()=0;  //! 收敛判据：返回 true 则停止迭代
};

template <typename T>
void FirstOrderOptimizer<T>::opt(){
    init();
    while(!stop_condition()){
        opt_step();
    }
    wrap_up();
}

#endif