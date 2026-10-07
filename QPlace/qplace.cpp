#include "qplace.h"
using namespace Eigen;

// ============================================================================
// QPPlacer —— 二次布局（Quadratic Placement），为 ePlace 生成初始布局
//
// 整体流程（Kraftwerk2 / ePlace 的经典做法）：
//   1. 把所有可动单元先堆到 core 区域中心（moveNodesCenterToCenter）；
//   2. 用 HPWL 的二次近似构造线性方程组 A·x = b，X / Y 两个方向各一套；
//      线长模型里只保留每个 net 的「边界 pin」（boundPinXmin/Xmax、Ymin/Ymax）
//      相关的 pin 对，其余 pin 对不贡献 HPWL，直接跳过以减小矩阵规模；
//   3. 用 BiCGSTAB 求解，把解写回单元坐标；
//   4. 位置变了 → 权重 w = 1/((p-1)·d) 也变了，于是重建矩阵再解，循环迭代。
//
// 重要前提：这里是「纯」二次布局，没有加任何扩散力 / 密度力（spreading），
// 可动单元只会朝所连 net 的加权重心聚拢，靠 terminal（固定 pin）把它撑开。
// 因此结果通常高度重叠，只作为初始解交给后续 ePlace 的全局布局去展开。
// ============================================================================

void QPPlacer::quadraticPlacement()
{
    // TESTING, one HPWL can be removed in the future
    db->moveNodesCenterToCenter(); //!!!!!! //! 第 1 步：所有可动单元搬到 core 中心，作为迭代起点

    double HPWL = db->calcHPWL(); //! 直接遍历所有 pin 求 HPWL

    double HPWL2 = db->calcNetBoundPins(); //! 顺带把每个 net 的 boundPin* 求出来（副作用），同时累加 HPWL

    assert(HPWL == HPWL2); //! 一致性校验：两种算法得到的线长必须相同（也保证 boundPin 已更新）

    int nodeCount = db->dbNodes.size(); //! 矩阵维度 = 可动单元数（不含 terminal）
    int maxIterationNumber;

    //! 迭代次数上限，命令行参数 -IPiteCount 指定（注意原名疑似拼写错误，应为 IPIterCount）
    if (gArg.CheckExist("IPiteCount"))
    {
        gArg.GetInt("IPiteCount", &maxIterationNumber);
    }
    else
    {
        maxIterationNumber = 20;
    }

    float xError, yError;                 //! BiCGSTAB 求解后返回的相对残差
    float target_error = 0.000001;        //! 收敛阈值

    printf("INFO:  initial HPWL: %.6lf\n", HPWL);
    cout << "INFO:  The matrix size is " << nodeCount << endl; //! matrix size: only nodes(movable modules)

    setNbThreads(8); //! a parameter(numThreads) fixed here //! 硬编码 8 线程，可考虑做成命令行参数

    // BCGSTAB settings
    //! A 是稀疏对称正定（半正定）矩阵；x 既是上一轮的解（作为初值 guess）也是本轮输出；b 是常数项
    SMatrix X_A(nodeCount, nodeCount), Y_A(nodeCount, nodeCount);            // A for Ax=b
    VectorXf X_x(nodeCount), Y_x(nodeCount), X_b(nodeCount), Y_b(nodeCount); // x and b for Ax=b
    double qp_time;

    //! 外层不动点迭代：解一次 → 更新位置 → 重算权重 → 再解，直到残差足够小或达到上限
    for (int i = 0;; i++)
    {

        time_start(&qp_time);
        //! 用当前坐标重建 A 与 b，同时把当前坐标填进 X_x / Y_x 当作求解初值
        createSparseMatrix(X_A, Y_A, X_x, Y_x, X_b, Y_b);

        //! 非对称/通用稀疏迭代求解器；这里用不上对称特性（A 理论上对称，但 BiCGSTAB 更稳）
        BiCGSTAB<SMatrix, IdentityPreconditioner> solver;
        solver.setMaxIterations(100); //! 内层迭代上限（注意：这是求解器内部迭代，和外层 for 的 i 不是一回事）

        //! X 方向
        solver.compute(X_A);
        X_x = solver.solveWithGuess(X_b, X_x); //! 带初值求解，初值 = 上一次的位置，收敛更快
        xError = solver.error();

        //! Y 方向
        solver.compute(Y_A);
        Y_x = solver.solveWithGuess(Y_b, Y_x);
        yError = solver.error();

        updateModuleLocation(X_x, Y_x); //! 把解向量写回每个 module 的中心坐标
        HPWL = db->calcNetBoundPins();  //! 重新统计线长（同时刷新 boundPin，供下一轮建矩阵用）

        if (gArg.CheckExist("debug") || gArg.CheckExist("fullPlot"))
        {
            PLOTTING::plotCurrentPlacement("Initial placement iteration " + to_string(i), db);
        }

        time_end(&qp_time);
        printf("INFO:  at iteration number %3d,  CG Error %.6lf,  HPWL %.6lf,  CPUtime %.2lf\n", i, max(xError, yError), HPWL, qp_time);

        //! 收敛判据：两个方向的残差都足够小。
        //! i > 4 是为了避开「假收敛」——首轮所有单元都在中心、pin 间距被 MIN_DISTANCE 夹住，
        //! 矩阵条件数特殊，残差可能一开始就很小，需要多跑几轮才算数。
        if (fabs(xError) < target_error && fabs(yError) < target_error && i > 4)
        {
            break;
        }
        if (i >= maxIterationNumber)
        {
            break;
        }
    }
}

// ----------------------------------------------------------------------------
// 构造 A·x = b：X / Y 两个方向各一套矩阵
//
// 数学形式（以 X 方向为例，pin 绝对坐标 = module 中心 + pin offset）：
//   对 net 上每个有效 pin 对 (i, j)，线长项贡献  w·(x_i + off_i − x_j − off_j)²
//   对 x_i 求偏导并令其为 0：
//       A[i][i] += w      A[j][j] += w
//       A[i][j] -= w      A[j][i] -= w
//       b[i]    += −w·(off_i − off_j)
//       b[j]    += −w·(off_j − off_i)
//   若其中一端是固定 terminal（例如 pin1 固定）：
//       A[2][2] += w
//       b[2]    += w·(X_fixed − off_2)      ← 固定端坐标并入常数项 b
//
//   权重 w = 1/((p−1)·d)，p 为该 net 的 pin 数，d 为两 pin 当前间距；
//   这是 Kraftwerk2 公式(8) 的形式，此处按 ePlace 的做法用 1 而不是 2 作分子。
//   间距 d 会被 MIN_DISTANCE 夹住下限，避免除零 / 权重爆炸。
// ----------------------------------------------------------------------------
void QPPlacer::createSparseMatrix(SMatrix &X_A, SMatrix &Y_A, VectorXf &X_x, VectorXf &Y_x, VectorXf &X_b, VectorXf &Y_b)
{
    vector<T> tripletListX, tripletListY; // create matrix A with triples, a triple：(i,j,v), means Aij=v
    //! 预留 1e7 个三元组：Triplet<float> 约 12 字节，两个 list 各占约 120MB，合计约 240MB。
    //! 在大设计上这是主要的内存开销，可按 net/pin 规模估算后动态调整。
    tripletListX.reserve(10000000);
    tripletListY.reserve(10000000);
    int nodeCount = db->dbNodes.size();

    //! 1) 初始化解向量与常数项：以当前坐标作为本轮求解的初值，b 清零
    Module *curNode = NULL;
    for (int i = 0; i < nodeCount; i++)
    {
        curNode = db->dbNodes[i];
        assert(curNode);
        assert(curNode->idx == i); //! dbNodes 的下标必须和 module->idx 一一对应，否则矩阵行列号就错了
        X_x(i) = curNode->getCenter().x; //!
        Y_x(i) = curNode->getCenter().y;
        X_b(i) = Y_b(i) = 0;
    }

    //! 2) 遍历所有 net，枚举 pin 对，往 A / b 里填贡献
    for (Net *curNet : db->dbNets)
    {
        assert(curNet);
        int pinCount = curNet->netPins.size();
        Pin *pin1;
        Pin *pin2;
        float constant1 = 1.0 / ((float)pinCount - 1.0); // see kraftwerk2 equation (8), here we follow eplace and use 1 instead of 2
        //! 注意：两点 net（pinCount==1，即单 pin net）会让 constant1 除零，实际数据里一般不会出现
        for (int j = 0; j < pinCount; j++)
        {
            pin1 = curNet->netPins[j];
            assert(pin1);
            for (int k = j + 1; k < pinCount; k++) //! 只枚举上三角，每个无序 pin 对处理一次
            {
                pin2 = curNet->netPins[k];
                assert(pin2);

                //! 同一个 module 上的两个 pin：距离为 0，对 HPWL 无贡献，跳过
                if (pin1->module == pin2->module)
                {
                    continue;
                }

                // ================= X 方向 =================
                //! 只有包含 boundPinXmin / boundPinXmax 的 pin 对才参与 X 方向线长项，
                //! 其余 pin 对不影响 HPWL（这是用二次项近似 HPWL 的关键简化）。
                //! boundPin* 由上一轮的 calcNetBoundPins() 刷新。
                if (pin1 == curNet->boundPinXmin || pin1 == curNet->boundPinXmax || pin2 == curNet->boundPinXmin || pin2 == curNet->boundPinXmax)
                {
                    float distanceX = fabs(pin1->getAbsolutePos().x - pin2->getAbsolutePos().x);

                    float weightX = 0.0f; // see kraftwerk2 equation (8), here we follow eplace and use 1 instead of 2
                    //! 间距夹下限：首轮所有单元叠在中心，d≈0 会让权重趋于无穷，
                    //! MIN_DISTANCE 同时起到「平滑」作用，让线长模型在早期更接近线性
                    if (float_greaterorequal(distanceX, MIN_DISTANCE))
                    {
                        weightX = constant1 / distanceX;
                    }
                    else
                    {
                        weightX = constant1 / MIN_DISTANCE;
                    }
                    if (weightX < 0)
                        printf("ERROR WEIGHT\n");

                    if (!pin1->module->isFixed && !pin2->module->isFixed) // both are movable modules
                    {
                        //! 两端都可动：标准的 2x2 块 [w −w; −w w] 摊到 A 的四个位置
                        tripletListX.push_back(T(pin1->module->idx, pin1->module->idx, weightX));
                        tripletListX.push_back(T(pin2->module->idx, pin2->module->idx, weightX));

                        tripletListX.push_back(
                            T(pin1->module->idx, pin2->module->idx, (-1.0) * weightX));
                        tripletListX.push_back(
                            T(pin2->module->idx, pin1->module->idx, (-1.0) * weightX));

                        //! offset 差异进 b：因为 A·x=b 里的 x 是 module 中心而非 pin 坐标
                        X_b(pin1->module->idx) += (-1.0) * weightX * ((pin1->offset.x) - (pin2->offset.x));
                        X_b(pin2->module->idx) += (-1.0) * weightX * ((pin2->offset.x) - (pin1->offset.x));
                    }

                    else if (pin1->module->isFixed && !pin2->module->isFixed) // 1 is terminal, 2 is movable
                    {
                        //! pin1 固定：只给可动端加对角项，固定端坐标并入 b（这就是把单元拉向 IO 的力）
                        tripletListX.push_back(T(pin2->module->idx, pin2->module->idx, weightX));
                        X_b(pin2->module->idx) += weightX * (pin1->getAbsolutePos().x - (pin2->offset.x));
                    }

                    else if (!pin1->module->isFixed && pin2->module->isFixed) // 2 is terminal, 1 is movable
                    {
                        //! 与上一种对称，仅固定端换成 pin2
                        tripletListX.push_back(T(pin1->module->idx, pin1->module->idx, weightX));
                        X_b(pin1->module->idx) += weightX * (pin2->getAbsolutePos().x - (pin1->offset.x));
                    }
                    //! 两端都固定的情况对求解无意义，不处理
                }

                // ================= Y 方向 =================
                //! 与 X 方向完全对称，只是换成 boundPinYmin / boundPinYmax 和 y 分量。
                //! ! 疑似 bug：条件里写了两遍 "pin2 == boundPinYmax"，却漏了 "pin1 == boundPinYmin"。
                //!   后果：当 boundPinYmin 恰好是 pin1（即它在 netPins 里下标更小）时，该 pin 对
                //!   不会进入 Y 方向线长项，Y 方向的矩阵会少一部分连接性。
                //!   X 方向的条件是四个都写全的，可以对照。此处保持原样，仅作标记。
                if (pin2 == curNet->boundPinYmin || pin1 == curNet->boundPinYmax || pin2 == curNet->boundPinYmin || pin2 == curNet->boundPinYmax)
                {
                    float distanceY = fabs(pin1->getAbsolutePos().y - pin2->getAbsolutePos().y);

                    float weightY = 0.0f; // see kraftwerk2 equation (8), here we follow eplace and use 1 instead of 2
                    if (float_greaterorequal(distanceY, MIN_DISTANCE))
                    {
                        weightY = constant1 / distanceY;
                    }
                    else
                    {
                        weightY = constant1 / MIN_DISTANCE;
                    }
                    if (weightY < 0)
                        printf("ERROR WEIGHT\n");

                    if (!pin1->module->isFixed && !pin2->module->isFixed) // both are movable modules
                    {
                        tripletListY.push_back(T(pin1->module->idx, pin1->module->idx, weightY));
                        tripletListY.push_back(T(pin2->module->idx, pin2->module->idx, weightY));

                        tripletListY.push_back(
                            T(pin1->module->idx, pin2->module->idx, (-1.0) * weightY));
                        tripletListY.push_back(
                            T(pin2->module->idx, pin1->module->idx, (-1.0) * weightY));

                        Y_b(pin1->module->idx) += (-1.0) * weightY * ((pin1->offset.y) - (pin2->offset.y));
                        Y_b(pin2->module->idx) += (-1.0) * weightY * ((pin2->offset.y) - (pin1->offset.y));
                    }

                    else if (pin1->module->isFixed && !pin2->module->isFixed) // 1 is terminal, 2 is movable
                    {
                        tripletListY.push_back(T(pin2->module->idx, pin2->module->idx, weightY));
                        Y_b(pin2->module->idx) += weightY * (pin1->getAbsolutePos().y - (pin2->offset.y));
                    }

                    else if (!pin1->module->isFixed && pin2->module->isFixed) // 2 is terminal, 1 is movable
                    {
                        tripletListY.push_back(T(pin1->module->idx, pin1->module->idx, weightY));
                        Y_b(pin1->module->idx) += weightY * (pin2->getAbsolutePos().y - (pin1->offset.y));
                    }
                }
            }
        }
    }

    //! 3) 三元组 → 稀疏矩阵。重复行列号的元素会被 Eigen 自动累加，所以上面可以放心 push 多个同位置项
    X_A.setFromTriplets(tripletListX.begin(), tripletListX.end());
    Y_A.setFromTriplets(tripletListY.begin(), tripletListY.end());

    //! 注意：若某个可动单元一个 net 都没连上（孤点），A 的对应行全为 0，
    //! 矩阵奇异，BiCGSTAB 解出的坐标可能是 NaN，需要靠数据保证每个单元都连线。
}

//! 把解向量 X_x / Y_x 写回 module 的中心坐标（只动 x、y，z 不变）
void QPPlacer::updateModuleLocation(VectorXf &X_x, VectorXf &Y_x)
{
    //! do indexs match?
    int nodeCount = db->dbNodes.size();
    Module *curNode;
    for (int i = 0; i < nodeCount; i++)
    {
        curNode = db->dbNodes[i];
        assert(curNode);
        assert(curNode->idx == i); //! 保证解向量的第 i 个分量就对应 dbNodes[i]
        db->setModuleCenter_2D(curNode, X_x(i), Y_x(i)); //! setModuleCenter_2D 会同步刷新该 module 上所有 pin 的绝对坐标
    }
}
