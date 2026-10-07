// ============================================================================
// 模块总览：DetailedPlacement（详细布局）— 本文件为头文件 / 数据结构声明
// ============================================================================
// 【职责】
//   DetailedPlacement 是整个布局流程（QPlace 二次布局 -> EPlace 全局布局 ->
//   Legalization 合法化）的最后一环。它的输入是一个**已经合法**的布局：
//   单元之间无重叠、都对齐到 site、都落在可放置区间（SiteRow 的 intervals）内。
//   它的目标是在**不破坏合法性**的前提下，通过局部的「移动 / 重排 / 交换」
//   进一步降低 HPWL（半周长线长）。因此它做的是小步长的局部搜索（local search），
//   而不是像 EPlace 那样的全局解析式（analytic）优化。
//
// 【整体流程】见 DetailedPlacer::detailedPlacement()：
//   1) initialization()：调用 placedb->removeBlockedSite()，把 macro / terminal
//      占用的 site 从 SiteRow 的 intervals 中剔除，得到「真正可放单元」的区间。
//      后续所有 row / segment / space 结构都以这份 intervals 为基础构建。
//   2) 外层循环 2 轮 { runLocalReordering(); runISM(); runGlobalSwap(); }
//      三种算子互补，反复迭代直到收益收敛。
//
// 【三种算子与各自的评估方式】
//   A. LocalReordering（局部重排，等价 branch-and-bound cell swap）
//      - 把每一行的每个可放置区间（LRSegment，类似 abacus 的一行）上的单元
//        按 x 排序，取一个大小为 windowSize 的滑动窗口（相邻窗口重叠 overlapSize）。
//      - 窗口内单元的**相对顺序**可以任意打乱，但总宽度不变、窗口跨度不变，
//        所以任意排列都是合法的：这保证了「搜索空间内全是合法解」。
//      - 用 DFS + 分支限界（LRSolver::depthFirstSolve）枚举排列，
//        代价用 net HPWL 增量式累加（LRSolution::solutionCost，对应 ntuplace3 的
//        m_bound）：只有当某条 net 上的窗口内单元**全部**放好之后，
//        这条 net 的 HPWL 才被计入（所以代价是单调不减的下界，可用于剪枝）。
//      - 原始顺序作为初始 bestSolution（上界），搜索结束后把窗口内单元
//        写回最佳排列的位置。
//   B. ISM（Independent Set Matching，独立集匹配）
//      - 借鉴 ntuplace3 的 grid_run()。用一个正方形的滑动窗口扫过整个芯片
//        （ISMSweep），窗口边长与重叠都按迭代轮次变化。
//      - 在每个窗口内挑出一组**互相没有连接关系**的单元（independent set），
//        关键在于：给它们各自分配一个**宽度都等于窗口内最大单元宽度 maxWIDTH**
//        的 slot。因为所有 slot 同宽，任意单元放到任意 slot 都不可能重叠 ->
//        **任意匹配结果都天然合法**，从而把「带约束的布局问题」转成
//        「无约束的二部图最小权匹配（LAP）」。
//      - 代价矩阵 cost[i][j] = 单元 i 放到 slot j 时的 HPWL，
//        用 lap2（Jonker-Volgenant 风格的 LAP 求解器）求全局最优指派。
//      - 独立集的意义：互相无连接的单元，其 HPWL 贡献彼此独立，
//        因此「求和最小 = 真实总代价最小」，匹配才是精确的。
//   C. GlobalSwap（全局交换）
//      - 论文 "An efficient and effective detailed placement algorithm" 的思路。
//      - 对每个单元计算它的 optimal region（由相连 net 的 bounding box 决定），
//        只在该区域内找候选位置：空白 slot，或**等宽**的单元做交换。
//      - 不找最优候选，只找「第一个有收益的」候选（first-fit），以省时间。
//      - 等宽交换同样是为了保证合法性：两个同宽单元对调，位置一定仍然合法。
//
// 【合法性是如何保证的】（本模块的核心设计）
//   三种算子都不是「先移动再修复」，而是**构造性地只产生合法解**：
//     - LocalReordering：只换顺序、不换总宽度 -> 不重叠。
//     - ISM：所有 slot 等宽（maxWIDTH）-> 任意指派都不重叠。
//     - GlobalSwap：只与空白或等宽单元交换 -> 不重叠。
//   所有 slot 的 x/y 都来自 rowSpaces（空闲区间）或已有单元的左下角坐标，
//   且都对齐到 site step / row height，因此也就自动满足对齐约束。
//
// 【核心数据结构】
//   - ISMRow：一行的抽象。rowSpaces（map<x起点, 长度>）维护空闲区间，
//     rowModules（map<x起点, Module*>）维护已放单元。ISM 与 GlobalSwap 共用。
//   - LRSegment：一行中的一个可放置区间（由 SiteRow::intervals 生成），
//     带 bottom / start / end 和 segModules（按 x 排序）。
//   - LRSolution：LocalReordering 的一个（部分/完整）解，
//     insertedCells / uninsertedCells / xLocations / netModuleCount / solutionCost。
//   - LRSolutionIterator：由部分解生成后继解（把下一个未放单元放到 currentX 处）。
//   - LRSolver：DFS + 分支限界的搜索器，维护 bestSolution / bestCost。
//   - lap2：线性指派问题（LAP）求解器，cost 矩阵 + assignment 结果。
//
// 【主要入口函数】
//   DetailedPlacer::detailedPlacement()        —— 顶层入口
//   ISMDP::ISMSweep() / ISMDP::ISMRun()        —— ISM 的滑动窗口与单次匹配
//   LocalReorderingDP::solve()                 —— 局部重排主循环
//   GlobalSwapDP::solve()                      —— 全局交换主循环
//   ISMRow::insertModule / removeModule / removeTail / insertSpace —— 区间维护原语
//   注意：本目录并未被根 CMakeLists.txt 的 add_subdirectory 包含，
//         因此当前不参与构建，代码可能与实际构建版本存在偏差。
// ============================================================================

#ifndef DETAILED_H
#define DETAILED_H
#include <global.h>
#include "placedb.h"
// 以下是本文件内部各组件的前向声明（它们的定义在下方依次展开）
class ISMDP;              // ISM：独立集匹配式详细布局
class ISMRow;             // 一行的空闲区间 / 已放单元的维护结构（ISM 与 GlobalSwap 共用）
class LocalReorderingDP;  // 局部重排（branch and bound cell swap）的驱动器
class LRSegment;          // 一行中的一个可放置区间，类似 abacus 的一行
class LRSolution;         // 局部重排中的一个（部分或完整）排列解
class LRSolutionIterator; // 由部分解枚举后继解的迭代器
class LRSolver;           // DFS + 分支限界求解器
class GlobalSwapDP;       // 全局交换（在 optimal region 内找空白 / 等宽单元交换）


// ============================================================================
// DetailedPlacer：详细布局的顶层驱动器
// ============================================================================
// 只负责「编排」：把三种局部搜索算子按固定顺序反复调用，自身不含算法细节。
// 三种算子的分工与互补关系：
//   - LocalReordering 处理**同一行内相邻少数单元**的顺序，粒度最细；
//   - ISM 处理**一个方形窗口内、跨若干行**的一批单元，用最优匹配重排；
//   - GlobalSwap 允许单元**跨较远距离**移动到它的 optimal region，
//     弥补前两者「只能在局部窗口内动」的局限。
// 因此先做细粒度的行内整理，再做窗口级匹配，最后做跨区域的全局交换，
// 并整体迭代若干轮，让三种算子的收益互相叠加。
class DetailedPlacer
{
public:
    DetailedPlacer()
    {
        Init();
    }
    DetailedPlacer(PlaceDB *_db)
    {
        Init();
        placedb = _db;
    }
    void Init()
    {
        placedb = NULL;
    }
    PlaceDB *placedb;

    //! 顶层入口：初始化后反复调用三种算子（详见类上方的总览说明）
    void detailedPlacement();

private:
    //! 准备工作：刷新 SiteRow 的 intervals，剔除被 macro / terminal 挡住的 site，
    //! 使后续所有 row / segment 都建立在「真正可放置」的区间上
    void initialization();

    //! 驱动 ISM（独立集匹配），内部带迭代收敛与时间上限控制
    void runISM();

    //! 驱动 GlobalSwap（在 optimal region 内做空白插入 / 等宽单元交换）
    void runGlobalSwap();

    void runLocalReordering(); // branch and bound cell swap = local reordering
};

// ============================================================================
// ISMDP：Independent Set Matching based Detailed Placement（独立集匹配）
// ============================================================================
// 核心思想（源自 ntuplace3 的 grid_run()）：
//   在一个正方形窗口内，选出一组**互不相连**的单元，并为每个单元准备一个
//   宽度相同（都等于窗口内最大单元宽度 maxWIDTH）的 slot。因为 slot 等宽，
//   「哪个单元放哪个 slot」就变成一个**无约束的二部图最小权匹配（LAP）**，
//   可以用 lap2 精确求最优解，且任意指派结果都自动合法（不会重叠）。
//   同时，因为单元间互不相连，各自对 HPWL 的贡献相互独立，
//   匹配最小化「各单元代价之和」才等价于最小化「窗口真实总线长」。
// 两个可切换的增强开关：
//   - doubleWindow      ：除当前窗口外，再向右取第 5~6 个窗口宽处的单元一起参与，
//                         扩大匹配规模（对应 ntuplace 的 double window 技巧）。
//   - independentCells  ：强制要求被选中的单元两两之间无共同 net（真正构成独立集）；
//                         关闭时匹配只是近似（代价不再可加），但能选到更多单元。
class ISMDP // ISMDP: Independent Set Matching based Detailed Placement
{
public:
    ISMDP()
    {
        Init();
    }
    ISMDP(PlaceDB *_db)
    {
        Init();
        placeDB = _db;
    }
    void Init()
    {
        placeDB = NULL;
        ISMRows.clear();
        maxModuleCount = 0;
        doubleWindow = false;
        independentCells = false;
    }
    void initialization();

    //! core function, see grid_run() in ntuplace3
    //! 用一个正方形窗口（边长 = windowSize 个行高）扫过整个芯片，
    //! 相邻窗口重叠 windowOverlap 个行高，避免窗口边界处的单元永远得不到优化
    void ISMSweep(int, int); // int windowSize, int windowOverlap

    PlaceDB *placeDB;
    vector<ISMRow> ISMRows;
    int maxWindow;      // see ntuplace class de_Detail.MAXWINDOW
    int maxModuleCount; // max number of modules in ISM window(s) (not one ISMRun()), see ntuplace class de_Detail.MAXMODULE
    //! 二者关系：maxWindow 控制一次匹配中 slot 的总数上界（矩阵规模、LAP 复杂度），
    //! maxModuleCount 控制从窗口里最多挑多少个单元参与匹配。
    //! 之所以分开：slot 数 = 被选中单元腾出的位置 + 额外加入的空白 slot，
    //! 而空白 slot 的多少取决于窗口内的空闲面积，需要单独限流。
    bool doubleWindow;
    bool independentCells;

private:
    //! 设置 maxWindow / maxModuleCount，并复位两个增强开关（默认都关闭）
    void initializeParams();

    //! 由 placedb 的 SiteRow::intervals 构建每一行的空闲区间，再填入已有单元
    void initializeISMRows();

    //! 对单个窗口执行一次完整的 ISM：选单元 -> 造 slot -> 填代价矩阵 -> LAP -> 写回
    void ISMRun(CRect);
};

// ============================================================================
// ISMRow：单行（placement row）的运行时视图，ISM 与 GlobalSwap 共用
// ============================================================================
// 为什么不直接用 placedb 里的 intervals？
//   因为匹配过程中单元会被反复「摘下来再放回去」，需要一个可以 O(log n)
//   增删、自动按 x 排序、并且能快速做区间合并 / 切分的容器。
//   std::map<double, ...> 以 x 坐标为键，天然有序且迭代器稳定，正好满足需求。
// 两个 map 共同描述同一行的占用情况，二者互补：
//   rowSpaces  —— 空闲区间：key = 区间起点 x，value = 区间长度
//   rowModules —— 已放单元：key = 单元左下角 x，value = 单元指针
// 不变式：两者覆盖的 x 区间互不重叠，合起来（再加 macro 等障碍）铺满整行。
class ISMRow // used for ISM and Global Swap
{
public:
    ISMRow()
    {
        Init();
    }
    //! _length 通常传 coreRegion.getWidth()，即整行长度；
    //! 行内真正可放的位置由 rowSpaces（来自 SiteRow::intervals）刻画，
    //! length 本身只作为行的几何范围记录，不参与合法性判断
    ISMRow(double _x = 0, double _y = 0, double _length = 0)
    {
        Init();
        ll.x = _x;
        ll.y = _y;
        length = _length;
    }
    void Init()
    {
        rowSpaces.clear();
        rowModules.clear();
    }

    //! 把一个单元放进本行：先找到能容下它的空闲区间，再把该区间「切掉」它占的部分。
    //! 返回 false 表示放不下（没有足够的连续空间 / 位置非法），调用方据此判失败。
    //! 这是**保证合法性的关键原语**：放不下就拒绝，绝不允许重叠。
    bool insertModule(double, Module *); // see ntuplace, detail.cpp insert_module()

    //! 摘走一个单元：把它占的 x 区间还回 rowSpaces（并做相邻区间合并），
    //! 同时从 rowModules 中删除。ISM 中被选中的单元会先被摘下，等匹配完再重新插入。
    void removeModule(Module *);

    //! 把 [x, x+width) 这段**原本空闲**的区间标记为「已被占用」（即从 rowSpaces 中扣掉），
    //! 用于给窄单元补齐到一个 maxWIDTH 宽的 slot。名字叫 removeTail 是因为
    //! 这段空间通常紧跟在单元右侧（tail）。与 insertSpace 互为逆操作。
    bool removeTail(double, double);  // double: x coordiante, double: width, see ntuplace, detail.cpp, remove_empty()

    //! 归还一段空闲区间 [x, x+width)，并与左右相邻的空闲区间做合并，
    //! 保证 rowSpaces 中不会出现两个首尾相接却分成两条的区间（否则会被误判成零散碎块）。
    bool insertSpace(double, double); // double: x coordiante, double: width, see ntuplace, detail.cpp add_empty()

    void showSpace();
    void showModule();

    POS_2D ll; // ll: lower left
    double length;
    map<double, double> rowSpaces;    // key: x-coordinate, value: length of the space of the , space is actually similar to interval, it is just inconvenient to operate directly the the intervals stored in placedb
    map<double, Module *> rowModules; // key: x-coordinate, value: module
    //! 注意：rowSpaces / rowModules 的下标都按「行内的绝对 x 坐标」，
    //! 而 ISMRow 之间用 ll.y 区分；模块跨多行时（macro）会在**每一行**各插一份，
    //! 因此 removeModule 只删一行里的那份记录，跨行单元需要调用方自己逐行处理。
    //? consider: use double or int as the key of the above 2 maps? I can use int for now, it should be easy to switch to double if necessary.
    //? consider: use double or int as the key of the above 2 maps? I can use int for now, it should be easy to switch to double if necessary.
    // use int as key is ok because the site step is usually integer, and x-coordinate is usually integer multiples of site step
};

// ============================================================================
// lap2：LAP（Linear Assignment Problem，线性指派问题）求解器
// ============================================================================
// 用途：ISM 里「哪个单元去哪个 slot」就是一个 LAP ——
//   行 = 单元（含 NULL 表示的空白），列 = slot，cost[i][j] = 单元 i 放到 slot j 的 HPWL，
//   求一个双射使总代价最小。因为所有 slot 等宽，任意双射都合法，
//   所以 LAP 的最优解就是窗口内这批单元的最优重排。
// 算法：Jonker-Volgenant 短增广路径（SAP）算法，O(n^3)。
//   先用列归约 + 行归约（reduction transfer）得到一个初始对偶可行解，
//   再用类 Dijkstra 的最短增广路径为尚未匹配的行逐个找增广路。
//   之所以手写成 int 版本：代价矩阵在 ISM 里被量化为整数（见 put()），
//   整数运算比 double 更快也避免浮点比较的抖动。
class lap2 // lap: linear assignment problem
{
public:
    lap2(int _deg)
    {
        degree = _deg;
        cost.resize(degree);
        for (unsigned int i = 0; i < cost.size(); ++i)
        {
            cost[i].resize(degree, 0);
        }
        assignment.resize(degree, 0);
        INF = INT_MAX;
        verbose = false;
    }
    ~lap2(void) {};

    //! 写入代价矩阵元素。注意这里把 double 线长**截断**成 int：
    //! 一是 LAP 用整数运算更快，二是线长本来就有小数噪声，截断误差可接受。
    //!//! 疑似问题：static_cast<int> 是向零截断而非四舍五入；若线长为负（不会）
    //!//! 或需要更高精度时会引入系统性偏差。另外 cost 为 int，若 HPWL 很大可能溢出。
    void put(const int &i, const int &j, double wl)
    {
        cost[i][j] = static_cast<int>(wl);
    }

    int INF;
    int degree;               // int m_deg;
    vector<vector<int>> cost; // vector<vector<int>> m_cost;
    vector<int> assignment;   // vector<int> m_assignment; assignment[i] = 第 i 行（单元）被指派的列（slot）

    //! 求解主函数，返回最优总代价；结果通过 getResult() 取回 assignment
    int lapSolve();

    void getResult(vector<int> &v)
    {
        v = this->assignment;
    }
    bool verbose;
};

// ============================================================================
// LocalReorderingDP：局部重排（= branch and bound cell swap）的驱动器
// ============================================================================
// 思路：把每一行的每个可放置区间（LRSegment）上的单元按 x 排序后，
//   用一个大小为 windowSize 的滑动窗口（相邻窗口重叠 overlapSize）取一小段连续单元，
//   在窗口内**穷举它们的排列**，找出 HPWL 最小的那个顺序。
// 为什么合法：窗口内单元的总宽度不变、占据的 x 跨度也不变（最左端对齐 currentX
//   依次紧凑排布），所以任何一种排列都只是「谁靠左谁靠右」的差别，
//   既不会互相重叠，也不会越出原区间，更不会破坏 site 对齐（宽度本身就是 site 的整数倍）。
// 为什么可行：窗口很小（默认 windowSize=3），排列数 n! 尚可枚举，
//   再配合分支限界（用「已完整放置的 net 的 HPWL 之和」作为单调下界）大幅剪枝。
class LocalReorderingDP // local reorder == bb cell swap
{
public:
    LocalReorderingDP()
    {
        Init();
    }
    LocalReorderingDP(PlaceDB *_db)
    {
        Init();
        placeDB = _db;
    }
    void Init()
    {
        placeDB = NULL;
        lrSegments.clear();
    }
    PlaceDB *placeDB;
    vector<LRSegment> lrSegments;

    //! 构建所有 LRSegment（来自 SiteRow::intervals）并把单元挂到对应 segment 上、按 x 排序
    void initialization();

    //! core function
    //! windowSize       ：一个窗口里放几个单元（默认 3，排列数 3! = 6，可枚举）
    //! overlapSize      ：相邻窗口重叠几个单元（默认 2，即窗口每次只右移 1 个单元），
    //!                    重叠是为了让跨窗口边界的单元也有机会被一起重排
    //! iterationNumber  ：整体扫描的轮数（默认 1）
    void solve(int, int, int); // int windowSize, int overlapSize, int iterationNumber

private:
    //! 由 placedb 的 SiteRow::intervals 生成 LRSegment，并把非 macro 单元二分查找后挂上去
    void initializeSegments();

    //! 对 [startModule, endModule) 这段相邻单元求最优排列，并把结果写回 placedb
    //! 返回 true 表示相比原顺序确有改善
    bool solveForBestOrder(vector<Module *>::iterator, vector<Module *>::iterator); // todo: implement me
};

// ============================================================================
// LRSegment：一行中的一个可放置区间（segment），概念上等价于 abacus 的一行
// ============================================================================
// 由 placedb->dbSiteRows 的 intervals 生成（intervals 已经剔除 macro / terminal 占用的 site）。
// 一个 segment 是「连续、中间没有障碍」的一段，因此其内部可以自由重排单元；
// 而 segment 之间的空隙不能放东西，所以重排绝不会跨 segment 进行。
class LRSegment // The placeable interval in each row is called a segment, pretty much like an abacus row. Initialized with placedb->dbSiteRows.intervals
{               // LR:local reordering
public:
    LRSegment()
    {
        Init();
    }
    //! 探针构造函数：用一个单元的 (y, x, x+width) 造一个「零宽度区间」式的 segment，
    //! 配合 compareLRSegment 做 lower_bound，从而定位该单元属于哪个 segment
    LRSegment(double _bottom, double _start, double _end)
    {
        Init();
        bottom = _bottom;
        start = _start;
        end = _end;
    }
    void Init()
    {
        bottom = -1;
        start = -1;
        end = -1;
        segModules.clear();
    }

    double bottom;               // The bottom y coordinate of this row of sites
    double start, end;           // The left and right coordinates of this row of sites
    vector<Module *> segModules; // Record the  modules on this segment
    //! segModules 在 initializeSegments() 末尾会按 x 升序排序，
    //! 这是 solve() 里「窗口 = 连续下标区间」这一假设成立的前提。
    void addModule(Module *module)
    {
        segModules.push_back(module);
    }
};

// ============================================================================
// compareLRSegment：LRSegment 的严格弱序比较器（供 lower_bound 二分查找使用）
// ============================================================================
// 排序键是 (bottom, start, end)：先按行（y）分，再按区间左端分，最后按右端分。
// 用途：给一个单元的 (y, x, x+width) 构造一个「探针 segment」，
//       用 lower_bound 在有序的 lrSegments 里定位它落在哪个 segment 中。
// //! 疑似问题：该比较器的语义是「区间包含/相邻判定」而非纯粹的字典序。
// //!   当 s.bottom == compSeg.bottom 且 s.start > compSeg.start 时直接返回 false，
// //!   而当 s.start <= compSeg.start 时又返回 s.end < compSeg.end，
// //!   这并不构成严格的字典序（缺少对 start 的直接比较），
// //!   可能导致 lower_bound 在存在重叠/相邻区间时定位到错误的 segment。
// //!   原作者在调用处也标注了 "//!!!! potential bug"，详见 detailed.cpp。
class compareLRSegment
{
public:
    bool operator()(const LRSegment &s, const LRSegment &compSeg)
    {
        if (s.bottom == compSeg.bottom)
        {
            if (s.start <= compSeg.start)
            {
                return s.end < compSeg.end;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return s.bottom < compSeg.bottom;
        }
    }
};

// ============================================================================
// LRSolution：局部重排窗口内的一个解（单元的一种排列）
// ============================================================================
// 这里的「解」指的是 DFS 搜索树上的一个节点，既可以是**完整解**（所有单元都排好，
// 位置确定，可以算出真实 HPWL），也可以是**部分解**（只有前缀若干单元排好）。
// 之所以要显式区分，是因为代价要**增量式**计算：
//   只有当某条 net 上的**窗口内所有单元**都排好之后，这条 net 的 bounding box
//   才最终确定，此时才把它的 HPWL 计入 solutionCost。
//   这样得到的 solutionCost 沿搜索路径**单调不减**，是一个合法的下界（bound），
//   可以用来做分支限界剪枝（见 LRSolver::depthFirstSolve）。
// 注意：为了做增量计算，单元会被**真实地移动**（setModuleLocation_2D），
//   因此搜索过程中 placedb 里的坐标是在不断被改写的，最终由 solveForBestOrder
//   统一按最佳解写回。
class LRSolution // one solution (ordering) of cells in a local reorder window
{
public:
    LRSolution(PlaceDB *_db)
    {
        Init();
        placedb = _db;
    }
    void Init()
    {
        placedb = NULL;
        solutionCost = 0;
        insertedCells.clear();
        uninsertedCells.clear();
        whiteSpaceWidth = 0;
        currentX = -1;
        xLocations.clear();
        netModuleCount.clear();
    }
    // LRSolution( const LRSolution& sol );

    //! 完整解 = 窗口内所有单元（含空白）都已排定位置
    bool isComplete() { return (uninsertedCells.size() == 0); } // insertedCells + uninsertedCells == all modules in the window
    double getCost()
    {
        return solutionCost;
    }
    //! 深拷贝一份解。LRSolver 保存最优解时必须 clone，
    //! 因为搜索过程中的临时解会被立刻 delete 掉（见 depthFirstSolve）。
    LRSolution *clone();
    void print();

    //! 统计窗口内每条 net 上**尚未排定**的单元个数。
    //! 每排定一个单元就把它所在 net 的计数减一，减到 0 说明这条 net 已完整，
    //! 此时才结算它的 HPWL —— 这是增量式代价计算的基础。
    void initializeNetModuleCount();

    void recalculateCost(); // Recalculate the total wirelength in insertedCells if this solution is complete

    PlaceDB *placedb;
    double solutionCost;            // For incrementally compute the solution cost, see m_bound in ntuplace3
    list<Module *> insertedCells;   // inserted modules (locations are determined in this solution)
    list<Module *> uninsertedCells; // to be inserted modules (locations not determined in this solution, NULL is whitespace)
    //! 窗口中所有空白被**合并成一个** NULL 单元参与排列
    //! （因为空白内部再怎么切分对线长无影响，合并可以显著减小搜索空间）。
    //! 代价：空白只能整体出现在某一处，无法拆成两段分别插在中间。
    double whiteSpaceWidth;         // The width of the whitespace
    //! 下一个待排单元应当放置的 x 坐标（左端），随排列推进不断右移
    double currentX;                // left most x coordinate of the whole segment of cells

    list<double> xLocations; // Record the x coordinate of each module in insertedCells

    //! key: net id, value: 该 net 上**窗口内**还没排定的单元数。
    //! 用 int 而非 Net* 作 key，是为了避开 std::map 用指针比较时的不确定性
    //! （见 initializeNetModuleCount() 里被注释掉的「find() fails when using pointer as key」）。
    map<int, int> netModuleCount; // number of modules in each related net of the modules in the window, key: net id, do not use Net* as key
};

// ============================================================================
// LRSolutionIterator：由一个（部分）解枚举出它的所有后继解
// ============================================================================
// 后继 = 「把某一个尚未排定的单元放到 currentX 处」。
// 因此一个含 k 个未排单元的部分解有 k 个后继，对应搜索树上 k 条分支；
// 对所有后继递归下去，就枚举出了窗口内单元的**全部排列**。
class LRSolutionIterator
{
public:
    LRSolutionIterator(const LRSolution *sol)
        : pointedSolution(sol)
    {
        reset();
    }
    void reset()
    {
        moduleIterator = pointedSolution->uninsertedCells.begin();
    }
    bool isDone()
    {
        return moduleIterator == pointedSolution->uninsertedCells.end();
    }
    //! 生成「把当前 moduleIterator 指向的单元放到 currentX」这一后继解，
    //! 同时增量更新 solutionCost / netModuleCount / xLocations
    LRSolution *createSuccessorSolution();

    list<Module *>::const_iterator moduleIterator; // The iterator points to current uninserted cell in m_sol
    const LRSolution *pointedSolution;
};

// ============================================================================
// LRSolver：局部重排的 DFS + 分支限界（branch and bound）求解器
// ============================================================================
// 搜索树：根节点 = 所有单元都未排定（initialSolution），
//   第 d 层的节点 = 已排定 d 个单元，边 = 选择下一个放哪个单元。
//   叶子 = 一个完整排列，其 solutionCost 就是该排列的真实 HPWL（窗口相关部分）。
// 剪枝：节点的 solutionCost 是单调不减的下界（只统计了已完整 net 的 HPWL），
//   所以一旦某节点的 cost 已经 >= bestCost，它的整棵子树都不可能更优，可直接剪掉。
// 初始上界：调用方先 setBestSolution(当前实际布局)，
//   这样搜索天然带一个「不比现状差」的保证（找不到更好的就保持原样）。
class LRSolver
{
public:
    LRSolver() : bestSolution(NULL), bestCost(DOUBLE_MAX) {}
    ~LRSolver() { delete bestSolution; }

    LRSolution *bestSolution;
    double bestCost;

    //! 入口：从 initialSolution 开始深搜，返回最优完整解（按值返回）
    LRSolution solve(LRSolution *);

    //! 设置初始上界（通常是当前的实际布局），只接受完整解
    void setBestSolution(LRSolution *);

    //! 若给定完整解更优则更新 bestSolution / bestCost
    void updateBestSolution(LRSolution *);

    //! 递归主体：完整解则更新最优，否则枚举后继并带剪枝地递归
    void depthFirstSolve(LRSolution *);
};

// ============================================================================
// GlobalSwapDP：全局交换（Global Swap）
// ============================================================================
// 论文：An efficient and effective detailed placement algorithm（Sun et al.）
// 与前两个算子不同，GlobalSwap 的搜索范围不是固定的几何窗口，
//   而是**每个单元各自的 optimal region** —— 由该单元相连的各条 net 的
//   bounding box 交集决定（placedb->getOptimialRegion），
//   直观上就是「这个单元放到哪里线长可能更短」的那一小片区域。
//   （理论依据：把单元移出其 net bounding box 之外只会让线长变大。）
// 候选位置的两种形式，都以「不破坏合法性」为前提：
//   1) 空白 slot：从 optimal region 内的空闲区间里切出若干等宽（= 本单元宽度）的位置；
//   2) 等宽单元：与 region 内**宽度相同**的单元对调位置（同宽 -> 对调后仍不重叠）。
// 效率取舍：论文明确指出，实际实现里**不找收益最大的候选，而找第一个有收益的候选**，
//   以大幅节省运行时间（见 solve() 中的 "looking for the first one ..." 注释）。
class GlobalSwapDP
{
public:
    GlobalSwapDP()
    {
        Init();
    }
    GlobalSwapDP(PlaceDB *_db)
    {
        Init();
        placeDB = _db;
    }
    void Init()
    {
        placeDB = NULL;
    }
    //! 主循环：遍历所有非 macro 单元，为它们各自在 optimal region 内找更好的位置
    void solve();
    PlaceDB *placeDB;
    //! 复用 ISMRow 作为行的运行时结构（空闲区间 rowSpaces + 已放单元 rowModules），
    //! 与 ISMDP::ISMRows 完全同构，只是由 GlobalSwap 自己维护
    vector<ISMRow> GSRows;
    int maxSlotCount; // max number of candidate slots chosen from the search region
    //! 注意 maxSlotCount 同时限制「空白 slot 数」和「空白 + 可交换单元数」，
    //! 是运行时间与优化质量之间的调节旋钮（默认 30）。

    //! 构建参数与 GSRows。//! 疑似问题：initialization() 只在 runGlobalSwap() 里调用一次，
    //! 而 solve() 会被调用 2 次；第二次复用的是已经被前一轮修改过的 GSRows。
    //! 若这是有意为之（增量维护）则没问题，但一旦中途有单元移动未被同步，
    //! GSRows 就会与 placedb 的真实坐标失配，且没有任何校验机制。
    void initialization();

private:
    void initializeParams();
    void initializeGSRows();
};
#endif