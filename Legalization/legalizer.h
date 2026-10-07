// ============================================================================
// 模块总览：Legalization（合法化）
// ============================================================================
// 【这个头文件/模块是干什么的】
//   上游 EPlace（全局布局）输出的结果只是「连续坐标上的密度均衡解」：
//   std cell 之间大量重叠、y 坐标不在任何一条 site row 的底边上、
//   x 坐标也没有落在 site 网格上。Legalization 负责把这份结果变成
//   物理上可制造的布局：消除重叠、x 对齐 site 网格、y 对齐 site row，
//   并且尽量让每个 cell 的移动距离最小（总位移越小，线长退化越少）。
//
// 【算法思路：为什么是 Abacus 而不是 Tetris】
//   模块主力是 AbacusLegalizer，严格按 Abacus 论文
//   (Spindler & Johannes, "Fast and Accurate Standard-Cell Based Legalization
//   with Consideration of Whitespace", DAC 2007) 实现。整体分三层：
//     ① 行分配（row assignment）
//        逐个 cell（按全局 x 递增处理）挑一条 subrow：以「最近行」为中心
//        向上、向下双向搜索，用「只算 y 位移」的下界做剪枝，
//        对少量候选行调用 placeRow(ABACUS_TRIAL) 在 subrow 的**拷贝**上试算真实代价，
//        取代价最小者再 placeRow(ABACUS_FINAL) 真正写入。
//     ② 行内合法化（in-row legalization）
//        每条 subrow 上维护一串 cluster。cell 放进来后若与前驱 cluster 重叠，
//        就把两者合并成一个刚体 cluster，再重新求该 cluster 的最优 x。
//        这是 Abacus 优于 Tetris 的关键：Tetris 只把 cell 往最近的空位里塞，
//        而 Abacus 会为了整体位移最小而「推着前面的一堆 cell 一起挪」。
//     ③ 坐标回写
//        placeRow 只决定 cluster 的位置，最后统一把 cluster 展开成各 cell 的 x
//        并写回 PlaceDB。
//
//   另一条支线是 SAMacroLegalizer（macro 合法化，对应 ePlace-MS 的 mLG）：
//   macro 不能像 std cell 那样塞进 site row，必须先单独合法化，
//   所以这里用「模拟退火 + 三类代价（线长 / 被压住的 cell 面积 / macro 间重叠）」
//   把 macro 推到互不重叠的位置；macro 定下来后，
//   PlaceDB::removeBlockedSite() 才能算出每条 site row 真正可用的 intervals。
//   因此调用顺序必须是：macro 合法化 → removeBlockedSite → std cell 的 Abacus 合法化。
//
// 【核心数据结构】
//   Obstacle                    —— 占住 site 的障碍物矩形（macro/terminal），已基本退场
//   Segtree / RectangleAreaSolution —— 扫描线 + 线段树，用来求「所有 macro 的并面积」，
//                                     进而得到 macro 之间的重叠量
//   SAMacroLegalizationBin_2D   —— mLG 的密度 bin，预先冻结 cell 密度以加速代价评估
//   AbacusCellCluster           —— Abacus 的 cluster，四个量 x/e/w/q 直接对应论文
//   AbacusRow                   —— subrow = 一条 site row 上一段连续可用区间 + 其 cluster 序列
//
// 【主要入口函数】
//   AbacusLegalizer::legalization()     —— std cell 合法化入口
//   SAMacroLegalizer::legalization()    —— macro 合法化入口
//   Legalization/ePlaceAbacus.cpp:main  —— 串联「读全局布局 .pl → 合法化 → 输出」的可执行入口
// ============================================================================

#ifndef LEGALIZER_H
#define LEGALIZER_H
#include "global.h"
#include "objects.h"
#include "placedb.h"
#include "plot.h"
#define ABACUS_TRIAL true
#define ABACUS_FINAL false
// ----------------------------------------------------------------------------
// Obstacle：占住 site 的障碍物，本质就是 CRect 的别名（矩形 ll/ur）
// 只服务于「把一条完整 site row 按障碍物切成若干可用 subrow」这一几何步骤，
// 并不参与后续的 Abacus 行内聚类。
//! 疑似问题：该类型是空壳（无任何新增成员），且与 placedb.cpp::removeBlockedSite()
//! 中就地使用的 CRect 障碍物重复。主线流程改为由 PlaceDB 直接生成
//! SiteRow::intervals 之后，本类型只剩已废弃的 initializeObstacles /
//! initializeSubrowsOld 还在引用，属于遗留死代码。
// ----------------------------------------------------------------------------
class Obstacle : public CRect
{
};

// ----------------------------------------------------------------------------
// Segtree：扫描线求「矩形并面积」时用到的线段树结点（见 RectangleAreaSolution）
//   cover      —— 当前结点代表的区间被矩形覆盖的次数；>0 表示整段被覆盖
//   length     —— 当前结点区间内**实际被覆盖**的总长度（由 pushup 自底向上维护）
//   max_length —— 当前结点区间的几何总长度，与覆盖无关，建树时一次性确定
// 为什么需要 max_length：当 cover>0 时该区间必然整段被覆盖，
// 直接取 max_length 即可，不必再递归累加子结点，从而把区间更新压到 O(log n)。
// ----------------------------------------------------------------------------
struct Segtree
{
    int cover;
    int length;
    int max_length;
};

// ----------------------------------------------------------------------------
// SAMacroLegalizationBin_2D：macro 合法化（mLG）使用的密度 bin
// 与 EPlace 全局布局的 bin 同构，但只服务于一件事：
//   快速估算「某个 macro 压住了多少 std cell 面积」（getCellAreaCoveredByMacro）。
// 手法是把每个 bin 的 nonMacroDensity = (cellArea + terminalArea + baseArea) / area
// 在初始化时算好并**冻结**。因为 mLG 期间 std cell 是固定不动的，
// 所以 macro 移动时只需拿「macro 与 bin 的重叠面积 × 该 bin 的密度」即可，
// 不必每次 SA 扰动都重算全芯片密度——这是 mLG 能在合理时间内跑完的关键近似。
// ----------------------------------------------------------------------------
class SAMacroLegalizationBin_2D
{
public:
    SAMacroLegalizationBin_2D()
    {
        init();
    }
    POS_2D center; // bin 中心，供后续按 bin 计算梯度/代价时使用
    POS_2D ll;     // bin 左下角（绝对坐标，已叠加 coreRegion.ll）
    POS_2D ur;     // bin 右上角
    float width;
    float height;
    float area; // bin 的几何面积；注意它会随 width/height 变化，由 getArea() 顺带刷新

    float cellArea;     // overlapped area between this bin and std cells
    float macroArea;    // overlapped area between this bin and macros
    float terminalArea; // overlapped area between this bin and terminals, should always be 0!!!!
    float baseArea;
    // todo: add virtual area?
    float nonMacroDensity; // = (cell area + terminal area + base area) / bin area

    //! 疑似问题：原注释称 terminalArea 恒为 0，但 initializeBins() 中确实会为
    //! 落在 coreRegion 内的 terminal 累加重叠面积，二者矛盾；
    //! 此外 macroArea 虽在每个 bin 上被累加，却从未参与任何代价计算
    //! （代价里的「macro 重叠」走的是 getMacroOverlapArea 的精确两两求交），属死变量。
    void init()
    {
        center.SetZero();
        ll.SetZero();
        ur.SetZero();
        cellArea = 0;
        macroArea = 0;
        terminalArea = 0;
        baseArea = 0;
        area = 0;
        width = 0;
        height = 0;
    }
    float getWidth() { return width; }
    float getHeight() { return height; }
    float getArea()
    {
        area = width * height;
        return width * height;
    }
};

// ----------------------------------------------------------------------------
// AbacusCellCluster：Abacus 算法中的「cluster（簇）」
//
// 为什么要 cluster：当若干 cell 在同一行内互相挤压时，它们必须紧贴成一段。
// 把这段 cell 捆成一个刚体、整体在行内滑动，就把「多个 cell 的最优摆放」
// 化简成「一个刚体的最优位置」——这正是 Abacus 比 Tetris 精度更高的根源：
// Tetris 只会把 cell 塞进最近的空位，Abacus 则会为了总位移最小而
// 推着前面的一整串 cell 一起移动。
//
// 四个量严格对应 Abacus 论文：
//   x —— cluster 左下角的 x 坐标（放置结果，必须对齐到 site step）
//   e —— cluster 内 cell 的个数
//   w —— cluster 内所有 cell 宽度之和，即 cluster 的总宽度
//   q —— Σ(cell_i 的理想 x − cluster 中排在它前面的 cell 宽度之和)。
//        直观理解：把每个 cell 的期望位置折算到「以 cluster 左端为原点」的坐标系后
//        再求和；于是 cluster 的理想左端位置就是 q/e，
//        它正好让所有 cell 的位移平方和（以及绝对位移和）最小。
// ----------------------------------------------------------------------------
class AbacusCellCluster
{
public:
    // cluster of std cells, used in abacus
    AbacusCellCluster()
    {
        Init();
    }
    void Init()
    {
        cells.clear();
        index = -1;
        x = 0.0;
        e = 0.0;
        w = 0.0;
        q = 0.0;
    }
    vector<Module *> cells; // cells should be ordered by x coordinate(non decreasing)
    int index;              // its index in the vector clusters of class AbacusRow
    //! 疑似问题：index 是 cluster 在 AbacusRow::clusters 中的下标缓存，
    //! 但 collapse() 合并后会执行 clusters.erase()，其后所有 cluster 的 index 均未回填，
    //! 于是 index 会整体偏大。placeRow() 里 addCell(clusters.back().index, cell)
    //! 正是拿这个陈旧 index 去索引，理论上可能越界（见 legalizer.cpp 中的同款标注）。
    // x,e,w,q of a cluster, see the abacus paper
    float x;
    float e;
    float w;
    float q;
};

// ----------------------------------------------------------------------------
// AbacusRow：一条 subrow
// subrow = 一条 site row 被 macro / terminal 切断后剩下的一段连续可用区间。
// 继承 SiteRow 得到 bottom / height / step / start / end 这些几何属性，
// 再额外维护：
//   clusters —— 行内已放置的 cluster 序列，按 x 递增排列（Abacus 合并只发生在相邻项之间）
//   width    —— 已放入 cell 的宽度之和
// width 只是一个「已占用宽度」的上界估计（不含 cluster 之间被推挤出来的空隙），
// 用途是快速判断某行是否还塞得下新 cell，从而跳过昂贵的 placeRow 试算。
// ----------------------------------------------------------------------------
class AbacusRow : public SiteRow
{
public:
    AbacusRow()
    {
        clusters.clear();
        width = 0;
        // lastClusterIndex=-1;
    }
    vector<AbacusCellCluster> clusters;
    double width;
    void addCell(int, Module *);      // 把一个 cell 追加进指定 cluster，并增量更新 e/q/w
    void addCluster(int, int);        // 把后一个 cluster 合并进前一个（Abacus 论文的 merge 步骤）
    void collapse(int);               // 重算 cluster 的最优 x，并在重叠时递归向前合并
};

// ----------------------------------------------------------------------------
// AbacusLegalizer：标准单元合法化器（本模块的主力）
//
// 算法脉络（严格遵循 Abacus 论文）：
//   ① 初始化 initialization()
//      - initializeCells()：收集所有非 macro 的 std cell，按全局 x 升序排序。
//        排序很重要：后续按此顺序插入，才能保证同一 cluster 内 cells 天然按 x 有序，
//        也保证 addCell 增量维护 q 的语义正确。
//      - initializeSubrows()：先调 PlaceDB::removeBlockedSite() 得到每条 site row
//        被 macro/terminal 切断后的 intervals，再把每个 interval 展开成一条
//        AbacusRow，最后按 bottom 升序排序（双向搜索 + 剪枝的前提）。
//   ② 行分配 legalization()
//      - 定位「最近行」closestRowIndex：cell 的 y 落在哪条 row 的高度范围内；
//      - 以最近行为中心，先向上、再向下搜索候选行：
//        用 costLowerBound = |cell.y − row.bottom| 作剪枝下界
//        （真实代价是含 x 位移的欧氏距离，必然 ≥ 该下界），
//        一旦下界已不小于当前最优 cost，更远的行只会更差，直接 break；
//      - 对候选行调 placeRow(..., ABACUS_TRIAL) 试算，试算在 subrow 的**拷贝**上进行。
//   ③ 行内合法化 placeRow()
//      - 若该 cell 可以独立成簇则新建 cluster，否则并入行内最后一个 cluster；
//      - collapse() 重算 cluster 的最优 x 并向下对齐到 site step、夹紧在 subrow 内，
//        若与前驱 cluster 重叠则合并成刚体并递归向前 collapse。
//   ④ 坐标回写：placeRow 只决定 cluster 的位置，最后统一把 cluster 展开成
//      各 cell 的实际 (x, row.bottom) 写回 PlaceDB。
//
// 为什么「先试算再落地」：cluster 合并是破坏性操作，把 cell 放进某行后再想换行的
// 代价极高；而对少量候选行各做一次拷贝试算的代价，远小于行分配错误带来的位移损失。
// ----------------------------------------------------------------------------
class AbacusLegalizer
{
    // strictly follow the abacus paper, read the paper for more details
public:
    AbacusLegalizer()
    {
        Init();
    }
    AbacusLegalizer(PlaceDB *_db)
    {
        Init();
        placeDB = _db;
    }
    void Init()
    {
        placeDB = NULL;
    }
    PlaceDB *placeDB;
    vector<AbacusRow> subrows;  // generated from dbSiteRows and Terminals, this is actually used during legalization
    vector<Obstacle> obstacles; // include macro and terminals, macros should be legalized first! this vector is used to obtain subrows, not necessary if we switch to initialize intervals for subrows in a member function of placedb 
    vector<Module *> dbCells;   // all std cells, sorted by x coordinate(non decreasing order)
    void legalization();

private:
    // set these functions private because they are called in other member functions only
    // follow the abacus paper
    void initialization();
    void initializeCells();
    void initializeObstacles(); // obstacles include macros and terminals, macros should be legalized first! 2024.08.26: this function is abandoned, and is integrated into removeBlockSites() in placedb.cpp
    void initializeSubrows();

    void initializeSubrowsOld();// abandoned

    double placeRow(Module *, int, bool); // trial=true 为试算（不写状态），false 为真正落地
};

// ----------------------------------------------------------------------------
// SAMacroLegalizer：macro 合法化器（对应 ePlace-MS 论文里的 mLG 阶段）
//
// 为什么 macro 要单独合法化：macro 尺寸大、不能像 std cell 那样塞进某条 site row，
// 而且它们会**占据并切断** site row。只有先把 macro 推到互不重叠的位置，
// PlaceDB::removeBlockedSite() 才能算出每条 site row 真正可用的 intervals，
// std cell 的 Abacus 合法化才有落脚点。所以 mLG 必须跑在 Abacus 之前。
//
// 算法思路（模拟退火，SA）：
//   代价函数（ePlace-MS 公式 30 的三项）：
//       cost = HPWL + miuD * (被 macro 压住的 std cell 面积) + miuO * (macro 之间的重叠面积)
//     第一项保证线长不退化；第二项避免 macro 压在 cell 密集区；
//     第三项直接驱动「消除 macro 重叠」这个唯一硬约束。
//   miuD / miuO 由初始量纲自动标定（见 initializeSAparams），
//   使三项在初始时刻处于同一数量级，避免某一项被淹没。
//   外层 j 循环（mLG iteration）：每轮重设 SA 温度与扰动半径 r 并退火，
//   同时按 beta 放大 miuO，让「消除重叠」的权重逐轮变强；
//   内层 k 循环（SA iteration）：温度与半径逐步衰减，每步随机挑一个 macro 做扰动。
//   终止条件：overlapFree，即「macro 总面积 − macro 并面积」降到 0。
//     注意这里用「并面积」而不是两两求交：并面积由扫描线+线段树 O(n log n) 求得，
//     比 O(n²) 的两两求交快得多，而且两者为 0 是等价的。
//
// 关键假设（原作者注释）：terminal 中不含 macro，且 MMS benchmark 的 terminal 面积为 0。
// ----------------------------------------------------------------------------
class SAMacroLegalizer
{
    // todo: implement a sa-based macro legalizer according to ePlace-MS
    // ! a crucial assumption: no macro in terminals!!!! terminals in the MMS benchmark should all have 0 area
public:
    SAMacroLegalizer()
    {
        Init();
    }
    SAMacroLegalizer(PlaceDB *_db)
    {
        Init();
        placeDB = _db;
    }
    void Init()
    {
        placeDB = NULL;
        dbMacros.clear();
        totalMacroArea = 0;
        totalCellArea = 0;
        bins.clear();
        binDimension.SetZero();
        binStep.SetZero();
        targetDensity = 1.0;

        totalHPWL = 0;
        totalCellAreaCovered = 0;
        totalMacroOverlap = 0;
        overlapFree = false;

        miuD = 0;
        miuO = 0;
        jLimit = 0;
        kLimit = 0;
        SAtemperature = 0;
        SAtemperatureCoef = 0;
        r.SetZero();
        u.SetZero();
        sa_r_stp.SetZero();
        beta = 0;
    }
    PlaceDB *placeDB;
    vector<Module *> dbMacros; // all std cells, sorted by x coordinate(non decreasing order)
    //! 疑似问题：注释与实际不符。dbMacros 里存的是 **macro**（非 std cell），
    //! 而且 initializeMacros() 之后并没有做任何排序，所以「sorted by x」也不成立。

    double totalMacroArea;// precision problem!!!! use int for total macro area causes bug
    double totalCellArea;

    float targetDensity; // equals to eplacer target density if mGP was ran

    //! bins 的第一维是 x、第二维是 y；元素用 new 分配且全程没有 delete，
    //! 属于可控的一次性泄漏（进程级工具），但若在同进程内重复初始化会累积泄漏。
    vector<vector<SAMacroLegalizationBin_2D *>> bins;
    VECTOR_2D_INT binDimension; // How many bins in X/Y direction
    VECTOR_2D binStep;          // length of a bin in X/Y direction

    // three things in the SA cost function!
    double totalHPWL;
    double totalCellAreaCovered;
    int totalMacroOverlap;
    //! 疑似问题：totalMacroOverlap 声明为 int，但 SAperturb() 里用
    //! 「+= 浮点差值」增量维护，每步都会截断取整，长迭代下误差会累积。
    //! 函数调用栈：SAMacroLegalizer::legalization() → SAMacroLegalization()
    //! → SAperturb() → getMacroCost()；由 ePlaceAbacus.cpp / eplace.cpp 按需触发。

    bool overlapFree; // key stop condition of macro legalization
    //! SA params
    //! 0. cost weights
    float miuD; // see ePlace-MS paper equation(30)
    float miuO; // see ePlace-MS paper equation(30)
    //! 1. iteration count upper limits
    int jLimit;
    int kLimit;
    //! 2. params related to SA computation
    float SAtemperature;
    float SAtemperatureCoef;
    VECTOR_2D r;                 // 当前扰动半径（x/y 方向的最大位移幅度），随退火逐步收缩
    VECTOR_2D u;                 // 见下
    VECTOR_2D sa_r_stp;          // name is same as in RePlAce
    //! 疑似问题：u 在 initializeSAparams/legalization 中被置为 (1.0, 1.0) 后
    //! 再未被修改，乘上它恒等于乘 1，是 RePlAce 遗产里的死变量。
    float beta;         // see ePlace-MS paper

    void legalization();

    void setTargetDensity(float _targetdensity)
    {
        targetDensity = _targetdensity;
    }

private:
    // set these functions private because they are called in other member functions only
    // follow the abacus paper
    void initialization();

    void initializeMacros();   // 收集 macro、累加面积，并把 macro 坐标离散化到整数/行高网格
    void initializeBins();     // 划分 bin 网格并预计算 terminal/base/cell 三类面积
    void initializeCost();     // 计算 HPWL / 被压 cell 面积 / macro 重叠三项初始代价
    void initializeSAparams(); // 由初始代价量纲反推 miuD / miuO，并设定迭代上限与 beta

    int getAreaCoveredByMacros(); // 用扫描线求所有 macro 的**并面积**（判断重叠是否为 0）
    double getMacrosOverlap();    // O(n^2) 两两求交的重叠面积，精度高但慢，仅调试用
    float getCellAreaCoveredByAllMacros();

    float getCellAreaCoveredByMacro(Module *); // 用冻结的 bin 密度近似估算被压住的 cell 面积
    int getMacroOverlapArea(Module *);         // 单个 macro 与其余 macro 的重叠面积（精确两两求交）
    float getMacroCost(Module *, float &, float &, float &);

    void SAMacroLegalization(); // 内层 SA 循环（k 循环）
    void SAperturb();           // 随机挑一个 macro 扰动一步，并按 Metropolis 准则接受/回滚
    bool acceptPerturb(float, float);
};

// ----------------------------------------------------------------------------
// RectangleAreaSolution：求「一组矩形的**并面积**」（矩形面积并，LeetCode 850 原题解法）
//
// 它在 legalization 里的用途（见 SAMacroLegalizer::getAreaCoveredByMacros）：
//   macro 之间的重叠量 = macro 总面积 − macro 并面积。
//   并面积为 0 的重叠量 ⇔ 所有 macro 互不重叠，这正是 mLG 的终止判据。
//   用扫描线在 O(n log n) 内求并面积，替代 O(n²) 的两两求交，
//   否则每轮 SA 都要重算，根本跑不动。
//
// 算法思路（扫描线 + 线段树）：
//   ① 把所有矩形的上下边界 y 收集起来排序去重（hbound），
//      相邻两个 y 之间形成一条「水平条带」，共 m-1 条；
//   ② 把每个矩形的左右边界做成 +1/-1 事件按 x 排序，从左向右扫描；
//   ③ 线段树维护「当前 x 处，被覆盖的 y 区间总长度」：
//      结点记录 cover（被覆盖次数）与 length（实际覆盖长度）；
//      cover>0 时整段取满 max_length，否则由子结点累加；
//   ④ 相邻两条扫描线之间的面积 = 当前覆盖总长度 × x 间距，累加即得并面积。
// 注意：这里的坐标必须是**整数**（macro 坐标已在 initializeMacros 中离散化），
// 因为线段树按离散 y 值建树，浮点会破坏去重与二分的正确性。
// ----------------------------------------------------------------------------
class RectangleAreaSolution
{
public:
    int rectangleArea(vector<vector<int>> &rectangles)
    {
        int n = rectangles.size();
        for (const auto &rect : rectangles)
        {
            // 下边界
            hbound.push_back(rect[1]);
            // 上边界
            hbound.push_back(rect[3]);
        }
        sort(hbound.begin(), hbound.end());
        hbound.erase(unique(hbound.begin(), hbound.end()), hbound.end());
        int m = hbound.size();
        // 线段树有 m-1 个叶子节点，对应着 m-1 个会被完整覆盖的线段，需要开辟 ~4m 大小的空间
        tree.resize(m * 4 + 1);
        init(1, 1, m - 1);

        vector<tuple<int, int, int>> sweep;
        for (int i = 0; i < n; ++i)
        {
            // 左边界
            sweep.emplace_back(rectangles[i][0], i, 1);
            // 右边界
            sweep.emplace_back(rectangles[i][2], i, -1);
        }
        sort(sweep.begin(), sweep.end());

        long long ans = 0;
        for (int i = 0; i < sweep.size(); ++i)
        {
            int j = i;
            // 把横坐标相同的所有事件归为一批：同一 x 上的覆盖变化必须一起生效，
            // 否则中间会出现「半更新」状态，导致条带长度算错
            while (j + 1 < sweep.size() && get<0>(sweep[i]) == get<0>(sweep[j + 1]))
            {
                ++j;
            }
            if (j + 1 == sweep.size())
            {
                // 最后一批事件之后已经没有下一个 x，无法再贡献面积，直接结束
                break;
            }
            // 一次性地处理掉一批横坐标相同的左右边界
            for (int k = i; k <= j; ++k)
            {
                auto &&[_, idx, diff] = sweep[k];
                // 使用二分查找得到完整覆盖的线段的编号范围
                int left = lower_bound(hbound.begin(), hbound.end(), rectangles[idx][1]) - hbound.begin() + 1;
                int right = lower_bound(hbound.begin(), hbound.end(), rectangles[idx][3]) - hbound.begin();
                update(1, 1, m - 1, left, right, diff);
            }
            ans += static_cast<long long>(tree[1].length) * (get<0>(sweep[j + 1]) - get<0>(sweep[j]));
            i = j;
        }
        //! 疑似问题：累加器 ans 是 long long，但函数返回类型是 int，
        //! 大规模设计里并面积可能溢出 int（>2^31）而被截断。
        //! 由于调用方只是拿它和 totalMacroArea 作差判 0，溢出会直接误判「无重叠」。
        return ans;
    }

    // 建树：确定每个结点区间的几何总长度 max_length（注意叶子 l 的长度取自 hbound[l]-hbound[l-1]，
    // 所以叶子编号从 1 开始，正好对应「第 l-1 条与第 l 条水平线之间的条带」）
    void init(int idx, int l, int r)
    {
        tree[idx].cover = tree[idx].length = 0;
        if (l == r)
        {
            tree[idx].max_length = hbound[l] - hbound[l - 1];
            return;
        }
        int mid = (l + r) / 2;
        init(idx * 2, l, mid);
        init(idx * 2 + 1, mid + 1, r);
        tree[idx].max_length = tree[idx * 2].max_length + tree[idx * 2 + 1].max_length;
    }

    // 区间加减：给 [ul, ur] 的覆盖次数加 diff（矩形左边界 +1、右边界 -1）
    void update(int idx, int l, int r, int ul, int ur, int diff)
    {
        if (l > ur || r < ul)
        {
            return;
        }
        if (ul <= l && r <= ur)
        {
            // 整结点被完全覆盖：只改 cover，无需下推到子结点（懒更新的简化版）
            tree[idx].cover += diff;
            pushup(idx, l, r);
            return;
        }
        int mid = (l + r) / 2;
        update(idx * 2, l, mid, ul, ur, diff);
        update(idx * 2 + 1, mid + 1, r, ul, ur, diff);
        pushup(idx, l, r);
    }

    // 自底向上维护 length：cover>0 → 整段满长；叶子的 cover==0 → 0；否则由子结点累加
    void pushup(int idx, int l, int r)
    {
        if (tree[idx].cover > 0)
        {
            tree[idx].length = tree[idx].max_length;
        }
        else if (l == r)
        {
            tree[idx].length = 0;
        }
        else
        {
            tree[idx].length = tree[idx * 2].length + tree[idx * 2 + 1].length;
        }
    }

private:
    vector<Segtree> tree;
    vector<int> hbound;
    // 作者：力扣官方题解
    // 链接：https://leetcode.cn/problems/rectangle-area-ii/solutions/1825859/ju-xing-mian-ji-ii-by-leetcode-solution-ulqz/
    // 来源：力扣（LeetCode）
    // 著作权归作者所有。商业转载请联系作者获得授权，非商业转载请注明出处。
};

#endif