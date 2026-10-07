#ifndef PLACEDB_H
#define PLACEDB_H
#include "objects.h"
// #include "plot.h"

// ============================================================================
// placedb.h —— PlaceDB：整个布局器的中心数据库（原子对象的「容器 + 管理器」）
//
// 职责：
//   PlaceDB 持有 PlaceCommon/objects.h 里定义的全部原子对象（Module / Net / Pin /
//   SiteRow / Tier / Interval），并对外提供坐标设置、线长统计、密度相关几何量、
//   布局行区间维护、结果输出等操作。Parser 负责往里填数据，QPlace（二次布局）、
//   EPlace（全局布局）、Legalization（合法化）、DetailedPlacement（详细布局）
//   都只通过它读写布局状态 —— 可以理解为「布局器的世界模型」。
//
// 核心数据成员与它们的关系：
//   dbNodes[]      可移动单元：std cell + 可移动 macro。约定 dbNodes[i]->idx == i，
//                  ePlace 的梯度向量也按下标与之一一对应
//   dbTerminals[]  固定 terminal（isFixed=true），含 NI 占位 terminal（isNI）
//   dbPins[]       所有 pin，是 Module 与 Net 之间的枢纽（双向索引）
//   dbNets[]       所有线网，线长（HPWL / WA / LSE）以 net 为单位计算
//   dbSiteRows[]   布局行；每行带 intervals，即被 macro / terminal 打断后剩余的可放置
//                  x 区间（由 removeBlockedSite() 计算，合法化的 subrow 就来自这里）
//   dbTiers[]      3D IC 的分层信息（当前 2D 主流程基本未使用）
//   moduleMap      名字 -> Module* 的反查表，供 .nets / .pl 按名字定位单元
//   coreRegion / chipRegion / commonRowHeight / totalRowArea
//                  由 setCoreRegion() / setChipRegion_2D() 从 dbSiteRows 推导出的全局
//                  几何量；bin 网格划分、密度与电场、合法化吸附全都以它们为基准
//
// 主要函数分组：
//   1. 建库与查表 ：addNode / addTerminal / addNet / addPin、allocate*Memory、
//                   getModuleFromName
//   2. 坐标设置   ：setModuleLocation_* / setModuleCenter_* / setModuleOrientation /
//                   moveModule_2D / randomPlacment / addNoise /
//                   saveNodesLocation / loadNodesLocation
//   3. 线长统计   ：calcHPWL / calcWA_Wirelength_2D / calcLSE_Wirelength_2D /
//                   calcNetBoundPins / calcModuleHPWL / calcModuleHPWLfast
//   4. 布局行维护 ：removeBlockedSite（重算各行 intervals）/ y2RowIndex
//   5. 输出与调试 ：showDBInfo / showRows / outputBookShelf
//                   （内部再调 outputAUX / outputNodes / outputPL / outputNets / outputSCL）
//
// 重要约定（!! 改动前务必确认）：
//   · Module 的 coor / center 是**私有**成员（PlaceDB 是其 friend），外部只能通过
//     setModuleCenter_* / setModuleLocation_* 修改坐标。原因是 pin 只存相对所属
//     module 的 offset，module 一动就必须同步刷新其上所有 pin 的 absolutePos；
//     把写坐标的入口收敛到 PlaceDB，才能保证这件事不被遗漏（否则线长与梯度读到
//     的是过期坐标，且这类错误很难定位）。
//   · 只有 !isFixed 的可动单元会被裁进 coreRegion，terminal / 固定 macro 不裁剪。
//   · 带 _2D 后缀的接口一律忽略 z 分量；3D 相关成员（layerCount、dbTiers、Tier）
//     目前是空壳，只有显式开启 3DIC 相关开关时才会真正生效。
// ============================================================================
class PlaceDB
{
public:
    //! 构造函数只做「清空到未初始化态」：数量类成员置 -1（表示尚未由 Parser 填写），
    //! 容器清空，coreRegion / chipRegion 归零。真正的填充由 Parser 调 add* 系列完成。
    PlaceDB()
    {
        layerCount = -1;
        moduleCount = -1;
        dbMacroCount = 0;
        commonRowHeight = -1;
        dbNodes.clear();
        dbTerminals.clear();
        dbPins.clear();
        dbNets.clear();
        dbSiteRows.clear();
        dbTiers.clear();
        coreRegion = CRect();
        chipRegion = CRect();
        totalRowArea = 0;
        nodesLocationRegister.clear();
    };
    int layerCount;  // how many layers? this is for 3dic
    int moduleCount; // number of modules
    int dbMacroCount;
    int netCount;
    int pinCount;
    double commonRowHeight; //! rowHeight that all(most of the times) rows share
    //! 疑似问题：netCount / pinCount 两个成员在构造函数里没有初始化（未置 -1），
    //!   而且全工程里似乎只有 showDBInfo() 中的同名**局部**变量在用它们，
    //!   即这两个成员既未初始化也可能从未被真正赋值，属于遗留的死成员。

    CRect coreRegion;   // In 2d placement the core region is just the rectangle that encloses all placement rows(see setCoreRegion), in 3d it might be shrunk. coreRegion should be smaller than the whole chip
    CRect chipRegion;   // Chip Region is obtained with coreRegion and all terminal locations. adapect1 is a good example. This should only be used for plot.
    float totalRowArea; //! area of all placement rows, equal or less than coreRegion area, usually equals coreRegion area. calculated in setCoreRegion

    //! dbXxs: vector for storing Xxs
    vector<Module *> dbNodes; // nodes include std cells and macros
    vector<Module *> dbTerminals;
    vector<Pin *> dbPins;
    vector<Net *> dbNets;
    vector<SiteRow> dbSiteRows;
    vector<Tier *> dbTiers;
    //! 注意 dbSiteRows 存的是 SiteRow **值**而非指针，所以 removeBlockedSite() 里
    //! 对 intervals 的修改是直接改容器内的元素；其余容器存的都是指针。

    map<string, Module *> moduleMap; // map module name to module pointer(module include nodes and terminals)

    //! 坐标快照（左下角 coor），用于「备份 → 试算 → 回滚」这类场景（见 save/loadNodesLocation）
    vector<POS_3D> nodesLocationRegister;

    //! —— 第 1 组：建库与区域初始化 ——
    void setCoreRegion();
    void init_tiers();
    //! 疑似问题：init_tiers() 在 placedb.cpp 里是**空实现**，且全工程未见调用点；
    //!   dbTiers 因此始终为空。3D IC 的分层逻辑目前并未真正接通。

    //! add* 系列只负责「按 idx 塞进对应槽位」，容器容量必须先由 allocate*Memory 预分配；
    //! 返回 / 记录的指针由 PlaceDB 持有，外部不要 delete。
    Module *addNode(int index, string name, float width, float height); // (frank) 2022-05-13 consider terminal_NI
    Module *addTerminal(int index, string name, float width, float height, bool isFixed, bool isNI);
    void addNet(Net *);
    int addPin(Module *, Net *, float, float);

    void allocateNodeMemory(int);
    void allocateTerminalMemory(int);
    void allocateNetMemory(int);
    void allocatePinMemory(int);

    //! 按名字反查 Module*；查不到返回 NULL（调用方必须判空，Parser 里曾因此崩过）
    Module *getModuleFromName(string);

    //! —— 第 2 组：坐标设置（唯一合法的写坐标入口）——
    //! 这些函数内部都会：① 必要时把坐标裁进 coreRegion ② 写 module 的 coor/center
    //! ③ 遍历 modulePins 刷新每个 pin 的 absolutePos。三步缺一不可。
    void setModuleLocation_2D(Module *, float, float);
    void setModuleLocation_2D(Module *, POS_3D);
    void setModuleCenter_2D(Module *, float, float);
    void setModuleCenter_2D(Module *, POS_3D);
    void setModuleCenter_2D(Module *, VECTOR_3D);
    void setModuleOrientation(Module *, int);
    void setModuleLocation_2D_random(Module *);
    void moveModule_2D(Module *, VECTOR_2D);
    void moveModule_2D(Module *, VECTOR_2D_INT);
    void randomPlacment(); //! 疑似问题：函数名拼写错误，应为 randomPlacement
    void addNoise();
    void saveNodesLocation();
    void loadNodesLocation();
    CRect getOptimialRegion(Module *); // for legalized results, normally called in detailed placement
    //! 疑似问题：函数名拼写错误，应为 getOptimalRegion。

    //! 只做 coreRegion 裁剪并返回合法中心坐标，**不**真的写回 module（只读版本）
    POS_3D getValidModuleCenter_2D(Module *module, float x, float y);

    //! —— 第 3 组：线长统计 ——
    double calcHPWL();
    double calcWA_Wirelength_2D(VECTOR_2D);
    double calcLSE_Wirelength_2D(VECTOR_2D);
    double calcNetBoundPins();
    double calcModuleHPWL(Module *);
    // double calcModuleHPWLunsafe(Module *);
    double calcModuleHPWLfast(Module *);

    void moveNodesCenterToCenter(); // used for initial 2D quadratic placement

    //! —— 第 4 组：布局行（site row）维护 ——
    void removeBlockedSite(); // calculate intervals of siterows considering macros and terminals that block sites, see void RemoveFixedBlockSite() and void RemoveMacroSite() in ntuplace
    //! 注意调用顺序：必须**先**做 macro 合法化（消除 macro 之间的重叠），再调本函数，
    //! 否则挖出来的可用区间会偏大（见 Legalization/legalizer.cpp 的说明）。

    void setChipRegion_2D();

    //! —— 第 5 组：输出与调试 ——
    void showDBInfo();
    void showRows();

    //! 输出 Bookshelf 格式结果；plOnly=true 时只输出 .pl（坐标变了、其余没变时用）
    //! 疑似问题：声明没有默认参数，但 Legalization/ePlaceAbacus.cpp 里存在
    //!   placedb->outputBookShelf() 的**无参**调用，那个 target 应当编译不过。
    void outputBookShelf(string, bool);

    int y2RowIndex(float);
    bool isConnected(Module *, Module *);

    // void plotCurrentPlacement(string);
private:
    //! 下面五个是 outputBookShelf 的内部实现，分别写一个 Bookshelf 文件：
    //! .aux（索引）/ .nodes（尺寸）/ .pl（坐标）/ .nets（连接关系）/ .scl（布局行）
    void outputAUX();
    void outputNodes();
    void outputPL();
    void outputNets();
    void outputSCL();
};
#endif