// ============================================================================
// 模块总览：DetailedPlacement（详细布局）— 本文件为算法实现
// ============================================================================
// 【职责】
//   布局流程最后一步。输入是**已合法**的布局（Legality：无重叠、对齐 site、
//   落在可放置区间内），输出仍是合法布局，但 HPWL 更低。
//   它只做局部搜索（local search），不做 EPlace 那样的全局解析式优化，
//   每次只改动少量单元的位置，靠反复迭代累积收益。
//
// 【算法思路：局部搜索的整体流程】
//   detailedPlacement() 里就是一个 2 轮的外层循环，每轮依次跑三个算子：
//     1. runLocalReordering() —— 细粒度：同一行内相邻若干个单元（默认窗口 3 个）
//        重新排列顺序，用 DFS + 分支限界求最优。
//     2. runISM()             —— 中粒度：用正方形滑动窗口扫全片，窗口内挑一批
//        互不相连的单元，构造等宽 slot，转成最小权二部匹配（LAP）精确求解。
//     3. runGlobalSwap()      —— 粗粒度：对每个单元，在它自己的 optimal region
//        （由相连 net 的 bounding box 决定）里找空白位或等宽单元做交换。
//   为什么是这个顺序：先做「行内小范围整理」，再做「跨若干行的窗口级重排」，
//   最后做「跨较大距离的搬迁」，粒度由细到粗；每种算子挪动单元的尺度不同，
//   互为补充，跑完一轮后某些区域会出现新的优化机会，所以要迭代多轮。
//
// 【如何评估增益（gain）】
//   三种算子评估的粒度不同，但都以 HPWL 为准：
//     - LocalReordering：以「窗口内相关 net 的 HPWL 之和」为代价，
//       且增量式计算 —— 只有一条 net 上的窗口内单元**全部**排定后才结算它的
//       HPWL（见 LRSolution::createSuccessorSolution / recalculateCost）。
//       这样代价沿搜索路径单调不减，可作为分支限界的下界。
//     - ISM：对每个「单元 i -> slot j」的组合，把单元**真的挪过去**再算
//       calcModuleHPWLfast(i)，填成代价矩阵 cost[i][j]，最后由 LAP 全局最小化总和。
//     - GlobalSwap：直接比较移动前后的 calcModuleHPWLfast 之和
//       （单元交换时同时算两个单元的线长），只要变小就接受。
//   注意：评估用的都是**局部**线长（只算被移动单元相关 net 的 bounding box），
//   因此是全局 HPWL 的一个代理指标，单次看可能有噪声，靠大量迭代收敛。
//
// 【如何保证合法（这是本模块最重要的设计）】
//   本模块**没有**「先移动、再合法化修复」这一步，而是构造性地只产生合法解：
//     - LocalReordering：窗口内单元只换顺序，总宽度和跨度都不变，
//       从最左端 currentX 起依次紧凑排布 -> 必然不重叠、不越界；
//     - ISM：所有 slot 的宽度都被统一成窗口内的最大单元宽度 maxWIDTH
//       （窄单元的右侧会补上一段 tail 空间凑成 maxWIDTH），
//       因此**任意**单元放到**任意** slot 都不会重叠 -> 匹配结果天然合法；
//     - GlobalSwap：只往空白 slot 里插，或只跟**等宽**单元对调 -> 也必然不重叠。
//   所有 slot 坐标都取自 rowSpaces（空闲区间）或已有单元的左下角，
//   而它们本身是对齐 site / row 的，所以对齐约束也自动满足。
//
// 【核心数据结构】
//   - ISMRow（rowSpaces + rowModules）：一行的空闲区间与已放单元，
//     ISM 与 GlobalSwap 共用，是本模块最基础的「合法性账本」；
//     四个原语 insertModule / removeModule / removeTail / insertSpace 负责维护它。
//   - LRSegment / LRSolution / LRSolutionIterator / LRSolver：LocalReordering 的
//     区间、解、后继枚举器与 DFS 搜索器。
//   - lap2：线性指派（LAP）求解器，ISM 的代价矩阵求解引擎。
//   它们都定义在 detailed.h，本文件是各自的实现。
//
// 【主要入口函数】
//   DetailedPlacer::detailedPlacement()   —— 顶层入口（编排三种算子）
//   ISMDP::ISMSweep() / ISMDP::ISMRun()   —— ISM 滑动窗口 / 单次匹配
//   LocalReorderingDP::solve()            —— 局部重排主循环
//   GlobalSwapDP::solve()                 —— 全局交换主循环
//   lap2::lapSolve()                      —— LAP（Jonker-Volgenant 短增广路）
//
// 【注意】本目录未被根 CMakeLists.txt 的 add_subdirectory 包含，当前不参与构建，
//        代码可能与主构建分支存在偏差；阅读时请以「设计意图」为主。
// ============================================================================

#include "detailed.h"

// ============================================================================
// 顶层入口：初始化 + 多轮 { 局部重排 -> ISM -> 全局交换 }
// ============================================================================
// 循环次数写死为 2：三种算子的收益在 2~3 轮后基本饱和，继续迭代性价比很低。
// 顺序上先 LR 再 ISM 再 GS：先收拾行内的局部乱序，再做窗口级最优匹配，
// 最后做跨区域的搬迁；每一轮结束后布局变化会为下一轮创造新的优化空间。
void DetailedPlacer::detailedPlacement()
{
    initialization();
    for (int i = 0; i < 2; i++)
    {
        runLocalReordering();
        runISM();
        runGlobalSwap();
    }
}
// 刷新 SiteRow 的 intervals：把 macro / terminal 挡住的 site 从可放置区间中剔除。
// 这一步必须放在最前面 —— ISM 的 rowSpaces、LR 的 segment、GS 的 rowSpaces
// 全都直接来自 intervals，如果 intervals 里还留有被占的 site，
// 后续就会把单元塞进 macro 身上，直接产生非法结果。
void DetailedPlacer::initialization()
{
    placedb->removeBlockedSite(); // intervals updated for detailed placement here
}

// ============================================================================
// ISM 的驱动循环：反复调用 ISMSweep，直到线长收益低于阈值或超时
// ============================================================================
// 收敛策略（沿用 ntuplace3）：
//   - 每一轮把窗口**逐渐变小**（windowSize 递减）、重叠**逐渐变大**，
//     即从「大范围粗调」过渡到「小范围精修」，类似模拟退火的降温过程；
//   - 当收益（相对上一轮的下降百分比）不再优于 stop 时，打开两个增强开关
//     （doubleWindow + independentCells），试图用更强的手段再榨出一点收益；
//   - 若打开增强开关后收益仍然很差（差于 stop*3），说明已收敛，直接退出。
void DetailedPlacer::runISM()
{
    ISMDP *ismDetailedPlacer = new ISMDP(placedb);
    ismDetailedPlacer->initialization();
    // ismDetailedPlacer->independentCells = false;

    double totalTime;
    time_start(&totalTime);

    double initialHPWL = placedb->calcHPWL();

    double previousHPWL = initialHPWL;

    double detailTimeStart = seconds();

    double totalDetailTime = 0;
    cout << "\nRunning ISM for detailed placement...\n";
    flush(cout);

    // placedb.SaveBlockLocation();

    // the following are 3 adjustable parameters
    double timeLimit = 28800;    // 28800= 8 hours, adjustable
    double stop = -0.2;          // adjustable parameter
    double windowSizeParam = 20; // in row height, adjustable, (20x row height for now)
    // assert(stop < 0);
    // 参数含义：
    //   timeLimit       —— 硬时间上限（默认 8 小时），防止在超大设计上跑飞
    //   stop            —— 「单轮收益阈值」。注意它是**负数**且比较的是
    //                      100*(cur/prev - 1)：线长**下降**时这个值为负，
    //                      所以 "> stop" 的真实语义是「这一轮的下降幅度不足 0.2%」，
    //                      即收益太小、需要换更强的手段（打开增强开关）。
    //   windowSizeParam —— 初始窗口边长（以行高为单位），随轮次递减

    totalDetailTime = seconds() - detailTimeStart;

    double accumulatedTimeStart = seconds(); // donnie 2006-03-13

    int index = 0;
    while (totalDetailTime < timeLimit)
    {
        // placedb.SaveBlockLocation();
        double oneIteTimeStart = seconds();

        // de.MAXWINDOW = param.de_MW; // MAXWINDOW reset here!
        // de.MAXMODULE = param.de_MM;

        // int run_para1 = param.de_window - i; // de_window = 20, this is window size
        // 窗口边长随轮次递减（20 -> 19 -> ...），窗口越小、一次匹配的单元越少，
        // 但覆盖得更密，相当于从全局粗调逐步过渡到局部精修。
        int windowSize = windowSizeParam - index;
        if (windowSize < 15)
        {
            // 缩到 15 以下就不再继续单调递减，而是在 15~19 之间循环，
            // 避免窗口太小导致每轮几乎无事可做（收益为 0 而空转）
            windowSize = 15 + index % 5;
        }

        // 相邻窗口的重叠量随轮次增大（但不超过 8），
        // 重叠越大 -> 窗口滑得越密 -> 边界处的单元被优化的机会越多，
        // 代价是同一个单元可能被重复处理，运行时间变长。
        int windowOverlap = 2 + (int)(index / 2);
        if (windowOverlap > 8)
        {
            windowOverlap = 8;
        }
        // cout << windowSize << " " << windowOverlap << endl;

        ismDetailedPlacer->ISMSweep(windowSize, windowOverlap);

        double currentHPWL = placedb->calcHPWL();
        // double wlx = placedb.GetHPWLp2p();
        double oneIteTime = double(seconds() - oneIteTimeStart);
        double accumulatedTime = seconds() - accumulatedTimeStart;

        printf(" iteration:%2d HPWL=%.0f (%.3f%%)(%.3f%%)   time: %d sec   all: %d sec\n",
               index, currentHPWL, 100.0 * (currentHPWL / previousHPWL - 1.0), 100.0 * (currentHPWL / initialHPWL - 1.0),
               (int)oneIteTime, (int)accumulatedTime);

        fflush(stdout);

        // 已经开了增强开关（independentCells）但收益仍然差于 stop*3
        // （即 0.6%），说明确实收敛了，退出主循环。
        // 要求 independentCells 已打开，是为了避免第一轮随机波动就误判收敛。
        if ((100.0 * (currentHPWL / previousHPWL - 1.0)) > stop * 3 && ismDetailedPlacer->independentCells)
        {
            break;
        }

        // 收益差于 stop（下降不足 0.2%）-> 打开增强开关，换更强的手段再试
        if ((100.0 * (currentHPWL / previousHPWL - 1.0)) > stop) // by donnie
        {
            if (ismDetailedPlacer->doubleWindow == false)
            {
                ismDetailedPlacer->doubleWindow = true;
                ismDetailedPlacer->independentCells = true;
                cout << "Consider independent cells and enabled double window now\n"; // ! start independent
            }
        }

        totalDetailTime = seconds() - detailTimeStart;
        index++;
        previousHPWL = currentHPWL;
    }
    time_end(&totalTime);
    double finalHPWL = placedb->calcHPWL();
    cout << "ISM result: Pin-to-pin HPWL= " << finalHPWL << " (" << 100.0 * (finalHPWL / initialHPWL - 1.0) << "%)\n";
    cout << "ISM total runtime:" << totalTime << "\n";
}

// ============================================================================
// 全局交换的驱动：建一次 GSRows，然后跑 2 轮 solve()
// ============================================================================
// 注意 initialization() 只调用一次而 solve() 调用两次：
// GSRows 是在 solve() 过程中**增量维护**的（每次移动都同步 insertModule/removeModule），
// 所以第二轮可以直接复用。这样省掉了重建全部行结构的开销。
void DetailedPlacer::runGlobalSwap()
{
    double GStime;
    GlobalSwapDP *GSDetailedPlacer = new GlobalSwapDP(placedb);
    GSDetailedPlacer->initialization();
    cout << "\nRunning Global Swap for detailed placement...\n";
    time_start(&GStime);
    for (int i = 0; i < 2; i++)
    {
        GSDetailedPlacer->solve();
    }
    time_end(&GStime);
    cout << "\nGlobal Swap finished in " << GStime << " seconds\n";
    //! 疑似问题：GSDetailedPlacer / ISM 的 ismDetailedPlacer / LR 的 lrDetailedPlacer
    //! 都用 new 创建却从未 delete，每轮循环都泄漏一次（对象本身不大，
    //! 但 ISMRows / GSRows 里的 map 会随行数累积，跑多轮会持续吃内存）。
}

// ============================================================================
// 局部重排的驱动：构建 segments 后跑 2 轮 solve(3, 2, 1)
// ============================================================================
// solve(3, 2, 1) = 窗口 3 个单元、相邻窗口重叠 2 个（即每次右移 1 个）、扫描 1 轮。
// 窗口取 3 是因为排列数 3! = 6，枚举代价可以接受；窗口再大就需要更强的剪枝。
void DetailedPlacer::runLocalReordering()
{
    LocalReorderingDP *lrDetailedPlacer = new LocalReorderingDP(placedb);
    lrDetailedPlacer->initialization();
    for (int i = 0; i < 2; i++)
    {
        cout << "\nRunning LR for detailed placement...\n";
        lrDetailedPlacer->solve(3, 2, 1);
        // lrDetailedPlacer->solve(3, 2, 1);
    }
    //! 疑似问题：与 runGlobalSwap 同理，lrDetailedPlacer 泄漏；
    //! 另外两轮之间没有重新 initializeSegments()，segModules 的顺序是
    //! solve() 结束后被 solveForBestOrder 就地更新的，所以顺序仍然有效。
}

// ============================================================================
// ISMSweep：用正方形窗口扫过整个芯片（对应 ntuplace3 的 grid_run()）
// ============================================================================
// 窗口是**正方形**的，边长 = windowSize 个行高（x 方向也用行高做单位，
// 这是因为标准单元的高宽比接近 1:1，用行高当尺子最自然）。
// 相邻窗口之间重叠 windowOverlap 个行高：位于窗口边界的单元如果只属于一个窗口，
// 它的可选 slot 就只有半边，优化效果会很差；重叠能让边界单元在下一个窗口里
// 处于中心位置，从而被充分优化。
void ISMDP::ISMSweep(int windowSize, int windowOverlap) // a square window, width=height=windowSize
{
    // 重叠量必须小于窗口边长，否则窗口步长 <= 0 会原地不动（死循环）。
    // 退化处理：干脆不重叠。
    if (windowOverlap >= windowSize)
    {
        windowOverlap = 0;
    }

    // fplan->CalcHPWL();
    // double wl1 = fplan->GetHPWLp2p();

    // 把芯片在 x / y 两个方向都按行高切成格子，窗口就是这些格子组成的方阵
    int x_num = (int)(placeDB->coreRegion.getWidth() / (placeDB->commonRowHeight));
    int y_num = (int)(placeDB->coreRegion.getHeight() / (placeDB->commonRowHeight));
    //	cout<<"\n total bins:"<<total<<"\n";
    //! 疑似问题：最后一个窗口可能超出 coreRegion（i + windowSize 可能 > y_num），
    //! 窗口右/上边界没有被裁剪到 coreRegion 之内。ISMRun() 内部对 y 方向做了裁剪
    //! （endRowBottom 与 coreRegion.ur.y 取 min），但 x 方向的 window.ur.x 没有裁剪，
    //! 可能导致 lower_bound(window.ur.x) 取到行末、把行尾单元也纳入匹配。
    for (int i = 0; i < y_num; i = i + windowSize - windowOverlap) // window slide across the chip region, overlap means overlap between windows
    {
        for (int j = 0; j < x_num; j = j + windowSize - windowOverlap)
        {
            CRect curWindow;

            curWindow.ll.x = placeDB->coreRegion.ll.x + j * placeDB->commonRowHeight;
            curWindow.ll.y = placeDB->coreRegion.ll.y + i * placeDB->commonRowHeight;
            curWindow.ur.x = curWindow.ll.x + windowSize * placeDB->commonRowHeight;
            curWindow.ur.y = curWindow.ll.y + windowSize * placeDB->commonRowHeight;

            // cout << curWindow;
            // ISMRun(fplan->m_coreRgn.left + j * placeDB->commonRowHeight, fplan->m_coreRgn.bottom + i * placeDB->commonRowHeight,
            //        window * placeDB->commonRowHeight, window * placeDB->commonRowHeight);
            ISMRun(curWindow);
        }
    }
}

// ============================================================================
// ISMRun：对一个窗口执行一次完整的独立集匹配（ISM 的核心）
// ============================================================================
// 五步流程（对应函数体内 1~5 的注释）：
//   1. 收集窗口内的单元，按宽度降序排好（moduleList）
//   2. 挑单元 + 造 slot：为每个被选中的单元准备一个**等宽（maxWIDTH）**的 slot；
//      若单元不足 maxModuleCount 个，再从窗口的空白区间里补一些空 slot
//   3. 填代价矩阵：cost[i][j] = 单元 i 放到 slot j 后的 HPWL
//   4. 跑 LAP（lap2::lapSolve）求最小代价的完美匹配
//   5. 按匹配结果把单元重新插入各行（insertModule），并写回 placedb 坐标
//
// 【为什么这样就能保证合法 —— 本函数最关键的一点】
//   slot 的宽度一律取窗口内最大单元的宽度 maxWIDTH：
//     - 宽度 == maxWIDTH 的单元：slot 就是它自己原来的位置；
//     - 宽度 <  maxWIDTH 的单元：slot 起点仍是它原来的位置，
//       但**占用**右侧紧邻的空白（tail），使 slot 总宽也凑成 maxWIDTH
//       （这一步要求右侧空白恰好紧邻，否则这个单元本轮不选）。
//   于是所有 slot 同宽，任何单元放进任何 slot 都不会越界、不会重叠，
//   「重排问题」因此退化成「无约束的指派问题」，LAP 的最优解可直接采用。
//   加空 slot 是为了让单元有机会搬到纯空白区域（相当于允许局部疏散）。
//
// 【为什么要独立集（independentCells）】
//   代价矩阵里每一项是**单个单元**的 HPWL，LAP 最小化的是它们的**和**。
//   只有当这些单元两两之间没有共享 net 时，各自线长才互不影响，
//   「求和最小」才真正等价于「窗口总 HPWL 最小」。否则匹配只是启发式近似。
//   所以 independentCells 打开后会逐个检查 isConnected() 来筛掉相连单元。
//
// 【关于 removeTail / insertSpace 的时序】
//   选出窄单元时立刻 removeTail 把它右侧那段空白标记为「已占用」，
//   这样第 2 步枚举空白 slot 时就不会把这段空间重复算成可用空白（否则会重复计数）；
//   等空白 slot 枚举完毕，再用 insertSpace 把这段空间还回去。
//   注意被选中单元腾出的位置是**在空白 slot 枚举之后**才 removeModule 的，
//   所以本轮里被选中单元的旧位置只以 slot 的形式出现一次，不会被重复计入空白。
void ISMDP::ISMRun(CRect window)
{

    // 1.build modules size set in the window
    // 2.build slot and modules vector, if it's length > maxModuleCount , cut it.
    // 3.calc cost matrix
    // 4.run bipartite matching
    // 5.save result

    // 窗口覆盖的行范围。注意 endRowBottom 减了一个行高：
    // y2RowIndex 算的是「y 落在第几行」，而我们要的是窗口**顶边所在的那一行**，
    // 若不加这个修正就会多算一行（让窗口实际比预期高一行）。
    int startRowIndex = max(placeDB->y2RowIndex(window.ll.y), 0); // lowest row

    double endRowBottom = min(double(window.ur.y), placeDB->coreRegion.ur.y - placeDB->commonRowHeight);
    int endRowIndex = placeDB->y2RowIndex(endRowBottom); // highest row

    // multimap 以宽度为 key，插入即自动按宽度升序排好；
    // 下面再用 rbegin()/rend() 反序遍历，就得到「宽度降序」的单元序列。
    // 之所以要从宽到窄处理：slot 宽度 = 最宽单元的宽度，
    // 先处理宽的能尽早确定 maxWIDTH，窄单元则靠补 tail 来凑宽度。
    multimap<double, Module *> moduleMap; // use multimap for auto sorting (by module width)

    list<Module *> moduleList; // list<module ID>

    // build module size map, with all modules in the window
    for (int i = startRowIndex; i <= endRowIndex; i++) // i< endRowIndex in ntuplace
    {
        // map<double, int>::iterator iter;// curModule
        // map<double, int>::iterator startModule, endModule;// begin and end in the original code

        auto startModule = ISMRows[i].rowModules.lower_bound(window.ll.x); //! rowModules initialized in initializeISMRows()
        auto endModule = ISMRows[i].rowModules.lower_bound(window.ur.x);

        for (auto curModuleIter = startModule; curModuleIter != endModule; curModuleIter++)
        {
            // if (curModuleIter->second->getHeight() == placeDB->commonRowHeight)
            if (!(curModuleIter->second->isMacro)) // macros are ignored
            {                                      //! adding all cells in the window
                moduleMap.insert(pair<double, Module *>(curModuleIter->second->getWidth(), curModuleIter->second));
            }
        }
        // ISMRows[i].showModule();
        // cout << "modulecountcp1: " << moduleMap.size() << endl;
    }

    // double window：除了当前窗口，再从右侧第 5~6 个窗口宽处捞一批单元进来一起匹配。
    // 目的：当本窗口收益变小时，扩大候选池（相当于把搜索范围跳到更远处），
    // 让匹配有更大的优化空间。这是 ntuplace 的一个提分技巧。
    //! 疑似问题：这里的循环上界是 i < endRowIndex（不含），而上面收集主窗口单元
    //! 用的是 i <= endRowIndex（含），两者不一致，最后一行不会被 double window 考虑。
    //! 另外被捞进来的单元其 slot 仍是它自己原来的位置，但那些位置**不在**本窗口内，
    //! 匹配后单元可能跨窗口移动很远，这里没有任何距离约束。
    if (doubleWindow == true) // enable double window, pick cells from 2 windows
    {
        for (int i = startRowIndex; i < endRowIndex; i++)
        {
            map<double, int>::iterator iter;
            // modified by Jin 20070726
            auto startModule = ISMRows[i].rowModules.lower_bound(window.ur.x + 5 * window.getWidth()); //?? double window?
            auto endModule = ISMRows[i].rowModules.lower_bound(window.ur.x + 6 * window.getWidth());   //??
            for (auto curModuleIter = startModule; curModuleIter != endModule; curModuleIter++)
            {
                // if (curModuleIter->second->getHeight() == placeDB->commonRowHeight)
                if (!(curModuleIter->second->isMacro)) // macros are ignored
                {
                    moduleMap.insert(pair<double, Module *>(curModuleIter->second->getWidth(), curModuleIter->second));
                }
            }
        }
    }

    //! moduleList stores ALL modules in current window(windows), sorted by module width, decreaingly (sorted increasingly in moduleMap)
    for (auto rIter = moduleMap.rbegin(); rIter != moduleMap.rend(); rIter++)
    {
        moduleList.push_back(rIter->second);
    }

    // 外层循环：一次匹配未必能消化掉窗口里所有单元，
    // 每轮从中挑出至多 maxModuleCount 个做一次 LAP，剩下的留到下一轮继续。
    // 循环终止条件：moduleList 被清空（每轮至少会选中 1 个单元 —— 最宽的那个
    // 总能被选中，因为它的 slot 就是自己的位置，一定放得下）。
    while (moduleList.size() != 0)
    {
        // cout<<"entered\n";
        // slots and modules: assign all modules to all slots, slots.size()==modules.size()?
        vector<POS_2D> slots;
        vector<Module *> modules;
        slots.reserve(maxModuleCount);
        modules.reserve(maxModuleCount);

        int insertedModuleCount = 0;
        // slot 的统一宽度 = 窗口内最宽单元的宽度。moduleList 按宽度降序排好，
        // 所以首元素就是最宽的。这是「任意指派都合法」的前提。
        const double maxWIDTH = (*moduleList.begin())->getWidth(); // WIDTH == largest width

        // The following 2 vectors are connected: removedRowIndexes[i] = the row that removedModules[i] is located in
        vector<Module *> removedModules;
        removedModules.reserve(maxModuleCount);

        vector<int> removedRowIndexes;
        removedRowIndexes.reserve(maxModuleCount);

        // The following 3 vectors are connected:
        // empty: space left after inserting a cell with width less than maxWIDTH, emptyX[i]: start of the empty,
        // emptyWidth[i]: length of the empty, emptyRowIndex[i]: row in which the ith empty located

        vector<double> emptyX;
        emptyX.reserve(maxModuleCount);

        vector<double> emptyWidth;
        emptyWidth.reserve(maxModuleCount);

        vector<int> emptyRowID;
        emptyRowID.reserve(maxModuleCount);

        for (auto curModuleIter = moduleList.begin(); curModuleIter != moduleList.end();)
        {
            Module *curModule = *(curModuleIter);
            POS_2D curModulePos = curModule->getLL_2D();
            bool inserted = false;
            int rowID = placeDB->y2RowIndex(curModulePos.y);

            // 情况 A：最宽的单元 —— 它自己的位置天然就是一个 maxWIDTH 宽的 slot
            if (curModule->getWidth() == maxWIDTH)
            {
                bool connection = false;
                if (independentCells == true) // independent set
                {
                    // 逐个检查它与已选中的单元是否共享 net；
                    // 只要与其中任何一个相连就不能同时选（否则代价不可加）
                    for (Module *storedModule : modules)
                    {
                        if (placeDB->isConnected(storedModule, curModule) == true)
                        {
                            connection = true;
                        }
                    }
                }
                if (connection == false)
                {
                    modules.push_back(curModule);
                    POS_2D slot; //?? no need to copy here, delete later
                    slot.x = curModulePos.x;
                    slot.y = curModulePos.y;
                    slots.push_back(slot);
                    // remove_module(*iter,rowID);
                    removedModules.push_back(curModule);
                    removedRowIndexes.push_back(rowID);
                    inserted = true;
                }
            }
            // 情况 B：较窄的单元 —— 需要它右侧**恰好紧邻**一段空白，
            // 用「自身宽度 + 右侧空白」凑出一个 maxWIDTH 宽的 slot
            else if (curModule->getWidth() < maxWIDTH) // width < WIDTH, need to check if there is space for a slot of maxWIDTH
            {
                auto spaceNextToCell = ISMRows[rowID].rowSpaces.upper_bound(curModulePos.x);
                if (spaceNextToCell != ISMRows[rowID].rowSpaces.end())
                {                                                                          // todo: understand when would the following 'if' happen, in eplace and ntuplace
                    if (spaceNextToCell->first == curModule->getWidth() + curModulePos.x)  //?? when would this happen for a legalized design???? if this won't happen, then ISM actually only consider cells with same size?
                    {                                                                      //? this happens only for the last cell of a group of abut cells? we can add the first cell too? what if we adopt chris chu's method: ignore overlap and do a legalization after swaps?
                        if ((curModule->getWidth() + spaceNextToCell->second) >= maxWIDTH) // the space occupied by the module + the space next to the module can hold a module of maximum size(maxWIDTH)
                        {
                            bool connection = false;
                            if (independentCells == true)
                            {
                                // for (int m = 0; m < (int)modules.size(); m++)
                                // {
                                //     if (placeDB->isConnected(modules[m], curModule) == true)
                                //     {
                                //         connection = true;
                                //     }
                                // }
                                for (Module *storedModule : modules)
                                {
                                    if (placeDB->isConnected(storedModule, curModule) == true)
                                    {
                                        connection = true;
                                    }
                                }
                            }
                            if (connection == false)
                            {
                                // tail = 单元右侧要被一并占用、用来把 slot 补到 maxWIDTH 的那段空白
                                double tailX = curModule->getWidth() + curModulePos.x;
                                double tailWidth = maxWIDTH - curModule->getWidth(); // space left EMPTY after removing the cell with width less than maxWidth

                                modules.push_back(curModule);
                                POS_2D slot;
                                slot.x = curModulePos.x;
                                slot.y = curModulePos.y;

                                slots.push_back(slot);
                                // remove_module(*iter,rowID);
                                removedModules.push_back(curModule);
                                removedRowIndexes.push_back(rowID);
                                emptyX.push_back(tailX); //! need to understand empty! 2024.8.21
                                emptyWidth.push_back(tailWidth);
                                emptyRowID.push_back(rowID);
                                // 先把 tail 从「可用空白」里扣掉，避免下面枚举空白 slot 时
                                // 把这段已经预订给本 slot 的空间又算成一份空白（重复计数）。
                                // 等空白 slot 枚举完后（见下方的 insertSpace 循环）再还回去。
                                ISMRows[rowID].removeTail(tailX, tailWidth); //!!! remove tail here, but add space later, why remove tail? for whitespace handling?
                                inserted = true;
                            }
                        }
                    }
                }
            }

            if (inserted == true)
            {
                auto moduleIter2 = curModuleIter;
                curModuleIter++;
                moduleList.erase(moduleIter2);
                insertedModuleCount++;
                if (insertedModuleCount >= maxModuleCount)
                {
                    break;
                }
            }
            else
            {
                curModuleIter++;
            }
        }
        // Handle whitespace and add extra slots
        // 若本轮选中的单元还不够多，就从窗口内的空白区间里再切一些空 slot 出来。
        // 空 slot 对应的 modules[i] == NULL，它的代价恒为 0（见下面代价矩阵），
        // 代表「这个位置空着不站人」—— 于是单元既可以在已占位置之间互换，
        // 也可以搬到纯空白区域，搜索空间更大。

        //! 疑似问题：这里的判据是 insertedModuleCount（**单元**数）与
        //! maxWindow（**slot** 总数上界，90）比较，而单元数的截断上界其实是
        //! maxModuleCount（128）。两个不同含义的量混用了。
        //! 另外下面内层 break 之后，外层两个 for 仍会继续跑并继续 push_back 新 slot，
        //! emptySlotCount 可以远超 maxWindow，起不到限流作用（break 只跳出最内层）。
        //! 还有一点：这段用的是 i < endRowIndex，与上面收集单元时的 i <= endRowIndex
        //! 不一致，窗口最上面那一行的空白被漏掉了。
        if (insertedModuleCount < maxWindow)
        {
            int emptySlotCount = 0;
            for (int i = startRowIndex; i < endRowIndex; i++)
            {
                for (auto curRowSpaceIter = ISMRows[i].rowSpaces.lower_bound(window.ll.x);
                     curRowSpaceIter != ISMRows[i].rowSpaces.lower_bound(window.ur.x);
                     curRowSpaceIter++)
                {
                    // 一段空白最多能切出 floor(空白长度 / maxWIDTH) 个等宽 slot
                    int slotCount = (int)(curRowSpaceIter->second / maxWIDTH);
                    for (int j = 0; j < slotCount; j++)
                    {
                        modules.push_back(NULL); // is this ok?
                        POS_2D slot;
                        slot.x = curRowSpaceIter->first + j * maxWIDTH;
                        slot.y = ISMRows[i].ll.y;

                        assert(placeDB->coreRegion.ll.y + i * placeDB->commonRowHeight == ISMRows[i].ll.y);

                        slots.push_back(slot); //! how was the vector 'position'used?
                        emptySlotCount++;
                        if (emptySlotCount + insertedModuleCount >= maxWindow)
                        {
                            break;
                        }
                    }
                }
            }
        }

        // 把上一步为凑 slot 而临时扣掉的 tail 空白还回去
        // （它们已经完成「占位」使命，避免被重复计为可用空白）
        for (int i = 0; i < emptyX.size(); i++)
        {
            ISMRows[emptyRowID[i]].insertSpace(emptyX[i], emptyWidth[i]); // reinsert tail
        }
        // remove all selected modules, and then assign all modules (stored in the vector modules) to all available positions, which are stored in the vector position
        // 把所有被选中的单元从行上摘下来。放在空白 slot 枚举**之后**是关键的时序：
        // 这样它们腾出的位置不会又被当成空白 slot 再生成一遍（否则同一个位置
        // 会同时出现在「单元 slot」和「空白 slot」里，匹配后可能产生重叠）。
        for (int i = 0; i < removedModules.size(); i++)
        {
            ISMRows[removedRowIndexes[i]].removeModule(removedModules[i]);
        }

        // 3. eatablish bimatching object, calc cost matrix
        // 到这里：modules 与 slots 一一对应，deg = slot 总数。
        // 注意 slots.size() == modules.size()，LAP 求的是它们之间的一个**双射**。
        int deg = modules.size(); // matrix degree
        // bimatching matrix(deg);
        lap2 matrix(deg);

        // calc cost matrix
        for (int i = 0; i < deg; i++)
        {
            // NULL 表示「这个 slot 空着」的虚拟单元：它放到哪里代价都是 0，
            // 于是 LAP 可以让真正的单元自由选择「留在原位」还是「搬到空白处」。
            if (modules[i] == NULL) // (donnie) whitespace?
            {
                for (int j = 0; j < deg; j++)
                {
                    matrix.put(i, j, 0);
                }
            }
            else
            {
                // bool set=false;
                POS_2D bestPos = slots[0];
                double bestWirelength = DOUBLE_MAX;
                // double ox=fplan->m_modules[modules[i]].m_x;
                // double oy=fplan->m_modules[modules[i]].m_y;
                // CNetLengthCalc netcalc(*fplan, modules[i]);
                // netcalc.init();

                // 代价矩阵的填充方式很「朴素」但很直观：
                // 把单元真的挪到 slot j 上，再算它相关 net 的 HPWL。
                // 因为所有单元都在**当前**的 placedb 坐标下评估，
                // 所以先算的单元会受后面单元位置的影响 —— 这正是为什么
                // 独立集（互不相连）很重要：相连单元的线长会互相耦合，
                // 矩阵里的和可加性就不成立了。
                for (int j = 0; j < deg; j++)
                {
                    placeDB->setModuleLocation_2D(modules[i], slots[j].x, slots[j].y);
                    double wl = placeDB->calcModuleHPWLfast(modules[i]);
                    // cerr<<" "<<wl;
                    // double wl=0;
                    // for(int k=0;k<(int)fplan->m_modules[modules[i]].m_netsId.size();k++)
                    //{
                    //   wl = wl+fplan->GetNetLength(fplan->m_modules[modules[i]].m_netsId[k]);
                    // }
                    // cerr<<","<<wl;
                    // if(wl!=wl2)
                    //	cerr<<" "<<wl<<","<<wl2;
                    // force module find empty space ^^  ==>FOOL!
                    // if( (modules[j]!=modules[i])&&(modules[j]!=-1))
                    //	wl=wl+100000;
                    if (wl < bestWirelength)
                    {
                        bestWirelength = wl;
                        bestPos = slots[j];
                    }
                    // matrix.costs.put(i,j,wl);
                    matrix.put(i, j, wl);
                }
                // fplan->SetModuleLocation(modules[i],ox,oy);
                // 扫完所有 slot 后，把单元暂时放在它自己「单看最优」的那个 slot 上，
                // 再继续算下一个单元的代价矩阵。
                //! 疑似问题：这样提前落子会让**后续**单元的 HPWL 计算基于
                //! 「前面单元已移到贪心最优位」的坐标，而最终 LAP 给出的指派
                //! 未必是这个贪心位置，于是代价矩阵与实际落地结果不一致
                //! （矩阵是在一组「中间态」坐标下评估出来的）。
                //! 更严谨的做法是先全部算完、最后统一按 LAP 结果落位。
                //! 另外若 deg==0 或其它异常导致 bestPos 未被更新，这里会用到 slots[0]。
                placeDB->setModuleLocation_2D(modules[i], bestPos.x, bestPos.y); //? is it necessary to set now? 2024.8.27
                                                                                 //? 2024.11.1 setting now would affect the bestWirelength calculation of the following cells?
            }
        }

        //	matrix.show();

        // 4.run matching
        // matrix.find();
        matrix.lapSolve();
        vector<int> result;
        matrix.getResult(result);

        // 5.save result
        // result[i] = 第 i 个单元被指派到的 slot 下标。
        // 因为所有 slot 等宽且互不重叠，这里直接插入即可，理论上不会失败；
        // insertModule 返回 false 只可能是 rowSpaces 账目与实际不一致（见下方提示）。
        for (int i = 0; i < deg; i++)
        {
            // NULL 单元代表空 slot，没有实体可放；
            // 又因为 NULL 是追加在 modules 尾部的，遇到第一个即可 break。
            if (modules[i] == NULL)
            {
                break;
            }

            int pos = placeDB->y2RowIndex(slots[result[i]].y);

            if (!ISMRows[pos].insertModule(slots[result[i]].x, modules[i])) //!!module reinserted here!!
            {
                // 插入失败 = 位置放不下，说明 rowSpaces 与实际占用已经失配。
                //! 疑似问题：这里失败后**不回滚**、也不撤销 setModuleLocation_2D，
                //! 会留下一个「placedb 里已挪走但行结构上没记录」的单元，
                //! 后续所有基于 rowSpaces 的判断都会继续出错（静默地破坏合法性）。
                if (!gArg.CheckExist("nocheck")) // (donnie) 2007-03-25
                {
                    // cout << "\ninsert" << pos << " x:" << position[result[i]].x
                    //      << " w:" << fplan->m_modules[modules[i]].m_width
                    //      << " mod:" << modules[i]
                    //      << " name:" << fplan->m_modules[modules[i]].m_name << "\n";
                    // m_de_row[pos].showspace();
                }
            }
            placeDB->setModuleLocation_2D(modules[i], slots[result[i]].x, slots[result[i]].y);
        }
        // for(int i=0;i<deg;i++)
        //{
        //	int pos=y2rowID(position[i].y);
        //	if(!m_de_row[pos].insert_module(position[i].x,fplan->m_modules[modules[i]].m_width,modules[i]))
        //	{
        //		cout<<"\ninsert"<<pos<<" x:"<<position[i].x<<" w:"<<fplan->m_modules[i].m_width<<" mod:"<<modules[i]<<" name:"<<fplan->m_modules[modules[i]].m_name<<"\n";
        //		m_de_row[pos].showspace();
        //	}
        //	fplan->SetModuleLocation(modules[i],position[result[i]].x,position[result[i]].y);
        // }
        // cout<<"lefted\n";
    }
}

// ISM 的初始化：先设参数，再建行结构。
// 注意这里**不**再调用 removeBlockedSite() —— 那个已经由 DetailedPlacer::initialization()
// 做过一次；重复调用会重复剔除区间，属于冗余（原作者把它注释掉了）。
void ISMDP::initialization()
{
    // placeDB->removeBlockedSite();
    initializeParams();
    initializeISMRows();
}

// 参数取值直接沿用 ntuplace3 的经验值。
// maxModuleCount=128 意味着一次 LAP 的规模上限约 128，O(n^3) 尚可接受；
// 再大匹配耗时会急剧上升，所以必须截断。
void ISMDP::initializeParams()
{
    // see ntuplace3, not sure if maxWindow=90, maxModuleCount=128, or maxWindow=maxModuleCount=64
    doubleWindow = false;
    independentCells = false;
    maxWindow = 90;
    maxModuleCount = 128;
}

// ============================================================================
// 构建 ISMRows：先把每行的可放置区间填进 rowSpaces，再把已有单元插进去
// ============================================================================
// 顺序很重要：必须**先**铺好全部空白区间，**后**插单元，
// 因为 insertModule() 的语义就是「从某个空白区间里切走一段」。
// 如果先插单元再铺空白，单元占用的部分就会被误当成可用空白。
void ISMDP::initializeISMRows()
{
    double ISMRowLength = placeDB->coreRegion.getWidth();
    // cout << "initialize ISM rows: \n";

    // 1.create ISMRows
    for (SiteRow curRow : placeDB->dbSiteRows)
    {
        ISMRow newRow(placeDB->coreRegion.ll.x, curRow.bottom, ISMRowLength);
        for (Interval curInterval : curRow.intervals)
        {
            // cout<<"father: "<<curInterval.start<<" "<<curInterval.getLength()<<endl;
            newRow.rowSpaces[curInterval.start] = curInterval.getLength();
            // cout<<"son: "<<newRow.rowSpaces[curInterval.start]<<endl;
        }
        ISMRows.push_back(newRow);
    }

    // 2. insert cells to ISMRows
    for (Module *curModule : placeDB->dbNodes) // all nodes, macros and cells
    {

        int coveredRowCount = 1;
        POS_2D curModulePos = curModule->getLL_2D();
        float moduleHeight = curModule->getHeight();
        // 高度超过一行高的必然是 macro（多倍行高），它会在**每一行**都占位置，
        // 所以要在它覆盖的每一行里都插一份记录。标准单元 coveredRowCount == 1。
        if (moduleHeight > placeDB->commonRowHeight) //?? macro included???
        {
            coveredRowCount = (int)(moduleHeight / placeDB->commonRowHeight);
            if (placeDB->commonRowHeight * coveredRowCount < moduleHeight)
            {
                coveredRowCount++; // for macros
            }
        }

        int pos = placeDB->y2RowIndex(curModulePos.y);

        //! 疑似问题：这里直接用 pos + i 访问 ISMRows，没有检查 pos + i 是否越界。
        //! 若 macro 顶部超出 coreRegion（或 y2RowIndex 的取整有偏差），会越界访问。
        //! 另外 macro 在 ISM / GlobalSwap 里后面又被 isMacro 过滤掉、从不参与匹配，
        //! 但它们在每一行都占掉了 rowSpaces，属于「只占位不参与」的合理设计。
        for (int i = 0; i < coveredRowCount; i++)
        {
            ISMRows[pos + i].insertModule(curModulePos.x, curModule);
        }
    }

    // placeDB->showRows();

    // for (ISMRow curRow : ISMRows)
    // {
    //     curRow.showSpace();
    //     curRow.showModule();
    // }
}

// ============================================================================
// ISMRow::insertModule —— 把一个单元放进本行（合法性的守门员）
// ============================================================================
// 做法：在 rowSpaces（按起点排序的空闲区间表）里找到「起点 <= x 的最后一个区间」，
// 检查它从 x 开始是否还剩下 >= moduleWidth 的长度；
// 放得下就把这段区间切出去（可能切成左右两截），并在 rowModules 里登记。
// 返回 false 表示放不下 —— 调用方据此判定失败，绝不允许强行插入造成重叠。
bool ISMRow::insertModule(double x, Module *curModule)
{
    //! both moduleX and moduleWidth should be integer for a legalized placement
    // 先 round 成整数：合法化之后所有坐标/宽度都应是 site step 的整数倍，
    // 取整可以消掉浮点累积误差导致的「差一点点放不下」这类假失败。
    double moduleX = round(x);
    double moduleWidth = round(curModule->getWidth());
    //@modified by Jin 20070727

    auto rowSpaceIter = rowSpaces.upper_bound(moduleX);
    if (rowSpaceIter == rowSpaces.begin())
    {
        // cout<<"wtf?\n";
        // upper_bound 返回 begin() 有两种含义：
        //   (a) 所有空闲区间的起点都 > x（x 太靠左）；
        //   (b) rowSpaces 为空（begin() == end()）。
        // 合法布局下都不该发生。
        return false; // x<all empty site's start point, shouldn't happen for a legalized placement?
    }
    else
    {
        --rowSpaceIter; // upper_bound returns the first element greater than key
        // 该区间从 x 处起算的剩余长度 = second - (x - first)，必须容得下 moduleWidth
        if ((rowSpaceIter->second - (moduleX - rowSpaceIter->first)) < moduleWidth)
        {
            // cout<<rowSpaceIter->second<<" "<<moduleX<<" "<<rowSpaceIter->first<<" "<<moduleWidth<<endl;
            return false; // space can't contain this cell
        }
        else
        {
            // l2 = 区间左端到单元左端的距离，即插入后被切下来的**左半段**长度
            double l2 = moduleX - rowSpaceIter->first;
            if (l2 == 0)
            {
                // 单元正好贴着区间左端放：左半段长度为 0，只需保留右半段（若有）
                if (((rowSpaceIter->first + rowSpaceIter->second) - (moduleX + moduleWidth)) != 0)
                {
                    rowSpaces[moduleX + moduleWidth] = (rowSpaceIter->first + rowSpaceIter->second) - (moduleX + moduleWidth);
                }
                rowSpaces.erase(rowSpaceIter->first); // space updated here. space initilization in the constructor didn't consider modules, here the spaces are modified considering cells
            }
            else
            {

                // 单元从区间中间放：右边可能还剩一截，单独登记；原区间缩短为左半段 l2
                // 注意「先登记右半段、再改左半段长度」的顺序：
                // 若先改长度，后面的计算就用不到原来的右端了。
                if ((moduleX + moduleWidth) < (rowSpaceIter->first + rowSpaceIter->second))
                {
                    rowSpaces[moduleX + moduleWidth] = (rowSpaceIter->first + rowSpaceIter->second) - (moduleX + moduleWidth);
                }
                rowSpaceIter->second = l2; // update length of space, which is l2.
            }
            rowModules[moduleX] = curModule;
            // cout<<"inserted\n"<<curModule->idx<<endl;
            return true;
        }
    }
}

// ============================================================================
// ISMRow::removeModule —— 摘走一个单元
// ============================================================================
// 两步：(1) 把它占的 x 区间还给 rowSpaces（并做相邻区间合并）；
//       (2) 从 rowModules 里删掉登记项。
// ISM 里被选中的单元会先全部摘下，等 LAP 出结果后再重新 insertModule，
// 所以「摘下」必须把空间完整归还，否则后面的 slot 会算错。
void ISMRow::removeModule(Module *module)
{

    //?SHOULD NOT use double as the index of map

    //! 疑似问题：insertSpace() 恒返回 true，其返回值赋给 flag 后**从未被使用**
    //! （死变量）。也就是说即便空间归还失败，这里也不会有任何提示。
    bool flag = insertSpace(module->getLL_2D().x, module->getWidth());
    double targetModuleX = round(module->getLL_2D().x);

    auto curRowModuleIter = rowModules.find(targetModuleX);
    if (curRowModuleIter == rowModules.end())
    {
        printf("Warning: mID %d not found\n", module->idx);
    }
    else if (curRowModuleIter->second != module)
    {
        printf("Warning: mID %d does not equal ite->second %d\n", module->idx, curRowModuleIter->second);
    }
    else
    {
        rowModules.erase(curRowModuleIter); // cell will be reinserted after ISM
    }
    //@modified by Jin 20070727
}

// ============================================================================
// ISMRow::removeTail —— 把一段**原本空闲**的区间标记为已占用
// ============================================================================
// 语义上等价于「从 rowSpaces 里挖走 [x, x+width)」，与 insertSpace 正好互为逆操作。
// 名字叫 removeTail 是因为这段空间通常是某个窄单元右侧（tail）的空白，
// ISM 用它与单元自身宽度拼出一个 maxWIDTH 的 slot。
// 逻辑与 insertModule 的「切区间」部分完全一致，只是不登记 rowModules。
bool ISMRow::removeTail(double x, double width)
{

    x = round(x);
    width = round(width);

    auto spaceIter = rowSpaces.upper_bound(x);
    if (spaceIter == rowSpaces.begin())
        return false; // x<all empty site's start point

    spaceIter--;
    // 区间从 x 起算的剩余长度不够 -> 这段 tail 并不真实存在，无需移除
    if ((spaceIter->second - (x - spaceIter->first)) < width) // (iter->first+iter->second)-(x+w) < 0, x+w is not convered by the empty indexed by iter, tail doesn't exist, no need to remove
        return false;

    double l2 = x - spaceIter->first;

    if (l2 == 0) // x == iter->first
    {
        // tail 贴着区间左端：挖走后只剩右半段
        if (((spaceIter->first + spaceIter->second) - (x + width)) != 0) // (iter->first+iter->second)-(x+w) > 0
        {
            rowSpaces[x + width] = (spaceIter->first + spaceIter->second) - (x + width);
        }
        rowSpaces.erase(spaceIter->first); // tail is removed, space after tail is inserted.
    }
    else
    { // what is this situation????: a long space is divided by the tail, and we remove tail (so the tail is no longer empty)
        // 这种情况是：tail 位于某个较长空白的**中间**，
        // 挖走后原空白被 tail 切成左右两截（左截 l2，右截另行登记）。
        if ((x + width) < (spaceIter->first + spaceIter->second))
        {

            rowSpaces[x + width] = (spaceIter->first + spaceIter->second) - (x + width);
        }
        spaceIter->second = l2;
    }

    return true;
}

// ============================================================================
// ISMRow::insertSpace —— 归还一段空闲区间，并与相邻空闲区间合并
// ============================================================================
// 为什么要合并：rowSpaces 里若出现两条首尾相接的区间，
// 「能切出几个 maxWIDTH 的 slot」这类按长度整除的判断就会把它们当成两块小碎片，
// 白白浪费本来可用的连续空间；同时也会让 upper_bound / lower_bound 的语义变复杂。
// 因此归还空间时必须左右各看一眼，能并就并。
// 四种组合：
//   front==true  && end==true  ：左右都不相接 -> 直接新插一条
//   front==true  && end==false ：只有右边相接 -> 并到右边（右侧区间改起点为 x）
//   front==false && end==true  ：只有左边相接 -> 并到左边（左侧区间加长）
//   front==false && end==false ：左右都相接   -> 三段并成一段（本段把左右粘起来）
bool ISMRow::insertSpace(double x, double width)
{
    // merge a spaces with existing spaces

    x = round(x);
    width = round(width);

    bool front; // true: no space in front of (and abut) the space to be added
    bool end;   // true: no space behind after (and abut) the space to be added

    auto curSpaceIter = rowSpaces.lower_bound(x);
    if (curSpaceIter == rowSpaces.begin())
    {
        front = true;
    }
    else
    {
        curSpaceIter--;
        if ((curSpaceIter->first + curSpaceIter->second) == x)
        {
            front = false;
        }
        else
        {
            front = true;
        }
    }
    curSpaceIter = rowSpaces.lower_bound(x);

    if (curSpaceIter == rowSpaces.end())
    {
        end = true;
    }
    else
    {
        if (curSpaceIter->first == (x + width))
        {
            end = false;
        }
        else
        {
            end = true;
        }
    }

    if ((front == true) && (end == true))
    {
        // 孤立的一段：直接插入
        rowSpaces[x] = width;
    }
    else if ((front == true) && (end == false))
    {
        // 与右边相接：新建一条以 x 为起点、长度为「本段 + 右段」的记录，再删掉右段
        curSpaceIter = rowSpaces.lower_bound(x);
        rowSpaces[x] = width + curSpaceIter->second;
        rowSpaces.erase(curSpaceIter->first); // or erase(curSpaceIter) ?
    }
    else if ((front == false) && (end == true))
    {
        // 与左边相接：左段原地加长即可
        curSpaceIter = rowSpaces.lower_bound(x);
        curSpaceIter--;
        rowSpaces[curSpaceIter->first] = curSpaceIter->second + width;
    }
    else // front==false && end == false
    {
        // 左右都相接：本段充当「粘合剂」，把左段、本段、右段合成一条
        // 注意必须先取出右段的 first/second，再 erase，否则迭代器失效后取不到值
        curSpaceIter = rowSpaces.lower_bound(x);
        double spaceBehindX = curSpaceIter->first;
        double spaceBehindWidth = curSpaceIter->second;
        curSpaceIter--;
        rowSpaces[curSpaceIter->first] = curSpaceIter->second + width + spaceBehindWidth;
        rowSpaces.erase(spaceBehindX);
    }
    return true;
}

// 调试用：打印本行的空闲区间（[起点, 长度] 序列）
void ISMRow::showSpace()
{
    cout << "\n=====ISM ROW SPACE ===\n";

    for (auto iter = rowSpaces.begin(); iter != rowSpaces.end(); iter++)
    {
        // modified by Jin 20070727
        printf("[%.10f,%.10f] ", iter->first, iter->second);
        // cout<<" ["<<iter->first<<","<<iter->second<<"] ";
        // modified by Jin 20070727
    }
    cout << '\n';
}

// 调试用：打印本行已放单元（[左下角 x, 单元 idx] 序列）
void ISMRow::showModule()
{
    cout << "\n=====ISM ROW Module ===\n";

    for (auto iter = rowModules.begin(); iter != rowModules.end(); iter++)
    {
        cout << " [" << iter->first << "," << iter->second->idx << "] ";
    }
    cout << '\n';
}

// ============================================================================
// lap2::lapSolve —— 线性指派问题（LAP）求解器
// ============================================================================
// 求解 min sum_i cost[i][assignment[i]]，assignment 是一个双射。
// 在 ISM 中：行 = 单元（NULL 表示空 slot），列 = slot，代价 = HPWL。
//
// 算法：Jonker-Volgenant 的短增广路（Shortest Augmenting Path）方法，O(n^3)。
// 它维护对偶变量 v[j]（列价格），使得归约代价 cost[i][j] - v[j] 满足
// 「每个已匹配行的最小归约代价为 0」这一互补松弛条件；然后对每个尚未匹配的行
// 用类 Dijkstra 的过程找一条最短增广路，把匹配数 +1，并更新对偶变量。
// 整体分为 6 段（见下方各段注释）：
//   1. COLUMN REDUCTION           列归约：为每列找最小代价行，得到初始匹配与 v
//   2. REDUCTION TRANSFER         归约转移：把只匹配一次的行做行归约，并入 v
//   3. AUGMENTING ROW REDUCTION   增广行归约：扫描自由行，尽量直接赋值（跑 2 遍）
//   4. AUGMENT SOLUTION           对剩余自由行逐个用 Dijkstra 找最短增广路
//   5. 更新列价格 v
//   6. 沿增广路翻转匹配，最后计算最优总代价
// 变量名（u/v/d/pred/free/collist/matches/colsol）沿用经典 LAPJV 实现，
// 便于和原始文献/参考代码对照。
int lap2::lapSolve()
{
    // input:
    // degree        - problem size
    // cost - cost matrix

    // output:
    // assignment     - column assigned to row in solution
    // assignment     - row assigned to column in solution
    // u          - dual variables, row reduction numbers
    // v          - dual variables, column reduction numbers
    //! 疑似问题：注释里提到输出对偶变量 u / v，但 v 是局部变量、u 虽声明却
    //! 全程未被使用（死变量）；两者都不会返回给调用方。
    bool unassignedfound;
    int i, imin, numfree = 0, prvnumfree, f, i0, k, freerow;
    int j, j1, j2, endofpath, last = 0, low, up;
    int min = 0, h, umin, usubmin, v2;
    j2 = 0;
    endofpath = 0;

    vector<int> free;
    free.resize(degree, 0);
    vector<int> collist;
    collist.resize(degree, 0);
    vector<int> matches;
    matches.resize(degree, 0);
    vector<int> d;
    d.resize(degree, 0);
    vector<int> pred;
    pred.resize(degree, 0);

    vector<int> colsol;
    colsol.resize(degree, 0);

    vector<int> u;
    u.resize(degree, 0);

    vector<int> v;
    v.resize(degree, 0);

    // ---- 1. COLUMN REDUCTION（列归约）----
    // 对每一列取最小代价作为该列的初始对偶变量 v[j]，
    // 若该最小代价所在行尚未被占用，就先把这一列指派给它（贪心建立初始匹配）。
    // 采用**逆序**遍历列是 LAPJV 的经验做法：这样先处理后面的列，
    // 冲突分布更均匀，初始匹配质量更好。
    for (j = degree - 1; j >= 0; j--) // reverse order gives better results.
    {
        // find minimum cost over rows.
        min = cost[0][j];
        imin = 0;
        for (i = 1; i < degree; i++)
            if (cost[i][j] < min)
            {
                min = cost[i][j];
                imin = i;
            }
        v[j] = min;

        if (++matches[imin] == 1)
        {
            // init assignment if minimum row assigned for first time.
            assignment[imin] = j;
            colsol[j] = imin;
        }
        else
            colsol[j] = -1; // row already assigned, column not assigned.
    }

    // ---- 2. REDUCTION TRANSFER（归约转移）----
    // 分成两类行：
    //   matches[i]==0 -> 一轮都没抢到列，放进 free 列表，留给后面的增广阶段处理；
    //   matches[i]==1 -> 恰好占了一列，对它做**行归约**：
    //      找出该行「除已占列之外」的最小归约代价 min，把 v[j1] 减去 min。
    //      这一步把行上的对偶信息转移到列上，为后续增广提供更好的初始价格。
    for (i = 0; i < degree; i++)
        if (matches[i] == 0) // fill list of unassigned 'free' rows.
            free[numfree++] = i;
        else if (matches[i] == 1) // transfer reduction from rows that are assigned once.
        {
            j1 = assignment[i];
            min = INF;
            for (j = 0; j < degree; j++)
                if (j != j1)
                    if (cost[i][j] - v[j] < min)
                        min = cost[i][j] - v[j];
            v[j1] = v[j1] - min;
        }

    // ---- 3. AUGMENTING ROW REDUCTION（增广行归约）----
    // 这一阶段试图「就地」把每个自由行安排掉：
    //   对自由行 i 找归约代价最小的列 j1 和次小的列 j2；
    //   若 j1 已被别的行 i0 占着，就让 i 抢占 j1，把被挤掉的 i0 拿回 free 列表
    //   （可能原地重试，也可能留到下一阶段）——这就是「增广」的雏形。
    // 为什么要跑两遍（loopcnt < 2）：
    //   第一遍中被挤掉的行会重新入队，原地重试往往又能安排掉一部分；
    //   实践中两遍就能把绝大多数行安排完，剩下的才交给代价更高的 Dijkstra 阶段。
    int loopcnt = 0; // do-loop to be done twice.
    do
    {
        loopcnt++;

        // scan all free rows.
        // in some cases, a free row may be replaced with another one to be scanned next.
        k = 0;
        prvnumfree = numfree;
        numfree = 0; // start list of rows still free after augmenting row reduction.
        while (k < prvnumfree)
        {
            i = free[k];
            k++;

            // find minimum and second minimum reduced cost over columns.
            umin = cost[i][0] - v[0];
            j1 = 0;
            usubmin = INF;
            for (j = 1; j < degree; j++)
            {
                h = cost[i][j] - v[j];
                if (h < usubmin)
                    if (h >= umin)
                    {
                        usubmin = h;
                        j2 = j;
                    }
                    else
                    {
                        usubmin = umin;
                        umin = h;
                        j2 = j1;
                        j1 = j;
                    }
            }

            // 若 j1 已被占，i0 就是当前占着 j1 的那一行
            i0 = colsol[j1];
            if (umin < usubmin)
                // 最小列严格优于次小列：抬高（对偶意义下降低）j1 的价格，
                // 使本行在 j1 上的归约代价升到次小值，为后面的增广留空间。
                // change the reduction of the minimum column to increase the minimum
                // reduced cost in the row to the subminimum.
                v[j1] = v[j1] - (usubmin - umin);
            else             // minimum and subminimum equal.
                if (i0 >= 0) // minimum column j1 is assigned.
                {
                    // swap columns j1 and j2, as j2 may be unassigned.
                    j1 = j2;
                    i0 = colsol[j2];
                }

            // (re-)assign i to j1, possibly de-assigning an i0.
            assignment[i] = j1;
            colsol[j1] = i;

            if (i0 >= 0) // minimum column j1 assigned earlier.
                if (umin < usubmin)
                    // put in current k, and go back to that k.
                    // continue augmenting path i - j1 with i0.
                    free[--k] = i0;
                else
                    // no further augmenting reduction possible.
                    // store i0 in list of free rows for next phase.
                    free[numfree++] = i0;
        }
    } while (loopcnt < 2); // repeat once.

    // ---- 4. AUGMENT SOLUTION（为每个剩下的自由行做增广）----
    // 第 3 阶段没能安排掉的行，必须真正找一条「增广路」：
    //   从自由行 freerow 出发，在「行-列」交替的图上找一条到某个**未匹配列**的
    //   最短路；沿这条路把匹配关系整体翻转，就能让匹配数 +1。
    // 这里用类 Dijkstra 的实现：d[j] 是到列 j 的最短距离，pred[j] 记录前驱行，
    // collist 是待扫描列的列表，low/up 把 collist 划成「已就绪 / 待扫描 / 未触及」三段。
    for (f = 0; f < numfree; f++)
    {
        freerow = free[f]; // start row of augmenting path.

        // Dijkstra shortest path algorithm.
        // runs until unassigned column added to shortest path tree.
        for (j = 0; j < degree; j++)
        {
            d[j] = cost[freerow][j] - v[j];
            pred[j] = freerow;
            collist[j] = j; // init column list.
        }

        // collist 三段划分（LAPJV 的经典技巧，避免真正的优先队列）：
        //   [0, low)      已确定最短距离、不再变动
        //   [low, up)     当前距离等于最小值的列，待扫描
        //   [up, degree)  尚未触及
        low = 0; // columns in 0..low-1 are ready, now none.
        up = 0;  // columns in low..up-1 are to be scanned for current minimum, now none.
                 // columns in up..degree-1 are to be considered later to find new minimum,
                 // at this stage the list simply contains all columns
        unassignedfound = false;
        do
        {
            if (up == low) // no more columns to be scanned for current minimum.
            {
                last = low - 1;

                // scan columns for up..degree-1 to find all indices for which new minimum occurs.
                // store these indices between low..up-1 (increasing up).
                min = d[collist[up++]];
                for (k = up; k < degree; k++)
                {
                    j = collist[k];
                    h = d[j];
                    if (h <= min)
                    {
                        if (h < min) // new minimum.
                        {
                            up = low; // restart list at index low.
                            min = h;
                        }
                        // new index with same minimum, put on undex up, and extend list.
                        collist[k] = collist[up];
                        collist[up++] = j;
                    }
                }

                // check if any of the minimum columns happens to be unassigned.
                // if so, we have an augmenting path right away.
                for (k = low; k < up; k++)
                    if (colsol[collist[k]] < 0)
                    {
                        endofpath = collist[k];
                        unassignedfound = true;
                        break;
                    }
            }

            if (!unassignedfound)
            {
                // update 'distances' between freerow and all unscanned columns, via next scanned column.
                j1 = collist[low];
                low++;
                i = colsol[j1];
                h = cost[i][j1] - v[j1] - min;

                for (k = up; k < degree; k++)
                {
                    j = collist[k];
                    v2 = cost[i][j] - v[j] - h;
                    if (v2 < d[j])
                    {
                        pred[j] = i;
                        if (v2 == min) // new column found at same minimum value
                            if (colsol[j] < 0)
                            {
                                // if unassigned, shortest augmenting path is complete.
                                endofpath = j;
                                unassignedfound = true;
                                break;
                            }
                            // else add to list to be scanned right away.
                            else
                            {
                                collist[k] = collist[up];
                                collist[up++] = j;
                            }
                        d[j] = v2;
                    }
                }
            }
        } while (!unassignedfound);

        // ---- 5. 更新列价格 ----
        // 对最短路树上已就绪的列，按最短距离与当前最小值的差调整对偶变量，
        // 使下一轮增广仍从对偶可行的状态出发。
        // update column prices.
        for (k = 0; k <= last; k++)
        {
            j1 = collist[k];
            v[j1] = v[j1] + d[j1] - min;
        }

        // ---- 6. 沿增广路翻转匹配 ----
        // 从终点 endofpath（那个未匹配的列）沿 pred 回溯到起点 freerow，
        // 把路径上「已匹配/未匹配」的关系整体取反，匹配数因此 +1。
        // reset row and column assignments along the alternating path.
        do
        {
            i = pred[endofpath];
            colsol[endofpath] = i;
            j1 = endofpath;
            endofpath = assignment[i];
            assignment[i] = j1;
        } while (i != freerow);
    }

    // calculate optimal cost.
    // 按最终 assignment 累加原始代价（注意不是归约代价），得到真实最优总代价。
    int lapcost = 0;
    for (i = 0; i < degree; i++)
    {
        j = assignment[i];
        lapcost = lapcost + cost[i][j];
    }
    //  this->assignment=assignment;
    //! 疑似问题：lapcost 用 int 累加 degree 个 int，代价很大时存在溢出风险；
    //! 且返回值目前被 ISMRun 忽略（只关心 assignment 本身）。

    return lapcost;
}

// 局部重排的初始化：只做一件事 —— 构建 segments（区间 + 挂在上面的有序单元）
void LocalReorderingDP::initialization()
{
    initializeSegments();
}

// ============================================================================
// LocalReorderingDP::solve —— 局部重排主循环
// ============================================================================
// 三层嵌套：
//   轮数 iterationNumber -> 遍历每个 segment -> 在 segment 内滑动窗口
// 窗口内单元按 x 升序排列，取连续 windowSize 个做最优排列；
// 相邻窗口重叠 overlapSize 个单元，保证跨窗口边界的单元也有机会被一起重排。
// 每个窗口求出的最优排列会**立即写回** placedb 和 segModules，
// 所以后续窗口看到的是更新后的坐标（类似 Gauss-Seidel 式的就地更新）。
void LocalReorderingDP::solve(int windowSize, int overlapSize, int iterationNumber)
{
    for (int i = 0; i < iterationNumber; i++)
    {
        for (auto curLRSegment = lrSegments.begin(); curLRSegment != lrSegments.end(); curLRSegment++)
        {
            // 1. No cells in this segment
            if (curLRSegment->segModules.size() < 1)
            {
                continue;
            }

            // 2. Only one cell in this segment: check the solution packing the cell to the left or right
            // 只有一个单元时没有「排列」可言，退化为三个候选位置里选线长最小的：
            // 原地不动 / 紧贴 segment 左端 / 紧贴 segment 右端。
            // 这三个位置都是合法的（唯一单元在本 segment 内怎么挪都不可能与别人重叠），
            // 且贴边往往能缩短与相邻 segment 上单元的距离。
            else if (curLRSegment->segModules.size() == 1)
            {
                Module *onlyModule = curLRSegment->segModules.front();
                POS_2D modulePos = onlyModule->getLL_2D();

                double originalX = modulePos.x;
                double originalWirelength = placeDB->calcModuleHPWLfast(onlyModule);

                double leftX = curLRSegment->start;
                placeDB->setModuleLocation_2D(onlyModule, leftX, modulePos.y);
                double leftWirelength = placeDB->calcModuleHPWLfast(onlyModule);

                double rightX = curLRSegment->end - onlyModule->getWidth();
                placeDB->setModuleLocation_2D(onlyModule, rightX, modulePos.y);
                double rightWirelength = placeDB->calcModuleHPWLfast(onlyModule);

                // OrigW is the smallest
                if (originalWirelength <= leftWirelength && originalWirelength <= rightWirelength)
                {
                    placeDB->setModuleLocation_2D(onlyModule, originalX, modulePos.y);
                }
                // leftW is the smallest
                else if (leftWirelength <= originalWirelength && leftWirelength <= rightWirelength)
                {
                    placeDB->setModuleLocation_2D(onlyModule, leftX, modulePos.y);
                }
                // rightW is the smallest
                else
                {
                    placeDB->setModuleLocation_2D(onlyModule, rightX, modulePos.y);
                }

                continue;
            }

            // 单元数多于一个窗口：滑动窗口逐段处理
            if (curLRSegment->segModules.size() > windowSize) // number of modules in the segment larger than window size
            {
                // 步长 = windowSize - overlapSize。默认 (3,2) 即每次右移 1 个单元，
                // 于是每个单元都会出现在 windowSize 个不同的窗口里，优化机会均等。
                //! 疑似问题：循环条件是 startModuleIte + windowSize <= end()，
                //! 因此**末尾不足 windowSize 个**的那几个单元永远不会被重排。
                for (auto startModuleIte = curLRSegment->segModules.begin(); startModuleIte + windowSize <= curLRSegment->segModules.end(); startModuleIte = startModuleIte + (windowSize - overlapSize))
                {
                    //! 疑似问题：bImprove 是**死变量** —— 赋值后从未被读取，
                    //! 既不影响滑动策略也不影响终止条件。
                    bool bImprove = solveForBestOrder(startModuleIte, startModuleIte + windowSize);
                }
            }
            else // curLRSegment->segModules.size() > windowSize
            {
                // 整个 segment 的单元数不超过窗口大小：一次性全部参与重排
                //! 疑似问题：这里的注释写反了（仍写 "> windowSize"），实际是 <= windowSize 的分支。
                bool bImprove = solveForBestOrder(curLRSegment->segModules.begin(), curLRSegment->segModules.end());
            }
        }
    }
}

// ============================================================================
// initializeSegments —— 由 SiteRow::intervals 建 segment，并把单元挂上去排好序
// ============================================================================
// 前提：placedb->removeBlockedSite() 已经把 macro / terminal 占用的 site 从
// intervals 中剔除，所以每个 interval 都是一段「连续、干净、可自由摆放」的空间。
// 之后 segModules 会按 x 升序排序 —— 这是 solve() 里「窗口 = 连续下标区间」
// 这一假设成立的前提，也是重排不会打乱整体左右关系的基础。
void LocalReorderingDP::initializeSegments()
{
    //! 1. create segments from intervals
    for (SiteRow curRow : placeDB->dbSiteRows)
    {
        for (Interval curInterval : curRow.intervals)
        {
            lrSegments.push_back(LRSegment(curRow.bottom, curInterval.start, curInterval.end));
        }
    }

    //! 2. add modules to segments
    //  printf("Init BB...");
    fflush(stdout);

    // Add all non-Macro modules into segments
    for (Module *curModule : placeDB->dbNodes)
    {
        if (curModule->isMacro) // Skip Macro modules
        {
            continue;
        }
        POS_2D curModulePos = curModule->getLL_2D();

        // Find the corresponding segment of this module, and insert the module id into the segment
        // 用单元的 (y, x, x+width) 造一个「探针 segment」，二分定位它落在哪个 segment。
        //! 疑似问题（原作者已标 //!!!! potential bug）：compareLRSegment 不是严格的
        //! 字典序（bottom 相同时只用 end 比较，没有直接比 start），
        //! lower_bound 可能定位到错误位置；定位失败时这里只打印警告并**跳过该单元**，
        //! 该单元于是永远不会参与局部重排（静默地少优化一批单元）。
        LRSegment compSeg(curModulePos.y, curModulePos.x, curModulePos.x + curModule->getWidth());
        auto iteFindSegment = lower_bound(lrSegments.begin(), lrSegments.end(), compSeg, compareLRSegment()); //!!!! potential bug
        if (iteFindSegment != lrSegments.end() && iteFindSegment->bottom == curModulePos.y)
        {
            iteFindSegment->addModule(curModule);
        }
        else
        {
            fprintf(stderr, "Warning: initializeSegments(), cannot find legal segment for "
                            "module '%s' at (%.2f, %.2f) w: %.2f h: %.2f\n",
                    curModule->name.c_str(),
                    curModulePos.x, curModulePos.y, curModule->getWidth(), curModule->getHeight());
            if (iteFindSegment != lrSegments.end())
            {
                fprintf(stderr, "   iteFindSegment  : bottom %.2f left: %.2f right: %.2f\n", iteFindSegment->bottom, iteFindSegment->start, iteFindSegment->end);
                fprintf(stderr, "   iteFindSegment-1: bottom %.2f left: %.2f right: %.2f\n", (iteFindSegment - 1)->bottom, (iteFindSegment - 1)->start, (iteFindSegment - 1)->end);
            }
        }

        // Warning: this module is not on any segments
        if (iteFindSegment == lrSegments.end())
        {
            cerr << "Warning: Module " << curModule->idx << " is not on any segments" << endl;
        }
    }

    // printf( "..s.." );
    fflush(stdout);
    // Sort the module id's by their x coordinates (increasingly)
    // 按 x 升序排序：之后「窗口」就可以简单地用连续下标区间 [i, i+windowSize) 表示，
    // 且窗口内的单元本来就是物理上相邻的一小段。
    //! 疑似问题：原注释提到 lambda 比较器有潜在 bug。此处用 float_less 做严格小于
    //! 比较，若两个单元 x 完全相等（例如同一位置的重叠，合法化后理论上不该发生）
    //! 会返回 false，排序结果不确定，但不算 UB。
    for (auto iteSegment = lrSegments.begin(); iteSegment != lrSegments.end(); iteSegment++) //!! potential bug regarding the lambda expression for comparing Module*
    {
        sort(iteSegment->segModules.begin(), iteSegment->segModules.end(), [=](Module *a, Module *b)
             { return float_less(a->getLL_2D().x, b->getLL_2D().x); });
    }

    // printf("done\n");
    fflush(stdout);
}

// ============================================================================
// solveForBestOrder —— 对 [startModule, endModule) 这段相邻单元求最优排列
// ============================================================================
// 步骤：
//   1) 构造 originalSolution：单元按**当前顺序**排好，算出它的真实 HPWL 作为初始上界。
//      （这一步很关键：保证搜索结果一定不比现状差，找不到更好的就保持原样。）
//   2) 构造 initialSolution：所有单元都还在 uninsertedCells 里，
//      并把窗口内所有空白合并成一个 NULL 单元一起参与排列，
//      currentX 设为窗口最左端。
//   3) LRSolver 做 DFS + 分支限界枚举全部排列（含空白的插入位置）。
//   4) 把最优排列写回 placedb 坐标，并同步更新 segModules 里的指针顺序。
// 合法性来源：所有排列都从 currentX 起依次紧凑排布，总宽度与窗口跨度都不变，
//   因此既不重叠也不越界；空白被当作一个整体单元参与排列，位置自然合法。
bool LocalReorderingDP::solveForBestOrder(vector<Module *>::iterator startModule, vector<Module *>::iterator endModule)
{
    // The original solution is the bound of CellSwap
    LRSolution originalSolution(placeDB);

    // // test code
    // double left_bound = m_pDB->m_modules[*startModule].m_x;
    // double right_bound = m_pDB->m_modules[*(endModule - 1)].m_x +
    //                      m_pDB->m_modules[*(endModule - 1)].m_width;
    // //@test code

    for (auto itePushItem = startModule; itePushItem != endModule; itePushItem++)
    {
        originalSolution.insertedCells.push_back(*itePushItem);
        originalSolution.xLocations.push_back((*itePushItem)->getLL_2D().x);
    }
    originalSolution.recalculateCost();

    LRSolution initialSolution(placeDB);
    initialSolution.solutionCost = 0;
    for (auto itePushItem = startModule; itePushItem != endModule; itePushItem++)
    {
        initialSolution.uninsertedCells.push_back(*itePushItem);
    }

    initialSolution.initializeNetModuleCount();

    float minX = FLOAT_MAX;
    float maxX = -FLOAT_MAX; //!
    float totalWidth = 0;

    // 求窗口的 x 跨度 [minX, maxX] 与单元总宽度。
    // 注意这里用的是**单元当前的实际坐标**，而不是 segment 的 start/end：
    // 窗口跨度只覆盖被选中单元本身占据的范围，重排严格限制在这个范围内，
    // 绝不会把单元挤出它本来所在的那块地盘。
    for (Module *uninsertedModule : initialSolution.uninsertedCells)
    {
        POS_2D curPos = uninsertedModule->getLL_2D();
        minX = min(minX, curPos.x);
        maxX = max(maxX, curPos.x + uninsertedModule->getWidth());
        totalWidth += uninsertedModule->getWidth();
    }

    // Set currentX
    // 排列从最左端开始，逐个往右紧凑码放
    initialSolution.currentX = minX;

    // If there is free space in current range, add a white space module
    // 窗口内若有零散空白，把它们**合并成一个** NULL 单元参与排列。
    // 合并的理由：空白内部再怎么切分，对线长没有任何影响，
    // 拆成多段只会让排列数爆炸（搜索空间变大）而收益为零。
    //! 疑似问题：合并成单块后，空白只能整体出现在排列的某一个位置，
    //! 无法「一部分在左、一部分在右」，某些最优布局会被漏掉（是剪枝近似，非精确）。
    if (totalWidth < maxX - minX)
    {
        initialSolution.whiteSpaceWidth = maxX - minX - totalWidth;
        initialSolution.uninsertedCells.push_back(NULL); // pack all free space into one single whitespace
    }

    LRSolver solver;
    // 当前实际布局作为初始上界：这样「找不到更好的」时结果就是原布局，绝不退化
    solver.setBestSolution(&originalSolution);
    // cout << "initial best: " << endl;
    originalSolution.print();
    LRSolution bestSol = solver.solve(&initialSolution);

    // Update the module coordinates
    // 按最优排列把每个单元写到它的新 x 上，并同步更新 segModules 的**指针顺序**
    // （iteUpdate 指向 segModules 中窗口对应的那段），
    // 这样下一轮滑动窗口看到的顺序仍然与 x 坐标一致。
    //! 疑似问题：DFS 过程中单元被真实移动过（createSuccessorSolution 里
    //! setModuleLocation_2D），而**被剪枝掉的分支不会回滚**这些移动；
    //! 虽然这里最后按 bestSol 统一写回，窗口内单元是一致的，
    //! 但被剪枝分支留下的「中间态坐标」会污染搜索期间的 HPWL 评估
    //! （calcNetHPWL 看的是 placedb 的当前坐标，包含那些没被采用的临时位置）。
    auto iteUpdate = startModule;
    auto iteList = bestSol.insertedCells.begin();
    auto iteLocation = bestSol.xLocations.begin();

    while (iteList != bestSol.insertedCells.end())
    {
        // NULL 是空白，没有实体可写，但它对应的 xLocations 也要一起跳过
        if (*iteList != NULL)
        {
            Module *curModule = *iteList;
            POS_2D curPos = curModule->getLL_2D();
            // Update coordiantes
            placeDB->setModuleLocation_2D(*iteList, *iteLocation, curPos.y);
            // Update module order
            *iteUpdate = *iteList;
            iteUpdate++;
        }

        iteList++;
        iteLocation++;
    }

    // 与初始上界（原布局）比较，报告本窗口是否真的有改善。
    // 注意：无论是否改善，坐标都已经按 bestSol 写回；
    // 若 bestSol 就是 originalSolution，则等价于原地不动。
    if (bestSol.solutionCost < originalSolution.solutionCost)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// ============================================================================
// LRSolver::solve —— 搜索入口：从 initialSolution 深搜，返回最优完整解
// ============================================================================
LRSolution LRSolver::solve(LRSolution *initialSolution)
{

    depthFirstSolve(initialSolution);
    // 只有当前上界是「原始布局」时才会出现 bestSolution 为空的情况
    // （即初始上界都没设置成功）；这里直接退出进程而非抛异常，
    // 属于早期的粗暴处理。
    if (bestSolution == NULL)
    {
        // throw domain_error ("no feasible solution found");
        printf("no feasible solution found");
        exit(0);
    }
    return *bestSolution;
}

// 设置初始上界。只接受**完整解**（部分解的代价不是真实 HPWL，不能当上界）。
void LRSolver::setBestSolution(LRSolution *curSolution)
{
    if (curSolution->isComplete())
    {
        bestCost = curSolution->getCost();
        // 必须 clone：调用方（solveForBestOrder）的 originalSolution 是栈对象，
        // 出作用域就析构了，直接存指针会悬空。
        bestSolution = curSolution->clone(); ///!!!!
    }
    //! 疑似问题：若传入的是部分解，这里会**静默什么都不做**，
    //! bestSolution 保持为 NULL，直到 solve() 里才被发现。
}

void LRSolver::updateBestSolution(LRSolution *curSolution)
{
    if (curSolution->getCost() < bestCost)
    {

        delete bestSolution;
        bestSolution = curSolution->clone();
        bestCost = curSolution->getCost();
    }
}

// ============================================================================
// LRSolver::depthFirstSolve —— DFS + 分支限界主体
// ============================================================================
// 递归结构：
//   若当前解已完整 -> 用它的代价更新最优解；
//   否则 -> 枚举每个「把某个未排单元放到 currentX」的后继，
//           只有后继代价仍严格小于当前上界时才递归进去（剪枝）。
// 剪枝为何正确：solutionCost 只累计「窗口内单元已全部排定的 net」的 HPWL，
//   随着排列推进只会往里加项，因此它单调不减，是完整解代价的一个合法下界；
//   一旦下界已经 >= bestCost，这棵子树里任何完整解都不可能更优。
void LRSolver::depthFirstSolve(LRSolution *curSolution)
{
    if (curSolution->isComplete())
    {
        updateBestSolution(curSolution);
        // cout << "one complete solution: " << endl;
        curSolution->print();
    }
    else
    {

        LRSolutionIterator *i = new LRSolutionIterator(curSolution);
        while (!i->isDone())
        {
            LRSolution *successor = i->createSuccessorSolution();

            // successor.Print();
            // 剪枝：代价下界已不优于当前最优 -> 整棵子树放弃
            //! 疑似问题：原注释对这里的剪枝条件存疑。注意 successor 是**部分解**，
            //! 它的 Cost 只是已完整 net 的线长之和，确实单调不减，可作下界；
            //! 但 bestCost 在搜索过程中会被更新，而已经被展开的兄弟分支
            //! 不会重新检查，属于正常行为。
            if (successor->getCost() < bestCost) // pruning here, but successor.Bound() is always less than bestObjective????
            {
                // cout << "adopted: " << endl;
                successor->print();
                depthFirstSolve(successor);
            }
            else
            {
                // cout << "not adopted: " << endl;
                successor->print();
            }

            //! 疑似问题：successor 被 delete，但它在构造时已经**真实移动过**
            //! placedb 里的单元坐标（见 createSuccessorSolution），
            //! 这些副作用不会被回滚。被剪枝分支留下的坐标会一直存在到
            //! 该窗口最终按 bestSol 写回为止，期间会污染其它 net 的 HPWL 计算。
            delete successor; //!!!
            i->moduleIterator++;
        }
        delete i; //!!!
    }
}

// 深拷贝一份解。LRSolver 保存最优解时必须 clone，
// 因为搜索中的临时解会被立刻 delete（见 depthFirstSolve）。
LRSolution *LRSolution::clone()
{
    return new LRSolution(*this);
}

// 调试打印。整个函数体被 if (0) 关掉了 —— 保留它是为了方便临时打开排查，
// 但注意 print() 在 depthFirstSolve 里被**无条件调用**（每个解都调），
// 一旦把 0 改成 1，输出量会随排列数爆炸。
void LRSolution::print()
{
    if (0)
    {
        cout << "Current: ";
        printf("(%d) ", uninsertedCells.size());
        for (Module *curModule : insertedCells)
        {
            if (curModule)
            {
                cout << curModule->name << " ";
            }
            else
            {
                cout << "OOO" << " ";
            }
        }
        cout << endl;

        cout << "Location: ";
        for (double curModuleLocation : xLocations)
        {
            cout << curModuleLocation << " ";
        }

        cout << endl;

        cout << "Cost: " << solutionCost << endl
             << endl;
    }
}

// ============================================================================
// initializeNetModuleCount —— 统计窗口内每条 net 上「还没排定」的单元数
// ============================================================================
// 这是增量式代价计算的基础设施。
//   每排定一个单元，就把它所在每条 net 的计数减一；
//   减到 0 说明这条 net 上属于本窗口的单元**全部**排定了，
//   它的 bounding box 才最终确定，此时才结算它的 HPWL。
// 这样得到的 solutionCost 沿搜索路径单调不减，可作分支限界的下界。
void LRSolution::initializeNetModuleCount()
{
    for (Module *uninsertedModule : uninsertedCells)
    {
        if (uninsertedModule == NULL)
        {
            continue;
        }

        for (Net *curNet : uninsertedModule->nets)
        {
            // if (netModuleCount.find(curNet) == netModuleCount.end()) // key not found
            // {
            //     //!!!!!! problem: find() fails when using pointer as key
            //     netModuleCount[curNet] = 0;
            // }
            // else
            // {
            // cout<<"ffff\n";
            // 用 net id（int）而非 Net* 作 key，避开指针作 key 时的比较不确定性
            // （见上面被注释掉的「find() fails when using pointer as key」）。
            // map 的 operator[] 会在 key 不存在时自动插入 0，所以「+1」即完成初始化。
            netModuleCount[curNet->idx] = netModuleCount[curNet->idx] + 1; // netModuleCount[curNet] will be inserted as 0 if not find
            // }
        }
    }

    // for (auto iteModule = uninsertedCells.begin(); iteModule != uninsertedCells.end(); iteModule++)
    // {
    //     if (*iteModule == NULL)
    //     {
    //         continue;
    //     }

    //     for (auto iteNet = (*iteModule)->nets.begin(); iteNet != (*iteModule)->nets.end(); iteNet++)
    //     {
    //         netModuleCount[*iteNet] = netModuleCount[*iteNet] + 1;
    //     }
    // }
}

// ============================================================================
// recalculateCost —— 完整解的重算：把窗口相关 net 的 HPWL 全部加总
// ============================================================================
// 与搜索中的增量计算不同，这里是「从零重算」，用在构造初始上界
// （originalSolution，即当前实际布局）时。
// 用 set 去重很关键：多个窗口内单元可能同属一条 net，
// 若不去重就会把同一条 net 的 HPWL 重复计入。
void LRSolution::recalculateCost()
{
    if (isComplete())
    {
        std::set<Net *> netSet; //! to avoid repeated calculation of net HPWL

        // Collet all involved net id
        for (Module *insertedModule : insertedCells)
        {
            // Skip Whitespace
            if (insertedModule == NULL)
            {
                continue;
            }

            for (Net *curNet : insertedModule->nets) // module->nets initialized in bookshelf parser
            {
                netSet.insert(curNet);
            }
        }

        // Calculate the involved wirelength
        // 注意 calcNetHPWL 读的是 placedb 里的**当前坐标**，
        // 所以调用它之前必须保证单元已经被摆到这个解对应的位置上。
        solutionCost = 0;
        for (Net *curNet : netSet)
        {
            solutionCost += curNet->calcNetHPWL();
        }
    }
    //! 疑似问题：若解不完整（isComplete() 为 false），这里**什么都不做**，
    //! solutionCost 保持旧值（构造时为 0），调用方无法区分「代价真的是 0」
    //! 和「没算」。目前只有 originalSolution 会调用它，恰好总是完整解，尚未暴露问题。
}

// ============================================================================
// createSuccessorSolution —— 生成「把当前未排单元放到 currentX」这一后继解
// ============================================================================
// 这是搜索树的「边」：一个含 k 个未排单元的部分解有 k 个后继，
// 分别对应「下一个放谁」，递归下去就枚举出了窗口内单元的**全部排列**。
// 后继解要复制父解的状态（已排列表、x 记录、net 计数、代价），
// 再把自己这一步的增量叠加上去。
LRSolution *LRSolutionIterator::createSuccessorSolution()
{
    LRSolution *succesorSolution = new LRSolution(pointedSolution->placedb);
    succesorSolution->whiteSpaceWidth = pointedSolution->whiteSpaceWidth;
    succesorSolution->xLocations = pointedSolution->xLocations;
    // succesorSolution->m_pDB = pointedSolution.m_pDB;

    // Push moduleIterator into list
    Module *curModule = *moduleIterator;
    succesorSolution->insertedCells = pointedSolution->insertedCells;
    succesorSolution->insertedCells.push_back(curModule);

    // Copy uninsertedCells except for moduleIterator
    for (Module *copyModule : pointedSolution->uninsertedCells)
    {
        if (copyModule != curModule)
        {
            succesorSolution->uninsertedCells.push_back(copyModule);
        }
    }

    // Move *moduleIterator module to new location
    // 注意这里**真实地**移动了 placedb 中的单元（不是只记在解里）：
    // 因为下一步要用 calcNetHPWL 算线长，而它读的就是 placedb 的当前坐标。
    //! 疑似问题：这个副作用在分支被剪枝 / successor 被 delete 时**不会回滚**，
    //! 因此搜索过程中 placedb 的坐标是「某条搜索路径的中间态」，
    //! 会污染后续 HPWL 的计算（详见 solveForBestOrder 处的说明）。
    if (curModule != NULL)
    {
        POS_2D curPos = curModule->getLL_2D();
        pointedSolution->placedb->setModuleLocation_2D(curModule, pointedSolution->currentX, curPos.y);
        succesorSolution->xLocations.push_back(pointedSolution->currentX);

        // Calculte the new currentX
        // 紧凑码放：下一个单元的左端紧贴本单元右端
        succesorSolution->currentX = pointedSolution->currentX + curModule->getWidth();
    }
    else
    {
        // NULL 代表空白：它占据 whiteSpaceWidth 的宽度，但不需要移动任何实体
        succesorSolution->xLocations.push_back(pointedSolution->currentX);

        // Calculte the new currentX
        succesorSolution->currentX = pointedSolution->currentX + pointedSolution->whiteSpaceWidth;
    }

    // Calculate the new m_bound
    succesorSolution->solutionCost = pointedSolution->solutionCost;
    succesorSolution->netModuleCount = pointedSolution->netModuleCount;

    // 增量结算：把刚排定的单元所在 net 的计数减一，
    // 减到 0 就说明这条 net 在窗口内的所有单元都已排定，可以结算它的 HPWL 了。
    if (curModule != NULL)
    {
        for (Net *curNet : curModule->nets)
        {
            succesorSolution->netModuleCount[curNet->idx]--;

            // No module in this net is unplaced
            // Calculate the net length
            assert(succesorSolution->netModuleCount[curNet->idx] >= 0);
            if (succesorSolution->netModuleCount[curNet->idx] == 0)
            {
                // cout << "added\n";
                // 只在「最后一个窗口内单元排定」时结算一次 -> 代价单调不减，可作下界
                succesorSolution->solutionCost += curNet->calcNetHPWL();
            }
            else
            {
                // cout << "value is: " << succesorSolution->netModuleCount[curNet] << endl;
            }
        }
    }
    return succesorSolution;
}

// ============================================================================
// GlobalSwapDP::solve —— 全局交换主循环
// ============================================================================
// 论文：An efficient and effective detailed placement algorithm
// 与前两个算子不同，这里**没有几何滑动窗口**，搜索范围是每个单元各自的
// optimal region —— 由它相连的各条 net 的 bounding box 交集决定
// （placedb->getOptimialRegion）。直观上就是「这个单元挪到哪儿线长可能更短」
// 的那一小片区域；理论依据是：把单元移出其 net bounding box 之外，
// 只会让相关 net 的包围盒变大，线长不可能更优。
//
// 对每个单元依次做四件事：
//   1. 求 optimal region，裁剪到 coreRegion 内，并对齐到 site step / 行高
//      （对齐是为了让候选位置天然满足对齐约束，避免产生非法位置）；
//   2. 在 region 内收集候选：空白 slot（切出等宽位置）+ 等宽单元（可安全对调）；
//   3. 先把单元从行结构上摘下，再依次试候选：
//      (a) 插进空白 slot —— 只要线长变小就接受；
//      (b) 与等宽单元对调 —— 只要两者线长之和变小就接受。
//      两者都是「找第一个有收益的」而非最优的（论文明确这么做以省时间）；
//   4. 都不行就放回原位，恢复行结构。
//
// 合法性保证：只在空白里插，或只与**等宽**单元对调 —— 两种情况都不会造成重叠；
// 且所有候选坐标都来自 rowSpaces 或已有单元的左下角，本身已对齐 site / row。
void GlobalSwapDP::solve()
{
    //! this function follows the algorithm introduced in the paper: An efficient and effective detailed placement algorithm
    int successCount = 0;
    for (Module *curCell : placeDB->dbNodes)
    {
        // macro 不参与：它尺寸大、通常已由专门的 macro legalization 处理，
        // 而且它跨多行，简单的等宽交换/空白插入都不适用
        if (curCell->isMacro)
        {
            continue;
        }
        float curCellWidth = curCell->getWidth();
        POS_2D curCellPos = curCell->getLL_2D();

        vector<POS_2D> candidateEmptySlots;
        vector<POS_2D> candidateModuleSlots;

        // 1. determine search region
        // optimal region = 该单元所有相连 net 的 bounding box 的交集
        CRect searchRegion = placeDB->getOptimialRegion(curCell);
        // possible extension: intersection of optimalRegion and other wanted regions (for example, regions with lower overflow or congestion).

        // boundary check
        // region 与 coreRegion 完全不相交 -> 没有任何合法位置可试，跳过。
        // 用 float_* 系列比较函数是为了规避浮点误差导致的边界抖动。
        if (float_greaterorequal(searchRegion.ll.y, placeDB->coreRegion.ur.y) || float_lessorequal(searchRegion.ur.y, placeDB->coreRegion.ll.y) || float_greaterorequal(searchRegion.ll.x, placeDB->coreRegion.ur.x) || float_lessorequal(searchRegion.ur.x, placeDB->coreRegion.ll.x))
        {
            continue;
        }

        searchRegion.ll.x = max(searchRegion.ll.x, placeDB->coreRegion.ll.x);
        searchRegion.ll.y = max(searchRegion.ll.y, placeDB->coreRegion.ll.y);
        searchRegion.ur.x = min(searchRegion.ur.x, placeDB->coreRegion.ur.x);
        searchRegion.ur.y = min(searchRegion.ur.y, placeDB->coreRegion.ur.y);

        // align search region to avoid bugs
        // 对齐的意义：后面所有候选 slot 都是从 region 里按 site 步长切出来的，
        // 如果 region 边界本身没对齐 site / 行，切出来的位置就可能落在非法的
        // 非 site 对齐点上，导致 insertModule 失败或产生未对齐的布局。
        // 1.1 align to sites
        // 左边界向下取整（floor）、右边界向上取整（ceil）—— 向外扩张而非收缩，
        // 保证不会因为取整而把整个 region 抹掉。
        double siteStep = placeDB->dbSiteRows.front().step;
        //? when use ceil for left and floor for right, all small search regions are ignored.
        //? How about ceil for right and floor for left, so we always have a search region?
        double newLeft = floor((searchRegion.ll.x - placeDB->coreRegion.ll.x) / siteStep) * siteStep + placeDB->coreRegion.ll.x;
        double newRight = ceil((searchRegion.ur.x - placeDB->coreRegion.ll.x) / siteStep) * siteStep + placeDB->coreRegion.ll.x;
        searchRegion.ll.x = newLeft;
        searchRegion.ur.x = newRight;

        // 1.2 align to placement rows
        double newBottom = floor((searchRegion.ll.y - placeDB->coreRegion.ll.y) / placeDB->commonRowHeight) * placeDB->commonRowHeight + placeDB->coreRegion.ll.y;
        double newTop = ceil((searchRegion.ur.y - placeDB->coreRegion.ll.y) / placeDB->commonRowHeight) * placeDB->commonRowHeight + placeDB->coreRegion.ll.y;
        searchRegion.ll.y = newBottom;
        searchRegion.ur.y = newTop;
        // watch out the search region is a point or a line
        // 退化成一个点或一条线 -> 宽度或高度为 0，放不下任何东西
        if (float_equal(searchRegion.ll.x, searchRegion.ur.x) || float_equal(searchRegion.ll.y, searchRegion.ur.y))
        {
            continue;
        }

        //! watch out when the cell is already in the region
        // 单元已经在自己的最优区域里了 -> 理论上线长已经接近最优，
        // 再挪动很可能只是白白折腾，直接跳过（也省时间）。
        //! 疑似问题：这个判断把「已在 optimal region 内」的单元全部排除，
        //! 但它们未必处在 region 内的**最优具体位置**（region 只保证 bounding box
        //! 层面的最优性），因此这里会漏掉一部分可优化单元。
        if (searchRegion.inside(curCellPos))
        {
            continue;
        }

        // 2. find empty slots and exchangable cells in the search region
        //  !similar to find slots and cells in ISM in implementation, need structures like ISMrows to manage spaces.
        //  currently we use a simplified version of the implementation introduced in the paper, no penalty, only exchange with cells with same size or empty slots
        //  'In the actual implementation, to saveruntime, for aselected cell,
        //  we do not pick the cell with the best benefit in its optimal region.
        //  Instead, we pick the first good cell that can give us certain benefit.'
        // 上面的取舍很关键：全局交换是 O(单元数 x 候选数) 的主循环，
        // 若对每个单元都要遍历全部候选挑最优，开销会大很多；
        // 而实验表明「取第一个有收益的」与「取最优的」质量差别不大。
        // 另外这里只处理「等宽交换 + 空白插入」，是论文完整算法的简化版
        // （原论文还有带惩罚项的更一般交换）。
        // cout << "ll: " << searchRegion.ll.y << endl;
        int startRowIndex = max(placeDB->y2RowIndex(searchRegion.ll.y), 0); // lowest row
        double endRowBottom = min(double(searchRegion.ur.y), placeDB->coreRegion.ur.y - placeDB->commonRowHeight);
        // cout << "ur: " << searchRegion.ur.y << endl;
        int endRowIndex = placeDB->y2RowIndex(endRowBottom); // highest row
        assert(endRowIndex >= startRowIndex);

        // 2.1 find empty slots
        // 从 region 内每一行的空闲区间里，按「本单元宽度」切出若干个 slot。
        // 这里用 curCellWidth 而不是 ISM 里的 maxWIDTH：因为只考虑
        // 「curCell 独占某个位置」这一种情形，不需要和其它单元等宽。
        //! 疑似问题：行范围是 i < endRowIndex（不含），而下面找可交换单元时
        //! 用的是 i <= endRowIndex（含），两者不一致，最高那一行的空白被漏掉。
        //! 另外 emptySlotCount 超过 maxSlotCount 时只 break 了最内层循环，
        //! 外层两个 for 仍会继续 push_back，限流并不严格。
        int emptySlotCount = 0;
        for (int i = startRowIndex; i < endRowIndex; i++)
        {
            for (auto curRowSpaceIter = GSRows[i].rowSpaces.lower_bound(searchRegion.ll.x);
                 curRowSpaceIter != GSRows[i].rowSpaces.lower_bound(searchRegion.ur.x);
                 curRowSpaceIter++)
            {
                int slotCount = (int)(curRowSpaceIter->second / curCellWidth);
                for (int j = 0; j < slotCount; j++)
                {
                    emptySlotCount++;
                    if (emptySlotCount > maxSlotCount)
                    {
                        break;
                    }
                    POS_2D slot;
                    slot.x = curRowSpaceIter->first + j * curCellWidth;
                    slot.y = GSRows[i].ll.y;

                    assert(placeDB->coreRegion.ll.y + i * placeDB->commonRowHeight == GSRows[i].ll.y);

                    candidateEmptySlots.push_back(slot);
                }
            }
        }

        // 2.2 find cells for swap
        // 只有在空白 slot 还不够多时才继续找可交换的单元，
        // 优先用空白（纯插入，不改变别人），其次才考虑交换（会动到另一个单元）。
        multimap<double, Module *> moduleMap; // use multimap for auto sorting (by module width)
        vector<Module *> candidateModules;    // correspond with the vector candidateModuleSlots
        //! 疑似问题：candidateModules 里可能包含 curCell 自己 ——
        //! 虽然上面已跳过「curCell 在 searchRegion 内」的情形，
        //! 但 searchRegion.inside() 用的是中心点判断，边界情形未必能排除。
        //! 若把 curCell 当成自己的交换对象，下面会先 removeModule 再 insertModule，
        //! 逻辑上是「原地不动」，不会致命，但会白算一次并被误计入 successCount。
        if (emptySlotCount < maxSlotCount)
        {
            // build module size map, with all modules in the window
            for (int i = startRowIndex; i <= endRowIndex; i++) // i< endRowIndex in ntuplace
            {
                auto startModule = GSRows[i].rowModules.lower_bound(searchRegion.ll.x); //! rowModules initialized in initializeISMRows()
                auto endModule = GSRows[i].rowModules.lower_bound(searchRegion.ur.x);

                for (auto curModuleIter = startModule; curModuleIter != endModule; curModuleIter++)
                {
                    // if (curModuleIter->second->getHeight() == placeDB->commonRowHeight)
                    if (!(curModuleIter->second->isMacro)) // macros are ignored
                    {
                        //! adding cells in the window with same sizes
                        // 只收**等宽**单元：两个同宽单元对调位置后，
                        // 各自占的宽度不变，因此绝不会与邻居重叠 —— 这是合法性的来源。
                        float curModuleWidth = curModuleIter->second->getWidth();
                        if (float_equal(curCellWidth, curModuleWidth))
                        {
                            moduleMap.insert(pair<double, Module *>(curModuleWidth, curModuleIter->second));
                        }
                    }
                }
            }

            //! moduleList stores ALL modules in current window(windows), sorted by module width, decreaingly (sorted increasingly in moduleMap)
            int swapModuleCount = 0;

            for (auto rIter = moduleMap.rbegin(); rIter != moduleMap.rend(); rIter++)
            {
                swapModuleCount++;
                if (emptySlotCount + swapModuleCount > maxSlotCount)
                {
                    break;
                }
                candidateModules.push_back(rIter->second);
                candidateModuleSlots.push_back(rIter->second->getLL_2D());
            }
        }

        // 3. exchange with empty slots if there are any, or exchange with a cell that give us certain benefit
        bool insertedFlag = false; // true if we find a candidate slot for current cell
        // cout << "CP2\n";
        // 先把 curCell 从它原来所在的行上摘下来：
        //   这样它腾出的空间会还给 rowSpaces，后面 insertModule 才有地方可放；
        //   同时 oldWL 必须在移动之前算好，作为比较的基准。
        int curCellRowID = placeDB->y2RowIndex(curCellPos.y); //! curCellPos is the original cell position before inserting
        GSRows[curCellRowID].removeModule(curCell);           // remove from curCellPos, not the potential new location
        double oldWL = placeDB->calcModuleHPWLfast(curCell);  // looking for the first one that gives us benefits, not the best one

        // 3.1 try to insert cell into the empty
        // 逐个试空白 slot：把单元真的挪过去再算线长，只要**变小**就立刻接受（first-fit）。
        for (POS_2D curSlot : candidateEmptySlots)
        {
            placeDB->setModuleLocation_2D(curCell, curSlot.x, curSlot.y);
            if (oldWL > placeDB->calcModuleHPWLfast(curCell))
            {
                insertedFlag = true;
                // cout << "CP3\n";
                int curSlotRowID = placeDB->y2RowIndex(curSlot.y);
                // update GSRows, only need insertModule and removeModule? remove tail actually means removespace, however it seems we only need to remove space of tails, so its called remove tail?
                GSRows[curSlotRowID].insertModule(curSlot.x, curCell);
                break; // looking for the first one that gives us benefits, not the best one
            }
        }

        if (!insertedFlag)
        {
            // 3.2 try to swap cells
            // 空白都不合适，就试着与某个等宽单元对调。
            // 判据是**两个单元线长之和**变小（而不是只看 curCell），
            // 否则可能「牺牲别人成全自己」，全局线长反而变差。
            //! 疑似问题：候选是按宽度降序（moduleMap 的 rbegin 顺序）取的，
            //! 宽度都相等时这个顺序等价于 multimap 内部的插入顺序，
            //! 与收益大小完全无关 —— 因此确实是「第一个有收益的」而非「最优的」。
            for (int i = 0; i < candidateModules.size(); i++)
            {
                double oldWL2 = oldWL + placeDB->calcModuleHPWLfast(candidateModules[i]); // save old wirelength
                // cout << "CP4\n";
                int curSlotRowID = placeDB->y2RowIndex(candidateModuleSlots[i].y);
                GSRows[curSlotRowID].removeModule(candidateModules[i]); // remove from original location
                // swap：curCell 去对方位置，对方来 curCell 的原位置（等宽 -> 一定放得下）
                placeDB->setModuleLocation_2D(curCell, candidateModuleSlots[i].x, candidateModuleSlots[i].y);
                placeDB->setModuleLocation_2D(candidateModules[i], curCellPos.x, curCellPos.y);

                if (oldWL2 > placeDB->calcModuleHPWLfast(curCell) + placeDB->calcModuleHPWLfast(candidateModules[i]))
                {
                    insertedFlag = true;
                    GSRows[curSlotRowID].insertModule(candidateModuleSlots[i].x, curCell);
                    GSRows[curCellRowID].insertModule(curCellPos.x, candidateModules[i]);
                    break;
                }
                else
                {
                    GSRows[curSlotRowID].insertModule(candidateModuleSlots[i].x, candidateModules[i]);
                    placeDB->setModuleLocation_2D(candidateModules[i], candidateModuleSlots[i].x, candidateModuleSlots[i].y); //! put it back!
                }
            }
        }

        if (!insertedFlag)
        {
            // no candidate is selected
            // 回滚：把单元放回原位，并把行结构恢复原样。
            // 这一步很重要 —— 前面已经 removeModule 腾出了空间，
            // 若不插回去，rowSpaces 就会凭空多出一块空白，账目与实际不一致。
            placeDB->setModuleLocation_2D(curCell, curCellPos.x, curCellPos.y);
            GSRows[curCellRowID].insertModule(curCellPos.x, curCell);
        }
        else
        {
            successCount++;
        }
    }
    cout << "Found position for " << successCount << " among all " << placeDB->dbNodes.size() - placeDB->dbMacroCount << " cells\n";

    // todo: when does global swap stop? Just iterate through all cells? record current HPWL decrease and decide stop criteria
}

void GlobalSwapDP::initialization()
{
    initializeParams();
    initializeGSRows();
}

// maxSlotCount 是全局交换唯一的性能旋钮：
// 它同时限制「空白候选数」和「空白 + 可交换单元总数」。
// 调大 -> 每个单元看的位置更多、质量更好，但主循环开销线性上升。
void GlobalSwapDP::initializeParams()
{
    maxSlotCount = 30; // affect runtime
}

// ============================================================================
// initializeGSRows —— 与 ISMDP::initializeISMRows() 完全同构
// ============================================================================
// 先用 SiteRow::intervals 铺好每行的空闲区间，再把已有单元（含 macro）插进去。
// 顺序不能颠倒：insertModule 的语义就是「从某个空闲区间里切走一段」。
void GlobalSwapDP::initializeGSRows()
{
    //! same as ISMDP::initializeISMRows()
    double ISMRowLength = placeDB->coreRegion.getWidth();
    // cout << "initialize ISM rows: \n";

    // 1.create ISMRows
    for (SiteRow curRow : placeDB->dbSiteRows)
    {
        ISMRow newRow(placeDB->coreRegion.ll.x, curRow.bottom, ISMRowLength);
        for (Interval curInterval : curRow.intervals)
        {
            // cout<<"father: "<<curInterval.start<<" "<<curInterval.getLength()<<endl;
            newRow.rowSpaces[curInterval.start] = curInterval.getLength();
            // cout<<"son: "<<newRow.rowSpaces[curInterval.start]<<endl;
        }
        GSRows.push_back(newRow);
    }

    // 2. insert cells to ISMRows
    for (Module *curModule : placeDB->dbNodes) // all nodes, macros and cells
    {

        int coveredRowCount = 1;
        POS_2D curModulePos = curModule->getLL_2D();
        float moduleHeight = curModule->getHeight();
        if (moduleHeight > placeDB->commonRowHeight) //?? macro included???
        {
            coveredRowCount = (int)(moduleHeight / placeDB->commonRowHeight);
            if (placeDB->commonRowHeight * coveredRowCount < moduleHeight)
            {
                coveredRowCount++; // for macros
            }
        }

        int pos = placeDB->y2RowIndex(curModulePos.y);

        //! 疑似问题：与 ISMDP::initializeISMRows() 一样，pos + i 没有做越界检查。
        //! 此外这里的整段逻辑与 ISMDP 版本几乎逐行重复，属于可抽取的公共代码。
        for (int i = 0; i < coveredRowCount; i++)
        {
            GSRows[pos + i].insertModule(curModulePos.x, curModule);
        }
    }
}
