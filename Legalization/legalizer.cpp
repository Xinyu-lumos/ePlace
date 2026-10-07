// ============================================================================
// 模块总览：Legalization 的实现（Abacus 标准单元合法化 + SA macro 合法化）
// ============================================================================
// 【文件职责】
//   把 EPlace 全局布局输出的「重叠、未对齐」结果，变成「无重叠、对齐 site、
//   位移尽量小」的合法布局。本文件实现两个合法化器：
//     · AbacusLegalizer    —— 标准单元（std cell）合法化，主力，严格遵循 Abacus 论文
//     · SAMacroLegalizer   —— macro 合法化（ePlace-MS 的 mLG 阶段），模拟退火
//
// 【为什么需要两段式（先 macro 后 std cell）】
//   macro 尺寸大、不能塞进 site row，而且会**占据并切断** site row。
//   只有先让 macro 互不重叠，PlaceDB::removeBlockedSite() 才能算出
//   每条 site row 真正可用的 intervals；std cell 才有地方可放。
//   所以正确顺序是：SAMacroLegalizer → removeBlockedSite → AbacusLegalizer。
//
// 【Abacus 的算法脉络】（AbacusLegalizer，本文件前半部分）
//   initialization()
//     ├─ initializeCells()   收集 std cell 并按全局 x 排序（保证行内插入顺序即 x 顺序）
//     └─ initializeSubrows() 由 site row 的 intervals 展开成 AbacusRow(subrow)，按 bottom 升序排
//   legalization()           逐个 cell 做行分配
//     ├─ 找「最近行」closestRowIndex
//     ├─ 向上 / 向下双向搜索候选行，用 |Δy| 下界剪枝
//     ├─ placeRow(ABACUS_TRIAL)  在 subrow 拷贝上试算代价，不污染真实状态
//     └─ placeRow(ABACUS_FINAL)  对最优行真正写入
//   placeRow() → AbacusRow::addCell / addCluster / collapse
//     行内合法化的核心是 cluster：cell 若与前驱 cluster 重叠就合并成刚体，
//     再按 q/e 求该刚体的最优 x 并对齐 site step —— 这就是 Abacus 优于 Tetris 的地方。
//
// 【mLG 的算法脉络】（SAMacroLegalizer，本文件后半部分）
//   代价 = HPWL + miuD * (被 macro 压住的 cell 面积) + miuO * (macro 间重叠面积)
//   外层 j 循环：重设 SA 温度/扰动半径，并按 beta 放大 miuO，逼迫重叠消失
//   内层 k 循环：温度与半径逐步衰减，SAperturb() 随机扰动一个 macro，
//               acceptPerturb() 按 Metropolis 准则决定接受或回滚
//   终止：overlapFree（macro 总面积 − macro 并面积 == 0）
//
// 【核心数据结构】见 legalizer.h：AbacusRow / AbacusCellCluster / SAMacroLegalizationBin_2D
// 【主要入口函数】
//   AbacusLegalizer::legalization()   —— std cell 合法化
//   SAMacroLegalizer::legalization()  —— macro 合法化
// ============================================================================

#include "legalizer.h"

// ============================================================================
// AbacusLegalizer::legalization —— 标准单元合法化的主入口
//
// 整体是「逐个 cell 做行分配 + 行内 Abacus 聚类」的贪心流程：
//   dbCells 已按全局 x 升序排好，因此**处理顺序即 x 顺序**，
//   这样同一行（以及同一 cluster）内插入的 cell 天然按 x 非递减排列，
//   Abacus 增量维护 q 的语义才成立。
//
// 为什么按 x 顺序处理：Abacus 的行内代价是在「已有 cluster 序列」上增量计算的，
// 若乱序插入，cluster 内的 cell 顺序就不再代表几何顺序，
// 算出来的最优 x 会把 cell 排成互相穿插的样子，位移被严重低估。
//
// 单 cell 的三步：找最近行 → 双向搜索候选行（带剪枝）→ 对最优行落地。
// ============================================================================
void AbacusLegalizer::legalization()
{
    initialization();
    int cellCount = dbCells.size();
    // cout << cellCount << endl;
    for (Module *curCell : dbCells)
    {
        // cout<<curCell->name<<endl;
        double cost = DOUBLE_MAX; // current minimum cost
        //! assume subrows are sorted by y coordinate(non decreasing)
        // 1. find closest row
        int subrowCount = subrows.size();
        // cout<<subrowCount<<endl;
        int bestRowIndex = -1;
        int closestRowIndex = 0;
        // 找到「最近行」：cell 的左下角 y 落在某条 row 的半高范围内即认为属于该行。
        // 用 0.5*height 作容差，等价于把 cell 的 y 吸附到最近的一条 row。
        //! 疑似问题：这里只比较 curRow.bottom 与 cell 的 **LL.y**（下边界），
        //! 没有减去/加上 cell 自身高度的一半，多倍高（multi-row-height）cell 的定位会偏。
        for (AbacusRow curRow : subrows)
        {
            if (float_lessorequal(fabs(curRow.bottom - curCell->getLL_2D().y), 0.5 * curRow.height)) //? double check here
            {
                break;
            }
            closestRowIndex++;
        }
        //! 疑似问题：若没有任何 row 满足条件，closestRowIndex 会停在 subrowCount（越界值）。
        //! 此时向上的循环不执行、向下的循环从 subrowCount-1 开始，侥幸不会越界，
        //! 但语义上应显式钳到 subrowCount-1。
        // cout <<closestRowIndex<<endl;
        //! 2. search in two directions(higher rows and lower rows) for the best row
        //? optimize here?
        // 双向搜索 + 剪枝的关键：
        //   真实代价是含 x 位移的欧氏距离 dist = sqrt(dx^2 + dy^2) >= |dy| = costLowerBound。
        //   由于 subrows 已按 bottom 升序排好，从最近行往外走 |dy| 单调增大，
        //   所以一旦 costLowerBound >= 当前最优 cost，更远的行必然更差，可以安全 break。
        //   这让「全行扫描」退化成「只看邻近若干行」，是合法化速度的关键。
        // first search in higher rows
        for (int i = closestRowIndex; i < subrowCount; i++)
        {
            // cout<<"higher "<<i<<endl;
            // 容量预筛：用「已占用宽度 width」粗略判断该行是否还塞得下，
            // 避免对明显放不下的行做昂贵的 placeRow 试算
            if (float_greater(subrows[i].width + curCell->getWidth(), subrows[i].end.x - subrows[i].start.x))
            {
                continue;
            }
            double costLowerBound = fabs(curCell->getLL_2D().y - subrows[i].bottom); // cost for moving cell in y direction only, //?quadratic or not?
            // determine cost to move a cell to this row!
            if (float_greaterorequal(costLowerBound, cost))
            {
                break;
            }
            // 试算：placeRow 内部会在 subrow 的拷贝上走一遍完整的 Abacus 插入 + collapse，
            // 因此不会污染真实状态，可以放心对多个候选行各试一次
            double newCost = placeRow(curCell, i, ABACUS_TRIAL);
            if (float_lessorequal(newCost, cost))
            {
                bestRowIndex = i;
                cost = newCost;
            }
        }
        // then search in lower rows
        for (int i = closestRowIndex - 1; i >= 0; i--)
        {
            // cout<<"lower "<<i<<endl;
            //! 疑似问题：i==0 时直接 break，导致 subrows[0] 在这一轮「向下搜索」中
            //! 永远不会被尝试（正确写法应是循环条件 i>=0，在 i<0 时自然结束）。
            //! 当 closestRowIndex>0 且第 0 行恰是最优行时，会选出次优行、白白增大位移。
            if (i == 0)
            {
                break;
            }
            if (float_greater(subrows[i].width + curCell->getWidth(), subrows[i].end.x - subrows[i].start.x))
            {
                continue;
            }
            double costLowerBound = fabs(curCell->getLL_2D().y - subrows[i].bottom); // cost for moving cell in y direction only
            if (float_greaterorequal(costLowerBound, cost))
            {
                break;
            }
            double newCost = placeRow(curCell, i, ABACUS_TRIAL);
            if (float_lessorequal(newCost, cost))
            {
                bestRowIndex = i;
                cost = newCost;
            }
        }
        //! 疑似问题：若所有 subrow 都因容量预筛被 continue 掉，bestRowIndex 会保持 -1。
        //! Debug 下由 assert 拦住，Release 下会以 subrows[-1] 越界访问。
        assert(bestRowIndex != -1);
        // 3. place cell to the best row
        // 只有这一步真正修改状态：累加行占用宽度、把 cell 写进 cluster 并 collapse
        placeRow(curCell, bestRowIndex, ABACUS_FINAL);
    }

    //! 4. write cell locations! In placeRow we only determine cluster's best position, without writing cell location

    // 回写阶段：cluster 只记录了整体左端 x，这里按「cluster 内 cell 的宽度顺序」
    // 依次把 x 展开给每个 cell，y 统一取所在 subrow 的 bottom（即对齐到 site row）。
    for (AbacusRow curSubrow : subrows)
    {
        for (AbacusCellCluster curCluster : curSubrow.clusters)
        {
            int x = curCluster.x;
            //! 疑似问题：x 声明为 int，而 curCluster.x 是 float、curCell->getWidth() 也是 float。
            //! 每次赋值/累加都会截断取整，误差沿 cluster 内 cell 逐个累积，
            //! 可能让 cell 偏移出 site 网格；末尾的 assert 也只在 Debug 下生效。
            for (Module *curCell : curCluster.cells)
            {
                placeDB->setModuleLocation_2D(curCell, x, curSubrow.bottom);
                x += curCell->getWidth();
                assert(x <= curSubrow.end.x);
            }
        }
    }
}

// ----------------------------------------------------------------------------
// AbacusLegalizer::initialization —— 合法化的准备工作
// 顺序不可交换：必须先有 cell 列表，再由（已合法化的）macro/terminal 算出可用 subrow。
// initializeObstacles() 已被注释掉——它的职责搬到了 PlaceDB::removeBlockedSite()，
// 因为「哪些 site 被挡住」是全局几何信息，放在 placedb 里才能被 legalizer
// 和 detailed placer 共用（详见 initializeSubrows 的注释）。
// ----------------------------------------------------------------------------
void AbacusLegalizer::initialization()
{
    // segmentFaultCP("cp2");
    initializeCells();
    // segmentFaultCP("cp3");
    // initializeObstacles();
    // segmentFaultCP("cp4");
    // initializeSubrows();
    // cout<<"\nusing new leal!\n";
    initializeSubrows();
}

// ----------------------------------------------------------------------------
// AbacusLegalizer::initializeCells —— 收集所有标准单元并按 x 升序排序
// 只收 !isMacro 的 std cell：macro 由 SAMacroLegalizer 单独处理，
// 它们在 Abacus 这里的角色是「挡住 site 的障碍物」，而不是要被摆放的对象。
// 排序是后续一切的前提（见 legalization 的注释）：处理顺序即 x 顺序。
// ----------------------------------------------------------------------------
void AbacusLegalizer::initializeCells()
{
    for (Module *curNode : placeDB->dbNodes)
    {
        assert(curNode);
        if (!curNode->isMacro)
        {
            dbCells.push_back(curNode);
        }
    }
    // segmentFaultCP("cp5");
    // float_less 是带 EPS 容差的比较，避免浮点抖动导致排序不稳定
    sort(dbCells.begin(), dbCells.end(), [=](Module *a, Module *b)
         { return float_less(a->getLL_2D().x, b->getLL_2D().x); });
    // segmentFaultCP("cp6");
    //? what's the potential bug when one cell is completely inside another? Should think about it.
}

// ----------------------------------------------------------------------------
// AbacusLegalizer::initializeObstacles —— 收集 macro / terminal 作为障碍物（已废弃）
// 原职责：把 coreRegion 内的 terminal 和所有 macro 收成 Obstacle 矩形并按 x 排序，
// 供 initializeSubrowsOld() 把完整 site row 切成可用 subrow。
//! 疑似问题：该函数已于 2024.08.26 被 PlaceDB::removeBlockedSite() 取代
//! （见 legalizer.h 的声明注释），initialization() 中也不再调用它，
//! 是一段仍然编译进产物的死代码；同时 AbacusLegalizer::obstacles 成员随之成为死变量。
// ----------------------------------------------------------------------------
void AbacusLegalizer::initializeObstacles()
{
    //! ignore terminals that are outside the core region
    //! overlap between macros should be eliminated first!
    for (Module *curTerminal : placeDB->dbTerminals)
    {
        if (float_greater(placeDB->coreRegion.ur.y, curTerminal->getLL_2D().y) && float_less(placeDB->coreRegion.ll.y, curTerminal->getUR_2D().y)) // overlap in y direction
        {
            if (float_greater(placeDB->coreRegion.ur.x, curTerminal->getLL_2D().x) && float_less(placeDB->coreRegion.ll.x, curTerminal->getUR_2D().x)) // overlap in x direction
            {
                Obstacle newObstacle;
                newObstacle.ll = curTerminal->getLL_2D();
                newObstacle.ur = curTerminal->getUR_2D();
                obstacles.push_back(newObstacle);
            }
        }
    }

    for (Module *curNode : placeDB->dbNodes)
    {
        if (curNode->isMacro)
        {
            Obstacle newObstacle;
            newObstacle.ll = curNode->getLL_2D();
            newObstacle.ur = curNode->getUR_2D();
            obstacles.push_back(newObstacle);
        }
    }
    //! sort obstacles by x coordinate first
    sort(obstacles.begin(), obstacles.end(), [=](Obstacle a, Obstacle b)
         {
            if (!float_equal(a.ll.x , b.ll.x))
            {
                return float_less(a.ll.x ,b.ll.x);
                
            }
            else if (!float_equal(a.ll.y , b.ll.y))
            {
                return float_less(a.ll.y ,b.ll.y);
            }
            else
            {
                // terminals/macros overlap!!
                cerr<<"TERMINALS/MACROS OVERLAP!\n";
                exit(0);
            } });
}

// ----------------------------------------------------------------------------
// AbacusLegalizer::initializeSubrowsOld —— 由障碍物切分 subrow（旧实现，已废弃）
// 思路：对每条完整 site row，按 x 从左到右扫过所有障碍物，
// 把「上一个障碍物的右边界 → 当前障碍物的左边界」之间的空档切成一条 subrow，
// 行尾剩下的空档再补一条；最后把所有 subrow 吸附到 site 网格并丢弃过窄的碎片。
//! 疑似问题：该函数已被 initializeSubrows()（走 PlaceDB::removeBlockedSite()）取代。
//! 区别在于本函数要求 obstacles 按 x 严格有序、且只做「全行扫描」，
//! 而 placedb 版本按 interval 做增量切分并顺带处理 site 对齐，逻辑更健壮。
// ----------------------------------------------------------------------------
void AbacusLegalizer::initializeSubrowsOld()
{
    // subrows actually equals intervals, so this function should be a member function of the placedb, which can be called by legalizer or detailed placer
    double siteStep = placeDB->dbSiteRows.front().step; //! assume step for all rows are identical!

    for (SiteRow curRow : placeDB->dbSiteRows) // generate subrows from dbRows
    {
        float currentX = curRow.start.x;
        float currentRowBottom = curRow.start.y; // start.y == curRow.bottom
        float currentRowTop = curRow.bottom + curRow.height;

        assert(curRow.bottom == curRow.start.y);
        for (Obstacle curObstacle : obstacles)
        {
            // assume all obstacles are inside the core region
            // ! and are ordered from left to right!

            // assuming row and obstacles are like this:
            // O for obstacle and --- for row
            // ---------------------
            //  OOO  OOOO  OO    OOOO  OOOOO
            //! curRow.end.x might be to the left of an obstacle's right boundary! although in almost all cases, curRow.end.x is to the right of all obstacles(only consider obstacles in the core region, see initializeObstacles)

            if (float_greaterorequal(curRow.start.x, curObstacle.ur.x)) // not overlap in x direction, following obstacles might overlap
            {
                continue;
            }
            if (float_greaterorequal(curObstacle.ll.x, curRow.end.x)) // not overlap in x direction, following obstacles will not overlap as well
            {
                break;
            }

            if (float_greater(currentRowTop, curObstacle.ll.y) && float_less(currentRowBottom, curObstacle.ur.y)) // overlap in y direction
            {
                // end.x might be less than start.x when:
                //    ------
                //   OOOOO
                AbacusRow newRow;
                newRow.start.x = currentX;
                newRow.start.y = currentRowBottom;

                newRow.end.x = curObstacle.ll.x; //!
                newRow.end.y = currentRowBottom;

                newRow.bottom = currentRowBottom;
                newRow.height = curRow.height;
                newRow.step = siteStep;

                currentX = curObstacle.ur.x;

                subrows.push_back(newRow);
            }
            if (float_greaterorequal(currentX, curRow.end.x))
            {
                break;
            }
        }
        // add the last subrow when row and obstacles are like this:
        // O for obstacle and --- for row
        // ---------------------
        //  OOO  OOOO  OOOO
        // or like this:
        // -------
        //           OOO
        if (float_less(currentX, curRow.end.x))
        {
            AbacusRow newRow;
            newRow.start.x = currentX;
            newRow.start.y = currentRowBottom;

            newRow.end.x = curRow.end.x; //!
            newRow.end.y = currentRowBottom;

            newRow.bottom = currentRowBottom;
            newRow.height = curRow.height;
            newRow.step = siteStep;

            subrows.push_back(newRow);
        }
    }

    // 后处理：把每条 subrow 的两端吸附到 site 网格，并丢弃连一个 site 都放不下的碎片。
    // 为什么必须吸附：cell 的 x 只能落在 site 边界上，若 subrow 端点不在网格上，
    // 贴着 subrow 端点放置的 cell 就会错位（这也是 removeBlockedSite 第 3 步做的事）。
    for (auto iter = subrows.begin(); iter != subrows.end();)
    {
        //! align subrows to sites after generating subrows, see ntuplace: FixFreeSiteBySiteStep(). Here we need to update the end and start of a site row, so end.x-start.x is an positive integer multiple of site step(site width)
        //? should start.x and end.x be integers too??? check ntuplace
        // 左端向上取整、右端向下取整（都相对 coreRegion.ll.x 对齐），
        // 保证吸附后的宽度一定是 site step 的正整数倍
        double subrowWidth = iter->end.x - iter->start.x;
        // subRowWidth might be less than 0 when:
        //    ------
        //   OOOOO
        if (float_less(subrowWidth, siteStep)) // assume iter->step > 0!
        {
            iter = subrows.erase(iter);
        }
        else
        {
            double newLeft = ceil((iter->start.x - placeDB->coreRegion.ll.x) / siteStep) * siteStep + placeDB->coreRegion.ll.x;
            double newRight = floor((iter->end.x - placeDB->coreRegion.ll.x) / siteStep) * siteStep + placeDB->coreRegion.ll.x;
            double newSubrowWidth = newRight - newLeft;
            assert(newSubrowWidth >= siteStep);
            if (float_greater(newSubrowWidth, 0.0))
            {
                iter->start.x = newLeft;
                iter->end.x = newRight;
            }
            else if (float_less(newSubrowWidth, 0.0))
            {
                // newSubrowWidth should >= 0.0
                cerr << "sub row new width < 0 when it should not\n";
                exit(0);
            }
            iter++;
        }
    }

    // sort subrows by y coordinate!
    sort(subrows.begin(), subrows.end(), [=](AbacusRow a, AbacusRow b)
         {
         if (!float_equal(a.bottom , b.bottom))
            {
                return float_less(a.bottom , b.bottom);
            }
            else if (!float_equal(a.start.x , b.start.x))
            {
                return float_less(a.start.x ,b.start.x);
            }
            else
            {
                // terminals/macros overlap!!
                cerr<<"SUBROWS OVERLAP!\n";
                exit(0);
            } });
}

// ----------------------------------------------------------------------------
// AbacusLegalizer::initializeSubrows —— 生成 Abacus 真正使用的 subrow
//
// 为什么先调 removeBlockedSite()：site row 会被 macro 和 terminal 切断，
// 只有先把这些被挡住的 site 挖掉，得到的 subrow 才是「真的能放 cell 的区间」。
// 这一步放在 PlaceDB 而不是本类，是因为 interval 是全局几何信息，
// legalizer 和后续 detailed placer 都要用（原作者注释也点明了这一点）。
//
// 前置条件（重要）：调用前 macro 必须已经合法化完毕。
// 若 macro 之间还有重叠，removeBlockedSite() 切出来的 interval 就会偏大，
// 导致 std cell 被放到 macro 底下的非法区域。
//
// 最后按 (bottom, start.x) 升序排序：双向搜索 + |Δy| 剪枝都依赖 bottom 单调。
// ----------------------------------------------------------------------------
void AbacusLegalizer::initializeSubrows()
{
    placeDB->removeBlockedSite();
    double siteStep = placeDB->dbSiteRows.front().step; //! 仍假设所有 row 的 site step 相同
    for (SiteRow curRow : placeDB->dbSiteRows) // generate subrows from dbRows
    {
        for (Interval curInterval : curRow.intervals)
        {
            AbacusRow newRow;
            newRow.start.x = curInterval.start;
            newRow.start.y = curRow.bottom;

            newRow.end.x = curInterval.end; //!
            newRow.end.y = curRow.bottom;

            newRow.bottom = curRow.bottom;
            newRow.height = curRow.height;
            newRow.step = siteStep;

            subrows.push_back(newRow);
        }
    }
    // sort subrows by y coordinate!
    sort(subrows.begin(), subrows.end(), [=](AbacusRow a, AbacusRow b)
         {
         if (!float_equal(a.bottom , b.bottom))
            {
                return float_less(a.bottom , b.bottom);
            }
            else if (!float_equal(a.start.x , b.start.x))
            {
                return float_less(a.start.x ,b.start.x);
            }
            else
            {
                // terminals/macros overlap!!
                cerr<<"SUBROWS OVERLAP!\n";
                exit(0);
            } });
}

// ============================================================================
// AbacusLegalizer::placeRow —— 把一个 cell 放到第 bestRow 条 subrow 上，并返回其位移代价
//
// trial == ABACUS_TRIAL（true）：**试算**。在 subrow 的一份拷贝上走完整流程，
//   不改动真实状态，因此可以对多个候选行各试一次再挑最优（legalization 的行分配就靠这个）。
// trial == ABACUS_FINAL（false）：**落地**。直接操作 subrows[bestRow] 本身，
//   并累加该行的已占用宽度 width（供后续 cell 的容量预筛使用）。
//
// 行内分两种情形（这就是 Abacus 的核心分支）：
//   · 行内还没有 cluster，或新 cell 的 x 已在最后一个 cluster 右端之外
//     → 它与已有 cell 不冲突，可以**独立成簇**；
//   · 否则它与最后一个 cluster 在 x 上重叠
//     → 直接并入该 cluster，再由 collapse() 决定是否需要继续向前合并。
// 注意两种情形最后都要 collapse()：原作者特地标注这一步是 Abacus 论文里没写、
// 但工程上必须补的（论文只说「插入后重算」，没说要处理越界夹紧与递归合并）。
// ============================================================================
double AbacusLegalizer::placeRow(Module *cell, int bestRow, bool trial)
{
    AbacusRow *tempRowPointer;
    AbacusRow tempRow;
    if (trial == ABACUS_TRIAL)
    {
        // cout<<"try "<< cell->name<<"at"<< bestRow<<endl;
        // 拷贝整条 subrow（含 clusters 深拷贝），后续所有修改只发生在这份副本上。
        // 代价是每次试算 O(行内 cluster 数)，换来的是「可安全地多方案比较」。
        tempRow = subrows[bestRow];
        tempRowPointer = &tempRow;
        // for trial place, would not actually modify subrows
    }
    else if (trial == ABACUS_FINAL)
    {
        tempRowPointer = &subrows[bestRow];
        // 只有落地时才累加占用宽度：这是 legalization 容量预筛的唯一依据
        tempRowPointer->width += cell->getWidth();
        // final place
    }
    else
    {
        cerr << "Not trial nor final\n";
        exit(0);
    }
    segmentFaultCP("K1"); // 调试用的分段打印，仅当命令行带 segDebug 时才真正输出
    // 判定「能否独立成簇」：行空，或新 cell 的左边界已在最后一个 cluster 右端（x+w）之后
    if (tempRowPointer->clusters.empty() || float_lessorequal(tempRowPointer->clusters.back().x + tempRowPointer->clusters.back().w, cell->getLL_2D().x))
    {
        AbacusCellCluster newCluster;
        newCluster.index = tempRowPointer->clusters.size();                                                //! 0 for the first cluster
        // 新 cluster 的初始位置：取 cell 的理想 x 并吸附到 site 网格
        //! 疑似问题：这个吸附表达式是恒等变换，并未真正对齐 site。
        //! int(x+0.5) 先四舍五入成整数，再 /step*step 只是把除法还原，
        //! 结果恒等于 int(x+0.5)；当 step != 1（例如 step=2、x=11）时会落在非 site 边界上。
        //! 正确写法应像 collapse() 那样先除后取整再乘：(int)(x / step) * step。
        //! 另外这里以绝对 0 为对齐基准，而 placedb::removeBlockedSite 是以 coreRegion.ll.x 为基准，
        //! 若 coreRegion.ll.x 不是 step 的整数倍，两者也会错位。
        newCluster.x = int((int(cell->getLL_2D().x + 0.5) / tempRowPointer->step) * tempRowPointer->step); // !! need to make sure that position of clusters align with sites!!!!
        //? site step: double instead of int, potential precision problems? although step is almost always 1 (never seen any case in which step != 1 so far)
        //! boundary check
        // 夹紧到 subrow 的可用区间内：左不过 start.x，右不过 end.x。
        // 这一步让「放在区间外」的 cell 被拉回区间内，是合法化的硬约束。
        if (float_less(newCluster.x, tempRowPointer->start.x))
        {
            newCluster.x = tempRowPointer->start.x;
        }
        if (float_greater(newCluster.x + cell->getWidth(), tempRowPointer->end.x))
        {
            newCluster.x = tempRowPointer->end.x - cell->getWidth();
        }

        tempRowPointer->clusters.push_back(newCluster);
        tempRowPointer->addCell(newCluster.index, cell); // addCell: actually always add cell to the last cluster, that is, cluster.back()
        tempRowPointer->collapse(newCluster.index);      //! this is missing in the abacus paper!
    }
    else
    {
        // cout<<cell->name<<"at here\n";
        // 与最后一个 cluster 重叠 → 并入它，由 collapse 决定要不要向前合并
        //! 疑似问题：这里用 clusters.back().index 作为下标。collapse() 合并后会 erase，
        //! 其后所有 cluster 的 index 字段都未回填，可能与实际下标不符（偏大）；
        //! 一旦发生合并，这里就可能索引越界。稳妥写法应直接用 clusters.size()-1。
        tempRowPointer->addCell(tempRowPointer->clusters.back().index, cell); // addCell: actually always add cell to the last cluster, that is, cluster.back()
        tempRowPointer->collapse(tempRowPointer->clusters.back().index);
    }
    //! operate with tempRowPointer
    // todo:calculate cell location and return cost!
    // here, cell must be in the last cluster
    float x = tempRowPointer->clusters.back().x; // !! need to make sure that position of clusters align with sites!!!! Align in collapse and when creating new clusters
    // 代价 = 该 cell 从全局布局位置到合法位置的**欧氏位移**。
    // 在 cluster 内按宽度顺序累加，找到 cell 自己的 x 后即可算出。
    // 用欧氏距离（而非平方）是为了与 legalization 中 |Δy| 的剪枝下界保持同量纲可比。
    for (Module *curCell : tempRowPointer->clusters.back().cells)
    {
        if (curCell == cell)
        {
            double cost = sqrt((curCell->getLL_2D().x - x) * (curCell->getLL_2D().x - x) + (curCell->getLL_2D().y - tempRowPointer->bottom) * (curCell->getLL_2D().y - tempRowPointer->bottom));
            return cost;
            //? sqrt or not?
        }
        else
        {
            x += curCell->getWidth();
        }
    }
    cerr << "cell not found in cluster??\n";
    exit(0);
}

// ============================================================================
// AbacusRow::addCell —— 把一个 cell 追加进指定 cluster，并**增量**维护 e / q / w
//
// 三个量的增量公式（Abacus 论文）：
//   e += 1                  cell 计数
//   q += cell.x - w(旧)     把该 cell 的理想 x 折算到「cluster 左端为原点」的坐标系后累加。
//                           减去旧 w 是因为该 cell 前面已经排了 w 那么宽的 cell。
//   w += cell.width        cluster 总宽度
// 于是 cluster 的理想左端 x = q / e（让所有 cell 的位移之和最小）。
// 之所以能增量维护：cell 是按全局 x 升序加入的，追加到尾部即保持几何顺序。
// ============================================================================
void AbacusRow::addCell(int clusterIndex, Module *cell)
{
    assert(clusterIndex < clusters.size());
    AbacusCellCluster &c = clusters[clusterIndex];
    c.cells.push_back(cell);
    c.e++;
    c.q += cell->getLL_2D().x - c.w;
    c.w += cell->getWidth();
}

// ============================================================================
// AbacusRow::addCluster —— 把后一个 cluster 合并进前一个（论文的 merge 步骤）
//
// 为什么需要合并：两个相邻 cluster 若在 x 上互相挤压，它们就必须贴在一起，
// 此后只能作为一个刚体整体移动，因此要合成一个 cluster 重新求最优位置。
//
// 关键在 q 的合并公式：q' += q - e * w_prime。
// 被并入的 cluster c 的 q 是以「c 的左端」为原点折算的，
// 而合并后 c 的左端变成了「cPrime 左端 + w_prime」，
// 所以每个 cell 的偏移量都要再减掉 w_prime，共 e 个 cell，即 - e * w_prime。
// ============================================================================
void AbacusRow::addCluster(int predecessorIndex, int clusterIndex)
{
    // cout<<predecessorIndex<<" fafa "<<clusterIndex<<endl;
    // cout<<clusters.size()<<endl;
    AbacusCellCluster &cPrime = clusters[predecessorIndex];
    AbacusCellCluster &c = clusters[clusterIndex];
    cPrime.cells.insert(cPrime.cells.end(), c.cells.begin(), c.cells.end());
    cPrime.e += c.e;
    cPrime.q += c.q - c.e * cPrime.w;
    cPrime.w += c.w;
    //! 疑似问题：合并后没有同步被删 cluster 之后各元素的 index 字段（见 legalizer.h 同款标注），
    //! 调用方 collapse() 下一行的 assert 也因此只在「恰好没发生过合并」时才有意义。
}

// ============================================================================
// AbacusRow::collapse —— Abacus 行内合法化的核心：重算 cluster 的 x 并在必要时递归向前合并
//
// 三步：
//   ① 按 x = (q/e) 求 cluster 的理想左端，并**向下吸附到 site 网格**
//      （先除 step、取整、再乘 step —— 这才是正确的吸附顺序）。
//   ② 夹紧到 subrow 的 [start.x, end.x]：保证 cell 不会溢出可用区间。
//      夹到左端/右端后，cluster 就可能被「推」到与前驱重叠。
//   ③ 若与前驱 cluster 重叠（前驱右端 > 自己左端），就合并成刚体，
//      再对合并后的 cluster 递归 collapse —— 因为合并会让总宽度变大，
//      可能继续顶到更前面的 cluster，所以必须一路向前传播。
//
// 递归终止：predecessorIndex 单调递减，最远到 0；且每次合并都会减少 cluster 总数。
// 这一步是原作者标注「Abacus 论文里缺失、但工程上必须补」的部分。
// ============================================================================
void AbacusRow::collapse(int clusterIndex)
{
    assert(clusterIndex < clusters.size());
    AbacusCellCluster &c = clusters[clusterIndex];
    // place c
    //! 先除后取整再乘，才是真正的向下吸附到 site 网格（对比 placeRow 里那处恒等变换）
    c.x = (int)(1.0 * c.q / c.e / step) * step;
    if (float_less(c.x, start.x))
    {
        c.x = start.x;
    }
    if (float_greater(c.x + c.w, end.x))
    {
        c.x = end.x - c.w;
    }
    int predecessorIndex = clusterIndex - 1;
    // cout<<"pindex here is: "<<predecessorIndex<<endl;
    if (predecessorIndex >= 0)
    {
        // 注意这里是**值拷贝**：比较用的是合并前的前驱状态，正好是我们想要的语义
        AbacusCellCluster cPrime = clusters[predecessorIndex];
        // 前驱右端 > 自己左端 ⇒ 两者重叠，必须合并成一个刚体
        //! 疑似问题：用裸的 > 比较浮点，未使用 float_greater 的 EPS 容差；
        //! 两者恰好相切（差值为 0）时不应合并，当前写法正确，但浮点误差下容易误判。
        if (cPrime.x + cPrime.w > c.x) // actually comparing int here
        {
            addCluster(predecessorIndex, clusterIndex);
            assert(clusters[clusterIndex].index == clusterIndex); //!
            clusters.erase(clusters.begin() + clusterIndex);

            // 合并后宽度变大，可能继续顶到更前面的 cluster，递归向前传播
            collapse(predecessorIndex);
        }
    }
}

// ============================================================================
// SAMacroLegalizer::legalization —— macro 合法化（mLG）的主入口
//
// 结构（对应 ePlace-MS 论文）：
//   外层 j 循环 = mLG iteration：每轮重新设定 SA 的初温与初始扰动半径，
//     跑完一整轮退火后，按 beta 放大 miuO（重叠惩罚权重），
//     让「消除 macro 重叠」这个硬约束逐轮变强，最终逼出无重叠解。
//   内层 = SAMacroLegalization() 里的 k 循环（SA iteration）。
//
// SA 参数的标定逻辑（这段基本照搬 RePlAce 的 macro.cpp）：
//   初温 T0 由「希望多大的代价恶化仍有 50% 概率被接受」反推：
//     由 exp(-Δ/T0) = 0.5 得 T0 = Δ / ln2，这里 Δ 取初始代价的 3%。
//   终温同理取 0.01%，再由 T0·coef^kLimit = T_last 反推出每步的温度衰减系数。
//   扰动半径 r 从「约 5% 的核区尺寸 / sqrt(macro 数)」线性衰减到 1，
//     即：早期大范围搜索、后期只做微调（模拟退火的经典 schedule）。
//
// 终止判据：overlapFree，即「macro 总面积 − macro 并面积」降到 0。
//   用并面积而非两两求交，是因为并面积可以用扫描线 O(n log n) 求出（见 RectangleAreaSolution）。
// ============================================================================
void SAMacroLegalizer::legalization()
{
    initialization();

    for (int j = 0; j < jLimit && !overlapFree; j++) // outter loop(mLG iteration, see ePlace-MS paper)
    {
        printf(
            "\
ITER mLG: %d\n\
    HPWL=%f\n\
    CELLCOVERED=%f\n\
    OVLP=%d\n\
",
            j + 1, totalHPWL, totalCellAreaCovered,
            totalMacroOverlap);

        //! 1. initialize parameters(t and r), this part of code is basically copied from RePlAce, macro.cpp

        // 每轮把「可接受恶化的幅度」也按 beta 放大，等价于让退火更激进
        float sa_init_neg_rate = 0.03 * pow(beta, (float)j);
        // 3% cost increase will be
        // accepted in 50% probability

        float sa_last_neg_rate = 0.0001 * pow(beta, (float)j);
        // 0.01% cost increase will be
        // accepted in 50% probability

        //! 疑似问题：注释写 0.01%，代码是 0.0001（即 0.01% 的 1/100，实为 0.01‰）。
        //! 不影响正确性，但与注释不符。
        float sa_init_t = sa_init_neg_rate / log(2.0); // based on the equation that
                                                       // exp(-1.0*sa_init_neg_rate/
                                                       // sa_init_t) = 0.5

        SAtemperature = sa_init_t;
        // 0.1 ; // 400 ; // from Howard's ICCCAS 13 paper

        SAtemperatureCoef = pow((sa_last_neg_rate / sa_init_neg_rate),
                                1.0 / (float)kLimit); // make sure sa_t(last) equals
                                                      // its expected value
        VECTOR_2D sa_n;
        VECTOR_2D sa_ncof;
        VECTOR_2D max_sa_r;
        // 半径按 sqrt(macro 数) 缩放：macro 越多，单个 macro 的合理活动范围越小
        sa_n.x = sa_n.y = sqrt(float(dbMacros.size()));

        sa_ncof.x = 0.05 * pow(beta, (float)j);
        sa_ncof.y = 0.05 * pow(beta, (float)j);

        max_sa_r.x = placeDB->coreRegion.getWidth() / sa_n.x * sa_ncof.x;
        max_sa_r.y = placeDB->coreRegion.getHeight() / sa_n.y * sa_ncof.y;

        r.x = max_sa_r.x;
        r.y = max_sa_r.y;

        u.x = u.y = 1.0;

        // 每步收缩量：kLimit 步后 r 从 max_sa_r 线性降到 1
        sa_r_stp.x = (max_sa_r.x - 1.0) / (float)kLimit;
        sa_r_stp.y = (max_sa_r.y - 1.0) / (float)kLimit;

        //! 2. SA iterations
        SAMacroLegalization();

        //! 3. update parameters
        // 只放大重叠权重 miuO：线长项 miuD 保持不变，
        // 这样后期会牺牲一定线长来优先消除重叠
        miuO *= beta;
        // no update for miuD
        // cout << getMacrosOverlap() << " fff " << totalMacroArea - getAreaCoveredByMacros() << endl;
        // 用「总面积 − 并面积」精确重算一次重叠量（增量维护的 totalMacroOverlap 有截断误差）
        totalMacroOverlap = totalMacroArea - getAreaCoveredByMacros();
        if (totalMacroOverlap <= 0)
        {
            overlapFree = 1;
        }
        if (gArg.CheckExist("fullPlot"))
        {
            PLOTTING::plotCurrentPlacement("mLG ite_" + to_string(j + 1), placeDB);
        }
    }

    totalMacroOverlap = totalMacroArea - getAreaCoveredByMacros();
    //! 疑似问题：多出来一个空语句 `;`（紧跟上一行的分号之后），无害但是笔误。
    ;
    printf(
        "\
    mLG done: \n\
    FINAL HPWL=%d\n\
    FINAL OVLP=%d\n\
",
        int(placeDB->calcHPWL()),
        totalMacroOverlap);
    // assert(totalMacroOverlap == 0);
}

// ============================================================================
// SAMacroLegalizer::SAMacroLegalization —— 内层 SA 循环（k 循环）
//
// 每个 k 步：
//   · 对**每个** macro 各做一次扰动尝试（内层 i 循环），
//     即「一轮 = 全体 macro 各动一次」，这样各 macro 的移动机会均等；
//   · 然后温度乘衰减系数、扰动半径线性收缩：早期大步搜索、后期局部微调。
// 一旦 overlapFree 被置起，两层循环都会立刻退出，不做无用功。
// ============================================================================
void SAMacroLegalizer::SAMacroLegalization()
{
    int innerLoopCount = dbMacros.size();
    // cout << "innerloopcount(macro count): " << innerLoopCount << endl;
    // cout << "kLimit: " << kLimit << endl;
    printf(
        "  -- ITER, TEMP, Rx, Ry, HPWL , DEN, OVLP\n");
    double ite100time;
    for (int k = 0; k < kLimit && !overlapFree; k++) // inner loop(SA iteration), see ePlace-MS paper
    {
        time_start(&ite100time);
        for (int i = 0; i < innerLoopCount && !overlapFree; i++) // there is one more for in RePlAce code, why???
        {
            SAperturb();
        }
        // update temperature and r
        SAtemperature *= SAtemperatureCoef;
        r.x -= sa_r_stp.x;
        r.y -= sa_r_stp.y;

        if (k % 100 == 99)
        {
            time_end(&ite100time);
            printf(
                "  -- %d, %.2e, %.2e, %.2e, %.8e, %.2e, "
                "\033[36m%d\033[0m\n",
                k / 100 + 1, SAtemperature, r.x, r.y,
                // tot_mac_hpwl ,
                totalHPWL, totalCellAreaCovered, totalMacroOverlap);
            cout << "ite100time: " << ite100time << endl;
        }
    }
}

// ============================================================================
// SAMacroLegalizer::initialization —— mLG 的准备工作，顺序有依赖
//   initializeMacros()   （先把 macro 坐标离散化，否则 bin 索引和并面积都无从算起）
//   → initializeBins()   （建 bin 网格并冻结 cell/terminal/base 密度）
//   → initializeCost()   （算三项初始代价，并判断是否一开始就无重叠）
//   → initializeSAparams()（用初始代价的量纲反推 miuD / miuO，必须最后做）
// ============================================================================
void SAMacroLegalizer::initialization()
{
    initializeMacros();
    initializeBins();
    initializeCost();
    initializeSAparams();
}

// ============================================================================
// SAMacroLegalizer::initializeMacros —— 收集 macro 并把它们的坐标离散化
//
// 为什么要离散化（照搬 RePlAce）：
//   · x 取整：后续「macro 并面积」是用扫描线 + 线段树按**整数坐标**求的
//     （见 RectangleAreaSolution），浮点会破坏 y 值去重与二分。
//   · y 对齐到行高网格：macro 的高度通常是行高的整数倍，
//     把 y 吸附到「coreRegion.ll.y + k * commonRowHeight」，
//     才能保证它正好压住整数条 site row，让 removeBlockedSite() 切得干净。
//
// 同时在这一遍里累加 totalMacroArea（macro 总面积）与 totalCellArea（std cell 总面积），
// 前者是判断重叠的基准，后者参与 miuD 标定。
//! 注意：totalMacroArea 用 double 而非 int —— 原作者注释指出用 int 会因精度/溢出出 bug。
// ============================================================================
void SAMacroLegalizer::initializeMacros()
{
    for (Module *curNode : placeDB->dbNodes)
    {
        assert(curNode);
        if (curNode->isMacro)
        {
            dbMacros.push_back(curNode);
            totalMacroArea += curNode->getArea(); // followed RePlAce code
            //! Discretization of macro coordinates, followed RePlAce code
            POS_2D legalLL;
            legalLL.x = (int)(curNode->getLL_2D().x + 0.5); // set macro center to integer here, but why according to row height?
            // y：先加半个行高再整除行高（即「四舍五入到最近的行边界」），
            // 基准是 coreRegion.ll.y 而不是绝对 0，最后再加回基准
            legalLL.y = ((int)((curNode->getLL_2D().y + 0.5 * placeDB->commonRowHeight - placeDB->coreRegion.ll.y /* 1.0 * ROW_Y0 */) / (double)placeDB->commonRowHeight)) * placeDB->commonRowHeight + (int)(placeDB->coreRegion.ll.y + 0.5) /* ROW_Y0 */;
            placeDB->setModuleLocation_2D(curNode, legalLL.x, legalLL.y);
        }
        else
        {
            totalCellArea += curNode->getArea();
        }
    }
    cout << "total macro area: " << fixed << totalMacroArea << endl;
}

// ============================================================================
// SAMacroLegalizer::initializeBins —— 建立 bin 网格并预先冻结各类面积
//
// 目的（对应 ePlace 论文的 bin 结构）：后续要反复问「某个 macro 压住了多少 cell 面积」，
// 逐个 bin × 逐个 cell 求交太慢，所以先把每个 bin 的
//   cellArea（与 std cell 的重叠面积）、terminalArea（与 terminal 的重叠）、
//   baseArea（不在任何 site row 上的「不可用」面积，按 targetDensity 折算）
// 一次性算好，合成 nonMacroDensity 冻结下来。
// 因为 mLG 期间 std cell 是**不动**的，所以这些密度只需算一次，
// 之后 macro 怎么移动都只要「重叠面积 × 密度」即可，这是 mLG 能跑起来的关键。
//
// 四步：① 定 bin 维数 → ② 分配 bin 并填几何 → ③ terminal 面积 → ④ base 面积 + cell 面积
// ============================================================================
void SAMacroLegalizer::initializeBins()
{
    segmentFaultCP("mLGbinInit");
    //! Bins are allocated on coreRegion! RePlAce code bin.cpp line 337 and bookshelfIO.cpp, also see ePlace paper
    ////////////////////////////////////////////////////////////////
    // calculate bin dimension and size
    ////////////////////////////////////////////////////////////////
    //! this code for bin dimension calculation is copied directly from RePlAce
    segmentFaultCP("mLGbinDim");

    int nodeCount = placeDB->dbNodes.size();
    float nodeArea = (totalCellArea + totalMacroArea);

    float coreRegionWidth = placeDB->coreRegion.getWidth();
    float coreRegionHeight = placeDB->coreRegion.getHeight();
    float coreRegionArea = float_mul(coreRegionWidth, coreRegionHeight);

    float averageNodeArea = 1.0 * nodeArea / nodeCount;
    float idealBinArea = averageNodeArea / targetDensity;
    // float idealBinArea = averageNodeArea / 1.0;

    int idealBinCount = INT_CONVERT(coreRegionArea / idealBinArea);

    bool isUpdate = false;
    // bin 维数取 2 的幂（便于后续用移位/整除定位 bin），并夹在 [4, 1024] 之间：
    // bin 太粗则密度近似失真，太细则代价评估变慢。
    // bin dimension upper bound: 1024 rather than 2048
    // for (int i = 1; i <= 10; i++)
    // { //! 4*4,8*8,16*16,32*32..., 2048*2048
    //! 疑似问题：候选维数从 (2<<1)=4 起跳，因此当 idealBinCount < 16 时
    //! 循环一次都命中不了，会直接落到下面的 1024 兜底——
    //! 也就是「cell 很少的小设计反而会得到最密的 1024×1024 网格」，与直觉相反。
    for (int i = 1; i < 10; i++)
    { //! 4*4,8*8,16*16,32*32..., 1024*1024
        if ((2 << i) * (2 << i) <= idealBinCount &&
            (2 << (i + 1)) * (2 << (i + 1)) > idealBinCount)
        {
            binDimension.x = binDimension.y = 2 << i;
            isUpdate = true;
            break;
        }
    }
    if (!isUpdate)
    {
        binDimension.x = binDimension.y = 1024; //!
    }

    cout << BLUE << "Macro legalizer bin dimension: " << binDimension << "\ncoreRegion width: " << coreRegionWidth << "\ncoreRegion height: " << coreRegionHeight << RESET << endl;

    binStep.x = float_div(coreRegionWidth, binDimension.x);
    binStep.y = float_div(coreRegionHeight, binDimension.y);

    cout << BLUE << "Macro legalizer bin step: " << binStep << RESET << endl;
    ////////////////////////////////////////////////////////////////
    // add bins
    ////////////////////////////////////////////////////////////////
    segmentFaultCP("addBin");
    // x stored in first dimension of vector
    //! bin index:
    //! # 5
    //! # 4
    //! # 3
    //! # 2
    //! # 1
    //! # 0
    //! #  0 1 2 3 4 5

    //?
    bins.resize(binDimension.x);

    double addBinTime;
    double terminalDensityTime;
    double baseDensityTime;

    time_start(&addBinTime);

    //! 疑似问题：所有 bin 都用 new 分配，但 SAMacroLegalizer 没有析构函数、也没有任何 delete，
    //! 属于内存泄漏（进程级工具可接受，但若在同一进程内重复初始化会持续累积）。
    for (int i = 0; i < binDimension.x; i++)
    {
        bins[i].resize(binDimension.y);
        for (int j = 0; j < binDimension.y; j++)
        {
            bins[i][j] = new SAMacroLegalizationBin_2D();
            // cout<<"adding bin "<<i<<","<<j<<endl;
            //!!! +db->coreRegion.ll.x to get coordinates!!
            bins[i][j]->ll.x = i * binStep.x + placeDB->coreRegion.ll.x;
            bins[i][j]->ll.y = j * binStep.y + placeDB->coreRegion.ll.y;

            bins[i][j]->width = binStep.x;
            bins[i][j]->height = binStep.y;

            bins[i][j]->ur.x = bins[i][j]->ll.x + bins[i][j]->width;
            bins[i][j]->ur.y = bins[i][j]->ll.y + bins[i][j]->height;

            bins[i][j]->area = binStep.x * binStep.y;

            bins[i][j]->center.x = bins[i][j]->ll.x + (float)0.5 * bins[i][j]->width;
            bins[i][j]->center.y = bins[i][j]->ll.y + (float)0.5 * bins[i][j]->height;
        }
    }

    time_end(&addBinTime);
    cout << "mLG bin add time: " << addBinTime << endl;
    ////////////////////////////////////////////////////////////////
    // terminal density calculation, calculate here because they are terminals and only needed to be considered once
    ////////////////////////////////////////////////////////////////
    segmentFaultCP("mLGterminalArea");
    VECTOR_2D_INT binStartIdx;
    VECTOR_2D_INT binEndIdx;

    time_start(&terminalDensityTime);

    // terminal 是**不可移动**的，它们占住的面积对所有 bin 而言是常量，
    // 因此只需在初始化时累加一次，之后 macro 怎么动都不用重算。
    for (Module *curTerminal : placeDB->dbTerminals)
    {
        //! only consider terminals inside the coreRegion
        //! assume no terminal would have a part inside the coreRegion and a part outside the coreRegion
        // 完全在 coreRegion 之外的 terminal 直接跳过（跨边界的部分被整体忽略，是上述假设的体现）
        if (curTerminal->getLL_2D().x < placeDB->coreRegion.ll.x || curTerminal->getLL_2D().x + curTerminal->getWidth() > placeDB->coreRegion.ur.x)
        {
            continue;
        }
        if (curTerminal->getLL_2D().y < placeDB->coreRegion.ll.y || curTerminal->getLL_2D().y + curTerminal->getHeight() > placeDB->coreRegion.ur.y)
        {
            continue;
        }

        // INT_DOWN 即 (int) 截断：把绝对坐标换算成「相对 coreRegion」后再除以 bin 步长
        //! 疑似问题：INT_DOWN 是向零截断而非向下取整，对负数会偏向 0；
        //! 这里靠上面「只处理完全在 coreRegion 内的对象」保证了被除数非负，
        //! 但宏/单元一旦越界（下面 cell 循环没有做同样的过滤），索引就会偏。
        binStartIdx.x = INT_DOWN((curTerminal->getLL_2D().x - placeDB->coreRegion.ll.x) / binStep.x);
        binEndIdx.x = INT_DOWN((curTerminal->getUR_2D().x - placeDB->coreRegion.ll.x) / binStep.x);

        binStartIdx.y = INT_DOWN((curTerminal->getLL_2D().y - placeDB->coreRegion.ll.y) / binStep.y);
        binEndIdx.y = INT_DOWN((curTerminal->getUR_2D().y - placeDB->coreRegion.ll.y) / binStep.y);

        // cout << curTerminal->getLL_2D() << curTerminal->getUR_2D() << endl;
        // cout << binStep << " " << db->coreRegion.ll << endl;
        assert(binStartIdx.x >= 0);
        assert(binEndIdx.x >= 0);
        assert(binStartIdx.y >= 0);
        assert(binEndIdx.y >= 0);

        // 右边界正好落在 bin 边界上时会算出 == binDimension 的索引，需钳回最后一个 bin
        if (binEndIdx.y >= binDimension.y)
        {
            binEndIdx.y = binDimension.y - 1;
        }

        if (binEndIdx.x >= binDimension.x)
        {
            binEndIdx.x = binDimension.x - 1;
        }

        // 把 terminal 与 bin 的重叠面积按 targetDensity 折算后累加进该 bin
        for (int i = binStartIdx.x; i <= binEndIdx.x; i++)
        {
            for (int j = binStartIdx.y; j <= binEndIdx.y; j++)
            {
                //! beware: density scaling!
                bins[i][j]->terminalArea += targetDensity * getOverlapArea_2D(bins[i][j]->ll, bins[i][j]->ur, curTerminal->getLL_2D(), curTerminal->getUR_2D());
            }
        }
    }

    time_end(&terminalDensityTime);
    cout << "Macro legalizer terminal density time: " << terminalDensityTime << endl;
    ////////////////////////////////////////////////////////////////
    // base density calculation, also only needed to be considered once
    ////////////////////////////////////////////////////////////////
    segmentFaultCP("mLGbaseArea");

    time_start(&baseDensityTime);

    for (int i = 0; i < binDimension.x; i++)
    {
        for (int j = 0; j < binDimension.y; j++)
        {
            // baseArea = bin 中「不在任何 site row 上」的面积，
            // 即天然放不了 cell 的死区（行间空隙、coreRegion 边缘等）。
            // 它和 terminal 一样是静态的，只需算一次。
            //! 疑似问题：这里是 O(bin 数 × site row 数) 的双重循环，
            //! 在 1024×1024 网格 + 上千条 row 的情况下会非常慢；
            //! 可以考虑按 y 先定位相关 row 再求交。
            float curBinAvailableArea = 0; // overlap area between current bin and placement rows
            for (SiteRow curRow : placeDB->dbSiteRows)
            {
                curBinAvailableArea += getOverlapArea_2D(bins[i][j]->ll, bins[i][j]->ur, curRow.getLL_2D(), curRow.getUR_2D());
            }
            // debugOutput("Bin area", bins[i][j]->area);
            // debugOutput("Available area", curBinAvailableArea);
            if (float_equal(bins[i][j]->area, curBinAvailableArea))
            {
                bins[i][j]->baseArea = 0;
            }
            else
            {
                // 同样按 targetDensity 折算，与 RePlAce 保持一致
                bins[i][j]->baseArea = targetDensity * (bins[i][j]->area - curBinAvailableArea); //! follow RePlAce bin.cpp line 433
            }
        }
    }
    time_end(&baseDensityTime);
    cout << "Base density time: " << baseDensityTime << endl;

    ////////////////////////////////////////////////////////////////
    // cell density calculation, calculate here because cells are fixed during macro legalization and only needed to be considered once
    ////////////////////////////////////////////////////////////////
    segmentFaultCP("cellArea");
    // std cell 在 mLG 期间固定不动，所以它们的 bin 密度同样只算一次。
    // 注意 macro 的面积**也**被累加进了 bins[][]->macroArea，
    // 但如 legalizer.h 中所标注，macroArea 从未参与任何代价计算。
    //! 疑似问题：这里没有像 terminal 那样过滤「部分落在 coreRegion 外」的 module，
    //! 一旦有 cell 越界，INT_DOWN 的向零截断会算出错误的 bin 索引（甚至 0），
    //! 与上面 terminal 的处理不一致。
    for (Module *curNode : placeDB->dbNodes) // ePlaceNodes: nodes and filler nodes
    {
        CRect rectForCurNode;
        rectForCurNode.ll = curNode->getLL_2D();
        rectForCurNode.ur = curNode->getUR_2D();

        VECTOR_2D_INT binStartIdx; // binStartIdx: index the index of the first bin that has overlap with a cell on X/Y direction
        VECTOR_2D_INT binEndIdx;
        binStartIdx.x = INT_DOWN((rectForCurNode.ll.x - placeDB->coreRegion.ll.x) / binStep.x);
        binEndIdx.x = INT_DOWN((rectForCurNode.ur.x - placeDB->coreRegion.ll.x) / binStep.x);

        binStartIdx.y = INT_DOWN((rectForCurNode.ll.y - placeDB->coreRegion.ll.y) / binStep.y);
        binEndIdx.y = INT_DOWN((rectForCurNode.ur.y - placeDB->coreRegion.ll.y) / binStep.y);

        if (!(binStartIdx.x >= 0))
        {
            cout << "Module pos: " << rectForCurNode.ll << " " << placeDB->coreRegion.ll << endl;
        }
        assert(binStartIdx.x >= 0);
        assert(binEndIdx.x >= 0);
        assert(binStartIdx.y >= 0);
        assert(binEndIdx.y >= 0);

        if (binEndIdx.y >= binDimension.y)
        {
            binEndIdx.y = binDimension.y - 1;
        }

        if (binEndIdx.x >= binDimension.x)
        {
            binEndIdx.x = binDimension.x - 1;
        }

        for (int i = binStartIdx.x; i <= binEndIdx.x; i++)
        {
            for (int j = binStartIdx.y; j <= binEndIdx.y; j++)
            {

                float overlapArea = getOverlapArea_2D(bins[i][j]->ll, bins[i][j]->ur, rectForCurNode.ll, rectForCurNode.ur);
                if (curNode->isMacro)
                {
                    bins[i][j]->macroArea += targetDensity * overlapArea; // macro area is never used??
                }
                else
                {
                    bins[i][j]->cellArea += overlapArea;
                }
            }
        }
    }

    // 汇总成每个 bin 的「非 macro 密度」：这是后续估算「macro 压住多少 cell」的唯一依据。
    // 命名 nonMacroDensity 是因为它刻意**不含** macroArea ——
    // 代价里的 macro 重叠由 getMacroOverlapArea() 精确两两求交负责，不走 bin 近似。
    for (int i = 0; i < binDimension.x; i++)
    {

        for (int j = 0; j < binDimension.y; j++)
        {
            bins[i][j]->nonMacroDensity = (bins[i][j]->cellArea + bins[i][j]->terminalArea + bins[i][j]->baseArea) / bins[i][j]->area;
        }
    }
}

// ============================================================================
// SAMacroLegalizer::initializeCost —— 计算三项代价的初始值
//   totalHPWL            —— 全芯片线长（顺带把所有 pin 的绝对位置初始化好）
//   totalCellAreaCovered —— 所有 macro 压住的 std cell 面积（bin 密度近似）
//   totalMacroOverlap    —— macro 总面积 − macro 并面积，即 macro 之间的重叠量
// 若一开始就没有重叠，直接置 overlapFree，外层循环会一次都不跑。
// ============================================================================
void SAMacroLegalizer::initializeCost()
{
    totalHPWL = placeDB->calcHPWL(); // all pin absolute position initialized here!

    totalCellAreaCovered = getCellAreaCoveredByAllMacros();

    totalMacroOverlap = totalMacroArea - getAreaCoveredByMacros(); // SA legalization terminates when total macro overlap == 0;
    if (totalMacroOverlap <= 0)
    {
        overlapFree = true;
    }
}

// ============================================================================
// SAMacroLegalizer::initializeSAparams —— 标定代价权重与迭代上限
//
// 为什么这样定 miuD / miuO：三项代价量纲完全不同（长度 / 面积 / 面积），
// 若直接用固定系数，某一项会被彻底淹没。这里的做法是**用初始状态反推**，
// 让三项在初始时刻大致等量：
//   miuD = HPWL / cellAreaCovered   → 使 miuD * cellAreaCovered ≈ HPWL
//   miuO = (HPWL + miuD*cellAreaCovered) / macroOverlap → 使 miuO * overlap ≈ 前两项之和
// 于是三项（线长、被压 cell、macro 重叠）在起点处于同一数量级。
// beta = 1.5：每轮外层迭代把 miuO 放大 1.5 倍，逐步把「消除重叠」变成主导目标。
// ============================================================================
void SAMacroLegalizer::initializeSAparams()
{
    jLimit = 1000;
    kLimit = 1000;

    //! 疑似问题：两处除法都没有防 0。若初始就没有 macro 重叠（totalMacroOverlap == 0），
    //! miuO 会是 inf/NaN；若没有 macro 压住任何 cell（totalCellAreaCovered == 0），
    //! miuD 同理。虽然这些情形下 overlapFree 往往已被置起、循环不会执行，
    //! 但参数本身已经被污染。
    miuD = totalHPWL / totalCellAreaCovered;
    miuO = (totalHPWL + totalCellAreaCovered * miuD) / (float)totalMacroOverlap;

    beta = 1.5;
}

// ============================================================================
// SAMacroLegalizer::getAreaCoveredByMacros —— 求所有 macro 的**并面积**
//
// 用途：macro 之间的重叠量 = 总面积 − 并面积，为 0 即表示已无重叠。
// 为什么用并面积而不是两两求交：两两求交是 O(n²)，每轮 SA 都要重算的话根本跑不动；
// 并面积用扫描线 + 线段树只要 O(n log n)（见 RectangleAreaSolution）。
// 前提是坐标必须是整数——所以 initializeMacros() 里先把 macro 离散化过。
// ============================================================================
int SAMacroLegalizer::getAreaCoveredByMacros()
{
    // todo: use segment tree to solve this problem(which is actually a classic OJ problem)

    RectangleAreaSolution solution;
    vector<vector<int>> rectangles;
    //! x and y should be integers here!
    for (Module *curMacro : dbMacros)
    {
        vector<int> llAndur;
        // 隐式把 float 截断成 int：依赖 initializeMacros() 已经做过离散化
        llAndur.push_back(curMacro->getLL_2D().x);
        llAndur.push_back(curMacro->getLL_2D().y);
        llAndur.push_back(curMacro->getUR_2D().x);
        llAndur.push_back(curMacro->getUR_2D().y);
        rectangles.push_back(llAndur);
    }
    return solution.rectangleArea(rectangles);
}

// ============================================================================
// SAMacroLegalizer::getMacrosOverlap —— O(n²) 精确两两求交的重叠面积
// 精度高于「总面积 − 并面积」（后者受整数离散化影响），但复杂度高，
// 因此只用于调试输出（当前唯一调用点已被注释掉）。
//! 疑似问题：实际上已无任何活调用点，是死代码；且返回 double 却命名为 overlap，
//! 语义上它统计的是「所有 macro 两两重叠面积之和」。
// ============================================================================
double SAMacroLegalizer::getMacrosOverlap()
{
    double ovlp = 0;

    for (int i = 0; i < dbMacros.size(); i++)
    {
        Module *mac1 = dbMacros[i];
        for (int j = i + 1; j < dbMacros.size(); j++)
        {
            Module *mac2 = dbMacros[j];
            ovlp += getOverlapArea_2D(mac2->getLL_2D(), mac2->getUR_2D(), mac1->getLL_2D(), mac1->getUR_2D());
        }
    }
    return ovlp;
}

// ============================================================================
// SAMacroLegalizer::getCellAreaCoveredByAllMacros —— 所有 macro 压住的 cell 面积之和
// 用 bin 密度近似而非精确求交（照搬 RePlAce）：
// 精确算需要「每个 macro × 每个 cell」求交，代价太高；
// 用 bin 的话只需「macro 与 bin 求交 × bin 密度」，bin 越小越精确。
// ============================================================================
float SAMacroLegalizer::getCellAreaCoveredByAllMacros()
{
    // follow the RePlAce code
    // an approximation of cell area covered by macros, because it would take too much time to get a precise value. This approximation can be very precise when bins are small enough
    //! 疑似问题：若两个 macro 互相重叠，被它们共同压住的 cell 面积会被重复计入两次，
    //! 从而高估「被压住的 cell 面积」。RePlAce 原实现同样如此。
    float cost = 0.0;
    for (Module *curMacro : dbMacros)
    {
        cost += getCellAreaCoveredByMacro(curMacro);
    }
    return cost;
}

// ============================================================================
// SAMacroLegalizer::getCellAreaCoveredByMacro —— 估算单个 macro 压住了多少 std cell 面积
//
// 做法：找出该 macro 覆盖到的所有 bin，把「macro 与 bin 的重叠面积 × bin 的
// nonMacroDensity」累加起来。密度是初始化时冻结好的，
// 所以这里的代价只与 macro 的包围盒有关，与 cell 数量无关——这正是 SA 能高频调用的前提。
// ============================================================================
float SAMacroLegalizer::getCellAreaCoveredByMacro(Module *curNode)
{
    float cost = 0;

    CRect rectForCurNode;
    rectForCurNode.ll = curNode->getLL_2D();
    rectForCurNode.ur = curNode->getUR_2D();

    VECTOR_2D_INT binStartIdx; // binStartIdx: index the index of the first bin that has overlap with a cell on X/Y direction
    VECTOR_2D_INT binEndIdx;
    binStartIdx.x = INT_DOWN((rectForCurNode.ll.x - placeDB->coreRegion.ll.x) / binStep.x);
    binEndIdx.x = INT_DOWN((rectForCurNode.ur.x - placeDB->coreRegion.ll.x) / binStep.x);

    binStartIdx.y = INT_DOWN((rectForCurNode.ll.y - placeDB->coreRegion.ll.y) / binStep.y);
    binEndIdx.y = INT_DOWN((rectForCurNode.ur.y - placeDB->coreRegion.ll.y) / binStep.y);

    if (!(binStartIdx.x >= 0))
    {
        cout << "Module pos: " << rectForCurNode.ll << " " << placeDB->coreRegion.ll << endl;
    }
    assert(binStartIdx.x >= 0);
    assert(binEndIdx.x >= 0);
    assert(binStartIdx.y >= 0);
    assert(binEndIdx.y >= 0);

    if (binEndIdx.y >= binDimension.y)
    {
        binEndIdx.y = binDimension.y - 1;
    }

    if (binEndIdx.x >= binDimension.x)
    {
        binEndIdx.x = binDimension.x - 1;
    }

    // 用冻结的 bin 密度做加权：重叠面积 × 密度 = 该 bin 里被压住的 cell 面积
    for (int i = binStartIdx.x; i <= binEndIdx.x; i++)
    {
        for (int j = binStartIdx.y; j <= binEndIdx.y; j++)
        {

            cost += getOverlapArea_2D(bins[i][j]->ll, bins[i][j]->ur, rectForCurNode.ll, rectForCurNode.ur) * bins[i][j]->nonMacroDensity;
        }
    }
    return cost;
}

// ============================================================================
// SAMacroLegalizer::getMacroOverlapArea —— 单个 macro 与其余所有 macro 的重叠面积（精确）
// 注意这里是**精确两两求交**，与「总面积 − 并面积」的全局判据不同：
// 单个 macro 的重叠量只涉及 n-1 次求交，代价可接受，且对 SA 的梯度方向更敏感。
//! 疑似问题：返回类型是 int，而 getOverlapArea_2D 返回浮点，每次 += 都会截断取整。
// ============================================================================
int SAMacroLegalizer::getMacroOverlapArea(Module *curMacro)
{
    int totalOverlap = 0;
    for (Module *otherMacro : dbMacros)
    {
        if (otherMacro == curMacro)
        {
            continue;
        }
        totalOverlap += getOverlapArea_2D(curMacro->getLL_2D(), curMacro->getUR_2D(), otherMacro->getLL_2D(), otherMacro->getUR_2D());
    }
    return totalOverlap;
}

// ============================================================================
// SAMacroLegalizer::getMacroCost —— 单个 macro 的 SA 代价（ePlace-MS 公式 30）
//   cost = HPWL(该 macro 的线长) + miuD * 被压住的 cell 面积 + miuO * 与其他 macro 的重叠面积
// 三个分量通过引用形参一并返回，供 SAperturb() 增量更新全局代价，
// 避免每接受一次扰动就重算全芯片。
// ============================================================================
float SAMacroLegalizer::getMacroCost(Module *curMacro, float &HPWL, float &cellCovered, float &macroOverlap)
{
    HPWL = placeDB->calcModuleHPWLfast(curMacro);
    cellCovered = getCellAreaCoveredByMacro(curMacro);
    macroOverlap = getMacroOverlapArea(curMacro);
    return HPWL + cellCovered * miuD + macroOverlap * miuO;
}

// ============================================================================
// SAMacroLegalizer::SAperturb —— SA 的一步：随机挑一个 macro 扰动，按 Metropolis 准则决定去留
//
// 流程：随机选 macro → 生成随机位移 delta → 移动 → 算新代价 → 接受则增量更新全局代价，
//       否则把 macro 移回原位（用中心点还原，比记住 delta 更稳妥）。
//
// 为什么 y 方向的位移要取整到 commonRowHeight 的倍数：
//   macro 通常占整数条 site row，只有让它的 y 始终落在行网格上，
//   后续 removeBlockedSite() 才能干净地把它挡住的 site 挖掉。
// ============================================================================
void SAMacroLegalizer::SAperturb()
{
    // this function randomly choose a macro and perturb its coordinate

    //! 1. choose a macro randomly
    int randomMacroIndex = 0;
    int rnd_idx = 0;
    double drnd_idx = 0;
    int mac_idx = 0;
    //! 疑似问题：mac_idx 声明后再未被使用，是死变量。

    rnd_idx = rand();
    drnd_idx = (double)rnd_idx / RAND_MAX;
    //! 疑似问题：drnd_idx 取值范围是 [0, 1]（含端点），因此 randomMacroIndex 可能等于
    //! dbMacros.size()，即越界一个。Debug 下有 assert 拦住，Release 下会读到野指针。
    randomMacroIndex = (int)(drnd_idx * (double)dbMacros.size());
    assert(randomMacroIndex < dbMacros.size());
    Module *chosenOne = dbMacros[randomMacroIndex];

    //! 2. move this macro randomly, first, decide the move vector(delta of coordinates)
    VECTOR_2D_INT delta;

    // 先记住原中心点：拒绝时需要用它把 macro 精确还原
    POS_3D curCenter = chosenOne->getCenter();

    VECTOR_2D_INT rnd;
    VECTOR_2D drnd;
    rnd.x = rand();
    rnd.y = rand();
    float RAND_MAX_INVERSE = (float)1.0 / RAND_MAX;
    // 归一化到 [-0.5, 0.5]：让位移以当前位置为中心、幅度为 r
    drnd.x = (float)rnd.x * RAND_MAX_INVERSE - 0.5;
    drnd.y = (float)rnd.y * RAND_MAX_INVERSE - 0.5;

    // x 随意取整；y 必须先按行高量化再取整，保证 macro 始终压在行网格上
    delta.x = round(drnd.x * r.x) * u.x;
    delta.y = round(drnd.y * r.y / placeDB->commonRowHeight) * placeDB->commonRowHeight * u.y;
    //! 疑似问题：delta 没有做 coreRegion 边界钳制。macro 可能被推出 coreRegion，
    //! 之后 getCellAreaCoveredByMacro() 里的 INT_DOWN 会得到负索引，触发 assert 或读越界；
    //! 且 rand() 未设种子，每次运行结果完全相同（可复现但不随机）。

    float curHPWLCost;
    float curCellCoveredCost;
    float curMacroOverlapCost;
    float curCost = getMacroCost(chosenOne, curHPWLCost, curCellCoveredCost, curMacroOverlapCost);

    placeDB->moveModule_2D(chosenOne, delta);
    // for (Pin *curModulePin : chosenOne->modulePins)
    // {
    //     // curModulePin->getAbsolutePos();
    //     curModulePin->absolutePos.x = curModulePin->module->getCenter().x + curModulePin->offset.x;
    //     curModulePin->absolutePos.y = curModulePin->module->getCenter().y + curModulePin->offset.y;
    //     curModulePin->absolutePos.z = curModulePin->module->getCenter().z;
    // }

    float newHPWLCost;
    float newCellCoveredCost;
    float newMacroOverlapCost;
    float newCost = getMacroCost(chosenOne, newHPWLCost, newCellCoveredCost, newMacroOverlapCost);

    bool accept = acceptPerturb(curCost, newCost);

    if (accept)
    {
        // 增量更新全局三项代价：只加减「被扰动 macro 这一个」的分量，
        // 避免每步都重算全芯片 HPWL / 密度
        totalHPWL += newHPWLCost - curHPWLCost;
        totalCellAreaCovered += newCellCoveredCost - curCellCoveredCost;
        //! 疑似问题：totalMacroOverlap 是 int，这里用浮点差值做 +=，每步都会截断取整。
        totalMacroOverlap += newMacroOverlapCost - curMacroOverlapCost;

        if (totalMacroOverlap <= 0 && newMacroOverlapCost <= 0 && curMacroOverlapCost > 0)
        {
            // 增量维护的重叠量刚跨过 0，用「总面积 − 并面积」精确复核一次，
            // 防止因为取整误差误判为「已无重叠」
            totalMacroOverlap = totalMacroArea - getAreaCoveredByMacros();
            // totalMacroOverlap = getMacrosOverlap();
            if (totalMacroOverlap <= 0)
            {
                overlapFree = true;
            }
        }
    }
    else
    {
        // 回滚：直接把中心点设回扰动前记录的值
        placeDB->setModuleCenter_2D(chosenOne, curCenter.x, curCenter.y);
        // for (Pin *curModulePin : chosenOne->modulePins)
        // {
        //     // curModulePin->getAbsolutePos();
        //     curModulePin->absolutePos.x = curModulePin->module->getCenter().x + curModulePin->offset.x;
        //     curModulePin->absolutePos.y = curModulePin->module->getCenter().y + curModulePin->offset.y;
        //     curModulePin->absolutePos.z = curModulePin->module->getCenter().z;
        // }
    }
}

// ============================================================================
// SAMacroLegalizer::acceptPerturb —— Metropolis 接受准则
//   代价变好（newCost < oldCost）⇒ costDelta < 0 ⇒ exp(正数) > 1 ⇒ 必然接受；
//   代价变差 ⇒ 按 exp(-Δ/T) 的概率接受，温度 T 越低越不容易接受劣解。
// 注意 costDelta 用的是**相对**变化量（除以 oldCost），
// 这样温度的物理含义就是「可接受多大的百分比恶化」，与代价的绝对量纲无关。
// ============================================================================
bool SAMacroLegalizer::acceptPerturb(float oldCost, float newCost)
{
    int rnd = 0;
    float drnd = 0;

    float expValue = 0;
    //! 疑似问题：oldCost 为 0 时会除以 0（理论上单个 macro 的 HPWL+各项至少有一项非零，
    //! 但对没有 pin 的 macro 确实可能出现 oldCost == 0）。
    float costDelta = (newCost - oldCost) / oldCost;

    rnd = rand();
    drnd = (float)rnd / RAND_MAX;

    expValue = std::exp(-1.0 * costDelta / SAtemperature);

    if (drnd < expValue)
    {
        return true;
    }
    else
    {
        return false;
    }
}
