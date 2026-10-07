#include "eplace.h"

// ============================================================================
// EPlacer_2D —— ePlace 全局布局（非线性静电场类比布局算法）
//
// 核心思想：把每个单元（含 filler）看成一团正电荷，用静电场互推来「铺开」单元。
// 目标函数：    f = WL(x) + λ · D(x)
//   · WL：线长（LSE 或 WA 平滑模型，用来逼近 HPWL 并可微）
//   · D ：密度代价 = Σ bin(电荷密度 × 电势)，由泊松方程 ∇²φ = −ρ 解出
//   · λ ：惩罚因子，迭代中自适应调整
//
// 单次迭代的梯度计算链路（totalGradientUpdate 是入口）：
//   binNodeDensityUpdate()      → 把单元/filler 面积摊到各 bin，得到 nodeDensity/fillerDensity
//   densityGradientUpdate()     → bin 密度送入 FFT 解泊松方程 → 得到电场 E 和电势 φ
//                                 再按「单元与 bin 的重叠面积 × E」加权求和得到密度梯度
//   wirelengthGradientUpdate()  → 先算溢出率 τ，由 τ 反推平滑系数 γ，再求线长梯度
//   → 合成 totalGradient，交给 Optimization/nesterov.hpp 做 Nesterov 一阶优化更新坐标
//
// 三个阶段（PLACEMENT_STAGE）：
//   mGP        宏观布局：所有单元（含 macro）+ filler 一起动
//   FILLERONLY 只重撒 filler，用于打破 mGP 后的局部拥塞
//   cGP        单元布局：macro 冻结，只有 std cell + filler 动
//
// 密度由四部分叠加（见 Bin_2D）：nodeDensity + fillerDensity + terminalDensity + baseDensity
// ============================================================================

void EPlacer_2D::setTargetDensity(float target)
{
    targetDensity = target;
    cout << padding << "Target density set at: " << targetDensity << padding << endl;
}

void EPlacer_2D::setPlacementStage(int stage)
{
    placementStage = stage;
}

//! 一次性初始化：filler 插入 → bin 网格构建 → 梯度容器分配 → 首轮梯度 → λ 初值
void EPlacer_2D::initialization()
{
    double binInitTime;
    fillerInitialization();
    time_start(&binInitTime);
    binInitialization();
    time_end(&binInitTime);
    cout << "Bin init time: " << binInitTime << endl;
    gradientVectorInitialization();
    totalGradientUpdate(); // first update
    penaltyFactorInitilization(); //! 注意：λ 的初始化依赖上面这次梯度，所以顺序不能颠倒
}

// ----------------------------------------------------------------------------
// filler（虚拟填充单元）插入
//
// filler 不连任何线网，只贡献面积/密度，作用是「占住空白区域」，
// 让密度力把真实单元往空白处推。总面积按 ePlace 论文公式(13)：
//     A_filler = A_whitespace · targetDensity − A_m
// 其中 A_whitespace = 布局行面积 − 与 terminal 的重叠面积，A_m 为单元面积（macro 按密度折算）。
// 单个 filler 尺寸取「中间 80% 单元的平均面积」，避免个别超大/超小单元拉偏。
// ----------------------------------------------------------------------------
void EPlacer_2D::fillerInitialization()
{
    segmentFaultCP("fillerInit");
    ////////////////////////////////////////////////////////////////
    // calculate whitespace area
    /////////////////////////////////////////////////////////////////

    //! whitespace: area of placement rows - overlap area between placement rows and terminals see RePlace opt.cpp whitespace_init() and ePlace Paper
    float whitespaceArea = 0;
    float totalOverLapArea = 0; // overlap area between placement rows and terminals

    //! 跳过 NI terminal（RePlAce 中的 "NI" / do-not-place 类型，如 IO filler、tapless 等占位）
    for (Module *curTerminal : db->dbTerminals)
    {
        if (curTerminal->isNI) //!
        {
            continue;
        }
        CRect terminal;
        terminal.ll = curTerminal->getLL_2D();
        terminal.ur = curTerminal->getUR_2D();

        for (SiteRow curRow : db->dbSiteRows)
        {
            CRect placementRow;
            placementRow.ll = curRow.getLL_2D();
            placementRow.ur = curRow.getUR_2D();
            totalOverLapArea += getOverlapArea_2D(terminal, placementRow);
        }
    }

    whitespaceArea = db->totalRowArea - totalOverLapArea;

    ////////////////////////////////////////////////////////////////
    // calculate node area
    ////////////////////////////////////////////////////////////////

    float nodeAreaScaled = 0; // node = std cells + movable macros
    float stdcellArea = 0;
    float macroArea = 0;

    for (Module *curNode : db->dbNodes)
    {
        assert(curNode->getArea() > 0);
        if (curNode->isMacro)
        {
            macroArea += curNode->getArea();
        }
        else
        {
            stdcellArea += curNode->getArea();
        }
    }

    ePlaceStdCellArea = stdcellArea;
    ePlaceMacroArea = macroArea;

    //! Am：式(13) 里的单元面积项。macro 面积乘 targetDensity 做折算，
    //! 因为 macro 通常不允许被塞满到 targetDensity，等效占用面积更小
    nodeAreaScaled = stdcellArea + macroArea * targetDensity; // see ePlace paper equation (13), Am is nodeArea here. or see RePlAce code opt.cpp line 86, total_modu_area equals nodeAreaScaled here
    //??? macro area should *= target density when calculating Am in(13)? see RePlAce code opt.cpp line 86 But terminal area wasn't *= target density when calculating Aws??? implement as this for now

    ////////////////////////////////////////////////////////////////
    // calculate filler area
    ////////////////////////////////////////////////////////////////

    float totalFillerArea = 0;                                         // area of all fillers, how to calculate: see ePlace paper (13)
    totalFillerArea = whitespaceArea * targetDensity - nodeAreaScaled; //! see ePlace paper (13)

    int nodeCount = db->dbNodes.size();

    vector<float> nodeArea; // sort node according to area with this because we don't want to sort dbNodes
    nodeArea.resize(nodeCount);

    for (int i = 0; i < nodeCount; i++)
    {
        nodeArea[i] = db->dbNodes[i]->getArea();
        // use this one below if we want no problems:
        // nodeArea[i]=db->dbNodes[i]->calcArea();
    }

    sort(nodeArea.begin(), nodeArea.end()); //! sort by area

    //! 取中间 80%（去掉最小 5% 和最大 5%）的平均面积作为单个 filler 的面积，
    //! 目的是避开极端值，让 filler 尺寸贴近「典型单元」大小
    float avg80TotalArea = 0;
    float avg80NodeArea = 0;
    int minIdx = (int)(0.05 * (float)nodeCount); //! for calculating average area of the middle 80% of all nodes(cells and macros)
    int maxIdx = (int)(0.95 * (float)nodeCount);

    for (int i = minIdx; i < maxIdx; i++)
    {
        avg80TotalArea += nodeArea[i];
    }

    avg80NodeArea = avg80TotalArea / ((float)(maxIdx - minIdx));

    //! filler 高度取行高，宽度由面积反推，保证 filler 和标准单元「同高」，能自然填满布局行
    float fillerArea = avg80NodeArea; //! use average area of the middle 90% of all nodes as filler area! see ePlace paper, filler insertion
    float fillerHeight = db->commonRowHeight;
    float fillerWidth = float_div(fillerArea, fillerHeight);

    ////////////////////////////////////////////////////////////////
    // add fillers and set filler locations randomly
    ////////////////////////////////////////////////////////////////

    int fillerCount = (int)(totalFillerArea / fillerArea + 0.5); //!

    ePlaceFillers.resize(fillerCount);

    float leftMost;
    float rightMost;

    //! filler 的 idx 从 nodeCount 开始往后排，保证 ePlaceNodesAndFillers 里 idx 连续
    for (int i = 0; i < fillerCount; i++)
    {
        string name = "f" + to_string(i);
        Module *curFiller = new Module(i + nodeCount, name, fillerWidth, fillerHeight, false, false);
        curFiller->isFiller = true; //!
        ePlaceFillers[i] = curFiller;

        db->setModuleLocation_2D_random(curFiller); //! 初始位置在 core 区域内随机撒点
    }

    //! 三个「参与优化的单元集合」，对应三个放置阶段：
    //!   ePlaceNodesAndFillers = 所有可动单元(含 macro) + filler  → mGP
    //!   ePlaceFillers         = 仅 filler                        → FILLERONLY
    //!   ePlaceCellsAndFillers = std cell + filler（不含 macro）   → cGP
    //! 这些 vector 里元素的先后顺序必须与各 gradient vector 的下标严格对应
    ePlaceNodesAndFillers = db->dbNodes;
    ePlaceNodesAndFillers.insert(ePlaceNodesAndFillers.end(), ePlaceFillers.begin(), ePlaceFillers.end()); //! fillers are stored after nodes in the vector

    for (Module *curCellOrFiller : ePlaceNodesAndFillers)
    {
        if (curCellOrFiller->isFiller)
        {
            ePlaceCellsAndFillers.push_back(curCellOrFiller);
        }
        else if (!curCellOrFiller->isMacro)
        {
            ePlaceCellsAndFillers.push_back(curCellOrFiller);
        }
    }
    // macro density scaling in density computation: RePlace bin.cpp line 1853
    // terminal(fixed macro) density scaling in density computation?: RePlace bin.cpp line 480

    // density scaling when calculating Aws?
}

// ----------------------------------------------------------------------------
// bin 网格初始化（划分整个 coreRegion）
//
// 1. 先确定 bin 数量：理想 bin 面积 = 平均单元面积 / targetDensity，
//    使「一个 bin 大致装一个单元」，再向上取到最近的 2 的幂（4×4 … 1024×1024），
//    这样 FFT 可以直接用基-2 算法，效率最高。
// 2. 逐个 bin 计算 ll/ur/center/area（坐标要加上 coreRegion.ll 偏移）。
// 3. terminalDensity：terminal 固定不动，重叠面积只算一次。
// 4. baseDensity：bin 内「不在任何布局行上」的不可放置区域（对应 RePlAce 的 virt_area），
//    同样固定不变，只算一次。
// ----------------------------------------------------------------------------
void EPlacer_2D::binInitialization()
{
    segmentFaultCP("binInit");
    //! Bins are allocated on coreRegion! RePlAce code bin.cpp line 337 and bookshelfIO.cpp, also see ePlace paper
    ////////////////////////////////////////////////////////////////
    // calculate bin dimension and size
    ////////////////////////////////////////////////////////////////
    //! this code for bin dimension calculation is copied directly from RePlAce
    segmentFaultCP("binDim");
    int nodeCount = db->dbNodes.size();
    float nodeArea = (ePlaceStdCellArea + ePlaceMacroArea);

    float coreRegionWidth = db->coreRegion.getWidth();
    float coreRegionHeight = db->coreRegion.getHeight();
    float coreRegionArea = float_mul(coreRegionWidth, coreRegionHeight);

    float averageNodeArea = 1.0 * nodeArea / nodeCount;
    float idealBinArea = averageNodeArea / targetDensity;

    int idealBinCount = INT_CONVERT(coreRegionArea / idealBinArea);

    bool isUpdate = false;
    // bin dimension upper bound: 1024 rather than 2048
    // for (int i = 1; i <= 10; i++)
    // { //! 4*4,8*8,16*16,32*32..., 2048*2048
    //! 在 4×4, 8×8, …, 1024×1024 里挑第一个「不超过 idealBinCount 的最大档」。
    //! 注意 2 << i 即 2^(i+1)：i=1 → 4，i=9 → 1024（故上限写成 i<10）
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

    cout << BLUE << "Bin dimension: " << binDimension << "\ncoreRegion width: " << coreRegionWidth << "\ncoreRegion height: " << coreRegionHeight << RESET << endl;

    binStep.x = float_div(coreRegionWidth, binDimension.x);
    binStep.y = float_div(coreRegionHeight, binDimension.y);

    cout << BLUE << "Bin step: " << binStep << RESET << endl;
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

    for (int i = 0; i < binDimension.x; i++)
    {
        bins[i].resize(binDimension.y);
        for (int j = 0; j < binDimension.y; j++)
        {
            bins[i][j] = new Bin_2D();
            // cout<<"adding bin "<<i<<","<<j<<endl;
            //!!! +db->coreRegion.ll.x to get coordinates!!
            bins[i][j]->ll.x = i * binStep.x + db->coreRegion.ll.x;
            bins[i][j]->ll.y = j * binStep.y + db->coreRegion.ll.y;

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
    cout << "Bin add time: " << addBinTime << endl;
    ////////////////////////////////////////////////////////////////
    // terminal density calculation, calculate here because they are terminals and only needed to be considered once
    ////////////////////////////////////////////////////////////////
    segmentFaultCP("terminalDensity");
    VECTOR_2D_INT binStartIdx;
    VECTOR_2D_INT binEndIdx;

    time_start(&terminalDensityTime);

    for (Module *curTerminal : db->dbTerminals)
    {
        //! only consider terminals inside the coreRegion
        //! assume no terminal would have a part inside the coreRegion and a part outside the coreRegion
        if (curTerminal->getLL_2D().x < db->coreRegion.ll.x || curTerminal->getLL_2D().x + curTerminal->getWidth() > db->coreRegion.ur.x)
        {
            continue;
        }
        if (curTerminal->getLL_2D().y < db->coreRegion.ll.y || curTerminal->getLL_2D().y + curTerminal->getHeight() > db->coreRegion.ur.y)
        {
            continue;
        }

        //! 把 terminal 的包围盒映射到 bin 下标区间 [binStartIdx, binEndIdx]，
        //! 之后只在这段区间内累加重叠面积（INT_DOWN 即向下取整）
        binStartIdx.x = INT_DOWN((curTerminal->getLL_2D().x - db->coreRegion.ll.x) / binStep.x);
        binEndIdx.x = INT_DOWN((curTerminal->getUR_2D().x - db->coreRegion.ll.x) / binStep.x);

        binStartIdx.y = INT_DOWN((curTerminal->getLL_2D().y - db->coreRegion.ll.y) / binStep.y);
        binEndIdx.y = INT_DOWN((curTerminal->getUR_2D().y - db->coreRegion.ll.y) / binStep.y);

        // cout << curTerminal->getLL_2D() << curTerminal->getUR_2D() << endl;
        // cout << binStep << " " << db->coreRegion.ll << endl;
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
                //! terminal 也按 targetDensity 折算，和 nodeDensity 的口径保持一致
                //! beware: density scaling!
                bins[i][j]->terminalDensity += targetDensity * getOverlapArea_2D(bins[i][j]->ll, bins[i][j]->ur, curTerminal->getLL_2D(), curTerminal->getUR_2D());
            }
        }
    }

    time_end(&terminalDensityTime);
    cout << "Terminal density time: " << terminalDensityTime << endl;
    ////////////////////////////////////////////////////////////////
    // base density calculation, also only needed to be considered once
    ////////////////////////////////////////////////////////////////
    segmentFaultCP("baseDensity");

    time_start(&baseDensityTime);

    for (int i = 0; i < binDimension.x; i++)
    {
        for (int j = 0; j < binDimension.y; j++)
        {
            //! bin 与所有布局行的重叠面积 = 该 bin 里「可放置」的面积；
            //! 差额 (bin面积 − 可放置面积) 就是不可放置区域，乘 targetDensity 记入 baseDensity。
            //! 完全落在布局行内的 bin 其 baseDensity 为 0。
            //! 注意这里是 O(binCount × rowCount) 的双层循环，bin 数多时是初始化热点。
            float curBinAvailableArea = 0; // overlap area between current bin and placement rows
            for (SiteRow curRow : db->dbSiteRows)
            {
                curBinAvailableArea += getOverlapArea_2D(bins[i][j]->ll, bins[i][j]->ur, curRow.getLL_2D(), curRow.getUR_2D());
            }
            // debugOutput("Bin area", bins[i][j]->area);
            // debugOutput("Available area", curBinAvailableArea);
            if (float_equal(bins[i][j]->area, curBinAvailableArea))
            {
                bins[i][j]->baseDensity = 0;
            }
            else
            {
                bins[i][j]->baseDensity = targetDensity * (bins[i][j]->area - curBinAvailableArea); //! follow RePlAce bin.cpp line 433
            }
        }
    }
    time_end(&baseDensityTime);
    cout << "Base density time: " << baseDensityTime << endl;
}

//! 分配各梯度容器。注意容量差异：
//!   wirelengthGradient 只有 dbNodes 个（filler 不连线，线长梯度恒为 0，不占空间）
//!   densityGradient / totalGradient 是 nodes + fillers 个
//!   cGPGradient 是 std cell + fillers 个，fillerGradient 只有 fillers 个
void EPlacer_2D::gradientVectorInitialization()
{
    wirelengthGradient.resize(db->dbNodes.size()); // fillers has no wirelength gradient
    densityGradient.resize(ePlaceNodesAndFillers.size());
    totalGradient.resize(ePlaceNodesAndFillers.size());

    cGPGradient.resize(ePlaceCellsAndFillers.size());
    fillerGradient.resize(ePlaceFillers.size());
}

// ----------------------------------------------------------------------------
// 更新全局溢出率 τ（globalDensityOverflow）
//
//   τ = Σ_bin max(0, ρ_bin − targetDensity) · binArea  /  Am
// 其中 ρ_bin = (nodeDensity + terminalDensity + baseDensity) / binArea。
//
// 注意：这里**故意不含 fillerDensity**——filler 是用来填白的，若计入会把 τ 人为抬高，
// 使算法误判为「还很拥塞」。τ 是 Nesterov 优化的收敛判据（见 nesterov.hpp stop_condition），
// 也是下面 γ 计算的输入。
// ----------------------------------------------------------------------------
void EPlacer_2D::densityOverflowUpdate()
{
    segmentFaultCP("densityOverflow");
    float globalOverflowArea = 0;
    float nodeAreaScaled = ePlaceStdCellArea + ePlaceMacroArea * targetDensity;
    float invertedBinArea = 1.0 / (binStep.x * binStep.y); // 1/bin area
    for (int i = 0; i < binDimension.x; i++)
    {
        for (int j = 0; j < binDimension.y; j++)
        {
            // nodeDensity+terminalDensity+baseDensity because filler density are not included
            globalOverflowArea += max(float(0.0), (bins[i][j]->nodeDensity + bins[i][j]->terminalDensity + bins[i][j]->baseDensity) * invertedBinArea - targetDensity) * bins[i][j]->area;
        }
    }
    globalDensityOverflow = globalOverflowArea / nodeAreaScaled; // see RePlAce code bin.cpp line 1183, opt.cpp line 86. And nodeAreaScaled in fillerInitialization() in this file
}

// ----------------------------------------------------------------------------
// 线长梯度更新
//
// Step1 先算平滑系数 γ（这里实际存的是 1/γ，跟随 RePlAce 的习惯）：
//   1/γ = (1 / (8·wb)) · 10^(−(k·τ + b))，k = 20/9，b = −11/9  （ePlace 论文式38）
//   τ 大（还很拥塞）→ γ 小 → 线长模型更接近线性、平滑更强，先顾着铺开；
//   τ 小（快收敛）  → γ 大 → 线长模型更接近真实 HPWL，开始精细优化。
//   代码里对 τ>1 和 τ<0.1 两个极端做了截断，避免指数项数值溢出/过小。
// Step2 用平滑线长模型（默认 WA，加 -LSE 参数则切到 LSE）逐 pin 求梯度并累加到单元上。
//   两种模型都依赖每个 net 的 X/Y 边界 pin，所以每一步都要先调 calcNetBoundPins() 刷新。
// ----------------------------------------------------------------------------
void EPlacer_2D::wirelengthGradientUpdate()
{
    ////////////////////////////////////////////////////////////////
    //! Step1: calculate gamma, see ePlace paper equation 38, RePlace code wlen.cpp line 141
    ////////////////////////////////////////////////////////////////
    // first, calculate tau(density overflow)
    densityOverflowUpdate();
    segmentFaultCP("wireLengthGradient");
    //! now calculate gamma with the updated tau, here we actually calculate 1/gamma for furthurer calculation
    VECTOR_2D baseWirelengthCoef;
    baseWirelengthCoef.x = baseWirelengthCoef.y = 0.125 /*0.5*/; // 0.125=1/8.0, 8.0:see ePlace paper equation 38. Notice that baseWirelngthCoef
                                                                 // is wcof00_org in RePlace code wlen.cpp,
                                                                 // and is tuned according to input benchmark in RePlAce main.cpp
    baseWirelengthCoef.x = baseWirelengthCoef.x / binStep.x;     // binStep: wb in ePlace paper equation 38, 1/8/wb=1/8.0wb
    baseWirelengthCoef.y = baseWirelengthCoef.y / binStep.y;

    if (globalDensityOverflow > 1.0)
    {
        baseWirelengthCoef.x *= 0.1;
        baseWirelengthCoef.y *= 0.1;
    }
    else if (globalDensityOverflow < 0.1)
    {
        baseWirelengthCoef.x *= 10.0;
        baseWirelengthCoef.y *= 10.0;
    }
    else
    {
        float temp;
        temp = 1.0 / pow(10.0, (globalDensityOverflow - 0.1) * 20 / 9.0 - 1.0); //! see eplace paper equation 38
        baseWirelengthCoef.x *= temp;                                           //!(1/8.0wb)*(1/10^(k*tau+b)), where k=20/9 and b=-11/9
        baseWirelengthCoef.y *= temp;
    }

    invertedGamma = baseWirelengthCoef;

    ////////////////////////////////////////////////////////////////
    //! Step2: calculate wirelength density for each nodes (not filler nodes)
    ////////////////////////////////////////////////////////////////
    // When using weighted-average wirelength model we would need X/Y/Z max and min in a net,
    // so update X/Y/Z max and min in all nets first, see ntuplace3D paper page 6: Stable Weighted-Average Wirelength Model
    // Also the numerators and denominators are pre-calculated for all nets for further use
    //! 先刷新 τ，因为 γ 是 τ 的函数；同时 calcNetBoundPins 会顺带更新每个 net 的边界 pin
    double HPWL = db->calcNetBoundPins();

    //! HPWL 在这里只是算出来备用（LSE 分支下未使用），真正用的是下面的梯度
    if (gArg.CheckExist("LSE"))
    {
        double LSE = db->calcLSE_Wirelength_2D(invertedGamma);
        int index = 0;
        for (Module *curNode : db->dbNodes) // use ePlaceNodesAndFillers?
        {
            assert(curNode->idx == index);
            wirelengthGradient[index].SetZero(); //! clear before updating
            for (Pin *curPin : curNode->modulePins)
            {
                VECTOR_2D gradient;
                gradient = curPin->net->getWirelengthGradientLSE_2D(invertedGamma, curPin);
                wirelengthGradient[index].x += gradient.x;
                wirelengthGradient[index].y += gradient.y;
                // get the wirelength gradient of this pin
            }
            index++;
        }
    }
    else
    {
        double WA = db->calcWA_Wirelength_2D(invertedGamma);
        int index = 0;
        for (Module *curNode : db->dbNodes) // use ePlaceNodesAndFillers?
        {
            assert(curNode->idx == index);
            wirelengthGradient[index].SetZero(); //! clear before updating
            for (Pin *curPin : curNode->modulePins)
            {
                VECTOR_2D gradient;
                gradient = curPin->net->getWirelengthGradientWA_2D(invertedGamma, curPin);
                wirelengthGradient[index].x += gradient.x;
                wirelengthGradient[index].y += gradient.y;
                // get the wirelength gradient of this pin
            }
            index++;
        }
    }
}

// ----------------------------------------------------------------------------
// 密度梯度更新（ePlace 的精髓：静电场）
//
// Step1 求电场：把每个 bin 的电荷密度 ρ = (node + base + filler + terminal) / binArea
//   送进 FFT，解泊松方程 ∇²φ = −ρ，得到每个 bin 的电势 φ 和电场 E = −∇φ。
//   （这里**包含** fillerDensity——filler 必须参与产生斥力，否则填白无从谈起。）
//   全局的 FFT 让这一步复杂度是 O(N log N) 而非 O(N²)。
// Step2 求梯度：对每个单元，把它覆盖到的 bin 上的 E 按「与 bin 的重叠面积」加权求和，
//   即 ∂D/∂x = Σ_bins overlapArea · E_bin。
//   local smooth：比 bin 还小的单元会被「撑大」到一整个 bin 宽再算重叠，
//   同时按 (实际尺寸 / bin尺寸) 缩放面积，避免小单元只在单个 bin 上产生脉冲式梯度。
// ----------------------------------------------------------------------------
void EPlacer_2D::densityGradientUpdate()
{
    ////////////////////////////////////////////////////////////////
    //! Step1: obtain electric field(e) through FFT
    ////////////////////////////////////////////////////////////////
    segmentFaultCP("densityGradient");
    replace::FFT_2D fft(binDimension.x, binDimension.y, binStep.x, binStep.y);
    float invertedBinArea = 1.0 / (binStep.x * binStep.y);
    for (int i = 0; i < binDimension.x; i++)
    {
        for (int j = 0; j < binDimension.y; j++)
        {
            float eDensity = bins[i][j]->nodeDensity + bins[i][j]->baseDensity + bins[i][j]->fillerDensity + bins[i][j]->terminalDensity; // consider filler area(density) here
            eDensity *= invertedBinArea;
            fft.updateDensity(i, j, eDensity);
        }
    }
    fft.doFFT();
    for (int i = 0; i < binDimension.x; i++)
    {
        for (int j = 0; j < binDimension.y; j++)
        {
            auto eForcePair = fft.getElectroForce(i, j);
            bins[i][j]->E.x = eForcePair.first;
            bins[i][j]->E.y = eForcePair.second;
            // std::cout<<"bin eforce x and y "<<Bins[i][j].Eforce.x<<" "<<Bins[i][j].Eforce.y<<std::endl;
            float electroPhi = fft.getElectroPhi(i, j);
            bins[i][j]->phi = electroPhi;
        }
        // sumPhi_ += electroPhi*static_cast<float>(bin->nonPlaceArea()+bin->instPlacedArea()+bin->fillerArea());
    }

    ////////////////////////////////////////////////////////////////
    //! Step2: calculate density(potential) gradient, see ePlace paper equation 16, RePlace charge.cpp line 451
    ////////////////////////////////////////////////////////////////
    int nodeCount = db->dbNodes.size();
    int index = 0;
    for (Module *curNode : ePlaceNodesAndFillers)
    {
        assert(index == curNode->idx);
        //! clear before updating
        densityGradient[index].SetZero();

        VECTOR_2D localSmoothLengthScale; // see ePlace paper page 15 or RePlace opt.cpp line 1460
        localSmoothLengthScale.x = 1.0;
        localSmoothLengthScale.y = 1.0;

        CRect rectForCurNode;
        rectForCurNode.ll = curNode->getLL_2D();
        rectForCurNode.ur = curNode->getUR_2D();

        //! beware: local smooth on x and y dimension, also needed here!
        //! binStart and binEnd should be calculated with inflated cell width and height, see replace charge.cpp line 408
        POS_3D cellCenter = curNode->getCenter();

        // if (!curNode->isMacro)
        // {
        //! 单元比 bin 小 → 把它在 X 方向撑到一整个 bin 宽，并记录缩放系数，
        //! 这样重叠面积算的是「撑大后」的，再乘 scale 还原成等效面积
        //! local smooth: not only for std cells because there may be small macros, like in MMS bigblue3
        if (float_less(curNode->getWidth(), binStep.x))
        {
            localSmoothLengthScale.x = curNode->getWidth() / binStep.x;
            rectForCurNode.ll.x = cellCenter.x - 0.5 * binStep.x;
            rectForCurNode.ur.x = cellCenter.x + 0.5 * binStep.x;
        }
        if (float_less(curNode->getHeight(), binStep.y))
        {
            localSmoothLengthScale.y = curNode->getHeight() / binStep.y;
            rectForCurNode.ll.y = cellCenter.y - 0.5 * binStep.y;
            rectForCurNode.ur.y = cellCenter.y + 0.5 * binStep.y;
        }
        // }

        VECTOR_2D_INT binStartIdx; // binStartIdx: index the index of the first bin that has overlap with a cell on X/Y direction
        VECTOR_2D_INT binEndIdx;
        binStartIdx.x = INT_DOWN((rectForCurNode.ll.x - db->coreRegion.ll.x) / binStep.x);
        binEndIdx.x = INT_DOWN((rectForCurNode.ur.x - db->coreRegion.ll.x) / binStep.x);

        binStartIdx.y = INT_DOWN((rectForCurNode.ll.y - db->coreRegion.ll.y) / binStep.y);
        binEndIdx.y = INT_DOWN((rectForCurNode.ur.y - db->coreRegion.ll.y) / binStep.y);

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

        //! beware: local smooth
        for (int i = binStartIdx.x; i <= binEndIdx.x; i++)
        {
            for (int j = binStartIdx.y; j <= binEndIdx.y; j++)
            {
                float overlapArea = localSmoothLengthScale.x * localSmoothLengthScale.y * getOverlapArea_2D(bins[i][j]->ll, bins[i][j]->ur, rectForCurNode.ll, rectForCurNode.ur);
                //! ????????????????????? watch out: do we need macro density scaling here? it seems RePlAce didn't do that
                densityGradient[index].x += overlapArea * bins[i][j]->E.x;
                densityGradient[index].y += overlapArea * bins[i][j]->E.y;
            }
        }

        index++;
    }
}

// ----------------------------------------------------------------------------
// 合成总梯度 —— 优化器每次迭代实际取用的就是这里的结果
//
//   可动单元：∇f = λ·∇D − ∇WL      （注意 ∇WL 取负号：要最小化线长，即往线长下降方向走）
//   filler  ：∇f = λ·∇D            （不连线，无线长项）
//
// 再乘一个 preconditioner = 1 / max(1, pin数 + λ·面积)，
// 让「连线多 / 面积大」的单元步长小一些，改善一阶优化的收敛性。
// 同时按阶段把结果分别填进 fillerGradient / cGPGradient，供 getGradient() 直接返回。
// ----------------------------------------------------------------------------
void EPlacer_2D::totalGradientUpdate()
{

    binNodeDensityUpdate();
    densityGradientUpdate();
    wirelengthGradientUpdate();

    segmentFaultCP("totalGradient");
    int index = 0;
    int cGPindex = 0;
    int fillerOnlyIndex = 0;
    for (Module *curNodeOrFiller : ePlaceNodesAndFillers) // branch in this for can be eliminated
    {
        totalGradient[index].SetZero();
        assert(index == curNodeOrFiller->idx); //! 下标必须和 module->idx 对齐，否则梯度会张冠李戴

        //! precondition
        float connectedNetNum = curNodeOrFiller->modulePins.size();
        // printf("Connected net num:%d, Pin count:%d\n",connectedNetNum,placer->ePlaceNodesAndFillers[idx]->modulePins.size());
        // assert(connectedNetNum == placer->ePlaceNodesAndFillers[idx]->modulePins.size());
        float charge = curNodeOrFiller->getArea();
        float preconditioner = 1 / max(1.0f, (connectedNetNum + lambda * charge));
        // preconditionedGradient[idx].x = preconditioner * placer->totalGradient[idx].x;
        // preconditionedGradient[idx].y = preconditioner * placer->totalGradient[idx].y;
        // calculate -gradient here
        if (curNodeOrFiller->isFiller)
        {
            // wirelength gradient of fillers should == 0
            //! filler 同时进 fillerGradient（FILLERONLY 阶段）和 cGPGradient（cGP 阶段）
            totalGradient[index].x = preconditioner * lambda * densityGradient[index].x;
            totalGradient[index].y = preconditioner * lambda * densityGradient[index].y;

            fillerGradient[fillerOnlyIndex] = totalGradient[index];
            fillerOnlyIndex++;

            cGPGradient[cGPindex] = totalGradient[index];
            cGPindex++;
        }
        else
        {
            totalGradient[index].x = preconditioner * (lambda * densityGradient[index].x - wirelengthGradient[index].x);
            totalGradient[index].y = preconditioner * (lambda * densityGradient[index].y - wirelengthGradient[index].y);
            if (!curNodeOrFiller->isMacro)
            {
                cGPGradient[cGPindex] = totalGradient[index];
                cGPindex++;
            }
        }

        index++;
    }
}

//! 按当前阶段返回优化器需要的梯度向量（三者的长度和顺序各不相同，见 gradientVectorInitialization）
//! ! 注意：placementStage 既不是这三个值时会走到函数末尾而没有 return，属未定义行为
vector<VECTOR_3D> EPlacer_2D::getGradient()
{
    if (placementStage == mGP)
    {
        return totalGradient;
    }
    else if (placementStage == FILLERONLY)
    {
        return fillerGradient;
    }
    else if (placementStage == cGP)
    {
        return cGPGradient;
    }
}

//! 返回参与当前阶段优化的单元的中心坐标，顺序与 getGradient() 一一对应
vector<VECTOR_3D> EPlacer_2D::getPosition()
{
    if (placementStage == mGP)
    {
        return getModulePositions(ePlaceNodesAndFillers);
    }
    else if (placementStage == FILLERONLY)
    {
        return getModulePositions(ePlaceFillers);
    }
    else if (placementStage == cGP)
    {
        return getModulePositions(ePlaceCellsAndFillers);
    }
}

//! 把优化器算出的新坐标写回数据库。setModuleCenter_2D 只改中心、不改尺寸，
//! 并会同步刷新该单元上所有 pin 的绝对坐标，供下一轮线长梯度使用
void EPlacer_2D::setPosition(vector<VECTOR_3D> modulePositions)
{
    int moduleCount;
    if (placementStage == mGP)
    {
        moduleCount = ePlaceNodesAndFillers.size();
        for (int i = 0; i < moduleCount; i++)
        {
            db->setModuleCenter_2D(ePlaceNodesAndFillers[i], modulePositions[i]);
        }
    }
    else if (placementStage == FILLERONLY)
    {
        moduleCount = ePlaceFillers.size();
        for (int i = 0; i < moduleCount; i++)
        {
            db->setModuleCenter_2D(ePlaceFillers[i], modulePositions[i]);
        }
    }
    else if (placementStage == cGP)
    {
        moduleCount = ePlaceCellsAndFillers.size();
        for (int i = 0; i < moduleCount; i++)
        {
            db->setModuleCenter_2D(ePlaceCellsAndFillers[i], modulePositions[i]);
        }
    }
}

// ----------------------------------------------------------------------------
// λ 初始化（ePlace 论文式 35）
//   λ₀ = Σ|∇WL| / Σ|∇D|
// 即让线长项和密度项的梯度量级相当，两项在优化中权重均衡。
// 分子只累加真实单元的线长梯度，分母要额外加上 filler 的密度梯度（下标 nodeCount 之后）。
// ----------------------------------------------------------------------------
void EPlacer_2D::penaltyFactorInitilization()
{
    lastHPWL = db->calcHPWL();
    float denominator = 0;
    float numerator = 0;

    int nodeCount = wirelengthGradient.size();
    int nodeAndFillerCount = densityGradient.size();

    for (int i = 0; i < nodeCount; i++)
    {
        numerator += fabs(wirelengthGradient[i].x);
        numerator += fabs(wirelengthGradient[i].y);
        denominator += fabs(densityGradient[i].x);
        denominator += fabs(densityGradient[i].y);
    }
    for (int i = nodeCount; i < nodeAndFillerCount; i++) //! filler 只有密度梯度，补到分母
    {
        denominator += fabs(densityGradient[i].x);
        denominator += fabs(densityGradient[i].y);
    }

    lambda = float_div(numerator, denominator);
}

// ----------------------------------------------------------------------------
// λ 自适应更新（ePlace 论文式 36）
//    multiplier = base^( −ΔHPWL / DELTA_HPWL_REF + 1 ) ，再夹到 [0.95, 1.05]
//    λ ← λ · multiplier
// 直观含义：
//    ΔHPWL > 0（线长变差，说明密度力推得太猛）→ 指数 < 1 → λ 缩小，放松密度约束；
//    ΔHPWL < 0（线长变好）                    → 直接用上界 1.05 → λ 增大，继续铺开。
// 这是让「铺开」和「缩短线长」两个目标逐步达到平衡的关键机制。
// ----------------------------------------------------------------------------
void EPlacer_2D::updatePenaltyFactor()
{
    // printf("penalty factor = %.10f\n", lambda);
    float curHPWL = db->calcHPWL();
    float multiplier;
    double deltaHPWL = curHPWL - lastHPWL;
    if ((deltaHPWL) < 0.0) //?? what if (curHPWL - lastHPWL)<0????? never considered before 2024.5.19
    {
        //! 线长改善了 → 直接取上界放大 λ，加大密度力继续铺开
        // cout << "multiplier is: " << pow(PENALTY_MULTIPLIER_BASE, (-(deltaHPWL) / DELTA_HPWL_REF + 1.0)) << " when < 0" << endl;
        multiplier = PENALTY_MULTIPLIER_UPPERBOUND;
    }
    else
    {
        multiplier = pow(PENALTY_MULTIPLIER_BASE, (-(deltaHPWL) / DELTA_HPWL_REF + 1.0)); // see ePlace-3D code opt.cpp line 1523
    }

    if (float_greater(multiplier, PENALTY_MULTIPLIER_UPPERBOUND))
    {
        multiplier = PENALTY_MULTIPLIER_UPPERBOUND;
    }
    if (float_less(multiplier, PENALTY_MULTIPLIER_LOWERBOUND))
    {
        multiplier = PENALTY_MULTIPLIER_LOWERBOUND;
    }
    lambda *= multiplier;
    //! 下面这段是曾经尝试过的 λ 硬上下限裁剪，现已注释掉——
    //! 目前完全靠 multiplier 的 [0.95, 1.05] 夹取来约束 λ 的变化幅度
    // if (penaltyFactor < 0.00001)
    // {
    //     penaltyFactor = 0.00001;
    // }
    // else if (penaltyFactor > 10.0)
    // {
    //     penaltyFactor = 10.0;
    // }
    lastHPWL = curHPWL;
}

//! mGP → FILLERONLY：把所有 filler 重新随机撒一遍，
//! 目的是打破 mGP 收敛后残留的局部拥塞，给后续 cGP 一个更好的初始分布
void EPlacer_2D::switch2FillerOnly()
{
    for (Module *curFiller : ePlaceFillers)
    {
        db->setModuleLocation_2D_random(curFiller);
    }
    placementStage = FILLERONLY;
}

//! mGP → cGP：macro 冻结，只剩 std cell + filler 继续优化
void EPlacer_2D::switch2cGP()
{
    //! 1. update lambda：mGP 跑得越久 λ 被抬得越高，这里按 1.1^(迭代数×0.1) 回调，
    //!    因为 cGP 阶段只剩单元，需要的密度力比 mGP 小
    lambda = lambda / pow(1.1, mGPIterationCount * 0.1);
    //! 2. update placement stage
    placementStage = cGP;
}

//! 打印当前迭代的关键指标（溢出率 τ、惩罚因子 λ、线长）
void EPlacer_2D::showInfo()
{
    cout << "Overflow: " << globalDensityOverflow << endl;
    cout << "penalty factor: " << fixed << lambda << endl;
    //! 注意打印的是 lastHPWL，即上一次 updatePenaltyFactor 时记录的线长，不是实时值
    cout << "HPWL: " << lastHPWL << endl
         << endl;
}

// 下面是一整段被注释掉的 plotCurrentPlacement：用 CImg 把当前布局画成 bmp（terminal 蓝、
// macro 橙、std cell 红、filler 绿）。绘图功能已迁到 Plot/plot.cpp（PLOTTING::plotCurrentPlacement），
// 这里保留作历史备份。
// void EPlacer_2D::plotCurrentPlacement(string imageName)
// {
//     string plotPath;
//     if (!gArg.GetString("plotPath", &plotPath))
//     {
//         plotPath = "./";
//     }

//     float chipRegionWidth = db->chipRegion.ur.x - db->chipRegion.ll.x;
//     float chipRegionHeight = db->chipRegion.ur.y - db->chipRegion.ll.y;

//     int minImgaeLength = 1000;

//     int imageHeight;
//     int imageWidth;

//     float opacity = 0.7;
//     int xMargin = 30, yMargin = 30;

//     if (chipRegionWidth < chipRegionHeight)
//     {
//         imageHeight = 1.0 * chipRegionHeight / (chipRegionWidth / minImgaeLength);
//         imageWidth = minImgaeLength;
//     }
//     else
//     {
//         imageWidth = 1.0 * chipRegionWidth / (chipRegionHeight / minImgaeLength);
//         imageHeight = minImgaeLength;
//     }

//     CImg<unsigned char> img(imageWidth + 2 * xMargin, imageHeight + 2 * yMargin, 1, 3, 255);

//     float unitX = imageWidth / chipRegionWidth,
//           unitY = imageHeight / chipRegionHeight;

//     for (Module *curTerminal : db->dbTerminals)
//     {
//         assert(curTerminal);
//         // ignore pin's location
//         if (curTerminal->isNI)
//         {
//             continue;
//         }
//         int x1 = getX(db->chipRegion.ll.x, curTerminal->getLL_2D().x, unitX) + xMargin;
//         int x2 = getX(db->chipRegion.ll.x, curTerminal->getUR_2D().x, unitX) + xMargin;
//         int y1 = getY(chipRegionHeight, db->chipRegion.ll.y, curTerminal->getLL_2D().y, unitY) + yMargin;
//         int y2 = getY(chipRegionHeight, db->chipRegion.ll.y, curTerminal->getUR_2D().y, unitY) + yMargin;
//         img.draw_rectangle(x1, y1, x2, y2, Blue, opacity);
//     }

//     for (Module *curNode : db->dbNodes)
//     {
//         assert(curNode);
//         int x1 = getX(db->chipRegion.ll.x, curNode->getLL_2D().x, unitX) + xMargin;
//         int x2 = getX(db->chipRegion.ll.x, curNode->getUR_2D().x, unitX) + xMargin;
//         int y1 = getY(chipRegionHeight, db->chipRegion.ll.y, curNode->getLL_2D().y, unitY) + yMargin;
//         int y2 = getY(chipRegionHeight, db->chipRegion.ll.y, curNode->getUR_2D().y, unitY) + yMargin;
//         if (curNode->isMacro)
//         {
//             img.draw_rectangle(x1, y1, x2, y2, Orange, opacity);
//         }
//         else
//         {
//             img.draw_rectangle(x1, y1, x2, y2, Red, opacity);
//         }
//     }

//     if ((gArg.CheckExist("debug") || placementStage == FILLERONLY))
//     {
//         for (Module *curNode : ePlaceFillers)
//         {
//             assert(curNode);
//             int x1 = getX(db->chipRegion.ll.x, curNode->getLL_2D().x, unitX) + xMargin;
//             int x2 = getX(db->chipRegion.ll.x, curNode->getUR_2D().x, unitX) + xMargin;
//             int y1 = getY(chipRegionHeight, db->chipRegion.ll.y, curNode->getLL_2D().y, unitY) + yMargin;
//             int y2 = getY(chipRegionHeight, db->chipRegion.ll.y, curNode->getUR_2D().y, unitY) + yMargin;
//             img.draw_rectangle(x1, y1, x2, y2, Green, opacity);
//         }
//     }

//     img.draw_text(50, 50, imageName.c_str(), Black, NULL, 1, 30);
//     img.save_bmp(string(plotPath + imageName + string(".bmp")).c_str());
//     cout << "INFO: BMP HAS BEEN SAVED: " << imageName + string(".bmp") << endl;
// }

//! 取出一组 module 的中心坐标（z 分量保留，实际 2D 布局里恒为 0）
vector<VECTOR_3D> EPlacer_2D::getModulePositions(vector<Module *> modules)
{
    int moduleCount = modules.size();
    vector<VECTOR_3D> res;
    res.resize(moduleCount);

    for (int i = 0; i < modules.size(); i++)
    {
        res[i] = modules[i]->getCenter();
    }
    return res;
}

// ----------------------------------------------------------------------------
// 把每个单元/filler 的面积摊进 bin（每次迭代都要重算，因为位置变了）
//
// 两个关键修正：
//   · local smooth：单元比 bin 还小时，把它撑大到一整个 bin 再算重叠，
//     否则小单元只在单个 bin 上留下尖峰，密度场和电场都会抖动。
//   · macro density scaling：macro 的面积乘 targetDensity 折算，
//     与 fillerInitialization 里 nodeAreaScaled 的口径保持一致。
// filler 的面积单独记进 fillerDensity，不混进 nodeDensity
//   （这样 densityOverflowUpdate 才能不计 filler 来算 τ）。
// ----------------------------------------------------------------------------
void EPlacer_2D::binNodeDensityUpdate()
{
    //!!!! clear nodeDensity for each bin before update!
    for (int i = 0; i < binDimension.x; i++)
    {
        for (int j = 0; j < binDimension.y; j++)
        {
            bins[i][j]->nodeDensity = 0;
            bins[i][j]->fillerDensity = 0;
        }
    }

    segmentFaultCP("nodeDensity");
    for (Module *curNode : ePlaceNodesAndFillers) // ePlaceNodes: nodes and filler nodes
    {
        bool localSmooth = false;         //! local smooth is applied only to std cells, is this right?
        bool macroDensityScaling = false; // density scaling, see ePlace paper
        //! 注：这两个 flag 下面只被赋值、从未被读取，实际判断走的是 isMacro / 尺寸比较，属于遗留的死变量

        VECTOR_2D localSmoothLengthScale; // see ePlace paper page 15 or RePlace opt.cpp line 1460
        localSmoothLengthScale.x = 1;
        localSmoothLengthScale.y = 1;

        CRect rectForCurNode;
        rectForCurNode.ll = curNode->getLL_2D();
        rectForCurNode.ur = curNode->getUR_2D();

        //! beware: local smooth on x and y dimension!
        //! binStart and binEnd should be calculated with inflated cell width and height, see replace bin.cpp line 1807
        //! local smooth: not only for std cells because there may be small macros, like in MMS bigblue3
        POS_3D cellCenter = curNode->getCenter();

        if (float_less(curNode->getWidth(), binStep.x))
        {
            localSmoothLengthScale.x = curNode->getWidth() / binStep.x;
            rectForCurNode.ll.x = cellCenter.x - 0.5 * binStep.x;
            rectForCurNode.ur.x = cellCenter.x + 0.5 * binStep.x;
        }
        if (float_less(curNode->getHeight(), binStep.y))
        {
            localSmoothLengthScale.y = curNode->getHeight() / binStep.y;
            rectForCurNode.ll.y = cellCenter.y - 0.5 * binStep.y;
            rectForCurNode.ur.y = cellCenter.y + 0.5 * binStep.y;
        }

        if (curNode->isMacro)
        {
            macroDensityScaling = true;
        }

        VECTOR_2D_INT binStartIdx; // binStartIdx: index the index of the first bin that has overlap with a cell on X/Y direction
        VECTOR_2D_INT binEndIdx;
        binStartIdx.x = INT_DOWN((rectForCurNode.ll.x - db->coreRegion.ll.x) / binStep.x);
        binEndIdx.x = INT_DOWN((rectForCurNode.ur.x - db->coreRegion.ll.x) / binStep.x);

        binStartIdx.y = INT_DOWN((rectForCurNode.ll.y - db->coreRegion.ll.y) / binStep.y);
        binEndIdx.y = INT_DOWN((rectForCurNode.ur.y - db->coreRegion.ll.y) / binStep.y);

        if (!(binStartIdx.x >= 0))
        {
            cout << "Module pos: " << rectForCurNode.ll << " " << db->coreRegion.ll << endl;
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

        //! beware: local smooth and density scaling!
        for (int i = binStartIdx.x; i <= binEndIdx.x; i++)
        {
            for (int j = binStartIdx.y; j <= binEndIdx.y; j++)
            {

                //! 三种口径分开累加：macro 折算 → nodeDensity；filler → fillerDensity；普通单元 → nodeDensity
                float overlapArea = getOverlapArea_2D(bins[i][j]->ll, bins[i][j]->ur, rectForCurNode.ll, rectForCurNode.ur);
                if (curNode->isMacro)
                {
                    bins[i][j]->nodeDensity += localSmoothLengthScale.x * localSmoothLengthScale.y * targetDensity * overlapArea;
                }
                else
                {
                    if (curNode->isFiller)
                    {
                        //? does filler need localSmooth?
                        bins[i][j]->fillerDensity += localSmoothLengthScale.x * localSmoothLengthScale.y * overlapArea;
                    }
                    else
                    {
                        bins[i][j]->nodeDensity += localSmoothLengthScale.x * localSmoothLengthScale.y * overlapArea;
                    }
                }
            }
        }
    }
}
