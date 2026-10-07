#include <random>

#include "placedb.h"
#include "global.h"

// ============================================================================
// placedb.cpp —— PlaceDB 的实现（类声明与分组见 placedb.h 头部的模块总览）
//
// 本文件按 placedb.h 的五组函数组织：
//   1. 建库与区域初始化 ：setCoreRegion / init_tiers / add* / allocate*Memory /
//                        getModuleFromName
//   2. 坐标设置         ：setModuleLocation_* / setModuleCenter_* / moveModule_2D /
//                        randomPlacment / addNoise / save,loadNodesLocation
//   3. 线长统计         ：calcHPWL / calcWA,LSE_Wirelength_2D / calcNetBoundPins /
//                        calcModuleHPWL(fast) / moveNodesCenterToCenter
//   4. 布局行维护       ：removeBlockedSite / setChipRegion_2D / y2RowIndex / isConnected
//   5. 输出与调试       ：showDBInfo / showRows / outputBookShelf + outputAUX/Nodes/
//                        PL/Nets/SCL
//
// 贯穿全文件的两条主线：
//
// 【主线一】坐标只能通过 setModuleCenter_* / setModuleLocation_* 修改。
//   因为 Pin 存的是「相对所属 module 的 offset」，pin 的绝对坐标 absolutePos 是
//   module 位置 + offset 的**缓存**。一旦有人绕过这两个入口直接改 module 的
//   coor / center（例如老代码里的 module->setCenter_2D），该 module 上所有 pin 的
//   absolutePos 就会停留在旧值，而 Net::calcNetHPWL / WA / LSE 都是直接读这个缓存的，
//   于是线长、梯度、收敛曲线都会静默地算错 —— 这类 bug 极难复现和定位。
//   所以这两个函数的固定套路是：
//       ①（仅对可动单元）把目标坐标裁进 coreRegion
//       ② 调 module->setLocation_2D / setCenter_2D 写 coor 与 center（二者互为反推）
//       ③ 遍历 module->modulePins 逐 pin 调 calculateAbsolutePos()
//
// 【主线二】dbSiteRows[].intervals 是「静态白区」，只在 removeBlockedSite() 里重算。
//   它表示每条布局行上不被 macro / terminal 占据、可以放 std cell 的 x 区间，
//   与 std cell 当前摆在哪儿无关。合法化的 subrow、详细布局的候选空位都由它派生。
//
// 关于 3D：文件里出现的 z / layerCount / dbTiers / Tier 相关代码都属于 3D IC 的预留
//   分支，当前 2D 主流程下 z 恒为 0、dbTiers 恒为空，这些分支不会生效。
// ============================================================================

// ----------------------------------------------------------------------------
// 由所有布局行推导 coreRegion（可放置区域）与 totalRowArea（布局行总面积）
//
// 做法：coreRegion 就是所有 site row 的包围盒 —— 下边界取首行的 bottom、上边界取末行的
// 顶边，左右边界取所有行 start.x / end.x 的最小 / 最大值。同时累加每行面积得到
// totalRowArea（ePlace 用它算 whitespace 面积，进而决定 filler 总面积）。
//
// 关键前提（!!）：dbSiteRows 必须已按 bottom **升序**排列，否则 front()/back() 取到的
// 并不是最下 / 最上的行，coreRegion 会算错。Parser 读 .scl 时是按文件顺序 push 的，
// 因此这里隐含依赖输入文件里行是从下往上排列的。
// ----------------------------------------------------------------------------
void PlaceDB::setCoreRegion()
{
    float bottom = dbSiteRows.front().bottom;
    float top = dbSiteRows.back().bottom + dbSiteRows.back().height;
    float left = dbSiteRows.front().start.x;
    float right = dbSiteRows.front().end.x;

    totalRowArea = 0;
    float curRowArea = 0;
    for (SiteRow curRow : dbSiteRows)
    {
        left = min(left, curRow.start.x);
        right = max(right, curRow.end.x);
        // printf( "right= %g\n", m_coreRgn.right );
        curRowArea = (curRow.end.x - curRow.start.x) * curRow.height;
        totalRowArea += curRowArea;
    }

    //! 注意：coreRegion 只由行的包围盒决定，完全不考虑 macro / terminal 占据的面积，
    //! 所以「core 面积」里是含被 macro 挡住的部分的（showDBInfo 里会另算 fixedAreaInCore）。
    coreRegion.ll = POS_2D(left, bottom);
    coreRegion.ur = POS_2D(right, top);

    cout << "Set core region from site info: ";
    coreRegion.Print();
}

// ----------------------------------------------------------------------------
// 3D IC 的分层初始化
// ----------------------------------------------------------------------------
void PlaceDB::init_tiers()
{
    //! 疑似问题：空实现，且全工程没有调用点 —— dbTiers 永远为空。
    //!   3D IC 相关（layerCount / Tier / z 坐标）目前只是预留的骨架。
}

// ----------------------------------------------------------------------------
// 建一个可移动单元（std cell 或可移动 macro），放进 dbNodes[index]
//
// 「是不是 macro」在这里判定：**高度 > commonRowHeight 即为 macro**（行高的整数倍）。
// 反过来，高度小于行高视为非法输入，直接报错退出 —— 因为 std cell 必须严丝合缝地
// 摆进布局行，矮一截就无法对齐 site。
// 调用前必须先 allocateNodeMemory() 预分配，且 commonRowHeight 已由 Parser 从 .scl 填好。
// ----------------------------------------------------------------------------
Module *PlaceDB::addNode(int index, string name, float width, float height)
{
    Module *node = new Module(index, name, width, height);
    assert(index < dbNodes.size());
    assert(commonRowHeight > EPS); //! potential float precision problem!!
    node->isMacro = false;
    if (float_greater(height, commonRowHeight))
    {
        node->isMacro = true;
        dbMacroCount++;
    }
    else if (float_greater(commonRowHeight, height))
    {
        printf("Cell %s height ERROR: %f\n", name.c_str(), height);
        exit(-1);
    }
    dbNodes[index] = node; // memory was previously allocated
    return node;
}

// ----------------------------------------------------------------------------
// 建一个 terminal（固定的 IO pad / 预放置单元），放进 dbTerminals[index]
//   isFixed：固定不动，不参与线长优化也不被裁进 coreRegion
//   isNI   ："NI" terminal，只占位、不参与密度与线长统计（输出时写作 terminal_NI）
// 与 addNode 不同，这里**不**检查高度是否等于行高 —— terminal 允许是任意形状。
// ----------------------------------------------------------------------------
Module *PlaceDB::addTerminal(int index, string name, float width, float height, bool isFixed, bool isNI)
{
    // assert(width != 0 && height != 0);
    Module *terminal = new Module(index, name, width, height, isFixed, isNI);
    assert(index < dbTerminals.size());
    dbTerminals[index] = terminal;
    return terminal;
}

//! 把 net 放进 dbNets[net->idx]；net 对象由 Parser 事先 new 好并填好 idx
void PlaceDB::addNet(Net *net)
{
    assert(net->idx < dbNets.size());
    dbNets[net->idx] = net;
}

// ----------------------------------------------------------------------------
// 建一个 pin，挂到所属 module 与所属 net 上
//
// (xOffset, yOffset) 是相对 module **中心**的偏移（不是相对左下角！见
// Pin::calculateAbsolutePos：absolutePos = module->getCenter() + offset）。
// Pin 构造函数里会做双向注册（module->modulePins 与 net->netPins 各 push 一次），
// 所以这里只需要在 dbPins 尾部追加并回写一个自增 id。
// dbPins 用的是 reserve，因此这里是 push_back 而不是按下标赋值（见 allocatePinMemory）。
// ----------------------------------------------------------------------------
int PlaceDB::addPin(Module *masterModule, Net *masterNet, float xOffset, float yOffset)
{
    dbPins.push_back(new Pin(masterModule, masterNet, xOffset, yOffset));
    int pinId = (int)dbPins.size() - 1;
    dbPins[pinId]->setId(pinId);
    return pinId;
}

//! 下面四个 allocate*Memory 供 Parser 在读入条目数后一次性预留容量，
//! 之后 add* 才能用 idx 直接定位（否则下标越界）。
void PlaceDB::allocateNodeMemory(int n)
{
    dbNodes.resize(n);
}

void PlaceDB::allocateTerminalMemory(int n)
{
    dbTerminals.resize(n);
}

void PlaceDB::allocateNetMemory(int n)
{
    dbNets.resize(n);
}

void PlaceDB::allocatePinMemory(int n) // difference between resize and reserve: with resize we can use index to index an element in a vector after allocation and before the element was instantiated
{
    dbPins.reserve(n); //! reserve instead of resize
    //! 用 reserve 而不是 resize 的原因：pin 是 addPin() 里 push_back 进去的，
    //! 不需要「先占座再填」。若用 resize 会先塞 n 个空指针，addPin 再 push_back
    //! 就会把 dbPins 撑成 2n 且前 n 个是野指针。
}

//! 按名字反查 Module*（nodes + terminals 都在这张表里）；查不到返回 NULL，调用方需判空
Module *PlaceDB::getModuleFromName(string name)
{
    map<string, Module *>::const_iterator ite = moduleMap.find(name);
    if (ite == moduleMap.end())
    {
        return NULL;
    }
    return ite->second;
}

// ----------------------------------------------------------------------------
// 按「左下角」定位一个 module —— 写坐标的核心入口之一（约定见文件头【主线一】）
//
// 三步：
//   ① 只对可动单元（!isFixed）做裁剪，把 [x, x+w] × [y, y+h] 夹进 coreRegion；
//      terminal 不裁剪，因为它们本来就允许压在 core 外面。
//      裁剪时故意多让出/缩进一个 EPS：保证裁完后的边界严格落在 coreRegion 内部，
//      避免 ePlace 用「(x - coreRegion.ll.x) / binStep」算 bin 下标时因为浮点误差
//      取到 -1 或越界的 bin（binNodeDensityUpdate 里就有 assert(binStartIdx >= 0)）。
//   ② module->setLocation_2D() 写 coor，并由 coor 反推出 center。
//   ③ 遍历该 module 上所有 pin 重算 absolutePos —— 这步是线长/梯度正确的前提。
// ----------------------------------------------------------------------------
void PlaceDB::setModuleLocation_2D(Module *module, float x, float y)
{
    //? use high precision comparison functions in global.h??
    //  check if x,y are legal(inside the chip)
    if (!module->isFixed)
    {
        if (x < coreRegion.ll.x)
        {
            x = coreRegion.ll.x;
        }
        if (x + module->width > coreRegion.ur.x)
        {
            x = coreRegion.ur.x - module->width - EPS;
        }
        if (y < coreRegion.ll.y)
        {
            y = coreRegion.ll.y;
        }
        if (y + module->height > coreRegion.ur.y)
        {
            y = coreRegion.ur.y - module->height - EPS;
        }
        //! 疑似问题：四条裁剪是「先按左下、再按右上」各判一次，彼此不联动。
        //!   若单元比 coreRegion 还宽/高，右上的裁剪会把 x 推到 ll.x 左边，
        //!   结果单元反而越出左边界（且不再有第二次检查兜住）。
    }

    module->setLocation_2D(x, y);
    for (Pin *curPin : module->modulePins)
    {
        curPin->calculateAbsolutePos();
    }
}

//! 重载：接受 POS_3D，只取 x / y（z 恒 0，2D 布局忽略）
void PlaceDB::setModuleLocation_2D(Module *module, POS_3D pos)
{
    setModuleLocation_2D(module, pos.x, pos.y);
}

// ----------------------------------------------------------------------------
// 按「中心」定位一个 module —— ePlace / Nesterov 优化器主用的入口
//
// 与 setModuleLocation_2D 的差别只在「参考点是中心还是左下角」：
// 裁剪判据换成 [x-w/2, x+w/2] × [y-h/2, y+h/2] 是否越界，最后调 setCenter_2D
// 由 center 反推 coor。刷新 pin 的那一步完全相同，也同样必不可少。
//
// 为什么优化器更爱用中心：梯度、线长模型、密度力都是围绕单元中心定义的，
// 用中心口径可以避免每次都做 ±w/2 的换算，也减少一次浮点截断误差。
// ----------------------------------------------------------------------------
void PlaceDB::setModuleCenter_2D(Module *module, float x, float y)
{
    //? use high precision comparison functions in global.h??
    //  check if x,y are legal(inside the chip)
    if (!module->isFixed)
    {
        if (x - 0.5 * module->width < coreRegion.ll.x)
        {
            // x = coreRegion.ll.x + 0.5 * module->width;
            x = coreRegion.ll.x + 0.5 * module->width + EPS;
        }
        if (x + 0.5 * module->width > coreRegion.ur.x)
        {
            // x = coreRegion.ur.x - 0.5 * module->width;
            x = coreRegion.ur.x - 0.5 * module->width - EPS;
        }
        if (y - 0.5 * module->height < coreRegion.ll.y)
        {
            // y = coreRegion.ll.y + 0.5 * module->height;
            y = coreRegion.ll.y + 0.5 * module->height + EPS;
        }
        if (y + 0.5 * module->height > coreRegion.ur.y)
        {
            // y = coreRegion.ur.y - 0.5 * module->height;
            y = coreRegion.ur.y - 0.5 * module->height - EPS;
        }
        //! 疑似问题：与 setModuleLocation_2D 同样的毛病 —— 两个方向的裁剪互不联动，
        //!   超大单元（宽/高大于 coreRegion）会被推到对侧边界之外。
    }

    module->setCenter_2D(x, y);
    for (Pin *curPin : module->modulePins)
    {
        curPin->calculateAbsolutePos();
    }
}

//! 下面两个重载分别接受 POS_3D 与 VECTOR_3D（后者来自优化器算出的位移/坐标），
//! 都只取 x / y 转发给上面的 float 版本
void PlaceDB::setModuleCenter_2D(Module *module, POS_3D pos)
{
    setModuleCenter_2D(module, pos.x, pos.y);
}

void PlaceDB::setModuleCenter_2D(Module *module, VECTOR_3D pos)
{
    setModuleCenter_2D(module, pos.x, pos.y);
}

// ----------------------------------------------------------------------------
// 只读版裁剪：返回「若把中心放到 (x,y)，裁剪后的合法中心」，不写回 module
// 供调用方先试探再决定是否落盘（例如合法化里先算目标位、再确认不越界）。
// ----------------------------------------------------------------------------
POS_3D PlaceDB::getValidModuleCenter_2D(Module *module, float x, float y)
{
    //! 疑似问题：这里**没有**像 setModuleCenter_2D 那样先判 !module->isFixed，
    //!   所以对 terminal 也会做 coreRegion 裁剪，与写坐标版本的行为不一致。
    if (x - 0.5 * module->width < coreRegion.ll.x)
    {
        // x = coreRegion.ll.x + 0.5 * module->width;
        x = coreRegion.ll.x + 0.5 * module->width + EPS;
    }
    if (x + 0.5 * module->width > coreRegion.ur.x)
    {
        // x = coreRegion.ur.x - 0.5 * module->width;
        x = coreRegion.ur.x - 0.5 * module->width - EPS;
    }
    if (y - 0.5 * module->height < coreRegion.ll.y)
    {
        // y = coreRegion.ll.y + 0.5 * module->height;
        y = coreRegion.ll.y + 0.5 * module->height + EPS;
    }
    if (y + 0.5 * module->height > coreRegion.ur.y)
    {
        // y = coreRegion.ur.y - 0.5 * module->height;
        y = coreRegion.ur.y - 0.5 * module->height - EPS;
    }
    POS_3D validPosition;
    validPosition.x = x;
    validPosition.y = y;
    //! 疑似问题：validPosition 的 z 分量沿用 POS_3D 默认构造的 0，没有从入参或
    //!   module 继承。2D 布局下无影响，3D 场景会丢 z。
    return validPosition;
}

//! 设置朝向（N / S / W / E 及其翻转，见 global.h 的 ORIENT）。
//! 注意：改朝向**不会**交换 width / height，也不会刷新 pin —— 因为 pin 的 offset
//! 是相对中心的，而 getUR_2D() 一律按未旋转计算（见 objects.cpp）。
void PlaceDB::setModuleOrientation(Module *module, int orientation)
{
    module->setOrientation(orientation);
}

// ----------------------------------------------------------------------------
// 把单元随机撒在 coreRegion 内的某个合法位置（按左下角定位）
//
// 做法：rand() 归一化到 [0,1) 后乘「可摆动范围」= coreRegion 尺寸 − 单元尺寸，
// 这样保证左下角取到任何随机数时单元都不会越出右上边界。
// 最后仍走 setModuleLocation_2D，由它再做一次裁剪并刷新 pin。
// ----------------------------------------------------------------------------
void PlaceDB::setModuleLocation_2D_random(Module *module)
{
    assert(module);
    float x = rand();
    float y = rand();

    assert(coreRegion.ur.x > coreRegion.ll.x);
    assert(coreRegion.ur.y > coreRegion.ll.y);

    float potentialRegionWidth = coreRegion.ur.x - coreRegion.ll.x; // potential region for randomly place module
    float potentialRegionHeight = coreRegion.ur.y - coreRegion.ll.y;

    potentialRegionHeight -= module->getHeight();
    potentialRegionWidth -= module->getWidth();
    //! 疑似问题：这里没有把生成的坐标再平移回 coreRegion.ll —— rand() 归一化后乘的是
    //!   「范围宽度」而不是「宽度后再加上 ll」，所以随机结果落在 [0, 范围] 而非
    //!   [ll, ur-w]，实际上完全依赖下面 setModuleLocation_2D 的裁剪把它拉回左下边界。
    //!   另外若单元比 coreRegion 大，potentialRegion* 会变成负数，随机坐标直接落到
    //!   coreRegion 左下方，同样只能靠裁剪兜底。
    //!   还有：用的是 rand() 且从未 srand()，每次运行的随机序列相同。

    float RAND_MAX_INVERSE = (float)1.0 / RAND_MAX;
    x = x * RAND_MAX_INVERSE * potentialRegionWidth;
    y = y * RAND_MAX_INVERSE * potentialRegionHeight;

    setModuleLocation_2D(module, x, y);
}

//! 按位移平移一个可动单元（浮点版）：取当前中心 → 加位移 → 走 setModuleCenter_2D
void PlaceDB::moveModule_2D(Module *module, VECTOR_2D delta)
{
    POS_3D curPos = module->getCenter();
    setModuleCenter_2D(module, curPos.x + delta.x, curPos.y + delta.y);
}

//! 按位移平移一个可动单元（整数版），给按 bin / 行数做整数位移的场合用
void PlaceDB::moveModule_2D(Module *module, VECTOR_2D_INT delta)
{
    VECTOR_2D_INT curPos;
    //! 疑似问题：先把 float 中心**截断**成 int，平移后又变回 float 写回，
    //!   会凭空引入最多 1 个单位的取整误差（单元被吸附到整数网格）。
    //!   若本意只是「加一个整数位移」，这里不该对 curPos 取整。
    curPos.x = module->getCenter().x;
    curPos.y = module->getCenter().y;
    setModuleCenter_2D(module, curPos.x + delta.x, curPos.y + delta.y);
}

//! 把所有可动单元随机撒开，用作「完全没有初始解」时的起点（正式流程走 QPlace）
void PlaceDB::randomPlacment()
{
    for (Module *curModule : dbNodes)
    {
        setModuleLocation_2D_random(curModule);
    }
}

// ----------------------------------------------------------------------------
// 下面四个是全设计的线长总和，遍历 dbNets 累加每个 net 的贡献即可。
//
// 三种线长模型的分工：
//   calcHPWL            精确半周线长（bounding box 半周长）。不可微，只用于**评估**
//                       与收敛判据，不进优化器。
//   calcWA/LSE_Wirelength_2D
//                       HPWL 的两种**平滑可微**近似（weighted-average / log-sum-exp），
//                       梯度下降真正优化的对象。invertedGamma = 1/γ，γ 越小越贴近 HPWL
//                       但也越"尖"，由 ePlace 在迭代中按溢出率自适应调整。
//
// 共同前提（!!）：这些计算直接读 pin 的 absolutePos 缓存，所以调用前所有 pin 的绝对
// 坐标必须是最新的 —— 也就是所有移动都必须经过 setModuleCenter*/setModuleLocation*。
// ----------------------------------------------------------------------------

double PlaceDB::calcHPWL() //! parallel this?
{
    double HPWL = 0;
    for (Net *curNet : dbNets)
    {
        HPWL += curNet->calcNetHPWL();
    }
    return HPWL;
}

double PlaceDB::calcWA_Wirelength_2D(VECTOR_2D invertedGamma)
{
    double WA = 0;
    for (Net *curNet : dbNets)
    {
        WA += curNet->calcWirelengthWA_2D(invertedGamma);
    }
    return WA;
}

double PlaceDB::calcLSE_Wirelength_2D(VECTOR_2D invertedGamma)
{
    double LSE = 0;
    for (Net *curNet : dbNets)
    {
        LSE += curNet->calcWirelengthLSE_2D(invertedGamma);
    }
    return LSE;
}

// ----------------------------------------------------------------------------
// 与 calcHPWL 数值等价，但走 Net::calcBoundPin()：顺带把每个 net 的
// boundPinXmin/Xmax/Ymin/Ymax(/Zmin/Zmax) 刷新出来。
//
// 这个「副作用」才是它存在的主要理由 —— QPlace 的二次布局只在边界 pin 之间建边，
// 必须先靠这一步把边界 pin 求出来（QPlace 里 assert(HPWL == HPWL2) 就是在校验
// 两种算法一致、间接确认 boundPin 已经更新）。
// ----------------------------------------------------------------------------
double PlaceDB::calcNetBoundPins()
{
    double HPWL = 0;
    for (Net *curNet : dbNets)
    {
        HPWL += curNet->calcBoundPin();
    }
    return HPWL;
}

// ----------------------------------------------------------------------------
// 单个 module 的「相关线长」：把它参与的所有 net 的 HPWL 加起来
//
// 语义：这不是该单元的线长，而是**移动它会影响到**的线长总量，用来评估移动的收益
//（详细布局里衡量「把这个 cell 挪走能省多少线长」）。
// ----------------------------------------------------------------------------
double PlaceDB::calcModuleHPWL(Module *curModule) //! assume there are no 2 pins in one module belong to a same net
{
    double HPWL = 0;
    for (Net *curNet : curModule->nets)
    {
        HPWL += curNet->calcNetHPWL();
    }
    return HPWL;
}

// double PlaceDB::calcModuleHPWLunsafe(Module *curModule)
// {
//     //!!!!this function is dangerous and is used for accelerating macro legalization only
//     double HPWL = 0;
//     for (Pin *curModulePin : curModule->modulePins)
//     {
//         // HPWL += curModulePin->net->calcNetHPWL();
//         float maxX = -FLOAT_MAX;
//         float minX = FLOAT_MAX;
//         // double maxY = DOUBLE_MIN;
//         float maxY = -FLOAT_MAX;
//         float minY = FLOAT_MAX;
//         // double maxZ = DOUBLE_MIN;
//         float maxZ = -FLOAT_MAX; // potential bug: double_min >0 so boundPinZmax might be null when all z == 0
//         float minZ = FLOAT_MAX;

//         POS_3D curPos;

//         for (Pin *curPin : curModulePin->net->netPins)
//         {
//             //!!! this is why this function is unsafe!!!
//             curPos = curPin->getAbsolutePos(); //!!!!!!!! must guarantee that the absoulte pos is up to date!!!!!! this is faster than use fetchAbsolutePos, probably because less function calling overhead?
//             // curPos = curPin->fetchAbsolutePos();
//             minX = min(minX, curPos.x);
//             maxX = max(maxX, curPos.x);
//             minY = min(minY, curPos.y);
//             maxY = max(maxY, curPos.y);
//             minZ = min(minZ, curPos.z);
//             maxZ = max(maxZ, curPos.z);
//         }
//         // if (!gArg.CheckExist("3DIC"))
//         // {
//         //     //? assert(maxZ == minZ == 0); this causes bug
//         //     assert(float_equal(maxZ, 0.0));
//         //     assert(float_equal(minZ, 0.0));
//         // }
//         HPWL += ((maxX - minX) + (maxY - minY) + (maxZ - minZ));
//     }
//     return HPWL;
// }

// ----------------------------------------------------------------------------
// calcModuleHPWL 的加速版（实测更快）：不再逐 net 调 calcNetHPWL()，而是
// 直接遍历本 module 的每个 pin、读缓存的 absolutePos，现场求该 net 的包围盒。
//
// 快在哪：省掉函数调用开销，且省掉 calcNetHPWL 里对 netPins 的重复遍历准备。
// 代价（!! 这也是它"不安全"的地方）：**直接读 curPin->absolutePos 缓存**，
// 因此调用前必须保证所有 pin 的绝对坐标已刷新，否则算的是过期值。
// ----------------------------------------------------------------------------
double PlaceDB::calcModuleHPWLfast(Module *curModule) // tested faster than calcModuleHPWL
{
    double HPWL = 0;
    for (Pin *curModulePin : curModule->modulePins)
    {
        // HPWL += curModulePin->net->calcNetHPWL();
        float maxX = -FLOAT_MAX;
        float minX = FLOAT_MAX;
        // double maxY = DOUBLE_MIN;
        float maxY = -FLOAT_MAX;
        float minY = FLOAT_MAX;
        // double maxZ = DOUBLE_MIN;
        //! 注意这里刻意用 -FLOAT_MAX 而不是 DOUBLE_MIN —— 因为 DOUBLE_MIN 是最小的
        //! **正**数（见 global.h），拿它当"负无穷大"初始化会让 maxZ 恒为正，
        //! 一旦所有 z 都为 0 就永远取不到正确的上界。
        float maxZ = -FLOAT_MAX; // potential bug: double_min >0 so boundPinZmax might be null when all z == 0
        float minZ = FLOAT_MAX;

        POS_3D curPos;

        for (Pin *curPin : curModulePin->net->netPins)
        {
            // curPos = curPin->getAbsolutePos();
            curPos = curPin->absolutePos;
            // curPos = curPin->fetchAbsolutePos();
            minX = min(minX, curPos.x);
            maxX = max(maxX, curPos.x);
            minY = min(minY, curPos.y);
            maxY = max(maxY, curPos.y);
            minZ = min(minZ, curPos.z);
            maxZ = max(maxZ, curPos.z);
        }
        //! 疑似问题：外层遍历的是本 module 的**每个 pin**，若同一 module 上有 2 个 pin
        //!   属于同一个 net，这条 net 的 HPWL 会被重复累加两次（原注释也点出了这个前提
        //!   "assume there are no 2 pins in one module belong to a same net"）。
        // if (!gArg.CheckExist("3DIC"))
        // {
        //     //? assert(maxZ == minZ == 0); this causes bug
        //     assert(float_equal(maxZ, 0.0));
        //     assert(float_equal(minZ, 0.0));
        // }
        //! 2D 布局下 z 恒为 0，(maxZ - minZ) 恒等于 0，这一项只有 3D IC 才有意义。
        HPWL += ((maxX - minX) + (maxY - minY) + (maxZ - minZ));
    }
    return HPWL;
}

// ----------------------------------------------------------------------------
// 把所有可动单元搬到 coreRegion 正中心 —— 二次布局（QPlace）的迭代起点
//
// 为什么先堆到中心：QPlace 求解的是纯二次线长模型，没有任何扩散/密度力，
// 若起点分散，解会被各自的初始位置"锚住"；全部堆到同一点后，单元只能靠所连 net
// 的拉力和 terminal 的支撑散开，解完全由连接关系决定，与初始随机性无关。
// ----------------------------------------------------------------------------
void PlaceDB::moveNodesCenterToCenter()
{
    POS_2D coreRegionCenter(0.5 * (coreRegion.ur.x + coreRegion.ll.x), 0.5 * (coreRegion.ur.y + coreRegion.ll.y));
    for (Module *curModule : dbNodes)
    {
        setModuleCenter_2D(curModule, coreRegionCenter.x, coreRegionCenter.y);
    }
}

// ============================================================================
// removeBlockedSite() —— 重算每条 site row 的 intervals（可放置 x 区间）
//
// 【它在做什么】
//   Parser 读入 .scl 时，每行只有 1 个覆盖整行的 interval。但真实版图里，行会被
//   macro（高于行高的硬宏）和压在 core 内的 terminal 横插一脚，那些位置根本放不了
//   std cell。本函数就是把这些「障碍物」从行上挖掉，得到若干段互不相交的可用区间：
//
//       行原本:  [==================================================]
//       挖掉后:  [======]        [============]        [============]
//                        ^macro^              ^terminal^
//
//   挖出来的 intervals 是**静态白区**（只取决于 macro/terminal 的位置，与 std cell
//   当前摆在哪儿无关），合法化的 subrow、详细布局的候选空位全部由它派生。
//
// 【三步流程】
//   1. 收集障碍物：与 coreRegion 有重叠的 terminal + 所有 macro，按 x（其次 y）排序。
//      排序是为了让挖除过程按从左到右推进，也让"重叠"能被立刻发现。
//   2. 逐障碍物挖除：先按 y 定位它压住了哪几行，再对每行的每个 interval 做 4 种
//      情形的区间减法（完全覆盖 / 切右端 / 切左端 / 中间挖洞一分为二）。
//   3. site 对齐：把区间端点吸附到 site 网格上（以 coreRegion.ll.x 为原点），
//      并丢弃宽度不足一个 site 的碎片 —— 对应 ntuplace 的 FixFreeSiteBySiteStep()。
//
// 【调用契约（!!）】
//   · 必须在 **macro 合法化之后** 调用：若 macro 之间还有重叠，本函数直接
//     "TERMINALS/MACROS OVERLAP!" 退出；即便不退出，切出来的区间也会偏大，
//     导致 std cell 被合法化到 macro 底下的非法区域。
//   · 只处理 dbSiteRows 的 y 维度；2D 流程下无 z 概念。
// ============================================================================
void PlaceDB::removeBlockedSite() // update intervals
{
    //! ignore terminals that are outside the core region
    //! overlap between macros should be eliminated first! (macro legalization)

    // 1. count all terminal and macros
    vector<CRect> obstacles;
    obstacles.clear();
    //! 只收与 coreRegion 有**实际重叠**的 terminal：完全在 core 外面的 pad 不影响
    //! 行内可用空间，收进来反而会污染下面那个"重叠即报错"的严格判据。
    for (Module *curTerminal : dbTerminals)
    {
        if (float_greater(coreRegion.ur.y, curTerminal->getLL_2D().y) && float_less(coreRegion.ll.y, curTerminal->getUR_2D().y)) // overlap in y direction
        {
            if (float_greater(coreRegion.ur.x, curTerminal->getLL_2D().x) && float_less(coreRegion.ll.x, curTerminal->getUR_2D().x)) // overlap in x direction
            {
                CRect newObstacle;
                newObstacle.ll = curTerminal->getLL_2D();
                newObstacle.ur = curTerminal->getUR_2D();
                obstacles.push_back(newObstacle);
            }
        }
    }

    //! macro 无条件收进来（不论是否与 core 重叠）：isMacro 的高度大于行高，
    //! 它压住的所有行都要挖掉。
    for (Module *curNode : dbNodes)
    {
        if (curNode->isMacro)
        {
            CRect newObstacle;
            newObstacle.ll = curNode->getLL_2D();
            newObstacle.ur = curNode->getUR_2D();
            obstacles.push_back(newObstacle);
        }
    }
    //! sort obstacles by x coordinate first
    sort(obstacles.begin(), obstacles.end(), [=](CRect a, CRect b)
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
                //! 左下角完全重合 = 必然互相重叠，无法正确做区间减法，直接终止。
                //! 注意这里只判了「左下角相同」，两个**部分**重叠（左下角不同）的
                //! 障碍物并不会被这条分支拦住 —— 那是调用方（macro 合法化）的责任。
                // terminals/macros overlap!!
                cerr<<"TERMINALS/MACROS OVERLAP!\n";
                exit(0);
            } });
            //! 疑似问题：lambda 用 [=] 捕获了 this（不需要），且参数按值传 CRect
            //!   （每行一个障碍物就要拷贝一次）；改成 [](const CRect &a, const CRect &b)
            //!   效率更好。另外这里的 exit(0) 属于"异常即杀进程"，上层无法降级处理。

    // 2. remove blocked sites and update intervals
    double siteStep = dbSiteRows.front().step; //! assume step for all rows are identical! so it's ok to use front()

    for (CRect curObstacle : obstacles)
    {
        vector<SiteRow>::iterator iteBeginRow, iteEndRow;
        // find the begin row and the end row (of all rows that are blocked by curObstacle)
        //!!!! important assumption here: dbSiteRows is sorted by bottom coordinate in an increasing order!
        //! 由于 dbSiteRows 按 bottom 升序，用两个线性扫描就能 O(行数) 定位
        //! 「被这个障碍物压住的第一行 / 最后一行」，无需二分、更无需全表比对。
        for (iteBeginRow = dbSiteRows.begin(); iteBeginRow < dbSiteRows.end(); iteBeginRow++)
        {
            if (iteBeginRow->bottom + iteBeginRow->height > curObstacle.ll.y)
            {
                break;
            }
        }

        for (iteEndRow = iteBeginRow; iteEndRow < dbSiteRows.end(); iteEndRow++)
        {
            if (iteEndRow->bottom + iteEndRow->height >= curObstacle.ur.y)
            {
                break;
            }
        }

        //! 障碍物顶边高于所有行的顶边（例如 terminal 压在最上面一行之上）时，
        //! 上面的循环会扫到 end()，这里回退一格把末行纳入处理范围。
        if (iteEndRow == dbSiteRows.end())
        {
            iteEndRow--;
        }
        assert(iteBeginRow != dbSiteRows.end());
        //! 疑似问题：dbSiteRows 为空时 dbSiteRows.end()-- 是未定义行为；
        //!   且 assert 在 NDEBUG（Release）下会被编译掉，届时会直接用野迭代器。

        Interval tempInterval;

        //! 对 [iteBeginRow, iteEndRow] 闭区间内的每一行做区间减法。
        //! 注意是**闭**区间，所以 iteEndRow 那一行也会被处理。
        for (vector<SiteRow>::iterator curRow = iteBeginRow; curRow <= iteEndRow; curRow++)
        {
            for (int i = 0; i < (signed)curRow->intervals.size(); i++)
            {
                tempInterval = curRow->intervals[i];

                if (tempInterval.start >= curObstacle.ur.x || tempInterval.end <= curObstacle.ll.x) // screen unnecessary checks
                {
                    continue;
                }

                //! 四种情形就是「区间 [start,end] 减去障碍物 [ll.x,ur.x]」的标准做法：
                if (tempInterval.start >= curObstacle.ll.x && tempInterval.end <= curObstacle.ur.x) // fully blocked
                {
                    //    ---
                    // MMMMMMMMM
                    //! ① 区间被完全吃掉 → 整个删掉。
                    //! 疑似问题：erase 之后没有 i--，原本位于 i+1 的区间会前移到 i，
                    //!   而循环末尾的 i++ 会**跳过**它，导致它不再与本障碍物比对。
                    //!   同一行内两个区间同时与本障碍物相交时才会触发（需要行内已存在缝隙），
                    //!   属于罕见的潜在漏处理。
                    curRow->intervals.erase(curRow->intervals.begin() + i);
                }
                else if (tempInterval.end > curObstacle.ur.x && tempInterval.start >= curObstacle.ll.x)
                {
                    // ------      -----
                    // MMM      MMMMM
                    //! ② 左端被吃 → 把 start 右移到障碍物右边界。
                    curRow->intervals[i].start = curObstacle.ur.x;
                }
                else if (tempInterval.start < curObstacle.ll.x && tempInterval.end <= curObstacle.ur.x)
                {
                    // ---------       -----
                    //     MMMMM          MMMMM
                    //! ③ 右端被吃 → 把 end 左移到障碍物左边界。
                    curRow->intervals[i].end = curObstacle.ll.x;
                }
                else if (tempInterval.start < curObstacle.ll.x && tempInterval.end > curObstacle.ur.x)
                {
                    // -----------
                    //    MMMM
                    //! ④ 障碍物在中间 → 一分为二：左半截收尾到障碍物左边界，
                    //!   再在 i+1 处插入右半截（从障碍物右边界到原 end）。
                    //!   插入后 intervals.size() 变大，外层 i 循环会一并扫到新增的那段
                    //!   （它由 start == ur.x 命中第 ① 个 continue，不会重复挖）。
                    curRow->intervals[i].end = curObstacle.ll.x;
                    curRow->intervals.insert(curRow->intervals.begin() + i + 1, Interval(curObstacle.ur.x, tempInterval.end));
                }
                else
                {
                    printf("Warning: Module Romoving Error\n");
                    // exit(-1);
                    //! 疑似问题：走到这里说明上面四种情形都没覆盖到（理论上不可能）。
                    //!   只打 warning 不退出，区间会保持"被挖之前"的错误状态；
                    //!   另外 "Romoving" 疑为 "Removing" 的拼写错误。
                }
            }
        }
    }

    //! 3. align intervals to sites after updating intervals, see ntuplace: FixFreeSiteBySiteStep(). Here we need to update the end and start of a site row, so end.x-start.x is an positive integer multiple of site step(site width)
    //! 第 3 步：把区间端点吸附到 site 网格，并丢掉放不下一个 site 的碎片。
    //!
    //! 为什么要对齐：合法化时 std cell 必须整 site 对齐地摆进 subrow。若区间端点
    //! 落在两个 site 中间，贴着端点放的 cell 就会错位半个 site，既压不住整条 site
    //! 又可能越出可用区。所以这里把左端**向上**取整（ceil）、右端**向下**取整（floor），
    //! 都吸附到「coreRegion.ll.x + k × siteStep」这个网格上，区间宽度自然成为
    //! siteStep 的正整数倍。
    //! 注意对齐基准是 coreRegion.ll.x 而不是绝对 0 —— 两者不同（见 legalizer.cpp 里
    //! 另一处以 0 为基准的对齐，两处口径不一致）。
    //!
    //! 为什么还要丢碎片：吸附后宽度不足一个 siteStep 的区间（例如被 macro 蹭掉一小
    //! 截的角落）实际上放不下任何 cell，留着只会让合法化徒劳地枚举空位。
    for (auto curRowIter = dbSiteRows.begin(); curRowIter != dbSiteRows.end(); curRowIter++)
    {
        //? should start.x and end.x be integers too??? check ntuplace

        for (auto curIntervalIter = curRowIter->intervals.begin(); curIntervalIter != curRowIter->intervals.end();)
        {
            double intervalWidth = curIntervalIter->end - curIntervalIter->start;
            // subRowWidth might be less than 0 when:
            //    ------
            //   OOOOO
            //! 宽度不足一个 site → 直接丢弃。注意 erase 返回下一个有效迭代器，
            //! 所以这个循环**不能**写 ++curIntervalIter（否则跳过一个且可能失效）。
            if (float_less(intervalWidth, siteStep)) // assume iter->step > 0!
            {
                curIntervalIter = curRowIter->intervals.erase(curIntervalIter);
            }
            else
            {
                //! 左端向右吸附（ceil），右端向左吸附（floor），都以 coreRegion.ll.x 为原点
                double newLeft = ceil((curIntervalIter->start - coreRegion.ll.x) / siteStep) * siteStep + coreRegion.ll.x;
                double newRight = floor((curIntervalIter->end - coreRegion.ll.x) / siteStep) * siteStep + coreRegion.ll.x;
                double newIntervalWidth = newRight - newLeft;
                // assert(newIntervalWidth >= siteStep);
                if (float_greater(newIntervalWidth, 0.0))
                {
                    curIntervalIter->start = newLeft;
                    curIntervalIter->end = newRight;
                    curIntervalIter++;
                }
                else
                {
                    //! 吸附后宽度为 0（或负）→ 这段已不可用，丢弃。
                    curIntervalIter = curRowIter->intervals.erase(curIntervalIter);
                    if (float_less(newIntervalWidth, 0.0))
                    {
                        // newIntervalWidth should >= 0.0
                        cerr << "sub row new width < 0 when it should not\n";
                        exit(0);
                        //! 疑似问题：newIntervalWidth < 0 说明左端吸附后跑到了右端右边，
                        //!   属于不该发生的状态；这里在 erase 之后才报错退出，
                        //!   而且 exit(0)（成功码）会让外部脚本误判为正常结束。
                    }
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------
// 计算 chipRegion —— 仅供画图/可视化使用的「芯片外框」
//
// 做法：从 coreRegion 出发，把所有 terminal（含压在 core 外的 pad）都包进去。
// 它和 coreRegion 的区别很重要：coreRegion 是可放置区域（布局、bin、密度都用它），
// chipRegion 只是「画出来的那张图有多大」。
// ----------------------------------------------------------------------------
void PlaceDB::setChipRegion_2D()
{
    chipRegion = coreRegion;
    for (Module *curTerminal : dbTerminals)
    {
        assert(curTerminal);
        chipRegion.ll.x = min(chipRegion.ll.x, curTerminal->getLL_2D().x);
        chipRegion.ll.y = min(chipRegion.ll.y, curTerminal->getLL_2D().y);
        // gmin.z = min(gmin.z, curTerminal->pmin.z);

        chipRegion.ur.x = max(chipRegion.ur.x, curTerminal->getUR_2D().x);
        chipRegion.ur.y = max(chipRegion.ur.y, curTerminal->getUR_2D().y);
        // gmax.z = max(gmax.z, curTerminal->pmax.z);
        //! 上面两行 z 相关代码被注释掉了：3D IC 下 chipRegion 是 CRect（2D），
        //! 装不下 z 方向的外扩，因此 z 分支目前不生效。
    }
    //!!!!!!!chipRegion.ll!=(0,0)
    //! 提醒：chipRegion 的原点通常不是 (0,0)（terminal 可能带负坐标），
    //! 画图时做坐标换算必须先减去 ll，不能直接用绝对坐标当像素坐标。
}

// ----------------------------------------------------------------------------
// 打印数据库摘要（<<<< DATABASE SUMMARIES >>>>），流程开始时调用一次
//
// 重点看三个指标：
//   Placement Util. = 可动面积 / (core 面积 − core 内固定面积)   —— 越大越难布
//   Core Density    = (可动 + core 内固定) / core 面积
//   Max net degree  —— 决定线长计算与矩阵规模的上界
//
// 注意 cellArea / macroArea 的划分就按 Module::isMacro（即高度是否大于行高）。
// ----------------------------------------------------------------------------
void PlaceDB::showDBInfo()
{
    double coreArea = coreRegion.getArea();
    double cellArea = 0;
    double macroArea = 0;
    double fixedArea = 0;
    double fixedAreaInCore = 0;
    double movableArea = 0;
    int macroCount = 0;
    int cellCount = 0;
    int terminalCount = dbTerminals.size();
    //! 疑似问题：下面两个局部变量与 PlaceDB 的同名**成员** netCount / pinCount 重名，
    //!   在 showDBInfo 内会遮蔽成员。成员因此既没被赋值也没被读到（见 placedb.h 的标注）。
    int netCount = dbNets.size();
    int pinCount = dbPins.size();
    int maxNetDegree = INT32_MIN;
    for (Module *curNode : dbNodes)
    {

        if (curNode->isMacro)
        {
            macroArea += curNode->getArea();
            macroCount++;
        }
        else
        {
            cellArea += curNode->getArea();
            cellCount++;
        }
    }
    movableArea = macroArea + cellArea;
    for (Module *curTerminal : dbTerminals)
    {
        fixedArea += curTerminal->getArea();
        fixedAreaInCore += getOverlapArea_2D(coreRegion.ll, coreRegion.ur, curTerminal->getLL_2D(), curTerminal->getUR_2D()); // notice that here we do not consider that some rows might have shorter width than the others
        //! 注：fixedAreaInCore 只是「terminal 与 coreRegion 矩形的重叠面积」，
        //! 没扣掉那些被 macro 挡住、本来也放不了 cell 的 site，所以是个偏乐观的估计。
    }
    int pin2 = 0, pin3 = 0, pin10 = 0, pin100 = 0;
    for (Net *curNet : dbNets)
    {
        int curPinCount = curNet->getPinCount();
        if (curPinCount > maxNetDegree)
        {
            maxNetDegree = curPinCount;
        }
        if (curPinCount == 2)
            pin2++;
        else if (curPinCount < 10)
            pin3++;
        else if (curPinCount < 100)
            pin10++;
        else
            pin100++;
        //! 疑似问题：pin3 这一档实际统计的是「度数 <10 且 !=2」，会把 1 度 net
        //!   （只有一个 pin，HPWL 恒为 0）也数进去，而打印文案写的是 "3-10"。
        //!   另外这里的分档用的是 <10 / <100，与打印的 "3-10"/"11-100" 边界略有出入。
    }

    printf("\n<<<< DATABASE SUMMARIES >>>>\n\n");
    printf("         Core region: ");
    coreRegion.Print();
    printf("   Row Height/Number: %.0f / %d (site step %f)\n", commonRowHeight, dbSiteRows.size(), dbSiteRows[0].step);
    printf("           Core Area: %.0f (%g)\n", coreArea, coreArea);
    printf("           Cell Area: %.0f (%.2f%%)\n", cellArea, 100.0 * cellArea / coreArea);
    if (macroCount > 0)
    {
        printf("          Macro Area: %.0f (%.2f%%)\n", macroArea, 100.0 * macroArea / coreArea);
        printf("  Macro/(Macro+Cell): %.2f%%\n", 100.0 * macroArea / (macroArea + cellArea));
    }
    printf("        Movable Area: %.0f (%.2f%%)\n", movableArea, 100.0 * movableArea / coreArea);
    if (terminalCount > 0)
    {
        printf("          Fixed Area: %.0f (%.2f%%)\n", fixedArea, 100.0 * fixedArea / coreArea);
        printf("  Fixed Area in Core: %.0f (%.2f%%)\n", fixedAreaInCore, 100.0 * fixedAreaInCore / coreArea);
    }
    // printf( "   (Macro+Cell)/Core: %.2f%%\n", 100.0*(macroArea+cellArea)/coreArea );
    printf("     Placement Util.: %.2f%% (=move/freeSites)\n", 100.0 * movableArea / (coreArea - fixedAreaInCore));
    printf("        Core Density: %.2f%% (=usedArea/core)\n", 100.0 * (movableArea + fixedAreaInCore) / coreArea);
    //! 疑似问题：两处的除数 (coreArea - fixedAreaInCore) 可能为 0 或负（terminal 把
    //!   core 全占满时），会得到 inf / nan 而不是报错。
    // printf( "           Site Area: %.0f (%.0f)", totalSiteArea, coreArea-fixedAreaInCore );
    printf("              Cell #: %d (=%dk)\n", cellCount, (cellCount / 1000));
    printf("            Object #: %d (=%dk) (fixed: %d) (macro: %d)\n", macroCount + cellCount + terminalCount, (macroCount + cellCount + terminalCount) / 1000, terminalCount, macroCount);
    //! 注：上面两处的 "=%dk" 用的是整数除法，例如 1500 个 cell 会显示 "=1k"（截断）。
    if (macroCount < 20)
    {
        for (Module *curNode : dbNodes)
            if (curNode->isMacro)
                printf(" Macro: %s\n", curNode->name.c_str());
    }
    printf("               Net #: %d (=%dk)\n", netCount, netCount / 1000);
    printf("               Max net degree=: %d\n", maxNetDegree);
    printf("                  Pin 2 (%d) 3-10 (%d) 11-100 (%d) 100- (%d)\n", pin2, pin3, pin10, pin100);
    printf("               Pin #: %d\n", pinCount);

    // printf( "               Pin #: %d (in: %d  out: %d  undefined: %d)\n", pinNum, inPinNum, outPinNum, undefPinNum );
    // double HPWL = calcHPWL();
    // printf("     Pin-to-Pin HPWL: %.0f (%g)\n", HPWL, HPWL);
}

// ----------------------------------------------------------------------------
// 逐行打印 intervals（调试用，配合 removeBlockedSite 检查挖除结果）
// ----------------------------------------------------------------------------
void PlaceDB::showRows()
{
    for (auto curRowIter = dbSiteRows.begin(); curRowIter != dbSiteRows.end(); curRowIter++)
    {
        cout << "\n=====DB ROW SPACE ===\n";

        for (auto iter = curRowIter->intervals.begin(); iter != curRowIter->intervals.end(); iter++)
        {
            // modified by Jin 20070727
            //! 注意打印的是 [起点, **长度**] 而不是 [起点, 终点]（第二个参数是
            //! getLength()），读日志时容易误读成区间右端点。
            printf("[%.10f,%.10f] ", iter->start, iter->getLength());
            // cout<<" ["<<iter->first<<","<<iter->second<<"] ";
            // modified by Jin 20070727
        }
        cout << '\n';
    }
}

// ----------------------------------------------------------------------------
// 输出 Bookshelf 格式的结果文件（交给 ntuplace3 之类的外部合法器，或用于评测）
//
//   plOnly = false：写全套 .aux / .nodes / .nets / .scl / .pl
//   plOnly = true ：只写 .pl（坐标变了、其余文件没变时的增量输出）
//
// 输出路径与文件名由命令行参数决定：
//   输出目录 <outputPath>（缺省 ./<benchmarkName>/），文件名 <benchmarkName>-<suffix>.xxx
// suffix 会被写回 gArg 的 "outputSuffix"，供后面几个 output* 函数各自读取。
// ----------------------------------------------------------------------------
void PlaceDB::outputBookShelf(string suffix, bool plOnly)
{
    string outputFilePath;
    string benchmarkName;
    gArg.GetString("benchmarkName", &benchmarkName);
    if (!gArg.GetString("outputPath", &outputFilePath))
    {
        outputFilePath = "./" + benchmarkName + "/";
    }

    gArg.Override("outputSuffix", suffix);

    if (!plOnly)
    {
        outputAUX();
        outputNodes();
        outputNets();
        outputSCL();
        //! 注：.aux 里声明了 .wts（线网权重文件），但这里并没有对应的 outputWts()，
        //!   也就是说 .wts 从未被真正生成过（外部工具通常允许缺该文件、按全 1 处理）。
    }

    outputPL();
}

//! 写 .aux：Bookshelf 的"目录"文件，列出本次结果由哪几个文件组成
void PlaceDB::outputAUX()
{
    string outputFilePath;
    gArg.GetString("outputPath", &outputFilePath);

    string benchmarkName;
    gArg.GetString("benchmarkName", &benchmarkName);

    string suffix;
    gArg.GetString("outputSuffix", &suffix);

    outputFilePath += benchmarkName;
    outputFilePath += "-" + suffix + ".aux";

    cout << "Output AUX file:" << outputFilePath << endl;

    //! 把最终路径写回 gArg，供后续流程（合法器）按 key 取用
    gArg.Override("outputAUX", outputFilePath);

    ofstream out(outputFilePath);
    if (!out)
    {
        cerr << "Cannot open output file\n";
        return;
    }

    out << "RowBasedPlacement : "
        << benchmarkName << "-" + suffix + ".nodes "
        << benchmarkName << "-" + suffix + ".nets "
        << benchmarkName << "-" + suffix + ".wts "
        << benchmarkName << "-" + suffix + ".pl "
        << benchmarkName << "-" + suffix + ".scl \n\n";
}

//! 写 .nodes：所有 module 的**名字与尺寸**（不含坐标，坐标在 .pl 里）
void PlaceDB::outputNodes()
{

    string outputFilePath;
    gArg.GetString("outputPath", &outputFilePath);

    string benchmarkName;
    gArg.GetString("benchmarkName", &benchmarkName);

    string suffix;
    gArg.GetString("outputSuffix", &suffix);

    //! 下面 4 个 output* 函数开头都是同一套"拼路径"样板：
    //!   <outputPath><benchmarkName>-<suffix>.<扩展名>
    //! suffix 由 outputBookShelf() 写进 gArg["outputSuffix"]，这里再读出来。
    outputFilePath += benchmarkName;
    outputFilePath += "-" + suffix + ".nodes";

    cout << "Output Nodes file:" << outputFilePath << endl;

    FILE *out;
    out = fopen(outputFilePath.c_str(), "w");
    if (!out)
    {
        cerr << "Cannot open output file\n";
        return;
    }

    fprintf(out, "UCLA nodes 1.0\n\n");
    fprintf(out, "NumNodes : %d\n", moduleCount);
    fprintf(out, "NumTerminals : %d\n\n", dbTerminals.size());
    //! 疑似问题：NumNodes 用的是成员 moduleCount，而 dbNodes 才是实际容器大小；
    //!   moduleCount 由 Parser 填写，若没填（构造函数里是 -1）这里就会写出 -1，
    //!   与下面真正列出的条目数不一致。改用 dbNodes.size() 更稳妥。

    // non-terminal
    for (Module *curNode : dbNodes)
    {
        double w = curNode->getWidth();
        double h = curNode->getHeight();

        fprintf(out, " %30s %10.0f %10.0f\n",
                curNode->name.c_str(),
                w,
                h);
    }

    // terminal
    for (Module *curTerminal : dbTerminals)
    {
        double w = curTerminal->getWidth();
        double h = curTerminal->getHeight();

        //! isNI 的 terminal 写作 terminal_NI：告诉下游合法器它只是占位，
        //! 不参与密度/线长统计，也不要当成真正的 IO pad 去对齐。
        if (curTerminal->isNI)
        {
            fprintf(out, " %10s %10.0f %10.0f terminal_NI\n",
                    curTerminal->name.c_str(),
                    w,
                    h);
        }
        else
        {
            fprintf(out, " %10s %10.0f %10.0f terminal\n",
                    curTerminal->name.c_str(),
                    w,
                    h);
        }
    }
    fprintf(out, "\n\n");
    fclose(out);
}

//! 写 .pl：所有 module 的**左下角坐标 + 朝向**（terminal 额外标 /FIXED 或 /FIXED_NI）
//! 这是唯一一个每次 outputBookShelf 都会写的文件（plOnly=true 时只写它）
void PlaceDB::outputPL()
{
    string outputFilePath;
    gArg.GetString("outputPath", &outputFilePath);

    string benchmarkName;
    gArg.GetString("benchmarkName", &benchmarkName);

    string suffix;
    gArg.GetString("outputSuffix", &suffix);

    outputFilePath += benchmarkName;
    outputFilePath += "-" + suffix + ".pl";

    //! 同样把 .pl 的路径回写 gArg（ePlaceAbacus 会读它把上一阶段坐标读回来）
    gArg.Override("outputPL", outputFilePath);

    cout << "Output PL file:" << outputFilePath << endl;

    // out "pl"
    FILE *out = fopen(outputFilePath.c_str(), "w");
    //! 疑似问题：这里 fopen 之后**没有**判空，与 outputNodes/Nets/SCL 不一致；
    //!   目录不存在时会直接往 NULL 里 fprintf 而崩。
    fprintf(out, "UCLA pl 1.0\n\n");

    char *orientN = "N";
    //! 注意：这里把所有单元的朝向一律写成 "N"，没有真实的朝向信息被输出
    //!   （Module::orientation 在本流程里基本没被维护）。

    for (Module *curNode : dbNodes)
    {
        fprintf(out, "%s\t%.0f\t%.0f : %s",
                curNode->name.c_str(),
                curNode->getLL_2D().x,
                curNode->getLL_2D().y,
                orientN);
        fprintf(out, "\n");
    }
    for (Module *curTerminal : dbTerminals)
    {
        fprintf(out, "%s\t%.0f\t%.0f : %s",
                curTerminal->name.c_str(),
                curTerminal->getLL_2D().x,
                curTerminal->getLL_2D().y,
                orientN);

        if (curTerminal->isNI)
            fprintf(out, " /FIXED_NI\n");
        else
            fprintf(out, " /FIXED\n");
    }
    fprintf(out, "\n\n");

    fclose(out);
}

//! 写 .nets：逐 net 写 NetDegree，再逐 pin 写「所属 module 名 + 方向 B + pin 相对偏移」
void PlaceDB::outputNets()
{
    string outputFilePath;
    gArg.GetString("outputPath", &outputFilePath);

    string benchmarkName;
    gArg.GetString("benchmarkName", &benchmarkName);

    string suffix;
    gArg.GetString("outputSuffix", &suffix);

    outputFilePath += benchmarkName;
    outputFilePath += "-" + suffix + ".nets";

    cout << "Output Nets file:" << outputFilePath << endl;

    FILE *out;
    out = fopen(outputFilePath.c_str(), "w");
    if (!out)
    {
        cerr << "Cannot open output file\n";
        return;
    }

    //! 写 .nets：逐 net 写 NetDegree，再逐 pin 写「所属 module 名 + 方向 + pin 相对偏移」。
    //! 注意写的是 pin->offset（相对 module 中心），不是绝对坐标。
    fprintf(out, "UCLA nets 1.0\n\n");
    fprintf(out, "NumNets : %d\n", dbNets.size());
    fprintf(out, "NumPins : %d\n", dbPins.size());
    //! 疑似问题：NumPins 用的是 dbPins.size()，但 Bookshelf 的 .nets 里 NumPins 指的
    //!   是「各 net 的 pin 条目总数」。由于一个物理 pin 在 dbPins 里只存一份，
    //!   两者在普通设计里相等；若存在被多个 net 共享的 pin 记录就会不一致。
    for (Net *curNet : dbNets)
    {
        fprintf(out, "NetDegree : %d\n", curNet->netPins.size());
        for (Pin *curPin : curNet->netPins)
        {
            fprintf(out, " %10s B : %.2f %.2f\n",

                    curPin->module->name.c_str(),
                    curPin->offset.x,
                    curPin->offset.y);
        }
    }
    fprintf(out, "\n");
    fclose(out);
}

//! 写 .scl：布局行信息（行的 y 坐标、行高、site 宽/间距、起始 x 与 site 数）
void PlaceDB::outputSCL()
{
    string outputFilePath;
    gArg.GetString("outputPath", &outputFilePath);

    string benchmarkName;
    gArg.GetString("benchmarkName", &benchmarkName);

    string suffix;
    gArg.GetString("outputSuffix", &suffix);

    outputFilePath += benchmarkName;
    outputFilePath += "-" + suffix + ".scl";

    cout << "Output SCL file:" << outputFilePath << endl;

    FILE *out;
    out = fopen(outputFilePath.c_str(), "w");
    if (!out)
    {
        cerr << "Cannot open output file\n";
        return;
    }

    fprintf(out, "UCLA scl 1.0\n");
    fprintf(out, "# Created       :\n");
    fprintf(out, "# User          :\n\n");
    fprintf(out, "NumRows : %d\n\n", dbSiteRows.size());

    char *ori[2] = {"N", "Y"};
    //! 疑似问题：ori 是**死变量**，从未被使用（下面的 Siteorient 被写死成 1）。
    //!   另外 "Y" 并不是合法的 Bookshelf 朝向值（应为 "N"/"S"/"W"/"E"/"FN" 等，
    //!   见 global.h 的 orientInt）。

    //! 每行输出一条 CoreRow。注意这里写的是行的**原始** start/end，
    //! 并没有把 removeBlockedSite() 算出来的 intervals 拆成多条 Subrow ——
    //! 也就是说被 macro 挖掉的区间信息不会体现在输出的 .scl 里。
    for (SiteRow curRow : dbSiteRows)
    {
        double step = curRow.step;
        if (step == 0)
            step = 1.0;
        fprintf(out, "CoreRow Horizontal\n");
        fprintf(out, " Coordinate    : %8.0f\n", curRow.bottom);
        fprintf(out, " Height        : %8.0f\n", curRow.height);
        fprintf(out, " Sitewidth     : %8.0f\n", step);
        fprintf(out, " Sitespacing   : %8.0f\n", step);
        fprintf(out, " Siteorient    : 1\n"); //%s\n", ori[i % 2] );
        fprintf(out, " Sitesymmetry  : 1\n");
        fprintf(out, " SubrowOrigin  : %8.0f Numsites : %8.0f\n",
                curRow.start.x, (curRow.end.x - curRow.start.x) / curRow.step);
        //! 疑似问题：上面 Numsites 用的是**未做 0 保护**的 curRow.step，
        //!   而 Sitewidth/Sitespacing 用的是修正后的局部 step。step == 0 时这里会除零。

        fprintf(out, "End\n");
    }
    fprintf(out, "\n");
    fclose(out);
}

// ----------------------------------------------------------------------------
// 下面整段是已被注释掉的 plotCurrentPlacement()：把当前布局画成 BMP。
// 它的职责已经搬到 Plot/plot.cpp（PLOTTING::plotCurrentPlacement），这里保留作参考。
// 画图的坐标基准是 chipRegion 而非 coreRegion —— 因为要把 core 外的 terminal 也画进去。
// ----------------------------------------------------------------------------
// void PlaceDB::plotCurrentPlacement(string imageName)
// {
//     string plotPath;
//     if (!gArg.GetString("plotPath", &plotPath))
//     {
//         plotPath = "./";
//     }

//     float chipRegionWidth = chipRegion.ur.x - chipRegion.ll.x;
//     float chipRegionHeight = chipRegion.ur.y - chipRegion.ll.y;

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

//     for (Module *curTerminal : dbTerminals)
//     {
//         assert(curTerminal);
//         // ignore pin's location
//         if (curTerminal->isNI)
//         {
//             continue;
//         }
//         int x1 = getX(chipRegion.ll.x, curTerminal->getLL_2D().x, unitX) + xMargin;
//         int x2 = getX(chipRegion.ll.x, curTerminal->getUR_2D().x, unitX) + xMargin;
//         int y1 = getY(chipRegionHeight, chipRegion.ll.y, curTerminal->getLL_2D().y, unitY) + yMargin;
//         int y2 = getY(chipRegionHeight, chipRegion.ll.y, curTerminal->getUR_2D().y, unitY) + yMargin;
//         img.draw_rectangle(x1, y1, x2, y2, Blue, opacity);
//     }

//     for (Module *curNode : dbNodes)
//     {
//         assert(curNode);
//         int x1 = getX(chipRegion.ll.x, curNode->getLL_2D().x, unitX) + xMargin;
//         int x2 = getX(chipRegion.ll.x, curNode->getUR_2D().x, unitX) + xMargin;
//         int y1 = getY(chipRegionHeight, chipRegion.ll.y, curNode->getLL_2D().y, unitY) + yMargin;
//         int y2 = getY(chipRegionHeight, chipRegion.ll.y, curNode->getUR_2D().y, unitY) + yMargin;
//         if (curNode->isMacro)
//         {
//             img.draw_rectangle(x1, y1, x2, y2, Orange, opacity);
//         }
//         else
//         {
//             img.draw_rectangle(x1, y1, x2, y2, Red, opacity);
//         }
//     }

//     img.draw_text(50, 50, imageName.c_str(), Black, NULL, 1, 30);
//     img.save_bmp(string(plotPath + imageName + string(".bmp")).c_str());
//     cout << "INFO: BMP HAS BEEN SAVED: " << imageName + string(".bmp") << endl;
// }

// ----------------------------------------------------------------------------
// 给所有可动单元的中心加一点随机扰动（-addNoise 时调用，在 mGP 之前）
//
// 目的：QPlace 出来的初始解里，大量单元会精确重合在同一个点上（见
// moveNodesCenterToCenter）。完全重合的单元在密度场里产生的力是对称的，
// ePlace 很难把它们分开（梯度互相抵消）。加一点噪声打破这种对称性。
//
// 扰动幅度：±(0.5 × 2.5%) × 单元尺寸，即约 ±1.25% 的宽/高 —— 刻意取得很小，
// 既能打破对称又不至于破坏初始解的线长质量。随机数用固定种子 0，便于复现调试。
// ----------------------------------------------------------------------------
void PlaceDB::addNoise()
{
    static std::mt19937 gen(0); // 0:random seed, fixed for debugging

    for (auto node : dbNodes)
    {
        VECTOR_2D range;
        range.x = 0.5 * 0.025 * node->getWidth();
        range.y = 0.5 * 0.025 * node->getHeight();
        //! 疑似问题：调用方注释写的是 "the noise range is [-avgbinStep, avgbinStep]"
        //!   （见 main/ePlace_main.cpp），但实际幅度是 ±1.25% 的**单元尺寸**，
        //!   与 bin 步长毫无关系，注释与实现不符。

        std::uniform_real_distribution<float> disx(-range.x, range.x);
        std::uniform_real_distribution<float> disy(-range.y, range.y);

        POS_3D pos = node->center;
        pos.x += disx(gen);
        pos.y += disy(gen);
        node->setCenter_2D(pos.x, pos.y);
        //! 疑似问题：这里**绕过**了 setModuleCenter_2D，直接调 module->setCenter_2D，
        //!   因此该单元上所有 pin 的 absolutePos 不会被刷新，线长/梯度会读到扰动前的
        //!   旧坐标（直到下一次经由 setModuleCenter* 的移动才修正）。
        //!   这正是文件头【主线一】强调"写坐标必须走 setModuleCenter*/setModuleLocation*"
        //!   的原因。另外这里也没有做 coreRegion 裁剪。
    }
}

//! 把所有可动单元的**左下角**坐标存进 nodesLocationRegister（快照/备份）
void PlaceDB::saveNodesLocation()
{
    int nodesCount = dbNodes.size();
    nodesLocationRegister.resize(nodesCount);

    for (int i = 0; i < nodesCount; i++)
    {
        nodesLocationRegister[i] = dbNodes[i]->getLocation();
    }
}

//! 用 saveNodesLocation() 存下的快照恢复所有可动单元的位置（回滚）
void PlaceDB::loadNodesLocation()
{
    int nodesCount = dbNodes.size();
    assert(nodesLocationRegister.size() == nodesCount);
    //! 注意 assert 在 NDEBUG 下会被编译掉，Release 构建里这个长度校验形同虚设；
    //!   而且若从未调用过 saveNodesLocation()，这里会直接断言失败。
    for (int i = 0; i < nodesCount; i++)
    {
        setModuleLocation_2D(dbNodes[i], nodesLocationRegister[i]);
    }
}

// ============================================================================
// getOptimialRegion() —— 求一个单元的「最优区域」（optimal region），详细布局用
//
// 出处：An Efficient and Effective Detailed Placement Algorithm（ISPD 的 DP 论文）。
// 思想：对单元参与的每条 net，先**排除本单元自己的 pin**，求其余 pin 的包围盒；
// 把这些包围盒的边界坐标全部收集起来，取中位数作为最优区域的四条边。
// 直观含义：「把这个单元挪到哪个矩形里，它所连的 net 最可能变短」——
// 于是详细布局只需要在这一小块区域里找空位，而不必搜索整个 core。
//
// 为什么取中位数而不是平均值：中位数对个别离得很远的 pin 不敏感，
// 得到的区域不会被一条长 net 拉得过大。
//
// 注意：返回值是**未**对齐到 site / 行高的原始矩形；对齐由调用方（detailed.cpp）
// 自行用 floor/ceil 处理（见那里对 searchRegion 的二次修正）。
// ============================================================================
CRect PlaceDB::getOptimialRegion(Module *module)
{
    // see the description of the 'optimal region' in the paper: An efficient and effective detailed placement algorithm
    //? Return a rectangle which has integer width and height. And its height is an integer multiple of the row height while its width align with row site
    //? is it necessary to align?
    vector<float> Xs;
    vector<float> Ys;

    for (Net *curNet : module->nets)
    {
        double maxX = -DOUBLE_MAX;
        double minX = DOUBLE_MAX;

        double maxY = -DOUBLE_MAX;
        double minY = DOUBLE_MAX;

        double curX;
        double curY;

        POS_3D curPos;
        double HPWL;
        //! 疑似问题：HPWL 是**死变量**，声明后从未被赋值也从未被使用（历史遗留）。

        for (Pin *curPin : curNet->netPins)
        {
            if (curPin->module == module)
            {
                continue;
            }
            curPos = curPin->absolutePos;
            //! 同样直接读缓存的 absolutePos，调用前需保证 pin 坐标是最新的。
            curX = curPos.x;
            curY = curPos.y;

            minX = min(minX, curX);
            maxX = max(maxX, curX);
            minY = min(minY, curY);
            maxY = max(maxY, curY);
        }
        //! 每条 net 贡献两个 x（min/max）和两个 y，所以 Xs / Ys 的长度一定是偶数
        //! —— 这正是下面能直接取「中间两个元素当中位数」的前提。
        //!
        //! 疑似问题：若某条 net 上除本 module 外**没有**其它 pin（net 的所有 pin 都在
        //!   这个 module 上，例如 1 度 net），minX/maxX 会停留在 ±DOUBLE_MAX 初值，
        //!   把垃圾值 push 进 Xs/Ys，中位数随之失真。
        Xs.push_back(minX);
        Xs.push_back(maxX);
        Ys.push_back(minY);
        Ys.push_back(maxY);
    }

    // find medians of Xs and Ys.
    // ! The size of Xs and Ys should be even.
    float left, right, bottom, up;

    // left = getKth(Xs, Xs.size() / 2 - 1); // left and right boundary of the optimal region
    // right = getKth(Xs, Xs.size() / 2);

    // bottom = getKth(Ys, Ys.size() / 2 - 1);
    // up = getKth(Ys, Ys.size() / 2);

    // or

    //! 上面被注释掉的 getKth 是 O(n) 的快速选择；这里改用完整 sort（O(n log n)）。
    //! 因为一个单元的 net 数通常很小（几条到几十条），两者差别不大，sort 更好维护。
    //! 疑似问题：getKth 的 partion() 会把 float 枢轴截成 int（见 global.h 的标注），
    //!   即便启用被注释的那条路径，结果也不可靠。
    sort(Xs.begin(), Xs.end());
    sort(Ys.begin(), Ys.end());

    left = Xs[Xs.size() / 2 - 1]; // left and right boundary of the optimal region
    right = Xs[Xs.size() / 2];

    bottom = Ys[Ys.size() / 2 - 1];
    up = Ys[Ys.size() / 2];
    //! 疑似问题：module 没有连任何 net（或 Xs 为空）时，Xs.size()/2 - 1 会下溢成
    //!   巨大的 size_t 值，越界访问。调用方需自行保证 module 至少有一条 net。

    // create CRect and return

    CRect result;
    result.ll.x = left;
    result.ll.y = bottom;
    result.ur.x = right;
    result.ur.y = up;
    //! 疑似问题：当 left > right 或 bottom > up（例如所有 net 都是 2 pin 且退化成一个点时）
    //!   会返回一个"反向"矩形，CRect::getWidth() 里的 assert(width > 0) 会炸。
    return result;
}

// ----------------------------------------------------------------------------
// y 坐标 → 所在 site row 的下标
//
// 直接按「离 coreRegion 下边界有几个行高」整除得到，是 O(1) 的。
// 前提（!!）：coreRegion.ll.y 必须正好等于第一行的 bottom，且所有行等高、行间无空隙
// —— 这由 setCoreRegion() 的推导方式保证（见那里的注释）。
// ----------------------------------------------------------------------------
int PlaceDB::y2RowIndex(float y)
{
    int index = (int)((y - coreRegion.ll.y) / commonRowHeight); //!!!!!!! assume coreRegion.ll.y == the bottom of the first row, check setCoreRegion()

    assert(index >= 0);
    assert(index < dbSiteRows.size());
    //! 疑似问题：(int) 是向零截断，所以 y 略小于 coreRegion.ll.y 时会得到 0 而不是 -1，
    //!   index >= 0 这条断言抓不住这种"擦边"情况；而 y 等于或超过 coreRegion.ur.y 时
    //!   会算出等于行数的下标，触发第二条断言（Release 下则是越界访问）。
    //!   调用方（detailed.cpp）需要自己先把 y 夹到 [ll.y, ur.y - rowHeight] 内。

    return index;
}

// ----------------------------------------------------------------------------
// 判断两个 module 是否共用至少一条 net（即"相连"）
//
// 详细布局的 independentCells 用它筛掉会互相影响的单元：只有互不相连的单元才能
// 在同一个窗口里独立地挪动而不互相干扰线长评估。
// ----------------------------------------------------------------------------
bool PlaceDB::isConnected(Module *module1, Module *module2)
{
    for (Net *module1Net : module1->nets)
    {
        for (Net *module2Net : module2->nets)
        {
            if (module1Net == module2Net)
            {
                return true;
            }
        }
    }
    //! 疑似问题：复杂度是 O(度数1 × 度数2) 的朴素两两比较，高扇出单元上会比较慢；
    //!   另外没有处理 module1 == module2 的情况（自己和自己必然"相连"，返回 true）。
    return false;
}
