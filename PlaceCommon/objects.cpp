#include <objects.h>

// ============================================================================
// objects.cpp —— Net / Pin / Module 等原子对象的实现
//
// 这里的核心是「线长模型」的三种实现（都在 Net 上）：
//   · calcNetHPWL      —— 精确半周长线长。遍历所有 pin 求 X/Y/Z 的极值，不可微，只用于评估打印
//   · calcBoundPin     —— 同样是 HPWL，但顺带把边界 pin 记到 boundPin* 上，供后续复用
//   · calcWirelength*  —— WA / LSE 两种**平滑可微**的近似线长，梯度下降实际用的就是它们
//
// 为什么需要平滑模型：HPWL 是 max/min 函数，在极值切换处不可导，
// 用带参数 γ 的 log-sum-exp（LSE）或加权平均（WA）来逼近，
// γ 越小越接近真实 HPWL 但越「陡」，γ 越大越平滑（见 EPlace 里 γ 随溢出率 τ 变化）。
// 注意：代码里传的都是 invertedGamma = 1/γ。
// ============================================================================

void Net::addPin(Pin *pin)
{
    netPins.push_back(pin);
}

int Net::getPinCount()
{
    return netPins.size();
}

void Net::allocateMemoryForPin(int n)
{
    netPins.reserve(n); //! 预分配，避免 push_back 反复扩容（解析大设计时影响明显）
}

//! 精确 HPWL = (maxX−minX) + (maxY−minY) + (maxZ−minZ)
//! 注意：这里直接读 curPin->absolutePos 缓存值（而不是调 getAbsolutePos），
//! 所以调用前必须保证所有 pin 的绝对坐标已刷新，否则算出来是旧值
double Net::calcNetHPWL()
{
    double maxX = -DOUBLE_MAX;
    double minX = DOUBLE_MAX;
    // double maxY = DOUBLE_MIN;
    double maxY = -DOUBLE_MAX;
    double minY = DOUBLE_MAX;
    // double maxZ = DOUBLE_MIN;
    double maxZ = -DOUBLE_MAX; // potential bug: double_min >0 so boundPinZmax might be null when all z == 0
    double minZ = DOUBLE_MAX;

    double curX;
    double curY;
    double curZ;
    POS_3D curPos;
    double HPWL;
    for (Pin *curPin : netPins)
    {
        // curPos = curPin->getAbsolutePos();
        curPos=curPin->absolutePos;
        curX = curPos.x;
        curY = curPos.y;
        curZ = curPos.z;
        minX = min(minX, curX);
        maxX = max(maxX, curX);
        minY = min(minY, curY);
        maxY = max(maxY, curY);
        minZ = min(minZ, curZ);
        maxZ = max(maxZ, curZ);
    }
    if (!gArg.CheckExist("3DIC"))
    {
        //? assert(maxZ == minZ == 0); this causes bug
        assert(float_equal(maxZ, 0.0));
        assert(float_equal(minZ, 0.0));
    }
    HPWL = ((maxX - minX) + (maxY - minY) + (maxZ - minZ));
    return HPWL;
}

double Net::calcBoundPin()
{
    // double maxX = DOUBLE_MIN;
    double maxX = -DOUBLE_MAX;
    double minX = DOUBLE_MAX;
    // double maxY = DOUBLE_MIN;
    double maxY = -DOUBLE_MAX;
    double minY = DOUBLE_MAX;
    // double maxZ = DOUBLE_MIN;
    double maxZ = -DOUBLE_MAX; // potential bug: double_min >0 so boundPinZmax might be null when all z == 0
    double minZ = DOUBLE_MAX;

    double curX;
    double curY;
    double curZ;
    POS_3D curPos;
    double HPWL;

    for (Pin *curPin : netPins)
    {
        curPos = curPin->getAbsolutePos();
        curX = curPos.x;
        curY = curPos.y;
        curZ = curPos.z;
        //!!!!! assume curX curY curZ always >= 0!!!
        assert(curZ == 0);
        if (curX < minX)
        {
            minX = curX;
            boundPinXmin = curPin;
        }

        if (curX > maxX)
        {
            maxX = curX;
            boundPinXmax = curPin;
        }

        if (curY < minY)
        {
            minY = curY;
            boundPinYmin = curPin;
        }

        if (curY > maxY)
        {
            maxY = curY;
            boundPinYmax = curPin;
        }

        if (curZ < minZ)
        {
            minZ = curZ;
            boundPinZmin = curPin;
        }

        if (curZ > maxZ)
        {
            maxZ = curZ;
            boundPinZmax = curPin;
        }
    }
    if (!gArg.CheckExist("3DIC"))
    {
        assert(float_equal(maxZ, 0.0));
        assert(float_equal(minZ, 0.0));
    }
    HPWL = ((maxX - minX) + (maxY - minY) + (maxZ - minZ));
    return HPWL;
}

void Net::clearBoundPins()
{
    boundPinXmax = NULL;
    boundPinXmin = NULL;
    boundPinYmax = NULL;
    boundPinYmin = NULL;
    boundPinZmax = NULL;
    boundPinZmin = NULL;
}

//! LSE（log-sum-exp）平滑线长：
//!   WL ≈ (Xmax − Xmin) + (1/α)·ln Σ e^[α(Xi−Xmax)] + (1/α)·ln Σ e^[α(Xmin−Xi)]
//! 其中 α = invertedGamma = 1/γ。α 越大越贴近 HPWL。
//! 每个 pin 的 e^[] 值和 sumMax/sumMin 都被缓存下来，梯度函数直接复用，避免重复计算。
//! expZeroFlg*：指数太小（< NEGATIVE_MAX_EXP）时 e^[] 下溢为 0，打标记让梯度端跳过该项。
double Net::calcWirelengthLSE_2D(VECTOR_2D invertedGamma)
{
    VECTOR_2D sumMax;
    VECTOR_2D sumMin;
    sumMax.SetZero();
    sumMin.SetZero();

    assert(boundPinXmax);
    assert(boundPinXmin);
    assert(boundPinYmax);
    assert(boundPinYmin);

    float pinMaxX = boundPinXmax->getAbsolutePos().x;
    float pinMaxY = boundPinYmax->getAbsolutePos().y;

    float pinMinX = boundPinXmin->getAbsolutePos().x;
    float pinMinY = boundPinYmin->getAbsolutePos().y;

    for (Pin *curPin : netPins)
    {
        POS_3D pinPosition = curPin->getAbsolutePos();
        VECTOR_2D expMax;
        VECTOR_2D expMin;
        expMax.x = (pinPosition.x - pinMaxX) * invertedGamma.x;
        expMin.x = (pinMinX - pinPosition.x) * invertedGamma.x;
        expMax.y = (pinPosition.y - pinMaxY) * invertedGamma.y;
        expMin.y = (pinMinY - pinPosition.y) * invertedGamma.y;

        if (expMax.x > NEGATIVE_MAX_EXP)
        {
            curPin->eMax_LSE.x = fastExp(expMax.x);
            sumMax.x += curPin->eMax_LSE.x;
            curPin->expZeroFlgMax_LSE.x = false;
        }
        else
        {
            curPin->expZeroFlgMax_LSE.x = true;
        }

        if (expMin.x > NEGATIVE_MAX_EXP)
        {
            curPin->eMin_LSE.x = fastExp(expMin.x);
            sumMin.x += curPin->eMin_LSE.x;
            curPin->expZeroFlgMin_LSE.x = false;
        }
        else
        {
            curPin->expZeroFlgMin_LSE.x = true;
        }

        if (expMax.y > NEGATIVE_MAX_EXP)
        {
            curPin->eMax_LSE.y = fastExp(expMax.y);
            sumMax.y += curPin->eMax_LSE.y;
            curPin->expZeroFlgMax_LSE.y = false;
        }
        else
        {
            curPin->expZeroFlgMax_LSE.y = true;
        }

        if (expMin.y > NEGATIVE_MAX_EXP)
        {
            curPin->eMin_LSE.y = fastExp(expMin.y);
            sumMin.y += curPin->eMin_LSE.y;
            curPin->expZeroFlgMin_LSE.y = false;
        }
        else
        {
            curPin->expZeroFlgMin_LSE.y = true;
        }
    }

    sumMax_LSE.x = sumMax.x;
    sumMax_LSE.y = sumMax.y;
    sumMin_LSE.x = sumMin.x;
    sumMin_LSE.y = sumMin.y;

    return (pinMaxX - pinMinX + log(sumMax.x) / invertedGamma.x + log(sumMin.x) / invertedGamma.x) +
           (pinMaxY - pinMinY + log(sumMax.y) / invertedGamma.y + log(sumMin.y) / invertedGamma.y);
}

//! WA（weighted-average，加权平均）平滑线长，见 NTUPlace3D 论文第 6 页：
//!   WL ≈ Σ(Xi·e^[α(Xi−Xmax)]) / Σ(e^[α(Xi−Xmax)])  −  Σ(Xi·e^[α(Xmin−Xi)]) / Σ(e^[α(Xmin−Xi)])
//! 即「用指数权重做加权平均」来替代 max 和 min，α = invertedGamma = 1/γ。
//! 相比 LSE 更稳定（分子分母都是有界加权平均），是本工程的默认模型。
//! 同样依赖 boundPin*，调用前必须先跑过 calcBoundPin()。
double Net::calcWirelengthWA_2D(VECTOR_2D invertedGamma)
{
    VECTOR_2D numeratorMax;
    VECTOR_2D denominatorMax;
    VECTOR_2D numeratorMin;
    VECTOR_2D denominatorMin;

    //! WA wirelength model, see NTUPlace 3D paper page 6 : Stable Weighted-Average Wirelength Model
    //! Here on X/Y dimension: WA wirelength = numeratorMax/denominatorMax - numeratorMin/denominatorMin, total wirelength = wirelength on X dimension + wirelength on Y dimension
    //! numerator and denominator are sum of the results of all pins, see the code below

    numeratorMax.SetZero();
    denominatorMax.SetZero();
    numeratorMin.SetZero();
    denominatorMin.SetZero();

    assert(boundPinXmax);
    assert(boundPinXmin);
    assert(boundPinYmax);
    assert(boundPinYmin);

    float pinMaxX = boundPinXmax->getAbsolutePos().x;
    float pinMaxY = boundPinYmax->getAbsolutePos().y;

    float pinMinX = boundPinXmin->getAbsolutePos().x;
    float pinMinY = boundPinYmin->getAbsolutePos().y;

    for (Pin *curPin : netPins)
    {
        POS_3D pinPosition = curPin->getAbsolutePos();
        VECTOR_2D expMax;                                       // (Xi-Xmax)/gamma in WA model (X/Y/Z)
        VECTOR_2D expMin;                                       // (Xmin-Xi)/gamma in WA model (X/Y/Z)
        expMax.x = (pinPosition.x - pinMaxX) * invertedGamma.x; //! wlen_cof is actually 1/gamma
        expMin.x = (pinMinX - pinPosition.x) * invertedGamma.x; //! wlen_cof used here!
        expMax.y = (pinPosition.y - pinMaxY) * invertedGamma.y;
        expMin.y = (pinMinY - pinPosition.y) * invertedGamma.y;
        // cout<<padding<<"expmax: "<<exp_max_x<<endl;
        if (expMax.x > NEGATIVE_MAX_EXP)
        {
            curPin->eMax_WA.x = fastExp(expMax.x);
            numeratorMax.x += pinPosition.x * curPin->eMax_WA.x;
            denominatorMax.x += curPin->eMax_WA.x;
            curPin->expZeroFlgMax_WA.x = false;
        }
        else
        {
            curPin->expZeroFlgMax_WA.x = true;
        }

        if (expMin.x > NEGATIVE_MAX_EXP)
        {
            curPin->eMin_WA.x = fastExp(expMin.x);
            numeratorMin.x += pinPosition.x * curPin->eMin_WA.x;
            denominatorMin.x += curPin->eMin_WA.x;
            curPin->expZeroFlgMin_WA.x = false;
        }
        else
        {
            curPin->expZeroFlgMin_WA.x = true;
        }

        if (expMax.y > NEGATIVE_MAX_EXP)
        {
            curPin->eMax_WA.y = fastExp(expMax.y);
            numeratorMax.y += pinPosition.y * curPin->eMax_WA.y;
            denominatorMax.y += curPin->eMax_WA.y;
            curPin->expZeroFlgMax_WA.y = false;
        }
        else
        {
            curPin->expZeroFlgMax_WA.y = true;
        }

        if (expMin.y > NEGATIVE_MAX_EXP)
        {
            curPin->eMin_WA.y = fastExp(expMin.y);
            numeratorMin.y += pinPosition.y * curPin->eMin_WA.y;
            denominatorMin.y += curPin->eMin_WA.y;
            curPin->expZeroFlgMin_WA.y = false;
        }
        else
        {
            curPin->expZeroFlgMin_WA.y = true;
        }
    }

    numeratorMax_WA.x = numeratorMax.x;
    numeratorMax_WA.y = numeratorMax.y;
    denominatorMax_WA.x = denominatorMax.x;
    denominatorMax_WA.y = denominatorMax.y;

    numeratorMin_WA.x = numeratorMin.x;
    numeratorMin_WA.y = numeratorMin.y;
    denominatorMin_WA.x = denominatorMin.x;
    denominatorMin_WA.y = denominatorMin.y;

    return (numeratorMax_WA.x / denominatorMax_WA.x - numeratorMin_WA.x / denominatorMin_WA.x) + (numeratorMax_WA.y / denominatorMax_WA.y - numeratorMin_WA.y / denominatorMin_WA.y);
}

// ----------------------------------------------------------------------------
// WA 模型对某个 pin 坐标的偏导
//
// WA 线长本身是「分式」：   f = numerator / denominator
// 所以对 Xi 求导要用商法则：∂f/∂Xi = (num'·den − den'·num) / den²
//   · 对 max 那一支：num' = e^[] + Xi·(α·e^[])，den' = α·e^[]
//   · 对 min 那一支：指数反向，所以 den' 的符号相反，商法则里变成 + 号
// 最后梯度 = ∂(max支)/∂Xi − ∂(min支)/∂Xi。
// 指数下溢（expZeroFlg）的项直接当 0 跳过。
// 注意：numerator/denominator 用的是 Net 上缓存的 *_WA 成员，
// 所以必须先调过 calcWirelengthWA_2D 才有正确值。
// ----------------------------------------------------------------------------
VECTOR_2D Net::getWirelengthGradientWA_2D(VECTOR_2D invertedGamma, Pin *curPin)
{
    assert(curPin);
    VECTOR_2D gradientOnCurrentPin = VECTOR_2D();
    VECTOR_2D gradientNumeratorMax = VECTOR_2D();
    VECTOR_2D gradientDenominatorMax = VECTOR_2D();
    VECTOR_2D gradientNumeratorMin = VECTOR_2D();
    VECTOR_2D gradientDenominatorMin = VECTOR_2D();
    VECTOR_2D gradientMax = VECTOR_2D();
    VECTOR_2D gradientMin = VECTOR_2D();
    // ? no SetZero here (called in default constructor)
    //? assert(gradientOnCurrentPin.x == gradientDenominatorMin.y == 0.0);
    assert(gradientOnCurrentPin.x == 0.0);
    assert(gradientDenominatorMin.y == 0.0);

    POS_3D curPinPosition = curPin->getAbsolutePos();

    if (!curPin->expZeroFlgMax_WA.x)
    { // if flg=0, assume grad=0
        gradientDenominatorMax.x = invertedGamma.x * curPin->eMax_WA.x;
        gradientNumeratorMax.x = curPin->eMax_WA.x + curPinPosition.x * gradientDenominatorMax.x;
        gradientMax.x =
            (gradientNumeratorMax.x * denominatorMax_WA.x - gradientDenominatorMax.x * numeratorMax_WA.x) /
            (denominatorMax_WA.x * denominatorMax_WA.x);
    }

    if (!curPin->expZeroFlgMax_WA.y)
    {
        gradientDenominatorMax.y = invertedGamma.y * curPin->eMax_WA.y;
        gradientNumeratorMax.y = curPin->eMax_WA.y + curPinPosition.y * gradientDenominatorMax.y;
        gradientMax.y =
            (gradientNumeratorMax.y * denominatorMax_WA.y - gradientDenominatorMax.y * numeratorMax_WA.y) /
            (denominatorMax_WA.y * denominatorMax_WA.y);
    }

    if (!curPin->expZeroFlgMin_WA.x)
    {
        gradientDenominatorMin.x = invertedGamma.x * curPin->eMin_WA.x;
        gradientNumeratorMin.x = curPin->eMin_WA.x - curPinPosition.x * gradientDenominatorMin.x;
        gradientMin.x =
            (gradientNumeratorMin.x * denominatorMin_WA.x + gradientDenominatorMin.x * numeratorMin_WA.x) /
            (denominatorMin_WA.x * denominatorMin_WA.x);
    }

    if (!curPin->expZeroFlgMin_WA.y)
    {
        gradientDenominatorMin.y = invertedGamma.y * curPin->eMin_WA.y;
        gradientNumeratorMin.y = curPin->eMin_WA.y - curPinPosition.y * gradientDenominatorMin.y;
        gradientMin.y =
            (gradientNumeratorMin.y * denominatorMin_WA.y + gradientDenominatorMin.y * numeratorMin_WA.y) /
            (denominatorMin_WA.y * denominatorMin_WA.y);
    }

    gradientOnCurrentPin.x = gradientMax.x - gradientMin.x;
    gradientOnCurrentPin.y = gradientMax.y - gradientMin.y;
    return gradientOnCurrentPin;
}

// ----------------------------------------------------------------------------
// LSE 模型对某个 pin 坐标的偏导（比 WA 简洁得多）
//   ∂/∂Xi [ (1/α)·ln Σ e^[α(Xi−Xmax)] ] = e^[α(Xi−Xmax)] / Σ e^[α(Xi−Xmax)]
// 即「自身指数项占总和的比例」，物理含义是 softmax 权重。
// min 那一支同理，最终梯度 = max支权重 − min支权重。
// ----------------------------------------------------------------------------
VECTOR_2D Net::getWirelengthGradientLSE_2D(VECTOR_2D invertedGamma, Pin *curPin)
{
    VECTOR_2D gradientOnCurrentPin = VECTOR_2D();
    VECTOR_2D gradientMax = VECTOR_2D(); // the gradient added by positive term
    VECTOR_2D gradientMin = VECTOR_2D(); // the gradient added by negative term

    gradientMax.x = (curPin->expZeroFlgMax_LSE.x ? 0 : curPin->eMax_LSE.x) / sumMax_LSE.x;
    gradientMin.x = (curPin->expZeroFlgMin_LSE.x ? 0 : curPin->eMin_LSE.x) / sumMin_LSE.x;
    gradientMax.y = (curPin->expZeroFlgMax_LSE.y ? 0 : curPin->eMax_LSE.y) / sumMax_LSE.y;
    gradientMin.y = (curPin->expZeroFlgMin_LSE.y ? 0 : curPin->eMin_LSE.y) / sumMin_LSE.y;

    gradientOnCurrentPin.x = gradientMax.x - gradientMin.x;
    gradientOnCurrentPin.y = gradientMax.y - gradientMin.y;
    return gradientOnCurrentPin;
}

//! 直接返回缓存的绝对坐标（不重算）。调用前需确保 module 移动后已触发 calculateAbsolutePos
POS_3D Pin::getAbsolutePos()
{
    // POS_3D absPos;
    // module->calcCenter();//?
    return absolutePos;
}

// POS_3D Pin::fetchAbsolutePos()
// {
//     return absolutePos;
// }

//! 绝对坐标 = 所属 module 的中心 + pin 的 offset。
//! z 分量直接取 module 的 z（offset 只有 x/y，因为 pin 都在单元表面同一层）
void Pin::calculateAbsolutePos()
{
    absolutePos.x = module->getCenter().x + offset.x;
    absolutePos.y = module->getCenter().y + offset.y;
    absolutePos.z = module->getCenter().z;
}

void Pin::setId(int index)
{
    idx = index;
}

void Pin::setNet(Net *_net)
{
    net = _net;
}

void Pin::setModule(Module *_module)
{
    module = _module;
}

void Pin::setDirection(int _direction)
{
    direction = _direction;
}

void Module::addPin(Pin *_pin)
{
    modulePins.push_back(_pin);
    nets.push_back(_pin->net); //! nets 是 modulePins 的派生索引，方便按单元遍历线网（如算单元 HPWL）
}

//! 左下角坐标（2D 投影，丢掉 z）
POS_2D Module::getLL_2D()
{
    POS_2D ll_2D;
    ll_2D.x = coor.x;
    ll_2D.y = coor.y;
    return ll_2D;
}

//! 右上角坐标 = 左下角 + 宽高。注意这里假设单元未旋转/翻转（orientation 不影响包围盒）
POS_2D Module::getUR_2D()
{
    POS_2D ur_2D;

    ur_2D.x = coor.x;
    ur_2D.y = coor.y;

    ur_2D.x += width;
    ur_2D.y += height;

    // assert(width != 0 && height != 0);
    return ur_2D;
}

void Module::setOrientation(int _oritentation)
{
    orientation = _oritentation;
}

//! need to check if coor is out side of the chip!!! but should be done in placeDB
//! 按「左下角」定位：先存 coor，再反推 center。私有函数，只能由 PlaceDB 调用
void Module::setLocation_2D(float _x, float _y, float _z)
{
    coor.x = _x;
    coor.y = _y;
    coor.z = _z;
    // update center
    center.x = coor.x + (float)0.5 * width; //! be careful of float problems
    center.y = coor.y + (float)0.5 * height;
    center.z = coor.z; //! z 上厚度为 0，中心 z 就等于底部 z
}

//! 按「中心」定位：先存 center，再反推左下角 coor。ePlace/优化器主要走这条路径
void Module::setCenter_2D(float _x, float _y, float _z)
{
    center.x = _x;
    center.y = _y;
    center.z = _z;
    // update coor
    coor.x = center.x - (float)0.5 * width; //! be careful of float problems
    coor.y = center.y - (float)0.5 * height;
    coor.z = center.z;
}

//! 布局行的左下角 = start
POS_2D SiteRow::getLL_2D()
{
    return start;
}

//! 布局行的右上角：end 存的是「右下角」，所以 y 要再加一个行高
POS_2D SiteRow::getUR_2D()
{
    POS_2D ur_2D = end;
    ur_2D.y += height;
    return ur_2D;
}
