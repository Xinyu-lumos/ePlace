#ifndef OBJECTS_H
#define OBJECTS_H
// ============================================================================
// objects.h —— 布局数据库的「原子对象」定义
//
// 层次关系： Net（线网） ⟷ Pin（引脚） ⟷ Module（单元/terminal/filler）
//           SiteRow（布局行）→ Interval（行内可用区间）
//           Tier（一个 2D 平面，3D IC 时有多层）
//
// 关键约定：
//   · Pin 同时挂在所属 Module（modulePins）和所属 Net（netPins）上，是双向索引的枢纽
//   · Pin 存的是相对所属 Module 的 offset，绝对坐标要加上 Module 的位置（见 getAbsolutePos）
//   · Module 的 coor / center 是私有的，只能通过 db->setModuleCenter/Location 修改，
//     以保证坐标变动时同步刷新其上所有 Pin 的绝对坐标
// ============================================================================
#include "global.h"
const int PIN_DIRECTION_OUT = 0;
const int PIN_DIRECTION_IN = 1;
const int PIN_DIRECTION_UNDEFINED = -1;
// todo: use maps for indexing by name
class Module;
class SiteRow;
class Row;
class Pin;
class Net;
class Tier;
class CRect;
class PlaceDB;
class Interval;

//! 线网：一组需要连在一起的 pin。线长模型（HPWL / WA / LSE）都以它为计算单位
class Net
{
public:
    Net()
    {
        init();
    }
    Net(int index)
    {
        init();
        idx = index;
    }
    int idx;
    vector<Pin *> netPins;
    //! boundPin pointers, used for quadratic placement utilizing bound2bound net model
    //! 边界 pin：X/Y/Z 方向上的最左/最右（最小/最大）pin。
    //! 只要求出这 6 个就能算 HPWL，二次布局也只在这些 pin 对之间建边（见 QPlace）
    Pin *boundPinXmin;
    Pin *boundPinXmax;
    Pin *boundPinYmin;
    Pin *boundPinYmax;
    Pin *boundPinZmin; // for 3D
    Pin *boundPinZmax; // for 3D

    //! WA（weighted-average）线长模型的中间量，预先累加好，求梯度时直接复用
    VECTOR_3D numeratorMax_WA;
    VECTOR_3D denominatorMax_WA;
    VECTOR_3D numeratorMin_WA;
    VECTOR_3D denominatorMin_WA;

    //! LSE（log-sum-exp）线长模型的中间量，同上
    VECTOR_3D sumMax_LSE;
    VECTOR_3D sumMin_LSE;

    void init()
    {
        idx = 0;
        netPins.clear();
        clearBoundPins();
        //! 疑似问题：这里把 numeratorMax_WA 和 denominatorMin_WA 各清零了两次，
        //!   而 numeratorMin_WA 和 denominatorMax_WA 从未被初始化（遗留的复制粘贴错误）
        numeratorMax_WA.SetZero();
        denominatorMin_WA.SetZero();
        numeratorMax_WA.SetZero();
        denominatorMin_WA.SetZero();
    }
    void addPin(Pin *);
    int getPinCount();
    void allocateMemoryForPin(int);
    double calcNetHPWL();   //! 直接遍历所有 pin 求 HPWL（精确但慢）
    double calcBoundPin();  //! 用边界 pin 求 HPWL，并顺带刷新 boundPin*（快，主流程用这个）
    void clearBoundPins();
    //! WA / LSE 两种平滑线长模型：可微，用来替代不可微的 HPWL 供梯度下降使用
    double calcWirelengthWA_2D(VECTOR_2D);
    double calcWirelengthLSE_2D(VECTOR_2D);
    //! 求某个 pin 上的线长梯度（对该 pin 所属 module 的坐标求导）
    VECTOR_2D getWirelengthGradientWA_2D(VECTOR_2D, Pin *);
    VECTOR_2D getWirelengthGradientLSE_2D(VECTOR_2D, Pin *);
};

//! 引脚：连接 Module 与 Net 的纽带。位置以「相对所属 module 的 offset」存储
class Pin
{
public:
    Pin()
    {
        init();
    }
    Pin(Module *masterModule, Net *masterNet, float x, float y)
        : direction(PIN_DIRECTION_UNDEFINED)
    {
        init();
        offset = POS_2D(x, y);
        setModule(masterModule);
        setNet(masterNet);
    }
    void init()
    {
        idx = -1;
        module = NULL;
        net = NULL;
        offset.SetZero();
        absolutePos.SetZero();
        direction = -1;
        eMin_WA.SetZero();
        eMax_WA.SetZero();
        expZeroFlgMax_WA.SetZero();
        expZeroFlgMin_WA.SetZero();
    }
    int idx;
    Module *module; //! 所属单元（反向指针）
    Net *net;       //! 所属线网（反向指针）
    POS_2D offset;  //! 相对 module 左下角的偏移，module 移动时它不变

    int direction; // 0 output  1 input  -1 not-define
    POS_3D getAbsolutePos(); //! 绝对坐标 = module 位置 + offset，用于线长/密度计算
    // POS_3D fetchAbsolutePos(); // currently for mLG only
    void calculateAbsolutePos(); //! 重算并缓存 absolutePos，由 module 移动后统一触发

    //! WA 模型里 e^[(Xmin−Xi)/γ] 一类的中间量；expZeroFlg* 标记指数是否已小到可视为 0（防下溢）
    VECTOR_3D eMin_WA;               // e^[(Xmin-Xi)/gamma] in WA model (X/Y/Z)
    VECTOR_3D eMax_WA;               // e^[(Xi-Xmax)/gamma] in WA model (X/Y/Z)
    VECTOR_3D_BOOL expZeroFlgMax_WA; // indicate if (Xi-Xmax)/gamma is too small that e^[(Xi-Xmax)/gamma] == 0
    VECTOR_3D_BOOL expZeroFlgMin_WA; // similar

    //! LSE 模型对应的中间量（当前代码路径只启用其中一种模型）
    VECTOR_3D eMin_LSE;               // e^[(Xmin-Xi)/gamma] in LSE model
    VECTOR_3D eMax_LSE;               // e^[(Xi-Xmax)/gamma] in LSE model
    VECTOR_3D_BOOL expZeroFlgMax_LSE; // indicate if (Xi-Xmax)/gamma is too small that e^[(Xi-Xmax)/gamma] == 0
    VECTOR_3D_BOOL expZeroFlgMin_LSE; // similar

    void setId(int);
    void setNet(Net *);
    void setModule(Module *);
    void setDirection(int);

// private:
    POS_3D absolutePos; //! 缓存的绝对坐标，避免每次重复计算
};

//! 单元：std cell / 可移动 macro / 固定 terminal / 虚拟 filler 通通都是 Module
class Module
{
public:
    friend class PlaceDB; //! PlaceDB 是唯一可以直接改 coor / center 的地方
    Module()
    {
        Init();
    }

    Module(int _index, string _name, float _width = 0, float _height = 0, bool _isFixed = false, bool _isNI = false)
    {
        Init();
        name = _name;
        width = _width;
        height = _height;
        area = float_mul(width, height); //! area calculated here!
        isFixed = _isFixed;
        isNI = _isNI; // 2022-05-13 (frank)
        idx = _index;
        assert(area >= 0);
    }
    int idx;
    Tier *tier; //! 所属层（3D IC 用；2D 布局里通常只有一个 tier）
    string name;
    float width;
    float height;
    float area;
    float orientation; //! 朝向（N/R90/S/R180…），合法化时用来判断能否翻转
    bool isMacro;      //! 硬宏（通常尺寸大、不参与 cGP、密度要缩放）
    bool isFixed;      //! 固定不动（terminal / 预放置单元）
    bool isNI;         //! "NI" terminal：占位性质，不参与密度与线长统计
    bool isFiller;     //! ePlace 的虚拟填充单元，只贡献密度、不连线
    vector<Pin *> modulePins; //! 该单元上的所有引脚
    vector<Net *> nets;       //! 该单元参与的所有线网（由 modulePins 派生，用于加速）
    void Init()
    {
        idx = -1;
        coor.SetZero();
        center.SetZero();
        width = 0;
        height = 0;
        area = 0;
        orientation = 0;
        isMacro = false;
        isFiller = false;
        isFixed = false;
        isNI = false;
        tier = NULL;
    }
    float calcArea()
    {
        area = float_mul(width, height);
        return area;
    }
    void addPin(Pin *);
    //! 下面这一组只读 accessor；注意 getCenter/getLocation 返回的是缓存值，
    //! 只有经过 db->setModuleCenter/Location 才会更新
    string getName() { return name; }
    float getWidth() { return width; }
    float getHeight() { return height; }
    POS_3D getLocation() { return coor; }
    POS_3D getCenter() { return center; }
    POS_2D getLL_2D();
    POS_2D getUR_2D();
    float getArea() { return area; }
    short int getOrientation() { return orientation; } //! 注意返回类型是 short int，而成员是 float，存在隐式截断
    void setOrientation(int);

private:
    //! these 2 functions should only be called in db->setModuleCenter/Location!!!
    //! 之所以设为私有：坐标一旦变动，必须同步刷新该 module 上所有 pin 的绝对坐标，
    //! 把入口收敛到 PlaceDB 才能保证不遗漏
    POS_3D coor;   // coor for coordinate //! 左下角坐标
    POS_3D center; // coor of center, be aware that center should be recalculated every time the module is moved, or before HPWL calculation
    void setLocation_2D(float, float, float = 0); //! 按左下角定位
    void setCenter_2D(float, float, float = 0);   //! 按中心定位
};

//! 抽象行（目前似乎只作为基类/占位存在）
class Row
{ // an abstract row
    //! 疑似问题：class 默认访问权限是 private，这两个构造函数没有 public 修饰，
    //!   外部根本无法构造 Row；且 Row 在全工程中未见被使用
    Row()
    {
        bottom = 0;
        height = 0;
        step = 0;
        start.SetZero();
        end.SetZero();
    }

    Row(double _bottom, double _height, double _step) : bottom(_bottom),
                                                        height(_height),
                                                        step(_step)
    {
        start.SetZero();
        end.SetZero();
    }
    double bottom;
    double height;
    double step;
    POS_2D start;
    POS_2D end;
};

//! 布局行（site row）：std cell 必须整行对齐地摆在这些行上
class SiteRow // a place row
{
public:
    SiteRow()
    {
        intervals.clear();
        bottom = 0;
        height = 0;
        step = 0;
        start.SetZero();
        end.SetZero();
        orientation = OR_N;
    }

    SiteRow(double _bottom, double _height, double _step) : bottom(_bottom),
                                                            height(_height),
                                                            step(_step),
                                                            orientation(OR_N)
    {
        intervals.clear();
        start.SetZero();
        end.SetZero();
    }

    double bottom;              // The bottom y coordinate of this SiteRow of sites //! 行的下边界 y
    double height;              // The height of this SiteRow of sites             //! 行高（= 单元高度）
    double step;                // The minimum x step of SiteRow.	by indark        //! x 方向最小步进（= site 宽度，合法化时用它做吸附）
    POS_2D start;               // lower left coor of this row;                    //! 行的左下角
    POS_2D end;                 //! lower right coor of this row;                  //! 行的右下角（注释写 lower right，实际是右上）
    ORIENT orientation;         // donnie 2006-04-23  N (0) or S (1)               //! 行朝向，影响单元能否翻转放置
    vector<Interval> intervals; //! 行内被 macro / 死区打断后剩余的可放置区间
    // double site_spacing;// site spacing in bookshelf format, equals to site width
    POS_2D getLL_2D();
    POS_2D getUR_2D();
    //! 下面两个是给 sort 用的比较器（按 bottom 排序），写成成员函数而非静态有点别扭但能用
    bool Lesser(SiteRow &r1, SiteRow &r2)
    {
        return (r1.bottom < r2.bottom);
    }
    bool Greater(SiteRow &r1, SiteRow &r2)
    {
        return (r1.bottom > r2.bottom);
    }
    bool isInside(const double &x, const double width)
    {
        //! 疑似问题：函数体被注释掉了，现在无条件 return false，等于没实现。
        //!   原本应遍历 intervals 判断 [x, x+width] 是否落在某个可用区间内
        // vector<double>::const_iterator ite;
        // for (ite = interval.begin(); ite != interval.end(); ite += 2)
        // {
        //     if (*ite > x)
        //         return false;
        //     // cout << "  sites(" << *ite << " " << *(ite+1) << ") ";
        //     if (*ite <= x && *(ite + 1) >= x)
        //         return true;
        // }
        return false;
    }
    friend inline std::ostream &operator<<(std::ostream &os, const SiteRow &row)
    {
        os << "SiteRow start: " << row.start << " SiteRow end:" << row.end;
        return os;
    }
};

//! 一个 2D 平面（3D IC 时每层一个 Tier；普通 2D 布局只有一个）
class Tier // one 2D plane(or one 2D chip)
{
public:
    Tier()
    {
        coreRegion.Init();
        layerNumber = -1;
        modules.clear();
        siteRows.clear();
    }
    CRect coreRegion;          //! 该层的可放置区域
    int layerNumber;           // which layer in a 3dic
    double rowHeight;          //! 该层的标准行高
    vector<Module *> modules;  //! 该层上的单元
    vector<SiteRow *> siteRows; //! 该层的布局行
    vector<Module *> terminals; // terminals: module that can't be moved //! 固定不动的 terminal
};

//! 布局行上一段「可放置」的 x 区间
class Interval
{
    // horizontal intervals of placement siterows, only store x coordinate, instead of POS_2D
    // there are intervals in rows because of macros or pre-defined dead zones.
    // interval is not an empty! empty is the available space in a row that is not covered by std cells
    // and interval is calculated without considering std cell locations
public:
    Interval()
    {
        SetZero();
    }
    Interval(float _start, float _end)
    {
        start = _start;
        end = _end;
    }
    inline void SetZero()
    {
        start = end = 0.0; //!! 0.0!!!!
    }
    float getLength()
    {
        return end - start;
    }
    float start;
    float end;
};
#endif