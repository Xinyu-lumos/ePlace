#ifndef EPLACE_H
#define EPLACE_H

#include "global.h"
#include "objects.h"
#include "placedb.h"
#include "fft.h"
// #include "plot.h"

//! 线长变化的参考尺度，用于把 ΔHPWL 归一化后算 λ 的乘子（见 updatePenaltyFactor）
#define DELTA_HPWL_REF 350000
#define PENALTY_MULTIPLIER_BASE 1.05 //1.05, equals PENALTY_MULTIPLIER_UPPERBOUND, follow Xplace(param_scheduler.py, step_density_weight()) and RePlAce
//! λ 每轮调整的上下界：单轮最多放大 / 缩小 5%，保证收敛平稳
#define PENALTY_MULTIPLIER_UPPERBOUND 1.05
#define PENALTY_MULTIPLIER_LOWERBOUND 0.95

//! 布局三阶段：
//!   mGP        宏观布局，所有可动单元（含 macro）+ filler 一起优化
//!   FILLERONLY 只重撒并优化 filler，用来打破 mGP 后的局部拥塞
//!   cGP        单元布局，macro 冻结，只优化 std cell + filler
enum PLACEMENT_STAGE
{
    mGP,
    FILLERONLY,
    cGP

};

class Bin_2D;
class Bin_3D; // call it cube?
class EPlacer_2D;

//! 一个 bin（bin 网格里的格子）：记录几何信息 + 四类密度 + 电场/电势
class Bin_2D
{
public:
    Bin_2D()
    {
        init();
    }
    POS_2D center;
    POS_2D ll;
    POS_2D ur;
    float width;
    float height;
    float area;

    //! here Density actually means Area(electric quantity) rather than electric density, the true density should be calculated as: area/bin area,
    //! and here the densities are actually the overlap area between each component and bin
    float nodeDensity;     //! density due to movable modules, use this to calculate overflow! see ePlace paper equation(37)
    float fillerDensity;   //! density due to fillers, use this and node density to do nesterov optimization
    float terminalDensity; // terminalDensity and baseDensity are calculated in binInitialization
    float baseDensity;     // see bin->virt_area in RePlAce, baseDensity: area inside a bin but not inside a placement row
    // todo: add virtual area?

    VECTOR_2D E; // electric field, see eplace paper
    float phi;   // electric potential, see eplace paper
    void init()
    {
        center.SetZero();
        ll.SetZero();
        ur.SetZero();
        nodeDensity = 0;
        fillerDensity = 0;
        terminalDensity = 0;
        baseDensity = 0;
        E.SetZero();
        phi = 0;
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

// ============================================================================
// EPlacer_2D —— ePlace 全局布局器（2D）
//
// 用法（由 Optimization/nesterov.hpp 的 EplaceNesterovOpt 驱动）：
//   初始化：setTargetDensity → initialization
//   每轮迭代：getPosition → 优化器算步长 → setPosition → updatePenaltyFactor
//            → totalGradientUpdate → getGradient → …
//   阶段切换：switch2FillerOnly / switch2cGP
//
// 下标约定（非常重要）：ePlaceNodesAndFillers / ePlaceCellsAndFillers / ePlaceFillers
// 三个容器里元素的下标，必须与 totalGradient / cGPGradient / fillerGradient 一一对应，
// 代码里到处都有 assert(index == module->idx) 在守这个约定。
// ============================================================================
class EPlacer_2D
{
public:
    EPlacer_2D()
    {
        init();
    }

    EPlacer_2D(PlaceDB *_db)
    {
        init();
        db = _db;
    }
    vector<vector<Bin_2D *>> bins; //! bins[x][y]，第一维是 X 方向；元素为裸指针，需外部释放

    float ePlaceStdCellArea; // calculated in fillerInitialization
    float ePlaceMacroArea;

    VECTOR_2D_INT binDimension; // How many bins in X/Y direction
    VECTOR_2D binStep;          // length of a bin in X/Y direction

    PlaceDB *db; //! 非拥有指针，生命周期由调用方管理

    vector<Module *> ePlaceFillers;         // for FILLERONLY placement stores all filler cells. I think filler cells shouldn't be stored in placedb. is this unnecessary when we have ePlaceNodesAndFillers?
    vector<Module *> ePlaceNodesAndFillers; //! contains nodes and fillers, for mGP
    vector<Module *> ePlaceCellsAndFillers; //! contains std cells and fillers, for cGP

    vector<VECTOR_3D> wirelengthGradient; // store wirelength gradient for nodes only(wirelength gradient for filler is always 0)
    vector<VECTOR_3D> densityGradient;    // store density gradient for fillers nodes (wirelength gradient for filler node is always 0)
    vector<VECTOR_3D> totalGradient;      // total gradient of objective function f including gradients of all components, for mGP
    vector<VECTOR_3D> fillerGradient;
    vector<VECTOR_3D> cGPGradient; // cell and fillers

    float targetDensity;         //!!!!!! so important //! 目标密度（如 0.9），密度缩放和溢出判据都以此为准
    float globalDensityOverflow; // !!!!! so important, tau //! 全局溢出率 τ，Nesterov 的收敛判据
    VECTOR_2D invertedGamma;     // gamma of the wa wirelength model,here we actually use 1/gamma, following RePlAce. gamma is different for different dimension

    float lambda; // penalty factor //! 惩罚因子，平衡线长与密度
    double lastHPWL; //! 上一次的线长，用于算 ΔHPWL 来调整 λ

    int placementStage;     //! 当前阶段，取值见 PLACEMENT_STAGE
    int mGPIterationCount;  //! mGP 阶段跑了多少轮，switch2cGP 回调 λ 时用到

    void init()
    {
        bins.clear();
        binDimension.SetZero();
        binStep.SetZero();

        targetDensity = 0;
        globalDensityOverflow = 0;
        invertedGamma.SetZero();
        lambda = 0.0;
        lastHPWL = 0.0;

        ePlaceStdCellArea = 0;
        ePlaceMacroArea = 0;

        db = NULL;
        ePlaceFillers.clear();
        ePlaceNodesAndFillers.clear();

        wirelengthGradient.clear();
        densityGradient.clear();

        placementStage = mGP;
        mGPIterationCount=0;
        
    }
    void setTargetDensity(float);
    void setPlacementStage(int);

    //! 一次性初始化：filler 插入 → bin 网格 → 梯度容器 → 首轮梯度 → λ 初值
    void initialization();

    void fillerInitialization(); //! 计算空白面积并按论文式(13)插入 filler
    void binInitialization();    //! 划分 bin 网格，并算好 terminalDensity / baseDensity（只算一次）
    void gradientVectorInitialization();

    void binNodeDensityUpdate();  //! only consider density from movable modules(nodes) in this function, because terminal density only needed to be calculated once, in binInitializaton()
    void densityOverflowUpdate(); // called in wirelenghGradientUpdate() //! 更新全局溢出率 τ（不计 filler）
    void wirelengthGradientUpdate(); //! 顺带更新 τ 与 γ，再求线长梯度
    void densityGradientUpdate();    //! FFT 解泊松方程求电场，再折算成密度梯度

    void totalGradientUpdate(); //! 合成 ∇f = λ∇D − ∇WL，并按阶段填充各梯度容器

    //! 以下三个是优化器的接口，返回值随 placementStage 变化
    vector<VECTOR_3D> getGradient();
    vector<VECTOR_3D> getPosition();

    void setPosition(vector<VECTOR_3D>);

    void penaltyFactorInitilization(); // initialize lambda 0, see ePlace paper equation 35 //! λ₀ = Σ|∇WL| / Σ|∇D|
    void updatePenaltyFactor();        //! 按 ΔHPWL 自适应缩放 λ，见论文式 36
    //! be aware of density scaling and local smooth in density calculation

    void switch2FillerOnly(); //! mGP → FILLERONLY，重撒 filler
    void switch2cGP();        //! mGP → cGP，冻结 macro 并回调 λ

    void showInfo();
    void showInfoFinal(); //! 注：只有声明、全工程找不到定义，调用会链接失败
    // void plotCurrentPlacement(string); //! 绘图已迁到 Plot/plot.cpp

private:
    vector<VECTOR_3D> getModulePositions(vector<Module *>);
};
#endif