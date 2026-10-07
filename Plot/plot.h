#ifndef PLOT_H
#define PLOT_H
// ============================================================================
// plot.h —— 布局可视化：用 CImg 把当前布局画成 BMP 图片
//
// 配色约定：terminal 蓝、macro 橙、std cell 红、ePlace filler 绿
// 用途：调试布局过程。加 -fullPlot 会逐轮出图，方便观察单元是如何被「推开」的。
// 依赖 CImg（以 cimg_display=0 编译，即纯文件输出、不需要 X11）。
// ============================================================================
#include "global.h"
#include "CImg.h"
#include "placedb.h"
#include "eplace.h"

using namespace cimg_library;
namespace PLOTTING
{
    const unsigned char Blue[] = {120, 200, 255},
                        Black[] = {0, 0, 0},
                        Green[] = {0, 150, 0},
                        Orange[] = {255, 165, 0},
                        Red[] = {255, 0, 0};
    //! 物理坐标 → 图像像素坐标的换算
    int getX(float, float, float);
    int getY(float, float, float, float); //! Y 轴需要翻转（图像 y 向下，布局 y 向上）

    void plotCurrentPlacement(string, PlaceDB *);   //! 画纯布局（terminal + 单元）
    void plotEPlace_2D(string, EPlacer_2D *);       //! 额外把 ePlace 的 filler 也画出来
}

#endif