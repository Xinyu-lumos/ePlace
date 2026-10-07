#include "plot.h"

// ============================================================================
// plot.cpp —— 布局可视化实现
//
// 两个绘图函数结构几乎相同（plotEPlace_2D 多了 filler 图层），流程都是：
//   1. 取 plotPath（命令行 -plotPath，缺省 "./"）
//   2. 按 chipRegion 长宽比算出图像尺寸，保证短边不小于 1000 像素
//   3. 逐个 module 换算成像素矩形画上去
//   4. 写标题文字、存 BMP
//
// 坐标换算要点：布局坐标系 y 轴向上、原点在 chipRegion.ll；
// 图像坐标系 y 轴向下、原点在左上角，所以 getY 要做一次镜像翻转。
// ============================================================================

using namespace PLOTTING;

//! 布局 x → 图像 x：平移到以 regionLLx 为原点，再按 unitX 缩放
int PLOTTING::getX(float regionLLx, float x, float unitX)
{
    return (x - regionLLx) * unitX;
}

// the Y-axis must be mirrored
//! 布局 y → 图像 y：先平移到以 regionLLy 为原点，再用 regionHeight 减去它完成上下翻转
int PLOTTING::getY(float regionHeight, float regionLLy, float y, float unitY)
{
    return (regionHeight - (y - regionLLy)) * unitY; //?
    // return (chipRegionHeight - y) * unitY;//?
}

//! 画「纯布局」图层：terminal（蓝）+ std cell（红）/ macro（橙）
void PLOTTING::plotCurrentPlacement(string imageName, PlaceDB *db)
{
    string plotPath;
    if (!gArg.GetString("plotPath", &plotPath))
    {
        plotPath = "./";
    }

    float chipRegionWidth = db->chipRegion.ur.x - db->chipRegion.ll.x;
    float chipRegionHeight = db->chipRegion.ur.y - db->chipRegion.ll.y;

    int minImgaeLength = 1000; //! 短边最小像素数，保证小芯片也能看清（原变量名 Imgae 拼写错误）

    int imageHeight;
    int imageWidth;

    float opacity = 0.7; //! 半透明，方便看清重叠区域
    int xMargin = 30, yMargin = 30; //! 四周留白，避免贴边的单元被裁掉

    //! 等比例缩放：让较长的一边正好等于 minImgaeLength，另一边按长宽比推算
    if (chipRegionWidth < chipRegionHeight)
    {
        imageHeight = 1.0 * chipRegionHeight / (chipRegionWidth / minImgaeLength);
        imageWidth = minImgaeLength;
    }
    else
    {
        imageWidth = 1.0 * chipRegionWidth / (chipRegionHeight / minImgaeLength);
        imageHeight = minImgaeLength;
    }

    CImg<unsigned char> img(imageWidth + 2 * xMargin, imageHeight + 2 * yMargin, 1, 3, 255); //! 3 通道 RGB，白底

    float unitX = imageWidth / chipRegionWidth,
          unitY = imageHeight / chipRegionHeight; //! 布局单位长度对应多少像素

    //! 先画 terminal（蓝）。NI 类 terminal 是占位性质，跳过不画
    for (Module *curTerminal : db->dbTerminals)
    {
        assert(curTerminal);
        // ignore pin's location
        if (curTerminal->isNI)
        {
            continue;
        }
        int x1 = getX(db->chipRegion.ll.x, curTerminal->getLL_2D().x, unitX) + xMargin;
        int x2 = getX(db->chipRegion.ll.x, curTerminal->getUR_2D().x, unitX) + xMargin;
        int y1 = getY(chipRegionHeight, db->chipRegion.ll.y, curTerminal->getLL_2D().y, unitY) + yMargin;
        int y2 = getY(chipRegionHeight, db->chipRegion.ll.y, curTerminal->getUR_2D().y, unitY) + yMargin;
        img.draw_rectangle(x1, y1, x2, y2, Blue, opacity);
    }

    //! 再画可动单元：macro 用橙色区分于红色 std cell
    for (Module *curNode : db->dbNodes)
    {
        assert(curNode);
        int x1 = getX(db->chipRegion.ll.x, curNode->getLL_2D().x, unitX) + xMargin;
        int x2 = getX(db->chipRegion.ll.x, curNode->getUR_2D().x, unitX) + xMargin;
        int y1 = getY(chipRegionHeight, db->chipRegion.ll.y, curNode->getLL_2D().y, unitY) + yMargin;
        int y2 = getY(chipRegionHeight, db->chipRegion.ll.y, curNode->getUR_2D().y, unitY) + yMargin;
        if (curNode->isMacro)
        {
            img.draw_rectangle(x1, y1, x2, y2, Orange, opacity);
        }
        else
        {
            img.draw_rectangle(x1, y1, x2, y2, Red, opacity);
        }
    }

    img.draw_text(50, 50, imageName.c_str(), Black, NULL, 1, 30);
    img.save_bmp(string(plotPath + imageName + string(".bmp")).c_str());
    cout << "INFO: BMP HAS BEEN SAVED: " << imageName + string(".bmp") << endl;
}

//! 画 ePlace 的布局：在 plotCurrentPlacement 的基础上，
//! 加 -debug 或处于 FILLERONLY 阶段时，额外把 filler 用绿色画出来
//! （filler 平时不画，是因为数量太多会盖住真实单元）
void PLOTTING::plotEPlace_2D(string imageName, EPlacer_2D *eplacer)
{
    string plotPath;
    if (!gArg.GetString("plotPath", &plotPath))
    {
        plotPath = "./";
    }

    float chipRegionWidth = eplacer->db->chipRegion.ur.x - eplacer->db->chipRegion.ll.x;
    float chipRegionHeight = eplacer->db->chipRegion.ur.y - eplacer->db->chipRegion.ll.y;

    int minImgaeLength = 1000;

    int imageHeight;
    int imageWidth;

    float opacity = 0.7;
    int xMargin = 30, yMargin = 30;

    if (chipRegionWidth < chipRegionHeight)
    {
        imageHeight = 1.0 * chipRegionHeight / (chipRegionWidth / minImgaeLength);
        imageWidth = minImgaeLength;
    }
    else
    {
        imageWidth = 1.0 * chipRegionWidth / (chipRegionHeight / minImgaeLength);
        imageHeight = minImgaeLength;
    }

    CImg<unsigned char> img(imageWidth + 2 * xMargin, imageHeight + 2 * yMargin, 1, 3, 255);

    float unitX = imageWidth / chipRegionWidth,
          unitY = imageHeight / chipRegionHeight;

    for (Module *curTerminal : eplacer->db->dbTerminals)
    {
        assert(curTerminal);
        // ignore pin's location
        if (curTerminal->isNI)
        {
            continue;
        }
        int x1 = getX(eplacer->db->chipRegion.ll.x, curTerminal->getLL_2D().x, unitX) + xMargin;
        int x2 = getX(eplacer->db->chipRegion.ll.x, curTerminal->getUR_2D().x, unitX) + xMargin;
        int y1 = getY(chipRegionHeight, eplacer->db->chipRegion.ll.y, curTerminal->getLL_2D().y, unitY) + yMargin;
        int y2 = getY(chipRegionHeight, eplacer->db->chipRegion.ll.y, curTerminal->getUR_2D().y, unitY) + yMargin;
        img.draw_rectangle(x1, y1, x2, y2, Blue, opacity);
    }

    for (Module *curNode : eplacer->db->dbNodes)
    {
        assert(curNode);
        int x1 = getX(eplacer->db->chipRegion.ll.x, curNode->getLL_2D().x, unitX) + xMargin;
        int x2 = getX(eplacer->db->chipRegion.ll.x, curNode->getUR_2D().x, unitX) + xMargin;
        int y1 = getY(chipRegionHeight, eplacer->db->chipRegion.ll.y, curNode->getLL_2D().y, unitY) + yMargin;
        int y2 = getY(chipRegionHeight, eplacer->db->chipRegion.ll.y, curNode->getUR_2D().y, unitY) + yMargin;
        if (curNode->isMacro)
        {
            img.draw_rectangle(x1, y1, x2, y2, Orange, opacity);
        }
        else
        {
            img.draw_rectangle(x1, y1, x2, y2, Red, opacity);
        }
    }

    if ((gArg.CheckExist("debug") || eplacer->placementStage == FILLERONLY))
    {
        //! filler 图层（绿）
        for (Module *curNode : eplacer->ePlaceFillers)
        {
            assert(curNode);
            int x1 = getX(eplacer->db->chipRegion.ll.x, curNode->getLL_2D().x, unitX) + xMargin;
            int x2 = getX(eplacer->db->chipRegion.ll.x, curNode->getUR_2D().x, unitX) + xMargin;
            int y1 = getY(chipRegionHeight, eplacer->db->chipRegion.ll.y, curNode->getLL_2D().y, unitY) + yMargin;
            int y2 = getY(chipRegionHeight, eplacer->db->chipRegion.ll.y, curNode->getUR_2D().y, unitY) + yMargin;
            img.draw_rectangle(x1, y1, x2, y2, Green, opacity);
        }
    }

    img.draw_text(50, 50, imageName.c_str(), Black, NULL, 1, 30);
    img.save_bmp(string(plotPath + imageName + string(".bmp")).c_str());
    cout << "INFO: BMP HAS BEEN SAVED: " << imageName + string(".bmp") << endl;
}
