// ============================================================================
// 文件/模块总览：Bookshelf 格式解析器的对外接口（实现见 parser.cpp）
// ----------------------------------------------------------------------------
// 职责：
//   把 EDA 全局布局（global placement）算法所用的 Bookshelf 格式 benchmark
//   读入 PlaceDB。本文件只负责声明接口。
//
// 支持的文件格式（UCLA / Bookshelf 系列 ASCII 文件）：
//   .aux   设计入口文件，一行列出该设计由哪些子文件组成，例如：
//          RowBasedPlacement : d.nodes d.nets d.wts d.pl d.scl
//   .nodes 单元/module 清单：名字、宽、高；第 4 列为 terminal / terminal_NI
//          时表示固定的 IO pad，否则是可移动 cell
//   .nets  线网清单：NetDegree 头部 + 每个 pin 所属的 module 名与 pin offset
//   .pl    初始位置文件：module 名、左下角坐标、方向
//   .scl   布局 site row 定义：行的底边坐标、行高、site 宽/间距、subrow 等
//
// 解析出的核心数据结构（定义在 PlaceDB/placedb.h 与 PlaceCommon/objects.h）：
//   dbNodes     -> vector<Module*>   可移动 std cell / macro
//   dbTerminals -> vector<Module*>   固定 IO pad
//   dbPins      -> vector<Pin*>      net 与 module 的连接关系（含 pin offset）
//   dbNets      -> vector<Net*>      线网
//   dbSiteRows  -> vector<SiteRow>   布局行，含 Interval 描述的可用 x 区间
//   moduleMap   -> map<string, Module*>  名字到 Module 指针的映射
//   commonRowHeight / coreRegion / chipRegion 由以上数据进一步推导
//
// 主要流程：
//   ReadFile(.aux) 解析入口文件 -> ReadSCLFile -> ReadNodesFile ->
//   ReadNetsFile -> ReadPLFile(init=true) -> setChipRegion_2D()
//   ReadPLFile 也可单独调用（init=false），用于把布局结果重新读回 PlaceDB。
// ============================================================================
#ifndef PARSER_H
#define PARSER_H
#include "placedb.h"
#include "global.h"
// BookshelfParser：Bookshelf 格式（.aux/.nodes/.nets/.pl/.scl）解析器。
// 调用方只需使用 ReadFile；ReadPLFile 单独暴露，便于重复载入 .pl。
class BookshelfParser
{
public:
    // 解析入口：读取 .aux 文件，按其列出的子文件名依次解析整个设计到 PlaceDB
    int ReadFile(string file, PlaceDB &db);
    // 读取 .pl 文件，设置每个 module 的初始坐标与方向。
    // init=true ：随 .aux 一起首次读取，此时文件名需要拼上 gArg("path") 前缀；
    // init=false：file 本身就是完整路径（例如读入上一轮布局输出的结果文件）。
    int ReadPLFile(string file, PlaceDB &db, bool init);

    // 以下三个子解析器属于内部实现。调用顺序不能随意调整：
    // .scl 要先读（确定行高），.nodes 次之（建立 moduleMap），
    // 最后才是 .nets（pin 需要通过名字反查 Module）。
private:
    // 解析 .scl：填充 db.dbSiteRows，并确定 db.commonRowHeight 与 coreRegion
    int ReadSCLFile(string file, PlaceDB &db);
    // 解析 .nodes：填充 db.dbNodes / db.dbTerminals 与 db.moduleMap
    int ReadNodesFile(string file, PlaceDB &db);
    // 解析 .nets：填充 db.dbNets / db.dbPins，并建立 net-pin-module 三方联系
    int ReadNetsFile(string file, PlaceDB &db);
};

// LEFDEFParser：为将来支持 LEF/DEF 工业格式预留的空壳类，目前没有任何实现。
class LEFDEFParser
{
    // maybe in the future...
};

#endif