// ============================================================================
// 模块总览：ePlaceAbacus —— 「只读入全局布局结果 → 做 Abacus 合法化 → 输出」的独立可执行入口
// ============================================================================
// 【这个文件的职责】
//   CMakeLists.txt 中名为 ePlaceAbacus 的可执行程序入口。
//   它**不跑** QPlace / EPlace 全局布局，而是假定磁盘上已经存在一份
//   全局布局的结果（<benchmarkName>-global.pl），只负责把这份重叠的结果合法化，
//   用于单独调试 Legalization 模块（对比之下，完整流程由 EPlace/eplace.cpp 驱动）。
//
// 【算法思路 / 主流程】
//   ① 解析 Bookshelf 输入（.aux 及其引用的 .nodes/.nets/.scl/.pl），建立 PlaceDB；
//      顺带从路径中推出 benchmark 名与绘图目录（先 rm -rf 再 mkdir，保证图是本次的）。
//   ② 读入上一轮全局布局的 .pl 坐标，覆盖掉 .nodes 里的初始位置；
//   ③ 画一张合法化前的图（global_RESULT）留作对照；
//   ④ AbacusLegalizer::legalization()：
//        initializeSubrows() 内部会先调 PlaceDB::removeBlockedSite()，
//        把被 macro / terminal 挡住的 site 挖掉，得到每条 site row 可用的 intervals，
//        再展开成 subrow；随后逐个 std cell（按 x 升序）做
//        「双向搜索候选行 + |Δy| 剪枝 + 行内 Abacus 聚类」，
//        消除重叠并把 cell 对齐到 site 网格与 site row，同时最小化位移。
//   ⑤ 画一张合法化后的图（LEGALIZED_RESULT），并把结果写回 Bookshelf 文件。
//
// 【核心数据结构】
//   BookshelfParser / PlaceDB（dbNodes、dbTerminals、dbSiteRows、coreRegion）
//   AbacusLegalizer（内部是 subrows → AbacusRow → AbacusCellCluster）
//
// 【主要入口函数】
//   main（本文件）；真正的算法入口见 legalizer.cpp 的 AbacusLegalizer::legalization()
// ============================================================================

#include "parser.h"
#include "arghandler.h"
#include "qplace.h"
#include "eplace.h"
#include "legalizer.h"
#include "plot.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
    BookshelfParser parser;
    PlaceDB *placedb = new PlaceDB();
    //! 疑似问题：placedb 用 new 分配但从未 delete（进程退出时由 OS 回收，属于可接受的工具级泄漏）。
    gArg.Init(argc, argv);

    if (argc < 2)
    {
        return 0;
    }
    // argv[1]+1 跳过开头的 '-'，于是 "-aux" 变成 "aux" 再比较
    if (strcmp(argv[1] + 1, "aux") == 0) // -aux, argv[1]=='-'
    {
        // bookshelf
        printf("Use BOOKSHELF placement format\n");

        string filename = argv[2];
        //! 疑似问题：没有检查 argc >= 3，若用户只给了 -aux 而没给文件名，argv[2] 会越界。
        string::size_type pos = filename.rfind("/");
        string benchmarkName;
        if (pos != string::npos)
        {
            printf("    Path = %s\n", filename.substr(0, pos + 1).c_str());
            gArg.Override("path", filename.substr(0, pos + 1));

            int length = filename.length();

            benchmarkName = filename.substr(pos + 1, length - pos);

            int len = benchmarkName.length();
            // 去掉 .aux 后缀，得到纯 benchmark 名，后续所有输出文件都以它为前缀
            if (benchmarkName.substr(len - 4, 4) == ".aux")
            {
                benchmarkName = benchmarkName.erase(len - 4, 4);
            }
            gArg.Override("benchmarkName", benchmarkName);
            cout << "    Benchmark: " << benchmarkName << endl;

            string plotPath;
            if (!gArg.GetString("plotPath", &plotPath))
            {
                plotPath = "./";
            }

            plotPath += "/" + benchmarkName + "/";

            // 清空上一次的绘图目录，避免新旧图片混在一起
            string cmd = "rm -rf " + plotPath;
            system(cmd.c_str());
            cmd = "mkdir -p " + plotPath;
            system(cmd.c_str());

            gArg.Override("plotPath", plotPath);

            cout << "    Plot path: " << plotPath << endl;
        }

        parser.ReadFile(argv[2], *placedb);
    }
    placedb->showDBInfo();

    string outputPath;
    string benchmarkName;
    //! 疑似问题：这里在读回 -global.pl **之前**就先 outputBookShelf() 了一次，
    //! 会把尚未合法化（甚至还是初始随机位置）的坐标写出，覆盖已有的输出文件。
    //! 若只是为了生成目录/文件骨架，代价是可能破坏上一轮的结果。
    placedb->outputBookShelf();
    gArg.GetString("outputPath", &outputPath);
    gArg.GetString("benchmarkName", &benchmarkName);

    // 读入上一阶段（EPlace 全局布局）产出的坐标，覆盖掉 .nodes 中的初始位置
    parser.ReadPLFile(outputPath + "/" + benchmarkName + "-global.pl", *placedb, false);
    PLOTTING::plotCurrentPlacement("global_RESULT",placedb);
    //! 疑似问题：这里**没有**先跑 SAMacroLegalizer 做 macro 合法化。
    //! 若设计中含 macro 且它们之间仍有重叠，removeBlockedSite() 切出的可用区间会偏大，
    //! std cell 就可能被放到 macro 底下的非法区域。完整流程里 mLG 应当在 Abacus 之前执行。
    AbacusLegalizer legalizer(placedb);
    legalizer.legalization();
    PLOTTING::plotCurrentPlacement("LEGALIZED_RESULT",placedb);
    placedb->outputBookShelf();
}
