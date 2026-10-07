// qplace 的独立调试入口：只读入 bookshelf 设计、跑一遍二次布局、画图看结果。
// 命令行：qplace_test -aux <path>/<benchmark>.aux [-IPiteCount N] [-fullPlot]
#include "parser.h"
#include "arghandler.h"
#include "qplace.h"
// #include "plot.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
    PlaceDB *placedb = new PlaceDB();
    gArg.Init(argc, argv); //! 把命令行参数解析进全局 gArg，后续 QPPlacer 里的 IPiteCount / debug / fullPlot 都从这里读

    if (argc < 2)
    {
        return 0;
    }
    if (strcmp(argv[1] + 1, "aux") == 0) // -aux, argv[1]=='-'
    {
        // bookshelf
        printf("Use BOOKSHELF placement format\n");

        string filename = argv[2];
        string::size_type pos = filename.rfind("/");
        string benchmarkName;
        if (pos != string::npos)
        {
            printf("    Path = %s\n", filename.substr(0, pos + 1).c_str());
            gArg.Override("path", filename.substr(0, pos + 1));

            int length = filename.length();

            benchmarkName = filename.substr(pos + 1, length - pos);

            int len = benchmarkName.length();
            if (benchmarkName.substr(len - 4, 4) == ".aux")
            {
                benchmarkName = benchmarkName.erase(len - 4, 4);
            }
            gArg.Override("benchmarkName", benchmarkName);
            cout << "    Benchmark: " << benchmarkName << endl;

            string plotPath;
            gArg.GetString("plotPath", &plotPath);

            plotPath += "/" + benchmarkName + "/";

            string cmd = "mkdir -p " + plotPath;
            system(cmd.c_str());

            gArg.Override("plotPath", plotPath);

            cout << "    Plot path: " << plotPath << endl;
        }
        BookshelfParser parser;
        parser.ReadFile(argv[2], *placedb);
    }
    placedb->showDBInfo(); //! 打印单元/线网/terminal 数量、core 区域等基本信息
    QPPlacer *qpplacer = new QPPlacer(placedb);
    qpplacer->quadraticPlacement(); //! 跑二次布局；加 -fullPlot 可逐轮输出布局图
}