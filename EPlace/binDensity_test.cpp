// ePlace 密度计算的性能/正确性测试入口：
//   先跑二次布局拿到初始解，再跑一遍 ePlace 的初始化和单轮梯度计算，打印耗时。
// 命令行：binDensity_test -aux <path>/<benchmark>.aux
#include "parser.h"
#include "arghandler.h"
#include "qplace.h"
#include "eplace.h"
// #include "plot.h"
#include <iostream>
using namespace std;

int main(int argc, char *argv[])
{
    PlaceDB *placedb = new PlaceDB();
    gArg.Init(argc, argv);

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

            string cmd = "rm -rf "+plotPath+";mkdir -p " + plotPath;
            system(cmd.c_str());

            gArg.Override("plotPath", plotPath);

            cout << "    Plot path: " << plotPath << endl;
        }
        BookshelfParser parser;
        parser.ReadFile(argv[2], *placedb);
    }
    placedb->showDBInfo();
    //! ePlace 需要一个初始解，所以先跑 QPlace 的二次布局
    QPPlacer *qpplacer = new QPPlacer(placedb);
    qpplacer->quadraticPlacement();

    double initializationTime;
    double iterationTime;
    EPlacer_2D *eplacer = new EPlacer_2D(placedb);
    eplacer->setTargetDensity(0.9); //! 目标密度，越高越紧凑但越难合法化

    time_start(&initializationTime);
    eplacer->initialization(); //! filler 插入 + bin 网格 + 首轮梯度 + λ 初值
    time_end(&initializationTime);

    //! 手动跑一轮梯度（正式流程里这些由 totalGradientUpdate 内部串起来调用，这里拆开是为了分别计时）
    time_start(&iterationTime);
    eplacer->binNodeDensityUpdate();
    eplacer->wirelengthGradientUpdate();
    eplacer->densityGradientUpdate();
    //! ! 注意：这里传了参数 1.0，但 eplace.h 里 totalGradientUpdate() 是无参的，
    //!   该测试文件已与当前接口脱节，无法编译（历史遗留）
    eplacer->totalGradientUpdate(1.0);
    time_end(&iterationTime);

    cout << "Initialization time: " << initializationTime << endl;
    cout << "Iteration time: " << iterationTime << endl;
}