// ============================================================================
// ePlace_main.cpp —— 完整布局流程的主入口
//
// 流水线（每个阶段都可用命令行参数跳过）：
//   1. 解析 Bookshelf 设计                     (-aux xxx.aux)
//   2. 二次布局 QPlace，产出重叠的初始解        (-noQP 跳过)
//   3. ePlace 全局布局 mGP                     (-nomGP 跳过)
//   4. macro 合法化 mLG（仅当设计里有 macro）   (-nomLG 跳过)
//   5. filler 重撒 FILLERONLY → 单元布局 cGP   (-nocGP 跳过)
//   6. 输出 Bookshelf 结果文件
//   7. 合法化：外部 ntuplace3 或内部 Abacus     (-noLegal 跳过)
//   8. 详细布局（局部重排/交换）                (-internalDP 启用)
//
// 其它常用参数：
//   -targetDensity 0.9   目标密度    -loadpl xxx.pl  读入已有布局（须配合 -noQP）
//   -addNoise            给初始解加扰动  -fullPlot     每轮出图
//   -bb / -bktrk         切换 Nesterov 步长策略
//   -legalizerPath <dir> 外部 ntuplace3 所在目录
//   -internalLegal <x>   用内部 AbacusLegalizer（目前只支持纯 std cell 设计）
// ============================================================================
#include "qplace.h"
#include "parser.h"
#include "arghandler.h"
#include "eplace.h"
#include "opt.hpp"
#include "nesterov.hpp"
#include "legalizer.h"
#include "detailed.h"
#include "plot.h"
#include <iostream>
#include <cstdio>
#include <unistd.h>
using namespace std;

int main(int argc, char *argv[])
{
    // Capture both C and C++ logging in the working directory.
    //! 把 stdout 重定向到 DUMP.txt，再把 stderr 也指向同一个文件，
    //! 这样 C 的 printf 和 C++ 的 cout 日志会合并到单个文件里，便于事后分析
    if (freopen("DUMP.txt", "w", stdout) == nullptr)
    {
        perror("Cannot open DUMP.txt");
        return 1;
    }
    if (dup2(fileno(stdout), fileno(stderr)) == -1)
    {
        perror("Cannot redirect stderr to DUMP.txt");
        return 1;
    }

    BookshelfParser parser;
    PlaceDB *placedb = new PlaceDB();
    gArg.Init(argc, argv); //! 命令行参数解析进全局 gArg，后续所有模块都从这里读配置

    if (argc < 2)
    {
        return 0;
    }
    if (strcmp(argv[1] + 1, "aux") == 0) // -aux, argv[1]=='-'
    {
        // bookshelf
        printf("Use BOOKSHELF placement format\n");

        //! 从 -aux 的路径里拆出「目录」和「benchmark 名」，分别写回 gArg 供后续输出使用
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

            string outputFilePath;
            if (!gArg.GetString("outputPath", &outputFilePath))
            {
                outputFilePath = "./";
            }
            outputFilePath = outputFilePath + "/" + benchmarkName + "/";
            gArg.Override("outputPath", outputFilePath);

            string plotPath = outputFilePath + "Graphs/";
            gArg.Override("plotPath", plotPath);

            //! 清空并重建输出目录（rm -rf + mkdir -p），保证每次运行的结果目录是干净的
            //! 注意：直接把路径拼进 system() 命令，路径含空格或特殊字符会有问题
            string cmd = "rm -rf " + outputFilePath;
            system(cmd.c_str());

            cmd = "mkdir -p " + outputFilePath;
            system(cmd.c_str());

            cmd = "rm -rf " + plotPath;
            system(cmd.c_str());

            cmd = "mkdir -p " + plotPath;
            system(cmd.c_str());

            cout << "    Plot path: " << plotPath << endl;
        }
        parser.ReadFile(argv[2], *placedb);
    }
    placedb->showDBInfo(); //! 打印设计的规模信息（单元数、线网数、行数、core 区域等）
    string plPath;
    if (gArg.GetString("loadpl", &plPath))
    {
        //! modules will be moved to center in QP, so if QP is not skipped, loading module locations from an existing pl file is meaningless
        parser.ReadPLFile(plPath, *placedb, false);
    }

    //! ---- 阶段 2：二次布局，产出初始解 ----
    QPPlacer *qpplacer = new QPPlacer(placedb);
    if (!gArg.CheckExist("noQP"))
    {
        qpplacer->quadraticPlacement();
    }

    //! 可选：给初始解加随机扰动（幅度为一个 bin 步长），常用于做鲁棒性实验
    if (gArg.CheckExist("addNoise"))
    {
        placedb->addNoise(); // the noise range is [-avgbinStep,avgbinStep]
    }

    double mGPTime;
    double mLGTime;
    double FILLERONLYtime;
    double cGPTime;

    //! ---- 阶段 3：ePlace 全局布局 ----
    EPlacer_2D *eplacer = new EPlacer_2D(placedb);

    //! 目标密度：越小越宽松但浪费面积，越大越紧凑但难合法化；缺省 1.0
    float targetDensity;
    if (!gArg.GetFloat("targetDensity", &targetDensity))
    {
        targetDensity = 1.0;
    }

    eplacer->setTargetDensity(targetDensity);
    eplacer->initialization(); //! filler 插入 + bin 网格 + 首轮梯度 + λ 初值

    //! 优化器用基类指针持有，主循环（init → 迭代 → wrap_up）在 opt.hpp 里
    FirstOrderOptimizer<VECTOR_3D> *opt = new EplaceNesterovOpt<VECTOR_3D>(eplacer);

    if (!gArg.CheckExist("nomGP"))
    {
        cout << "mGP started!\n";

        time_start(&mGPTime);
        opt->opt(); //! 跑到 stop_condition 成立（τ 达标或超 MAX_ITERATION）才返回
        time_end(&mGPTime);

        cout << "mGP finished!\n";
        cout << "Final HPWL: " << int(placedb->calcHPWL()) << endl;
        cout << "mGP time: " << mGPTime << endl;
        PLOTTING::plotCurrentPlacement("mGP result", placedb);
    }

    ///////////////////////////////////////////////////
    // legalization and detailed placement
    ///////////////////////////////////////////////////

    //! ---- 阶段 4~5：macro 合法化，然后重撒 filler 再跑 cGP ----
    //! 注意：只有设计里确实有 macro 才走这条路
    if (placedb->dbMacroCount > 0 && !gArg.CheckExist("nomLG"))
    {
        //! mLG：模拟退火式的 macro 合法化，把重叠的 macro 拉开并对齐
        SAMacroLegalizer *macroLegalizer = new SAMacroLegalizer(placedb);
        macroLegalizer->setTargetDensity(targetDensity);
        cout << "Start mLG, total macro count: " << placedb->dbMacroCount << endl;
        time_start(&mLGTime);
        macroLegalizer->legalization();
        time_end(&mLGTime);

        PLOTTING::plotCurrentPlacement("mLG result", placedb);
        cout << "mLG finished. HPWL after mLG: " << int(placedb->calcHPWL()) << endl;
        cout << "mLG time: " << mLGTime << endl;
        // exit(0);

        if (!gArg.CheckExist("nocGP"))
        {
            //! FILLERONLY：macro 已固定，重撒 filler 打破 mGP 残留的局部拥塞
            eplacer->switch2FillerOnly();
            cout << "filler placement started!\n";

            time_start(&FILLERONLYtime);
            opt->opt();
            time_end(&FILLERONLYtime);

            cout << "filler placement finished!\n";
            //! 疑似问题：这里打印的是 mGPTime，应改为 FILLERONLYtime
            cout << "FILLERONLY time: " << mGPTime << endl;
            PLOTTING::plotCurrentPlacement("FILLERONLY result", placedb);

            //! cGP：macro 冻结，只优化 std cell + filler
            eplacer->switch2cGP();
            cout << "cGP started!\n";

            time_start(&cGPTime);
            opt->opt();
            time_end(&cGPTime);

            cout << "cGP finished!\n";
            cout << "cGP Final HPWL: " << int(placedb->calcHPWL()) << endl;
            //! 疑似问题：这里打印的是 mGPTime，应改为 cGPTime
            cout << "cGP time: " << mGPTime << endl;
            PLOTTING::plotCurrentPlacement("cGP result", placedb);
        }
    }

    //! ---- 阶段 6：输出全局布局结果（供外部合法器使用）----
    placedb->outputBookShelf("eGP",false); // output, files will be used for legalizers such as ntuplace3

    //! ---- 阶段 7：合法化（消除重叠、对齐 site row）----
    if (!gArg.CheckExist("noLegal"))
    {
        string legalizerPath;
        if (gArg.GetString("legalizerPath", &legalizerPath))
        {
            //! 路线 A：调用外部的 ntuplace3 可执行程序做合法化 + 详细布局
            string outputAUXPath;
            string outputPLPath;
            string outputPath;
            string benchmarkName;

            gArg.GetString("outputAUX", &outputAUXPath);
            gArg.GetString("outputPL", &outputPLPath);
            gArg.GetString("outputPath", &outputPath);
            gArg.GetString("benchmarkName", &benchmarkName);

            //! -noglobal 表示只做合法化/详细布局，不再做全局布局
            string cmd = legalizerPath + "/ntuplace3" + " -aux " + outputAUXPath + " -loadpl " + outputPLPath + " -noglobal" + " -out " + outputPath + benchmarkName + "-ntu" + " > " + outputPath + "ntuplace3-log.txt";

            cout << RED << "Running legalizer and detailed placer: " << cmd << RESET << endl;
            system(cmd.c_str());
        }
        else if (gArg.GetString("internalLegal", &legalizerPath))
        {
            //! 路线 B：用内置的 AbacusLegalizer（目前只支持纯 std cell 设计，没有 macro 合法器）
            // for std cell design only, since we don't have macro legalizer for now
            cout << "Calling internal legalizer: " << endl;
            cout << "Global HPWL: " << int(placedb->calcHPWL()) << endl;

            AbacusLegalizer *legalizer = new AbacusLegalizer(placedb);
            legalizer->legalization();

            cout << "Legal HPWL: " << int(placedb->calcHPWL()) << endl;
            PLOTTING::plotCurrentPlacement("Cell legalized result", placedb);

            placedb->outputBookShelf("eLG",true);
        }
        else
        {
            cout<<"Legalization not done!!!\n";
            exit(0);
        }
    }

    //! ---- 阶段 8：内置详细布局（局部重排 / 交换，进一步压线长）----
    if (gArg.CheckExist("internalDP"))// currently only works after internal legalization
    {
        cout << "Calling internal detailed placement: " << endl;
        cout << "HPWL before detailed placement: " << int(placedb->calcHPWL()) << endl;

        DetailedPlacer *detailedPlacer = new DetailedPlacer(placedb);
        detailedPlacer->detailedPlacement();

        cout << "HPWL after detailed placement: " << int(placedb->calcHPWL()) << endl;
        PLOTTING::plotCurrentPlacement("Detailed placement result", placedb);

        placedb->outputBookShelf("eDP",true);
    }
}
