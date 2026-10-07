// ============================================================================
// ePlace_test.cpp —— 完整布局流程的测试入口
//
// 与 main/ePlace_main.cpp 的流水线**完全相同**，唯一区别是它不把日志重定向到 DUMP.txt，
// 而是直接打到终端，方便调试时实时观察。
//
// 流水线：解析 → QPlace 二次布局 → ePlace mGP → mLG(macro 合法化)
//        → FILLERONLY → cGP → 输出 → 合法化 → 详细布局
// 各阶段可用 -noQP / -nomGP / -nomLG / -nocGP / -noLegal 跳过，
// -internalLegal 走内置 Abacus 合法器，-internalDP 启用内置详细布局。
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
using namespace std;

int main(int argc, char *argv[])
{
    BookshelfParser parser;
    PlaceDB *placedb = new PlaceDB();
    gArg.Init(argc, argv); //! 命令行参数解析进全局 gArg

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

            string outputFilePath;
            if (!gArg.GetString("outputPath", &outputFilePath))
            {
                outputFilePath = "./";
            }
            outputFilePath = outputFilePath + "/" + benchmarkName + "/";
            gArg.Override("outputPath", outputFilePath);

            string plotPath = outputFilePath + "Graphs/";
            gArg.Override("plotPath", plotPath);

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
    placedb->showDBInfo();
    string plPath;
    if (gArg.GetString("loadpl", &plPath))
    {
        //! modules will be moved to center in QP, so if QP is not skipped, loading module locations from an existing pl file is meaningless
        parser.ReadPLFile(plPath, *placedb, false);
    }

    //! ---- 阶段 2：二次布局产出初始解（加 -noQP 可跳过）----
    QPPlacer *qpplacer = new QPPlacer(placedb);
    if (!gArg.CheckExist("noQP"))
    {
        qpplacer->quadraticPlacement();
    }

    //! 可选：给初始解加一个 bin 步长的随机扰动
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

    //! 目标密度，缺省 1.0
    float targetDensity;
    if (!gArg.GetFloat("targetDensity", &targetDensity))
    {
        targetDensity = 1.0;
    }

    eplacer->setTargetDensity(targetDensity);
    eplacer->initialization(); //! filler 插入 + bin 网格 + 首轮梯度 + λ 初值

    //! 优化器用基类指针持有，主循环在 opt.hpp 的 opt() 里
    FirstOrderOptimizer<VECTOR_3D> *opt = new EplaceNesterovOpt<VECTOR_3D>(eplacer);

    if (!gArg.CheckExist("nomGP"))
    {
        cout << "mGP started!\n";

        time_start(&mGPTime);
        opt->opt(); //! 跑到 τ 达标或超过 MAX_ITERATION 才返回
        time_end(&mGPTime);

        cout << "mGP finished!\n";
        cout << "Final HPWL: " << int(placedb->calcHPWL()) << endl;
        cout << "mGP time: " << mGPTime << endl;
        PLOTTING::plotCurrentPlacement("mGP result", placedb);
    }

    ///////////////////////////////////////////////////
    // legalization and detailed placement
    ///////////////////////////////////////////////////

    //! ---- 阶段 4~5：macro 合法化 → 重撒 filler → cGP（仅当设计含 macro）----
    if (placedb->dbMacroCount > 0 && !gArg.CheckExist("nomLG"))
    {
        SAMacroLegalizer *macroLegalizer = new SAMacroLegalizer(placedb); //! mLG：模拟退火式 macro 合法化
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
            //! FILLERONLY：macro 已固定，重撒 filler 打破残留拥塞
            eplacer->switch2FillerOnly();
            cout << "filler placement started!\n";

            time_start(&FILLERONLYtime);
            opt->opt();
            time_end(&FILLERONLYtime);

            cout << "filler placement finished!\n";
            //! 疑似问题：打印的是 mGPTime，应改为 FILLERONLYtime
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
            //! 疑似问题：打印的是 mGPTime，应改为 cGPTime
            cout << "cGP time: " << mGPTime << endl;
            PLOTTING::plotCurrentPlacement("cGP result", placedb);
        }
    }

    //! ---- 阶段 6：输出全局布局结果 ----
    placedb->outputBookShelf("eGP",false); // output, files will be used for legalizers such as ntuplace3

    //! ---- 阶段 7：合法化：外部 ntuplace3（-legalizerPath）或内部 Abacus（-internalLegal）----
    if (!gArg.CheckExist("noLegal"))
    {
        string legalizerPath;
        if (gArg.GetString("legalizerPath", &legalizerPath))
        {
            string outputAUXPath;
            string outputPLPath;
            string outputPath;
            string benchmarkName;

            gArg.GetString("outputAUX", &outputAUXPath);
            gArg.GetString("outputPL", &outputPLPath);
            gArg.GetString("outputPath", &outputPath);
            gArg.GetString("benchmarkName", &benchmarkName);

            string cmd = legalizerPath + "/ntuplace3" + " -aux " + outputAUXPath + " -loadpl " + outputPLPath + " -noglobal" + " -out " + outputPath + benchmarkName + "-ntu" + " > " + outputPath + "ntuplace3-log.txt";

            cout << RED << "Running legalizer and detailed placer: " << cmd << RESET << endl;
            system(cmd.c_str());
        }
        else if (gArg.GetString("internalLegal", &legalizerPath))
        {
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

    //! ---- 阶段 8：内置详细布局 ----
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