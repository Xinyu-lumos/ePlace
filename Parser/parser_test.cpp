// ============================================================================
// 文件/模块总览：Parser 模块的独立测试 / 冒烟入口
// ----------------------------------------------------------------------------
// 职责：
//   构造一个 PlaceDB，用 BookshelfParser 读入命令行指定的 .aux 设计，
//   用来单独验证 Parser 能否正确解析，不牵涉后续的任何布局算法。
//
// 用法：
//   parser_test -aux <design>.aux
//
// 支持格式与核心数据结构：见 parser.h / parser.cpp 开头的总览说明
// （Bookshelf：.aux -> .nodes / .nets / .wts / .pl / .scl -> PlaceDB）。
//
// 主要流程：
//   1) 检查命令行参数是否为 -aux；
//   2) 从 .aux 的路径里截取出目录，用 gArg.Override("path", ...) 登记，
//      这样后面各个子解析器拼接子文件名时才能找到它们；
//   3) 调用 BookshelfParser::ReadFile 完成整份设计的解析。
// ============================================================================
#include "parser.h"
#include "arghandler.h"
#include <iostream>
using namespace std;

    // 程序入口：仅用于冒烟测试 Parser 本身
int main(int argc, char *argv[])
{
    // 所有解析结果都写入这个 PlaceDB 实例
    PlaceDB placedb;

    // 至少要给出 -aux 选项。
    //! 疑似问题：argc < 2 的判据不够——真正被访问的是 argv[2]，
    //! 应当判断 argc < 3。目前传入 "-aux" 而不带文件名时，argv[2] 越界读取。
    //! 另外参数不合法时在不开任何提示的情况下直接返回 0，不易排查。
    if (argc < 2)
    {
        return 0;
    }
    // argv[1] 形如 "-aux"，argv[1] + 1 跳过前导的短横杠，
    // 再与 "aux" 比较，等价于同时容忍 "-aux" 这种写法；
    // 只有匹配时才会继续解析（其它参数直接被忽略）。
    if (strcmp(argv[1] + 1, "aux") == 0) // -aux, argv[1]=='-'
    {
        // bookshelf
        printf("Use BOOKSHELF placement format\n");

        // 取出 .aux 所在目录并登记为全局变量 path，
        // 供 ReadSCLFile / ReadNodesFile / ReadNetsFile / ReadPLFile 拼路径使用。
        // rfind 找的是最后一个斜杠；若命令行给的是裸文件名（不含斜杠），
        // 则 path 保持默认，子文件会从当前工作目录下查找。
        string filename = argv[2];
        string::size_type pos = filename.rfind("/");
        if (pos != string::npos)
        {
            printf("    Path = %s\n", filename.substr(0, pos + 1).c_str());
            gArg.Override("path", filename.substr(0, pos + 1));
        }

        // 创建解析器并把整个设计读入 placedb
        BookshelfParser parser;
        parser.ReadFile(argv[2], placedb);
    }
}