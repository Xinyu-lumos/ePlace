// ============================================================================
// arghandler.cpp —— 命令行参数解析实现
//
// 每个 Get* 都是同一套逻辑：
//   1. 先查 m_override（运行时写入的值），命中就直接返回 —— override 优先级最高
//   2. 没命中就用 FindCaptionIndex 在 argv 里找 caption 的下标
//   3. 值就是 caption 后面那一个 token（m_argv[index + 1]）
// ============================================================================
#include "arghandler.h"
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <map>

CArgHandler gArg; // global variable //! 全工程共享的唯一实例

CArgHandler::CArgHandler()
{
    m_argc = 0;
    m_argv = NULL;
    //! 疑似问题：m_isDev 在构造函数里没初始化（只在 Init() 里赋值），
    //!   若在 Init 之前调用 IsDev() 会读到未初始化的值
    m_debugLevel = 0;
}

CArgHandler::~CArgHandler()
{
}

void CArgHandler::Init(const int argc, char *argv[])
{
    m_isDev = false;
    m_argc = argc;
    m_argv = argv;
    m_debugLevel = GetCount("v"); //! -v 的个数当作调试等级（-vvv → 3）

    // check "devdev"
    //! 开发者模式判定：直接按下标逐个比对字符 "devdev"（等价于 argv[i] == "-devdev"）
    //! 疑似问题：没有先检查 argv[i] 的长度和首字符，若参数短于 7 个字符会越界读
    for (int i = 1; i < argc; i++)
    {
        if (argv[i][1] == 'd' && argv[i][2] == 'e' && argv[i][3] == 'v' &&
            argv[i][4] == 'd' && argv[i][5] == 'e' && argv[i][6] == 'v')
        {
            m_isDev = true;
        }
    }
}

//! 取整数参数：override 优先，其次命令行
//! 疑似问题：`if (index + 1 > m_argc)` 应为 `>=`（m_argc 是元素个数，最大合法下标是 m_argc-1），
//!   当 caption 恰好是最后一个参数时会越界读一个。GetDouble/GetFloat/GetString 同样问题
bool CArgHandler::GetInt(const string caption, int *variable)
{
    map<string, string>::const_iterator ite;
    ite = m_override.find(caption);
    if (ite != m_override.end())
    {
        *variable = atoi(ite->second.c_str());
        return true;
    }

    int index = FindCaptionIndex(caption);

    if (index < 0)
        return false;

    if (index + 1 > m_argc)
        return false;

    *variable = atoi(m_argv[index + 1]);
    //return false;
    return true;
}

bool CArgHandler::GetDouble(const string caption, double *variable)
{
    map<string, string>::const_iterator ite;
    ite = m_override.find(caption);
    if (ite != m_override.end())
    {
        *variable = atof(ite->second.c_str());
        return true;
    }

    int index = FindCaptionIndex(caption);

    if (index < 0)
        return false;

    if (index + 1 > m_argc)
        return false;

    *variable = atof(m_argv[index + 1]);
    //return false;
    return true;
}

// (kaie)
//! 取 float 参数，逻辑与 GetInt / GetDouble 完全一致
bool CArgHandler::GetFloat(const string caption, float *variable)
{
    map<string, string>::const_iterator ite;
    ite = m_override.find(caption);
    if (ite != m_override.end())
    {
        *variable = (float)atof(ite->second.c_str());
        return true;
    }

    int index = FindCaptionIndex(caption);

    if (index < 0)
        return false;

    if (index + 1 > m_argc)
        return false;

    *variable = (float)atof(m_argv[index + 1]);
    //return false;
    return true;
}
// @(kaie)

// 2007-02-13 (donnie)
//! const char* 版本只是转发到 string 版本，方便调用方直接写字面量
bool CArgHandler::GetString(const char *caption, string *variable)
{
    return GetString(string(caption), variable);
}

//! 取字符串参数：override 优先，其次命令行
bool CArgHandler::GetString(const string caption, string *variable)
{
    // look for caption in the map
    map<string, string>::const_iterator ite;
    ite = m_override.find(caption);
    if (ite != m_override.end()) //! 命中 override：直接返回写入过的值
    {

        *variable = ite->second;
        return true;
    }
    // not found in map, look for caption in input arguments(argv)
    int index = FindCaptionIndex(caption);

    if (index < 0)
        return false;

    if (index + 1 > m_argc)
        return false;

    *variable = m_argv[index + 1];
    //return false;
    return true;
}

//! 开关型参数是否存在（如 -fullPlot、-noQP）。override 里只要写入过就算存在
bool CArgHandler::CheckExist(const string caption)
{
    map<string, string>::const_iterator ite;
    ite = m_override.find(caption);
    if (ite != m_override.end())
        return true;

    if (FindCaptionIndex(caption) > 0)
        return true;
    return false;
}

//! 统计某个参数在命令行里出现了几次（用于 -v 这类可重复的参数）
int CArgHandler::GetCount(const string caption)
{
    int count = 0;
    //! 从 i=1 开始跳过程序名；m_argv[i]+1 是跳过前导的 '-'
    for (int i = 1; i < m_argc; i++)
    {
        if (strcmp(m_argv[i] + 1, caption.c_str()) == 0)
            count++;
    }
    return count;
}

//! 在 argv 里查找 caption 的下标（跳过 '-' 后做全串比较），找不到返回 -1
int CArgHandler::FindCaptionIndex(const string caption)
{
    for (int i = 1; i < m_argc; i++)
    {
        //! 疑似问题：若 argv[i] 不是以 '-' 开头（比如用户误传了裸参数），
        //!   +1 会把首字符当参数名的一部分来比，可能误命中
        if (strcmp(m_argv[i] + 1, caption.c_str()) == 0)
            return i;
    }
    return -1;
}

//! 取消某个参数的 override，使其恢复读取命令行原值；原本没有 override 则返回 false
bool CArgHandler::RemoveOverride(const string caption)
{
    map<string, string>::iterator ite;
    ite = m_override.find(caption);
    if (ite == m_override.end())
        return false;
    m_override.erase(ite);
    return true;
}

//! 运行时写入/覆盖一个参数值。本工程里 main 靠它把 benchmark 名、输出/绘图路径
//! 写回 gArg，后续模块就能直接 GetString 取用，不必层层传递
void CArgHandler::Override(const string caption, const string value)
{
    m_override[caption] = value;// add caption to argument map and set value
}
