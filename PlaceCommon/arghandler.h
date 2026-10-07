#ifndef ARGHANDLER_H
#define ARGHANDLER_H

// ============================================================================
// arghandler.h —— 命令行参数解析器
//
// 支持的写法：以 '-' 开头的 caption，后面跟 0 个或多个值，例如
//     ePlace -aux design.aux -targetDensity 0.9 -fullPlot -noQP
//   · -aux design.aux     带一个值 → GetString("aux", &s)
//   · -fullPlot / -noQP   不带值   → CheckExist("fullPlot")（开关型参数）
//
// 除了解析，还支持 Override()：在运行时**写入/覆盖**某个参数的值。
// 这个功能在本工程里很关键——main 里解析出 benchmark 名和输出目录后，
// 会 Override 回 gArg，供后面所有模块直接读取，不用层层传参。
//
// 全局单例 gArg 在 arghandler.cpp 里定义，全工程共享。
// ============================================================================
#include <string>
#include <cstring>
#include <map>
using namespace std;
/**
	@author Indark <indark@eda.ee.ntu.edu.tw>
	2006-09-22  Updated by donnie
*/

class CArgHandler
{
public:
	CArgHandler();
	~CArgHandler();

	//! 保存 main 的 argc/argv，并解析 -v 数量作为 debugLevel
	void Init(const int argc, char *argv[]);

	//! 取带值的参数，返回 false 表示该参数不存在（此时不改动 *variable）
	bool GetInt(const string caption, int *variable);
	bool GetDouble(const string caption, double *variable);
	bool GetFloat(const string caption, float *variable); // kaie
	bool GetString(const string caption, string *value);
	bool GetString(const char *caption, string *value);
	int GetCount(const string caption); //! 统计某参数出现了几次（如 -v -v -v → 3，用作 debug 等级）

	//! 运行时写入/覆盖参数值（覆盖后 Get* 优先读 override 表）
	void Override(const string caption, const string value);
	bool RemoveOverride(const string caption); //! 取消覆盖，恢复读命令行原值

	bool CheckExist(string caption); //! 开关型参数是否存在
	bool IsDev() { return m_isDev; };          //! 是否为开发者模式（连续两个 -devdev）
	int GetDebugLevel() { return m_debugLevel; }; //! -v 的个数

private:
	char **m_argv;
	int m_argc;
	int m_debugLevel;
	bool m_isDev;
	map<string, string> m_override; //! Override() 写入的值存在这里，优先于命令行

private:
	int FindCaptionIndex(const string caption); //! 在 argv 里找 caption 的下标，找不到返回 -1
};

extern CArgHandler gArg; // global variable //! 全工程共享的全局单例

#endif
