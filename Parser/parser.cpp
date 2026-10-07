// ============================================================================
// 文件/模块总览：Bookshelf 格式解析器实现
// ----------------------------------------------------------------------------
// 职责：
//   实现 BookshelfParser 的四个解析函数，把 .aux 指向的 benchmark 文件
//   （.scl / .nodes / .nets / .pl）逐步读入 PlaceDB。
//
// 支持的文件格式与典型内容（可参考 examples/smoke/ 下的样例）：
//   .aux   RowBasedPlacement : smoke.nodes smoke.nets smoke.wts smoke.pl smoke.scl
//   .scl   UCLA scl 1.0 / NumRows : 20 / CoreRow Horizontal / Coordinate : 0 /
//          Height : 1 / Sitewidth : 1 / Sitespacing : 1 / Siteorient : N /
//          SubrowOrigin : 0 Numsites : 20 / End
//   .nodes UCLA nodes 1.0 / NumNodes : 104 / NumTerminals : 4 / c0 1 1 /
//          p0 1 1 terminal
//   .nets  UCLA nets 1.0 / NumNets : 184 / NumPins : 368 / NetDegree : 2 n0 /
//          c0 B : 0 0   （module 名、方向、冒号后的 pin offset）
//   .pl    UCLA pl 1.0 / c0 0 0 : N
//
// 解析出的核心数据结构：
//   PlaceDB::dbNodes     可移动 Module（node）
//   PlaceDB::dbTerminals 固定 Module（terminal / terminal_NI）
//   PlaceDB::dbPins      Pin，记录所属 Module、所属 Net 及相对偏移
//   PlaceDB::dbNets      Net，记录其所有 Pin
//   PlaceDB::dbSiteRows  SiteRow，含 <bottom,height,step,start,end,intervals>
//   PlaceDB::moduleMap   module 名 -> Module*，供 nets/pl 按名字反查
//
// 主要流程：
//   ReadFile(.aux)  -> ReadSCLFile   先建立 site row 才能知道 row height
//                   -> ReadNodesFile 读入 module/terminal 并建立名字索引
//                   -> ReadNetsFile  读入 net 与 pin（依赖上面的索引）
//                   -> ReadPLFile    读入初始位置
//                   -> setChipRegion_2D() 生成用于可视化的 chip region
// ============================================================================
#include "parser.h"
#include "arghandler.h"

// ============================================================================
// 函数：BookshelfParser::ReadFile
// 职责：解析总入口。读取 .aux 文件，得到该设计由哪些 .nodes/.nets/.wts/
//       .pl/.scl 子文件组成，再依次调用各子解析器把设计灌入 PlaceDB。
// 参数：file 为 .aux 的完整路径（注意不拼 gArg("path")，见下方说明）
// 返回：恒为 0；打开文件失败时直接 exit(-1)。
// ============================================================================
int BookshelfParser::ReadFile(string file, PlaceDB &db)
{
	cout << "Reading AUX file: " << file << endl;

	// ifstream 变量名叫 out（实为输入流）。这里是唯一一处直接用传入路径打开文件的地方，
	// 其余子解析器都会先拼上 gArg("path") 前缀再由自己去打开。
	ifstream out(file);
	if (!out)
	{
		cout << "\tCannot open aux file\n";
		exit(-1);
	}

	// .aux 的有效内容其实只有第一行，形如：
	//   RowBasedPlacement : smoke.nodes smoke.nets smoke.wts smoke.pl smoke.scl
	char tmp[10000];
	out.getline(tmp, 10000);

	char tmp1[500];
	char tmp2[500];
	char file_nodes[500];
	char file_nets[500];
	char file_wts[500];
	char file_pl[500];
	char file_scl[500];

	// 按空白切出 7 个 token：tmp1 = "RowBasedPlacement"，tmp2 = "："，
	// 其后依次是 nodes / nets / wts / pl / scl 五个子文件名。
	//! 疑似问题：没有检查 sscanf 的返回值是否为 7。若某 benchmark 的 .aux 首行 token
	//! 数不足（例如缺少 .wts 或 .scl），未被填充的文件名缓冲是未初始化的随机内容，
	//! 后续 open 会以一个错误的名字失败。建议用返回值判断并给出明确报错。
	sscanf(tmp, "%s %s %s %s %s %s %s",
		   tmp1, tmp2, file_nodes, file_nets, file_wts, file_pl, file_scl);

	// 打印解析到的子文件名。注意 file_wts 只在这里被打印，本解析器从未真正读取 .wts
	// （线网权重文件），所有 net 的权重在解析阶段一律视为 1。
	printf("\t\t%s %s %s %s %s\n", file_nodes, file_nets, file_wts, file_pl, file_scl);

	// 调用顺序不能乱：
	//  1) scl   -> 得到 dbSiteRows 与 commonRowHeight，后续很多步骤依赖行高
	//  2) nodes -> 得到 dbNodes/dbTerminals，并建立 moduleMap（名字 -> Module*）
	//  3) nets  -> pin 需要靠 moduleMap 反查所属的 Module
	//  4) pl    -> 读初始坐标，setModuleLocation_2D 会按 coreRegion 做越界裁剪
	// Read SCL first to know the row height.
	ReadSCLFile(file_scl, db);	   // read core-row information
	ReadNodesFile(file_nodes, db); // blocks & cell width/height
	ReadNetsFile(file_nets, db);   // read net file
	ReadPLFile(file_pl, db, true); // initial module locations
	// coreRegion 只覆盖 site row；再把它和所有 terminal 的外框合并成 chipRegion，
	// chipRegion 一般只用于画图 / 可视化（见 PlaceDB::setChipRegion_2D）。
	db.setChipRegion_2D();
	return 0;
}

// ============================================================================
// 函数：BookshelfParser::ReadSCLFile
// 职责：解析 .scl，为每个 CoreRow 生成一条 SiteRow 写入 db.dbSiteRows，
//       并在结尾导出 db.commonRowHeight 与 coreRegion。
//
// 本函数关心每个 CoreRow 的：底边 y（bottom）、行高（height）、
// site 宽度（同时也是 x 方向最小步进 step）、行朝向 N/S，
// 以及每个 subrow 的可用区间 [起点 x, 终点 x)（存成 Interval）。
// 其余关键字（NumRows、Sitespacing、Sitesymmetry、End 等）被忽略或读作临时量。
//
// 字符串处理的边界情况：
//   - 每行只切前 6 个 token；"<值>" 所在的列位置随关键字不同而不同，
//     例如 "Coordinate : 0" 的值落在第 3 个 token；
//   - "SubrowOrigin : 0 Numsites : 20" 一行里塞了两个 key:value，
//     所以 Numsites 这个名字出现在第 4 个 token（tmp4），值在第 6 个 token。
// ============================================================================
int BookshelfParser::ReadSCLFile(string file, PlaceDB &db)
{
	// In this version, we don't care the pre-placed location.
	// Read Numrows, Height, and Numsites only.

	// 除 .aux 之外，子文件名都要拼上 gArg("path")（.aux 所在目录）才是完整路径
	string path;
	gArg.GetString("path", &path);
	path += file;
	ifstream in(path.c_str());
	if (!in)
	{
		cerr << "Cannot open .scl file: " << file << endl;
		exit(-1);
	}

	// nRows / height / lineNumber 这三个量目前读出来后并未真正参与任何逻辑：
	// nRows 没有和实际读到的 CoreRow 条数做校验，height 只被反复赋值，
	// lineNumber 只自增从未被读取（原本打算用于出错定位）。
	int nRows, height;
	nRows = height = -1;

	int lineNumber = 0;

	// 首行一般是版本头 "UCLA scl 1.0"，这里先把它消费掉；下面 #if 0 里
	// 原本有严格校验，被关掉后不带版本头的文件也能被接受。
	// check file format string
	char tmp[10000];
	in.getline(tmp, 10000);
	lineNumber++;

#if 0
	if( strcmp( "UCLA scl 1.0", tmp ) != 0 )
	{
		cerr << "SCL file header error (UCLA scl 1.0)\n";
		exit(-1);
	}
#endif

	// 每行最多取前 6 个 token 就够用
	char tmp1[500], tmp2[500], tmp3[500], tmp4[500], tmp5[500], tmp6[500];

	// vSites 直接引用 db 里的容器，push_back / back() 的写库操作即生效
	vector<SiteRow> &vSites = db.dbSiteRows;
	// Bookshelf 格式下 site 的默认宽度取 1（沿用作者的默认值）
	double SiteWidth = 1; // default value is 1 (2005/2/14 donnie)
	// 行内 site 的默认朝向为 N。orient 定义在循环外，会跨 CoreRow 保留上一次的取值
	ORIENT orient = OR_N; // 2006-04-23 (donnie)
	// 逐行做 "关键字 -> 填当前行字段" 的解析，状态就是 vSites.back()。
	//! 疑似问题：所有分支都直接用 vSites.back()，但没有先判断容器非空。
	//! 若在第一个 CoreRow 之前就出现 Coordinate / Height / Numsites 等行，
	//! 对空 vector 调用 back() 属于未定义行为。
	while (in.getline(tmp, 10000))
	{
		// 关键一步：sscanf 遇到短行不会把 6 个 %s 全部写入，未被写入的缓冲会残留
		// 上一行的内容；这里先把每个 token 串置为空串，之后才能用非空与否来判断
		// "该 token 是否存在"。否则像 "Sitespacing : 1" 这样的短行会误命中
		// 依赖靠后 token 的分支。
		tmp1[0] = tmp2[0] = tmp3[0] = tmp4[0] = tmp5[0] = tmp6[0] = '\0';
		sscanf(tmp, "%s %s %s %s %s %s", tmp1, tmp2, tmp3, tmp4, tmp5, tmp6);

		// Numrows / NumRows 两种大小写拼写都做兼容，值在第 3 个 token
		if (strcmp(tmp1, "Numrows") == 0 || strcmp(tmp1, "NumRows") == 0)
		{
			nRows = atoi(tmp3);
			////test code
			// printf("get numrows %f\n", atoi( tmp3 ) );
			////@test code
		}
		// 遇到 "CoreRow Horizontal" 表示新的一行开始：追加一个空 SiteRow
		// 作为后续几行 Coordinate/Height/... 的写入目标。
		//! 疑似问题：只接受 Horizontal 行，"CoreRow Vertical" 会被静默忽略，
		//! 既不报错也没有任何提示。
		else if (strcmp(tmp1, "CoreRow") == 0 && strcmp(tmp2, "Horizontal") == 0) // start of a row
		{
			vSites.push_back(SiteRow());
			////test code
			// printf("get corerow, site size: %d\n", vSites.size() );
			////@test code
		}
		// "Coordinate : y" 给出这一行的底边 y 坐标
		else if (strcmp(tmp1, "Coordinate") == 0)
		{
			vSites.back().bottom = atof(tmp3);
			////test code
			// printf("get coordinate %f\n", atof( tmp3 ) );
			////@test code
		}
		// "Height : h" 给出这一行的行高（即 site 的高度）
		else if (strcmp(tmp1, "Height") == 0)
		{
			vSites.back().height = atof(tmp3);
			height = atof(tmp3);
			////test code
			// printf("get height %f\n", atof( tmp3 ) );
			////@test code
		}
		// Bookshelf 里 site 宽度与 site 间距相等，因此 row.step 直接用 Sitewidth 填充。
		//! 疑似问题：文件里的 "Sitespacing" 行被完全忽略。若某 benchmark 的
		//! Sitespacing 与 Sitewidth 不相等，这里得到的 step 就是错的。
		else if (strcmp(tmp1, "Sitewidth") == 0) //! site spacing equals to site width in bookshelf format
		{
			SiteWidth = atof(tmp3);
			vSites.back().step = atof(tmp3); // step equals sitewidth in bookshelf format
											 ////test code
											 // printf("get sitewidth %f\n", atof( tmp3 ) );
											 ////@test code
		}
		// "SubrowOrigin : x Numsites : n" 一行同时含两个 key:value，Numsites 落在 tmp4，
		// 起点取 tmp3（subrow origin），终点 = 起点 + 站点数 * Sitewidth。
		// 一个 CoreRow 可以有多条 subrow（被 macro/terminal 打断），每次追加一个 Interval。
		//! 疑似问题：start / end 每次都被无条件覆盖，多 subrow 时只留下最后一条
		//! subrow 的范围，而不是所有 subrow 的最小左端与最大右端。
		else if (strcmp(tmp4, "Numsites") == 0 || strcmp(tmp4, "NumSites") == 0)
		{
			double subOrigin = atof(tmp3); //! subrowOrigin
			double numSites = atof(tmp6);
			// initialize interval
			vSites.back().intervals.push_back(Interval(subOrigin, (numSites * SiteWidth) + subOrigin));

			vSites.back().start = POS_2D(subOrigin, vSites.back().bottom);
			vSites.back().end = POS_2D((numSites * SiteWidth) + subOrigin, vSites.back().bottom);
			// printf("get numsites: %f %f\n", atof(tmp3), ( atof( tmp6 )*SiteWidth ) + atof( tmp3 ) );
		}
		// Bookshelf 只在 N / S 两种行朝向里选一：S 与 FS 都归到 OR_S，其余一律 OR_N
		// （也就是不会用到 W / E）。
		else if (strcmp(tmp1, "Siteorient") == 0 || strcmp(tmp1, "SiteOrient") == 0) // donnie 2006-04-23
		{
			if (strcmp(tmp3, "S") == 0 || strcmp(tmp3, "FS") == 0)
				orient = OR_S;
			else
				orient = OR_N;
			vSites.back().orientation = orient;
		}
	}

	// cout << "     Numrows: " << nRows << "\n";
	// cout << "      Height: " << height << "\n";
	// cout << "    Numsites: " << nSites << "\n";
	// cout << " Core region: (" << left << "," << bottom << ")-("
	//	 << left+nSites << "," << bottom + nRows * height << ")\n";

	// 假定整个设计的行高统一，取最后一行的高度作为 db.commonRowHeight；
	// 再由 setCoreRegion() 用所有行的包围盒算出 coreRegion 并累加 totalRowArea。
	//! 疑似问题：直接取 vSites.back()，既没有校验各行高度是否一致（不一致时会
	//! 默默用最后一行的值覆盖所有行），也没有防范 .scl 里一条 row 都没有的情况
	//! （空容器上调用 back() 是未定义行为）。
	db.commonRowHeight = vSites.back().height; //!
	db.setCoreRegion();

	// 固定 module（terminal / macro）所占用的 site 需要被剔除，
	// 这部分逻辑已移交给 PlaceDB::removeBlockedSite()，本函数不再处理。
	//! 疑似问题：紧接着的原注释里 "Romove" 疑为 "Remove" 的拼写错误。
	// Romove the sites occupied by the fixed module

	// included in CPlaceDB
	//@Romove the sites occupied by the fixed module

	return 0;
}

// ============================================================================
// 函数：BookshelfParser::ReadNodesFile
// 职责：解析 .nodes，把所有 module 分成两类写入 PlaceDB：
//        可移动 node      -> db.dbNodes
//        固定 IO pad      -> db.dbTerminals（terminal / terminal_NI）
//       同时建立 db.moduleMap（名字 -> Module*），供 .nets 与 .pl 按名字反查。
//
// 格式：头部给出 NumNodes / NumTerminals，正文每行一个 module：
//        c0 1 1            名字 宽 高（没有第 4 列 => 可移动 node）
//        p0 1 1 terminal   第 4 列为 terminal => 固定 IO pad
// 由 nNodes = NumNodes - NumTerminals 得到可移动单元数。
// ============================================================================
int BookshelfParser::ReadNodesFile(string file, PlaceDB &db)
{
	string path;
	gArg.GetString("path", &path);
	path += file;
	ifstream in(path.c_str());
	if (!in)
	{
		cerr << "\tCannot open nodes file: " << file << endl;
		exit(-1);
	}

	// nModules = 总 module 数；nNodes = 可移动数；nTerminals = IO pad 数
	int nModules, nNodes, nTerminals;
	nModules = nNodes = nTerminals = -1;

	int lineNumber = 0;

	// 首行版本头 "UCLA nodes 1.0" 被直接消费掉；严格版本校验写在 #if 0 里，
	// 因此不带版本头的文件同样能被接受。
	// check file format string
	char tmp[10000], tmp2[10000], tmp3[10000];
	in.getline(tmp, 10000);
	lineNumber++;

#if 0
	if( strcmp( "UCLA nodes 1.0", tmp ) != 0 )
	{
		cerr << "Nodes file header error (not UCLA nodes 1.0)\n";
		return -1;
	}
#endif

	// check file header
	// checkFormat 统计头部两个关键字的命中情况，凑齐 2 个才算头部合法
	int checkFormat = 0;
	while (in.getline(tmp, 10000))
	{
		lineNumber++;

		// cout << tmp << endl;
		if (tmp[0] == '#')
			continue;
		// 头部行形如 "NumNodes : 104"：用 strrchr 找最后一个冒号，
		// 再 atoi 冒号之后的数字。
		if (strncmp("NumNodes", tmp, 8) == 0)
		{
			char *pNumber = strrchr(tmp, ':');
			nModules = atoi(pNumber + 1);
			checkFormat++;
		}
		else if (strncmp("NumTerminals", tmp, 12) == 0)
		{
			char *pNumber = strrchr(tmp, ':');
			nTerminals = atoi(pNumber + 1);
			checkFormat++;
		}

		// 头部信息凑齐即可跳出；后面的正文行留给下一个循环去解析
		if (checkFormat == 2)
			break;
	}
	// 可移动单元数 = 总 module 数 - terminal 数
	nNodes = nModules - nTerminals;

	//! 疑似问题：头部不合法时只打印一句警告就继续往下执行，随后
	//! allocateNodeMemory(nNodes) 会拿到 -1 或未初始化的值，vector::resize 到
	//! 一个极大容量，轻则内存耗尽重则崩溃。这里应当像 .pl 那样直接退出。
	if (checkFormat != 2)
	{
		cerr << "** Block file header error (miss NumNodes or NumTerminals)\n";
	}
	// 打印规模信息；可移动单元数超过 1000 时顺便换算成 k 便于阅读
	cout << "    NumModules: " << nModules << endl;
	cout << "    NumNodes: " << nNodes;
	if (nNodes > 1000)

		cout << " (= " << nNodes / 1000 << "k)";
	cout << endl;
	cout << "    Terminals: " << nTerminals << endl;

	// 按头部声明的数量为 dbNodes / dbTerminals 预开槽位，之后按下标写入
	db.allocateNodeMemory(nNodes);
	db.allocateTerminalMemory(nTerminals);

	// Read modules and terminals.
	char name[10000];
	char type[10000];
	double width, height;
	int nodeIndex = 0;
	int terminalIndex = 0;
	int index;
	Module *curModule = NULL;
	// 第二遍扫描：逐行读入 module 正文
	while (in.getline(tmp, 10000))
	{

		lineNumber++;

		// 跳过空行。
		//! 疑似问题：这里没有像头部循环那样跳过以 # 开头的注释行，
		//! .nodes 正文里若出现注释行会被当成一个 module 去解析。
		if (tmp[0] == '\0')
			continue;
		// 每行按 "名字 宽 高 [类型]" 读取，第 4 列可选，所以每轮先把 type 清空。
		//! 疑似问题：tmp2 / tmp3（宽、高）没有同样清空。若某行只有名字（或格式异常），
		//! sscanf 不会写入它们，width / height 就沿用了上一行的残值。
		type[0] = '\0';
		sscanf(tmp, "%s %s %s %s",
			   name, tmp2, tmp3, type);

		// 宽高先以字符串读入，再转成浮点数
		width = atof(tmp2);
		height = atof(tmp3);

		// 下面这段 #if 0 原本用于拒绝行高不一致的混合尺寸（mixed-size）benchmark
#if 0
        if( oldH != -1 && h != oldH && strcmp( type, "terminal" ) != 0)
        {
            cerr << "The program cannot handle mixed-size benchmark currently.";
            exit(0);
        }
        oldH = h;
#endif

		// terminal：固定 IO pad。
		// addTerminal 的最后两个 bool 分别是 isFixed = true（不可移动）、isNI = false
		if (strcmp(type, "terminal") == 0)
		{
			index = terminalIndex;
			curModule = db.addTerminal(terminalIndex, name, width, height, true, false);
			terminalIndex++;
		}
		// terminal_NI：同样是固定 pad，但标记为 NI（do-not-place 类占位，
		// 如 IO filler 等），在密度计算与画图时会被跳过。
		else if (strcmp(type, "terminal_NI") == 0) // (frank) 2022-05-13 consider terminal_NI
		{
			index = terminalIndex;
			curModule = db.addTerminal(terminalIndex, name, width, height, true, true);
			terminalIndex++;
		}
		// 其余情况都当成可移动 node 处理
		else
		{
			index = nodeIndex;
			curModule = db.addNode(nodeIndex, name, width, height);
			nodeIndex++;
		}
		// 登记到名字索引，.nets 里的 pin 名称与 .pl 里的坐标都靠它反查
		//! 疑似问题：若 .nodes 中有重名 module，后出现的会覆盖前者且不报错，
		//! 调用方不会察觉 module 数量与 map 大小不一致。
		db.moduleMap[name] = curModule; // recorded in map
	}

	// 解析完做一致性检查：正文实际数量要与头部声明的 NumNodes / NumTerminals 吻合
	// check if modules number and terminal number match
	if (terminalIndex + nodeIndex != nModules)
	{
		cerr << "Error: There are " << terminalIndex + nodeIndex << " modules in the file\n";
		exit(-1);
	}
	if (terminalIndex != nTerminals)
	{
		cerr << "Error: There are " << terminalIndex << " terminals in the file\n";
		exit(-1);
	}

	// db.moduleCount 记录含 terminal 在内的 module 总数
	db.moduleCount = nNodes + nTerminals;
	return 0;
}

// ============================================================================
// 函数：BookshelfParser::ReadNetsFile
// 职责：解析 .nets，建立 Net 与 Pin，并把它们与 Module 三方关联起来，
//       最后写入 db.dbNets / db.dbPins，并登记 db.netCount / db.pinCount。
//
// 格式（见 examples/smoke/smoke.nets）：
//        NetDegree : 2 n0     <- net 头部：行数/度数，后面的 net 名被忽略
//        c0 B : 0 0           <- pin：module 名、方向(B/I/O)、冒号后的 x/y offset
//
// 处理要点：
//   - 先读头部 NumNets / NumPins 决定 net 与 pin 的容器容量；
//   - 每个 net 先申请 degree 个 pin 的存储，再逐行读入它的 pin；
//   - pin 的 offset 缺省补 0，最终 pin 的绝对位置 = module 位置 + offset。
// ============================================================================
int BookshelfParser::ReadNetsFile(string file, PlaceDB &db)
{
	string path;
	gArg.GetString("path", &path);
	path += file;
	ifstream in(path.c_str());
	if (!in)
	{
		cerr << "Cannot open net file: " << file << endl;
		exit(-1);
	}

	// nNets / nPins 来自头部声明，最后会与实际读到的数量做校验
	int nNets, nPins;
	nNets = nPins = -1;

	int lineNumber = 0;

	// 首行版本头 "UCLA nets 1.0" 直接消费掉，严格校验同样被 #if 0 关掉
	// check file format string
	char tmp[10000];
	in.getline(tmp, 10000);
	lineNumber++;

#if 0
	if( strcmp( "UCLA nets 1.0", tmp ) != 0 )
	{
		cerr << "Nets file header error (UCLA nets 1.0)\n";
		exit(-1);
	}
#endif

	// check file header
	// 与 .nodes 相同：凑齐 NumNets + NumPins 两个关键字才认为头部合法
	int checkFormat = 0;
	while (in.getline(tmp, 10000))
	{
		lineNumber++;

		// cout << tmp << endl;
		if (tmp[0] == '#')
			continue;
		// "NumNets : 184"：取最后一个冒号之后的数字
		if (strncmp("NumNets", tmp, 7) == 0)
		{
			char *pNumber = strrchr(tmp, ':');
			nNets = atoi(pNumber + 1);
			// dbNets 用 resize 预开槽位，之后 dbNets[net->idx] 才能按下标直接赋值
			db.allocateNetMemory(nNets);
			checkFormat++;
		}
		// "NumPins : 368"
		else if (strncmp("NumPins", tmp, 7) == 0)
		{
			char *pNumber = strrchr(tmp, ':');
			nPins = atoi(pNumber + 1);
			// dbPins 用 reserve 而不是 resize：pin 是逐个 push_back 进去的
			db.allocatePinMemory(nPins);
			checkFormat++;
		}

		if (checkFormat == 2)
			break;
	}

	//! 疑似问题：与 .nodes 同样的问题——头部不合法只打印警告而不退出，
	//! 后续 allocateNetMemory / allocatePinMemory 会拿到非法容量。
	if (checkFormat != 2)
	{
		cerr << "** Net file header error\n";
	}

	cout << "         Nets: " << nNets << endl;
	cout << "         Pins: " << nPins << endl;

	// 正文解析用到的 token 缓冲：tmp1/tmp2 用于头两个 token，tmp3/tmp4 用于 pin offset
	char tmp1[2000], tmp2[2000], tmp3[2000], tmp4[2000];
	int maxDegree = 0;
	//! 疑似问题：degree 声明时没有初始化。若某行 sscanf 没能读到第三个字段，
	//! degree 将是不确定值，直接用它做循环次数与扩容规模会很危险。
	int degree;
	int pinIndex = 0;
	int netIndex = 0;
	while (in.getline(tmp, 10000))
	{
		lineNumber++;

		if (tmp[0] == '\0')
			continue;

		// net 头部形如 "NetDegree : 2 n0"：这里只读前 3 个 token，
		// 末尾的 net 名字（n0）被忽略——net 的编号由 netIndex 顺序生成。
		sscanf(tmp, "%s %s %d", tmp1, tmp2, &degree);
		// 头部必须是 NetDegree 且度数非负，否则认为语法不支持并报错返回。
		//! 疑似问题：返回值写成了 01。八进制 01 恰好等于 1，结果没错但容易误解，
		//! 疑是笔误。另外这里 return 时不会释放已 new 出来的 Net（内存泄漏）。
		if (strcmp(tmp1, "NetDegree") != 0 || degree < 0)
		{
			cerr << "Syntax unsupport in line " << lineNumber << ": "
				 << tmp1 << endl;
			return 01;
		}
		// 每个 net 新建一个 Net 对象，编号从 0 开始递增
		Net *net = new Net(netIndex); // default constructer
		Module *module;
		int vCount;
		int pinId;
		double xOffset, yOffset;
		// 顺便统计最大线网度数（maxDegree），可用于估算 net 建模的开销
		if (degree > maxDegree)
			maxDegree = degree;

		// 为该 net 预留 degree 个 pin 的存储空间
		net->allocateMemoryForPin(degree);

		// 累计全局 pin 计数，最后与头部 NumPins 对账
		pinIndex += degree; // will read "degree" pins

		// 紧接着的 degree 行就是该 net 的每个 pin
		for (int j = 0; j < degree; j++)
		{
			// 逐行读 pin。注意这里在内层直接 getline，若文件在 net 中途结束，
			// 读失败后 tmp 会残留上一行内容而不会被察觉。
			in.getline(tmp, 10000);
			lineNumber++;
			tmp3[0] = '\0';
			tmp4[0] = '\0';
			// pin 行形如 "c0 B : 0 0"：tmp1 = module 名，tmp2 = 方向，
			// 冒号后的 tmp3 / tmp4 是 pin 相对 module 左下角的 x / y offset。
			// scanf 格式串里的 " : " 可以匹配任意空白（含零个），所以写成
			// "c0 B:0 0" 或 "c0 B : 0 0" 都能正确切分。
			//! 疑似问题：tmp2 里的 pin 方向（B/I/O）被读出来后丢弃，没有写入
			//! Pin::direction，所有 pin 的 direction 保持为构造时的 -1（未定义）。
			vCount = sscanf(tmp, "%s %s : %s %s", tmp1, tmp2, tmp3, tmp4);

			// 上面已把 tmp3 / tmp4 清空：offset 缺省补 0（等价于落在 module 左下角）
			if (tmp3[0] != '\0')
				xOffset = atof(tmp3);
			else
				xOffset = 0;
			if (tmp4[0] != '\0')
				yOffset = atof(tmp4);
			else
				yOffset = 0;

			// 按名字反查所属 module
			//! 疑似问题：这里没有判空。若 .nets 引用了 .nodes 中不存在的名字，
			//! getModuleFromName 返回 NULL，随后会在 addPin / addPin 里解引用崩溃
			//! （ReadPLFile 中同样的情况是做了检查并报错退出的）。
			module = db.getModuleFromName(tmp1);
			// 新建 Pin，同时把 <module, net, offset> 三者的关联记下来，返回 pin 下标
			pinId = db.addPin(module, net, xOffset, yOffset);
			// module 侧也挂上这个 pin；Module::addPin 会顺带把 net 追加到
			// module->nets（因此不必再手工维护 module 与 net 的关系表）。
			module->addPin(db.dbPins[pinId]);
			net->addPin(db.dbPins[pinId]);

		// 原作者这里留过一段被注释掉的代码，本意是往 module 的 net 表里去重；
		// 由于 Module::addPin 已经自动在维护 nets 列表，这段代码已不再必要。
			// 2005/2/2 (donnie)
			// TODO: Remove duplicate netsIds
			// 2007/3/9 (indark)
			// remove duplicated netsIds
			//? is this necessary?
			// bool found = false;
			// for (unsigned int z = 0; z < db.m_modules[moduleId].m_netsId.size(); z++)
			// {
			// 	if (nReadNets == db.m_modules[moduleId].m_netsId[z])
			// 	{
			// 		found = true;
			// 		break;
			// 	}
			// }
			// if (!found)
			// 	db.m_modules[moduleId].m_netsId.push_back(nReadNets);
		}
		// 把 net 放入 db.dbNets[net->idx]
		db.addNet(net);
		netIndex++;

		// if (nReadNets % stepNet == 0 && nNets > stepNet)
		// 	printf("#%d...\n", nReadNets);
	}

	// TODO: if nReadNets > nNets may have memory problem.
	// check if modules number and terminal number match
	// 一致性检查：实际读入的 net 数 / pin 数要与头部声明一致
	if (nNets != netIndex)
	{
		cerr << "Error: There are " << netIndex << " nets in the file\n";
		exit(-1);
	}
	if (pinIndex != nPins)
	{
		cerr << "Error: There are " << pinIndex << " pins in the file\n";
		exit(-1);
	}

	// 写入 db 的 pin / net 总数，供后续 HPWL 计算与密度模型使用
	db.pinCount = nPins;
	db.netCount = nNets;

	// 打印最大线网度数（#if 1 恒定开启的调试输出）
#if 1
	cout << "Max net degree= " << maxDegree << endl;
#endif
	return 0;
}

// ============================================================================
// 函数：BookshelfParser::ReadPLFile
// 职责：解析 .pl，逐个设置 module 的初始位置（左下角坐标）与方向。
//
// 格式：每行一个 module，例如 "c0 0 0 : N"（名字、x、y、方向），
//       也兼容带括号写成 "c0 (0,0) : N" 的写法。
//
// 参数 init：
//   true  -> 随 .aux 首次读取，file 只是文件名，需要拼上 gArg("path") 前缀；
//   false -> file 已是完整路径，用于把上一轮的布局结果重新读回 PlaceDB。
//
// 字符串处理要点：
//   先把 ( ) , : = \r 等分隔符全部替换成空格，再用 "%s %f %f %s" 取值，
//   这样可以同时兼容 "N : 1 2" 与 "N (1,2)" 两种写法，并处理 Windows 换行。
// ============================================================================
int BookshelfParser::ReadPLFile(string file, PlaceDB &db, bool init)
{
	string path;
	// init 为真时要把文件名拼到 gArg("path") 之后
	if (init)
	{
		gArg.GetString("path", &path);
		path += file;
		cout << "Initialize module position with file: " << file << "\n";
	}
	else
	{
		path = file;
		cout << "Setting module position with file: " << file << "\n";
	}

	ifstream in(path.c_str());
	if (!in)
	{
		cerr << "\tCannot open PL file: " << file << endl;
		exit(-1);
	}

	int lineNumber = 0;

	// 首行是版本头 "UCLA pl 1.0"，这里只是消费掉；
	// 紧跟着被注释掉的一段原本会校验版本头并返回 -1。
	// check file format string
	char tmp[10000];
	in.getline(tmp, 10000);
	lineNumber++;
	// if( strcmp( "UCLA pl 1.0", tmp ) != 0 )
	//{
	//	cerr << "PL file header format error (UCLA pl 1.0)\n";
	//	return -1;
	// }

	char name[10000];
	char orientation[1000];
	float x, y;
	// 逐行解析每个 module 的位置
	while (in.getline(tmp, 10000))
	{
		lineNumber++;

		// cout << tmp << endl;
		// 跳过注释行与空行
		if (tmp[0] == '#')
			continue;
		if (tmp[0] == '\0')
			continue;

		// 把各种分隔符统一替换成空格：既处理 "(x,y)" 的括号与逗号，
		// 也处理方向前的冒号、"x = 1" 的等号，以及 Windows 的回车符 \r。
		for (int i = 0; i < (int)strlen(tmp); i++)
		{
			if (tmp[i] == '(')
				tmp[i] = ' ';
			if (tmp[i] == ',')
				tmp[i] = ' ';
			if (tmp[i] == ')')
				tmp[i] = ' ';
			if (tmp[i] == ':')
				tmp[i] = ' ';
			if (tmp[i] == '=')
				tmp[i] = ' ';
			if (tmp[i] == '\r')
				tmp[i] = ' ';
		}

		// 每轮先把名字缓冲清空，便于后续按返回值判断该行是否成功读到内容
		name[0] = '\0';
		// 期望取到 4 个字段：<名字> <x> <y> <方向>，ret 为实际匹配成功的字段数
		int ret = sscanf(tmp, "%s %f %f %s ", name, &x, &y, orientation);

		// 一个字段都没读到（例如纯空白行），直接跳过
		if (ret <= 0)
			continue;

		// 只接受 3 个字段（无方向）或 4 个字段（带方向）。
		//! 疑似问题：出错时只打印一句提示然后 continue 跳过该行，并不退出；
		//! 被跳过的 module 会保留默认位置 (0,0)，后续布局结果悄悄失真。
		if (ret != 4 && ret != 3)
		{
			// cerr << "Error in the PL file: <" << tmp << ">\n";
			// exit(-1);
			printf("Syntax (may) error in line %d. Please check. (ret = %d)\n",
				   lineNumber, ret);
			continue; // skip this line...
		}

		// 没有方向信息时默认按 N（north）处理
		if (ret == 3)
		{
			// printf( "Block %s does not has orientation.\n", name );
			orientation[0] = 'N';
			orientation[1] = '\0';
		}

		// 按名字找 module；这里的 assert 被注释掉了，改用下面的显式判空报错
		Module *module = db.getModuleFromName(name);
		// assert(module);
			// .pl 里出现了 .nodes 中不存在的名字，属于文件不一致，直接退出
		if (!module)
		{
			cerr << "Error: module name " << name << " not found in line "
				 << lineNumber << " file: " << file << endl;
			exit(-1);
		}

		// 设置左下角坐标。注意 setModuleLocation_2D 会把非固定的 module 裁剪到
		// coreRegion 之内；terminal 因为 isFixed 为真不受裁剪限制。
		db.setModuleLocation_2D(module, x, y);
		// orientInt 把方向字符串转成 enum ORIENT 的整数值。
		//! 疑似问题：orientInt 内部 assert(strlen <= 2) 且 switch 没有 default 分支，
		//! 一旦 .pl 的方向列出现非预期字符串（例如第二个字符是错的），可能走到
		//! 无返回值的路径上，导致返回值不确定。建议这里先校验再传入。
		//! 疑似问题：解析结束后没有检查是否所有 module 都被赋值过位置；
		//! .pl 里缺失的 module 会静默停留在默认位置 (0,0)。
		db.setModuleOrientation(module, orientInt(orientation));
	}

	return 0;
}