// ============================================================================
// 模块总览：2D DCT / DST 封装（Ooura fftsg2d 包的 C++ 移植版）
// ----------------------------------------------------------------------------
// 【职责】
//   为 FFT/fft.cpp 的 FFT_2D 提供三个 2D 变换：
//     ddct2d  ：X 方向 cos 基 + Y 方向 cos 基（密度 ρ 的正变换、电势 φ 的反变换）
//     ddsct2d ：X 方向 sin 基 + Y 方向 cos 基（电场 E_x 的反变换）
//     ddcst2d ：X 方向 cos 基 + Y 方向 sin 基（电场 E_y 的反变换）
//   上层语义（泊松方程频域求解）见 FFT/fft.h 与 FFT/fft.cpp，本文件只管「怎么算」。
//
// 【2D 如何复用 1D】—— 行列分离法（row-column）
//   2D 变换是可分离的，可以拆成两轮 1D 变换，完全不需要写 2D 的蝶形核：
//     第一轮：对每一行 a[i]（长度 n2，即 Y 方向）调用 1D 的 ddct / ddst；
//     第二轮：把「列」抽到临时缓冲 t 里（长度 n1，即 X 方向）再做 1D 变换，然后写回。
//   这就是本文件全部四个函数做的事：三个入口函数各自先跑第一轮（行），
//   再调用 ddxt2d_sub() 跑第二轮（列），靠参数 ics 决定第二轮用 cos 基还是 sin 基。
//
// 【关键参数】
//   n1 / n2：两个维度的长度（本项目即 binCntX / binCntY），Ooura 要求为 2 的幂
//   isgn  ：−1 正变换、+1 反变换
//   ics   ：第二轮（n1 方向）用哪种基，0 = DCT（cos）、1 = DST（sin）
//   a     ：a[n1][n2] 二维数组，原地（in-place）输入输出
//   t     ：列变换用的临时缓冲，传 NULL 时由库内部分配并在返回前释放
//   ip / w：位反转工作表与 cos/sin 三角函数表（由 FFT_2D 持有，跨调用复用）
//
// 【命名规则】dd[X][Y]2d，第一个字母对应 n1（X），第二个字母对应 n2（Y），
//   c = cos 基（DCT）、s = sin 基（DST）。之所以要区分 sin/cos，是因为
//   ePlace 里对电势求导得到电场时，「在哪个方向求导，那个方向的基就变成 sin」。
//
//! 疑似问题：若编译时定义了 USE_FFT2D_PTHREADS / USE_FFT2D_WINTHREADS，
//! 大数据量分支会调用 ddxt2d0_subth() 与 ddxt2d_subth()，但这两个函数在本目录的
//! 任何文件里都没有定义（Ooura 原版定义在 fftsg2d.c 中，这里移植时被裁掉了），
//! 会得到链接错误。本项目默认不开启这两个宏，走的是下面的单线程分支。
// ============================================================================

#include <stdio.h>
#include <stdlib.h>
#include "fft.h"
#define fft2d_alloc_error_check(p)                        \
  {                                                       \
    if((p) == NULL) {                                     \
      fprintf(stderr, "fft2d memory allocation error\n"); \
      exit(1);                                            \
    }                                                     \
  }

#ifdef USE_FFT2D_PTHREADS
#define USE_FFT2D_THREADS
#ifndef FFT2D_MAX_THREADS
#define FFT2D_MAX_THREADS 4
#endif
#ifndef FFT2D_THREADS_BEGIN_N
#define FFT2D_THREADS_BEGIN_N 65536
#endif
#include <pthread.h>
#define fft2d_thread_t pthread_t
#define fft2d_thread_create(thp, func, argp)                   \
  {                                                            \
    if(pthread_create(thp, NULL, func, (void *)(argp)) != 0) { \
      fprintf(stderr, "fft2d thread error\n");                 \
      exit(1);                                                 \
    }                                                          \
  }
#define fft2d_thread_wait(th)                  \
  {                                            \
    if(pthread_join(th, NULL) != 0) {          \
      fprintf(stderr, "fft2d thread error\n"); \
      exit(1);                                 \
    }                                          \
  }
#endif /* USE_FFT2D_PTHREADS */

#ifdef USE_FFT2D_WINTHREADS
#define USE_FFT2D_THREADS
#ifndef FFT2D_MAX_THREADS
#define FFT2D_MAX_THREADS 4
#endif
#ifndef FFT2D_THREADS_BEGIN_N
#define FFT2D_THREADS_BEGIN_N 131072
#endif
#include <windows.h>
#define fft2d_thread_t HANDLE
#define fft2d_thread_create(thp, func, argp)                       \
  {                                                                \
    DWORD thid;                                                    \
    *(thp) = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)(func), \
                          (LPVOID)(argp), 0, &thid);               \
    if(*(thp) == 0) {                                              \
      fprintf(stderr, "fft2d thread error\n");                     \
      exit(1);                                                     \
    }                                                              \
  }
#define fft2d_thread_wait(th)          \
  {                                    \
    WaitForSingleObject(th, INFINITE); \
    CloseHandle(th);                   \
  }
#endif /* USE_FFT2D_WINTHREADS */

namespace replace {


// ddcst2d：X(n1) 方向用 cos 基、Y(n2) 方向用 sin 基的 2D 变换
//   （ePlace 中用于把 Ê_y 反变换成空域电场 E_y）
void ddcst2d(int n1, int n2, int isgn, float **a, float *t, int *ip, float *w) {
#ifdef USE_FFT2D_THREADS
  void ddxt2d0_subth(int n1, int n2, int ics, int isgn, float **a, int *ip,
                     float *w);
  void ddxt2d_subth(int n1, int n2, int ics, int isgn, float **a, float *t,
                    int *ip, float *w);
#endif 
  int n, nw, nc, itnull, nthread, nt, i;

  // 两个入口函数共用的准备工作（ddsct2d / ddct2d 中的这段完全相同）：
  // 取两维的最大长度 n，作为三角函数表的规模依据
  n = n1;
  if(n < n2) {
    n = n2;
  }
  // ip[0] 缓存「已生成的旋转因子表规模 nw」，不够大就重建（makewt）
  nw = ip[0];
  if(n > (nw << 2)) {
    nw = n >> 2;
    makewt(nw, ip, w);
  }
  // ip[1] 缓存「已生成的余弦表规模 nc」，不够大就重建（makect，存在 w 的后半段）
  nc = ip[1];
  if(n > nc) {
    nc = n;
    makect(nc, ip, w + nw);
  }
  // t 为列变换的临时缓冲；传 NULL（上层 FFT_2D 就是这么传的）时在这里分配，
  // 并用 itnull 标记「返回前需要由本函数释放」
  itnull = 0;
  if(t == NULL) {
    itnull = 1;
    nthread = 1;
#ifdef USE_FFT2D_THREADS
    nthread = FFT2D_MAX_THREADS;
#endif 
    nt = 4 * nthread * n1;
    if(n2 == 2 * nthread) {
      nt >>= 1;
    }
    else if(n2 < 2 * nthread) {
      nt >>= 2;
    }
    t = (float *)malloc(sizeof(float) * nt);
    fft2d_alloc_error_check(t);
  }
#ifdef USE_FFT2D_THREADS
  if((float)n1 * n2 >= (float)FFT2D_THREADS_BEGIN_N) {
    ddxt2d0_subth(n1, n2, 1, isgn, a, ip, w);
    ddxt2d_subth(n1, n2, 0, isgn, a, t, ip, w);
  }
  else
#endif 
  {
    // 第一轮：对每一行（Y 方向，长度 n2）做 1D 的 DST —— 对应函数名里的第二个字母 s
    for(i = 0; i < n1; i++) {
      ddst(n2, isgn, a[i], ip, w);
    }
    // 第二轮：对每一列（X 方向，长度 n1）做 1D 变换，ics = 0 表示用 DCT（cos 基）
    ddxt2d_sub(n1, n2, 0, isgn, a, t, ip, w);
  }
  if(itnull != 0) {
    free(t);
  }
}

//1
// ddsct2d：X(n1) 方向用 sin 基、Y(n2) 方向用 cos 基的 2D 变换
//   （ePlace 中用于把 Ê_x 反变换成空域电场 E_x；与 ddcst2d 的区别只是两个方向互换了）
//   准备工作的注释见 ddcst2d，此处不再重复
void ddsct2d(int n1, int n2, int isgn, float **a, float *t, int *ip, float *w) {
  #ifdef USE_FFT2D_THREADS
  void ddxt2d0_subth(int n1, int n2, int ics, int isgn, float **a, int *ip,
                     float *w);
  void ddxt2d_subth(int n1, int n2, int ics, int isgn, float **a, float *t,
                    int *ip, float *w);
#endif /* USE_FFT2D_THREADS */
  int n, nw, nc, itnull, nthread, nt, i;

  n = n1;
  if(n < n2) {
    n = n2;
  }
  nw = ip[0];
  if(n > (nw << 2)) {
    nw = n >> 2;
    makewt(nw, ip, w);
  }
  nc = ip[1];
  if(n > nc) {
    nc = n;
    makect(nc, ip, w + nw);
  }
  itnull = 0;
  if(t == NULL) {
    itnull = 1;
    nthread = 1;
#ifdef USE_FFT2D_THREADS
    nthread = FFT2D_MAX_THREADS;
#endif /* USE_FFT2D_THREADS */
    nt = 4 * nthread * n1;
    if(n2 == 2 * nthread) {
      nt >>= 1;
    }
    else if(n2 < 2 * nthread) {
      nt >>= 2;
    }
    t = (float *)malloc(sizeof(float) * nt);
    fft2d_alloc_error_check(t);
  }
#ifdef USE_FFT2D_THREADS
  if((float)n1 * n2 >= (float)FFT2D_THREADS_BEGIN_N) {
    ddxt2d0_subth(n1, n2, 0, isgn, a, ip, w);
    ddxt2d_subth(n1, n2, 1, isgn, a, t, ip, w);
  }
  else
#endif /* USE_FFT2D_THREADS */
  {
    // 第一轮：每一行（Y 方向）用 DCT（cos 基）
    for(i = 0; i < n1; i++) {
      ddct(n2, isgn, a[i], ip, w);
    }
    // 第二轮：每一列（X 方向）用 DST（sin 基），ics = 1
    ddxt2d_sub(n1, n2, 1, isgn, a, t, ip, w);
  }
  if(itnull != 0) {
    free(t);
  }
}


//1
// ddct2d：两个方向都用 cos 基的 2D DCT
//   （ePlace 中用于密度 ρ 的正变换与电势 φ 的反变换，是本模块用得最多的一个）
void ddct2d(int n1, int n2, int isgn, float **a, float *t, int *ip, float *w) {
#ifdef USE_FFT2D_THREADS
  void ddxt2d0_subth(int n1, int n2, int ics, int isgn, float **a, int *ip,
                     float *w);
  void ddxt2d_subth(int n1, int n2, int ics, int isgn, float **a, float *t,
                    int *ip, float *w);
#endif /* USE_FFT2D_THREADS */
  int n, nw, nc, itnull, nthread, nt, i;

  n = n1;
  if(n < n2) {
    n = n2;
  }
  nw = ip[0];
  if(n > (nw << 2)) {
    nw = n >> 2;
    makewt(nw, ip, w);
  }
  nc = ip[1];
  if(n > nc) {
    nc = n;
    makect(nc, ip, w + nw);
  }
  itnull = 0;
  if(t == NULL) {
    itnull = 1;
    nthread = 1;
#ifdef USE_FFT2D_THREADS
    nthread = FFT2D_MAX_THREADS;
#endif /* USE_FFT2D_THREADS */
    nt = 4 * nthread * n1;
    if(n2 == 2 * nthread) {
      nt >>= 1;
    }
    else if(n2 < 2 * nthread) {
      nt >>= 2;
    }
    t = (float *)malloc(sizeof(float) * nt);
    fft2d_alloc_error_check(t);
  }
#ifdef USE_FFT2D_THREADS
  if((float)n1 * n2 >= (float)FFT2D_THREADS_BEGIN_N) {
    ddxt2d0_subth(n1, n2, 0, isgn, a, ip, w);
    ddxt2d_subth(n1, n2, 0, isgn, a, t, ip, w);
  }
  else
#endif /* USE_FFT2D_THREADS */
  {
    // 第一轮：每一行（Y 方向）用 DCT
    for(i = 0; i < n1; i++) {
      ddct(n2, isgn, a[i], ip, w);
    }
    // 第二轮：每一列（X 方向）也用 DCT，ics = 0
    ddxt2d_sub(n1, n2, 0, isgn, a, t, ip, w);
  }
  if(itnull != 0) {
    free(t);
  }
}

//1
// ddxt2d_sub：行列分离法的第二轮 —— 沿 n1（X）方向做 1D 变换，是三个 2D 入口的公共尾巴
//   ics = 0 → 该方向用 DCT（cos 基）；ics = 1 → 该方向用 DST（sin 基）
//   做法：把 n1 长度的「列」从 a 中抽到连续缓冲 t 再做 1D 变换（1D 变换要求连续内存），
//         每次同时取 4 列一起搬，既摊薄搬运开销，也让 4 次 1D 变换的数据都在缓存里
//   t 的布局：[ 第 j 列 | 第 j+1 列 | 第 j+2 列 | 第 j+3 列 ]，每段长度 n1
void ddxt2d_sub(int n1, int n2, int ics, int isgn, float **a, float *t, int *ip,
                float *w) {
  int i, j;

  //! 疑似问题：下面的循环按 j += 4 步进，隐含要求 n2 是 4 的倍数（或恰好为 2）。
  //! Ooura 原库要求各维长度为 2 的幂，ePlace 的 bin 数也满足该条件；
  //! 但若传入其它 n2（如 3、5、6），最后一次迭代会访问 a[i][j+1..j+3] 而越界。
  if(n2 > 2) {
    for(j = 0; j < n2; j += 4) {
      for(i = 0; i < n1; i++) {
        t[i] = a[i][j];
        t[n1 + i] = a[i][j + 1];
        t[2 * n1 + i] = a[i][j + 2];
        t[3 * n1 + i] = a[i][j + 3];
      }
      // 对这 4 列各做一次 1D 变换（长度 n1，沿 X 方向）；ics 决定用 cos 基还是 sin 基
      if(ics == 0) {
        ddct(n1, isgn, t, ip, w);
        ddct(n1, isgn, &t[n1], ip, w);
        ddct(n1, isgn, &t[2 * n1], ip, w);
        ddct(n1, isgn, &t[3 * n1], ip, w);
      }
      else {
        ddst(n1, isgn, t, ip, w);
        ddst(n1, isgn, &t[n1], ip, w);
        ddst(n1, isgn, &t[2 * n1], ip, w);
        ddst(n1, isgn, &t[3 * n1], ip, w);
      }
      // 把变换后的 4 列写回 a 的对应列
      for(i = 0; i < n1; i++) {
        a[i][j] = t[i];
        a[i][j + 1] = t[n1 + i];
        a[i][j + 2] = t[2 * n1 + i];
        a[i][j + 3] = t[3 * n1 + i];
      }
    }
  }
  // n2 == 2 的特例：只有两列，单独处理（无法凑满 4 列）
  else if(n2 == 2) {
    for(i = 0; i < n1; i++) {
      t[i] = a[i][0];
      t[n1 + i] = a[i][1];
    }
    if(ics == 0) {
      ddct(n1, isgn, t, ip, w);
      ddct(n1, isgn, &t[n1], ip, w);
    }
    else {
      ddst(n1, isgn, t, ip, w);
      ddst(n1, isgn, &t[n1], ip, w);
    }
    for(i = 0; i < n1; i++) {
      a[i][0] = t[i];
      a[i][1] = t[n1 + i];
    }
  }
}
}
