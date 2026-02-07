#ifndef EXDQLM_OMP_COMPAT_H
#define EXDQLM_OMP_COMPAT_H

// If OpenMP is available, include it; otherwise provide minimal stubs.
#ifdef _OPENMP
  #include <omp.h>
#else
  inline int  omp_get_max_threads() { return 1; }
  inline int  omp_get_num_threads() { return 1; }
  inline int  omp_get_thread_num()  { return 0; }
  inline void omp_set_num_threads(int) {}
  // If you use omp_get_wtime() anywhere, also add:
  // #include <chrono>
  // inline double omp_get_wtime() {
  //   using clk = std::chrono::steady_clock;
  //   static const auto t0 = clk::now();
  //   return std::chrono::duration<double>(clk::now() - t0).count();
  // }
#endif

#endif // EXDQLM_OMP_COMPAT_H
