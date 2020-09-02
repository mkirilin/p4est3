dnl   P4EST_SIMD_GCC_CPU_SUPPORTS(INSTRUCTION-SET,
dnl   PREFIX, ACTION-IF-FOUND, ACTION-IF-NOT-FOUND)
dnl   This macro performs the checks for the avx2 instruction sets
dnl
AC_DEFUN([P4EST_SIMD_GCC_CPU_SUPPORTS],
  [AC_REQUIRE([AC_PROG_CC])
   AC_LANG_PUSH([C])
   AS_VAR_PUSHDEF([simd_cv_gcc_check_cpu_init],
        [AS_TR_SH([simd_cv_gcc_cpu_supports_$1])])
   BACKUP_CFLAGS="${CFLAGS}"
   CFLAGS="${BACKUP_CFLAGS} -m$1"
   AC_CACHE_CHECK([for $1 instruction support],
     [simd_cv_gcc_check_cpu_init],
     [AC_RUN_IFELSE(
       [AC_LANG_PROGRAM( [#include <immintrin.h> ],
       [ __m256i a;
         a = _mm256_set1_epi32 (1);
         _mm256_abs_epi32 (a);
         return 0;
        ])],
        [simd_cv_gcc_check_cpu_init=yes],
        [simd_cv_gcc_check_cpu_init=no],
        [simd_cv_gcc_check_cpu_init=no]
     )]
   )
   CFLAGS="$BACKUP_CFLAGS"
   AC_LANG_POP([C])
   AS_VAR_IF([simd_cv_gcc_check_cpu_init], [yes], [$3], [$4])
   AS_VAR_POPDEF([simd_cv_gcc_check_cpu_init])
])

dnl P4EST_ARG_DISABLE_AVX2(NAME, COMMENT)
dnl Check configure line for --enable/disable-avx2, default enabled.
dnl If enabled, run a compile test with an AVX2 program.
dnl Only if that succeeds, shell variable P4EST_ENABLE_AVX2 stays "yes."
dnl On success, we add -mavx2 to the CFLAGS and also
dnl set the #define and AM_CONDITIONAL P4EST_ENABLE_AVX2.
dnl
AC_DEFUN([P4EST_ARG_DISABLE_AVX2],
[
AC_ARG_ENABLE([$1],
              [AS_HELP_STRING([--disable-$1], [$2])],,
              [enableval=yes])
if test "x$enableval" != xno ; then
  P4EST_SIMD_GCC_CPU_SUPPORTS([avx2], [P4EST], [], [enableval=no])
fi
if test "x$enableval" != xno ; then
  AC_DEFINE([AVX2], 1, [DEPRECATED (use P4EST_ENABLE_AVX2 instead)])
  AC_DEFINE([ENABLE_AVX2], 1, [Undefine if: $2])
  SIMD_FEATURE_CFLAGS="-m[]avx2"
  CFLAGS="$CFLAGS $SIMD_FEATURE_CFLAGS"
fi
AM_CONDITIONAL([P4EST_ENABLE_AVX2], [test "x$enableval" != xno])
P4EST_ENABLE_AVX2="$enableval"
])
