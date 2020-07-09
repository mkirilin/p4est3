dnl  P4EST_CHECK_AVX2([PREFIX])
dnl  If --disable-avx2 feature is not set:
dnl  chech if target CPU supports SIMD avx2 intrinsics.
dnl  If it is supported, set conditional PREFIX_HAVE_AVX2.
dnl  If it is supported, set valid right-hand side for a
dnl  C #define PREFIX_HAVE_AVX2_INSTRUCTIONS as well.
dnl
dnl   See also P4EST_SIMD_GCC_CPU_SUPPORTS(INSTRUCTION-SET,
dnl   [PREFIX], [ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
dnl   which is the actual macro that perform
dnl   the checks for the instruction sets.

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
   AS_VAR_IF([simd_cv_gcc_check_cpu_init],[yes],
         [AC_DEFINE(
           AS_TR_CPP([$2_HAVE_$$1_INSTRUCTIONS]),
           [1],
           [Define if $1 instructions are supported])
          $3],
          [$4]
         )
   AS_VAR_POPDEF([simd_cv_gcc_check_cpu_init])
])

AC_DEFUN([P4EST_CHECK_AVX2],
 [ AM_CONDITIONAL([$1_HAVE_AVX2], [test "xyes" != xno])
   AM_COND_IF([P4EST_ENABLE_AVX2],[
     P4EST_SIMD_GCC_CPU_SUPPORTS(avx2, [$1],
      [SIMD_FEATURE_CFLAGS="-m[]avx2"
       $1_HAVE_AVX2="yes"],
      [$1_HAVE_AVX2="no"]
     )
     AC_SUBST([SIMD_FEATURE_CFLAGS])
     [CFLAGS="$CFLAGS $SIMD_FEATURE_CFLAGS"]
     ],[])
])
