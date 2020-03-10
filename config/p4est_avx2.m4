dnld  P4EST_CHECK_AVX2([ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
dnl
dnl   See also SIMD_GCC_CPU_SUPPORTS(INSTRUCTION-SET,
dnl   [ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND])
dnl   which is the actual macro that perform
dnl   the checks for the instruction sets.

AC_DEFUN([SIMD_GCC_CPU_SUPPORTS],
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
           AS_TR_CPP([HAVE_$1_INSTRUCTIONS]),
           [1],
           [Define if $1 instructions are supported])
          $2],
          [$3]
         )
   AS_VAR_POPDEF([simd_cv_gcc_check_cpu_init])
])

AC_DEFUN([P4EST_CHECK_AVX2],
 [m4_foreach_w(
   [simd_feature],
   [avx2],
   [SIMD_GCC_CPU_SUPPORTS(simd_feature,
     [SIMD_FEATURE_CFLAGS="$SIMD_FEATURE_CFLAGS -m[]simd_feature"],
     [])
  ])
  AC_SUBST([SIMD_FEATURE_CFLAGS])
  m4_ifval([$1],[$1],
    [CFLAGS="$CFLAGS $SIMD_FEATURE_CFLAGS"])
])
