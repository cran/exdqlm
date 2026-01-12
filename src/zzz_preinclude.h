#ifndef EXDQLM_ZZZ_PREINCLUDE_H
#define EXDQLM_ZZZ_PREINCLUDE_H

// Forced into every C++ translation unit via:
//   PKG_CPPFLAGS = ... -include zzz_preinclude.h
//
// Purpose:
// Alma/RHEL hardening may define _GLIBCXX_ASSERTIONS (e.g. via -Wp,-D_GLIBCXX_ASSERTIONS).
// When enabled, libstdc++ headers compile in assertion helpers that can reference
// abort() and printf(), which triggers R CMD check --as-cran compiled-code warnings.
//
// Fix:
// Undefine _GLIBCXX_ASSERTIONS before any system/libstdc++ headers are processed.

#ifdef _GLIBCXX_ASSERTIONS
#  undef _GLIBCXX_ASSERTIONS
#endif

// Optional belt-and-suspenders: ensure <cassert> assert() compiles out.
// Keep this if your toolchain ever injects -UNDEBUG (debug builds).
#ifndef NDEBUG
#  define NDEBUG 1
#endif

#ifndef BOOST_DISABLE_ASSERTS
#  define BOOST_DISABLE_ASSERTS
#endif

#endif // EXDQLM_ZZZ_PREINCLUDE_H
