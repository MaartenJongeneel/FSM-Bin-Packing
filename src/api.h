#pragma once

#if defined _WIN32 || defined __CYGWIN__
#  define FSMSequenceOne_DLLIMPORT __declspec(dllimport)
#  define FSMSequenceOne_DLLEXPORT __declspec(dllexport)
#  define FSMSequenceOne_DLLLOCAL
#else
// On Linux, for GCC >= 4, tag symbols using GCC extension.
#  if __GNUC__ >= 4
#    define FSMSequenceOne_DLLIMPORT __attribute__((visibility("default")))
#    define FSMSequenceOne_DLLEXPORT __attribute__((visibility("default")))
#    define FSMSequenceOne_DLLLOCAL __attribute__((visibility("hidden")))
#  else
// Otherwise (GCC < 4 or another compiler is used), export everything.
#    define FSMSequenceOne_DLLIMPORT
#    define FSMSequenceOne_DLLEXPORT
#    define FSMSequenceOne_DLLLOCAL
#  endif // __GNUC__ >= 4
#endif // defined _WIN32 || defined __CYGWIN__

#ifdef FSMSequenceOne_STATIC
// If one is using the library statically, get rid of
// extra information.
#  define FSMSequenceOne_DLLAPI
#  define FSMSequenceOne_LOCAL
#else
// Depending on whether one is building or using the
// library define DLLAPI to import or export.
#  ifdef FSMSequenceOne_EXPORTS
#    define FSMSequenceOne_DLLAPI FSMSequenceOne_DLLEXPORT
#  else
#    define FSMSequenceOne_DLLAPI FSMSequenceOne_DLLIMPORT
#  endif // FSMSequenceOne_EXPORTS
#  define FSMSequenceOne_LOCAL FSMSequenceOne_DLLLOCAL
#endif // FSMSequenceOne_STATIC