#ifndef LAS_INLINE_H_
#define LAS_INLINE_H_
// check both clang and intel before gnu
#if defined(__clang__)
#define LAS_INLINE inline
#elif defined(__INTEL_COMPILER)
#define LAS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__GNUG__)
#define LAS_INLINE __attribute__((always_inline))
#elif defined(_MSC_VER)
#define LAS_INLINE inline
#endif
#endif
