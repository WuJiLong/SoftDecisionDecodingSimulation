#ifndef DEFINE_HPP
#define DEFINE_HPP
#define PI 3.14159265358979323846
#define EXP 2.71828182845904523536
#define PolynomialSize 16
#define MASK 0xf
#define T 2
#define _m_ 4
#if MASK == 0x7 && _m_ == 3
#define GFX 0xd
#elif MASK == 0xf && _m_ == 4
#define GFX 0x19
#elif MASK == 0x1f && _m_ == 5
#define GFX 0x25
#elif MASK == 0x3f && _m_ == 6
#define GFX 0x43
#elif MASK == 0x7f && _m_ == 7
#define GFX 0x83
#elif MASK == 0xff && _m_ == 8
#define GFX 0x11D
#else
#define GFX 0
#error "define value setting error!"
#endif

typedef unsigned long long Code;
#define ABS(x) ((x)>0?(x):-(x))
#endif
