#ifndef Vector_h
#define Vector_h

#include <cmath>
#include <intrin.h>
#define SSE_SHUFFLE(x,y,z,w) _MM_SHUFFLE(w,z,y,x)


struct float3
{
    union {
        struct { float x, y, z; };
        float v[3];
    };
};

struct float4
{
    union {
        struct { float x, y, z, w; };
        float v[4];
    };
};

inline float3 operator+(float3 a, float3 b) { float3 r = { a.x + b.x, a.y + b.y, a.z + b.z }; return r; }
inline float3 operator-(float3 a, float3 b) { float3 r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
inline float3 operator*(float3 a, float3 b) { float3 r = { a.x * b.x, a.y * b.y, a.z * b.z }; return r; }
inline float3 operator/(float3 a, float3 b) { float3 r = { a.x / b.x, a.y / b.y, a.z / b.z }; return r; }
inline float3 operator+(float3 a, float b) { float3 r = { a.x + b, a.y + b, a.z + b }; return r; }
inline float3 operator-(float3 a, float b) { float3 r = { a.x - b, a.y - b, a.z - b }; return r; }
inline float3 operator*(float3 a, float b) { float3 r = { a.x * b, a.y * b, a.z * b }; return r; }
inline float3 operator/(float3 a, float b) { float3 r = { a.x / b, a.y / b, a.z / b }; return r; }
inline float3 operator+(float a, float3 b) { float3 r = { a + b.x, a + b.y, a + b.z }; return r; }
inline float3 operator-(float a, float3 b) { float3 r = { a - b.x, a - b.y, a - b.z }; return r; }
inline float3 operator*(float a, float3 b) { float3 r = { a * b.x, a * b.y, a * b.z }; return r; }
inline float3 operator/(float a, float3 b) { float3 r = { a / b.x, a / b.y, a / b.z }; return r; }

inline float dot(float3 a, float3 b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline float length_sq(float3 v)
{
    return dot(v, v);
}
inline float length(float3 v)
{
    return sqrt(length_sq(v));
}



typedef __m128 simdfloat4;

struct soa43
{
    union {
        simdfloat4 v[3];
        struct { simdfloat4 x, y, z; };
    };
};

struct soa44
{
    union {
        simdfloat4 v[4];
        struct { simdfloat4 x, y, z, w; };
    };
};

inline soa43 soa_transpose34(const simdfloat4 &v0, const simdfloat4 &v1, const simdfloat4 &v2, const simdfloat4 &v3)
{
    simdfloat4 r1 = _mm_unpacklo_ps(v0, v1);
    simdfloat4 r2 = _mm_unpacklo_ps(v2, v3);
    simdfloat4 r3 = _mm_unpackhi_ps(v0, v1);
    simdfloat4 r4 = _mm_unpackhi_ps(v2, v3);
    soa43 ret = {
        _mm_shuffle_ps(r1, r2, SSE_SHUFFLE(0, 1, 0, 1)),
        _mm_shuffle_ps(r1, r2, SSE_SHUFFLE(2, 3, 2, 3)),
        _mm_shuffle_ps(r3, r4, SSE_SHUFFLE(0, 1, 0, 1)) };
    return ret;
}

inline soa44 soa_transpose44(const simdfloat4 &v0, const simdfloat4 &v1, const simdfloat4 &v2)
{
    simdfloat4 zero = _mm_set_ps1(0.0f);
    simdfloat4 r1 = _mm_unpacklo_ps(v0, v1);
    simdfloat4 r2 = _mm_unpacklo_ps(v2, zero);
    simdfloat4 r3 = _mm_unpackhi_ps(v0, v1);
    simdfloat4 r4 = _mm_unpackhi_ps(v2, zero);
    soa44 ret = {
        _mm_shuffle_ps(r1, r2, SSE_SHUFFLE(0, 1, 0, 1)),
        _mm_shuffle_ps(r1, r2, SSE_SHUFFLE(2, 3, 2, 3)),
        _mm_shuffle_ps(r3, r4, SSE_SHUFFLE(0, 1, 0, 1)),
        _mm_shuffle_ps(r3, r4, SSE_SHUFFLE(2, 3, 2, 3)) };
    return ret;
}

#endif // Vector_h
