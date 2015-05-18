#ifndef Vector_ispc_h
#define Vector_ispc_h

struct float3 { float x, y, z; };

static inline float3 operator+(float3 a, float3 b) { float3 r = { a.x + b.x, a.y + b.y, a.z + b.z }; return r; }
static inline float3 operator-(float3 a, float3 b) { float3 r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static inline float3 operator*(float3 a, float3 b) { float3 r = { a.x * b.x, a.y * b.y, a.z * b.z }; return r; }
static inline float3 operator/(float3 a, float3 b) { float3 r = { a.x / b.x, a.y / b.y, a.z / b.z }; return r; }
static inline float3 operator+(float3 a, float b) { float3 r = { a.x + b, a.y + b, a.z + b }; return r; }
static inline float3 operator-(float3 a, float b) { float3 r = { a.x - b, a.y - b, a.z - b }; return r; }
static inline float3 operator*(float3 a, float b) { float3 r = { a.x * b, a.y * b, a.z * b }; return r; }
static inline float3 operator/(float3 a, float b) { float3 r = { a.x / b, a.y / b, a.z / b }; return r; }
static inline float3 operator+(float a, float3 b) { float3 r = { a + b.x, a + b.y, a + b.z }; return r; }
static inline float3 operator-(float a, float3 b) { float3 r = { a - b.x, a - b.y, a - b.z }; return r; }
static inline float3 operator*(float a, float3 b) { float3 r = { a * b.x, a * b.y, a * b.z }; return r; }
static inline float3 operator/(float a, float3 b) { float3 r = { a / b.x, a / b.y, a / b.z }; return r; }
static inline uniform float3 operator+(uniform float3 a, uniform float3 b) { uniform float3 r = { a.x + b.x, a.y + b.y, a.z + b.z }; return r; }
static inline uniform float3 operator-(uniform float3 a, uniform float3 b) { uniform float3 r = { a.x - b.x, a.y - b.y, a.z - b.z }; return r; }
static inline uniform float3 operator*(uniform float3 a, uniform float3 b) { uniform float3 r = { a.x * b.x, a.y * b.y, a.z * b.z }; return r; }
static inline uniform float3 operator/(uniform float3 a, uniform float3 b) { uniform float3 r = { a.x / b.x, a.y / b.y, a.z / b.z }; return r; }
static inline uniform float3 operator+(uniform float3 a, uniform float b) { uniform float3 r = { a.x + b, a.y + b, a.z + b }; return r; }
static inline uniform float3 operator-(uniform float3 a, uniform float b) { uniform float3 r = { a.x - b, a.y - b, a.z - b }; return r; }
static inline uniform float3 operator*(uniform float3 a, uniform float b) { uniform float3 r = { a.x * b, a.y * b, a.z * b }; return r; }
static inline uniform float3 operator/(uniform float3 a, uniform float b) { uniform float3 r = { a.x / b, a.y / b, a.z / b }; return r; }
static inline uniform float3 operator+(uniform float a, uniform float3 b) { uniform float3 r = { a + b.x, a + b.y, a + b.z }; return r; }
static inline uniform float3 operator-(uniform float a, uniform float3 b) { uniform float3 r = { a - b.x, a - b.y, a - b.z }; return r; }
static inline uniform float3 operator*(uniform float a, uniform float3 b) { uniform float3 r = { a * b.x, a * b.y, a * b.z }; return r; }
static inline uniform float3 operator/(uniform float a, uniform float3 b) { uniform float3 r = { a / b.x, a / b.y, a / b.z }; return r; }

static inline uniform float3 reduce_add(float3 v)
{
    uniform float3 r = { reduce_add(v.x), reduce_add(v.y), reduce_add(v.z) };
    return r;
}

static inline float dot(float3 a, float3 b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}
static inline uniform float dot(uniform float3 a, uniform float3 b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

static inline float length_sq(float3 v)
{
    return dot(v, v);
}
static inline uniform float length_sq(uniform float3 v)
{
    return dot(v, v);
}
static inline float length(float3 v)
{
    return sqrt(length_sq(v));
}
static inline uniform float length(uniform float3 v)
{
    return sqrt(length_sq(v));
}

#endif // Vector_ispc_h
