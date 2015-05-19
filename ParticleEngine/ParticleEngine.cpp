#include <cstdlib>
#include <algorithm>
#include <random>
#include <ppl.h>
#include "ParticleEngine.h"
#include "ParticleCore_ispc.h"

using namespace concurrency;

void* peAlignedAlloc(size_t size, size_t align)
{
#ifdef _MSC_VER
    return _aligned_malloc(size, align);
#elif defined(__APPLE__)
    return malloc(size);
#else  // _MSC_VER
    return memalign(align, size);
#endif // _MSC_VER
}

void peAlignedFree(void *p)
{
#ifdef _MSC_VER
    _aligned_free(p);
#else  // _MSC_VER
    free(p);
#endif // _MSC_VER
}




peCopyToTextureBase *g_CopyToTexture;

class peContext
{
public:
    peContext(int n);
    ~peContext();

    void setUpdateRoutine(peUpdateRoutine v);
    void enbaleMultiThreading(bool v);
    void setParticleSize(float v);

    void update(float dt);
    void copyDataToTexture(void *texture, int width, int height);

private:
    void updateVelocity_Plain(float dt, int begin ,int end);
    void updatePosition_Plain(float dt, int begin, int end);
    void update_Plain(float dt);

    void updateVelocity_SIMD(float dt, int begin, int end);
    void updatePosition_SIMD(float dt, int begin, int end);
    void update_SIMD(float dt);

    void updateVelocity_SIMDSoA(float dt, int begin, int end);
    void updatePosition_SIMDSoA(float dt, int begin, int end);
    void update_SIMDSoA(float dt);

    void update_ISPC(float dt);

private:
    peParticle *m_particles;
    peParticleSoA m_soa;
    int m_particle_count;

    float m_particle_size;
    float m_pressure_stiffness;

    peUpdateRoutine m_routine;
    bool m_multi_threading;
    int m_task_granularity;
};


peContext::peContext(int n)
    : m_particle_count(n)
    , m_particle_size(0.2f)
    , m_pressure_stiffness(500.0f)
    , m_routine(peE_ISPC)
    , m_multi_threading(true)
    , m_task_granularity(256)
{
    m_particles = (peParticle*)peAlignedAlloc(sizeof(peParticle)*n);
    m_soa.pos_x = (float*)peAlignedAlloc(sizeof(float)*n);
    m_soa.pos_y = (float*)peAlignedAlloc(sizeof(float)*n);
    m_soa.pos_z = (float*)peAlignedAlloc(sizeof(float)*n);
    m_soa.vel_x = (float*)peAlignedAlloc(sizeof(float)*n);
    m_soa.vel_y = (float*)peAlignedAlloc(sizeof(float)*n);
    m_soa.vel_z = (float*)peAlignedAlloc(sizeof(float)*n);

    // set random position
    std::mt19937 rand;
    std::uniform_real_distribution<float> dist(-5.0f, 5.0f);
    const float4 zero4 = {0.0f};
    for (int i = 0; i < n; ++i) {
        float4 pos = { dist(rand), dist(rand) + 5.0f, dist(rand), 0.0f };
        m_particles[i].position = pos;
        m_particles[i].velocity = zero4;
    }
}

peContext::~peContext()
{
    peAlignedFree(m_particles);
    peAlignedFree(m_soa.pos_x);
    peAlignedFree(m_soa.pos_y);
    peAlignedFree(m_soa.pos_z);
    peAlignedFree(m_soa.vel_x);
    peAlignedFree(m_soa.vel_y);
    peAlignedFree(m_soa.vel_z);
}


void peContext::setUpdateRoutine(peUpdateRoutine v)
{
    m_routine = v;
}

void peContext::enbaleMultiThreading(bool v)
{
    m_multi_threading = v;
}

void peContext::setParticleSize(float v)
{
    m_particle_size = v;
}

void peContext::copyDataToTexture(void *texture, int width, int height)
{
    if (g_CopyToTexture && texture) {
        g_CopyToTexture->copy(texture, width, height, m_particles, sizeof(peParticle)*m_particle_count);
    }
}

void peContext::update(float dt)
{
    switch (m_routine)
    {
    case peE_Plain:     update_Plain(dt);   break;
    case peE_SIMD:      update_SIMD(dt);    break;
    case peE_SIMDSoA:   update_SIMDSoA(dt); break;
    case peE_ISPC:      update_ISPC(dt);    break;
    case peE_CSharp: break;
    default: break;
    }
}


// Plain C++ implementation

void peContext::updateVelocity_Plain(float dt, int begin, int end)
{
    const float particle_size2 = m_particle_size * 2.0f;
    const float rcp_particle_size2 = 1.0f / (m_particle_size * 2.0f);

    // パーティクル同士の押し返し
    for (int i = begin; i < end; ++i) {
        float3 &pos1 = (float3&)m_particles[i].position;
        float3 accel = { 0.0f, 0.0f, 0.0f };
        for (int j = 0; j < m_particle_count; ++j) {
            float3 &pos2 = (float3&)m_particles[j].position;
            float3 diff = pos2 - pos1;
            float3 dir = diff * rcp_particle_size2;
            float dist = length(diff);
            if (dist > 0.0f) {
                float3 a = dir * (std::min<float>(0.0f, dist - particle_size2) * m_pressure_stiffness);
                accel = accel + a;
            }
        }

        float3 &vel = (float3&)m_particles[i].velocity;
        vel = vel + accel * dt;
    }

    // 床との衝突
    const float3 floor_normal = { 0.0f, 1.0f, 0.0f };
    const float floor_distance = -m_particle_size;
    for (int i = begin; i < end; ++i) {
        float3 &pos = (float3&)m_particles[i].position;
        float d = dot(pos, floor_normal) + floor_distance;
        float3 accel = floor_normal * (-std::min<float>(0.0f, d) * m_pressure_stiffness);
        float3 &vel = (float3&)m_particles[i].velocity;
        vel = vel + accel * dt;
    }

    // 重力加速
    const float3 gravity_direction = { 0.0f, -1.0f, 0.0f };
    const float gravity_strength = 5.0f;
    for (int i = begin; i < end; ++i) {
        float3 accel = gravity_direction * gravity_strength;
        float3 &vel = (float3&)m_particles[i].velocity;
        vel = vel + accel * dt;
    }
}

void peContext::updatePosition_Plain(float dt, int begin, int end)
{
    for (int i = begin; i < end; ++i) {
        float3 &pos = (float3&)m_particles[i].position;
        float3 &vel = (float3&)m_particles[i].velocity;
        pos = pos + (vel * dt);
    }
}

void peContext::update_Plain(float dt)
{
    if (m_multi_threading) {
        parallel_for(0, m_particle_count, m_task_granularity, [&](int begin){
            int end = std::min<int>(m_particle_count, begin + m_task_granularity);
            updateVelocity_Plain(dt, begin, end);
        });
        parallel_for(0, m_particle_count, m_task_granularity, [&](int begin){
            int end = std::min<int>(m_particle_count, begin + m_task_granularity);
            updatePosition_Plain(dt, begin, end);
        });
    }
    else {
        updateVelocity_Plain(dt, 0, m_particle_count);
        updatePosition_Plain(dt, 0, m_particle_count);
    }
}



// SIMD implementation

inline __m128 select(__m128 v1, __m128 v2, __m128 control)
{
    __m128 t1 = _mm_andnot_ps(control, v1);
    __m128 t2 = _mm_and_ps(v2, control);
    return _mm_or_ps(t1, t2);
}

inline __m128 dot(__m128 v1, __m128 v2)
{
    __m128 d = _mm_mul_ps(v1, v2);
    __m128 t = _mm_shuffle_ps(d, d, _MM_SHUFFLE(2, 1, 2, 1));
    d = _mm_add_ss(d, t);
    t = _mm_shuffle_ps(t, t, _MM_SHUFFLE(1, 1, 1, 1));
    d = _mm_add_ss(d, t);
    return _mm_shuffle_ps(d, d, _MM_SHUFFLE(0, 0, 0, 0));
}

inline __m128 length(__m128 v)
{
    __m128 lensq = _mm_mul_ps(v, v);
    __m128 t = _mm_shuffle_ps(lensq, lensq, _MM_SHUFFLE(1, 2, 1, 2));
    lensq = _mm_add_ss(lensq, t);
    t = _mm_shuffle_ps(t, t, _MM_SHUFFLE(1, 1, 1, 1));
    lensq = _mm_add_ss(lensq, t);
    lensq = _mm_shuffle_ps(lensq, lensq, _MM_SHUFFLE(0, 0, 0, 0));
    lensq = _mm_sqrt_ps(lensq);
    return lensq;
}

void peContext::updateVelocity_SIMD(float dt_, int begin, int end)
{
    __m128 particle_size2 = _mm_set_ps1(m_particle_size*2.0f);
    __m128 rcp_particle_size2 = _mm_set_ps1(1.0f / (m_particle_size*2.0f));
    __m128 pressure_stiffness = _mm_set_ps1(m_pressure_stiffness);
    __m128 dt = _mm_set_ps1(dt_);
    __m128 zero = _mm_setzero_ps();
    __m128 negone = _mm_set_ps1(-1.0f);

    // パーティクル同士の押し返し
    for (int i = begin; i < end; ++i) {
        __m128 pos1 = _mm_load_ps(m_particles[i].position.v);
        __m128 accel = zero;
        for (int j = 0; j < m_particle_count; ++j) {
            __m128 pos2 = _mm_load_ps(m_particles[j].position.v);
            __m128 diff = _mm_sub_ps(pos2, pos1);
            __m128 dir = _mm_mul_ps(diff, rcp_particle_size2);
            __m128 dist = length(diff);
            __m128 overlap = _mm_min_ps(zero, _mm_mul_ps(_mm_sub_ps(dist, particle_size2), pressure_stiffness));
            __m128 a = _mm_mul_ps(dir, overlap);
            accel = _mm_add_ps(accel, select(zero, a, _mm_cmpgt_ps(dist, zero)));
        }

        __m128 vel = _mm_load_ps(m_particles[i].velocity.v);
        vel = _mm_add_ps(vel, _mm_mul_ps(accel, dt));
        _mm_store_ps(m_particles[i].velocity.v, vel);
    }

    // 床との衝突
    const __m128 floor_normal = _mm_set_ps(0.0f, 0.0f, 1.0f, 0.0f);
    const __m128 floor_distance = _mm_set1_ps(-m_particle_size);
    for (int i = begin; i < end; ++i) {
        __m128 pos = _mm_load_ps(m_particles[i].position.v);
        __m128 d = _mm_add_ps(dot(pos, floor_normal), floor_distance);

        __m128 vel = _mm_load_ps(m_particles[i].velocity.v);
        __m128 accel = _mm_mul_ps(floor_normal, _mm_mul_ps(_mm_mul_ps(_mm_min_ps(zero, d), negone), pressure_stiffness));
        vel = _mm_add_ps(vel, _mm_mul_ps(accel, dt));
        _mm_store_ps(m_particles[i].velocity.v, vel);
    }

    // 重力加速
    const __m128 gravity_direction = _mm_set_ps(0.0f, 0.0f, -1.0f, 0.0f);
    const __m128 gravity_strength = _mm_set1_ps(5.0f);
    for (int i = begin; i < end; ++i) {
        __m128 vel = _mm_load_ps(m_particles[i].velocity.v);
        __m128 accel = _mm_mul_ps(gravity_direction, gravity_strength);
        vel = _mm_add_ps(vel, _mm_mul_ps(accel, dt));
        _mm_store_ps(m_particles[i].velocity.v, vel);
    }
}

void peContext::updatePosition_SIMD(float dt_, int begin, int end)
{
    __m128 dt = _mm_set_ps1(dt_);
    for (int i = begin; i < end; ++i) {
        __m128 pos = _mm_load_ps(m_particles[i].position.v);
        __m128 vel = _mm_load_ps(m_particles[i].velocity.v);
        pos = _mm_add_ps(pos, _mm_mul_ps(vel, dt));
        _mm_store_ps(m_particles[i].position.v, pos);
    }
}

void peContext::update_SIMD(float dt)
{
    if (m_multi_threading) {
        parallel_for(0, m_particle_count, m_task_granularity, [&](int begin){
            int end = std::min<int>(m_particle_count, begin + m_task_granularity);
            updateVelocity_SIMD(dt, begin, end);
        });
        parallel_for(0, m_particle_count, m_task_granularity, [&](int begin){
            int end = std::min<int>(m_particle_count, begin + m_task_granularity);
            updatePosition_SIMD(dt, begin, end);
        });
    }
    else {
        updateVelocity_SIMD(dt, 0, m_particle_count);
        updatePosition_SIMD(dt, 0, m_particle_count);
    }
}



// SIMD SoA implementation

void peSoAnize(peParticleSoA &dst, const peParticle *src, int begin, int end)
{
    for (size_t i = begin; i < end; i += 4)
    {
        soa43 tpos = soa_transpose34(
            (simdfloat4&)src[i + 0].position,
            (simdfloat4&)src[i + 1].position,
            (simdfloat4&)src[i + 2].position,
            (simdfloat4&)src[i + 3].position);
        ((simdfloat4&)dst.pos_x[i]) = tpos.x;
        ((simdfloat4&)dst.pos_y[i]) = tpos.y;
        ((simdfloat4&)dst.pos_z[i]) = tpos.z;

        soa43 tvel = soa_transpose34(
            (simdfloat4&)src[i + 0].velocity,
            (simdfloat4&)src[i + 1].velocity,
            (simdfloat4&)src[i + 2].velocity,
            (simdfloat4&)src[i + 3].velocity);
        ((simdfloat4&)dst.vel_x[i]) = tvel.x;
        ((simdfloat4&)dst.vel_y[i]) = tvel.y;
        ((simdfloat4&)dst.vel_z[i]) = tvel.z;
    }
}

void peAoSnize(peParticle *dst, const peParticleSoA &src, int begin, int end)
{
    for (size_t i = begin; i < end; i += 4)
    {
        soa44 tpos = soa_transpose44(
            ((simdfloat4&)src.pos_x[i]),
            ((simdfloat4&)src.pos_y[i]),
            ((simdfloat4&)src.pos_z[i]));
        dst[i + 0].position = (float4&)tpos.v[0];
        dst[i + 1].position = (float4&)tpos.v[1];
        dst[i + 2].position = (float4&)tpos.v[2];
        dst[i + 3].position = (float4&)tpos.v[3];

        soa44 tvel = soa_transpose44(
            ((simdfloat4&)src.vel_x[i]),
            ((simdfloat4&)src.vel_y[i]),
            ((simdfloat4&)src.vel_z[i]));
        dst[i + 0].velocity = (float4&)tvel.v[0];
        dst[i + 1].velocity = (float4&)tvel.v[1];
        dst[i + 2].velocity = (float4&)tvel.v[2];
        dst[i + 3].velocity = (float4&)tvel.v[3];
    }
}



inline __m128 soa_dot(__m128 x1, __m128 y1, __m128 z1, __m128 x2, __m128 y2, __m128 z2)
{
    __m128 x = _mm_mul_ps(x1, x2);
    __m128 y = _mm_mul_ps(y1, y2);
    __m128 z = _mm_mul_ps(z1, z2);
    return _mm_add_ps(_mm_add_ps(x, y), z);
}

inline __m128 soa_length(__m128 x, __m128 y, __m128 z)
{
    return _mm_sqrt_ps(soa_dot(x,y,z, x,y,z));
}

inline __m128 reduce_add(__m128 v1)
{
    __m128 v2 = _mm_movehl_ps(v1, v1);                      // z,w,z,w
    v1 = _mm_add_ps(v1, v2);                                // xz,yw,zz,ww
    v2 = _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(2, 2, 2, 2));   // yw, yw, yw, yw
    v1 = _mm_add_ss(v1, v2);                                // xzyw, ...
    return v1;
}

void peContext::updateVelocity_SIMDSoA(float dt_, int begin, int end)
{
    __m128 particle_size2 = _mm_set_ps1(m_particle_size*2.0f);
    __m128 rcp_particle_size2 = _mm_set_ps1(1.0f / (m_particle_size*2.0f));
    __m128 pressure_stiffness = _mm_set_ps1(m_pressure_stiffness);
    __m128 dt = _mm_set_ps1(dt_);
    __m128 zero = _mm_setzero_ps();
    __m128 negone = _mm_set_ps1(-1.0f);

    // パーティクル同士の押し返し
    for (int i = begin; i < end; ++i) {
        __m128 pos1x = _mm_set_ps1(m_soa.pos_x[i]);
        __m128 pos1y = _mm_set_ps1(m_soa.pos_y[i]);
        __m128 pos1z = _mm_set_ps1(m_soa.pos_z[i]);
        __m128 accelx = zero;
        __m128 accely = zero;
        __m128 accelz = zero;

        for (int j = 0; j < m_particle_count; j += 4) {
            __m128 pos2x = _mm_load_ps(&m_soa.pos_x[j]);
            __m128 pos2y = _mm_load_ps(&m_soa.pos_y[j]);
            __m128 pos2z = _mm_load_ps(&m_soa.pos_z[j]);

            __m128 diffx = _mm_sub_ps(pos2x, pos1x);
            __m128 diffy = _mm_sub_ps(pos2y, pos1y);
            __m128 diffz = _mm_sub_ps(pos2z, pos1z);

            __m128 dirx = _mm_mul_ps(diffx, rcp_particle_size2);
            __m128 diry = _mm_mul_ps(diffy, rcp_particle_size2);
            __m128 dirz = _mm_mul_ps(diffz, rcp_particle_size2);

            __m128 dist = soa_length(diffx, diffy, diffz);
            __m128 overlap = _mm_min_ps(zero, _mm_mul_ps(_mm_sub_ps(dist, particle_size2), pressure_stiffness));

            __m128 ax = _mm_mul_ps(dirx, overlap);
            __m128 ay = _mm_mul_ps(diry, overlap);
            __m128 az = _mm_mul_ps(dirz, overlap);
            __m128 gt = _mm_cmpgt_ps(dist, zero);
            accelx = _mm_add_ps(accelx, select(zero, ax, gt));
            accely = _mm_add_ps(accely, select(zero, ay, gt));
            accelz = _mm_add_ps(accelz, select(zero, az, gt));
        }

        __m128 velx = _mm_set_ss(m_soa.vel_x[i]);
        __m128 vely = _mm_set_ss(m_soa.vel_y[i]);
        __m128 velz = _mm_set_ss(m_soa.vel_z[i]);
        velx = _mm_add_ss(velx, _mm_mul_ss(reduce_add(accelx), dt));
        vely = _mm_add_ss(vely, _mm_mul_ss(reduce_add(accely), dt));
        velz = _mm_add_ss(velz, _mm_mul_ss(reduce_add(accelz), dt));
        _mm_store_ss(&m_soa.vel_x[i], velx);
        _mm_store_ss(&m_soa.vel_y[i], vely);
        _mm_store_ss(&m_soa.vel_z[i], velz);
    }

    // 床との衝突
    const __m128 floor_normal_x = _mm_set1_ps(0.0f);
    const __m128 floor_normal_y = _mm_set1_ps(1.0f);
    const __m128 floor_normal_z = _mm_set1_ps(0.0f);
    const __m128 floor_distance = _mm_set1_ps(-m_particle_size);
    for (int i = begin; i < end; i += 4) {
        __m128 posx = _mm_load_ps(&m_soa.pos_x[i]);
        __m128 posy = _mm_load_ps(&m_soa.pos_y[i]);
        __m128 posz = _mm_load_ps(&m_soa.pos_z[i]);
        __m128 d = _mm_add_ps(soa_dot(posx, posy, posz, floor_normal_x, floor_normal_y, floor_normal_z), floor_distance);
        __m128 a = _mm_mul_ps(_mm_mul_ps(_mm_min_ps(zero, d), negone), pressure_stiffness);

        __m128 velx = _mm_load_ps(&m_soa.vel_x[i]);
        __m128 vely = _mm_load_ps(&m_soa.vel_y[i]);
        __m128 velz = _mm_load_ps(&m_soa.vel_z[i]);

        __m128 accelx = _mm_mul_ps(floor_normal_x, a);
        __m128 accely = _mm_mul_ps(floor_normal_y, a);
        __m128 accelz = _mm_mul_ps(floor_normal_z, a);
        velx = _mm_add_ps(velx, _mm_mul_ps(accelx, dt));
        vely = _mm_add_ps(vely, _mm_mul_ps(accely, dt));
        velz = _mm_add_ps(velz, _mm_mul_ps(accelz, dt));

        _mm_store_ps(&m_soa.vel_x[i], velx);
        _mm_store_ps(&m_soa.vel_y[i], vely);
        _mm_store_ps(&m_soa.vel_z[i], velz);
    }

    // 重力加速
    const __m128 gravity_direction_x = _mm_set1_ps(0.0f);
    const __m128 gravity_direction_y = _mm_set1_ps(-1.0f);
    const __m128 gravity_direction_z = _mm_set1_ps(0.0f);
    const __m128 gravity_strength = _mm_set1_ps(5.0f);
    for (int i = begin; i < end; i += 4) {
        __m128 velx = _mm_load_ps(&m_soa.vel_x[i]);
        __m128 vely = _mm_load_ps(&m_soa.vel_y[i]);
        __m128 velz = _mm_load_ps(&m_soa.vel_z[i]);
        velx = _mm_add_ps(velx, _mm_mul_ps(_mm_mul_ps(gravity_direction_x, gravity_strength), dt));
        vely = _mm_add_ps(vely, _mm_mul_ps(_mm_mul_ps(gravity_direction_y, gravity_strength), dt));
        velz = _mm_add_ps(velz, _mm_mul_ps(_mm_mul_ps(gravity_direction_z, gravity_strength), dt));
        _mm_store_ps(&m_soa.vel_x[i], velx);
        _mm_store_ps(&m_soa.vel_y[i], vely);
        _mm_store_ps(&m_soa.vel_z[i], velz);
    }
}

void peContext::updatePosition_SIMDSoA(float dt_, int begin, int end)
{
    __m128 dt = _mm_set_ps1(dt_);
    for (int i = begin; i < end; i += 4) {
        __m128 posx = _mm_load_ps(&m_soa.pos_x[i]);
        __m128 posy = _mm_load_ps(&m_soa.pos_y[i]);
        __m128 posz = _mm_load_ps(&m_soa.pos_z[i]);

        __m128 velx = _mm_load_ps(&m_soa.vel_x[i]);
        __m128 vely = _mm_load_ps(&m_soa.vel_y[i]);
        __m128 velz = _mm_load_ps(&m_soa.vel_z[i]);

        posx = _mm_add_ps(posx, _mm_mul_ps(velx, dt));
        posy = _mm_add_ps(posy, _mm_mul_ps(vely, dt));
        posz = _mm_add_ps(posz, _mm_mul_ps(velz, dt));

        _mm_store_ps(&m_soa.pos_x[i], posx);
        _mm_store_ps(&m_soa.pos_y[i], posy);
        _mm_store_ps(&m_soa.pos_z[i], posz);
    }
}

void peContext::update_SIMDSoA(float dt)
{
    if (m_multi_threading)
    {
        parallel_for(0, m_particle_count, m_task_granularity, [&](int begin){
            int end = std::min<int>(m_particle_count, begin + m_task_granularity);
            peSoAnize(m_soa, m_particles, begin, end);
        });
        parallel_for(0, m_particle_count, m_task_granularity, [&](int begin){
            int end = std::min<int>(m_particle_count, begin + m_task_granularity);
            updateVelocity_SIMDSoA(dt, begin, end);
        });
        parallel_for(0, m_particle_count, m_task_granularity, [&](int begin){
            int end = std::min<int>(m_particle_count, begin + m_task_granularity);
            updatePosition_SIMDSoA(dt, begin, end);
        });
        parallel_for(0, m_particle_count, m_task_granularity, [&](int begin){
            int end = std::min<int>(m_particle_count, begin + m_task_granularity);
            peAoSnize(m_particles, m_soa, begin, end);
        });
    }
    else
    {
        peSoAnize(m_soa, m_particles, 0, m_particle_count);
        updateVelocity_SIMDSoA(dt, 0, m_particle_count);
        updatePosition_SIMDSoA(dt, 0, m_particle_count);
        peAoSnize(m_particles, m_soa, 0, m_particle_count);
    }
}



//  ISPC implementation

void peContext::update_ISPC(float dt)
{
    ispc::Context ic;
    ic.pos_x = m_soa.pos_x;
    ic.pos_y = m_soa.pos_y;
    ic.pos_z = m_soa.pos_z;
    ic.vel_x = m_soa.vel_x;
    ic.vel_y = m_soa.vel_y;
    ic.vel_z = m_soa.vel_z;
    ic.particle_count = m_particle_count;
    ic.pressure_stiffness = m_pressure_stiffness;
    ic.particle_size = m_particle_size;
    ic.rcp_particle_size2 = 1.0f / (m_particle_size * 2.0f);
    ic.timestep = dt;


    if (m_multi_threading)
    {
        parallel_for(0, m_particle_count, m_task_granularity, [&](int begin){
            int end = std::min<int>(m_particle_count, begin + m_task_granularity);
            peSoAnize(m_soa, m_particles, begin, end);
        });
        parallel_for(0, m_particle_count, m_task_granularity, [&](int begin){
            int end = std::min<int>(m_particle_count, begin + m_task_granularity);
            ispc::UpdateVelocity(ic, begin, end);
        });
        parallel_for(0, m_particle_count, m_task_granularity, [&](int begin){
            int end = std::min<int>(m_particle_count, begin + m_task_granularity);
            ispc::UpdatePosition(ic, begin, end);
        });
        parallel_for(0, m_particle_count, m_task_granularity, [&](int begin){
            int end = std::min<int>(m_particle_count, begin + m_task_granularity);
            peAoSnize(m_particles, m_soa, begin, end);
        });
    }
    else
    {
        peSoAnize(m_soa, m_particles, 0, m_particle_count);
        ispc::UpdateVelocity(ic, 0, m_particle_count);
        ispc::UpdatePosition(ic, 0, m_particle_count);
        peAoSnize(m_particles, m_soa, 0, m_particle_count);
    }
}



peCLinkage peExport peContext* peCreateContext(int n) { return new peContext(n); }
peCLinkage peExport void peDestroyContext(peContext *ctx) { delete ctx; }

peCLinkage peExport void peSetParticleSize(peContext *ctx, float v)
{
    ctx->setParticleSize(v);
}

peCLinkage peExport void peSetUpdateRoutine(peContext *ctx, peUpdateRoutine v)
{
    ctx->setUpdateRoutine(v);
}

peCLinkage peExport void peEnableMultiThreading(peContext *ctx, bool v)
{
    ctx->enbaleMultiThreading(v);
}

peCLinkage peExport void peUpdate(peContext *ctx, float dt)
{
    ctx->update(dt);
}

peCLinkage peExport void peCopyDataToTexture(peContext *ctx, void *texture, int width, int height)
{
    if (ctx) {
        ctx->copyDataToTexture(texture, width, height);
    }
}





// Graphics device identifiers in Unity
enum GfxDeviceRenderer
{
    kGfxRendererOpenGL = 0,          // OpenGL
    kGfxRendererD3D9,                // Direct3D 9
    kGfxRendererD3D11,               // Direct3D 11
    kGfxRendererGCM,                 // Sony PlayStation 3 GCM
    kGfxRendererNull,                // "null" device (used in batch mode)
    kGfxRendererHollywood,           // Nintendo Wii
    kGfxRendererXenon,               // Xbox 360
    kGfxRendererOpenGLES,            // OpenGL ES 1.1
    kGfxRendererOpenGLES20Mobile,    // OpenGL ES 2.0 mobile variant
    kGfxRendererMolehill,            // Flash 11 Stage3D
    kGfxRendererOpenGLES20Desktop,   // OpenGL ES 2.0 desktop variant (i.e. NaCl)
    kGfxRendererCount
};

// Event types for UnitySetGraphicsDevice
enum GfxDeviceEventType {
    kGfxDeviceEventInitialize = 0,
    kGfxDeviceEventShutdown,
    kGfxDeviceEventBeforeReset,
    kGfxDeviceEventAfterReset,
};

peCopyToTextureBase* peCreateCopyToTextureD3D9(void *device);
peCopyToTextureBase* peCreateCopyToTextureD3D11(void *device);
peCopyToTextureBase* peCreateCopyToTextureOpenGL(void *device);

#define peSupportD3D11

peCLinkage peExport void UnitySetGraphicsDevice(void* device, int deviceType, int eventType)
{
    if (eventType == kGfxDeviceEventInitialize) {
#ifdef peSupportD3D9
        if (deviceType == kGfxRendererD3D9)
        {
            // todo
        }
#endif // peSupportD3D9
#ifdef peSupportD3D11
        if (deviceType == kGfxRendererD3D11)
        {
            g_CopyToTexture = peCreateCopyToTextureD3D11(device);
        }
#endif // peSupportD3D11
#ifdef peSupportOpenGL
        if (deviceType == kGfxRendererOpenGL)
        {
            g_CopyToTexture = peCreateCopyToTextureOpenGL(device);
        }
#endif // peSupportOpenGL
    }

    if (eventType == kGfxDeviceEventShutdown) {
        delete g_CopyToTexture;
        g_CopyToTexture = nullptr;
    }
}

peCLinkage peExport void UnityRenderEvent(int eventID)
{
}
