#include <cstdlib>
#include <algorithm>
#include <random>
#include "ParticleEngine.h"
#include "ParticleCore_ispc.h"

#include "DirectXMath.h"
using namespace DirectX;


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


class peContext
{
public:
    peContext(int n);
    ~peContext();

    void update_Plain(float dt);
    void update_SIMD(float dt);
    void update_SIMDSoA(float dt);
    void update_ISPC(float dt);

private:
    peParticle *m_particles;
    peParticleSoA m_soa;
    int m_particle_count;

    float m_particle_size;
    float m_pressure_stiffness;
};


peContext::peContext(int n)
    : m_particle_count(n)
    , m_particle_size(0.1f)
    , m_pressure_stiffness(1000.0f)
{
    m_particles = (peParticle*)peAlignedAlloc(sizeof(peParticle)*n);
    m_soa.pos_x = (float*)peAlignedAlloc(sizeof(float)*n);
    m_soa.pos_y = (float*)peAlignedAlloc(sizeof(float)*n);
    m_soa.pos_z = (float*)peAlignedAlloc(sizeof(float)*n);
    m_soa.vel_x = (float*)peAlignedAlloc(sizeof(float)*n);
    m_soa.vel_y = (float*)peAlignedAlloc(sizeof(float)*n);
    m_soa.vel_z = (float*)peAlignedAlloc(sizeof(float)*n);

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

void peContext::update_Plain(float dt)
{
    const float rcp_particle_size2 = 1.0f / (m_particle_size*2.0f);

    // パーティクル同士の押し返し
    for (int i = 0; i < m_particle_count; ++i) {
        float3 &pos1 = (float3&)m_particles[i].position;
        float3 accel = { 0.0f, 0.0f, 0.0f };
        for (int j = 0; j < m_particle_count; ++j) {
            float3 &pos2 = (float3&)m_particles[j].position;
            float3 diff = pos2 - pos1;
            float3 dir = diff * rcp_particle_size2;
            float dist = length(diff);
            float3 a = dir * (std::min<float>(0.0f, dist - (m_particle_size*2.0f)) * m_pressure_stiffness);
            accel = accel + (a * (dist > 0.0f, 1.0f, 0.0f));
        }

        float3 &vel = (float3&)m_particles[i].velocity;
        vel = vel + accel * dt;
    }

    // 床との衝突
    const float3 floor_normal = { 0.0f, 1.0f, 0.0f };
    const float floor_distance = 5.0f;
    for (int i = 0; i < m_particle_count; ++i) {
        float3 &pos = (float3&)m_particles[i].position;
        float d = dot(pos, floor_normal) + floor_distance;

        float3 &vel = (float3&)m_particles[i].velocity;
        vel = vel + (floor_normal * std::min<float>(0.0f, -d * m_pressure_stiffness) * dt);
    }

    // 重力加速
    const float3 gravity_direction = { 0.0f, -1.0f, 0.0f };
    const float gravity_strength = 5.0f;
    for (int i = 0; i < m_particle_count; ++i) {
        float3 &vel = (float3&)m_particles[i].velocity;
        vel = vel + (gravity_direction * gravity_strength * dt);
    }

    for (int i = 0; i < m_particle_count; ++i) {
        float3 &pos = (float3&)m_particles[i].position;
        float3 &vel = (float3&)m_particles[i].velocity;
        pos = pos + (vel * dt);
    }
}

void peContext::update_SIMD(float dt_)
{
    XMVECTOR particle_size2 = XMVectorReplicate(m_particle_size*2.0f);
    XMVECTOR rcp_particle_size2 = XMVectorReplicate(1.0f / (m_particle_size*2.0f));
    XMVECTOR pressure_stiffness = XMVectorReplicate(m_pressure_stiffness);
    XMVECTOR dt = XMVectorReplicate(dt_);
    XMVECTOR zero = XMVectorZero();

    for (int i = 0; i < m_particle_count; ++i) {
        XMVECTOR pos1 = XMLoadFloat4A((XMFLOAT4A*)&m_particles[i].position);
        XMVECTOR accel = zero;
        for (int j = 0; j < m_particle_count; ++j) {
            XMVECTOR pos2 = XMLoadFloat4A((XMFLOAT4A*)&m_particles[j].position);
            XMVECTOR diff = XMVectorSubtract(pos2, pos1);
            XMVECTOR dir = XMVectorMultiply(diff, rcp_particle_size2);
            XMVECTOR dist = XMVector3Length(diff);
            XMVECTOR a = XMVectorMultiply(dir,
                XMVectorMin(zero, XMVectorMultiply(XMVectorSubtract(dist, particle_size2), pressure_stiffness)));
            accel = XMVectorAdd(accel, XMVectorSelect(a, zero, XMVectorGreater(dist, zero)));
        }

        XMVECTOR vel = XMLoadFloat4A((XMFLOAT4A*)&m_particles[i].velocity);
        vel = XMVectorAdd(vel, XMVectorMultiply(accel, dt));
        XMStoreFloat4A((XMFLOAT4A*)&m_particles[i].velocity, vel);
    }
}




void peSoAnize(peParticleSoA &dst, const peParticle *src, size_t n)
{
    for (size_t i = 0; i < n; i += 4)
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

void peAoSnize(peParticle *dst, const peParticleSoA &src, size_t n)
{
    for (size_t i = 0; i < n; i += 4)
    {
        soa44 tpos = soa_transpose44(
            ((simdfloat4&)src.pos_x[i]),
            ((simdfloat4&)src.pos_x[i]),
            ((simdfloat4&)src.pos_x[i]));
        dst[i + 0].position = (float4&)tpos.v[0];
        dst[i + 1].position = (float4&)tpos.v[1];
        dst[i + 2].position = (float4&)tpos.v[2];
        dst[i + 3].position = (float4&)tpos.v[3];

        soa44 tvel = soa_transpose44(
            ((simdfloat4&)src.vel_x[i]),
            ((simdfloat4&)src.vel_x[i]),
            ((simdfloat4&)src.vel_x[i]));
        dst[i + 0].velocity = (float4&)tpos.v[0];
        dst[i + 1].velocity = (float4&)tpos.v[1];
        dst[i + 2].velocity = (float4&)tpos.v[2];
        dst[i + 3].velocity = (float4&)tpos.v[3];
    }
}


inline XMVECTOR XMVectorSoALength(XMVECTOR x, XMVECTOR y, XMVECTOR z)
{
    XMVECTOR x2 = XMVectorMultiply(x, x);
    XMVECTOR y2 = XMVectorMultiply(y, y);
    XMVECTOR z2 = XMVectorMultiply(z, z);
    return XMVectorSqrt(XMVectorAdd(XMVectorAdd(x2, y2), z2));
}

inline XMVECTOR XMVectorReduceAdd(XMVECTOR v1)
{
    XMVECTOR v2 = _mm_movehl_ps(v1, v1);                    // v2: z,w,z,w
    v1 = _mm_add_ps(v1, v2);                                // v1: xz,yw,zz,ww
    v2 = _mm_shuffle_ps(v1, v1, _MM_SHUFFLE(2, 2, 2, 2));   // v2: yw, yw, yw, yw
    v1 = _mm_add_ss(v1, v2);                                // v1: xzyw, ...
    return v1;
}

void peContext::update_SIMDSoA(float dt_)
{
    XMVECTOR particle_size2 = XMVectorReplicate(m_particle_size*2.0f);
    XMVECTOR rcp_particle_size2 = XMVectorReplicate(1.0f / (m_particle_size*2.0f));
    XMVECTOR pressure_stiffness = XMVectorReplicate(m_pressure_stiffness);
    XMVECTOR dt = XMVectorReplicate(dt_);
    XMVECTOR zero = XMVectorZero();

    peSoAnize(m_soa, m_particles, m_particle_count);
    for (int i = 0; i < m_particle_count; ++i) {
        XMVECTOR pos1x = XMVectorReplicate(m_soa.pos_x[i]);
        XMVECTOR pos1y = XMVectorReplicate(m_soa.pos_y[i]);
        XMVECTOR pos1z = XMVectorReplicate(m_soa.pos_z[i]);
        XMVECTOR accelx = zero;
        XMVECTOR accely = zero;
        XMVECTOR accelz = zero;

        for (int j = 0; j < m_particle_count; j += 4) {
            XMVECTOR pos2x = XMLoadFloat4A((XMFLOAT4A*)&m_soa.pos_x[j]);
            XMVECTOR pos2y = XMLoadFloat4A((XMFLOAT4A*)&m_soa.pos_y[j]);
            XMVECTOR pos2z = XMLoadFloat4A((XMFLOAT4A*)&m_soa.pos_z[j]);

            XMVECTOR diffx = XMVectorSubtract(pos2x, pos1x);
            XMVECTOR diffy = XMVectorSubtract(pos2y, pos1y);
            XMVECTOR diffz = XMVectorSubtract(pos2z, pos1z);

            XMVECTOR dirx = XMVectorMultiply(diffx, rcp_particle_size2);
            XMVECTOR diry = XMVectorMultiply(diffy, rcp_particle_size2);
            XMVECTOR dirz = XMVectorMultiply(diffz, rcp_particle_size2);

            XMVECTOR dist = XMVectorSoALength(diffx, diffy, diffz);
            XMVECTOR penetration = XMVectorMin(zero, XMVectorMultiply(XMVectorSubtract(dist, particle_size2), pressure_stiffness));

            XMVECTOR ax = XMVectorMultiply(dirx, penetration);
            XMVECTOR ay = XMVectorMultiply(diry, penetration);
            XMVECTOR az = XMVectorMultiply(dirz, penetration);
            accelx = XMVectorAdd(accelx, XMVectorSelect(ax, zero, XMVectorGreater(dist, zero)));
            accely = XMVectorAdd(accely, XMVectorSelect(ax, zero, XMVectorGreater(dist, zero)));
            accelz = XMVectorAdd(accelz, XMVectorSelect(ax, zero, XMVectorGreater(dist, zero)));
        }

        XMVECTOR velx = XMVectorReplicate(m_soa.vel_x[i]);
        XMVECTOR vely = XMVectorReplicate(m_soa.vel_y[i]);
        XMVECTOR velz = XMVectorReplicate(m_soa.vel_z[i]);
        velx = XMVectorAdd(velx, XMVectorMultiply(XMVectorReduceAdd(accelx), dt));
        vely = XMVectorAdd(vely, XMVectorMultiply(XMVectorReduceAdd(accely), dt));
        velz = XMVectorAdd(velz, XMVectorMultiply(XMVectorReduceAdd(accelz), dt));
        XMStoreFloat(&m_soa.vel_x[i], velx);
        XMStoreFloat(&m_soa.vel_y[i], vely);
        XMStoreFloat(&m_soa.vel_z[i], velz);
    }
    peAoSnize(m_particles, m_soa, m_particle_count);
}

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

    peSoAnize(m_soa, m_particles, m_particle_count);
    ispc::UpdatePressure(ic);
    ispc::Integrate(ic);
    peAoSnize(m_particles, m_soa, m_particle_count);
}



peCLinkage peExport peContext* peCreateContext(int n) { return new peContext(n); }
peCLinkage peExport void peDestroyContext(peContext *ctx) { delete ctx; }

peCLinkage peExport void peUpdate_Plain(peContext *ctx, float dt) { ctx->update_Plain(dt); }
peCLinkage peExport void peUpdate_SIMD(peContext *ctx, float dt) { ctx->update_SIMD(dt); }
peCLinkage peExport void peUpdate_SIMDSoA(peContext *ctx, float dt) { ctx->update_SIMDSoA(dt); }
peCLinkage peExport void peUpdate_ISPC(peContext *ctx, float dt) { ctx->update_ISPC(dt); }
