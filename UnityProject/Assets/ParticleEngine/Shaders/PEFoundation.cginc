#ifndef PEFoundation_h
#define PEFoundation_h

#include "../../BatchRenderer/Shaders/Math.cginc"

int         g_batch_begin;
int         g_num_instances;
float       g_size;
float       g_fade_time;
float       g_spin;
float4      g_instance_data_size;
sampler2D   g_instance_data;
int         g_use_buffer;

#if defined(SHADER_API_D3D11)
struct Particle
{
    float4 position;
    float4 velocity;
};
StructuredBuffer<Particle> g_instance_buffer;
#endif



void GetParticleParams(int iid, out float4 o_pos, out float4 o_vel)
{
#if defined(SHADER_API_D3D11)
    if(g_use_buffer) {
        o_pos   = g_instance_buffer[iid].position;
        o_vel   = g_instance_buffer[iid].velocity;
    }
    else
#endif
    {
        float i = iid * 2;
        float4 t = float4(
            g_instance_data_size.xy * float2(fmod(i, g_instance_data_size.z) + 0.5, floor(i/g_instance_data_size.z) + 0.5),
            0.0, 0.0);
        float4 pitch = float4(g_instance_data_size.x, 0.0, 0.0, 0.0);
        o_pos = tex2Dlod(g_instance_data, t + pitch*0.0);
        o_vel = tex2Dlod(g_instance_data, t + pitch*1.0);
    }
}

void ParticleTransform(inout appdata_full v, out float4 pos, out float4 vel)
{
    int iid = v.texcoord1.x + g_batch_begin;
    GetParticleParams(iid, pos, vel);

    v.vertex.xyz *= g_size;
    v.vertex.xyz += pos.xyz;
    if(iid >= g_num_instances) {
        v.vertex.xyz = 0.0;
    }
}


#endif // PEFoundation_h
