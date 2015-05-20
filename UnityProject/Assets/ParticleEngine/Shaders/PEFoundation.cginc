#ifndef PEFoundation_h
#define PEFoundation_h

#include "../../BatchRenderer/Shaders/Math.cginc"

int         g_batch_begin;
int         g_num_instances;
float       g_size;
float       g_fade_time;
float       g_spin;
float4      g_instance_data_size;
#if defined(SHADER_API_D3D11)
Texture2D<float4> g_instance_data;
#endif


float3 iq_rand( float3 p )
{
    p = float3( dot(p,float3(127.1,311.7,311.7)), dot(p,float3(269.5,183.3,183.3)), dot(p,float3(269.5,183.3,183.3)) );
    return frac(sin(p)*43758.5453)*2.0-1.0;
}


void GetParticleParams(int iid, out float4 o_pos, out float4 o_vel)
{
    float i = iid * 2;
    float4 t = float4(
        g_instance_data_size.xy * float2(fmod(i, g_instance_data_size.z) + 0.5, floor(i/g_instance_data_size.z) + 0.5),
        0.0, 0.0);
    float4 pitch = float4(g_instance_data_size.x, 0.0, 0.0, 0.0);
#if defined(SHADER_API_D3D11)
    uint x = iid % 128;
    uint y = iid / 128;
    o_pos   = g_instance_data[uint2(x*2+0, y)];
    o_vel   = g_instance_data[uint2(x*2+1, y)];
#endif
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
