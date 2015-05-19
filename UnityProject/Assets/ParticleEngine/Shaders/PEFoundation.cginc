#ifndef PEFoundation_h
#define PEFoundation_h

#include "../../BatchRenderer/Shaders/Math.cginc"

int         g_batch_begin;
sampler2D   g_instance_data;
float       g_size;
float       g_fade_time;
float       g_spin;
float4      g_instance_data_size;


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
    o_pos   = tex2Dlod(g_instance_data, t + pitch*0.0);
    o_vel   = tex2Dlod(g_instance_data, t + pitch*1.0);
}

void ParticleTransform(inout appdata_full v, out float4 pos, out float4 vel)
{
    int iid = v.texcoord1.x + g_batch_begin;
    GetParticleParams(iid, pos, vel);

    v.vertex.xyz *= g_size;
    v.vertex.xyz += pos.xyz;
}


#endif // PEFoundation_h
