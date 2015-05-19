#ifndef PESurface_h
#define PESurface_h

#define PE_ENABLE_SPIN
#define PE_ENABLE_HEAT_EMISSION

#ifdef PE_SHADOW_COLLECTOR
#   define SHADOW_COLLECTOR_PASS
#endif // PE_SHADOW_COLLECTOR

#include "UnityCG.cginc"
#include "PEFoundation.cginc"


#ifdef PE_DEPTH_PREPASS
    #define PE_SHADOW_CASTER
#endif // PE_DEPTH_PREPASS



#if defined(PE_STANDARD) || defined(PE_SURFACE)
    sampler2D _MainTex;
    fixed4 _Color;

    struct Input {
        float2 uv_MainTex;
    #ifdef PE_ENABLE_HEAT_EMISSION
        float4 velocity;
    #endif // PE_ENABLE_HEAT_EMISSION
    };

    void vert(inout appdata_full v, out Input data)
    {
        UNITY_INITIALIZE_OUTPUT(Input,data);

        float4 pos;
        float4 vel;
        ParticleTransform(v, pos, vel);

    #ifdef PE_ENABLE_HEAT_EMISSION
        data.velocity = vel;
    #endif // PE_ENABLE_HEAT_EMISSION
    }
#endif // defined(PE_STANDARD) || defined(PE_SURFACE)



// legacy surface shader
#ifdef PE_SURFACE
    void surf(Input data, inout SurfaceOutput o)
    {
        o.Albedo = _Color * tex2D(_MainTex, data.uv_MainTex);

    #ifdef PE_ENABLE_HEAT_EMISSION
        float speed = data.velocity.w;
        float ei = max(speed-2.0, 0.0) * 1.0;
        o.Emission = float3(0.25, 0.05, 0.025)*ei;
    #endif // PE_ENABLE_HEAT_EMISSION
    }
#endif // PE_SURFACE



// standard shader
#ifdef PE_STANDARD
    half _Glossiness;
    half _Metallic;

    void surf(Input IN, inout SurfaceOutputStandard o)
    {
        fixed4 c = tex2D(_MainTex, IN.uv_MainTex) * _Color;
        o.Albedo = c.rgb;
        o.Metallic = _Metallic;
        o.Smoothness = _Glossiness;
        o.Alpha = c.a;

    #ifdef PE_ENABLE_HEAT_EMISSION
        float speed = IN.velocity.w;
        float ei = max(speed-2.0, 0.0) * 1.0;
        o.Emission += float3(0.25, 0.05, 0.025)*ei;
    #endif // PE_ENABLE_HEAT_EMISSION
    }
#endif // PE_STANDARD



// shadow caster
#ifdef PE_SHADOW_CASTER
    struct v2f
    { 
        V2F_SHADOW_CASTER;
    };

    v2f vert( appdata_full v )
    {
        float4 pos;
        float4 vel;
        ParticleTransform(v, pos, vel);

        v2f o;
        TRANSFER_SHADOW_CASTER(o)
        return o;
    }

    float4 frag( v2f i ) : SV_Target
    {
        SHADOW_CASTER_FRAGMENT(i)
    }
#endif // PE_SHADOW_CASTER



// legacy shadow collector
#ifdef PE_SHADOW_COLLECTOR
    struct v2f { 
        V2F_SHADOW_COLLECTOR;
    };

    v2f vert( appdata_full v )
    {
        float4 pos;
        float4 vel;
        ParticleTransform(v, pos, vel);

        v2f o;
        TRANSFER_SHADOW_COLLECTOR(o)
        return o;
    }
        
    fixed4 frag (v2f i) : SV_Target
    {
        SHADOW_COLLECTOR_FRAGMENT(i)
    }
#endif // PE_SHADOW_COLLECTOR


#endif // PESurface_h
