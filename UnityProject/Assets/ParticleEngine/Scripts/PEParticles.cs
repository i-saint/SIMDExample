using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;
#if UNITY_EDITOR
using UnityEditor;
#endif


public class PEParticles : MonoBehaviour
{
    [StructLayout(LayoutKind.Explicit)]
    public struct peParticle
    {
        [FieldOffset(0)]  public Vector3 position;
        [FieldOffset(0)]  public Vector4 position4;
        [FieldOffset(16)] public Vector3 velocity;
        [FieldOffset(16)] public Vector4 velocity4;
    }
    public delegate void CSUpdateRoutine(float dt, int begin, int end);

    public enum peUpdateRoutine
    {
        Plain,
        SIMD,
        SIMDSoA,
        ISPC,
        CSharp,
    }

    [DllImport ("ParticleEngine")] public static extern IntPtr  peCreateContext(int num);
    [DllImport ("ParticleEngine")] public static extern void    peDestroyContext(IntPtr ctx);

    [DllImport ("ParticleEngine")] public static extern void    peSetUpdateRoutine(IntPtr ctx, peUpdateRoutine v);
    [DllImport ("ParticleEngine")] public static extern void    peEnableMultiThreading(IntPtr ctx, bool v);
    [DllImport ("ParticleEngine")] public static extern void    peSetParticleSize(IntPtr ctx, float sv);
    [DllImport ("ParticleEngine")] public static extern void    peSetPressureStiffness(IntPtr ctx, float v);
    [DllImport ("ParticleEngine")] public static extern void    peSetWallStiffness(IntPtr ctx, float v);
    [DllImport ("ParticleEngine")] public static extern void    peSetCSUpdateRoutine(IntPtr ctx, CSUpdateRoutine vel, CSUpdateRoutine pos);
    [DllImport ("ParticleEngine")] public static extern void    peUpdate(IntPtr ctx, float dt);
    [DllImport ("ParticleEngine")] public static extern void    peCopyDataToTexture(IntPtr ctx, IntPtr texture, int width, int height);
    [DllImport ("ParticleEngine")] unsafe public static extern peParticle* peGetParticles(IntPtr ctx);

    [DllImport ("ParticleEngine")] public static extern IntPtr  peBenchmark(IntPtr ctx, int loop_count);



    public const int DataTextureWidth = 128;

    public peUpdateRoutine m_routine = peUpdateRoutine.ISPC;
    public bool m_multi_threading = true;
    public int m_particle_count = 1024*4;
    public float m_particle_size = 0.1f;
    public float m_pressure_stiffness = 500.0f;
    public float m_wall_stiffness = 1500.0f;

    IntPtr m_ctx;

    public void SetUpdateRoutibe(int v)
    {
        m_routine = (peUpdateRoutine)v;
    }

    public void SetParticleCount(int v)
    {
        m_particle_count = v;
        ResetContext();
    }

    public void EnableMultiThreading(bool v)
    {
        m_multi_threading = v;
    }

    public void CopyDataToTexture(Texture tex)
    {
        peCopyDataToTexture(m_ctx, tex.GetNativeTexturePtr(), tex.width, tex.height);
    }



    void ResetContext()
    {
        if (m_ctx != IntPtr.Zero)
        {
            peDestroyContext(m_ctx);
            m_ctx = IntPtr.Zero;
        }
        m_particle_count = Mathf.Max(DataTextureWidth, m_particle_count);
        m_ctx = peCreateContext(m_particle_count);
    }


    void OnEnable()
    {
        ResetContext();
    }

    void OnDisable()
    {
        peDestroyContext(m_ctx);
        m_ctx = IntPtr.Zero;
    }

    void Update()
    {
        peSetUpdateRoutine(m_ctx, m_routine);
        peEnableMultiThreading(m_ctx, m_multi_threading);
        peSetParticleSize(m_ctx, m_particle_size);
        peSetPressureStiffness(m_ctx, m_pressure_stiffness);
        peSetWallStiffness(m_ctx, m_wall_stiffness);
        peSetCSUpdateRoutine(m_ctx, Update_Velocity, Update_Position);

        peUpdate(m_ctx, Time.deltaTime);
    }


    unsafe void Update_Velocity(float dt, int begin, int end)
    {
        peParticle *particles = peGetParticles(m_ctx);
        float particle_size2 = m_particle_size * 2.0f;
        float rcp_particle_size2 = 1.0f / (m_particle_size * 2.0f);

        // パーティクル同士の押し返し
        for (int i = begin; i < end; ++i) {
            Vector3 pos1 = particles[i].position;
            Vector3 accel = Vector3.zero;
            for (int j = 0; j < m_particle_count; ++j) {
                Vector3 pos2 = particles[j].position;
                Vector3 diff = pos2 - pos1;
                Vector3 dir = diff * rcp_particle_size2;
                float dist = diff.magnitude;
                if (dist > 0.0f) {
                    Vector3 a = dir * (Mathf.Min(0.0f, dist - particle_size2) * m_pressure_stiffness);
                    accel = accel + a;
                }
            }

            Vector3 vel = particles[i].velocity;
            vel = vel + accel * dt;
            particles[i].velocity = vel;
        }

        // 床との衝突
        Vector3 floor_normal = new Vector3(0.0f, 1.0f, 0.0f);
        float floor_distance = -m_particle_size;
        for (int i = begin; i < end; ++i) {
            Vector3 pos = particles[i].position;
            float d = Vector3.Dot(pos, floor_normal) + floor_distance;
            Vector3 accel = floor_normal * (-Mathf.Min(0.0f, d) * m_wall_stiffness);
            Vector3 vel = particles[i].velocity;
            vel = vel + accel * dt;
            particles[i].velocity = vel;
        }

        // 重力加速
        Vector3 gravity_direction = new Vector3(0.0f, -1.0f, 0.0f);
        float gravity_strength = 5.0f;
        for (int i = begin; i < end; ++i) {
            Vector3 accel = gravity_direction * gravity_strength;
            Vector3 vel = particles[i].velocity;
            vel = vel + accel * dt;
            particles[i].velocity = vel;
        }
    }

    unsafe void Update_Position(float dt, int begin, int end)
    {
        peParticle *particles = peGetParticles(m_ctx);
        for (int i = begin; i < end; ++i) {
            Vector3 pos = particles[i].position;
            Vector3 vel = particles[i].velocity;
            pos = pos + (vel * dt);
            particles[i].position = pos;
        }
    }

    public void Benchmark()
    {
        peSetCSUpdateRoutine(m_ctx, Update_Velocity, Update_Position);
        IntPtr str = peBenchmark(m_ctx, 1);
        string result = Marshal.PtrToStringAnsi(str);
        Debug.Log(result);
    }
}
