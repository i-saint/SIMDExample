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
    [DllImport ("ParticleEngine")] public static extern void    peUpdate(IntPtr ctx, float dt);
    [DllImport ("ParticleEngine")] public static extern void    peCopyDataToTexture(IntPtr ctx, IntPtr texture, int width, int height);


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
        peDestroyContext(m_ctx);
        m_particle_count = v;
        m_ctx = peCreateContext(m_particle_count);
    }

    public void EnableMultiThreading(bool v)
    {
        m_multi_threading = v;
    }

    public void CopyDataToTexture(Texture tex)
    {
        peCopyDataToTexture(m_ctx, tex.GetNativeTexturePtr(), tex.width, tex.height);
    }


    void OnEnable()
    {
        m_particle_count = Mathf.Max(DataTextureWidth, m_particle_count);
        m_ctx = peCreateContext(m_particle_count);
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
        peUpdate(m_ctx, Time.deltaTime);
    }
}
