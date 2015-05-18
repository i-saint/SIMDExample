using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;
#if UNITY_EDITOR
using UnityEditor;
#endif


[AddComponentMenu("MassParticle/Renderer")]
public class PERenderer : BatchRendererBase
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

    [DllImport ("ParticleEngine")] public static extern void    peSetUpdateRoutine(IntPtr ctx, peUpdateRoutine r);
    [DllImport ("ParticleEngine")] public static extern void    peUpdate(IntPtr ctx, float dt);
    [DllImport ("ParticleEngine")] public static extern void    peCopyDataToTexture(IntPtr ctx, IntPtr texture, int width, int height);



    public int m_particle_count = 1024*4;
    public float m_particle_size = 0.1f;
    public peUpdateRoutine m_routine = peUpdateRoutine.ISPC;
    IntPtr m_ctx;
    RenderTexture m_data_texture;
    Bounds m_bounds;


    public RenderTexture GetInstanceTexture() { return m_data_texture; }

#if UNITY_EDITOR
    void Reset()
    {
        //m_mesh = AssetDatabase.LoadAssetAtPath("Assets/BatchRenderer/Meshes/cube.asset", typeof(Mesh)) as Mesh;
        //m_material = AssetDatabase.LoadAssetAtPath("Assets/MassParticle/Materials/MPStandard.mat", typeof(Material)) as Material;
        m_bounds_size = Vector3.one * 2.0f;
    }
#endif


    public override Material CloneMaterial(Material src, int nth)
    {
        Material m = new Material(src);
        m.SetInt("g_batch_begin", nth * m_instances_par_batch);
        m.SetTexture("g_instance_data", m_data_texture);

        Vector4 ts = new Vector4(
            1.0f / m_data_texture.width,
            1.0f / m_data_texture.height,
            m_data_texture.width,
            m_data_texture.height);
        m.SetVector("g_instance_data_size", ts);

        // fix rendering order for transparent objects
        if (m.renderQueue >= 3000)
        {
            m.renderQueue = m.renderQueue + (nth + 1);
        }
        return m;
    }


    public virtual void ReleaseGPUResources()
    {
        if (m_data_texture != null)
        {
            m_data_texture.Release();
            m_data_texture = null;
        }
        ClearMaterials();
    }

    public virtual void ResetGPUResoures()
    {
        ReleaseGPUResources();

        m_data_texture = new RenderTexture(1024*2, 64, 0, RenderTextureFormat.ARGBFloat, RenderTextureReadWrite.Default);
        m_data_texture.filterMode = FilterMode.Point;
        m_data_texture.Create();

        UpdateGPUResources();
    }

    public override void UpdateGPUResources()
    {
        peCopyDataToTexture(m_ctx, m_data_texture.GetNativeTexturePtr(), m_data_texture.width, m_data_texture.height);

        ForEachEveryMaterials((v) =>
        {
            v.SetInt("g_num_max_instances", m_max_instances);
            v.SetInt("g_num_instances", m_instance_count);
        });
    }


    public override void OnEnable()
    {
        m_ctx = peCreateContext(m_particle_count);

        base.OnEnable();
        ResetGPUResoures();
    }

    public override void OnDisable()
    {
        base.OnDisable();
        ReleaseGPUResources();
        peDestroyContext(m_ctx);
    }

    public override void LateUpdate()
    {
        base.LateUpdate();
    }

    public override void OnDrawGizmos()
    {
    }
}
