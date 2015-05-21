using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;
#if UNITY_EDITOR
using UnityEditor;
#endif


[RequireComponent(typeof(PEParticles))]
public class PERenderer : BatchRendererBase
{
    PEParticles m_particles;
    public RenderTexture m_data_texture;


    public RenderTexture GetInstanceTexture() { return m_data_texture; }

#if UNITY_EDITOR
    void Reset()
    {
        m_mesh = AssetDatabase.LoadAssetAtPath("Assets/BatchRenderer/Meshes/cube.asset", typeof(Mesh)) as Mesh;
        m_material = AssetDatabase.LoadAssetAtPath("Assets/ParticleEngine/Materials/PEStandard.mat", typeof(Material)) as Material;
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

        m_data_texture = new RenderTexture(PEParticles.DataTextureWidth * 2, 128, 0, RenderTextureFormat.ARGBFloat);
        m_data_texture.filterMode = FilterMode.Point;
        m_data_texture.enableRandomWrite = true;
        m_data_texture.generateMips = false;
        m_data_texture.Create();

        UpdateGPUResources();
    }

    public override void UpdateGPUResources()
    {
        bool use_buffer = m_particles.m_routine == PEParticles.peUpdateRoutine.ComputeShader && m_particles.m_cb_particles != null;
        if (!use_buffer)
        {
            m_particles.CopyDataToTexture(m_data_texture);
        }

        ForEachEveryMaterials((v) =>
        {
            v.SetInt("g_num_max_instances", m_max_instances);
            v.SetInt("g_num_instances", m_instance_count);
            v.SetInt("g_use_buffer", use_buffer ? 1 : 0);
            if (use_buffer)
            {
                v.SetBuffer("g_instance_buffer", m_particles.m_cb_particles);
            }
        });
    }


    public override void OnEnable()
    {
        m_particles = GetComponent<PEParticles>();
        m_max_instances = 16384;

        base.OnEnable();
        ResetGPUResoures();
    }

    public override void OnDisable()
    {
        base.OnDisable();
        ReleaseGPUResources();
    }

    public override void LateUpdate()
    {
        m_instance_count = m_particles.m_particle_count;
        base.LateUpdate();
    }

    public override void OnDrawGizmos()
    {
    }
}
