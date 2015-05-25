using UnityEngine;
using System.Collections;

[RequireComponent(typeof(UnityEngine.UI.Text))]
public class FrameRateCounter : MonoBehaviour
{
    public PEParticles m_particles;
    public float m_update_interval = 0.5f;
    private float m_last_time;
    private float m_accum = 0.0f; // FPS accumulated over the interval
    private int m_frames = 0; // Frames drawn over the interval
    private float m_time_left; // Left time for current interval
    private float m_fps;
    private float m_average_particle_update_time;
    float m_accum_p;

    void Start()
    {
        m_time_left = m_update_interval;
        m_last_time = Time.realtimeSinceStartup;
    }

    void Update()
    {
        if (m_particles != null)
        {
            m_accum_p += m_particles.m_update_time;
        }
        float now = Time.realtimeSinceStartup;
        float delta = now - m_last_time;
        m_last_time = now;
        m_time_left -= delta;
        m_accum += 1.0f / delta;
        ++m_frames;

        // Interval ended - update result
        if (m_time_left <= 0.0)
        {
            m_fps = m_accum / m_frames;
            m_average_particle_update_time = m_accum_p / m_frames * 1000.0f;
            string t = m_fps.ToString("f2") + " FPS\n";
            t += m_average_particle_update_time.ToString("f2") + " ms";
            GetComponent<UnityEngine.UI.Text>().text = t;
            m_time_left = m_update_interval;
            m_accum = 0.0f;
            m_accum_p = 0.0f;
            m_frames = 0;
        }
    }
}
