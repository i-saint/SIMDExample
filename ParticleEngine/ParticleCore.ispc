#include "Vector_ispc.h"



struct Context
{
   float *pos_x;
   float *pos_y;
   float *pos_z;
   float *vel_x;
   float *vel_y;
   float *vel_z;

   int particle_count;
   float particle_size;
   float rcp_particle_size2;
   float pressure_stiffness;
   float wall_stiffness;
   float timestep;
};


export void UpdateVelocity(uniform Context &ctx, uniform int begin, uniform int end)
{
    // パーティクル同士の押し返し
    for(uniform int i=begin; i<end; ++i) {
        uniform float3 pos1 = {ctx.pos_x[i], ctx.pos_y[i], ctx.pos_z[i] };
        float3 accel = {0.0f, 0.0f, 0.0f};
        for(uniform int j=0; j<ctx.particle_count; j+=programCount) {
            int js = j + programIndex;
            float3 pos2 = {ctx.pos_x[js], ctx.pos_y[js], ctx.pos_z[js] };
            float3 diff = pos2 - pos1;
            float distance = length(diff);
            float3 dir = diff * ctx.rcp_particle_size2;
            float3 a = dir * (min(0.0f, distance-(ctx.particle_size*2.0f)) * ctx.pressure_stiffness);
            accel = accel + (a * select(distance > 0.0f, 1.0f, 0.0f));
        }

        uniform float3 vel = {ctx.vel_x[i], ctx.vel_y[i], ctx.vel_z[i] };
        vel = vel + reduce_add(accel) * ctx.timestep;
        ctx.vel_x[i] = vel.x;
        ctx.vel_y[i] = vel.y;
        ctx.vel_z[i] = vel.z;
    }

    // 床との衝突
    uniform float3 floor_normal = {0.0f, 1.0f, 0.0f};
    uniform float floor_distance = -ctx.particle_size;
    for (uniform int i=begin; i<end; i+=programCount) {
        int is = i + programIndex;
        float3 pos = {ctx.pos_x[is], ctx.pos_y[is], ctx.pos_z[is] };

        float d = dot(pos, floor_normal) + floor_distance;
        float3 accel = floor_normal * (-min(0.0f, d) * ctx.wall_stiffness);

        float3 vel = {ctx.vel_x[is], ctx.vel_y[is], ctx.vel_z[is] };
        vel = vel + accel * ctx.timestep;
        ctx.vel_x[is] = vel.x;
        ctx.vel_y[is] = vel.y;
        ctx.vel_z[is] = vel.z;
    }

    // 重力加速
    uniform float3 gravity_direction = {0.0f, -1.0f, 0.0f};
    uniform float gravity_strength = 5.0f;
    for (uniform int i=begin; i<end; i+=programCount) {
        int is = i + programIndex;
        float3 accel = gravity_direction * gravity_strength;
        float3 vel = {ctx.vel_x[is], ctx.vel_y[is], ctx.vel_z[is] };
        vel = vel + accel * ctx.timestep;
        ctx.vel_x[is] = vel.x;
        ctx.vel_y[is] = vel.y;
        ctx.vel_z[is] = vel.z;
    }
}

export void UpdatePosition(uniform Context &ctx, uniform int begin, uniform int end)
{
    for(uniform int i=begin; i<end; i+=programCount) {
        int is = i + programIndex;
        float3 pos = {ctx.pos_x[is], ctx.pos_y[is], ctx.pos_z[is] };
        float3 vel = {ctx.vel_x[is], ctx.vel_y[is], ctx.vel_z[is] };

        pos = pos + (vel * ctx.timestep);
        ctx.pos_x[is] = pos.x;
        ctx.pos_y[is] = pos.y;
        ctx.pos_z[is] = pos.z;
    }
}
