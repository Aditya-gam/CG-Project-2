#include "shaders.h"

// Lookup maps to access a shader by name.
std::map<std::string,shader_v> vertex_shader_map;
std::map<std::string,shader_f> fragment_shader_map;

// Simplest useful vertex shader; just copies over the positions.
void vertex_shader_trivial(const data_vertex& in, data_geometry& out,
    const float * uniform_data)
{
    const vertex_p& v = *(const vertex_p*)in.data;
    out.gl_Position = vec4(v.position,1);
}

// This shader reorders the input.
void vertex_shader_reorder(const data_vertex& in, data_geometry& out,
    const float * uniform_data)
{
    const vertex_r& v = *(const vertex_r*)in.data;
    out.gl_Position = vec4(v.position,1);
    std::swap(out.gl_Position[0],out.gl_Position[1]);
}

// This shader (1) transforms vertices and (2) copies over colors.
void vertex_shader_color(const data_vertex& in, data_geometry& out,
    const float * uniform_data)
{
    const vertex_pc& v = *(const vertex_pc*)in.data;
    vertex_pc& o = *(vertex_pc*)out.data;
    mat4& xform = *(mat4*)uniform_data;
    out.gl_Position = xform * vec4(v.position,1);
    o.color = v.color;
}

// This shader (1) transforms vertices and (2) copies over colors.
void vertex_shader_color2(const data_vertex& in, data_geometry& out, const float * uniform_data)
{
    const vertex_pc2& v = *(const vertex_pc2*)in.data;
    
    // Allocate memory for the output attributes
    out.data = new float[MAX_FLOATS_PER_VERTEX]();

    // Apply transformation
    mat4& xform = *(mat4*)uniform_data;
    out.gl_Position = xform * vec4(v.position,1);

    // Copy color data
    out.data[0] = v.color[0];  // R
    out.data[1] = v.color[1];  // G
    out.data[2] = v.color[2];  // B
}


// This shader just transforms vertices
void vertex_shader_transform(const data_vertex& in, data_geometry& out,
    const float * uniform_data)
{
    // Extract the transformation matrix from uniform_data
    const mat4& xform = *(const mat4*)uniform_data;

    // Read the input vertex position
    const vertex_p& v = *(const vertex_p*)in.data;

    // Apply the transformation matrix (Projection * ModelView * Vertex)
    out.gl_Position = xform * vec4(v.position, 1.0f);
    
    // Perform **homogeneous division** (perspective divide)
    if (out.gl_Position[3] != 0.0f) {
        out.gl_Position[0] /= out.gl_Position[3];
        out.gl_Position[1] /= out.gl_Position[3];
        out.gl_Position[2] /= out.gl_Position[3];  // Depth normalization
    }
}


// Simple fragment shader: set the fragment to red
void fragment_shader_red(const data_fragment& in, data_output& out,
    const float * uniform_data)
{
    out.output_color = vec4(1,0,0,0);
}

// Simple fragment shader: set the fragment to green
void fragment_shader_green(const data_fragment& in, data_output& out,
    const float * uniform_data)
{
    out.output_color = vec4(0,1,0,0);
}

// Simple fragment shader: set the fragment to blue
void fragment_shader_blue(const data_fragment& in, data_output& out,
    const float * uniform_data)
{
    out.output_color = vec4(0,0,1,0);
}

// Simple fragment shader: set the fragment to white
void fragment_shader_white(const data_fragment& in, data_output& out,
    const float * uniform_data)
{
    out.output_color = vec4(1,1,1,0);
}

// Simple fragment shader: pass through globally constant color
void fragment_shader_uniform(const data_fragment& in, data_output& out,
    const float * uniform_data)
{
    transform_color& tc = *(transform_color*)uniform_data;
    out.output_color = vec4(tc.color,0);
}

// Simple fragment shader: pass through interpolated per-vertex color
void fragment_shader_gouraud(const data_fragment& in, data_output& out,
    const float * uniform_data)
{
    vertex_pc& v = *(vertex_pc*)in.data;
    out.output_color = vec4(v.color,0);
}

// Simple fragment shader: pass through interpolated per-vertex color
void fragment_shader_gouraud2(const data_fragment& in, data_output& out, const float * uniform_data)
{
    if (!in.data) {
        out.output_color = vec4(1, 1, 1, 1);  // Default to white if no data
        return;
    }

    vec3 color;
    color[0] = in.data[0];  // R
    color[1] = in.data[1];  // G
    color[2] = in.data[2];  // B

    out.output_color = vec4(color, 1);
}

// Assign shaders to the maps so they can be accessed by name.
void register_named_shaders()
{
    vertex_shader_map["trivial"]=vertex_shader_trivial;
    vertex_shader_map["reorder"]=vertex_shader_reorder;
    vertex_shader_map["transform"]=vertex_shader_transform;
    vertex_shader_map["color"]=vertex_shader_color;
    vertex_shader_map["color2"]=vertex_shader_color2;
    fragment_shader_map["red"]=fragment_shader_red;
    fragment_shader_map["green"]=fragment_shader_green;
    fragment_shader_map["blue"]=fragment_shader_blue;
    fragment_shader_map["white"]=fragment_shader_white;
    fragment_shader_map["gouraud"]=fragment_shader_gouraud;
    fragment_shader_map["gouraud2"]=fragment_shader_gouraud2;
    fragment_shader_map["uniform"]=fragment_shader_uniform;
}
