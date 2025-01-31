#include "driver_state.h"
#include <cstring>
#include <vector>  // Ensure this is included
#include <algorithm>  // For std::min, std::max



driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width = width;
    state.image_height = height;

    // Allocate memory for image_color and image_depth
    state.image_color = new pixel[width * height];
    state.image_depth = new float[width * height];

    // Initialize image_color to black (0,0,0)
    for (int i = 0; i < width * height; i++)
    {
        state.image_color[i] = make_pixel(0, 0, 0);  // Black color
    }

    // Initialize image_depth to maximum depth value (1.0 for a z-buffer)
    std::fill_n(state.image_depth, width * height, 1.0f);
}


// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    if (state.num_vertices == 0 || state.floats_per_vertex == 0) {
        std::cout << "No vertex data available for rendering." << std::endl;
        return;
    }

    switch (type) {
        case render_type::triangle:
            for (int i = 0; i < state.num_vertices; i += 3) {
                data_geometry v0, v1, v2;
                data_vertex in0, in1, in2;

                // Allocate memory for attributes
                v0.data = new float[state.floats_per_vertex];
                v1.data = new float[state.floats_per_vertex];
                v2.data = new float[state.floats_per_vertex];

                // Set up vertex data pointers
                in0.data = &state.vertex_data[i * state.floats_per_vertex];
                in1.data = &state.vertex_data[(i + 1) * state.floats_per_vertex];
                in2.data = &state.vertex_data[(i + 2) * state.floats_per_vertex];

                // Copy attributes to geometry data
                std::memcpy(v0.data, in0.data, sizeof(float) * state.floats_per_vertex);
                std::memcpy(v1.data, in1.data, sizeof(float) * state.floats_per_vertex);
                std::memcpy(v2.data, in2.data, sizeof(float) * state.floats_per_vertex);

                // Call vertex shader
                state.vertex_shader(in0, v0, state.uniform_data);
                state.vertex_shader(in1, v1, state.uniform_data);
                state.vertex_shader(in2, v2, state.uniform_data);

                // Clip and rasterize
                clip_triangle(state, v0, v1, v2);

                // Free allocated memory
                delete[] v0.data;
                delete[] v1.data;
                delete[] v2.data;
            }
            break;

        case render_type::indexed:
            for (int i = 0; i < state.num_triangles * 3; i += 3) {
                data_geometry v0, v1, v2;
                data_vertex in0, in1, in2;

                // Get vertex indices from the index buffer
                int idx0 = state.index_data[i];
                int idx1 = state.index_data[i + 1];
                int idx2 = state.index_data[i + 2];

                // Allocate memory for attributes
                v0.data = new float[state.floats_per_vertex];
                v1.data = new float[state.floats_per_vertex];
                v2.data = new float[state.floats_per_vertex];

                // Set up vertex data pointers using indices
                in0.data = &state.vertex_data[idx0 * state.floats_per_vertex];
                in1.data = &state.vertex_data[idx1 * state.floats_per_vertex];
                in2.data = &state.vertex_data[idx2 * state.floats_per_vertex];

                // Copy attributes to geometry data
                std::memcpy(v0.data, in0.data, sizeof(float) * state.floats_per_vertex);
                std::memcpy(v1.data, in1.data, sizeof(float) * state.floats_per_vertex);
                std::memcpy(v2.data, in2.data, sizeof(float) * state.floats_per_vertex);

                // Call vertex shader
                state.vertex_shader(in0, v0, state.uniform_data);
                state.vertex_shader(in1, v1, state.uniform_data);
                state.vertex_shader(in2, v2, state.uniform_data);

                // Clip and rasterize
                clip_triangle(state, v0, v1, v2);

                // Free allocated memory
                delete[] v0.data;
                delete[] v1.data;
                delete[] v2.data;
            }
            break;

        case render_type::fan:
            if (state.num_vertices < 3) return;
            for (int i = 1; i < state.num_vertices - 1; ++i) {
                data_geometry v0, v1, v2;
                data_vertex in0, in1, in2;

                v0.data = new float[state.floats_per_vertex];
                v1.data = new float[state.floats_per_vertex];
                v2.data = new float[state.floats_per_vertex];

                in0.data = &state.vertex_data[0];
                in1.data = &state.vertex_data[i * state.floats_per_vertex];
                in2.data = &state.vertex_data[(i + 1) * state.floats_per_vertex];

                std::memcpy(v0.data, in0.data, sizeof(float) * state.floats_per_vertex);
                std::memcpy(v1.data, in1.data, sizeof(float) * state.floats_per_vertex);
                std::memcpy(v2.data, in2.data, sizeof(float) * state.floats_per_vertex);

                state.vertex_shader(in0, v0, state.uniform_data);
                state.vertex_shader(in1, v1, state.uniform_data);
                state.vertex_shader(in2, v2, state.uniform_data);

                clip_triangle(state, v0, v1, v2);

                delete[] v0.data;
                delete[] v1.data;
                delete[] v2.data;
            }
            break;

        case render_type::strip:
            if (state.num_vertices < 3) return;
            for (int i = 0; i < state.num_vertices - 2; ++i) {
                data_geometry v0, v1, v2;
                data_vertex in0, in1, in2;

                v0.data = new float[state.floats_per_vertex];
                v1.data = new float[state.floats_per_vertex];
                v2.data = new float[state.floats_per_vertex];

                in0.data = &state.vertex_data[i * state.floats_per_vertex];
                in1.data = &state.vertex_data[(i + 1) * state.floats_per_vertex];
                in2.data = &state.vertex_data[(i + 2) * state.floats_per_vertex];

                std::memcpy(v0.data, in0.data, sizeof(float) * state.floats_per_vertex);
                std::memcpy(v1.data, in1.data, sizeof(float) * state.floats_per_vertex);
                std::memcpy(v2.data, in2.data, sizeof(float) * state.floats_per_vertex);

                state.vertex_shader(in0, v0, state.uniform_data);
                state.vertex_shader(in1, v1, state.uniform_data);
                state.vertex_shader(in2, v2, state.uniform_data);

                clip_triangle(state, v0, v1, v2);

                delete[] v0.data;
                delete[] v1.data;
                delete[] v2.data;
            }
            break;

        default:
            std::cout << "Invalid render type!" << std::endl;
            return;
    }
}




data_geometry perspective_interpolate(const data_geometry& in1, const data_geometry& in2, float w1, float w2, const interp_type interp_rules[MAX_FLOATS_PER_VERTEX])
{
    data_geometry result;

    // Compute interpolation factor `t` based on homogeneous w values
    float d1 = in1.gl_Position[3];  // w_A
    float d2 = in2.gl_Position[3];  // w_B

    float t = d1 / (d1 - d2);  // Correct interpolation factor

    // Interpolate gl_Position
    result.gl_Position = in1.gl_Position * (1 - t) + in2.gl_Position * t;

    // Allocate new memory for attributes
    result.data = new float[MAX_FLOATS_PER_VERTEX];

    // Compute α', β', γ' in screen space (unnormalized barycentric weights)
    float alpha_prime = (1 - t) * d1;
    float beta_prime = t * d2;

    float k = (alpha_prime / d1) + (beta_prime / d2);  // Normalization factor

    float alpha = (alpha_prime / d1) / k;
    float beta = (beta_prime / d2) / k;

    for (int i = 0; i < MAX_FLOATS_PER_VERTEX; i++) {
        if (interp_rules[i] == interp_type::noperspective) {
            // No perspective correction: Interpolate linearly
            result.data[i] = in1.data[i] * (1 - t) + in2.data[i] * t;
        } else {
            // Perspective correction: Interpolate using corrected weights
            float attr1 = in1.data[i] / d1;
            float attr2 = in2.data[i] / d2;
            float interpolated_attr = (alpha * attr1) + (beta * attr2);
            result.data[i] = interpolated_attr * result.gl_Position[3];  // Convert back to perspective space
        }
    }

    return result;
}









data_geometry interpolate(const data_geometry& a, const data_geometry& b, float t) {
    data_geometry result;
    result.gl_Position = a.gl_Position * (1 - t) + b.gl_Position * t;
    result.data = new float[MAX_FLOATS_PER_VERTEX];
    for (int i = 0; i < MAX_FLOATS_PER_VERTEX; i++) {
        result.data[i] = a.data[i] * (1 - t) + b.data[i] * t;
    }
    return result;
}



void clip_triangle(driver_state& state, const data_geometry& v0,
                   const data_geometry& v1, const data_geometry& v2, int face)
{
    if (face == 6) {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }

    std::vector<data_geometry> inside, outside;
    std::vector<float> inside_w, outside_w;

    // Clipping planes: Near (z = -w) and Far (z = w)
    vec4 planes[6] = {
        vec4(1, 0, 0, 1),  vec4(-1, 0, 0, 1),  // Left and Right
        vec4(0, 1, 0, 1),  vec4(0, -1, 0, 1), // Top and Bottom
        vec4(0, 0, 1, 1),  vec4(0, 0, -1, 1)  // Near and Far
    };

    vec4 plane = planes[face];

    const data_geometry* vertices[3] = { &v0, &v1, &v2 };

    for (int i = 0; i < 3; i++) {
        float d = dot(vertices[i]->gl_Position, plane);
        if (d >= 0) {
            inside.push_back(*vertices[i]);
            inside_w.push_back(vertices[i]->gl_Position[3]);
        } else {
            outside.push_back(*vertices[i]);
            outside_w.push_back(vertices[i]->gl_Position[3]);
        }
    }

    if (inside.size() == 3) {
        clip_triangle(state, v0, v1, v2, face + 1);
    } 
    else if (inside.size() == 2 && outside.size() == 1) {
        data_geometry new_v1 = perspective_interpolate(inside[0], outside[0], inside_w[0], outside_w[0], state.interp_rules);
        data_geometry new_v2 = perspective_interpolate(inside[1], outside[0], inside_w[1], outside_w[0], state.interp_rules);

        clip_triangle(state, inside[0], inside[1], new_v1, face + 1);
        clip_triangle(state, inside[1], new_v1, new_v2, face + 1);
    } 
    else if (inside.size() == 1 && outside.size() == 2) {
        data_geometry new_v1 = perspective_interpolate(inside[0], outside[0], inside_w[0], outside_w[0], state.interp_rules);
        data_geometry new_v2 = perspective_interpolate(inside[0], outside[1], inside_w[0], outside_w[1], state.interp_rules);

        clip_triangle(state, inside[0], new_v1, new_v2, face + 1);
    }
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
vec3 to_screen_space(const vec4& ndc, int width, int height) {
    return vec3(
        (ndc[0] / ndc[3] * 0.5f + 0.5f) * width,
        (ndc[1] / ndc[3] * 0.5f + 0.5f) * height,
        ndc[2] / ndc[3]
    );
}

vec3 compute_barycentric(const vec3& A, const vec3& B, const vec3& C, const vec3& P) {
    vec3 v0 = B - A, v1 = C - A, v2 = P - A;
    float d00 = dot(v0, v0);
    float d01 = dot(v0, v1);
    float d11 = dot(v1, v1);
    float d20 = dot(v2, v0);
    float d21 = dot(v2, v1);
    float denom = d00 * d11 - d01 * d01;
    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;
    return vec3(u, v, w);
}

void rasterize_triangle(driver_state& state, const data_geometry& v0,
                        const data_geometry& v1, const data_geometry& v2)
{
    // Convert NDC coordinates to screen space
    auto ndc_to_screen = [&](const vec4& v) -> vec3 {
        return vec3(
            (v[0] * 0.5f + 0.5f) * state.image_width,
            (v[1] * 0.5f + 0.5f) * state.image_height,
            v[2]
        );
    };

    vec3 p0 = ndc_to_screen(v0.gl_Position);
    vec3 p1 = ndc_to_screen(v1.gl_Position);
    vec3 p2 = ndc_to_screen(v2.gl_Position);

    int min_x = std::max(0, (int)std::min({p0[0], p1[0], p2[0]}));
    int max_x = std::min(state.image_width - 1, (int)std::max({p0[0], p1[0], p2[0]}));
    int min_y = std::max(0, (int)std::min({p0[1], p1[1], p2[1]}));
    int max_y = std::min(state.image_height - 1, (int)std::max({p0[1], p1[1], p2[1]}));

    auto edge_function = [](const vec3& a, const vec3& b, const vec3& c) -> float {
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    };

    float area = edge_function(p0, p1, p2);
    if (area == 0) return; // Skip degenerate triangles

    for (int y = min_y; y <= max_y; y++) {
        for (int x = min_x; x <= max_x; x++) {
            vec3 p(x + 0.5f, y + 0.5f, 0);

            float w0 = edge_function(p1, p2, p) / area;
            float w1 = edge_function(p2, p0, p) / area;
            float w2 = edge_function(p0, p1, p) / area;

            if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
                int pixel_index = y * state.image_width + x;

                // **Fix Perspective Depth Interpolation**
                float z0 = v0.gl_Position[2] / v0.gl_Position[3];
                float z1 = v1.gl_Position[2] / v1.gl_Position[3];
                float z2 = v2.gl_Position[2] / v2.gl_Position[3];

                float depth = w0 * z0 + w1 * z1 + w2 * z2;

                if (depth < state.image_depth[pixel_index]) {
                    state.image_depth[pixel_index] = depth;

                    // Allocate memory for interpolated attributes
                    data_fragment fragment;
                    fragment.data = new float[state.floats_per_vertex];

                    // **Corrected Attribute Interpolation**
                    for (int i = 0; i < state.floats_per_vertex; i++) {
                        if (state.interp_rules[i] == interp_type::smooth) {
                            // Perspective-correct interpolation
                            float w_div_z = w0 / v0.gl_Position[3] + w1 / v1.gl_Position[3] + w2 / v2.gl_Position[3];
                            fragment.data[i] = (w0 * v0.data[i] / v0.gl_Position[3] +
                                                w1 * v1.data[i] / v1.gl_Position[3] +
                                                w2 * v2.data[i] / v2.gl_Position[3]) / w_div_z;
                        } else if (state.interp_rules[i] == interp_type::noperspective) {
                            // Image-space barycentric interpolation
                            fragment.data[i] = w0 * v0.data[i] + w1 * v1.data[i] + w2 * v2.data[i];
                        } else if (state.interp_rules[i] == interp_type::flat) {
                            // Flat interpolation (first vertex attribute)
                            fragment.data[i] = v0.data[i];
                        }
                    }

                    // Process the fragment
                    data_output output;
                    state.fragment_shader(fragment, output, state.uniform_data);

                    int r = std::min(255, static_cast<int>(output.output_color[0] * 255));
                    int g = std::min(255, static_cast<int>(output.output_color[1] * 255));
                    int b = std::min(255, static_cast<int>(output.output_color[2] * 255));

                    state.image_color[pixel_index] = make_pixel(r, g, b);

                    // Cleanup allocated memory
                    delete[] fragment.data;
                }
            }
        }
    }
}
