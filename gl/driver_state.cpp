#include "driver_state.h"
#include <cstring>
#include <vector>    // For std::vector
#include <algorithm> // For std::min, std::max

// -----------------------------------------------------------------------------
// driver_state Constructor and Destructor
// -----------------------------------------------------------------------------
driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// -----------------------------------------------------------------------------
// initialize_render: Allocate and initialize the render buffers.
// -----------------------------------------------------------------------------
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width = width;
    state.image_height = height;

    // Allocate memory for the color and depth buffers.
    state.image_color = new pixel[width * height];
    state.image_depth = new float[width * height];

    // Initialize the color buffer to black.
    for (int i = 0; i < width * height; i++) {
        state.image_color[i] = make_pixel(0, 0, 0);
    }

    // Initialize the depth buffer to 1.0 (the farthest depth).
    std::fill_n(state.image_depth, width * height, 1.0f);
}

// -----------------------------------------------------------------------------
// render: Render the geometry stored in the driver_state.
// -----------------------------------------------------------------------------
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

                // Allocate attribute arrays for each vertex.
                v0.data = new float[state.floats_per_vertex];
                v1.data = new float[state.floats_per_vertex];
                v2.data = new float[state.floats_per_vertex];

                // Set pointers into the vertex_data array.
                in0.data = &state.vertex_data[i * state.floats_per_vertex];
                in1.data = &state.vertex_data[(i + 1) * state.floats_per_vertex];
                in2.data = &state.vertex_data[(i + 2) * state.floats_per_vertex];

                // Copy raw vertex data to geometry data.
                std::memcpy(v0.data, in0.data, sizeof(float) * state.floats_per_vertex);
                std::memcpy(v1.data, in1.data, sizeof(float) * state.floats_per_vertex);
                std::memcpy(v2.data, in2.data, sizeof(float) * state.floats_per_vertex);

                // Run the vertex shader on each vertex.
                state.vertex_shader(in0, v0, state.uniform_data);
                state.vertex_shader(in1, v1, state.uniform_data);
                state.vertex_shader(in2, v2, state.uniform_data);

                // Clip (which may generate new vertices) and then rasterize.
                clip_triangle(state, v0, v1, v2);

                // Free the allocated attribute arrays.
                delete[] v0.data;
                delete[] v1.data;
                delete[] v2.data;
            }
            break;

        case render_type::indexed:
            for (int i = 0; i < state.num_triangles * 3; i += 3) {
                data_geometry v0, v1, v2;
                data_vertex in0, in1, in2;

                // Retrieve vertex indices from the index buffer.
                int idx0 = state.index_data[i];
                int idx1 = state.index_data[i + 1];
                int idx2 = state.index_data[i + 2];

                // Allocate attribute arrays.
                v0.data = new float[state.floats_per_vertex];
                v1.data = new float[state.floats_per_vertex];
                v2.data = new float[state.floats_per_vertex];

                // Set up vertex pointers using the indices.
                in0.data = &state.vertex_data[idx0 * state.floats_per_vertex];
                in1.data = &state.vertex_data[idx1 * state.floats_per_vertex];
                in2.data = &state.vertex_data[idx2 * state.floats_per_vertex];

                // Copy attribute data.
                std::memcpy(v0.data, in0.data, sizeof(float) * state.floats_per_vertex);
                std::memcpy(v1.data, in1.data, sizeof(float) * state.floats_per_vertex);
                std::memcpy(v2.data, in2.data, sizeof(float) * state.floats_per_vertex);

                // Run the vertex shader.
                state.vertex_shader(in0, v0, state.uniform_data);
                state.vertex_shader(in1, v1, state.uniform_data);
                state.vertex_shader(in2, v2, state.uniform_data);

                // Clip and rasterize.
                clip_triangle(state, v0, v1, v2);

                // Free allocated attribute arrays.
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


// -----------------------------------------------------------------------------
// perspective_interpolate: Compute the intersection vertex along an edge crossing a clip plane.
// -----------------------------------------------------------------------------
/**
 * perspective_interpolate
 *
 * Given two geometry vertices (forming an edge that crosses a clipping plane),
 * this function computes the intersection vertex.
 *
 * For attributes flagged as "smooth", perspective–correct interpolation is performed
 * using the clip–space interpolation factor t_clip (computed from the clip–space distances).
 *
 * For attributes flagged as "noperspective", an interpolation factor t_ndc is computed
 * in normalized device coordinates (NDC) using the appropriate coordinate for the current
 * clip plane (determined by the parameter 'face'). This ensures that non–perspective
 * interpolation remains linear in screen space.
 *
 * For flat interpolation the attribute from the first vertex is used.
 *
 * Parameters:
 *   in1, in2       - The two geometry vertices defining the edge.
 *   d1, d2         - The signed distances (dot(gl_Position, plane)) for in1 and in2.
 *   interp_rules   - The interpolation rules (one per attribute).
 *   face           - The index of the current clip plane:
 *                      0: Near, 1: Far, 2: Left, 3: Right, 4: Bottom, 5: Top.
 *
 * Returns:
 *   A new geometry vertex at the intersection of the edge with the clip plane.
 */
data_geometry perspective_interpolate(const data_geometry& in1, const data_geometry& in2,
    float d1, float d2, const interp_type interp_rules[MAX_FLOATS_PER_VERTEX], int face)
{
    data_geometry result;

    // Compute the clip–space interpolation factor.
    float t_clip = d1 / (d1 - d2);

    // For noperspective interpolation, compute the interpolation factor in NDC.
    float t_ndc = t_clip; // Default value.
    {
        // Compute normalized device coordinates for both vertices.
        vec4 ndc1 = in1.gl_Position / in1.gl_Position[3];
        vec4 ndc2 = in2.gl_Position / in2.gl_Position[3];
        switch(face) {
            case 0: // Near plane: use z coordinate, expected ndc.z >= -1
                t_ndc = (ndc1[2] - (-1)) / (ndc1[2] - ndc2[2]);
                break;
            case 1: // Far plane: ndc.z <= 1
                t_ndc = (1 - ndc1[2]) / (ndc2[2] - ndc1[2]);
                break;
            case 2: // Left plane: ndc.x >= -1
                t_ndc = (ndc1[0] - (-1)) / (ndc1[0] - ndc2[0]);
                break;
            case 3: // Right plane: ndc.x <= 1
                t_ndc = (1 - ndc1[0]) / (ndc2[0] - ndc1[0]);
                break;
            case 4: // Bottom plane: ndc.y >= -1
                t_ndc = (ndc1[1] - (-1)) / (ndc1[1] - ndc2[1]);
                break;
            case 5: // Top plane: ndc.y <= 1
                t_ndc = (1 - ndc1[1]) / (ndc2[1] - ndc1[1]);
                break;
            default:
                t_ndc = t_clip;
                break;
        }
    }

    // Interpolate the clip–space position using t_clip.
    result.gl_Position = in1.gl_Position * (1 - t_clip) + in2.gl_Position * t_clip;
    // Allocate memory for the interpolated attribute data.
    result.data = new float[MAX_FLOATS_PER_VERTEX];
    for (int i = 0; i < MAX_FLOATS_PER_VERTEX; i++) {
        switch(interp_rules[i]) {
            case interp_type::flat:
                // For flat interpolation, always use the attribute from the provoking vertex.
                // (In our implementation the provoking vertex is the first vertex of the original triangle.)
                result.data[i] = in1.data[i];
                break;
            case interp_type::noperspective:
                // Noperspective interpolation: use the NDC–based interpolation factor.
                result.data[i] = in1.data[i] * (1 - t_ndc) + in2.data[i] * t_ndc;
                break;
            case interp_type::smooth:
                {
                    // Perspective–correct (smooth) interpolation.
                    float w1 = in1.gl_Position[3];
                    float w2 = in2.gl_Position[3];
                    float attr1 = in1.data[i] / w1;
                    float attr2 = in2.data[i] / w2;
                    float interp_a = attr1 * (1 - t_clip) + attr2 * t_clip;
                    float recip_w = (1 - t_clip) / w1 + t_clip / w2;
                    result.data[i] = interp_a / recip_w;
                }
                break;
            default:
                // Default to linear interpolation.
                result.data[i] = in1.data[i] * (1 - t_clip) + in2.data[i] * t_clip;
                break;
        }
    }
    return result;
}


// -----------------------------------------------------------------------------
// clip_triangle: Recursively clip a triangle against the canonical view volume.
// -----------------------------------------------------------------------------
/**
 * clip_triangle
 *
 * Recursively clips a triangle (defined by vertices v0, v1, and v2) against the
 * six clipping planes of the canonical view volume.
 *
 * When an edge crosses a clip plane the new vertex is computed using
 * perspective_interpolate(). (For attributes flagged as noperspective the
 * interpolation factor is computed in NDC.)
 *
 * When face == 6 (i.e. all clip planes have been processed) the triangle is sent
 * to rasterize_triangle().
 *
 * NOTE: The clip plane order has been changed so that the near and far planes are
 * processed first. This ensures that the proper homogeneous coordinates are used
 * for perspective–correct interpolation.
 */
void clip_triangle(driver_state& state, const data_geometry& v0,
                   const data_geometry& v1, const data_geometry& v2, int face)
{
    if (face == 6) {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }

    // --- MODIFIED CLIP PLANE ORDER ---
    // New order: 0: Near, 1: Far, 2: Left, 3: Right, 4: Bottom, 5: Top.
    vec4 planes[6] = {
        vec4(0, 0, 1, 1),    // Near:   z + w >= 0
        vec4(0, 0, -1, 1),   // Far:   -z + w >= 0
        vec4(1, 0, 0, 1),    // Left:   x + w >= 0
        vec4(-1, 0, 0, 1),   // Right: -x + w >= 0
        vec4(0, 1, 0, 1),    // Bottom: y + w >= 0
        vec4(0, -1, 0, 1)    // Top:    -y + w >= 0
    };
    // --------------------------------------------------------------------------

    vec4 plane = planes[face];

    std::vector<data_geometry> inside, outside;
    std::vector<float> inside_d, outside_d;
    const data_geometry* vertices[3] = { &v0, &v1, &v2 };
    for (int i = 0; i < 3; i++) {
        float d = dot(vertices[i]->gl_Position, plane);
        if (d >= 0) {
            inside.push_back(*vertices[i]);
            inside_d.push_back(d);
        } else {
            outside.push_back(*vertices[i]);
            outside_d.push_back(d);
        }
    }

    if (inside.size() == 3) {
        // All vertices are inside; proceed with clipping against the next plane.
        clip_triangle(state, v0, v1, v2, face + 1);
    }
    else if (inside.size() == 2 && outside.size() == 1) {
        // TWO VERTICES INSIDE; the proper re‐triangulation must be performed.
        // Let:
        //    inside[0] = v0, inside[1] = v1, and outside[0] = v2.
        // Compute the intersections on the edges (v0,v2) and (v1,v2).
        data_geometry i0 = perspective_interpolate(inside[0], outside[0],
                                                   inside_d[0], outside_d[0],
                                                   state.interp_rules, face);
        data_geometry i1 = perspective_interpolate(inside[1], outside[0],
                                                   inside_d[1], outside_d[0],
                                                   state.interp_rules, face);
        // The clipped polygon is now a quadrilateral with vertices:
        //    [inside[0], inside[1], i1, i0]
        // We triangulate it into two triangles:
        //    Triangle 1: (inside[0], inside[1], i1)
        //    Triangle 2: (inside[0], i1, i0)
        clip_triangle(state, inside[0], inside[1], i1, face + 1);
        clip_triangle(state, inside[0], i1, i0, face + 1);
    }
    else if (inside.size() == 1 && outside.size() == 2) {
        // ONE VERTEX INSIDE; compute intersections with both outside vertices.
        data_geometry new_v1 = perspective_interpolate(inside[0], outside[0],
                                                       inside_d[0], outside_d[0],
                                                       state.interp_rules, face);
        data_geometry new_v2 = perspective_interpolate(inside[0], outside[1],
                                                       inside_d[0], outside_d[1],
                                                       state.interp_rules, face);
        clip_triangle(state, inside[0], new_v1, new_v2, face + 1);
    }
}


// -----------------------------------------------------------------------------
// interpolate: Simple linear interpolation between two geometry vertices.
// -----------------------------------------------------------------------------
data_geometry interpolate(const data_geometry& a, const data_geometry& b, float t) {
    data_geometry result;
    result.gl_Position = a.gl_Position * (1 - t) + b.gl_Position * t;
    result.data = new float[MAX_FLOATS_PER_VERTEX];
    for (int i = 0; i < MAX_FLOATS_PER_VERTEX; i++) {
        result.data[i] = a.data[i] * (1 - t) + b.data[i] * t;
    }
    return result;
}

// -----------------------------------------------------------------------------
// to_screen_space: Convert clip-space coordinates to screen-space by performing perspective division.
// -----------------------------------------------------------------------------
vec3 to_screen_space(const vec4& ndc, int width, int height) {
    return vec3(
        (ndc[0] / ndc[3] * 0.5f + 0.5f) * width,
        (ndc[1] / ndc[3] * 0.5f + 0.5f) * height,
        ndc[2] / ndc[3]
    );
}

// -----------------------------------------------------------------------------
// compute_barycentric: Compute barycentric coordinates for point P relative to triangle ABC.
// -----------------------------------------------------------------------------
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

// -----------------------------------------------------------------------------
// rasterize_triangle: Rasterize a triangle defined by three geometry vertices.
// -----------------------------------------------------------------------------
void rasterize_triangle(driver_state& state, const data_geometry& v0,
                        const data_geometry& v1, const data_geometry& v2)
{
    auto ndc_to_screen = [&](const vec4& v) -> vec3 {
        return to_screen_space(v, state.image_width, state.image_height);
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
    if (area == 0) return;
    for (int y = min_y; y <= max_y; y++) {
        for (int x = min_x; x <= max_x; x++) {
            vec3 p(x + 0.5f, y + 0.5f, 0);
            float w0 = edge_function(p1, p2, p) / area;
            float w1 = edge_function(p2, p0, p) / area;
            float w2 = edge_function(p0, p1, p) / area;
            if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
                int pixel_index = y * state.image_width + x;
                float z0 = v0.gl_Position[2] / v0.gl_Position[3];
                float z1 = v1.gl_Position[2] / v1.gl_Position[3];
                float z2 = v2.gl_Position[2] / v2.gl_Position[3];
                float depth = w0 * z0 + w1 * z1 + w2 * z2;
                if (depth < state.image_depth[pixel_index]) {
                    state.image_depth[pixel_index] = depth;
                    data_fragment fragment;
                    fragment.data = new float[state.floats_per_vertex];
                    for (int i = 0; i < state.floats_per_vertex; i++) {
                        if (state.interp_rules[i] == interp_type::smooth) {
                            float w_div_z = w0 / v0.gl_Position[3] + w1 / v1.gl_Position[3] + w2 / v2.gl_Position[3];
                            fragment.data[i] = (w0 * v0.data[i] / v0.gl_Position[3] +
                                                w1 * v1.data[i] / v1.gl_Position[3] +
                                                w2 * v2.data[i] / v2.gl_Position[3]) / w_div_z;
                        } else if (state.interp_rules[i] == interp_type::noperspective) {
                            fragment.data[i] = w0 * v0.data[i] + w1 * v1.data[i] + w2 * v2.data[i];
                        } else if (state.interp_rules[i] == interp_type::flat) {
                            fragment.data[i] = v0.data[i];
                        }
                    }
                    data_output output;
                    state.fragment_shader(fragment, output, state.uniform_data);
                    int r = std::min(255, static_cast<int>(output.output_color[0] * 255));
                    int g = std::min(255, static_cast<int>(output.output_color[1] * 255));
                    int b = std::min(255, static_cast<int>(output.output_color[2] * 255));
                    state.image_color[pixel_index] = make_pixel(r, g, b);
                    delete[] fragment.data;
                }
            }
        }
    }
}