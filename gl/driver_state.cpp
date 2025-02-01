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
    for (int i = 0; i < width * height; i++)
    {
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
    // (The existing implementation remains unchanged.)
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


/**
 * perspective_interpolate
 *
 * Given two geometry vertices (along an edge crossing the clip boundary),
 * compute the intersection vertex.
 *
 * @param in1         The first vertex (inside the clipping plane).
 * @param in2         The second vertex (outside the clipping plane).
 * @param d1          The signed distance (dot product with the plane) for in1.
 * @param d2          The signed distance for in2.
 * @param interp_rules  The interpolation rules for each attribute.
 *
 * @return The interpolated vertex (in clip space) along the edge.
 */
data_geometry perspective_interpolate(const data_geometry& in1, const data_geometry& in2,
    float d1, float d2, const interp_type interp_rules[MAX_FLOATS_PER_VERTEX])
{
    data_geometry result;

    // Compute the interpolation factor t using the plane distances.
    float t = d1 / (d1 - d2);

    // Linearly interpolate the clip-space position.
    result.gl_Position = in1.gl_Position * (1 - t) + in2.gl_Position * t;

    // Allocate a new attribute array.
    result.data = new float[MAX_FLOATS_PER_VERTEX];

    // Interpolate each vertex attribute.
    for (int i = 0; i < MAX_FLOATS_PER_VERTEX; i++) {
        if (interp_rules[i] == interp_type::noperspective || interp_rules[i] == interp_type::flat) {
            // Simple linear interpolation for non–perspective or flat.
            result.data[i] = in1.data[i] * (1 - t) + in2.data[i] * t;
        } else { // interp_type::smooth: perspective–correct interpolation.
            float w1 = in1.gl_Position[3];
            float w2 = in2.gl_Position[3];
            float attr1 = in1.data[i] / w1;
            float attr2 = in2.data[i] / w2;
            float interp_a = attr1 * (1 - t) + attr2 * t;
            float recip_w = (1 - t) / w1 + t / w2;
            result.data[i] = interp_a / recip_w;
        }
    }

    return result;
}

/**
 * clip_triangle
 *
 * Recursively clips a triangle (specified by vertices v0, v1, and v2)
 * against the six clip planes. When face == 6, the triangle is passed on
 * for rasterization.
 *
 * For each clip plane, we compute the signed distance d = dot(gl_Position, plane)
 * for each vertex.
 *
 * @param state  The driver state.
 * @param v0     The first vertex.
 * @param v1     The second vertex.
 * @param v2     The third vertex.
 * @param face   The current clipping plane (0 to 5).
 */
void clip_triangle(driver_state& state, const data_geometry& v0,
                   const data_geometry& v1, const data_geometry& v2, int face)
{
    if (face == 6) {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }

    std::vector<data_geometry> inside, outside;
    std::vector<float> inside_d, outside_d;

    // Define the six clip planes. Each is represented as (A, B, C, D) so that
    // the plane equation is: A*x + B*y + C*z + D*w >= 0.
    vec4 planes[6] = {
        vec4(1, 0, 0, 1),   // Left:    x + w >= 0
        vec4(-1, 0, 0, 1),  // Right:  -x + w >= 0
        vec4(0, 1, 0, 1),   // Bottom:  y + w >= 0
        vec4(0, -1, 0, 1),  // Top:    -y + w >= 0
        vec4(0, 0, 1, 1),   // Near:   z + w >= 0
        vec4(0, 0, -1, 1)   // Far:   -z + w >= 0
    };

    vec4 plane = planes[face];

    // Compute signed distances for each vertex.
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
        // All vertices are inside; clip against the next plane.
        clip_triangle(state, v0, v1, v2, face + 1);
    }
    else if (inside.size() == 2 && outside.size() == 1) {
        // Two vertices inside; compute intersections along the two edges.
        data_geometry new_v1 = perspective_interpolate(inside[0], outside[0], inside_d[0], outside_d[0], state.interp_rules);
        data_geometry new_v2 = perspective_interpolate(inside[1], outside[0], inside_d[1], outside_d[0], state.interp_rules);

        // Form two new triangles.
        clip_triangle(state, inside[0], inside[1], new_v1, face + 1);
        clip_triangle(state, inside[1], new_v1, new_v2, face + 1);
    }
    else if (inside.size() == 1 && outside.size() == 2) {
        // One vertex inside; compute intersections with both outside vertices.
        data_geometry new_v1 = perspective_interpolate(inside[0], outside[0], inside_d[0], outside_d[0], state.interp_rules);
        data_geometry new_v2 = perspective_interpolate(inside[0], outside[1], inside_d[0], outside_d[1], state.interp_rules);

        // The clipped triangle consists of the inside vertex and the two intersection points.
        clip_triangle(state, inside[0], new_v1, new_v2, face + 1);
    }
}

/**
 * Simple linear interpolation between two geometry vertices.
 */
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
// to_screen_space: Convert clip-space coordinates to screen-space, performing
// perspective division.
// -----------------------------------------------------------------------------

vec3 to_screen_space(const vec4& ndc, int width, int height) {
    return vec3(
        (ndc[0] / ndc[3] * 0.5f + 0.5f) * width,
        (ndc[1] / ndc[3] * 0.5f + 0.5f) * height,
        ndc[2] / ndc[3]
    );
}

// Compute barycentric coordinates (used in rasterization).
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
// This function performs barycentric interpolation, calls the fragment shader,
// and performs z-buffering.
// -----------------------------------------------------------------------------

void rasterize_triangle(driver_state& state, const data_geometry& v0,
                        const data_geometry& v1, const data_geometry& v2)
{
    // --- FIX: Use proper perspective division for screen-space conversion ---
    // Instead of using a lambda that does not divide by w, we call to_screen_space,
    // which performs the perspective divide. This change fixes test 05.txt while
    // leaving test 10.txt (which already has w = 1) unchanged.
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
    if (area == 0) return; // Skip degenerate triangles

    for (int y = min_y; y <= max_y; y++) {
        for (int x = min_x; x <= max_x; x++) {
            vec3 p(x + 0.5f, y + 0.5f, 0);

            float w0 = edge_function(p1, p2, p) / area;
            float w1 = edge_function(p2, p0, p) / area;
            float w2 = edge_function(p0, p1, p) / area;

            if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
                int pixel_index = y * state.image_width + x;

                // Perspective-correct depth interpolation.
                float z0 = v0.gl_Position[2] / v0.gl_Position[3];
                float z1 = v1.gl_Position[2] / v1.gl_Position[3];
                float z2 = v2.gl_Position[2] / v2.gl_Position[3];
                float depth = w0 * z0 + w1 * z1 + w2 * z2;

                if (depth < state.image_depth[pixel_index]) {
                    state.image_depth[pixel_index] = depth;

                    // Allocate memory for interpolated attributes.
                    data_fragment fragment;
                    fragment.data = new float[state.floats_per_vertex];

                    // Interpolate attributes according to the interpolation rules.
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

                    // Process the fragment through the fragment shader.
                    data_output output;
                    state.fragment_shader(fragment, output, state.uniform_data);

                    int r = std::min(255, static_cast<int>(output.output_color[0] * 255));
                    int g = std::min(255, static_cast<int>(output.output_color[1] * 255));
                    int b = std::min(255, static_cast<int>(output.output_color[2] * 255));

                    state.image_color[pixel_index] = make_pixel(r, g, b);

                    // Cleanup allocated memory.
                    delete[] fragment.data;
                }
            }
        }
    }
}
