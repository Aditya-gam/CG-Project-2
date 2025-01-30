#include "driver_state.h"
#include <cstring>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <array>

driver_state::driver_state() {}

driver_state::~driver_state() {
    delete[] image_color;
    delete[] image_depth;
}

// Initialize color and depth buffers
void initialize_render(driver_state& state, int width, int height) {
    state.image_width = width;
    state.image_height = height;

    // Allocate memory for color buffer
    state.image_color = new pixel[width * height];
    std::fill_n(state.image_color, width * height, make_pixel(0, 0, 0)); // Initialize to black

    // Allocate memory for depth buffer
    state.image_depth = new float[width * height];
    std::fill_n(state.image_depth, width * height, 1.0f); // Initialize depth buffer to farthest
}

// Render function to process vertices and invoke rasterization
void render(driver_state& state, render_type type) {
    if (type != render_type::triangle) {
        std::cerr << "Only render_type::triangle is implemented." << std::endl;
        return;
    }

    // Iterate over each triangle
    for (int i = 0; i < state.num_vertices; i += 3) {
        data_geometry v[3];

        for (int j = 0; j < 3; j++) {
            data_vertex vertex;
            vertex.data = &state.vertex_data[(i + j) * state.floats_per_vertex];

            // Invoke vertex shader to transform the vertex
            state.vertex_shader(vertex, v[j], state.uniform_data);
        }

        // Clip the triangle
        clip_triangle(state, v[0], v[1], v[2], 0);
    }
}

// Basic clipping function
void clip_triangle(driver_state& state, const data_geometry& v0,
                   const data_geometry& v1, const data_geometry& v2, int face) {
    if (face == 6) {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }

    // For now, bypass clipping and simply proceed to the next face
    clip_triangle(state, v0, v1, v2, face + 1);
}

// Barycentric coordinate method for rasterization
void rasterize_triangle(driver_state& state, const data_geometry& v0,
                        const data_geometry& v1, const data_geometry& v2) {
    auto ndc_to_screen = [&](vec4 v) -> vec2 {
        return vec2((v[0] * 0.5f + 0.5f) * state.image_width,
                    (v[1] * 0.5f + 0.5f) * state.image_height);
    };

    vec2 p0 = ndc_to_screen(v0.gl_Position);
    vec2 p1 = ndc_to_screen(v1.gl_Position);
    vec2 p2 = ndc_to_screen(v2.gl_Position);

    std::array<float, 3> x_values = {p0[0], p1[0], p2[0]};
    std::array<float, 3> y_values = {p0[1], p1[1], p2[1]};

    int min_x = std::max(0, static_cast<int>(std::floor(*std::min_element(x_values.begin(), x_values.end()))));
    int max_x = std::min(state.image_width - 1, static_cast<int>(std::ceil(*std::max_element(x_values.begin(), x_values.end()))));
    int min_y = std::max(0, static_cast<int>(std::floor(*std::min_element(y_values.begin(), y_values.end()))));
    int max_y = std::min(state.image_height - 1, static_cast<int>(std::ceil(*std::max_element(y_values.begin(), y_values.end()))));

    auto edge_function = [](const vec2& a, const vec2& b, const vec2& c) -> float {
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    };

    float area = edge_function(p0, p1, p2);
    if (area == 0) return; // Skip degenerate triangles

    for (int y = min_y; y <= max_y; y++) {
        for (int x = min_x; x <= max_x; x++) {
            vec2 p = vec2(x + 0.5f, y + 0.5f);

            float w0 = edge_function(p1, p2, p) / area;
            float w1 = edge_function(p2, p0, p) / area;
            float w2 = edge_function(p0, p1, p) / area;

            if (w0 >= 0 && w1 >= 0 && w2 >= 0) {
                int pixel_index = y * state.image_width + x;

                // **Corrected Depth Interpolation**
                float z0 = v0.gl_Position[2] / v0.gl_Position[3];
                float z1 = v1.gl_Position[2] / v1.gl_Position[3];
                float z2 = v2.gl_Position[2] / v2.gl_Position[3];

                float depth = w0 * z0 + w1 * z1 + w2 * z2;

                if (depth < state.image_depth[pixel_index]) {
                    state.image_depth[pixel_index] = depth;

                    // Allocate memory for interpolated attributes
                    data_fragment fragment;
                    fragment.data = new float[MAX_FLOATS_PER_VERTEX]();

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
