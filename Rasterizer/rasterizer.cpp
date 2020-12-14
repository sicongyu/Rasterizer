#include "rasterizer.h"
#include <memory>
#include <cmath>

Rasterizer::Rasterizer(int width, int height, const char* path) :
	width(width), height(height), model(path) {
	framebuffer = new TGAColor[width * height];
	memset(framebuffer, 0, sizeof(TGAColor) * width * height);
	light_dir = glm::vec3(0.0f, 0.0f, -1.0f);
}

Rasterizer::~Rasterizer() {
	delete[] framebuffer;
}

ScanLine::ScanLine(int width, int height, const char* path) :
	Rasterizer(width, height, path) {
	polygon_table = new std::vector<PolygonEntry>[height];
	edge_table = new std::vector<EdgeEntry>[height];
	z_buffer = new float[height * width];
	memset(z_buffer, 0, sizeof(float) * width * height);
}

ScanLine::~ScanLine() {
	delete[] polygon_table;
	delete[] edge_table;
}

void ScanLine::draw() {

	// create polygon_table and edge_table
	for (int i = 0; i < model.num_faces; i++) {
		// for every vertices in the triangle
		glm::vec3 vertices[3];
		// map for vertices: ymax-0, ymid-1, ymin-2
		int y_map[3];
		float y_max = -INFINITY;
		float y_min = INFINITY;
		for (int j = 0; j < 3; j++) {
			// convert into screen space
			glm::vec3 vtx = (model.vertices[3 * i + j] - model.min) / (model.max - model.min);
			vtx.x *= width;
			vtx.y *= height;

			if (vtx.y > y_max) {
				y_max = vtx.y;
				y_map[0] = j;
			}
			else if (vtx.y < y_min) {
				y_min = vtx.y;
				y_map[2] = j;
			}
			vertices[j] = vtx;
		}
		y_map[1] = 3 - y_map[0] - y_map[2];

		// put entry into polygon_table
		glm::vec3 n = model.normals[i];
		float intensity = glm::dot(n, light_dir);
		PolygonEntry polygon_entry = { model.normals[i], i, (int)y_max - (int)y_min, intensity * white};
		polygon_table[int(y_max)].push_back(polygon_entry);

		// combination of (0, 1), (0, 2), (1, 2), put entries into edge_table
		for (int m = 0; m < 2; m++) {
			for (int n = m + 1; n < 3; n++) {
				glm::vec3 vtx0 = vertices[y_map[m]];
				glm::vec3 vtx1 = vertices[y_map[n]];
				glm::vec3 edge = vtx1 - vtx0;
				EdgeEntry edge_entry = { vtx0.x, //x of upper vertex
								  vtx0.z, // z of upper vertex
								  -(edge.x / edge.y), // dx (edge.y is negative)
								  i,
								  (int)vtx0.y - (int)vtx0.y,
								  };
				edge_table[int(vtx0.y)].push_back(edge_entry);
			}
		}
	}

	// scan process: from top to bottom
	for (int i = height; i >= 0; i--) {
		// add new active polygons and edges
		if (!polygon_table[i].empty()) {
			for (auto polygon : polygon_table[i]) {
				active_polygon_table.push_back(polygon);
				EdgeEntry* edge0;
				EdgeEntry* edge1;
				for (auto edge : edge_table[i]) {
					if (edge.faceID == polygon.faceID) {
						if (!edge0) {
							edge0 = &edge;
						}
						else if (!edge1) {
							edge1 = &edge;
						}
					}
				}
				// ensure that edge0 is the leftside edge
				if (edge0->dx > edge1->dx) {
					std::swap(edge0, edge1);
				}
				float c = polygon.plane.z == 0 ? 0.0001 : polygon.plane.z;
				ActiveEdgeEntry active_edge_entry = {
					edge0->x_at_ymax, edge0->dx, edge0->dy,
					edge1->x_at_ymax, edge1->dx, edge1->dy,
					edge0->z_at_ymax, -polygon.plane.x / c, polygon.plane.y / c };
			}
		}

		// update pixel for each scanline.
			
	}
}