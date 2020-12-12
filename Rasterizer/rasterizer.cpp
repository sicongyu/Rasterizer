#include "rasterizer.h"
#include <memory>
#include <cmath>

Rasterizer::Rasterizer(int width, int height, const char* path) :
	width(width), height(height), model(path) {
	framebuffer = new TGAColor[width * height];
	memset(framebuffer, 0, sizeof(TGAColor) * width * height);
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
		glm::vec3 v[3];
		float y_max = -INFINITY;
		for (int j = 0; j < 3; j++) {
			glm::vec3 vtx = (model.vertices[3 * i + j] - model.min) / (model.max - model.min);
			vtx.x *= width;
			vtx.y *= height;
			if (vtx.y > y_max)
				y_max = vtx.y;
			v[j] = vtx;
		}
		
	}
}