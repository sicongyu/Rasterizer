#pragma once

#include "objloader.hpp"
#include "scanlines.h"
#include "tgaimage.h"
#include <vector>

class Rasterizer{

public:

	Rasterizer(int width, int height, const char* path);

	~Rasterizer();

	virtual void draw() = 0;

protected:

	int width;

	int height;

	Model model;

	glm::vec3 light_dir;

	TGAColor* framebuffer;
};

class ScanLine: public Rasterizer{
public:

	ScanLine(int width, int height, const char* path);

	~ScanLine();

	void draw();

private:
	std::vector<PolygonEntry> * polygon_table;
	std::vector<EdgeEntry> * edge_table;
	std::vector<PolygonEntry> active_polygon_table;
	std::vector<ActiveEdgeEntry> active_edge_table;
	float* z_buffer;
	
};