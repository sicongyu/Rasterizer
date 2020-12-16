#pragma once

#include "objloader.hpp"
#include "scanlines.h"
#include "tgaimage.h"
#include <vector>
#include <unordered_map>

class Rasterizer{

public:

	Rasterizer(int width, int height, const char* path);

	~Rasterizer();

	virtual void draw() = 0;

	TGAImage* getFramebuffer() { return framebuffer; }

protected:

	int width;

	int height;

	Model model;

	glm::vec3 light_dir;

	TGAImage* framebuffer;
};

class ScanLine: public Rasterizer{
public:

	ScanLine(int width, int height, const char* path);

	~ScanLine();

	void draw();

private:
	std::unordered_map<int, PolygonEntry> * polygon_table;
	// an polygon could have multiple edges indexed by the same ID
	std::unordered_multimap<int, EdgeEntry> * edge_table;
	std::unordered_map<int, PolygonEntry> active_polygon_table;
	std::unordered_map<int, ActiveEdgeEntry> active_edge_table;
	//std::vector<PolygonEntry> * polygon_table;
	//std::vector<EdgeEntry> * edge_table;
	//std::vector<PolygonEntry> active_polygon_table;
	//std::vector<ActiveEdgeEntry> active_edge_table;
	float* z_buffer;
	
};