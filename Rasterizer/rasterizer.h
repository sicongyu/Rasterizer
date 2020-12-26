#pragma once

#include "objloader.hpp"
#include "scanlines.h"
#include "tgaimage.h"
#include <glm/glm.hpp>
#include <vector>
#include <unordered_map>

#define HIERACHY_ZBUFFER 1

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

	std::vector<std::unordered_map<int, PolygonEntry>> polygon_table;

	// an polygon could have multiple edges indexed by the same ID
	std::vector<std::unordered_multimap<int, EdgeEntry>> edge_table;

	std::unordered_map<int, PolygonEntry> active_polygon_table;

	std::unordered_map<int, ActiveEdgeEntry> active_edge_table;
};

class ScanLine: public Rasterizer{
public:

	ScanLine(int width, int height, const char* path);

	~ScanLine();

	void draw();

private:
#if !HIERACHY_ZBUFFER
	float* z_buffer;
#else
	float** z_buffer;

	int _lod;
#endif
};

//class HiZScanLine : public Rasterizer {
//public:
//
//	HiZScanLine(int width, int height, const char* path);
//
//	HiZScanLine();
//
//	void draw();
//
//private:
//
//	float** z_buffer;
//};