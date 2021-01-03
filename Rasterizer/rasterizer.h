#pragma once

#include "objloader.hpp"
#include "scanlines.h"
#include "tgaimage.h"
#include "octree.h"
#include <glm/glm.hpp>
#include <vector>
#include <unordered_map>

#define HIERACHY_ZBUFFER 1
#define OCTREE 1
#define VISUALIZE_OCTREE 0

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

	TGAImage* framebuffer;

	glm::vec3 light_dir;

	glm::mat4 view;

	glm::mat4 ortho;

	TGAColor white = TGAColor(255, 255, 255, 255);
};

class ScanLine: public Rasterizer{
public:

	ScanLine(int width, int height, const char* path);

	~ScanLine();

	virtual void draw();

protected:
#if HIERACHY_ZBUFFER || OCTREE
	float** z_buffer;

	int _lod;

	bool TraverseZBuffer(float targetZ, int mip_x, int mip_y, int lod, int xmin, int xmax, int ymin, int ymax);

	bool UpdateZBuffer(float z, int x, int y);
#else
	float* z_buffer;
#endif

private:
	std::vector<std::unordered_map<int, PolygonEntry>> polygon_table;

	// an polygon could have multiple edges indexed by the same ID
	std::vector<std::unordered_multimap<int, EdgeEntry>> edge_table;

	std::unordered_map<int, PolygonEntry> active_polygon_table;

	std::unordered_map<int, ActiveEdgeEntry> active_edge_table;
};

class OctreeZBuffer : public ScanLine {
public:
	OctreeZBuffer(int width, int height, const char* path);

	~OctreeZBuffer();

	void draw();

	friend class Octree;

private:
	std::vector<glm::vec3> _verticesScreenspace;

	Octree _octree;
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