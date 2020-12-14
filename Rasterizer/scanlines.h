#pragma once
#include <glm/glm.hpp>
#include "tgaimage.h"

struct PolygonEntry {
	glm::vec3 plane;
	int faceID;
	int dy;
	TGAColor color;
};

struct EdgeEntry {
	float x_at_ymax;
	float z_at_ymax;
	float dx;
	int faceID;
	int dy;
};

struct ActiveEdgeEntry {
	float xl;
	float dxl;
	int dyl;
	float xr;
	float dxr;
	int dyr;
	float zl;
	float dzx;
	float dzy;
	int id;
};