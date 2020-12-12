#pragma once
#include <glm/glm.hpp>

struct PolygonEntry {
	glm::vec4 plane;
	int faceID;
	int dy;
	glm::vec3 color;
};

struct EdgeEntry {
	float x_at_ymax;
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