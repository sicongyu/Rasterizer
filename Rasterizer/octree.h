#pragma once
#include "objloader.hpp"
#include "scanlines.h"
#include <list>
#include <glm/glm.hpp>

#define MAX_DEPTH 10
#define FLT_MAX 3.402823466e+38F

class OctreeZBuffer;

struct AABB {
	glm::vec3 max;
	glm::vec3 min;
};

class OctreeNode
{
public:
	// Methods
	OctreeNode(const std::list<int>& faceIDs, const glm::vec3& center, float length, int depth);

	~OctreeNode();

	friend class Octree;

	friend class OctreeZBuffer;

private:
	// Fields
	OctreeNode* _subNodes[8];

	std::list<int> _faceIDs;

	bool _isLeaf;

	glm::vec3 _center;

	float _length;

	static Model* _model;

	static std::vector<AABB> _aabbs;
};

class Octree
{
public:
	// Methods
	Octree(Model* model, OctreeZBuffer& rasterizer);

	~Octree();

	OctreeNode* Root() {
		return _root;
	}

private:
	// Fields
	OctreeNode* _root;
};