#include "octree.h"
#include "rasterizer.h"
#include <numeric>

Model* OctreeNode::_model = nullptr;

std::vector<AABB> OctreeNode::_aabbs;

OctreeNode::OctreeNode(const std::list<int>& faceIDs, const glm::vec3& center, float length, int depth) : 
	_faceIDs(faceIDs), _center(center), _length(length)
{
	_isLeaf = true;
	if (_faceIDs.size() <= 2 || depth == MAX_DEPTH) {
		return;
	}

	std::list<int> subLists [8];
	for (auto iter = _faceIDs.begin(); iter != _faceIDs.end();) {
		auto aabb = this->_aabbs[*iter];
		auto bitmap = [this](glm::vec3 corner) { return (0x4 & (corner.x < _center.x) << 2) | (0x2 & (corner.y < _center.y) << 1) | (0x1 & (corner.z < _center.z)); };
		// bounding box doesn't intersect with any partition plane
		int map;
		if (bitmap(aabb.min) == (map = bitmap(aabb.max))) {
			// 1, 3; 5, 7
			// 0, 2; 4, 6
			subLists[map].push_back(*iter);
			iter = _faceIDs.erase(iter);
			_isLeaf = false;
		}
		else {
			iter++;
		}
	}

	if (_isLeaf = false) {
		auto signedLength = [this](bool positive) {
			return positive ? _length / 4 : -_length / 4;
		};
		auto center = [this, signedLength] (int map) {
			return glm::vec3(
				_center.x + signedLength((map & 0x4) >> 2),
				_center.y + signedLength((map & 0x2) >> 1),
				_center.z + signedLength(map & 0x1)); 
		};
		for (int i = 0; i < 8; i++) {
			_subNodes[i] = new OctreeNode(subLists[i], center(i), _length / 2, depth + 1);
		}
	}
}

OctreeNode::~OctreeNode()
{
	for (auto subNodePtr : _subNodes) {
		if (subNodePtr) {
			delete subNodePtr;
			subNodePtr = nullptr;
		}
	}
}

Octree::Octree(Model* model, OctreeZBuffer& rasterizer)
{
	std::list<int> faceIDs(model->num_faces);
	std::iota(std::begin(faceIDs), std::end(faceIDs), 0);
	// do view/proj transformation for all vertices HERE
	// the bounding box for the entire mesh
	glm::vec3 model_max = glm::vec3(-FLT_MAX);
	glm::vec3 model_min = glm::vec3(FLT_MAX);
	for (int i = 0; i < model->num_faces; i++) {
		// the bounding box for each triangle
		glm::vec3 max = glm::vec3(-FLT_MAX);
		glm::vec3 min = glm::vec3(FLT_MAX);
		for (int j = 0; j < 3; j++) {
			glm::vec3 vtx = model->vertices[3 * i + j];
			vtx = glm::vec3(rasterizer.ortho * rasterizer.view * glm::vec4(vtx, 1.0f));
			vtx.x = (vtx.x + 1) * rasterizer.width / 2;
			vtx.y = (vtx.y + 1) * rasterizer.height / 2;
			rasterizer._verticesScreenspace.push_back(vtx);
			max = glm::max(max, vtx);
			min = glm::min(min, vtx);
			model_max = glm::max(model_max, max);
			model_min = glm::min(model_min, min);
		}
		OctreeNode::_aabbs.push_back(AABB{ max, min });
	}
	float length = box_max(model_max.x - model_min.x, box_max(model_max.y - model_min.y, model_max.z - model_min.z));
	_root = new OctreeNode(faceIDs, (model_max + model_min) / 2.0f, length, 1);
	OctreeNode::_model = model;
}

Octree::~Octree()
{

}