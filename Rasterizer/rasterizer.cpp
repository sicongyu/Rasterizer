#include "rasterizer.h"
#include <memory>
#include <cmath>
#include <functional>
#include "timer.h"
#include <glm/gtc/matrix_transform.hpp>

Rasterizer::Rasterizer(int width, int height, const char* path) :
	width(width), height(height), model(path) {
	framebuffer = new TGAImage(width, height, TGAImage::RGB);
	light_dir = glm::vec3(0.0f, 0.0f, 1.0f);
	// the lookat matrix for cerberus
	glm::vec3 cameraFront = glm::vec3(0.2f, 0.0f, -1.0f);
	glm::vec3 eye = glm::vec3(0.0, 0.0, 3.0);
	view = glm::lookAt(
		eye,
		eye + cameraFront,
		glm::vec3(0.0, 1.0, 0.0) // up
	);
	ortho = glm::ortho(-1.2f, 0.7f, -1.0f, 1.0f);
}

Rasterizer::~Rasterizer() {
	delete framebuffer;
}

ScanLine::ScanLine(int width, int height, const char* path) :
	Rasterizer(width, height, path) {
	polygon_table = std::vector<std::unordered_map<int, PolygonEntry>>(height);
	edge_table = std::vector<std::unordered_multimap<int, EdgeEntry>>(height);
#if !HIERACHY_ZBUFFER
	z_buffer = new float[height * width];
	std::fill_n(z_buffer, height * width, -FLT_MAX);
#else
	// Hierachy z-buffer: the smaller lod is, the blurer the picture
	_lod = log2(width) + 1; // Assume height == width
	z_buffer = new float*[_lod];
	for (int i = 0; i < _lod; i++) {
		int buffer_length = pow(4, i);
		z_buffer[i] = new float[buffer_length];
		std::fill_n(z_buffer[i], buffer_length, -FLT_MAX);
	}
#endif
}

ScanLine::~ScanLine() {
#if !HIERACHY_ZBUFFER
	delete[] z_buffer;
#else
	for (int i = 0; i < _lod; i++) {
		delete[] z_buffer[i];
	}
	delete[] z_buffer;
#endif // !HIERACHY_ZBUFFER

}

//HiZScanLine::HiZScanLine(int width, int height, const char* path) :
//	Rasterizer(width, height, path) {
//	z_buffer = new float*[HIERACHY_ZBUFFER];
//}

void ScanLine::draw() {
	TIMING_BEGIN("Start Rasterization")
	// create polygon_table and edge_table
	for (int i = 0; i < model.num_faces; i++) {
		// for every vertices in the triangle
		glm::vec3 vertices[3];
		// map for vertices: ymax-0, ymid-1, ymin-2
		int y_map[3];
		float y_max = -FLT_MAX;
		float y_min = FLT_MAX;
		float x_max = -FLT_MAX;
		float x_min = FLT_MAX;
		float z_max = -FLT_MAX;
		glm::vec2 abs_max = glm::vec2(box_max(abs(model.max.x), abs(model.min.x)), box_max(abs(model.max.y), abs(model.min.y)));
		glm::vec3 world_normal = model.normals[3 * i]; // every three vtxs share the same normal;
		for (int j = 0; j < 3; j++) {
			// convert into screen space
			//glm::vec3 vtx = model.vertices[3 * i + j];
			//glm::vec2 vtx2d = (glm::vec2(vtx.x, vtx.y) / (1.5f * abs_max) + glm::vec2(1.0f)) * glm::vec2(width / 2, height / 2);
			//vtx = glm::vec3(vtx2d.x, vtx2d.y, vtx.z);

			glm::vec3 vtx = glm::vec3(ortho * view * glm::vec4(model.vertices[3 * i + j], 1.0));
			vtx.x = (vtx.x + 1) * width / 2;
			vtx.y = (vtx.y + 1) * height / 2;

			if (vtx.y > y_max) {
				y_max = vtx.y;
				y_map[0] = j;
			}
			if (vtx.y < y_min) {
				y_min = vtx.y;
				y_map[2] = j;
			}
			x_max = box_max(vtx.x, x_max);
			x_min = box_min(vtx.x, x_min);
			z_max = box_max(vtx.z, z_max);
			vertices[j] = vtx;
		}

		int polygon_dy = (int)y_max - (int)y_min + 1;

		// Occasionally when three vertices share the same y value, i.e., in the xoz plane
		if (y_map[0] == y_map[2]) {
			// assign an oreder manually
			y_map[2] = 2;
		}
		y_map[1] = 3 - y_map[0] - y_map[2];
		if (!(y_map[1] >= 0 && y_map[1] <= 2)) {
			throw "Error: invalid y_map index";
		}

		float intensity = box_max(glm::dot(world_normal, -light_dir), 0.0f);

		glm::vec3 edges[3];

		// combination of (0, 1), (0, 2), (1, 2), put entries into edge_table
		for (int m = 0; m < 2; m++) {
			for (int n = m + 1; n < 3; n++) {
				glm::vec3 vtx0 = vertices[y_map[m]];
				glm::vec3 vtx1 = vertices[y_map[n]];

				glm::vec3 edge = vtx1 - vtx0;
				edges[m + n - 1] = edge;

				int edge_dy = (int)vtx0.y - (int)vtx1.y + 1;
				if (edge_dy == 1 && m + n == 1) {
					continue; // discard edges that covers only one row
				}

				EdgeEntry edge_entry = { 
					vtx0.x, //x of upper vertex
					vtx0.z, // z of upper vertex
					-(edge.x / edge.y), // dx (edge.y is negative)
					i, // faceID
					edge_dy // dy, notice that the height of vtx0 must be greater than vtx1
				};
				//edge_table[int(vtx0.y)].push_back(edge_entry);
				if (int(vtx0.y) < 0 || int(vtx0.y) >= height) {
					continue;
				}
				edge_table[int(vtx0.y)].emplace(edge_entry.faceID, edge_entry);
			}
		}

		glm::vec3 screen_normal = glm::normalize(glm::cross(edges[0], edges[1]));
		screen_normal = glm::dot(screen_normal, world_normal) > 0 ? screen_normal : -screen_normal;

		world_normal = (world_normal + glm::vec3(1)) / 2.0f;

		// put entry into polygon_table
		PolygonEntry polygon_entry = {
			screen_normal, // plane
			i, //faceID
			polygon_dy, // dy
			white * intensity, // color
			//white, // debug purpose
			//TGAColor(i * (255.0 / model.num_faces), 0, 0, 255),
			//TGAColor(255 * world_normal.x, 255 * world_normal.y, 255 * world_normal.z, 255),
			x_min, // xl
			x_max, // xr
			z_max // znear
		};
		//polygon_table[int(y_max)].push_back(polygon_entry);

		if (int(y_max) < 0 || int(y_max) >= height) {
			continue;
		}
		polygon_table[int(y_max)].emplace(polygon_entry.faceID, polygon_entry);
	}

	// scan process: from top to bottom
	for (int y = height-1; y >= 0; y--) {
		// add new active polygons and edges
		if (!polygon_table[y].empty()) {
			for (auto polygon_table_iter : polygon_table[y]) {
				int faceID = polygon_table_iter.first;
				auto polygon = polygon_table_iter.second;


#if HIERACHY_ZBUFFER
				// do znear check before rasterization a polygon
				if (!TraverseZBuffer(polygon.znear, 0, 0, 0, polygon.xl, polygon.xr, y - polygon.dy, y)) {
					continue;
				}
#endif // HIERACHY_ZBUFFER

				// add entries to active polygon table
				//active_polygon_table.push_back(polygon);
				active_polygon_table.emplace(faceID, polygon);
				
				if (faceID == 4092) {
					int err = 0;
				}

				// serach for corresponding edges of the polygon: there must to be at least two edges, and we do not consider the three-edge case
				std::vector<EdgeEntry> edges;
				auto cnt = edge_table[y].count(faceID);
				if (!(cnt == 2 || cnt == 3)) {
					throw "Invalid number of edges with the same ID in a single scan";
				}
				auto edge = edge_table[y].find(faceID);
				for (;cnt > 0; cnt--, edge++) {
					edges.push_back(edge->second);
				}

				// ensure that edge0 is the leftside edge: two cases
				// 1. /\ (most cases) 
				// 2. \/ (speical one: the upper edge is horizontal
				if ((edges[0].x_at_ymax == edges[1].x_at_ymax && edges[0].dx > edges[1].dx) || (edges[0].x_at_ymax > edges[1].x_at_ymax)) {
					std::reverse(edges.begin(), edges.end());
				}

				// add entries to active edge table
				//float c = polygon.plane.z == 0 ? 0.0001 : polygon.plane.z;
				float c = polygon.plane.z;
				ActiveEdgeEntry active_edge_entry = {
					edges[0].x_at_ymax, edges[0].dx, edges[0].dy,
					edges[1].x_at_ymax, edges[1].dx, edges[1].dy,
					edges[0].z_at_ymax, -polygon.plane.x / c, polygon.plane.y / c,
					faceID
				};
				active_edge_table.emplace(faceID, active_edge_entry);
				//active_edge_table.push_back(active_edge_entry);
			}
		}

		for (auto active_edge_pair_iter = active_edge_table.begin(); active_edge_pair_iter != active_edge_table.end();) {
			/*--- consider this case will not happen anymore ---*/
			//if (active_edge->second.dyl == 0 && active_edge->second.dyr == 0) {
			//	active_edge_table.erase(active_edge);
			//}
			auto& active_edge = active_edge_pair_iter->second;
			auto& active_polygon = active_polygon_table[active_edge_pair_iter->first];

			if (active_edge.faceID == 4092) {
				int err = 0;
			}

			// update pixels from left to right
			for (int x = active_edge.xl + 0.5f; x <= active_edge.xr + 0.5f; x++) {
				// update buffer contents
				float z = active_edge.zl + (x - active_edge.xl) * active_edge.dzx;
#if !HIERACHY_ZBUFFER
				if (z > z_buffer[y * width + x]) {
					z_buffer[y * width + x] = z;
					framebuffer->set(x, y, active_polygon.color);
				}
#else			
				if (z > z_buffer[_lod - 1][y * width + x]) {
					z_buffer[_lod - 1][y * width + x] = z;
					framebuffer->set(x, y, active_polygon.color);
					// update the hierachy z-buffer
					for (int i = _lod - 2; i >= 0; i--) {
						int mip_x = x / (int)pow(2, _lod - (1 + i));
						int mip_y = y / (int)pow(2, _lod - (1 + i));
						auto& z_in_buffer = z_buffer[i][mip_y * (int)pow(2, i) + mip_x];
						for (int j = 0; j < 4; j++) {
							int sub_x = 2 * mip_x + j % 2;
							int sub_y = 2 * mip_y + j / 2;
							auto& z_in_sub = z_buffer[i + 1][sub_y * (int)pow(2, i + 1) + sub_x];
							z_in_buffer = box_min(z_in_sub, z_in_buffer);
						}
					}
				}
#endif // !HIERACHY_ZBUFFER

				//debug
				if (y == 308 && x == 395) {
					int err = 0;
				}
			}

			{
				active_edge.xl += active_edge.dxl;
				active_edge.xl = box_max(active_polygon.xl, active_edge.xl);
				active_edge.xr += active_edge.dxr;
				active_edge.xr = box_min(active_polygon.xr, active_edge.xr);
				active_edge.zl += active_edge.dzx * active_edge.dxl + active_edge.dzy;
			}

			// update active edge
			active_edge.dyl--;
			active_edge.dyr--;
			auto polygon_dy = --active_polygon.dy;

			if (active_edge.dyl == 0 && active_edge.dyr == 0) {
				if (polygon_dy != 0) {
					throw "Error: Unfinished Polygon Rasterization";
				}
				active_polygon_table.erase(active_edge.faceID);
				active_edge_pair_iter = active_edge_table.erase(active_edge_pair_iter);
				continue;
			}
			else if (active_edge.dyl == 0 || active_edge.dyr == 0) {
				auto new_edge = edge_table[y].find(active_edge.faceID)->second;
				// dzx, dzy, zl and xl won't change
				if (active_edge.dyl == 0) {
					active_edge.dxl = new_edge.dx;
					active_edge.dyl = new_edge.dy - 1;
				}
				else {
					active_edge.dxr = new_edge.dx;
					active_edge.dyr = new_edge.dy - 1;
				}
			}
			++active_edge_pair_iter;
		}
	}
	TIMING_END("Rasterization Done")
}

#if HIERACHY_ZBUFFER || OCTREE
bool ScanLine::TraverseZBuffer(float targetZ, int mip_x, int mip_y, int lod, int xmin, int xmax, int ymin, int ymax) {
	// if the polygon has znear closer than current z: try to draw it!
	if (targetZ > z_buffer[lod][mip_y * (int)pow(2, lod) + mip_x]) {
		// if has reached base z-buffer: draw it!
		if (lod == this->_lod - 1) {
			return true;
		}
		// recursively subdivide the quad
		else {
			for (int i = 0; i < 4; i++) {
				int sub_x = 2 * mip_x + i % 2;
				int sub_y = 2 * mip_y + i / 2;
				int sub_node_length = pow(2, _lod - lod - 2);
				if (sub_x < xmin / sub_node_length || sub_x > xmax / sub_node_length ||
					sub_y < ymin / sub_node_length || sub_y > ymax / sub_node_length) {
					continue;
				}
				if (TraverseZBuffer(targetZ, sub_x, sub_y, lod + 1, xmin, xmax, ymin, ymax)) {
					return true;
				}
			}
		}
	}
	return false;
}

#endif

OctreeZBuffer::OctreeZBuffer(int width, int height, const char* path) : 
	ScanLine(width, height, path),
	_octree(&model, *this){
}

OctreeZBuffer::~OctreeZBuffer() {}

void OctreeZBuffer::draw() {
	TIMING_BEGIN("Start Octree Rasterization")

	// Traverse through the octree
	std::list<OctreeNode*> stage;
	stage.push_back(_octree.Root());
	while (!stage.empty()) {
		OctreeNode* node = stage.back();
		// visibility test of the node
		float halfLength = node->_length / 2;
		if (!TraverseZBuffer(node->_center.z + halfLength, 0, 0, 0,
			(int)(node->_center.x - halfLength),
			(int)(node->_center.x + halfLength),
			(int)(node->_center.y - halfLength),
			(int)(node->_center.y + halfLength))) 
		{
			continue;
		}
		// if there are triangles in this ndoe: rasterize it!
		for (auto i : node->_faceIDs) {
			glm::vec3 world_normal = model.normals[3 * i]; // every three vtxs share the same normal;
			// for every vertices in the triangle
			glm::vec3 vertices[3];
			int y_map[3];
			float y_max = -FLT_MAX;
			float y_min = FLT_MAX;
			float x_max = -FLT_MAX;
			float x_min = FLT_MAX;
			float z_max = -FLT_MAX;
			for (int j = 0; j < 3; j++) {
				// convert into screen space
				//glm::vec3 vtx = model.vertices[3 * i + j];
				//glm::vec2 vtx2d = (glm::vec2(vtx.x, vtx.y) / (1.5f * abs_max) + glm::vec2(1.0f)) * glm::vec2(width / 2, height / 2);
				//vtx = glm::vec3(vtx2d.x, vtx2d.y, vtx.z);

				glm::vec3 vtx = _verticesScreenspace[3 * i + j];
				if (vtx.y > y_max) {
					y_max = vtx.y;
					y_map[0] = j;
				}
				if (vtx.y < y_min) {
					y_min = vtx.y;
					y_map[2] = j;
				}
				x_max = box_max(vtx.x, x_max);
				x_min = box_min(vtx.x, x_min);
				z_max = box_max(vtx.z, z_max);
				vertices[j] = vtx;
			}

			int polygon_dy = (int)y_max - (int)y_min + 1;

			// Occasionally when three vertices share the same y value, i.e., in the xoz plane
			if (y_map[0] == y_map[2]) {
				// assign an oreder manually
				y_map[2] = 2;
			}
			y_map[1] = 3 - y_map[0] - y_map[2];
			if (!(y_map[1] >= 0 && y_map[1] <= 2)) {
				throw "Error: invalid y_map index";
			}

			float intensity = box_max(glm::dot(world_normal, -light_dir), 0.0f);

			glm::vec3 edges[3]; // raw data of edges are kept for computing screen-space normal
			std::vector<EdgeEntry> edgeEntries(3); // Edge entries are kept for creating active edge entry
			for (int m = 0; m < 2; m++) {
				for (int n = m + 1; n < 3; n++) {
					glm::vec3 vtx0 = vertices[y_map[m]];
					glm::vec3 vtx1 = vertices[y_map[n]];
					glm::vec3 edge = vtx1 - vtx0;
					edges[m + n - 1] = edge;

					int edge_dy = (int)vtx0.y - (int)vtx1.y + 1;
					if (edge_dy == 1 && m + n == 1) {
						continue; // discard edges that covers only one row
					}

					EdgeEntry edge_entry = {
						vtx0.x, //x of upper vertex
						vtx0.z, // z of upper vertex
						-(edge.x / edge.y), // dx (edge.y is negative)
						i, // faceID
						edge_dy // dy, notice that the height of vtx0 must be greater than vtx1
					};
					if (int(vtx0.y) < 0 || int(vtx0.y) >= height) {
						continue;
					}
					edgeEntries[m + n - 1] = edge_entry;
				}
			}

			glm::vec3 screen_normal = glm::normalize(glm::cross(edges[0], edges[1]));
			screen_normal = glm::dot(screen_normal, world_normal) > 0 ? screen_normal : -screen_normal;

			//world_normal = (world_normal + glm::vec3(1)) / 2.0f;

			PolygonEntry polygon = {
				screen_normal, // plane
				i, //faceID
				polygon_dy, // dy
				white * intensity, // color
				//white, // debug purpose
				//TGAColor(i * (255.0 / model.num_faces), 0, 0, 255),
				//TGAColor(255 * world_normal.x, 255 * world_normal.y, 255 * world_normal.z, 255),
				x_min, // xl
				x_max, // xr
				z_max // znear
			};

			ActiveEdgeEntry active_edge = {
				edgeEntries[0].x_at_ymax, edgeEntries[0].dx, edgeEntries[0].dy,
				edgeEntries[1].x_at_ymax, edgeEntries[1].dx, edgeEntries[1].dy,
				edgeEntries[0].z_at_ymax, -polygon.plane.x / polygon.plane.z, polygon.plane.y / polygon.plane.z,
				i
			};

			for (int y = y_max; y > y_max - polygon_dy; y--) {
				// update pixels from left to right
				for (int x = active_edge.xl + 0.5f; x <= active_edge.xr + 0.5f; x++) {
					// update buffer contents
					float z = active_edge.zl + (x - active_edge.xl) * active_edge.dzx;	
					if (z > z_buffer[_lod - 1][y * width + x]) {
						z_buffer[_lod - 1][y * width + x] = z;
						framebuffer->set(x, y, polygon.color);
						// update the hierachy z-buffer
						for (int i = _lod - 2; i >= 0; i--) {
							int mip_x = x / (int)pow(2, _lod - (1 + i));
							int mip_y = y / (int)pow(2, _lod - (1 + i));
							auto& z_in_buffer = z_buffer[i][mip_y * (int)pow(2, i) + mip_x];
							for (int j = 0; j < 4; j++) {
								int sub_x = 2 * mip_x + j % 2;
								int sub_y = 2 * mip_y + j / 2;
								auto& z_in_sub = z_buffer[i + 1][sub_y * (int)pow(2, i + 1) + sub_x];
								z_in_buffer = box_min(z_in_sub, z_in_buffer);
							}
						}
					}
				}
				active_edge.xl += active_edge.dxl;
				active_edge.xl = box_max(polygon.xl, active_edge.xl);
				active_edge.xr += active_edge.dxr;
				active_edge.xr = box_min(polygon.xr, active_edge.xr);
				active_edge.zl += active_edge.dzx * active_edge.dxl + active_edge.dzy;

				// update active edge
				active_edge.dyl--;
				active_edge.dyr--;

				if (active_edge.dyl == 0 || active_edge.dyr == 0) {
					auto new_edge = edgeEntries[2];
					// dzx, dzy, zl and xl won't change
					if (active_edge.dyl == 0) {
						active_edge.dxl = new_edge.dx;
						active_edge.dyl = new_edge.dy - 1;
					}
					else {
						active_edge.dxr = new_edge.dx;
						active_edge.dyr = new_edge.dy - 1;
					}
				}
			}
		}
		stage.pop_back();
		if (!node->_isLeaf) {
			for (int i = 7; i >= 0; i--) {
				stage.push_back(node->_subNodes[i]);
			}
		}
	}
	TIMING_END("Rasterization Done")
}