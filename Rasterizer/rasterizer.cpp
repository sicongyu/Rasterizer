#include "rasterizer.h"
#include <memory>
#include <cmath>

#define FLT_MAX 3.402823466e+38F

Rasterizer::Rasterizer(int width, int height, const char* path) :
	width(width), height(height), model(path) {
	//framebuffer = new TGAColor[width * height];
	//memset(framebuffer, 0, sizeof(TGAColor) * width * height);
	framebuffer = new TGAImage(width, height, TGAImage::RGB);
	light_dir = glm::vec3(0.0f, 0.0f, 1.0f);
}

Rasterizer::~Rasterizer() {
	delete framebuffer;
}

ScanLine::ScanLine(int width, int height, const char* path) :
	Rasterizer(width, height, path) {
	//polygon_table = new std::vector<PolygonEntry>[height];
	//edge_table = new std::vector<EdgeEntry>[height];
	polygon_table = std::vector<std::unordered_map<int, PolygonEntry>>(height);
	edge_table = std::vector<std::unordered_multimap<int, EdgeEntry>>(height);
	z_buffer = new float[height * width];
	std::fill_n(z_buffer, height * width, -FLT_MAX);
	//memset(z_buffer, 128, sizeof(float) * width * height);
}

ScanLine::~ScanLine() {
	//delete[] polygon_table;
	//delete[] edge_table;
	delete[] z_buffer;
}

void ScanLine::draw() {
	TGAColor white = TGAColor(255, 255, 255, 255);
	// create polygon_table and edge_table
	for (int i = 0; i < model.num_faces; i++) {
		// for every vertices in the triangle
		glm::vec3 vertices[3];
		// map for vertices: ymax-0, ymid-1, ymin-2
		int y_map[3];
		float y_max = -FLT_MAX;
		float y_min = FLT_MAX;
		glm::vec2 abs_max = glm::vec2(box_max(abs(model.max.x), abs(model.min.x)), box_max(abs(model.max.y), abs(model.min.y)));
		glm::vec3 normal;
		for (int j = 0; j < 3; j++) {
			// convert into screen space
			//glm::vec3 vtx = (model.vertices[3 * i + j] - model.min) / (1.5f * (model.max - model.min));
			glm::vec3 vtx = model.vertices[3 * i + j];
			glm::vec2 vtx2d = (glm::vec2(vtx.x, vtx.y) / (1.5f * abs_max) + glm::vec2(1.0f)) * glm::vec2(width / 2, height / 2);
			vtx = glm::vec3(vtx2d.x, vtx2d.y, vtx.z);
			//vtx.x *= width;
			//vtx.y *= height;

			if (vtx.y > y_max) {
				y_max = vtx.y;
				y_map[0] = j;
			}
			if (vtx.y < y_min) {
				y_min = vtx.y;
				y_map[2] = j;
			}
			vertices[j] = vtx;
			normal = model.normals[3 * i]; // every three vtxs share the same normal
		}

		int polygon_dy = (int)y_max - (int)y_min + 1;
		if (polygon_dy == 1) {
			continue; // skip triangles that covers only one row;
		}

		// Occasionally when three vertices share the same y value, i.e., in the xoz plane
		if (y_map[0] == y_map[2]) {
			// assign an oreder manually
			y_map[2] = 2;
		}
		y_map[1] = 3 - y_map[0] - y_map[2];
		if (!(y_map[1] >= 0 && y_map[1] <= 2)) {
			throw "Error: invalid y_map index";
		}

		// put entry into polygon_table
		float intensity = glm::dot(normal, light_dir);
		PolygonEntry polygon_entry = {
			normal,
			i,
			polygon_dy,
			//white * intensity,
			//white, // debug purpose
			TGAColor(i * (255.0 / model.num_faces), 0, 0, 255)
			//TGAColor(255 * normal.x, 255 * normal.y, 255 * normal.z, 255)
		};
		//polygon_table[int(y_max)].push_back(polygon_entry);
		polygon_table[int(y_max)].emplace(polygon_entry.faceID, polygon_entry);

		// combination of (0, 1), (0, 2), (1, 2), put entries into edge_table
		for (int m = 0; m < 2; m++) {
			for (int n = m + 1; n < 3; n++) {
				glm::vec3 vtx0 = vertices[y_map[m]];
				glm::vec3 vtx1 = vertices[y_map[n]];

				int edge_dy = (int)vtx0.y - (int)vtx1.y + 1;
				if (edge_dy == 1) {
					continue; // discard edges that covers only one row
				}

				glm::vec3 edge = vtx1 - vtx0;
				EdgeEntry edge_entry = { 
					vtx0.x, //x of upper vertex
					vtx0.z, // z of upper vertex
					-(edge.x / edge.y), // dx (edge.y is negative)
					i, // faceID
					edge_dy // dy, notice that the height of vtx0 must be greater than vtx1
				};
				//edge_table[int(vtx0.y)].push_back(edge_entry);
				edge_table[int(vtx0.y)].emplace(edge_entry.faceID, edge_entry);
			}
		}
	}

	// scan process: from top to bottom
	for (int y = height-1; y >= 0; y--) {
		if (y == 460) {
			int err = 0;
		}
		// add new active polygons and edges
		if (!polygon_table[y].empty()) {
			for (auto polygon_table_iter : polygon_table[y]) {
				int faceID = polygon_table_iter.first;
				auto polygon = polygon_table_iter.second;
				// add entries to active polygon table
				//active_polygon_table.push_back(polygon);
				active_polygon_table.emplace(faceID, polygon);
				
				if (faceID == 13) {
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
				//for (edge : edge_table[y]) {
				//	if (edge.faceID == polygon.faceID) {
				//		if (!edge0) {
				//			edge0 = &edge;
				//		}
				//		else if (!edge1) {
				//			edge1 = &edge;
				//		}
				//	}
				//}

				// ensure that edge0 is the leftside edge
				if (edges[0].x_at_ymax > edges[1].x_at_ymax) {
					std::reverse(edges.begin(), edges.end());
				}
				//if (edge0->dx > edge1->dx) {
				//	std::swap(edge0, edge1);
				//}

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

			// update pixels from left to right
			for (int x = active_edge.xl; x < active_edge.xr; x++) {
				// update buffer contents
				float z = active_edge.zl + x * active_edge.dzx;
				if (z > z_buffer[y * width + x]) {
					z_buffer[y * width + x] = z;
					framebuffer->set(x, y, active_polygon_table[active_edge.faceID].color);
				}
			}

			// update active edge
			active_edge.dyl--;
			active_edge.dyr--;
			auto polygon_dy = --active_polygon_table[active_edge.faceID].dy;

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

			{
				active_edge.xl += active_edge.dxl;
				active_edge.xr += active_edge.dxr;
				active_edge.zl += active_edge.dzx * active_edge.dxl + active_edge.dzy;
			}
			++active_edge_pair_iter;
		}

		//// check if any polygon ends here (newly added edges may also end immediately -- may not handled properly)
		//for (auto active_polygon_iter = active_polygon_table.begin(); active_polygon_iter != active_polygon_table.end(); active_polygon_iter++) {
		//	active_polygon_iter->second.dy--;
		//	if (active_polygon_iter->second.dy == 0) {
		//		// remove the polygon itself
		//		active_polygon_table.erase(active_polygon_iter);
		//	}
		//}
	}
}