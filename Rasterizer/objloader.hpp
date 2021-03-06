#ifndef OBJLOADER_H
#define OBJLOADER_H

#include <vector>
#include <glm/glm.hpp>

#define box_max(a,b) (a>b?a:b)
#define box_min(a,b) (a>b?b:a)

class Model {

public:

	std::vector<glm::vec3> vertices;
	std::vector<glm::vec2> uvs;
	std::vector<glm::vec3> normals;
	glm::vec3 max;
	glm::vec3 min;
	unsigned int num_faces;

	Model();
	Model(const char* path);

private:
	bool loadOBJ(
		const char * path,
		std::vector<glm::vec3> & out_vertices,
		std::vector<glm::vec2> & out_uvs,
		std::vector<glm::vec3> & out_normals,
		glm::vec3& max,
		glm::vec3& min
	);

};

#endif