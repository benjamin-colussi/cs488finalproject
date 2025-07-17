///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////   CS 488/688 base code             /////////////////////////////////////////////////////////////////////////
//////////   (written by Toshiya Hachisuka)   /////////////////////////////////////////////////////////////////////////
//////////   (extended by Benjamin Colussi)   /////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// hello world
#pragma once
#define NOMINMAX

// linear algebra 
#include "linalg.h"
using namespace linalg::aliases;

// standard libraries
#include <iostream>
#include <cfloat>
#include <vector>

// global constants
constexpr float PI = 3.14159265358979f;
constexpr float PI_OVER_TWO = PI / 2.0f;
constexpr float PI_OVER_FOUR = PI / 4.0f;
constexpr float ONE_OVER_PI = 1.0f / PI;
constexpr float ONE_OVER_TWO_PI = 1.0f / (2.0f * PI);
constexpr float ONE_OVER_FOUR_PI = 1.0f / (4.0f * PI);
constexpr float Deg_To_Rad = PI / 180.0f;
constexpr float Rad_To_Deg = 180.0f / PI;
constexpr float EPSILON = 5e-5f;
constexpr float REFR_INDEX_AIR = 1.00029f;

// window size and resolution
constexpr int globalWidth = 512;
constexpr int globalHeight = 384;

// fixed camera parameters
constexpr float globalAspectRatio = static_cast<float>(globalWidth) / static_cast<float>(globalHeight);
constexpr float globalFOV = 45.0f;
constexpr float globalFilmSize = 0.032f;
const float globalDistanceToFilm = globalFilmSize / (2.0f * tan(globalFOV * Deg_To_Rad * 0.5f));

// dynamic camera parameters
constexpr float3 globalEye(0.0f, 0.0f, 1.5f);
constexpr float3 globalLookat(0.0f, 0.0f, 0.0f);
constexpr float3 globalUp(0.0f, 1.0f, 0.0f);
const float3 globalViewDir = normalize(globalLookat - globalEye);
const float3 globalRight = normalize(cross(globalViewDir, globalUp));

// fast random number generator based pcg32_fast
#include <stdint.h>
namespace PCG32 {
	static uint64_t mcg_state = 0xcafef00dd15ea5e5u; // must be odd
	static uint64_t const multiplier = 6364136223846793005u;
	uint32_t pcg32_fast(void) {
		uint64_t x = mcg_state;
		const unsigned count = (unsigned)(x >> 61);
		mcg_state = x * multiplier;
		x ^= x >> 22;
		return (uint32_t)(x >> (22 + count));
	}
	float rand() {
		return float(double(pcg32_fast()) / 4294967296.0);
	}
}

// random normal generator
#include <random>
constexpr float mu = 0.0f;
constexpr float sigma = 1.0f;
std::random_device rd; // random seed
std::mt19937 gen(rd()); // mersenne twister engine
std::normal_distribution<float> std_norm(mu, sigma);

// global constants and variables
constexpr int MAX_REFLECTION_RECURSION_DEPTH = 4;
constexpr int MAX_REFRACTION_RECURSION_DEPTH = 4;
constexpr int MAXIMUM_PATH_LENGTH = 0;

// switches
constexpr bool SAMPLE_UNIFORM_HEMISPHERE = false;
constexpr bool SAMPLE_COS_WEIGHTED_HEMISPHERE = true;
constexpr bool SAMPLE_UNIFORM_HEMISPHERE_LIKE_A_NOOB = false;



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// image
class Image {

	std::vector<float3> pixels;
	int width = 0, height = 0;

	public:

		void resize(const int newWidth, const int newHeight) {
			this->pixels.resize(newWidth * newHeight);
			this->width = newWidth;
			this->height = newHeight;
		}

		void clear() {
			for (int j = 0; j < height; j++) {
				for (int i = 0; i < width; i++) {
					this->pixel(i, j) = float3(0.0f);
				}
			}
		}

		explicit Image(const int newWidth = 0, const int newHeight = 0) {
			this->resize(newWidth, newHeight);
			this->clear();
		}

		float3& pixel(const int i, const int j) {
			return this->pixels[i + j * width];
		}

};

// final image to be computed
Image globalImage(globalWidth, globalHeight);



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// material type
enum MaterialType {
	LAMBERTIAN,
	METAL,
	GLASS,
	LIGHT
};

// uber material
class Material {

	public:

		// fields
		std::string name;
		MaterialType type = LAMBERTIAN;
		float eta = 1.0f; // index of refraction
		float glossiness = 1.0f;
		float3 Ka = float3(0.0f); // ambient colour
		float3 Kd = float3(0.9f); // diffuse colour
		float3 Ks = float3(0.0f); // specular colour
		float Ns = 0.0; // specular exponent
		float3 emission;

		// support 8-bit texture
		bool isTextured = false;
		unsigned char* texture = nullptr;
		int textureWidth = 0;
		int textureHeight = 0;

		// functions
		Material() = default;
		~Material() = default;
		void setReflectance(const float3& c) {
			if (type == LAMBERTIAN) {
				Kd = c;
			}
			else if (type == METAL) {
				// empty
			}
			else if (type == GLASS) {
				// empty
			}
		}
		float3 fetchTexture(const float2& tex) const {
			int x = int(tex.x * textureWidth) % textureWidth;
			int y = int(tex.y * textureHeight) % textureHeight;
			if (x < 0) x += textureWidth;
			if (y < 0) y += textureHeight;
			int pix = (x + y * textureWidth) * 3;
			const unsigned char r = texture[pix + 0];
			const unsigned char g = texture[pix + 1];
			const unsigned char b = texture[pix + 2];
			return float3(r, g, b) / 255.0f;
		}

		// material spectrum
		float3 spectrum() const {

			// perfect diffuse
			if (type == LAMBERTIAN) return Kd * ONE_OVER_PI;

			// perfect specular reflection
			if (type == METAL) {}

			// default
			return float3(0.0f);
		}

		// sample direction
		float3 sampleDirection(const float3& wo, const float3& n) const {

			// perfect diffuse
			if (type == LAMBERTIAN) {

				// using normalized 3d standard normal
				if (SAMPLE_UNIFORM_HEMISPHERE) {
					float3 d = float3(std_norm(gen), std_norm(gen), std_norm(gen));
					while (d.x == 0 && d.y == 0 && d.z == 0) d = float3(std_norm(gen), std_norm(gen), std_norm(gen));
					d = normalize(d);
					if (dot(n, d) < 0) d *= -1;
					return d;
				}

				// using malley's method
				else if (SAMPLE_COS_WEIGHTED_HEMISPHERE) {

					// sample cos weighted positive unit hemisphere
					float x = 2 * PCG32::rand() - 1;
					float y = 2 * PCG32::rand() - 1;
					if (x == 0 && y == 0) return float3(x, y, 1);
					float theta, r;
					if (abs(x) > abs(y)) {
						r = x;
						theta = PI_OVER_FOUR * (y / x);
					}
					else {
						r = y;
						theta = PI_OVER_TWO - PI_OVER_FOUR * (x / y);
					}
					float3 d = float3(r * cos(theta), r * sin(theta), sqrtf(1 - x * x - y * y));

					// build orthonormal basis with normal
					float sign = copysignf(1, n.z);
					const float a = -1 / (sign + n.z);
					const float b = n.x * n.y * a;
					const float3 b1 = float3(1 + sign * n.x * n.x * a, sign * b, -sign * n.x);
					const float3 b2 = float3(b, sign + n.y * n.y * a, -n.y);

					// centre about normal
					return d.x * b1 + d.y * b2 + d.z * n;
				}

				// checking if this works ...
				else if (SAMPLE_UNIFORM_HEMISPHERE_LIKE_A_NOOB) {
					
					// using spherical coordinates
					const float z = PCG32::rand();
					const float r = sqrtf(1 - z * z);
					const float phi = 2 * PI * PCG32::rand();
					const float3 d = float3(r * std::cos(phi), r * std::sin(phi), z);

					// build orthonormal basis with normal
					float sign = copysignf(1, n.z);
					const float a = -1 / (sign + n.z);
					const float b = n.x * n.y * a;
					const float3 b1 = float3(1 + sign * n.x * n.x * a, sign * b, -sign * n.x);
					const float3 b2 = float3(b, sign + n.y * n.y * a, -n.y);

					// centre about normal
					return d.x * b1 + d.y * b2 + d.z * n;
				}
			}

			// perfect specular reflection
			if (type == METAL) {
				float3 reflection = wo - 2 * dot(wo, n) * n;
				// if (dot(reflection, n) < 0) reflection -= 2 * dot(reflection, n) * n; // unecessary? ...
				return normalize(reflection);
			}

			// default
			return float3(0.0f);
		}

		// sample pdf
		float pdf(const float3& n, const float3& wi) const {

			// perfect diffuse
			if (type == LAMBERTIAN) {

				// hemisphere subtends 2 Pi steradians
				if (SAMPLE_UNIFORM_HEMISPHERE) return ONE_OVER_TWO_PI;

				// inputs must be normalized
				else if (SAMPLE_COS_WEIGHTED_HEMISPHERE) return dot(n, wi) * ONE_OVER_PI;

				// checking if this works ...
				else if (SAMPLE_UNIFORM_HEMISPHERE_LIKE_A_NOOB) return ONE_OVER_TWO_PI;
			}

			// perfect specular reflection
			if (type == METAL) {}

			// default
			return 0.0f;
		}

};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// ray
class Ray {

	public:

		float3 o, d;
		Ray() : o(), d(float3(0.0f, 0.0f, 1.0f)) {}
		Ray(const float3& o, const float3& d) : o(o), d(d) {}

};

// hit info
class HitInfo {

	public:

		float t; // distance
		float3 P; // location
		float3 N; // shading normal vector
		float2 T; // texture coordinate
		const Material* material; // material of object hit

		float3 G; // geometric normal
		bool front; // hit front or back

};

// axis-aligned bounding box
class AABB {

	private:

		float3 minp, maxp, size;

	public:

		float3 get_minp() { return minp; }
		float3 get_maxp() { return maxp; }
		float3 get_size() { return size; }

		AABB() {
			minp = float3(FLT_MAX);
			maxp = float3(-FLT_MAX);
			size = float3(0.0f);
		}

		void reset() {
			minp = float3(FLT_MAX);
			maxp = float3(-FLT_MAX);
			size = float3(0.0f);
		}

		int getLargestAxis() const {
			if ((size.x > size.y) && (size.x > size.z)) return 0;
			if (size.y > size.z) return 1;
			return 2;
		}

		void fit(const float3& v) {
			if (minp.x > v.x) minp.x = v.x;
			if (minp.y > v.y) minp.y = v.y;
			if (minp.z > v.z) minp.z = v.z;
			if (maxp.x < v.x) maxp.x = v.x;
			if (maxp.y < v.y) maxp.y = v.y;
			if (maxp.z < v.z) maxp.z = v.z;
			size = maxp - minp;
		}

		float area() const { return 2.0f * (size.x * size.y + size.y * size.z + size.z * size.x); }

		bool intersect(HitInfo& minHit, const Ray& ray) const {
			// set minHit.t as the distance to the intersection point
			// return true/false if the ray hits or not
			float tx1 = (minp.x - ray.o.x) / ray.d.x;
			float ty1 = (minp.y - ray.o.y) / ray.d.y;
			float tz1 = (minp.z - ray.o.z) / ray.d.z;

			float tx2 = (maxp.x - ray.o.x) / ray.d.x;
			float ty2 = (maxp.y - ray.o.y) / ray.d.y;
			float tz2 = (maxp.z - ray.o.z) / ray.d.z;

			if (tx1 > tx2) {
				const float temp = tx1;
				tx1 = tx2;
				tx2 = temp;
			}

			if (ty1 > ty2) {
				const float temp = ty1;
				ty1 = ty2;
				ty2 = temp;
			}

			if (tz1 > tz2) {
				const float temp = tz1;
				tz1 = tz2;
				tz2 = temp;
			}

			float t1 = tx1; if (t1 < ty1) t1 = ty1; if (t1 < tz1) t1 = tz1;
			float t2 = tx2; if (t2 > ty2) t2 = ty2; if (t2 > tz2) t2 = tz2;

			if (t1 > t2) return false;
			if ((t1 < 0.0) && (t2 < 0.0)) return false;

			minHit.t = t1;
			return true;
		}

};

// triangle
struct Triangle {
	float3 positions[3];
	float3 normals[3];
	float2 texcoords[3];
	int idMaterial = 0;
	AABB bbox;
	float3 center;
};

// triangle mesh
class TriangleMesh {

	public:

		// fields
		std::vector<Triangle> triangles;
		std::vector<Material> materials;
		AABB bbox;

		// ray trace against triangle mesh
		bool raytraceTriangle(HitInfo& result, const Ray& ray, const Triangle& tri, float tMin, float tMax) const {

			// ray-triangle intersection, Cramer's rule, A = [a-b  a-c  d], x = [beta  gamma  t], b = [a-o]
			const float3 ba = tri.positions[0] - tri.positions[1];
			const float3 ca = tri.positions[0] - tri.positions[2];
			const float detA = dot(cross(ba, ca), ray.d);
			if (detA == 0) return false;

			const float3 oa = tri.positions[0] - ray.o;
			const float detAbeta = dot(cross(oa, ca), ray.d);
			const float detAgamma = dot(cross(ba, oa), ray.d);
			const float detAt = dot(cross(ba, ca), oa);

			if (detA > 0) {
				if (detAbeta < 0 || detA < detAbeta) return false;
				if (detAgamma < 0 || detA < detAgamma) return false;
				if (detA - detAbeta - detAgamma < 0 || detA < detA - detAbeta - detAgamma) return false;
				if (detAt < tMin * detA || tMax * detA < detAt) return false;
			}
			else if (detA < 0) {
				if (detAbeta > 0 || detA > detAbeta) return false;
				if (detAgamma > 0 || detA > detAgamma) return false;
				if (detA - detAbeta - detAgamma > 0 || detA > detA - detAbeta - detAgamma) return false;
				if (detAt > tMin * detA || tMax * detA > detAt) return false;
			}

			// only after checking
			const float beta = detAbeta / detA;
			const float gamma = detAgamma / detA;
			const float t = detAt / detA;
			const float alpha = 1.0f - beta - gamma; // might want to check these sum to less than 1

			// fill in "result" when there is an intersection
			result.t = t;
			result.P = alpha * tri.positions[0] + beta * tri.positions[1] + gamma * tri.positions[2];
			result.N = normalize(alpha * tri.normals[0] + beta * tri.normals[1] + gamma * tri.normals[2]);
			result.T = alpha * tri.texcoords[0] + beta * tri.texcoords[1] + gamma * tri.texcoords[2];
			result.material = &materials[tri.idMaterial];

			// calculate geometric normal to front of triangle
			float3 g = normalize(cross(ba, ca));
			// assuming the interpolated shading normal is pointing towards the front of the triangle
			if (dot(g, result.N) < 0) g *= -1;
			result.G = g;

			// check if hitting front of triangle
			if (dot(-ray.d, result.G) > 0) result.front = true;
			else result.front = false;

			// return true/false if there is an intersection or not
			return true;
		}

		// some precalculation for bounding boxes
		void preCalc() {
			bbox.reset();
			for (int i = 0, _n = (int)triangles.size(); i < _n; i++) {
				this->triangles[i].bbox.reset();
				this->triangles[i].bbox.fit(this->triangles[i].positions[0]);
				this->triangles[i].bbox.fit(this->triangles[i].positions[1]);
				this->triangles[i].bbox.fit(this->triangles[i].positions[2]);
				this->triangles[i].center = (this->triangles[i].positions[0] + this->triangles[i].positions[1] + this->triangles[i].positions[2]) * (1.0f / 3.0f);
				this->bbox.fit(this->triangles[i].positions[0]);
				this->bbox.fit(this->triangles[i].positions[1]);
				this->bbox.fit(this->triangles[i].positions[2]);
			}
		}

		// load .obj file
		bool load(const char* filename, const float4x4& ctm = linalg::identity) {
			int nVertices = 0;
			float* vertices;
			float* normals;
			float* texcoords;
			int nIndices;
			int* indices;
			int* matid = nullptr;

			printf("Loading \"%s\"...\n", filename);
			ParseOBJ(filename, nVertices, &vertices, &normals, &texcoords, nIndices, &indices, &matid);
			if (nVertices == 0) return false;
			this->triangles.resize(nIndices / 3);

			if (matid != nullptr) {
				for (unsigned int i = 0; i < materials.size(); i++) {
					// convert .mlt data into BSDF definitions
					// you may change the followings in the final project if you want
					materials[i].type = LAMBERTIAN;
					if (materials[i].Ns == 100.0f) {
						materials[i].type = METAL;
					}
					if (materials[i].name.compare(0, 5, "glass", 0, 5) == 0) {
						materials[i].type = GLASS;
						materials[i].eta = 1.5f;
					}
				}
			} else {
				// use default Lambertian
				this->materials.resize(1);
			}

			for (unsigned int i = 0; i < this->triangles.size(); i++) {
				const int v0 = indices[i * 3 + 0];
				const int v1 = indices[i * 3 + 1];
				const int v2 = indices[i * 3 + 2];

				this->triangles[i].positions[0] = float3(vertices[v0 * 3 + 0], vertices[v0 * 3 + 1], vertices[v0 * 3 + 2]);
				this->triangles[i].positions[1] = float3(vertices[v1 * 3 + 0], vertices[v1 * 3 + 1], vertices[v1 * 3 + 2]);
				this->triangles[i].positions[2] = float3(vertices[v2 * 3 + 0], vertices[v2 * 3 + 1], vertices[v2 * 3 + 2]);

				if (normals != nullptr) {
					this->triangles[i].normals[0] = float3(normals[v0 * 3 + 0], normals[v0 * 3 + 1], normals[v0 * 3 + 2]);
					this->triangles[i].normals[1] = float3(normals[v1 * 3 + 0], normals[v1 * 3 + 1], normals[v1 * 3 + 2]);
					this->triangles[i].normals[2] = float3(normals[v2 * 3 + 0], normals[v2 * 3 + 1], normals[v2 * 3 + 2]);
				} else {
					// no normal data, calculate the normal for a polygon
					const float3 e0 = this->triangles[i].positions[1] - this->triangles[i].positions[0];
					const float3 e1 = this->triangles[i].positions[2] - this->triangles[i].positions[0];
					const float3 n = normalize(cross(e0, e1));

					this->triangles[i].normals[0] = n;
					this->triangles[i].normals[1] = n;
					this->triangles[i].normals[2] = n;
				}

				// material id
				this->triangles[i].idMaterial = 0;
				if (matid != nullptr) {
					// read texture coordinates
					if ((texcoords != nullptr) && materials[matid[i]].isTextured) {
						this->triangles[i].texcoords[0] = float2(texcoords[v0 * 2 + 0], texcoords[v0 * 2 + 1]);
						this->triangles[i].texcoords[1] = float2(texcoords[v1 * 2 + 0], texcoords[v1 * 2 + 1]);
						this->triangles[i].texcoords[2] = float2(texcoords[v2 * 2 + 0], texcoords[v2 * 2 + 1]);
					} else {
						this->triangles[i].texcoords[0] = float2(0.0f);
						this->triangles[i].texcoords[1] = float2(0.0f);
						this->triangles[i].texcoords[2] = float2(0.0f);
					}
					this->triangles[i].idMaterial = matid[i];
				} else {
					this->triangles[i].texcoords[0] = float2(0.0f);
					this->triangles[i].texcoords[1] = float2(0.0f);
					this->triangles[i].texcoords[2] = float2(0.0f);
				}
			}
			printf("Loaded \"%s\" with %d triangles.\n", filename, int(triangles.size()));

			delete[] vertices;
			delete[] normals;
			delete[] texcoords;
			delete[] indices;
			delete[] matid;

			return true;
		}

		// destructonator
		~TriangleMesh() {
			materials.clear();
			triangles.clear();
		}

		// single triangle by default
		void createSingleTriangle() {
			triangles.resize(1);
			materials.resize(1);
			triangles[0].idMaterial = 0;
			triangles[0].positions[0] = float3(-0.5f, -0.5f, 0.0f);
			triangles[0].positions[1] = float3(0.5f, -0.5f, 0.0f);
			triangles[0].positions[2] = float3(0.0f, 0.5f, 0.0f);
			const float3 e0 = this->triangles[0].positions[1] - this->triangles[0].positions[0];
			const float3 e1 = this->triangles[0].positions[2] - this->triangles[0].positions[0];
			const float3 n = normalize(cross(e0, e1));
			triangles[0].normals[0] = n;
			triangles[0].normals[1] = n;
			triangles[0].normals[2] = n;
			triangles[0].texcoords[0] = float2(0.0f, 0.0f);
			triangles[0].texcoords[1] = float2(0.0f, 1.0f);
			triangles[0].texcoords[2] = float2(1.0f, 0.0f);
		}

	private:

		std::string GetBaseDir(const std::string& filepath) {
			if (filepath.find_last_of("/\\") != std::string::npos) return filepath.substr(0, filepath.find_last_of("/\\"));
			return "";
		}
		std::string base_dir;

		void LoadMTL(const std::string fileName) {
			FILE* fp = fopen(fileName.c_str(), "r");

			Material mtl;
			mtl.texture = nullptr;
			char line[81];
			while (fgets(line, 80, fp) != nullptr) {
				float r, g, b, s;
				std::string lineStr;
				lineStr = line;
				int i = int(materials.size());

				if (lineStr.compare(0, 6, "newmtl", 0, 6) == 0) {
					lineStr.erase(0, 7);
					mtl.name = lineStr;
					mtl.isTextured = false;
				} else if (lineStr.compare(0, 2, "Ka", 0, 2) == 0) {
					lineStr.erase(0, 3);
					sscanf(lineStr.c_str(), "%f %f %f\n", &r, &g, &b);
					mtl.Ka = float3(r, g, b);
				} else if (lineStr.compare(0, 2, "Kd", 0, 2) == 0) {
					lineStr.erase(0, 3);
					sscanf(lineStr.c_str(), "%f %f %f\n", &r, &g, &b);
					mtl.Kd = float3(r, g, b);
				} else if (lineStr.compare(0, 2, "Ks", 0, 2) == 0) {
					lineStr.erase(0, 3);
					sscanf(lineStr.c_str(), "%f %f %f\n", &r, &g, &b);
					mtl.Ks = float3(r, g, b);
				} else if (lineStr.compare(0, 2, "Ns", 0, 2) == 0) {
					lineStr.erase(0, 3);
					sscanf(lineStr.c_str(), "%f\n", &s);
					mtl.Ns = s;
					mtl.texture = nullptr;
					materials.push_back(mtl);
				} else if (lineStr.compare(0, 6, "map_Kd", 0, 6) == 0) {
					lineStr.erase(0, 7);
					lineStr.erase(lineStr.size() - 1, 1);
					materials[i - 1].isTextured = true;
				}
			}

			fclose(fp);
		}

		void ParseOBJ(const char* fileName, int& nVertices, float** vertices, float** normals, float** texcoords, int& nIndices, int** indices, int** materialids) {
		// local function in C++...
		struct {
			void operator()(char* word, int* vindex, int* tindex, int* nindex) {
				const char* null = " ";
				char* ptr;
				const char* tp;
				const char* np;

				// by default, the texture and normal pointers are set to the null string
				tp = null;
				np = null;

				// replace slashes with null characters and cause tp and np to point
				// to character immediately following the first or second slash
				for (ptr = word; *ptr != '\0'; ptr++) {
					if (*ptr == '/') {
						if (tp == null) {
							tp = ptr + 1;
						} else {
							np = ptr + 1;
						}

						*ptr = '\0';
					}
				}

				*vindex = atoi(word);
				*tindex = atoi(tp);
				*nindex = atoi(np);
			}
		} get_indices;

		base_dir = GetBaseDir(fileName);
		#ifdef _WIN32
			base_dir += "\\";
		#else
			base_dir += "/";
		#endif

		FILE* fp = fopen(fileName, "r");
		int nv = 0, nn = 0, nf = 0, nt = 0;
		char line[81];
		if (!fp) {
			printf("Cannot open \"%s\" for reading\n", fileName);
			return;
		}

		while (fgets(line, 80, fp) != NULL) {
			std::string lineStr;
			lineStr = line;

			if (lineStr.compare(0, 6, "mtllib", 0, 6) == 0) {
				lineStr.erase(0, 7);
				lineStr.erase(lineStr.size() - 1, 1);
				LoadMTL(base_dir + lineStr);
			}

			if (line[0] == 'v') {
				if (line[1] == 'n') {
					nn++;
				} else if (line[1] == 't') {
					nt++;
				} else {
					nv++;
				}
			} else if (line[0] == 'f') {
				nf++;
			}
		}
		fseek(fp, 0, 0);

		float* n = new float[3 * (nn > nf ? nn : nf)];
		float* v = new float[3 * nv];
		float* t = new float[2 * nt];

		int* vInd = new int[3 * nf];
		int* nInd = new int[3 * nf];
		int* tInd = new int[3 * nf];
		int* mInd = new int[nf];

		int nvertices = 0;
		int nnormals = 0;
		int ntexcoords = 0;
		int nindices = 0;
		int ntriangles = 0;
		bool noNormals = false;
		bool noTexCoords = false;
		bool noMaterials = true;
		int cmaterial = 0;

		while (fgets(line, 80, fp) != NULL) {
			std::string lineStr;
			lineStr = line;

			if (line[0] == 'v') {
				if (line[1] == 'n') {
					float x, y, z;
					sscanf(&line[2], "%f %f %f\n", &x, &y, &z);
					float l = sqrt(x * x + y * y + z * z);
					x = x / l;
					y = y / l;
					z = z / l;
					n[nnormals] = x;
					nnormals++;
					n[nnormals] = y;
					nnormals++;
					n[nnormals] = z;
					nnormals++;
				} else if (line[1] == 't') {
					float u, v;
					sscanf(&line[2], "%f %f\n", &u, &v);
					t[ntexcoords] = u;
					ntexcoords++;
					t[ntexcoords] = v;
					ntexcoords++;
				} else {
					float x, y, z;
					sscanf(&line[1], "%f %f %f\n", &x, &y, &z);
					v[nvertices] = x;
					nvertices++;
					v[nvertices] = y;
					nvertices++;
					v[nvertices] = z;
					nvertices++;
				}
			}
			if (lineStr.compare(0, 6, "usemtl", 0, 6) == 0) {
				lineStr.erase(0, 7);
				if (materials.size() != 0) {
					for (unsigned int i = 0; i < materials.size(); i++) {
						if (lineStr.compare(materials[i].name) == 0) {
							cmaterial = i;
							noMaterials = false;
							break;
						}
					}
				}

			} else if (line[0] == 'f') {
				char s1[32], s2[32], s3[32];
				int vI, tI, nI;
				sscanf(&line[1], "%s %s %s\n", s1, s2, s3);

				mInd[ntriangles] = cmaterial;

				// indices for first vertex
				get_indices(s1, &vI, &tI, &nI);
				vInd[nindices] = vI - 1;
				if (nI) {
					nInd[nindices] = nI - 1;
				} else {
					noNormals = true;
				}

				if (tI) {
					tInd[nindices] = tI - 1;
				} else {
					noTexCoords = true;
				}
				nindices++;

				// indices for second vertex
				get_indices(s2, &vI, &tI, &nI);
				vInd[nindices] = vI - 1;
				if (nI) {
					nInd[nindices] = nI - 1;
				} else {
					noNormals = true;
				}

				if (tI) {
					tInd[nindices] = tI - 1;
				} else {
					noTexCoords = true;
				}
				nindices++;

				// indices for third vertex
				get_indices(s3, &vI, &tI, &nI);
				vInd[nindices] = vI - 1;
				if (nI) {
					nInd[nindices] = nI - 1;
				} else {
					noNormals = true;
				}

				if (tI) {
					tInd[nindices] = tI - 1;
				} else {
					noTexCoords = true;
				}
				nindices++;

				ntriangles++;
			}
		}

		*vertices = new float[ntriangles * 9];
		if (!noNormals) {
			*normals = new float[ntriangles * 9];
		} else {
			*normals = 0;
		}

		if (!noTexCoords) {
			*texcoords = new float[ntriangles * 6];
		} else {
			*texcoords = 0;
		}

		if (!noMaterials) {
			*materialids = new int[ntriangles];
		} else {
			*materialids = 0;
		}

		*indices = new int[ntriangles * 3];
		nVertices = ntriangles * 3;
		nIndices = ntriangles * 3;

		for (int i = 0; i < ntriangles; i++) {
			if (!noMaterials) {
				(*materialids)[i] = mInd[i];
			}

			(*indices)[3 * i] = 3 * i;
			(*indices)[3 * i + 1] = 3 * i + 1;
			(*indices)[3 * i + 2] = 3 * i + 2;

			(*vertices)[9 * i] = v[3 * vInd[3 * i]];
			(*vertices)[9 * i + 1] = v[3 * vInd[3 * i] + 1];
			(*vertices)[9 * i + 2] = v[3 * vInd[3 * i] + 2];

			(*vertices)[9 * i + 3] = v[3 * vInd[3 * i + 1]];
			(*vertices)[9 * i + 4] = v[3 * vInd[3 * i + 1] + 1];
			(*vertices)[9 * i + 5] = v[3 * vInd[3 * i + 1] + 2];

			(*vertices)[9 * i + 6] = v[3 * vInd[3 * i + 2]];
			(*vertices)[9 * i + 7] = v[3 * vInd[3 * i + 2] + 1];
			(*vertices)[9 * i + 8] = v[3 * vInd[3 * i + 2] + 2];

			if (!noNormals) {
				(*normals)[9 * i] = n[3 * nInd[3 * i]];
				(*normals)[9 * i + 1] = n[3 * nInd[3 * i] + 1];
				(*normals)[9 * i + 2] = n[3 * nInd[3 * i] + 2];

				(*normals)[9 * i + 3] = n[3 * nInd[3 * i + 1]];
				(*normals)[9 * i + 4] = n[3 * nInd[3 * i + 1] + 1];
				(*normals)[9 * i + 5] = n[3 * nInd[3 * i + 1] + 2];

				(*normals)[9 * i + 6] = n[3 * nInd[3 * i + 2]];
				(*normals)[9 * i + 7] = n[3 * nInd[3 * i + 2] + 1];
				(*normals)[9 * i + 8] = n[3 * nInd[3 * i + 2] + 2];
			}

			if (!noTexCoords) {
				(*texcoords)[6 * i] = t[2 * tInd[3 * i]];
				(*texcoords)[6 * i + 1] = t[2 * tInd[3 * i] + 1];

				(*texcoords)[6 * i + 2] = t[2 * tInd[3 * i + 1]];
				(*texcoords)[6 * i + 3] = t[2 * tInd[3 * i + 1] + 1];

				(*texcoords)[6 * i + 4] = t[2 * tInd[3 * i + 2]];
				(*texcoords)[6 * i + 5] = t[2 * tInd[3 * i + 2] + 1];
			}

		}
		fclose(fp);

		delete[] n;
		delete[] v;
		delete[] t;
		delete[] nInd;
		delete[] vInd;
		delete[] tInd;
		delete[] mInd;
	}

};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// BVH node
class BVHNode {
public:
	bool isLeaf; // is the node a leaf
	int idLeft, idRight; // left and right children
	int triListNum; // number of triangles
	int* triList; // triangle list
	AABB bbox; // axis-aligned bounding box
};

// bounding volume hierarchy
class BVH {

	public:

		// fields
		const TriangleMesh* triangleMesh = nullptr;
		BVHNode* node = nullptr;
		const float costBBox = 1.0f;
		const float costTri = 1.0f;
		int leafNum = 0;
		int nodeNum = 0;

		// default ctor
		BVH() {}

		// builds the BVH based on the triangle mesh object
		void build(const TriangleMesh* mesh);

		// checks for closest hit by traversing the BVH nodes
		bool intersect(HitInfo& result, const Ray& ray, float tMin = 0.0f, float tMax = FLT_MAX) const {
			bool hit = false;
			HitInfo tempMinHit;
			result.t = FLT_MAX;
			if (this->node[0].bbox.intersect(tempMinHit, ray)) hit = traverse(result, ray, 0, tMin, tMax);
			if (result.t != FLT_MAX) hit = true;
			return hit;
		}

		// traverses the BVH nodes
		bool traverse(HitInfo& result, const Ray& ray, int node_id, float tMin, float tMax) const;

	private:

		void sortAxis(int* obj_index, const char axis, const int li, const int ri) const;
		int splitBVH(int* obj_index, const int obj_num, const AABB& bbox);

};

// sort bounding boxes (in case you want to build SAH-BVH)
void BVH::sortAxis(int* obj_index, const char axis, const int li, const int ri) const {
	int i, j;
	float pivot;
	int temp;
	i = li;
	j = ri;
	pivot = triangleMesh->triangles[obj_index[(li + ri) / 2]].center[axis];
	while (true) {
		while (triangleMesh->triangles[obj_index[i]].center[axis] < pivot) ++i;
		while (triangleMesh->triangles[obj_index[j]].center[axis] > pivot) --j;
		if (i >= j) break;
		temp = obj_index[i];
		obj_index[i] = obj_index[j];
		obj_index[j] = temp;
		++i;
		--j;
	}
	if (li < (i - 1)) sortAxis(obj_index, axis, li, i - 1);
	if ((j + 1) < ri) sortAxis(obj_index, axis, j + 1, ri);
}

// #define SAHBVH // use this in once you have SAH-BVH
int BVH::splitBVH(int* obj_index, const int obj_num, const AABB& bbox) {

#ifndef SAHBVH

	int bestAxis, bestIndex;
	AABB bboxL, bboxR, bestbboxL, bestbboxR;
	int* sorted_obj_index = new int[obj_num];

	// split along the largest axis
	bestAxis = bbox.getLargestAxis();

	// sorting along the axis
	this->sortAxis(obj_index, bestAxis, 0, obj_num - 1);
	for (int i = 0; i < obj_num; ++i) sorted_obj_index[i] = obj_index[i];

	// split in the middle
	bestIndex = obj_num / 2 - 1;

	bboxL.reset();
	for (int i = 0; i <= bestIndex; ++i) {
		const Triangle& tri = triangleMesh->triangles[obj_index[i]];
		bboxL.fit(tri.positions[0]);
		bboxL.fit(tri.positions[1]);
		bboxL.fit(tri.positions[2]);
	}

	bboxR.reset();
	for (int i = bestIndex + 1; i < obj_num; ++i) {
		const Triangle& tri = triangleMesh->triangles[obj_index[i]];
		bboxR.fit(tri.positions[0]);
		bboxR.fit(tri.positions[1]);
		bboxR.fit(tri.positions[2]);
	}

	bestbboxL = bboxL;
	bestbboxR = bboxR;

#else // implement SAH-BVH here

	// x-axis
	AABB* bboxL = new AABB[obj_num];
	AABB* bboxR = new AABB[obj_num];
	int* sorted_obj_index_x = new int[obj_num];
	this->sortAxis(obj_index, 0, 0, obj_num - 1);
	for (int i = 0; i < obj_num; ++i) sorted_obj_index_x[i] = obj_index[i];

	// sweep min to max
	bboxL[0].reset();
	bboxL[0].fit(triangleMesh->triangles[obj_index[0]].positions[0]);
	bboxL[0].fit(triangleMesh->triangles[obj_index[0]].positions[1]);
	bboxL[0].fit(triangleMesh->triangles[obj_index[0]].positions[2]);
	for (int i = 1; i < obj_num; ++i) {
		bboxL[i].reset();
		bboxL[i] = bboxL[i - 1];
		const Triangle& tri = triangleMesh->triangles[obj_index[i]];
		bboxL[i].fit(tri.positions[0]);
		bboxL[i].fit(tri.positions[1]);
		bboxL[i].fit(tri.positions[2]);
	}

	// sweep max to min
	bboxR[0].reset();
	bboxR[0].fit(triangleMesh->triangles[obj_index[obj_num - 1]].positions[0]);
	bboxR[0].fit(triangleMesh->triangles[obj_index[obj_num - 1]].positions[1]);
	bboxR[0].fit(triangleMesh->triangles[obj_index[obj_num - 1]].positions[2]);
	for (int i = 1; i < obj_num; ++i) {
		bboxR[i].reset();
		bboxR[i] = bboxR[i - 1];
		const Triangle& tri = triangleMesh->triangles[obj_index[obj_num - i - 1]];
		bboxR[i].fit(tri.positions[0]);
		bboxR[i].fit(tri.positions[1]);
		bboxR[i].fit(tri.positions[2]);
	}

	// calculate SAH
	float SAHcostminX = FLT_MAX;
	int bestIndexX = INT_MAX;
	for (int i = 0; i < obj_num; ++i) {
		float SAHcost = 2 * bbox.area() + (bboxL[i].area() * (i + 1) + bboxR[obj_num - i - 1].area() * (obj_num - i));
		if (SAHcost < SAHcostminX) {
			SAHcostminX = SAHcost;
			bestIndexX = i;
		}
	}

	// y-axis
	AABB* bboyL = new AABB[obj_num];
	AABB* bboyR = new AABB[obj_num];
	int* sorted_obj_index_y = new int[obj_num];
	this->sortAxis(obj_index, 1, 0, obj_num - 1);
	for (int i = 0; i < obj_num; ++i) sorted_obj_index_y[i] = obj_index[i];

	// sweep min to max
	bboyL[0].reset();
	bboyL[0].fit(triangleMesh->triangles[obj_index[0]].positions[0]);
	bboyL[0].fit(triangleMesh->triangles[obj_index[0]].positions[1]);
	bboyL[0].fit(triangleMesh->triangles[obj_index[0]].positions[2]);
	for (int i = 1; i < obj_num; ++i) {
		bboyL[i].reset();
		bboyL[i] = bboyL[i - 1];
		const Triangle& tri = triangleMesh->triangles[obj_index[i]];
		bboyL[i].fit(tri.positions[0]);
		bboyL[i].fit(tri.positions[1]);
		bboyL[i].fit(tri.positions[2]);
	}

	// sweep max to min
	bboyR[0].reset();
	bboyR[0].fit(triangleMesh->triangles[obj_index[obj_num - 1]].positions[0]);
	bboyR[0].fit(triangleMesh->triangles[obj_index[obj_num - 1]].positions[1]);
	bboyR[0].fit(triangleMesh->triangles[obj_index[obj_num - 1]].positions[2]);
	for (int i = 1; i < obj_num; ++i) {
		bboyR[i].reset();
		bboyR[i] = bboyR[i - 1];
		const Triangle& tri = triangleMesh->triangles[obj_index[obj_num - i - 1]];
		bboyR[i].fit(tri.positions[0]);
		bboyR[i].fit(tri.positions[1]);
		bboyR[i].fit(tri.positions[2]);
	}

	// calculate SAH
	float SAHcostminY = FLT_MAX;
	int bestIndexY = INT_MAX;
	for (int i = 0; i < obj_num; ++i) {
		float SAHcost = 2 * bbox.area() + (bboyL[i].area() * (i + 1) + bboyR[obj_num - i - 1].area() * (obj_num - i));
		if (SAHcost < SAHcostminY) {
			SAHcostminY = SAHcost;
			bestIndexY = i;
		}
	}

	// z-axis
	AABB* bbozL = new AABB[obj_num];
	AABB* bbozR = new AABB[obj_num];
	int* sorted_obj_index_z = new int[obj_num];
	this->sortAxis(obj_index, 2, 0, obj_num - 1);
	for (int i = 0; i < obj_num; ++i) sorted_obj_index_z[i] = obj_index[i];

	// sweep min to max
	bbozL[0].reset();
	bbozL[0].fit(triangleMesh->triangles[obj_index[0]].positions[0]);
	bbozL[0].fit(triangleMesh->triangles[obj_index[0]].positions[1]);
	bbozL[0].fit(triangleMesh->triangles[obj_index[0]].positions[2]);
	for (int i = 1; i < obj_num; ++i) {
		bbozL[i].reset();
		bbozL[i] = bbozL[i - 1];
		const Triangle& tri = triangleMesh->triangles[obj_index[i]];
		bbozL[i].fit(tri.positions[0]);
		bbozL[i].fit(tri.positions[1]);
		bbozL[i].fit(tri.positions[2]);
	}

	// sweep max to min
	bbozR[0].reset();
	bbozR[0].fit(triangleMesh->triangles[obj_index[obj_num - 1]].positions[0]);
	bbozR[0].fit(triangleMesh->triangles[obj_index[obj_num - 1]].positions[1]);
	bbozR[0].fit(triangleMesh->triangles[obj_index[obj_num - 1]].positions[2]);
	for (int i = 1; i < obj_num; ++i) {
		bbozR[i].reset();
		bbozR[i] = bbozR[i - 1];
		const Triangle& tri = triangleMesh->triangles[obj_index[obj_num - i - 1]];
		bbozR[i].fit(tri.positions[0]);
		bbozR[i].fit(tri.positions[1]);
		bbozR[i].fit(tri.positions[2]);
	}

	// calculate SAH
	float SAHcostminZ = FLT_MAX;
	int bestIndexZ = INT_MAX;
	for (int i = 0; i < obj_num; ++i) {
		float SAHcost = 2 * bbox.area() + (bbozL[i].area() * (i + 1) + bbozR[obj_num - i - 1].area() * (obj_num - i));
		if (SAHcost < SAHcostminZ) {
			SAHcostminZ = SAHcost;
			bestIndexZ = i;
		}
	}

	// calculate minimum SAH
	float SAHcostmin = FLT_MAX;
	int bestIndex = INT_MAX;
	AABB bestbboxL, bestbboxR;
	int* sorted_obj_index = new int[obj_num];
	if (SAHcostminX < SAHcostminY && SAHcostminX < SAHcostminZ) {
		SAHcostmin = SAHcostminX;
		bestIndex = bestIndexX;
		bestbboxL = bboxL[bestIndex];
		bestbboxR = bboxR[obj_num - bestIndex - 1];
		for (int i = 0; i < obj_num; ++i) sorted_obj_index[i] = sorted_obj_index_x[i];
	}
	else if (SAHcostminY < SAHcostminZ) {
		SAHcostmin = SAHcostminY;
		bestIndex = bestIndexY;
		bestbboxL = bboyL[bestIndex];
		bestbboxR = bboyR[obj_num - bestIndex - 1];
		for (int i = 0; i < obj_num; ++i) sorted_obj_index[i] = sorted_obj_index_y[i];
	}
	else {
		SAHcostmin = SAHcostminZ;
		bestIndex = bestIndexZ;
		bestbboxL = bbozL[bestIndex];
		bestbboxR = bbozR[obj_num - bestIndex - 1];
		for (int i = 0; i < obj_num; ++i) sorted_obj_index[i] = sorted_obj_index_z[i];
	}

	delete[] bboxL;
	delete[] bboxR;
	delete[] bboyL;
	delete[] bboyR;
	delete[] bbozL;
	delete[] bbozR;

	// stop subdividing
	if (obj_num * bbox.area() < SAHcostmin) {
		delete[] sorted_obj_index;
		this->nodeNum++;
		this->node[this->nodeNum - 1].bbox = bbox;
		this->node[this->nodeNum - 1].isLeaf = true;
		this->node[this->nodeNum - 1].triListNum = obj_num;
		this->node[this->nodeNum - 1].triList = new int[obj_num];
		for (int i = 0; i < obj_num; i++) this->node[this->nodeNum - 1].triList[i] = obj_index[i];
		int temp_id;
		temp_id = this->nodeNum - 1;
		this->leafNum++;
		return temp_id;
	}

#endif

	// leaf node
	if (obj_num <= 4) {
		delete[] sorted_obj_index;
		this->nodeNum++;
		this->node[this->nodeNum - 1].bbox = bbox;
		this->node[this->nodeNum - 1].isLeaf = true;
		this->node[this->nodeNum - 1].triListNum = obj_num;
		this->node[this->nodeNum - 1].triList = new int[obj_num];
		for (int i = 0; i < obj_num; i++) this->node[this->nodeNum - 1].triList[i] = obj_index[i];
		int temp_id;
		temp_id = this->nodeNum - 1;
		this->leafNum++;
		return temp_id;
	}

	// split obj_index into two
	int* obj_indexL = new int[bestIndex + 1];
	int* obj_indexR = new int[obj_num - (bestIndex + 1)];
	for (int i = 0; i <= bestIndex; ++i) obj_indexL[i] = sorted_obj_index[i];
	for (int i = bestIndex + 1; i < obj_num; ++i) obj_indexR[i - (bestIndex + 1)] = sorted_obj_index[i];
	delete[] sorted_obj_index;
	int obj_numL = bestIndex + 1;
	int obj_numR = obj_num - (bestIndex + 1);

	// recursive call to build a tree
	this->nodeNum++;
	int temp_id;
	temp_id = this->nodeNum - 1;
	this->node[temp_id].bbox = bbox;
	this->node[temp_id].isLeaf = false;
	this->node[temp_id].idLeft = splitBVH(obj_indexL, obj_numL, bestbboxL);
	this->node[temp_id].idRight = splitBVH(obj_indexR, obj_numR, bestbboxR);

	// clean up and return
	delete[] obj_indexL;
	delete[] obj_indexR;
	return temp_id;
}

// you may keep this part as-is
void BVH::build(const TriangleMesh* mesh) {
	triangleMesh = mesh;

	// construct the bounding volume hierarchy
	const int obj_num = (int)(triangleMesh->triangles.size());
	int* obj_index = new int[obj_num];
	for (int i = 0; i < obj_num; ++i) obj_index[i] = i;
	this->nodeNum = 0;
	this->node = new BVHNode[obj_num * 2];
	this->leafNum = 0;

	// calculate a scene bounding box
	AABB bbox;
	for (int i = 0; i < obj_num; i++) {
		const Triangle& tri = triangleMesh->triangles[obj_index[i]];
		bbox.fit(tri.positions[0]);
		bbox.fit(tri.positions[1]);
		bbox.fit(tri.positions[2]);
	}

	// ---------- building BVH ----------
	printf("Building BVH...\n");
	splitBVH(obj_index, obj_num, bbox);
	printf("Done.\n");

	delete[] obj_index;
}

// you may keep this part as-is
bool BVH::traverse(HitInfo& minHit, const Ray& ray, int node_id, float tMin, float tMax) const {
	bool hit = false;
	HitInfo tempMinHit, tempMinHitL, tempMinHitR;
	bool hit1, hit2;

	// check if node is a leaf
	if (this->node[node_id].isLeaf) {

		// for all triangles in the node's triangle list
		for (int i = 0; i < (this->node[node_id].triListNum); ++i) {

			// check if intersects with triangle i
			if (triangleMesh->raytraceTriangle(tempMinHit, ray, triangleMesh->triangles[this->node[node_id].triList[i]], tMin, tMax)) {
				hit = true;
				if (tempMinHit.t < minHit.t) minHit = tempMinHit;
			}
		}
	}

	// node is not a leaf
	else {

		// check for intersection with left and right children BVH nodes
		hit1 = this->node[this->node[node_id].idLeft].bbox.intersect(tempMinHitL, ray);
		hit2 = this->node[this->node[node_id].idRight].bbox.intersect(tempMinHitR, ray);

		// only if the hit is closer than the previous one
		hit1 = hit1 && (tempMinHitL.t < minHit.t);
		hit2 = hit2 && (tempMinHitR.t < minHit.t);

		// traversal of children nodes order based on which hit is closest
		if (hit1 && hit2) {
			if (tempMinHitL.t < tempMinHitR.t) {
				hit = traverse(minHit, ray, this->node[node_id].idLeft, tMin, tMax);
				hit = traverse(minHit, ray, this->node[node_id].idRight, tMin, tMax);
			} else {
				hit = traverse(minHit, ray, this->node[node_id].idRight, tMin, tMax);
				hit = traverse(minHit, ray, this->node[node_id].idLeft, tMin, tMax);
			}
		}
		else if (hit1) hit = traverse(minHit, ray, this->node[node_id].idLeft, tMin, tMax);
		else if (hit2) hit = traverse(minHit, ray, this->node[node_id].idRight, tMin, tMax);
	}

	// returns the closest hit
	return hit;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// spherical light source
class Sphere {

	public:

		float3 centre;
		float radius;
		float radiusInv;
		Material material;

		// constructor
		Sphere(float3 c, float r, Material m): centre(c), radius(r), radiusInv(1 / r), material(m) {}

		// check if ray intersects sphere
		bool intersect(HitInfo& hitInfo, const Ray& ray, float tMin, float tMax) const {

			// solve quadratic vector equation
			float3 oc = centre - ray.o;
			float b = dot(ray.d, oc);
			float c = dot(oc, oc) - radius * radius;
			float discriminant = b * b - c;

			// no intersection
			if (discriminant < 0) return false;

			// intersection
			discriminant = sqrtf(discriminant);

			// find minimum t
			float t = b - discriminant;
			if (t < tMin || tMax < t) {
				t = b + discriminant;
				if (t < tMin || tMax < t) {
					return false;
				}
			}

			// set hit info and return
			hitInfo.t = t;
			hitInfo.P = ray.o + hitInfo.t * ray.d;
			hitInfo.G = normalize(hitInfo.P - centre);
			hitInfo.material = &material;
			return true;
		}

		// generate random point on the sphere
		float3 sampleSurface(const float3& n) const {

			// using normalized 3d standard normal
			if (SAMPLE_UNIFORM_HEMISPHERE) {
				float3 d = float3(std_norm(gen), std_norm(gen), std_norm(gen));
				while (d.x == 0 && d.y == 0 && d.z == 0) d = float3(std_norm(gen), std_norm(gen), std_norm(gen));
				d = normalize(d);
				if (dot(n, d) < 0) d *= -1;
				return centre + radius * d;
			}

			// using malley's method
			else if (SAMPLE_COS_WEIGHTED_HEMISPHERE) {

				// sample cos weighted positive unit hemisphere
				float x = 2 * PCG32::rand() - 1;
				float y = 2 * PCG32::rand() - 1;
				if (x == 0 && y == 0) return float3(x, y, 1);
				float theta, r;
				if (abs(x) > abs(y)) {
					r = x;
					theta = PI_OVER_FOUR * (y / x);
				}
				else {
					r = y;
					theta = PI_OVER_TWO - PI_OVER_FOUR * (x / y);
				}
				float3 d = float3(r * cos(theta), r * sin(theta), sqrtf(1 - x * x - y * y));

				// build orthonormal basis with normal
				float sign = copysignf(1, n.z);
				const float a = -1 / (sign + n.z);
				const float b = n.x * n.y * a;
				const float3 b1 = float3(1 + sign * n.x * n.x * a, sign * b, -sign * n.x);
				const float3 b2 = float3(b, sign + n.y * n.y * a, -n.y);

				// centre about normal
				return centre + radius * (d.x * b1 + d.y * b2 + d.z * n);
			}

			// checking if this works ...
			else if (SAMPLE_UNIFORM_HEMISPHERE_LIKE_A_NOOB) {
				
				// using spherical coordinates
				const float z = PCG32::rand();
				const float r = sqrtf(1 - z * z);
				const float phi = 2 * PI * PCG32::rand();
				const float3 d = float3(r * std::cos(phi), r * std::sin(phi), z);

				// build orthonormal basis with normal
				float sign = copysignf(1, n.z);
				const float a = -1 / (sign + n.z);
				const float b = n.x * n.y * a;
				const float3 b1 = float3(1 + sign * n.x * n.x * a, sign * b, -sign * n.x);
				const float3 b2 = float3(b, sign + n.y * n.y * a, -n.y);

				// centre about normal
				return centre + radius * (d.x * b1 + d.y * b2 + d.z * n);
			}
		}

		// sample pdf
		float pdf(const float3& n, const float3& wi) const {

			// hemisphere subtends 2 Pi steradians
			if (SAMPLE_UNIFORM_HEMISPHERE) return ONE_OVER_TWO_PI * radiusInv * radiusInv;

			// inputs must be normalized
			else if (SAMPLE_COS_WEIGHTED_HEMISPHERE) return dot(n, wi) * ONE_OVER_PI * radiusInv * radiusInv;

			// checking if this works ...
			else if (SAMPLE_UNIFORM_HEMISPHERE_LIKE_A_NOOB) return ONE_OVER_TWO_PI * radiusInv * radiusInv;
		}

};

// forward declaration
static float3 pathShader(Ray ray);

// scene definition
class Scene {

	public:

		// scene components
		std::vector<TriangleMesh*> objects;
		std::vector<Sphere*> lights; // change this to list of geometry pointers
		std::vector<BVH> bvhs;

		// add object triangle mesh to scene
		void addObject(TriangleMesh* pObj) { objects.push_back(pObj); }

		// add spherical light source to scene
		void addLight(Sphere* pObj) { lights.push_back(pObj); }

		// compute BVH
		void preCalc() {
			bvhs.resize(objects.size());
			for (int i = 0, i_n = static_cast<int>(objects.size()); i < i_n; ++i) {
				objects[i]->preCalc();
				bvhs[i].build(objects[i]);
			}
		}

		// eye ray generation
		static Ray eyeRay(const int x, const int y) {

			// compute the camera coordinate system
			const float3 wDir = normalize(-globalViewDir); // back
			const float3 uDir = normalize(cross(globalUp, wDir)); // right
			const float3 vDir = cross(wDir, uDir); // up

			// compute the pixel location in world coordinate space using camera coordinate space
			const float imPlaneUPos = (static_cast<float>(x) + 0.5f) / static_cast<float>(globalWidth) - 0.5f;
			const float imPlaneVPos = (static_cast<float>(y) + 0.5f) / static_cast<float>(globalHeight) - 0.5f;
			const float3 pixelPos = globalEye + globalAspectRatio * globalFilmSize * imPlaneUPos * uDir + globalFilmSize * imPlaneVPos * vDir - globalDistanceToFilm * wDir;

			// trace ray through centre of each pixel
			return Ray{globalEye, normalize(pixelPos - globalEye)};
		}

		// ray-scene intersection
		bool intersect(HitInfo& minHit, const Ray& ray, const float tMin = 0.0f, const float tMax = FLT_MAX) const {
			bool hit = false;
			HitInfo tempMinHit;
			minHit.t = FLT_MAX;
			for (int i = 0, i_n = static_cast<int>(objects.size()); i < i_n; ++i) {
				if (bvhs[i].intersect(tempMinHit, ray, tMin, tMax)) {
					if (tempMinHit.t < minHit.t) {
						hit = true;
						minHit = tempMinHit;
					}
				}
			}
			for (int i = 0, i_n = static_cast<int>(lights.size()); i < i_n; ++i) {
				if (lights[i]->intersect(tempMinHit, ray, tMin, tMax)) {
					if (tempMinHit.t < minHit.t) {
						hit = true;
						minHit = tempMinHit;
					}
				}
			}
			return hit;
		}

		// path tracing
		void pathTrace(const int rootStrata) const {

			// division
			const float invStrata = 1.0f / static_cast<float>(rootStrata * rootStrata);
			const float strataHeight = 1.0f / static_cast<float>(globalHeight * rootStrata);
			const float strataWidth = 1.0f / static_cast<float>(globalWidth * rootStrata);

			// compute the camera coordinate system
			const float3 wDir = normalize(-globalViewDir); // back
			const float3 uDir = normalize(cross(globalUp, wDir)); // right
			const float3 vDir = cross(wDir, uDir); // up

			// accumulate shade
			float3 shade;

			// OpenMP
			#pragma omp parallel for schedule(dynamic, 1) private(shade)

			// loop over rows
			for (int j = 0; j < globalHeight; ++j) {

				// report progress
				fprintf(stderr, "\rRendering %5.2f%%", 100.0 * j / (globalHeight - 1));

				// loop over columns
				for (int i = 0; i < globalWidth; ++i) {

					// reset shade
					shade = float3(0.0f, 0.0f, 0.0f);

					// stratified subpixel sampling
					for (int y = 0; y < rootStrata; ++y) {
						for (int x = 0; x < rootStrata; ++x) {

							// generate random jitter
							const float Randy = PCG32::rand() * strataHeight;
							const float Randolf = PCG32::rand() * strataWidth;

							// compute coordinate of jittered sample
							const float v = static_cast<float>(j * rootStrata + y) + Randy;
							const float u = static_cast<float>(i * rootStrata + x) + Randolf;

							// compute sample location in world space
							const float imPlaneVPos = v * strataHeight - 0.5f;
							const float imPlaneUPos = u * strataWidth - 0.5f;
							const float3 pixelPos = globalEye + globalAspectRatio * globalFilmSize * imPlaneUPos * uDir + globalFilmSize * imPlaneVPos * vDir - globalDistanceToFilm * wDir;

							// trace ray through sample location
							shade += pathShader(Ray(globalEye, normalize(pixelPos - globalEye)));
						}
					}

					// average over strata samples
					globalImage.pixel(i, j) = shade * invStrata;
				}
			}
		}

};

// global scene object
static Scene globalScene;



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// path tracing shading
static float3 pathShader(Ray ray) {

	// radiance, throughput
	float3 radiance(0.0f), throughput(1.0f);

	// multiple importance sampling
	float probBRDF = 0.0f;
	float cosThetaMax;

	// hit specular material
	bool specular = false;

	// trace path
	int pathLength = 0;
	while (true) {

		// check intersection
		HitInfo hitInfo;
		if (globalScene.intersect(hitInfo, ray) == false) break;
		const float3 hitPoint = hitInfo.P + hitInfo.G * EPSILON;
		const float3 wo = -ray.d;
		++pathLength;

		// next event estimation
		const int k = 0;
		const Sphere* light = globalScene.lights[k];
		// const float3 lightNormal = normalize(hitPoint - light->centre);

		// hit emissive surface
		if (hitInfo.material->type == LIGHT) {

			// camera ray intersection or hit specular material
			if (pathLength == 1 || specular) radiance += throughput * hitInfo.material->emission;

			// multiple importance sampling
			else {
				// throughput *= probBRDF / (probBRDF + light->pdf(lightNormal, hitInfo.G));
				throughput *= probBRDF / (probBRDF + (1 / (2 * PI * (1 - cosThetaMax))));
				radiance += throughput * hitInfo.material->emission;
			}
		}





		// distance from centre of spherical light
		// float3 lightNormal = hitPoint - light->centre;
		float3 lightNormal = light->centre - hitPoint;
		float oneOverDistanceSquared = 1 / dot(lightNormal, lightNormal);
		lightNormal *= sqrtf(oneOverDistanceSquared);

		// random friends
		const float Bertrand = PCG32::rand();
		const float Randolf = PCG32::rand();

		// random direction in cone of spherical light
		// const float sinThetaMax2 = light->radius * light->radius * oneOverDistanceSquared;
		// const float sinThetaMax = sqrtf(sinThetaMax2);
		cosThetaMax = sqrtf(std::max(0.0f, 1 - light->radius * light->radius * oneOverDistanceSquared));
		const float cosTheta = 1 + (cosThetaMax - 1) * Bertrand;
		// const float sinTheta2 = 1 - cosTheta * cosTheta;
		const float sinTheta = sqrtf(1 - cosTheta * cosTheta);
		// const float cosAlpha = sinTheta2 / sinThetaMax + cosTheta * sqrtf(1 - sinTheta2 / sinThetaMax2);
		// const float sinAlpha = sqrtf(1 - cosAlpha * cosAlpha);
		const float phi = 2 * PI * Randolf;

		// const float cosThetaMax = sqrtf(1 - light->radius * light->radius * oneOverDistanceSquared);
		// const float cosTheta = 1 - Bertrand + Randolf * cosThetaMax;
		// const float sinTheta = sqrtf(1 - cosTheta * cosTheta);
		// const float phi = 2 * PI * Randolf;

		// float cos_theta = std::lerp(r1, cos_theta_max, 1.f);
		// float sin_theta = std::sqrt(1.f - cos_theta * cos_theta);
		// float phi = 2 * M_PI * r2;
		// return std::cos(phi) * sin_theta * x + std::sin(phi) * sin_theta * y + cos_theta * z;

		// build orthonormal basis
		float sign = copysignf(1, lightNormal.z);
		const float a = -1 / (sign + lightNormal.z);
		const float b = lightNormal.x * lightNormal.y * a;
		const float3 b1 = float3(1 + sign * lightNormal.x * lightNormal.x * a, sign * b, -sign * lightNormal.x);
		const float3 b2 = float3(b, sign + lightNormal.y * lightNormal.y * a, -lightNormal.y);

		// centre about normal
		// const float3 lightPointNormal = cosAlpha * lightNormal + sinAlpha * std::cos(phi) * b1 + sinAlpha * std::sin(phi) * b2;
		// const float3 lightPoint = light->centre + light->radius * lightPointNormal;
		// float3 wi = lightPoint - hitPoint;
		// oneOverDistanceSquared = 1 / dot(wi, wi);
		// wi *= sqrtf(oneOverDistanceSquared);

		const float3 wi = cosTheta * lightNormal + sinTheta * std::cos(phi) * b1 + sinTheta * std::sin(phi) * b2;

		// trace shadow ray from hit to light
		HitInfo shadowHitInfo;
		if (globalScene.intersect(shadowHitInfo, Ray(hitPoint, wi)) && shadowHitInfo.material->type == LIGHT) {
			
			const float probLight = 1 / (2 * PI * (1 - cosThetaMax));
			const float weight = probLight / (probLight + hitInfo.material->pdf(hitInfo.G, wi));

			// const float geometry = dot(hitInfo.G, wi) * dot(lightPointNormal, -wi) * oneOverDistanceSquared;
			// radiance += weight * throughput * hitInfo.material->spectrum() * light->material.emission * geometry / probLight;

			radiance += weight * throughput * hitInfo.material->spectrum() * light->material.emission * dot(hitInfo.G, wi) / probLight; // looks much better without geometry term
		}



/*
		// sample surface of light source
		float3 lightPoint = light->sampleSurface(lightNormal);
		float3 hitToLight = lightPoint - hitPoint;
		const float oneOverDistanceSquared = 1 / dot(hitToLight, hitToLight);
		const float oneOverDistance = sqrtf(oneOverDistanceSquared);
		hitToLight *= oneOverDistance;

		// trace shadow ray from hit to light
		HitInfo shadowHitInfo;
		if (globalScene.intersect(shadowHitInfo, Ray(hitPoint, hitToLight)) && shadowHitInfo.material->type == LIGHT) {

			// check if intersection with the sampled point
			if (dot(hitToLight, normalize(lightPoint - light->centre))) {
				const float probLight = light->pdf(lightNormal, shadowHitInfo.G);
				const float weight = probLight / (probLight + hitInfo.material->pdf(hitInfo.G, hitToLight));
				const float geometry = dot(hitInfo.G, hitToLight) * dot(shadowHitInfo.G, -hitToLight) * oneOverDistanceSquared;
				radiance += weight * throughput * hitInfo.material->spectrum() * light->material.emission * geometry / probLight;
			}
		}
*/



		// continue path
		ray = Ray(hitPoint, hitInfo.material->sampleDirection(wo, hitInfo.G));

		// update throughput
		const float probBRDF = hitInfo.material->pdf(hitInfo.G, ray.d);
		throughput *= hitInfo.material->spectrum() * dot(hitInfo.G, ray.d) / probBRDF;

		// russian roulette
		if (pathLength > MAXIMUM_PATH_LENGTH) {
			float probabilityOfContinuing = std::max(throughput.x, std::max(throughput.y, throughput.z));
			if (PCG32::rand() < probabilityOfContinuing) throughput /= probabilityOfContinuing;
			else break;
		}
	}

	// return radiance
	return radiance;
}



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
