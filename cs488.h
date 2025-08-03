///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////   CS 488 base code                 /////////////////////////////////////////////////////////////////////////
//////////   (written by Toshiya Hachisuka)   /////////////////////////////////////////////////////////////////////////
//////////   CS 488 final project             /////////////////////////////////////////////////////////////////////////
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
#include <cmath>
#include <climits>

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

// window size and resolution
constexpr int globalWidth = 640; // 512
constexpr int globalHeight = 480; // 384

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

// global constants
constexpr int MINIMUM_PATH_LENGTH = 0;

// refractive indices
constexpr float REFR_IND_AIR = 1.00029f;
constexpr float REFR_IND_GLASS = 1.458f;
constexpr float AIR_TO_GLASS = REFR_IND_AIR / REFR_IND_GLASS;
constexpr float GLASS_TO_AIR = REFR_IND_GLASS / REFR_IND_AIR;
constexpr float AIR_GLASS_R = (REFR_IND_AIR - REFR_IND_GLASS) * (REFR_IND_AIR - REFR_IND_GLASS) / ((REFR_IND_AIR + REFR_IND_GLASS) * (REFR_IND_AIR + REFR_IND_GLASS));

// homogeneous volumetric scattering
constexpr float ABSORPTION = 0.005f; // cornellbox: 0.0, mirrodin: 0.01
constexpr float SCATTERING = 0.025f; // cornellbox: 1.0, mirrodin: 0.05
constexpr float ATTENUATION = ABSORPTION + SCATTERING;
constexpr float ONE_OVER_ATTENUATION = 1.0f / ATTENUATION;

// switches
constexpr bool ATMOSPHERIC_SCATTERING = true;

// choose a scene to render
constexpr bool RENDER_CORNELL_BOX = false;
constexpr bool RENDER_GLASS_DIORAMA = false;
constexpr bool RENDER_MICROFACET_DIORAMA = false;
constexpr bool RENDER_FINAL_SCENE = false;

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
	LIGHT,
	LAMBERTIAN,
	METAL,
	GLASS,
	MICROFACET
};

// uber material
class Material {

	public:

		// fields
		std::string name;
		MaterialType type = LAMBERTIAN;
		float3 emission = float3(0.0f);

		// colour
		float3 Ka = float3(0.0f); // ambient colour
		float3 Kd = float3(0.9f); // diffuse colour
		float3 Ks = float3(0.0f); // specular colour
		float Ns = 0.0f; // specular exponent

		// microfacet
		float roughness = 0.25f;
		float alpha = roughness * roughness;
		float alpha2 = alpha * alpha;
		float k = alpha / 2.0f;
		void setRoughness(const float& r) {
			roughness = r;
			alpha = roughness * roughness;
			alpha2 = alpha * alpha;
			k = alpha / 2.0f;
		}

		// ctor & dtor
		Material() = default;
		~Material() = default;

		// 8-bit texture
		bool isTextured = false;
		unsigned char* texture = nullptr;
		int textureWidth = 0;
		int textureHeight = 0;
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

		// brdf value
		float3 brdf(const float3& wo, const float3& n, const float3& wi) const {

			// perfect diffuse
			if (type == LIGHT || type == LAMBERTIAN) return Kd * ONE_OVER_PI;

			// microfacet reflection
			else if (type == MICROFACET) {

				// calculate half direction
				const float3 m = normalize(wo + wi);

				// microfacet distribution function
				const float nom = std::max(0.0f, dot(n, m));
				const float base = nom * nom * (alpha2 - 1) + 1;
				const float base2 = base * base;
				float D;
				if (base2 > 0) D = alpha2 * ONE_OVER_PI / base2;
				else return float3(0.0f);

				// evaluate fresnel term
				const float c = 1 - dot(m, wo);
				float3 F = Kd + (1 - Kd) * c * c * c * c * c;

				// G1 wo
				const float noo = std::max(0.0f, dot(n, wo));
				const float denomo = (noo * (1 - k) + k);
				float G1o;
				if (denomo > 0) G1o = noo / denomo;
				else return float3(0.0f);

				// G1 wi
				const float noi = std::max(0.0f, dot(n, wi));
				const float denomi = (noi * (1 - k) + k);
				float G1i;
				if (denomi > 0) G1i = noi / denomi;
				else return float3(0.0f);

				// masking shadowing function
				const float G = G1o * G1i;

				// macrosurface brdf
				const float denominator = noo * noi;
				if (denominator > 0) return D * F * G / (4 * denominator);
				else return float3(0.0f);
			}

			// default
			else return float3(0.0f);
		}

		// sample direction
		float3 sample(const float3& wo, const float3& n) const {

			// perfect diffuse reflection
			if (type == LIGHT || type == LAMBERTIAN) {

				// sample unit disk
				float x = 2 * PCG32::rand() - 1;
				float y = 2 * PCG32::rand() - 1;
				if (x == 0 && y == 0) return float3(x, y, 1);
				float theta, r;
				if (std::abs(x) > std::abs(y)) {
					r = x;
					theta = PI_OVER_FOUR * (y / x);
				}
				else {
					r = y;
					theta = PI_OVER_TWO - PI_OVER_FOUR * (x / y);
				}
				x = r * std::cos(theta);
				y = r * std::sin(theta);

				// project to hemisphere safely
				float z = 1 - x * x - y * y;
				if (z < EPSILON) z = 0;
				else z = sqrtf(z);

				// build orthonormal basis with normal
				const float sign = copysignf(1, n.z);
				const float a = -1 / (sign + n.z);
				const float b = n.x * n.y * a;
				const float3 b1 = float3(1 + sign * n.x * n.x * a, sign * b, -sign * n.x);
				const float3 b2 = float3(b, sign + n.y * n.y * a, -n.y);

				// centre about normal
				return x * b1 + y * b2 + z * n;
			}

			// microfacet reflection
			else if (type == MICROFACET) {

				// random friends
				const float Bertrand = PCG32::rand();
				const float Randolf = PCG32::rand();

				// sample microfacet normal
				const float theta = std::atan((alpha * sqrtf(Bertrand)) / sqrtf(1 - Bertrand));
				const float phi = 2 * PI * Randolf;
				const float x = std::sin(theta) * std::cos(phi);
				const float y = std::sin(theta) * std::sin(phi);
				const float z = std::cos(theta);

				// build orthonormal basis with normal
				const float sign = copysignf(1, n.z);
				const float a = -1 / (sign + n.z);
				const float b = n.x * n.y * a;
				const float3 b1 = float3(1 + sign * n.x * n.x * a, sign * b, -sign * n.x);
				const float3 b2 = float3(b, sign + n.y * n.y * a, -n.y);

				// centre about macrosurface normal
				const float3 m = x * b1 + y * b2 + z * n;

				// return reflected direction
				return 2 * dot(m, wo) * m - wo;
			}

			// default
			else return float3(0.0f);
		}

		// pdf value at sample
		float pdf(const float3& wo, const float3& n, const float3& wi) const {

			// perfect diffuse reflection
			if (type == LIGHT || type == LAMBERTIAN) return std::max(0.0f, dot(n, wi)) * ONE_OVER_PI;

			// microfacet reflection
			else if (type == MICROFACET) {

				// calculate half direction
				const float3 m = normalize(wo + wi);

				// microfacet distribution function
				const float nom = std::max(0.0f, dot(n, m));
				const float base = nom * nom * (alpha2 - 1) + 1;
				const float base2 = base * base;
				float D;
				if (base2 > 0) D = alpha2 * ONE_OVER_PI / base2;
				else return 0.0f;

				// return pdf value
				const float moi = std::max(0.0f, dot(m, wi));
				if (moi > 0) return D * nom / (4 * dot(m, wi));
				else return 0.0f;
			}

			// default
			else return 0.0f;
		}

		// throughput weight
		float3 weight(const float3& wo, const float3& n, const float3& wi) const {

			// perfect diffuse
			if (type == LIGHT || type == LAMBERTIAN) return Kd;

			// microfacet reflection
			else if (type == MICROFACET) {

				// calculate half direction
				const float3 m = normalize(wo + wi);
				const float nom = std::max(0.0f, dot(n, m));
				const float moi = std::max(0.0f, dot(m, wi));

				// evaluate fresnel term
				// const float c = 1 - dot(m, wo);
				// float3 F;
				// if (metallic) F = Kd + (1 - Kd) * c * c * c * c * c;
				// else {
				// 	const float F0 = 0.16f * (reflectance * reflectance);
				// 	F = float3(F0 + (1 - F0) * c * c * c * c * c);
				// }

				// G1 wo
				const float noo = std::max(0.0f, dot(n, wo));
				const float denomo = (noo * (1 - k) + k);
				float G1o;
				if (denomo > 0) G1o = noo / denomo;
				else return float3(0.0f);

				// G1 wi
				const float noi = std::max(0.0f, dot(n, wi));
				const float denomi = (noi * (1 - k) + k);
				float G1i;
				if (denomi > 0) G1i = noi / denomi;
				else return float3(0.0f);

				// masking shadowing function
				const float G = G1o * G1i;

				// macrosurface brdf
				const float denominator = nom * noi;
				if (denominator > 0) return Kd * ONE_OVER_PI + G * moi / denominator;
				// if (denominator > 0) return float3(1.0f) * G * moi / denominator;
				// if (denominator > 0) return (1 - F) * Kd * ONE_OVER_PI + G * moi / denominator;
				else return float3(0.0f);
			}

			// default
			else return float3(0.0f);
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

		float t;
		float3 P;
		float3 G;
		float3 N;
		float2 T;
		const Material* material;
		bool inside;

};

// axis aligned bounding box
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

			// interpolate only after checking
			const float beta = detAbeta / detA;
			const float gamma = detAgamma / detA;
			const float t = detAt / detA;
			const float alpha = 1.0f - beta - gamma;

			// fill in result for hit
			result.t = t;
			result.P = alpha * tri.positions[0] + beta * tri.positions[1] + gamma * tri.positions[2];
			result.N = normalize(alpha * tri.normals[0] + beta * tri.normals[1] + gamma * tri.normals[2]);
			result.T = alpha * tri.texcoords[0] + beta * tri.texcoords[1] + gamma * tri.texcoords[2];
			result.material = &materials[tri.idMaterial];

			// geometric normal for reflection
			// float3 g = normalize(cross(ba, ca));
			// if (dot(g, -ray.d) < 0) g *= -1;
			// result.G = g;

			// geometric normal for refraction
			float3 g = normalize(cross(ba, ca));
			if (dot(result.N, g) < 0) g *= -1;
			result.G = g;

			// inside object for refraction
			if (dot(result.G, -ray.d) < 0) result.inside = true;
			else result.inside = false;

			// // calculate geometric normal to front of triangle
			// float3 g = normalize(cross(ba, ca));
			// // assuming the interpolated shading normal is pointing towards the front of the triangle
			// if (dot(g, result.N) < 0) g *= -1;
			// result.G = g;

			// // check if hitting front of triangle
			// if (dot(-ray.d, result.G) > 0) result.front = true;
			// else result.front = false;

			// return true for hit
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

			printf("Loading \"%s\" ...\n", filename);
			ParseOBJ(filename, nVertices, &vertices, &normals, &texcoords, nIndices, &indices, &matid);
			if (nVertices == 0) return false;
			this->triangles.resize(nIndices / 3);



			/*

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
			}
			else {
				// use default Lambertian
				this->materials.resize(1);
			}

			*/

			// convert .mlt data into BSDF definitions
			if (matid != nullptr) {

				// for all materials
				for (unsigned int i = 0; i < materials.size(); i++) {

					// default
					materials[i].type = LAMBERTIAN;

					// metal
					if (materials[i].Ns == 100.0f) materials[i].type = METAL;

					// glass
					if (materials[i].name.compare(0, 5, "glass", 0, 5) == 0) materials[i].type = GLASS;

					// microfacet
					if (materials[i].name.compare(0, 8, "titanium", 0, 8) == 0) {
						materials[i].type = MICROFACET;
						materials[i].setRoughness(0.1f);
					}
					if (materials[i].name.compare(0, 6, "chrome", 0, 6) == 0) {
						materials[i].type = MICROFACET;
						materials[i].setRoughness(0.1f);
					}
					if (materials[i].name.compare(0, 4, "iron", 0, 4) == 0) {
						materials[i].type = MICROFACET;
						materials[i].setRoughness(0.1f);
					}
					if (materials[i].name.compare(0, 6, "nickel", 0, 6) == 0) {
						materials[i].type = MICROFACET;
						materials[i].setRoughness(0.1f);
					}
					if (materials[i].name.compare(0, 8, "platinum", 0, 8) == 0) {
						materials[i].type = MICROFACET;
						materials[i].setRoughness(0.1f);
					}
					if (materials[i].name.compare(0, 6, "copper", 0, 6) == 0) {
						materials[i].type = MICROFACET;
						materials[i].setRoughness(0.1f);
					}
					if (materials[i].name.compare(0, 9, "palladium", 0, 9) == 0) {
						materials[i].type = MICROFACET;
						materials[i].setRoughness(0.1f);
					}
					if (materials[i].name.compare(0, 4, "zinc", 0, 4) == 0) {
						materials[i].type = MICROFACET;
						materials[i].setRoughness(0.1f);
					}
					if (materials[i].name.compare(0, 4, "gold", 0, 4) == 0) {
						materials[i].type = MICROFACET;
						materials[i].setRoughness(0.1f);
					}
					if (materials[i].name.compare(0, 8, "aluminum", 0, 8) == 0) {
						materials[i].type = MICROFACET;
						materials[i].setRoughness(0.1f);
					}
					if (materials[i].name.compare(0, 6, "silver", 0, 6) == 0) {
						materials[i].type = MICROFACET;
						materials[i].setRoughness(0.1f);
					}
				}
			}

			// default
			else this->materials.resize(1);



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
			printf("Cannot open \"%s\" for reading.\n", fileName);
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

#define SAHBVH // use this in once you have SAH-BVH
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
	printf("Building BVH ...\n");
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



// sphere
class Sphere {

	public:

		float3 centre;
		float radius;
		Material material;

		// constructor
		Sphere(float3 c, float r, Material m): centre(c), radius(r), material(m) {}

		// check if ray intersects sphere
		bool intersect(HitInfo& hitInfo, const Ray& ray, float tMin, float tMax) const {

			// find distance to intersection
			const float3 oc = ray.o - centre;
			const float b = dot(ray.d, oc);
			const float c = dot(oc, oc) - radius * radius;
			float discriminant = b * b - c;
			if (discriminant < 0) return false;
			discriminant = sqrtf(discriminant);
			const float t1 = -b + discriminant;
			const float t2 = -b - discriminant;
			float t;
			if (t1 <= 0 && t2 <= 0) return false;
			else if (t1 <= 0 || t2 <= 0) t = std::max(t1, t2);
			else t = std::min(t1, t2);

			// hit info
			hitInfo.t = t;
			hitInfo.P = ray.o + hitInfo.t * ray.d;
			hitInfo.G = normalize(hitInfo.P - centre);
			hitInfo.material = &material;
			if (dot(hitInfo.G, -ray.d) < 0) hitInfo.inside = true;
			else hitInfo.inside = false;
			return true;
		}

};

// forward declaration
static float3 pathShader(Ray ray);
static float3 volumeShader(Ray ray);

// scene
class Scene {

	public:

		// scene components
		std::vector<Sphere*> balls;
		std::vector<Sphere*> lights;
		std::vector<TriangleMesh*> objects;
		std::vector<BVH> bvhs;

		// add spherical light source to scene
		void addLight(Sphere* pObj) { lights.push_back(pObj); }

		// add ball to scene
		void addBall(Sphere* pObj) { balls.push_back(pObj); }

		// add object triangle mesh to scene
		void addObject(TriangleMesh* pObj) { objects.push_back(pObj); }

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
			for (int i = 0, i_n = static_cast<int>(lights.size()); i < i_n; ++i) {
				if (lights[i]->intersect(tempMinHit, ray, tMin, tMax)) {
					if (tempMinHit.t < minHit.t) {
						hit = true;
						minHit = tempMinHit;
					}
				}
			}
			for (int i = 0, i_n = static_cast<int>(balls.size()); i < i_n; ++i) {
				if (balls[i]->intersect(tempMinHit, ray, tMin, tMax)) {
					if (tempMinHit.t < minHit.t) {
						hit = true;
						minHit = tempMinHit;
					}
				}
			}
			for (int i = 0, i_n = static_cast<int>(objects.size()); i < i_n; ++i) {
				if (bvhs[i].intersect(tempMinHit, ray, tMin, tMax)) {
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
							if (ATMOSPHERIC_SCATTERING) shade += volumeShader(Ray(globalEye, normalize(pixelPos - globalEye)));
							else shade += pathShader(Ray(globalEye, normalize(pixelPos - globalEye)));
						}
					}

					// average over strata samples
					const float3 pixelValue = shade * invStrata;
					if (std::isnan(pixelValue.x) || std::isnan(pixelValue.y) || std::isnan(pixelValue.z)) {
						globalImage.pixel(i, j) = float3(1.0f, 0.0f, 1.0f); // return pink
						std::cout << "we got NaaN @ pixel (" << i << ", " << j << ") !! :^(" << std::endl;
					}
					else globalImage.pixel(i, j) = pixelValue;
				}
			}
		}

};

// global scene object
static Scene globalScene;

// path tracing with surface scattering
static float3 pathShader(Ray ray) {

	// radiance, throughput
	float3 radiance(0.0f), throughput(1.0f);

	// multiple importance sampling
	const int numberOfLights = static_cast<int>(globalScene.lights.size());
	const float probLight = 1.0f / static_cast<float>(numberOfLights);
	float probBRDF(0.0f), probNEE(0.0f);

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

		// hit emissive material
		if (hitInfo.material->type == LIGHT) {

			// camera ray intersection or hit specular material
			if (pathLength == 1 || specular) {
				if (RENDER_FINAL_SCENE) radiance += throughput * hitInfo.material->Kd * dot(hitInfo.G, wo);
				else radiance += throughput * hitInfo.material->emission * dot(hitInfo.G, wo);
			}

			// multiple importance sampling
			else radiance += (probBRDF / (probBRDF + probNEE)) * throughput * hitInfo.material->emission * dot(hitInfo.G, wo);
		}

		// specular metal
		else if (hitInfo.material->type == METAL) {
			const float3 reflection = normalize(-wo + 2 * dot(hitInfo.G, wo) * hitInfo.G);
			ray = Ray(hitPoint, reflection);
			throughput *= hitInfo.material->Ks;
			specular = true;
			continue;
		}

		// specular glass
		else if (hitInfo.material->type == GLASS) {

			// continue path
			float3 origin, direction;
			const float Randerson = PCG32::rand();

			// outside air to glass
			if (!hitInfo.inside) {

				// reflection coefficient
				const float c = 1 - dot(hitInfo.G, -ray.d);
				const float R = AIR_GLASS_R + (1 - AIR_GLASS_R) * c * c * c * c * c;

				// reflect
				if (Randerson < R) {
					origin = hitInfo.P + EPSILON * hitInfo.G;
					direction = normalize(ray.d - 2 * dot(ray.d, hitInfo.G) * hitInfo.G);
					if (dot(hitInfo.G, direction) < 0) direction -= 2 * dot(hitInfo.G, direction) * hitInfo.G;
				}

				// refract
				else {
					origin = hitInfo.P - EPSILON * hitInfo.G;
					const float radicand = 1 - AIR_TO_GLASS * AIR_TO_GLASS * (1 - dot(ray.d, hitInfo.G) * dot(ray.d, hitInfo.G));
					direction = normalize(AIR_TO_GLASS * (ray.d - dot(ray.d, hitInfo.G) * hitInfo.G) - sqrtf(radicand) * hitInfo.G);
					if (dot(hitInfo.G, direction) > 0) direction -= 2 * dot(hitInfo.G, direction) * hitInfo.G;
				}
			}

			// inside glass to air
			else {

				// reflection coefficient
				const float c = 1 + dot(hitInfo.G, -ray.d);
				const float R = AIR_GLASS_R + (1 - AIR_GLASS_R) * c * c * c * c * c;

				// check for total internal reflection
				const float radicand = 1 - GLASS_TO_AIR * GLASS_TO_AIR * (1 - dot(ray.d, hitInfo.G) * dot(ray.d, hitInfo.G));

				// total internal reflection
				if (radicand < 0) {
					origin = hitInfo.P - EPSILON * hitInfo.G;
					direction = normalize(ray.d - 2 * dot(ray.d, hitInfo.G) * hitInfo.G);
					if (dot(hitInfo.G, direction) > 0) direction -= 2 * dot(hitInfo.G, direction) * hitInfo.G;
				}

				// reflect
				else if (Randerson < R) {
					origin = hitInfo.P - EPSILON * hitInfo.G;
					direction = normalize(ray.d - 2 * dot(ray.d, hitInfo.G) * hitInfo.G);
					if (dot(hitInfo.G, direction) > 0) direction -= 2 * dot(hitInfo.G, direction) * hitInfo.G;
				}

				// refract
				else {
					origin = hitInfo.P + EPSILON * hitInfo.G;
					direction = normalize(GLASS_TO_AIR * (ray.d - dot(ray.d, hitInfo.G) * hitInfo.G) - sqrtf(radicand) * hitInfo.G);
					if (dot(hitInfo.G, direction) < 0) direction -= 2 * dot(hitInfo.G, direction) * hitInfo.G;
				}
			}

			// continue path
			ray = Ray(origin, direction);
			specular = true;
			continue;
		}

		// choose light source with equal probability
		int lightIndex = std::floorf(PCG32::rand() * numberOfLights);
		if (lightIndex < 0) lightIndex = 0;
		else if (numberOfLights - 1 < lightIndex) lightIndex = numberOfLights - 1;
		const Sphere* light = globalScene.lights[lightIndex];
		
		// direction and distance to centre of light
		float3 hitToLight = light->centre - hitPoint;
		const float oneOverDistanceSquared = 1 / dot(hitToLight, hitToLight);
		hitToLight *= sqrtf(oneOverDistanceSquared);

		// calculate random direction towards visible spherical cap
		const float cosThetaMax = sqrtf(std::max(0.0f, 1 - light->radius * light->radius * oneOverDistanceSquared));
		const float cosTheta = 1 + (cosThetaMax - 1) * PCG32::rand();
		const float sinTheta = sqrtf(1 - cosTheta * cosTheta);
		const float phi = 2 * PI * PCG32::rand();

		// build orthonormal basis
		const float sign = copysignf(1, hitToLight.z);
		const float a = -1 / (sign + hitToLight.z);
		const float b = hitToLight.x * hitToLight.y * a;
		const float3 b1 = float3(1 + sign * hitToLight.x * hitToLight.x * a, sign * b, -sign * hitToLight.x);
		const float3 b2 = float3(b, sign + hitToLight.y * hitToLight.y * a, -hitToLight.y);

		// centre about direction to centre of light
		const float3 wi = cosTheta * hitToLight + sinTheta * std::cos(phi) * b1 + sinTheta * std::sin(phi) * b2;

		// calculate probability
		probNEE = probLight / (2 * PI * (1 - cosThetaMax));

		// check visibility
		HitInfo shadowHitInfo;
		if (globalScene.intersect(shadowHitInfo, Ray(hitPoint, wi)) && shadowHitInfo.material->type == LIGHT) {
			probBRDF = hitInfo.material->pdf(wo, hitInfo.G, wi);
			const float weight = 1 / (probNEE + probBRDF);
			if (probBRDF > 0) radiance += weight * throughput * hitInfo.material->brdf(wo, hitInfo.G, wi) * light->material.emission * dot(hitInfo.G, wi);
		}

		// continue path and update throughput
		ray = Ray(hitPoint, hitInfo.material->sample(wo, hitInfo.G));
		probBRDF = hitInfo.material->pdf(wo, hitInfo.G, ray.d);
		if (probBRDF > 0) throughput *= hitInfo.material->brdf(wo, hitInfo.G, ray.d) * dot(hitInfo.G, ray.d) / probBRDF;
		specular = false;

		// russian roulette
		if (pathLength > MINIMUM_PATH_LENGTH) {
			float probabilityOfContinuing = std::max(throughput.x, std::max(throughput.y, throughput.z));
			if (probabilityOfContinuing < 1) {
				probabilityOfContinuing = std::max(0.25f, probabilityOfContinuing);
				if (PCG32::rand() < probabilityOfContinuing) throughput /= probabilityOfContinuing;
				else break;
			}
		}
	}

	// return radiance
	return radiance;
}

// path tracing with volume scattering
static float3 volumeShader(Ray ray) {

	// radiance, throughput
	float3 radiance(0.0f), throughput(1.0f);

	// multiple importance sampling
	const int numberOfLights = static_cast<int>(globalScene.lights.size());
	const float probLight = 1.0f / static_cast<float>(numberOfLights);
	float probBRDF(0.0f), probNEE(0.0f);

	// hit specular material
	bool specular = false;

	// trace path
	int pathLength = 0;
	while (true) {

		// free flight sampling
		const float s = -std::log(1 - PCG32::rand()) * ONE_OVER_ATTENUATION;

		// check intersection
		HitInfo hitInfo;

		// into the void
		if (!globalScene.intersect(hitInfo, ray) && RENDER_FINAL_SCENE) {
			hitInfo.t = FLT_MAX;
			hitInfo.inside = false;
		}
		else if (!globalScene.intersect(hitInfo, ray)) break;

		// volume interaction
		if (s < hitInfo.t && !hitInfo.inside) {

			// point within volume
			const float3 pointInVolume = ray.o + s * ray.d;
			++pathLength;

			// weight throughput
			throughput *= ONE_OVER_ATTENUATION * SCATTERING;

			// choose light source with equal probability
			int lightIndex = std::floorf(PCG32::rand() * numberOfLights);
			if (lightIndex < 0) lightIndex = 0;
			else if (numberOfLights - 1 < lightIndex) lightIndex = numberOfLights - 1;
			const Sphere* light = globalScene.lights[lightIndex];
			
			// direction and distance to centre of light
			float3 hitToLight = light->centre - pointInVolume;
			const float oneOverDistanceSquared = 1 / dot(hitToLight, hitToLight);
			hitToLight *= sqrtf(oneOverDistanceSquared);

			// calculate random direction towards visible spherical cap
			const float cosThetaMax = sqrtf(std::max(0.0f, 1 - light->radius * light->radius * oneOverDistanceSquared));
			const float cosTheta = 1 + (cosThetaMax - 1) * PCG32::rand();
			const float sinTheta = sqrtf(1 - cosTheta * cosTheta);
			const float phi = 2 * PI * PCG32::rand();

			// build orthonormal basis
			const float sign = copysignf(1, hitToLight.z);
			const float a = -1 / (sign + hitToLight.z);
			const float b = hitToLight.x * hitToLight.y * a;
			const float3 b1 = float3(1 + sign * hitToLight.x * hitToLight.x * a, sign * b, -sign * hitToLight.x);
			const float3 b2 = float3(b, sign + hitToLight.y * hitToLight.y * a, -hitToLight.y);

			// centre about direction to centre of light
			const float3 wi = cosTheta * hitToLight + sinTheta * std::cos(phi) * b1 + sinTheta * std::sin(phi) * b2;

			// calculate probability
			probNEE = probLight / (2 * PI * (1 - cosThetaMax));

			// check visibility
			HitInfo shadowHitInfo;
			if (globalScene.intersect(shadowHitInfo, Ray(pointInVolume, wi)) && shadowHitInfo.material->type == LIGHT) {

				// account for possible scattering
				const float transmittance = std::exp(-ATTENUATION * shadowHitInfo.t);

				// multiple importance sampling
				const float probDistance = transmittance;
				probBRDF = probDistance * ONE_OVER_FOUR_PI;
				const float weight = 1 / (probNEE + probBRDF);

				// accumulate radiance
				radiance += weight * throughput * transmittance * SCATTERING * ONE_OVER_FOUR_PI * light->material.emission;
			}

			// random scattered direction on unit sphere
			const float z = PCG32::rand();
			const float r = sqrtf(1 - z * z);
			const float butt = 2 * PI * PCG32::rand();
			ray = Ray(pointInVolume, float3(r * std::cos(butt), r * std::sin(butt), z));
			probBRDF = ONE_OVER_FOUR_PI;
		}

		// surface interaction
		else {

			// doink the hit point
			const float3 hitPoint = hitInfo.P + hitInfo.G * EPSILON;
			const float3 wo = -ray.d;
			++pathLength;

			// hit emissive material
			if (hitInfo.material->type == LIGHT) {

				// camera ray intersection or hit specular material
				if (pathLength == 1 || specular) {
					if (RENDER_FINAL_SCENE) radiance += throughput * hitInfo.material->Kd * dot(hitInfo.G, wo);
					else radiance += throughput * hitInfo.material->emission * dot(hitInfo.G, wo);
				}

				// multiple importance sampling
				else {

					// outside
					if (!hitInfo.inside) {
						const float probDistance = std::exp(-ATTENUATION * hitInfo.t);
						probBRDF *= probDistance;
					}
					
					// accumulate radiance
					radiance += (probBRDF / (probBRDF + probNEE)) * throughput * hitInfo.material->emission * dot(hitInfo.G, wo);
				}
			}

			// specular metal
			else if (hitInfo.material->type == METAL) {
				const float3 reflection = normalize(-wo + 2 * dot(hitInfo.G, wo) * hitInfo.G);
				ray = Ray(hitPoint, reflection);
				throughput *= hitInfo.material->Ks;
				specular = true;
				continue;
			}

			// specular glass
			else if (hitInfo.material->type == GLASS) {

				// continue path
				float3 origin, direction;
				const float Randerson = PCG32::rand();

				// outside air to glass
				if (!hitInfo.inside) {

					// reflection coefficient
					const float c = 1 - dot(hitInfo.G, -ray.d);
					const float R = AIR_GLASS_R + (1 - AIR_GLASS_R) * c * c * c * c * c;

					// reflect
					if (Randerson < R) {
						origin = hitInfo.P + EPSILON * hitInfo.G;
						direction = normalize(ray.d - 2 * dot(ray.d, hitInfo.G) * hitInfo.G);
						if (dot(hitInfo.G, direction) < 0) direction -= 2 * dot(hitInfo.G, direction) * hitInfo.G;
					}

					// refract
					else {
						origin = hitInfo.P - EPSILON * hitInfo.G;
						const float radicand = 1 - AIR_TO_GLASS * AIR_TO_GLASS * (1 - dot(ray.d, hitInfo.G) * dot(ray.d, hitInfo.G));
						direction = normalize(AIR_TO_GLASS * (ray.d - dot(ray.d, hitInfo.G) * hitInfo.G) - sqrtf(radicand) * hitInfo.G);
						if (dot(hitInfo.G, direction) > 0) direction -= 2 * dot(hitInfo.G, direction) * hitInfo.G;
					}
				}

				// inside glass to air
				else {

					// reflection coefficient
					const float c = 1 + dot(hitInfo.G, -ray.d);
					const float R = AIR_GLASS_R + (1 - AIR_GLASS_R) * c * c * c * c * c;

					// check for total internal reflection
					const float radicand = 1 - GLASS_TO_AIR * GLASS_TO_AIR * (1 - dot(ray.d, hitInfo.G) * dot(ray.d, hitInfo.G));

					// total internal reflection
					if (radicand < 0) {
						origin = hitInfo.P - EPSILON * hitInfo.G;
						direction = normalize(ray.d - 2 * dot(ray.d, hitInfo.G) * hitInfo.G);
						if (dot(hitInfo.G, direction) > 0) direction -= 2 * dot(hitInfo.G, direction) * hitInfo.G;
					}

					// reflect
					else if (Randerson < R) {
						origin = hitInfo.P - EPSILON * hitInfo.G;
						direction = normalize(ray.d - 2 * dot(ray.d, hitInfo.G) * hitInfo.G);
						if (dot(hitInfo.G, direction) > 0) direction -= 2 * dot(hitInfo.G, direction) * hitInfo.G;
					}

					// refract
					else {
						origin = hitInfo.P + EPSILON * hitInfo.G;
						direction = normalize(GLASS_TO_AIR * (ray.d - dot(ray.d, hitInfo.G) * hitInfo.G) - sqrtf(radicand) * hitInfo.G);
						if (dot(hitInfo.G, direction) < 0) direction -= 2 * dot(hitInfo.G, direction) * hitInfo.G;
					}
				}

				// continue path
				ray = Ray(origin, direction);
				specular = true;
				continue;
			}

			// choose light source with equal probability
			int lightIndex = std::floorf(PCG32::rand() * numberOfLights);
			if (lightIndex < 0) lightIndex = 0;
			else if (numberOfLights - 1 < lightIndex) lightIndex = numberOfLights - 1;
			const Sphere* light = globalScene.lights[lightIndex];
			
			// direction and distance to centre of light
			float3 hitToLight = light->centre - hitPoint;
			const float oneOverDistanceSquared = 1 / dot(hitToLight, hitToLight);
			hitToLight *= sqrtf(oneOverDistanceSquared);

			// calculate random direction towards visible spherical cap
			const float cosThetaMax = sqrtf(std::max(0.0f, 1 - light->radius * light->radius * oneOverDistanceSquared));
			const float cosTheta = 1 + (cosThetaMax - 1) * PCG32::rand();
			const float sinTheta = sqrtf(1 - cosTheta * cosTheta);
			const float phi = 2 * PI * PCG32::rand();

			// build orthonormal basis
			const float sign = copysignf(1, hitToLight.z);
			const float a = -1 / (sign + hitToLight.z);
			const float b = hitToLight.x * hitToLight.y * a;
			const float3 b1 = float3(1 + sign * hitToLight.x * hitToLight.x * a, sign * b, -sign * hitToLight.x);
			const float3 b2 = float3(b, sign + hitToLight.y * hitToLight.y * a, -hitToLight.y);

			// centre about direction to centre of light
			const float3 wi = cosTheta * hitToLight + sinTheta * std::cos(phi) * b1 + sinTheta * std::sin(phi) * b2;

			// calculate probability
			probNEE = probLight / (2 * PI * (1 - cosThetaMax));

			// check visibility
			HitInfo shadowHitInfo;
			if (globalScene.intersect(shadowHitInfo, Ray(hitPoint, wi)) && shadowHitInfo.material->type == LIGHT) {

				// outside
				float transmittance(1.0f), probDistance(1.0f);
				if (!hitInfo.inside) {
					transmittance = std::exp(-ATTENUATION * shadowHitInfo.t);
					probDistance = transmittance;
				}

				// multiple importance sampling
				probBRDF = hitInfo.material->pdf(wo, hitInfo.G, wi) * probDistance;
				const float weight = 1 / (probNEE + probBRDF);

				// accumulate radiance
				if (probBRDF > 0) radiance += weight * throughput * transmittance * hitInfo.material->brdf(wo, hitInfo.G, wi) * light->material.emission * dot(hitInfo.G, wi);
			}

			// continue path and update throughput
			ray = Ray(hitPoint, hitInfo.material->sample(wo, hitInfo.G));
			probBRDF = hitInfo.material->pdf(wo, hitInfo.G, ray.d);
			if (probBRDF > 0) throughput *= hitInfo.material->brdf(wo, hitInfo.G, ray.d) * dot(hitInfo.G, ray.d) / probBRDF;
			specular = false;
		}





		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		if (std::isnan(radiance.x) || std::isnan(radiance.y) || std::isnan(radiance.z)) std::cout << "radiance" << std::endl;
		if (std::isnan(throughput.x) || std::isnan(throughput.y) || std::isnan(throughput.z)) std::cout << "throughput" << std::endl;
		if (std::isnan(probBRDF)) std::cout << "probBRDF" << std::endl;
		if (std::isnan(probLight)) std::cout << "probLight" << std::endl;
		if (std::isnan(probNEE)) std::cout << "probNEE" << std::endl;
		if (std::isnan(ray.d.x) || std::isnan(ray.d.y) || std::isnan(ray.d.z)) std::cout << "ray.d" << std::endl;
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





		// russian roulette
		if (pathLength > MINIMUM_PATH_LENGTH) {
			float probabilityOfContinuing = std::max(throughput.x, std::max(throughput.y, throughput.z));
			if (probabilityOfContinuing < 1) {
				probabilityOfContinuing = std::max(0.25f, probabilityOfContinuing);
				if (PCG32::rand() < probabilityOfContinuing) throughput /= probabilityOfContinuing;
				else break;
			}
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