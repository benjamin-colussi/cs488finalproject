// smallpt, a Path Tracer by Kevin Beason, 2009
// Make : g++ -O3 -fopenmp explicit.cpp -o explicit
//        Remove "-fopenmp" for g++ version < 4.2
// Usage: time ./explicit 16 && xv image.ppm

#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>

struct Vec {
	double x, y, z;
	Vec(double x_ = 0, double y_ = 0, double z_ = 0) { x = x_; y = y_; z = z_; }
	Vec operator+(const Vec& b) const { return Vec(x + b.x, y + b.y, z + b.z); }
	Vec operator-(const Vec& b) const { return Vec(x - b.x, y - b.y, z - b.z); }
	Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
	Vec mult(const Vec& b) const { return Vec(x * b.x, y * b.y, z * b.z); }
	Vec& norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
	double dot(const Vec& b) const { return x * b.x + y * b.y + z * b.z; }
	Vec operator%(Vec& b) {return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
};

struct Ray {
	Vec o, d;
	Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};

enum Refl_t {
	DIFFUSE, SPECULAR, REFRACTIVE
};

struct Sphere {

	double radius;
	Vec position, emission, colour;
	Refl_t refl;

	// ctor
	Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_):
		radius(rad_), position(p_), emission(e_), colour(c_), refl(refl_) {}

	// returns distance, 0 if no hit
	double intersect(const Ray &r) const {

		// solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		Vec op = position - r.o;
		double t;
		double eps = 1e-4;
		double b = op.dot(r.d);
		double det = b * b - op.dot(op) + radius * radius;
		if (det < 0) return 0;
		det = sqrt(det);
		return (t = b - det) > eps ? t : (t = b + det) > eps ? t : 0;
	}

};

// scene: radius, position, emission, colour, material
Sphere spheres[] = {

	// left, right, back, front, bottom, top
	Sphere(1e5, Vec( 1e5+1,40.8,81.6),  Vec(), Vec(.75,.25,.25), DIFFUSE),
	Sphere(1e5, Vec(-1e5+99,40.8,81.6), Vec(), Vec(.25,.25,.75), DIFFUSE),
	Sphere(1e5, Vec(50,40.8, 1e5),      Vec(), Vec(.75,.75,.75), DIFFUSE),
	Sphere(1e5, Vec(50,40.8,-1e5+170),  Vec(), Vec(),            DIFFUSE),
	Sphere(1e5, Vec(50, 1e5, 81.6),     Vec(), Vec(.75,.75,.75), DIFFUSE),
	Sphere(1e5, Vec(50,-1e5+81.6,81.6), Vec(), Vec(.75,.75,.75), DIFFUSE),

	// balls
	Sphere(16.5, Vec(27,16.5,47), Vec(), Vec(1,1,1) * .999, DIFFUSE),
	Sphere(16.5, Vec(73,16.5,78), Vec(), Vec(1,1,1) * .999, DIFFUSE),

	// light
	Sphere(1.5, Vec(50,81.6-16.5,81.6), Vec(4,4,4) * 100, Vec(), DIFFUSE),
};



// global variables
int numSpheres = sizeof(spheres) / sizeof(Sphere);
inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1 : x; }
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }
inline bool intersect(const Ray& r, double& t, int& id) {
	double n = sizeof(spheres) / sizeof(Sphere);
	double d;
	double inf = t = 1e20;
	for(int i = int(n); i--; ) {
		if((d = spheres[i].intersect(r)) && d < t) {
			t = d;
			id = i;
		}
	}
	return t < inf;
}



// shade
Vec radiance(const Ray& r, int depth, unsigned short* Xi, int E = 1) {



	// hit info
	double t; // distance to intersection
	int id = 0; // id of intersected object
	if (!intersect(r, t, id)) return Vec(); // if miss return black
	const Sphere& obj = spheres[id]; // the hit object

	// surface properties
	Vec x = r.o + r.d * t; // ray intersection point
	Vec n = (x - obj.position).norm(); // sphere normal
	Vec nl = n.dot(r.d) < 0 ? n : n * -1; // properly oriented surface normal
	Vec f = obj.colour; // object colour, BRDF modulator

	// max reflection using russian roulette
	double prob = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;
	if (++depth > 5 || !prob) {
		if (erand48(Xi) < prob) f = f * (1 / prob);
		else return obj.emission * E;
	}
	// if (++depth > 5) return obj.emission * E;



	// ideal diffuse reflection
	if (obj.refl == DIFFUSE) {

		// loop over explicit lights
		Vec e;
		for (int i = 0; i < numSpheres; i++) {

			// skip non-lights
			const Sphere& s = spheres[i];
			if (s.emission.x <= 0 && s.emission.y <= 0 && s.emission.z <= 0) continue;

			// generate random direction towards light spheres, sampling a point on each light source
			Vec sw = s.position - x;
			Vec su = ((fabs(sw.x) > 0.1 ? Vec(0, 1) : Vec(1)) % sw).norm();
			Vec sv = sw % su;
			double cos_a_max = sqrt(1 - s.radius * s.radius / (x - s.position).dot(x - s.position));
			double eps1 = erand48(Xi), eps2 = erand48(Xi);
			double cos_a = 1 - eps1 + eps1 * cos_a_max;
			double sin_a = sqrt(1 - cos_a * cos_a);
			double phi = 2 * M_PI * eps2;
			Vec l = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;
			l.norm();

			// trace shadow ray, check if the connection is not blocked
			if (intersect(Ray(x, l), t, id) && id == i) {
				double omega = 2 * M_PI * (1 - cos_a_max);
				e = e + f.mult(s.emission * l.dot(nl) * omega) * M_1_PI;  // 1 / pi for brdf
			}
		}

		// importance sampling with cosine, orthonormal basis transformation
		double r1 = 2 * M_PI * erand48(Xi); // angle around
		double r2 = erand48(Xi), r2s = sqrt(r2); // distance from centre
		Vec w = nl; // w = normal
		Vec u = ((fabs(w.x) > 0.1 ? Vec(0,1) : Vec(1)) % w).norm(); // u is perp to w
		Vec v = w % u; // v is perp to w and u
		Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm(); // random reflected ray

		// recursive call of light transport equation
		return obj.emission * E + e + f.mult(radiance(Ray(x, d), depth, Xi, 0));
	}



	// ideal specular reflection
	if (obj.refl == SPECULAR) {
		return obj.emission + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));
	}





	// ideal dielectric refraction
	Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // reflected ray
	bool into = n.dot(nl) > 0; // ray from outside going in?
	double nc = 1; // refraction index of air
	double nt = 1.5; // refraction index of glass
	double nnt = into ? nc / nt : nt / nc;
	double ddn = r.d.dot(nl); // dot(ray.d, geometric normal)
	double cos2t;

	// total internal reflection
	if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) {
		return obj.emission + f.mult(radiance(reflRay, depth, Xi));
	}

	// no total internal reflection
	Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();

	double a = nt - nc;
	double b = nt + nc;
	double R0 = a * a / (b * b);

	double c = 1 - (into ? -ddn : tdir.dot(n));
	double Re = R0 + (1 - R0) * c * c * c * c * c; // R
	double Tr = 1 - Re; // T

	double P = 0.25 + 0.5 * Re; // clamping prob
	double RP = Re / P; // R / P
	double TP = Tr / (1 - P); // T / (1 - P)



	

	// russian roulette
	return obj.emission + f.mult(depth > 2 ? (erand48(Xi) < P ?
		radiance(reflRay, depth, Xi) * RP :
		radiance(Ray(x, tdir), depth, Xi) * TP) :
		radiance(reflRay, depth, Xi) * Re + radiance(Ray(x, tdir), depth, Xi) * Tr);
}



// path tracer
int main(int argc, char *argv[]) {

    // image size
    int w = 1024, h = 768;

    // number of samples
    int samps = argc == 2 ? atoi(argv[1]) / 4 : 1;

    // camera position and direction
    Ray cam(Vec(50,52,295.6), Vec(0,-0.042612,-1).norm());
    Vec cx = Vec(w * .5135 / h);
    Vec cy = (cx % cam.d).norm() * .5135;

	// accumulated radiance
    Vec r;

	// final image
    Vec* c = new Vec[w * h];

    // loop over rows
    for (int y = 0; y < h; y++) {

        // output progress
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));

        // loop over columns
		const unsigned short y3 = y * y * y;
    	unsigned short Xi[3] = {0, 0, y3};
        for (unsigned short x = 0; x < w; x++) {

            // 2x2 subpixel rows
            for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++) {

                // 2x2 subpixel columns
                for (int sx = 0; sx < 2; sx++, r = Vec()) {

                    // samples
                    for (int s = 0; s < samps; s++) {

                        double r1 = 2 * erand48(Xi);
                        double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);

                        double r2 = 2 * erand48(Xi);
                        double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);

                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) + cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;

                        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi) * (1.0 / samps);

                    } // camera rays are pushed ^^^^^ forward to start in interior

                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * 0.25;
                }
            }
        }
    }

    // write to .ppm file
    // FILE *f = fopen("image.ppm", "w");
	const std::string s = "image_" + std::to_string(samps * 4) + "_spp.ppm";
	FILE *f = fopen(s.c_str(), "w");
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++) {
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
    }

}




















// SMALL PT METHOD FOR GETTING RANDOM DIRECTION TO SPHERE LIGHT

// const float cosThetaMax = sqrtf(1 - light->radius * light->radius * oneOverDistanceSquared);
// const float cosTheta = 1 - Bertrand + Randolf * cosThetaMax;
// const float sinTheta = sqrtf(1 - cosTheta * cosTheta);
// const float phi = 2 * PI * Randolf;

// // build orthonormal basis with normal
// float sign = copysignf(1, lightNormal.z);
// const float a = -1 / (sign + lightNormal.z);
// const float b = lightNormal.x * lightNormal.y * a;
// const float3 b1 = float3(1 + sign * lightNormal.x * lightNormal.x * a, sign * b, -sign * lightNormal.x);
// const float3 b2 = float3(b, sign + lightNormal.y * lightNormal.y * a, -lightNormal.y);

// // centre about normal
// const float3 d = normalize(cosTheta * lightNormal + sinTheta * std::cos(phi) * b1 + sinTheta * std::sin(phi) * b2);










// SAMPLE COSINE WEIGHTED HEMISPHERE

// // using malley's method
// else if (SAMPLE_COS_WEIGHTED_HEMISPHERE) {

// 	// sample cos weighted positive unit hemisphere
// 	// float x = 2 * PCG32::rand() - 1;
// 	// float y = 2 * PCG32::rand() - 1;
// 	// if (x == 0 && y == 0) return float3(x, y, 1);
// 	// float theta, r;
// 	// if (abs(x) > abs(y)) {
// 	// 	r = x;
// 	// 	theta = PI_OVER_FOUR * (y / x);
// 	// }
// 	// else {
// 	// 	r = y;
// 	// 	theta = PI_OVER_TWO - PI_OVER_FOUR * (x / y);
// 	// }
// 	// float3 d = float3(r * cos(theta), r * sin(theta), sqrtf(1 - x * x - y * y));

// 	// next attempt
// 	const float r = sqrtf(PCG32::rand());
// 	const float theta = 2 * PI * PCG32::rand();
// 	const float x = r * std::cos(theta);
// 	const float y = r * std::cos(theta);
// 	float3 d = float3(x, y, sqrtf(1 - x * x - y * y));

// 	// build orthonormal basis with normal
// 	float sign = copysignf(1, n.z);
// 	const float a = -1 / (sign + n.z);
// 	const float b = n.x * n.y * a;
// 	const float3 b1 = float3(1 + sign * n.x * n.x * a, sign * b, -sign * n.x);
// 	const float3 b2 = float3(b, sign + n.y * n.y * a, -n.y);

// 	// centre about normal
// 	return d.x * b1 + d.y * b2 + d.z * n;
// }








// FIRST ATTEMP AT NEE

// // calculate vector from hit to light
// const float3 lightNormal = normalize(hitPoint - light->centre);
// float3 lightPoint = light->sampleSurface(lightNormal);
// float3 hitToLight = lightPoint - hitPoint;
// const float oneOverDistanceSquared = 1 / dot(hitToLight, hitToLight);
// const float oneOverDistance = sqrtf(oneOverDistanceSquared);
// hitToLight *= oneOverDistance;

// // trace shadow ray from hit to light
// // maybe also confirm we are hitting the same light source that we are sampling ...
// // under construction ...
// HitInfo shadowHitInfo;
// if (globalScene.intersect(shadowHitInfo, Ray(hitPoint, hitToLight)) && shadowHitInfo.material->type == LIGHT) {

// 	// check if intersection with the sampled point
// 	if ((shadowHitInfo.t - EPSILON) * oneOverDistance < 1 && 1 < (shadowHitInfo.t + EPSILON) * oneOverDistance) {
// 		const float geometry = dot(hitInfo.G, hitToLight) * dot(shadowHitInfo.G, -hitToLight) * oneOverDistanceSquared;
// 		const float probLight = light->pdf(shadowHitInfo.G, -hitToLight);
// 		// const float weight = probLight / (probLight + hitInfo.material->pdf(hitInfo.G, hitToLight));
// 		// radiance += throughput * hitInfo.material->spectrum() * light->material.emission * geometry / probLight;
// 		radiance += 0.5f * throughput * hitInfo.material->spectrum() * light->material.emission * geometry / probLight;
// 	}
// }












// BUILDING CONIC SAMPLE TOWARDS SURFACE OF SPHERE LIGHT

// // distance from centre of spherical light
// // float3 lightNormal = hitPoint - light->centre;
// float3 lightNormal = light->centre - hitPoint;
// float oneOverDistanceSquared = 1 / dot(lightNormal, lightNormal);
// lightNormal *= sqrtf(oneOverDistanceSquared);

// // random direction in cone of spherical light
// // const float sinThetaMax2 = light->radius * light->radius * oneOverDistanceSquared;
// // const float sinThetaMax = sqrtf(sinThetaMax2);
// // const float cosThetaMax = sqrtf(std::max(0.0f, 1 - sinThetaMax2));
// // const float cosTheta = 1 + (cosThetaMax - 1) * PCG32::rand();
// // const float sinTheta2 = 1 - cosTheta * cosTheta;
// // const float cosAlpha = sinTheta2 / sinThetaMax + cosTheta * sqrtf(1 - sinTheta2 / sinThetaMax2);
// // const float sinAlpha = sqrtf(1 - cosAlpha * cosAlpha);
// // const float phi = 2 * PI * PCG32::rand();

// const float Bertrand = PCG32::rand();
// const float Randolf = PCG32::rand();
// const float cosThetaMax = sqrtf(1 - light->radius * light->radius * oneOverDistanceSquared);
// const float cosTheta = 1 - Bertrand + Randolf * cosThetaMax;
// const float sinTheta = sqrtf(1 - cosTheta * cosTheta);
// const float phi = 2 * PI * Randolf;

// // build orthonormal basis
// float sign = copysignf(1, lightNormal.z);
// const float a = -1 / (sign + lightNormal.z);
// const float b = lightNormal.x * lightNormal.y * a;
// const float3 b1 = float3(1 + sign * lightNormal.x * lightNormal.x * a, sign * b, -sign * lightNormal.x);
// const float3 b2 = float3(b, sign + lightNormal.y * lightNormal.y * a, -lightNormal.y);

// // centre about normal
// // const float3 lightPointNormal = cosAlpha * lightNormal + sinAlpha * std::cos(phi) * b1 + sinAlpha * std::sin(phi) * b2;
// // const float3 lightPoint = light->centre + light->radius * lightPointNormal;
// // float3 wi = lightPoint - hitPoint;
// // oneOverDistanceSquared = 1 / dot(wi, wi);
// // wi *= sqrtf(oneOverDistanceSquared);

// // const float tmax = length(d); // use this to clamp intersection routine
// const float3 wi = cosTheta * lightNormal + sinTheta * std::cos(phi) * b1 + sinTheta * std::sin(phi) * b2;

// // trace shadow ray from hit to light
// HitInfo shadowHitInfo;
// if (globalScene.intersect(shadowHitInfo, Ray(hitPoint, wi)) && shadowHitInfo.material->type == LIGHT) {
// 	const float pdf = 2 * PI * (1 - cosThetaMax);
// 	// const float geometry = dot(hitInfo.G, wi) * dot(lightPointNormal, -wi) * oneOverDistanceSquared;
// 	const float geometry = dot(hitInfo.G, wi) * dot(shadowHitInfo.G, -wi) / (shadowHitInfo.t * shadowHitInfo.t);
// 	radiance += 0.5f * throughput * hitInfo.material->spectrum() * light->material.emission * geometry * pdf;
// }








// SAMPLE UNIT DISK - BUT NOT CONCENTRICALLY

// next attempt
// const float r = sqrtf(PCG32::rand());
// const float theta = 2 * PI * PCG32::rand();
// const float x = r * std::cos(theta);
// const float y = r * std::cos(theta);
// float3 d = float3(x, y, sqrtf(1 - x * x - y * y));










/*

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

*/











/*

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

*/