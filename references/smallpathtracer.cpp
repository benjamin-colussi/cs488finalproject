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
	Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25), DIFFUSE),
	Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75), DIFFUSE),
	Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75), DIFFUSE),
	Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(), DIFFUSE),
	Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75), DIFFUSE),
	Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75), DIFFUSE),

	// balls
	Sphere(16.5,Vec(27,16.5,47),Vec(),Vec(1,1,1)*.999, SPECULAR),
	Sphere(16.5,Vec(73,16.5,78),Vec(),Vec(1,1,1)*.999, SPECULAR),

	// light
	Sphere(1.5, Vec(50,81.6-16.5,81.6),Vec(4,4,4)*100,Vec(), DIFFUSE),
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
	double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; // takes max value
	if (++depth > 5 || !p) {
		if (erand48(Xi) < p) f = f * (1 / p);
		else return obj.emission * E;
	}



	// ideal diffuse reflection
	if (obj.refl == DIFFUSE) {

		// importance sampling with cosine, orthonormal basis transformation
		double r1 = 2 * M_PI * erand48(Xi); // angle around
		double r2 = erand48(Xi), r2s = sqrt(r2); // distance from centre
		Vec w = nl; // w = normal
		Vec u = ((fabs(w.x) > 0.1 ? Vec(0,1) : Vec(1)) % w).norm(); // u is perp to w
		Vec v = w % u; // v is perp to w and u
		Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm(); // random reflected ray

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

		// recursive call of light transport equation
		return obj.emission * E + e + f.mult(radiance(Ray(x, d), depth, Xi, 0));
	}



	// ideal specular reflection
	if (obj.refl == SPECULAR) {
		return obj.emission + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));
	}



	// ideal dielectric refraction
	Ray reflRay(x, r.d-n*2*n.dot(r.d));
	bool into = n.dot(nl)>0;                // Ray from outside going in?
	double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;

	// total internal reflection
	if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0) return obj.emission + f.mult(radiance(reflRay,depth,Xi));

	// no total internal reflection
	Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
	double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
	double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
	return obj.emission + f.mult(depth>2 ? (erand48(Xi)<P ?   // Russian roulette
	radiance(reflRay,depth,Xi)*RP:radiance(Ray(x,tdir),depth,Xi)*TP) :
	radiance(reflRay,depth,Xi)*Re+radiance(Ray(x,tdir),depth,Xi)*Tr);
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
