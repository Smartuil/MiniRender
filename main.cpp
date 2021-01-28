#pragma once


#include "model.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
const TGAColor blue = TGAColor(0, 0, 255, 255);
const TGAColor yellow = TGAColor(255, 255, 0, 255);

Model* model = NULL;
const int width = 800;
const int height = 800;
const int depth = 255;
int* zbuffer = NULL;
Vec3f light_dir(0, 0, -1);
Vec3f camera(0, 0, 3);

//Lesson 1
void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
	for (float t = 0.; t < 1.; t += .1) {
		int x = x0 + (x1 - x0) * t;
		int y = y0 + (y1 - y0) * t;
		image.set(x, y, color);
	}
}

void line2(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
	for (int x = x0; x <= x1; x++) {
		float t = (x - x0) / (float)(x1 - x0);
		int y = y0 * (1. - t) + y1 * t;
		image.set(x, y, color);
	}
}

void line3(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
	bool steep = false;
	if (std::abs(x0 - x1) < std::abs(y0 - y1)) { // if the line is steep, we transpose the image 
		std::swap(x0, y0);
		std::swap(x1, y1);
		steep = true;
	}
	if (x0 > x1) { // make it left−to−right 
		std::swap(x0, x1);
		std::swap(y0, y1);
	}
	for (int x = x0; x <= x1; x++) {
		float t = (x - x0) / (float)(x1 - x0);
		int y = y0 * (1. - t) + y1 * t;
		if (steep) {
			image.set(y, x, color); // if transposed, de−transpose 
		}
		else {
			image.set(x, y, color);
		}
	}
}

void line4(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
	bool steep = false;
	if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
		std::swap(x0, y0);
		std::swap(x1, y1);
		steep = true;
	}
	if (x0 > x1) {
		std::swap(x0, x1);
		std::swap(y0, y1);
	}
	int dx = x1 - x0;
	int dy = y1 - y0;
	float derror = std::abs(dy / float(dx));
	float error = 0;
	int y = y0;
	for (int x = x0; x <= x1; x++) {
		if (steep) {
			image.set(y, x, color);
		}
		else {
			image.set(x, y, color);
		}
		error += derror;
		if (error > .5) {
			y += (y1 > y0 ? 1 : -1);
			error -= 1.;
		}
	}
}

void line5(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
	bool steep = false;
	if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
		std::swap(x0, y0);
		std::swap(x1, y1);
		steep = true;
	}
	if (x0 > x1) {
		std::swap(x0, x1);
		std::swap(y0, y1);
	}
	int dx = x1 - x0;
	int dy = y1 - y0;
	int derror2 = std::abs(dy) * 2;
	int error2 = 0;
	int y = y0;
	for (int x = x0; x <= x1; x++) {
		if (steep) {
			image.set(y, x, color);
		}
		else {
			image.set(x, y, color);
		}
		error2 += derror2;
		if (error2 > dx) {
			y += (y1 > y0 ? 1 : -1);
			error2 -= dx * 2;
		}
	}
}

//Lesson 2
void line6(Vec2i p0, Vec2i p1, TGAImage& image, TGAColor color) {
	bool steep = false;
	if (std::abs(p0.x - p1.x) < std::abs(p0.y - p1.y)) {
		std::swap(p0.x, p0.y);
		std::swap(p1.x, p1.y);
		steep = true;
	}
	if (p0.x > p1.x) {
		std::swap(p0, p1);
	}

	for (int x = p0.x; x <= p1.x; x++) {
		float t = (x - p0.x) / (float)(p1.x - p0.x);
		int y = p0.y * (1. - t) + p1.y * t;
		if (steep) {
			image.set(y, x, color);
		}
		else {
			image.set(x, y, color);
		}
	}
}

void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color) {
	line6(t0, t1, image, color);
	line6(t1, t2, image, color);
	line6(t2, t0, image, color);
}

void triangle2(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color) {
	// sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!) 
	if (t0.y > t1.y) std::swap(t0, t1);
	if (t0.y > t2.y) std::swap(t0, t2);
	if (t1.y > t2.y) std::swap(t1, t2);
	line6(t0, t1, image, green);
	line6(t1, t2, image, green);
	line6(t2, t0, image, red);
}

void triangle3(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color) {
	// sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!) 
	if (t0.y > t1.y) std::swap(t0, t1);
	if (t0.y > t2.y) std::swap(t0, t2);
	if (t1.y > t2.y) std::swap(t1, t2);
	int total_height = t2.y - t0.y;
	for (int y = t0.y; y <= t1.y; y++) {
		int segment_height = t1.y - t0.y + 1;
		float alpha = (float)(y - t0.y) / total_height;
		float beta = (float)(y - t0.y) / segment_height; // be careful with divisions by zero 
		Vec2i A = t0 + (t2 - t0) * alpha;
		Vec2i B = t0 + (t1 - t0) * beta;
		image.set(A.x, y, red);
		image.set(B.x, y, green);
	}
}

void triangle4(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color) {
	// sort the vertices, t0, t1, t2 lower−to−upper (bubblesort yay!) 
	if (t0.y > t1.y) std::swap(t0, t1);
	if (t0.y > t2.y) std::swap(t0, t2);
	if (t1.y > t2.y) std::swap(t1, t2);
	int total_height = t2.y - t0.y;
	for (int y = t0.y; y <= t1.y; y++) {
		int segment_height = t1.y - t0.y + 1;
		float alpha = (float)(y - t0.y) / total_height;
		float beta = (float)(y - t0.y) / segment_height; // be careful with divisions by zero 
		Vec2i A = t0 + (t2 - t0) * alpha;
		Vec2i B = t0 + (t1 - t0) * beta;
		if (A.x > B.x) std::swap(A, B);
		for (int j = A.x; j <= B.x; j++) {
			image.set(j, y, color); // attention, due to int casts t0.y+i != A.y 
		}
	}
	for (int y = t1.y; y <= t2.y; y++) {
		int segment_height = t2.y - t1.y + 1;
		float alpha = (float)(y - t0.y) / total_height;
		float beta = (float)(y - t1.y) / segment_height; // be careful with divisions by zero 
		Vec2i A = t0 + (t2 - t0) * alpha;
		Vec2i B = t1 + (t2 - t1) * beta;
		if (A.x > B.x) std::swap(A, B);
		for (int j = A.x; j <= B.x; j++) {
			image.set(j, y, color); // attention, due to int casts t0.y+i != A.y 
		}
	}
}

void triangle5(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage& image, TGAColor color) {
	if (t0.y == t1.y && t0.y == t2.y) return; // i dont care about degenerate triangles
	if (t0.y > t1.y) std::swap(t0, t1);
	if (t0.y > t2.y) std::swap(t0, t2);
	if (t1.y > t2.y) std::swap(t1, t2);
	int total_height = t2.y - t0.y;
	for (int i = 0; i < total_height; i++) {
		bool second_half = i > t1.y - t0.y || t1.y == t0.y;
		int seg_height = second_half ? t2.y - t1.y : t1.y - t0.y;
		float alpha = (float)i / total_height;
		float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / seg_height; // be careful: with above conditions no division by zero here
		Vec2i A = t0 + (t2 - t0) * alpha;
		Vec2i B = second_half ? t1 + (t2 - t1) * beta : t0 + (t1 - t0) * beta;
		if (A.x > B.x) std::swap(A, B);
		for (int j = A.x; j <= B.x; j++) {
			image.set(j, t0.y + i, color); // attention, due to int casts t0.y+i != A.y
		}
	}
}

Vec3f cross(Vec3f v1, Vec3f v2) {
	return Vec3f(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

Vec3f barycentric(Vec2i* pts, Vec2i P) {
	Vec3f u = cross(Vec3f(pts[2][0] - pts[0][0], pts[1][0] - pts[0][0], pts[0][0] - P[0]), Vec3f(pts[2][1] - pts[0][1], pts[1][1] - pts[0][1], pts[0][1] - P[1]));
	/* `pts` and `P` has integer value as coordinates
	   so `abs(u[2])` < 1 means `u[2]` is 0, that means
	   triangle is degenerate, in this case return something with negative coordinates */
	if (std::abs(u[2]) < 1) return Vec3f(-1, 1, 1);
	return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
}

void triangle6(Vec2i* pts, TGAImage& image, TGAColor color) {
	Vec2i bboxmin(image.get_width() - 1, image.get_height() - 1);
	Vec2i bboxmax(0, 0);
	Vec2i clamp(image.get_width() - 1, image.get_height() - 1);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			bboxmin[j] = std::max(0, std::min(bboxmin[j], pts[i][j]));
			bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
		}
	}
	Vec2i P;
	for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
		for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
			Vec3f bc_screen = barycentric(pts, P);
			if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
			image.set(P.x, P.y, color);
		}
	}
}

//Lesson3
void line7(Vec2i p0, Vec2i p1, TGAImage& image, TGAColor color) {
	bool steep = false;
	if (std::abs(p0.x - p1.x) < std::abs(p0.y - p1.y)) {
		std::swap(p0.x, p0.y);
		std::swap(p1.x, p1.y);
		steep = true;
	}
	if (p0.x > p1.x) {
		std::swap(p0, p1);
	}

	for (int x = p0.x; x <= p1.x; x++) {
		float t = (x - p0.x) / (float)(p1.x - p0.x);
		int y = p0.y * (1. - t) + p1.y * t + .5;
		if (steep) {
			image.set(y, x, color);
		}
		else {
			image.set(x, y, color);
		}
	}
}

void rasterize(Vec2i p0, Vec2i p1, TGAImage& image, TGAColor color, int ybuffer[]) {
	if (p0.x > p1.x) {
		std::swap(p0, p1);
	}
	for (int x = p0.x; x <= p1.x; x++) {
		float t = (x - p0.x) / (float)(p1.x - p0.x);
		int y = p0.y * (1. - t) + p1.y * t + .5;
		if (ybuffer[x] < y) {
			ybuffer[x] = y;
			image.set(x, 0, color);
		}
	}
}

Vec3f barycentric2(Vec3f A, Vec3f B, Vec3f C, Vec3f P) {
	Vec3f s[2];
	//for (int i = 2; i--; ) {
	//	std::cout << C[i] << " " << A[i];
	//	s[i][0] = C[i] - A[i];
	//	s[i][1] = B[i] - A[i];
	//	s[i][2] = A[i] - P[i];
	//}
	s[0] = Vec3f(C[0] - A[0], B[0] - A[0], A[0] - P[0]);
	s[1] = Vec3f(C[1] - A[1], B[1] - A[1], A[1] - P[1]);
	Vec3f u = cross(s[0], s[1]);
	if (std::abs(u[2]) > 1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
		return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
	return Vec3f(-1, 1, 1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

void triangle7(Vec3f* pts, float* zbuffer, TGAImage& image, TGAColor color) {
	//包围盒
	Vec2f bboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
	Vec2f clamp(image.get_width() - 1, image.get_height() - 1);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			bboxmin[j] = std::max(0.f, std::min(bboxmin[j], pts[i][j]));
			bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
		}
	}

	Vec3f P;
	for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
		for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
			Vec3f bc_screen = barycentric2(pts[0], pts[1], pts[2], P);
			if (bc_screen.x < 0 || bc_screen.y < 0 || bc_screen.z < 0) continue;
			P.z = 0;
			for (int i = 0; i < 3; i++) P.z += pts[i][2] * bc_screen[i];//？？？
			if (zbuffer[int(P.x + P.y * width)] < P.z) {
				zbuffer[int(P.x + P.y * width)] = P.z;
				image.set(P.x, P.y, color);
			}
		}
	}
}

Vec3f world2screen(Vec3f v) {
	return Vec3f(
		int((v.x + 1.) * width / 2. + .5),
		int((v.y + 1.) * height / 2. + .5),
		v.z);
}

void triangle8(Vec3i t0, Vec3i t1, Vec3i t2, TGAImage& image, TGAColor color, int* zbuffer) {
	if (t0.y == t1.y && t0.y == t2.y) return; // i dont care about degenerate triangles
	if (t0.y > t1.y) std::swap(t0, t1);
	if (t0.y > t2.y) std::swap(t0, t2);
	if (t1.y > t2.y) std::swap(t1, t2);

	int total_height = t2.y - t0.y;
	for (int i = 0; i < total_height; i++) {
		bool second_half = i > t1.y - t0.y || t1.y == t0.y;
		int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
		float alpha = (float)i / total_height;
		float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height; // be careful: with above conditions no division by zero here
		Vec3i A = t0 + (t2 - t0) * alpha;
		Vec3i B = second_half ? t1 + (t2 - t1) * beta : t0 + (t1 - t0) * beta;
		if (A.x > B.x) std::swap(A, B);
		for (int j = A.x; j <= B.x; j++) {
			float phi = B.x == A.x ? 1. : (float)(j - A.x) / (float)(B.x - A.x);
			Vec3i P = A + (B - A) * phi;
			P.x = j; P.y = t0.y + i; // a hack to fill holes (due to int cast precision problems)
			int idx = j + (t0.y + i) * width;
			if (zbuffer[idx] < P.z) {
				zbuffer[idx] = P.z;
				image.set(P.x, P.y, color); // attention, due to int casts t0.y+i != A.y
			}
		}
	}
}

void triangle9(Vec3i t0, Vec3i t1, Vec3i t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, TGAImage& image, float intensity, int* zbuffer) {
	if (t0.y == t1.y && t0.y == t2.y) return; // i dont care about degenerate triangles
	if (t0.y > t1.y) { std::swap(t0, t1); std::swap(uv0, uv1); }
	if (t0.y > t2.y) { std::swap(t0, t2); std::swap(uv0, uv2); }
	if (t1.y > t2.y) { std::swap(t1, t2); std::swap(uv1, uv2); }

	int total_height = t2.y - t0.y;
	for (int i = 0; i < total_height; i++) {
		bool second_half = i > t1.y - t0.y || t1.y == t0.y;
		int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
		float alpha = (float)i / total_height;
		float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height; // be careful: with above conditions no division by zero here
		Vec3i A = t0 + Vec3f(t2 - t0) * alpha;
		Vec3i B = second_half ? t1 + Vec3f(t2 - t1) * beta : t0 + Vec3f(t1 - t0) * beta;
		Vec2i uvA = uv0 + (uv2 - uv0) * alpha;
		Vec2i uvB = second_half ? uv1 + (uv2 - uv1) * beta : uv0 + (uv1 - uv0) * beta;
		if (A.x > B.x) { std::swap(A, B); std::swap(uvA, uvB); }
		for (int j = A.x; j <= B.x; j++) {
			float phi = B.x == A.x ? 1. : (float)(j - A.x) / (B.x - A.x);
			Vec3i   P = Vec3f(A) + Vec3f(B - A) * phi;
			Vec2i uvP = uvA + (uvB - uvA) * phi;
			int idx = P.x + P.y * width;
			if (zbuffer[idx] < P.z) {
				zbuffer[idx] = P.z;
				//TGAColor color = model->diffuse(uvP);
				//image.set(P.x, P.y, TGAColor(color.r * intensity, color.g * intensity, color.b * intensity));
				image.set(P.x, P.y, model->diffuse(uvP) * intensity);
			}
		}
	}
}

//Lesson4
void line8(Vec3i p0, Vec3i p1, TGAImage& image, TGAColor color) {
	bool steep = false;
	if (std::abs(p0.x - p1.x) < std::abs(p0.y - p1.y)) {
		std::swap(p0.x, p0.y);
		std::swap(p1.x, p1.y);
		steep = true;
	}
	if (p0.x > p1.x) {
		std::swap(p0, p1);
	}

	for (int x = p0.x; x <= p1.x; x++) {
		float t = (x - p0.x) / (float)(p1.x - p0.x);
		int y = p0.y * (1. - t) + p1.y * t + .5;
		if (steep) {
			image.set(y, x, color);
		}
		else {
			image.set(x, y, color);
		}
	}
}

Vec3f m2v(Matrix m) {
	return Vec3f(m[0][0] / m[3][0], m[1][0] / m[3][0], m[2][0] / m[3][0]);
}

Matrix v2m(Vec3f v) {
	Matrix m(4, 1);
	m[0][0] = v.x;
	m[1][0] = v.y;
	m[2][0] = v.z;
	m[3][0] = 1.f;
	return m;
}

Matrix viewport(int x, int y, int w, int h) {
	Matrix m = Matrix::identity(4);
	m[0][3] = x + w / 2.f;
	m[1][3] = y + h / 2.f;
	m[2][3] = depth / 2.f;

	m[0][0] = w / 2.f;
	m[1][1] = h / 2.f;
	m[2][2] = depth / 2.f;
	return m;
}

Matrix translation(Vec3f v) {
	Matrix Tr = Matrix::identity(4);
	Tr[0][3] = v.x;
	Tr[1][3] = v.y;
	Tr[2][3] = v.z;
	return Tr;
}

Matrix zoom(float factor) {
	Matrix Z = Matrix::identity(4);
	Z[0][0] = Z[1][1] = Z[2][2] = factor;
	return Z;
}

Matrix rotation_x(float cosangle, float sinangle) {
	Matrix R = Matrix::identity(4);
	R[1][1] = R[2][2] = cosangle;
	R[1][2] = -sinangle;
	R[2][1] = sinangle;
	return R;
}

Matrix rotation_y(float cosangle, float sinangle) {
	Matrix R = Matrix::identity(4);
	R[0][0] = R[2][2] = cosangle;
	R[0][2] = sinangle;
	R[2][0] = -sinangle;
	return R;
}

Matrix rotation_z(float cosangle, float sinangle) {
	Matrix R = Matrix::identity(4);
	R[0][0] = R[1][1] = cosangle;
	R[0][1] = -sinangle;
	R[1][0] = sinangle;
	return R;
}

void triangle10(Vec3i t0, Vec3i t1, Vec3i t2, Vec2i uv0, Vec2i uv1, Vec2i uv2, TGAImage& image, float intensity, int* zbuffer) {
	if (t0.y == t1.y && t0.y == t2.y) return; // i dont care about degenerate triangles
	if (t0.y > t1.y) { std::swap(t0, t1); std::swap(uv0, uv1); }
	if (t0.y > t2.y) { std::swap(t0, t2); std::swap(uv0, uv2); }
	if (t1.y > t2.y) { std::swap(t1, t2); std::swap(uv1, uv2); }

	int total_height = t2.y - t0.y;
	for (int i = 0; i < total_height; i++) {
		bool second_half = i > t1.y - t0.y || t1.y == t0.y;
		int segment_height = second_half ? t2.y - t1.y : t1.y - t0.y;
		float alpha = (float)i / total_height;
		float beta = (float)(i - (second_half ? t1.y - t0.y : 0)) / segment_height; // be careful: with above conditions no division by zero here
		Vec3i A = t0 + Vec3f(t2 - t0) * alpha;
		Vec3i B = second_half ? t1 + Vec3f(t2 - t1) * beta : t0 + Vec3f(t1 - t0) * beta;
		Vec2i uvA = uv0 + (uv2 - uv0) * alpha;
		Vec2i uvB = second_half ? uv1 + (uv2 - uv1) * beta : uv0 + (uv1 - uv0) * beta;
		if (A.x > B.x) { std::swap(A, B); std::swap(uvA, uvB); }
		for (int j = A.x; j <= B.x; j++) {
			float phi = B.x == A.x ? 1. : (float)(j - A.x) / (B.x - A.x);
			Vec3i   P = Vec3f(A) + Vec3f(B - A) * phi;
			Vec2i uvP = uvA + (uvB - uvA) * phi;
			int idx = P.x + P.y * width;
			if (zbuffer[idx] < P.z) {
				zbuffer[idx] = P.z;
				image.set(P.x, P.y, model->diffuse(uvP) * intensity);
			}
		}
	}
}

int main(int argc, char** argv) {

	/*{
		//Lesson 0
		TGAImage image(100, 100, TGAImage::RGB);
		image.set(52, 41, red);
		image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
		image.write_tga_file("output.tga");
		return 0;

		//Lesson 2 First attempt
		TGAImage image(100, 100, TGAImage::RGB);
		//line2(13, 20, 80, 40, image, white);
		line3(13, 20, 80, 40, image, white);
		line3(20, 13, 40, 80, image, red);
		line3(80, 40, 13, 20, image, red);
		image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
		image.write_tga_file("output.tga");
		return 0;

		if (2 == argc) {
			model = new Model(argv[1]);
		}
		else {
			model = new Model("obj/african_head/african_head.obj");
		}

		TGAImage image(width, height, TGAImage::RGB);
		for (int i = 0; i < model->nfaces(); i++) {
			std::vector<int> face = model->face(i);
			for (int j = 0; j < 3; j++) {
				Vec3f v0 = model->vert(face[j]);
				Vec3f v1 = model->vert(face[(j + 1) % 3]);
				int x0 = (v0.x + 1.) * width / 2.;
				int y0 = (v0.y + 1.) * height / 2.;
				int x1 = (v1.x + 1.) * width / 2.;
				int y1 = (v1.y + 1.) * height / 2.;
				line3(x0, y0, x1, y1, image, white);
			}
		}

		//image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
		image.flip_horizontally();
		image.write_tga_file("output.tga");
		delete model;
		return 0;
	}*/

	/*{
		//Lesson 2
		TGAImage image(width, height, TGAImage::RGB);

		Vec2i t0[3] = { Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80) };
		Vec2i t1[3] = { Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180) };
		Vec2i t2[3] = { Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180) };

		triangle5(t0[0], t0[1], t0[2], image, red);
		triangle5(t1[0], t1[1], t1[2], image, white);
		triangle5(t2[0], t2[1], t2[2], image, green);


		//image.flip_horizontally(); // i want to have the origin at the left bottom corner of the image
		//image.flip_vertically();
		image.write_tga_file("output.tga");
		return 0;

		TGAImage frame(200, 200, TGAImage::RGB);
		Vec2i pts[3] = { Vec2i(10,10), Vec2i(100, 30), Vec2i(190, 160) };
		triangle6(pts, frame, TGAColor(255, 0, 0));
		frame.write_tga_file("framebuffer.tga");
		return 0;

		if (2 == argc) {
			model = new Model(argv[1]);
		}
		else {
			model = new Model("obj/african_head/african_head.obj");
		}
		TGAImage image(width, height, TGAImage::RGB);
		Vec3f light_dir(0, 0, -1);
		//for (int i = 0; i < model->nfaces(); i++) {
		//	std::vector<int> face = model->face(i);
		//	Vec2i screen_coords[3];
		//	for (int j = 0; j < 3; j++) {
		//		Vec3f world_coords = model->vert(face[j]);
		//		screen_coords[j] = Vec2i((world_coords.x + 1.) * width / 2., (world_coords.y + 1.) * height / 2.);
		//	}
		//	triangle6(screen_coords, image, TGAColor(rand() % 255, rand() % 255, rand() % 255, 255));
		//}

		for (int i = 0; i < model->nfaces(); i++) {
			std::vector<int> face = model->face(i);
			Vec2i screen_coords[3];
			Vec3f world_coords[3];
			for (int j = 0; j < 3; j++) {
				Vec3f v = model->vert(face[j]);
				screen_coords[j] = Vec2i((v.x + 1.) * width / 2., (v.y + 1.) * height / 2.);
				world_coords[j] = v;
			}
			Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
			n.normalize();
			float intensity = n * light_dir;
			if (intensity > 0) {
				triangle6(screen_coords, image, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
			}
		}

		//image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
		image.write_tga_file("output.tga");
		delete model;
		return 0;
	}*/

	/*{
		//Lesson 3
		//{
		//	// just dumping the 2d scene (yay we have enough dimensions!)
		//	TGAImage scene(width, height, TGAImage::RGB);

		//	// scene "2d mesh"
		//	line7(Vec2i(20, 34), Vec2i(744, 400), scene, red);
		//	line7(Vec2i(120, 434), Vec2i(444, 400), scene, green);
		//	line7(Vec2i(330, 463), Vec2i(594, 200), scene, blue);

		//	// screen line
		//	line7(Vec2i(10, 10), Vec2i(790, 10), scene, white);

		//	//scene.flip_vertically(); // i want to have the origin at the left bottom corner of the image
		//	scene.write_tga_file("scene.tga");
		//}

		//{
		//	TGAImage render(width, 16, TGAImage::RGB);
		//	int ybuffer[width];
		//	for (int i = 0; i < width; i++) {
		//		ybuffer[i] = std::numeric_limits<int>::min();
		//	}
		//	rasterize(Vec2i(20, 34), Vec2i(744, 400), render, red, ybuffer);
		//	rasterize(Vec2i(120, 434), Vec2i(444, 400), render, green, ybuffer);
		//	rasterize(Vec2i(330, 463), Vec2i(594, 200), render, blue, ybuffer);

		//	// 1-pixel wide image is bad for eyes, lets widen it
		//	for (int i = 0; i < width; i++) {
		//		for (int j = 1; j < 16; j++) {
		//			render.set(i, j, render.get(i, 0));
		//		}
		//	}
		//	render.flip_vertically(); // i want to have the origin at the left bottom corner of the image
		//	render.write_tga_file("render.tga");
		//}

		{
			//和下边的方法结果一样
			if (2 == argc) {
				model = new Model(argv[1]);
			}
			else {
				model = new Model("obj/african_head/african_head.obj");
			}
			float* zbuffer = new float[width * height];
			for (int i = width * height; i--; zbuffer[i] = -std::numeric_limits<float>::max());

			TGAImage image(width, height, TGAImage::RGB);
			Vec3f light_dir(0, 0, -1);

			for (int i = 0; i < model->nfaces(); i++) {
				std::vector<int> face = model->face(i);
				Vec2i screen_coords[3];
				Vec3f world_coords[3];
				for (int j = 0; j < 3; j++) {
					Vec3f v = model->vert(face[j]);
					screen_coords[j] = Vec2i((v.x + 1.) * width / 2., (v.y + 1.) * height / 2.);
					world_coords[j] = v;
				}
				Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
				n.normalize();
				float intensity = n * light_dir;
				if (intensity > 0) {
					Vec3f pts[3];
					for (int i = 0; i < 3; i++) pts[i] = Vec3f(screen_coords[i].x, screen_coords[i].y, world_coords[i].z);
					triangle7(pts, zbuffer, image, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255));
				}
			}

			image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
			image.write_tga_file("output.tga");//生成的图片稍微有点亮（和网站图片相比）
			delete model;
		}

		{
			if (2 == argc) {
				model = new Model(argv[1]);
			}
			else {
				model = new Model("obj/african_head/african_head.obj");
			}

			//初始化zbuffer
			zbuffer = new int[width * height];
			for (int i = 0; i < width * height; i++) {
				zbuffer[i] = std::numeric_limits<int>::min();
			}

			{
				// draw the model
				TGAImage image(width, height, TGAImage::RGB);

				for (int i = 0; i < model->nfaces(); i++) {
					std::vector<int> face = model->face(i);
					Vec3i screen_coords[3];
					Vec3f world_coords[3];
					for (int j = 0; j < 3; j++) {
						Vec3f v = model->vert(face[j]);
						screen_coords[j] = Vec3i((v.x + 1.) * width / 2., (v.y + 1.) * height / 2., (v.z + 1.) * depth / 2.);
						world_coords[j] = v;
					}
					Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
					n.normalize();
					float intensity = n * light_dir;
					if (intensity > 0) {
						triangle8(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(intensity * 255, intensity * 255, intensity * 255, 255), zbuffer);
					}
				}
				//image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
				image.write_tga_file("output.tga");
			}

			{
				// dump z-buffer (debugging purposes only)
				TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
				for (int i = 0; i < width; i++) {
					for (int j = 0; j < height; j++) {
						zbimage.set(i, j, TGAColor(zbuffer[i + j * width], 1));
					}
				}
				//zbimage.flip_vertically(); // i want to have the origin at the left bottom corner of the image
				zbimage.write_tga_file("zbuffer.tga");
			}
			delete model;
			delete[] zbuffer;
		}

	}

	{
		//作业
		if (2 == argc) {
			model = new Model(argv[1]);
		}
		else {
			model = new Model("obj/african_head/african_head.obj");
		}

		zbuffer = new int[width * height];
		for (int i = 0; i < width * height; i++) {
			zbuffer[i] = std::numeric_limits<int>::min();
		}

		{
			// draw the model
			TGAImage image(width, height, TGAImage::RGB);
			for (int i = 0; i < model->nfaces(); i++) {
				std::vector<int> face = model->face(i);
				Vec3i screen_coords[3];
				Vec3f world_coords[3];
				for (int j = 0; j < 3; j++) {
					Vec3f v = model->vert(face[j]);
					screen_coords[j] = Vec3i((v.x + 1.) * width / 2., (v.y + 1.) * height / 2., (v.z + 1.) * depth / 2.);
					world_coords[j] = v;
				}
				Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
				n.normalize();
				float intensity = n * light_dir;
				if (intensity > 0) {
					Vec2i uv[3];
					for (int k = 0; k < 3; k++) {
						uv[k] = model->uv(i, k);
					}
					triangle9(screen_coords[0], screen_coords[1], screen_coords[2], uv[0], uv[1], uv[2], image, intensity, zbuffer);
				}
			}

			image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
			image.write_tga_file("output.tga");
		}

		{
			// dump z-buffer (debugging purposes only)
			TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					zbimage.set(i, j, TGAColor(zbuffer[i + j * width], 1));
				}
			}
			zbimage.flip_vertically(); // i want to have the origin at the left bottom corner of the image
			zbimage.write_tga_file("zbuffer.tga");
		}
		delete model;
		delete[] zbuffer;
	}*/

	/*{
		//Lesson4
		{
			//if (2 == argc) {
			//	model = new Model(argv[1]);
			//}
			//else {
			//	model = new Model("obj/african_head/african_head.obj");
			//}

			//TGAImage image(width, height, TGAImage::RGB);
			//Matrix VP = viewport(width / 4, width / 4, width / 2, height / 2);

			//{ // draw the axes
			//	Vec3f x(1.f, 0.f, 0.f), y(0.f, 1.f, 0.f), o(0.f, 0.f, 0.f);
			//	o = m2v(VP * v2m(o));
			//	x = m2v(VP * v2m(x));
			//	y = m2v(VP * v2m(y));
			//	line8(o, x, image, red);
			//	line8(o, y, image, green);
			//}


			//for (int i = 0; i < model->nfaces(); i++) {
			//	std::vector<int> face = model->face(i);
			//	for (int j = 0; j < (int)face.size(); j++) {
			//		Vec3f wp0 = model->vert(face[j]);
			//		Vec3f wp1 = model->vert(face[(j + 1) % face.size()]);

			//		{ // draw the original model
			//			Vec3f sp0 = m2v(VP * v2m(wp0));
			//			Vec3f sp1 = m2v(VP * v2m(wp1));
			//			line8(sp0, sp1, image, white);
			//		}

			//		{ // draw the deformed model
			//			Matrix T = zoom(1.5);
			//			//                  Matrix T = Matrix::identity(4);
			//			//                  T[0][1] = 0.333;
			//			//                Matrix T = translation(Vec3f(.33, .5, 0))*rotation_z(cos(10.*M_PI/180.), sin(10.*M_PI/180.));
			//			Vec3f sp0 = m2v(VP * T * v2m(wp0));
			//			Vec3f sp1 = m2v(VP * T * v2m(wp1));
			//			line8(sp0, sp1, image, yellow);
			//		}
			//	}
			//	break;
			//}


			//image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
			//image.write_tga_file("output.tga");
			//delete model;
		}

		{
			if (2 == argc) {
				model = new Model(argv[1]);
			}
			else {
				model = new Model("obj/african_head/african_head.obj");
			}

			zbuffer = new int[width * height];
			for (int i = 0; i < width * height; i++) {
				zbuffer[i] = std::numeric_limits<int>::min();
			}

			{ // draw the model
				Matrix Projection = Matrix::identity(4);
				Matrix ViewPort = viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
				Projection[3][2] = -1.f / camera.z;

				TGAImage image(width, height, TGAImage::RGB);
				for (int i = 0; i < model->nfaces(); i++) {
					std::vector<int> face = model->face(i);
					Vec3i screen_coords[3];
					Vec3f world_coords[3];
					for (int j = 0; j < 3; j++) {
						Vec3f v = model->vert(face[j]);
						screen_coords[j] = Vec3f(ViewPort*Projection*Matrix(v));
						world_coords[j] = v;
					}
					Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
					n.normalize();
					float intensity = n * light_dir;
					if (intensity > 0) {
						Vec2i uv[3];
						for (int k = 0; k < 3; k++) {
							uv[k] = model->uv(i, k);
						}
						triangle10(screen_coords[0], screen_coords[1], screen_coords[2], uv[0], uv[1], uv[2], image, intensity, zbuffer);
					}
				}

				image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
				image.write_tga_file("output.tga");
			}

			{ // dump z-buffer (debugging purposes only)
				TGAImage zbimage(width, height, TGAImage::GRAYSCALE);
				for (int i = 0; i < width; i++) {
					for (int j = 0; j < height; j++) {
						zbimage.set(i, j, TGAColor(zbuffer[i + j * width], 1));
					}
				}
				zbimage.flip_vertically(); // i want to have the origin at the left bottom corner of the image
				zbimage.write_tga_file("zbuffer.tga");
			}
			delete model;
			delete[] zbuffer;
		}
	}*/
	return 0;
}