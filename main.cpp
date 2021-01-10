#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);
const TGAColor green = TGAColor(0, 255, 0, 255);
Model* model = NULL;
const int width = 800;
const int height = 800;

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

int main(int argc, char** argv) {
	//Lesson 0
	//TGAImage image(100, 100, TGAImage::RGB);
	//image.set(52, 41, red);
	//image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	//image.write_tga_file("output.tga");
	//return 0;

	//Lesson 2 First attempt
	//TGAImage image(100, 100, TGAImage::RGB);
	////line2(13, 20, 80, 40, image, white);
	//line3(13, 20, 80, 40, image, white);
	//line3(20, 13, 40, 80, image, red);
	//line3(80, 40, 13, 20, image, red);
	//image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	//image.write_tga_file("output.tga");
	//return 0;

	//if (2 == argc) {
	//	model = new Model(argv[1]);
	//}
	//else {
	//	model = new Model("obj/african_head/african_head.obj");
	//}

	//TGAImage image(width, height, TGAImage::RGB);
	//for (int i = 0; i < model->nfaces(); i++) {
	//	std::vector<int> face = model->face(i);
	//	for (int j = 0; j < 3; j++) {
	//		Vec3f v0 = model->vert(face[j]);
	//		Vec3f v1 = model->vert(face[(j + 1) % 3]);
	//		int x0 = (v0.x + 1.) * width / 2.;
	//		int y0 = (v0.y + 1.) * height / 2.;
	//		int x1 = (v1.x + 1.) * width / 2.;
	//		int y1 = (v1.y + 1.) * height / 2.;
	//		line3(x0, y0, x1, y1, image, white);
	//	}
	//}

	////image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	//image.flip_horizontally();
	//image.write_tga_file("output.tga");
	//delete model;
	//return 0;

	//Lesson 2
	//TGAImage image(width, height, TGAImage::RGB);

	//Vec2i t0[3] = { Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80) };
	//Vec2i t1[3] = { Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180) };
	//Vec2i t2[3] = { Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180) };

	//triangle5(t0[0], t0[1], t0[2], image, red);
	//triangle5(t1[0], t1[1], t1[2], image, white);
	//triangle5(t2[0], t2[1], t2[2], image, green);


	////image.flip_horizontally(); // i want to have the origin at the left bottom corner of the image
	////image.flip_vertically();
	//image.write_tga_file("output.tga");
	//return 0;

	//TGAImage frame(200, 200, TGAImage::RGB);
	//Vec2i pts[3] = { Vec2i(10,10), Vec2i(100, 30), Vec2i(190, 160) };
	//triangle6(pts, frame, TGAColor(255, 0, 0));
	//frame.write_tga_file("framebuffer.tga");
	//return 0;

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
	//	triangle6(screen_coords/*[0], screen_coords[1], screen_coords[2]*/, image, TGAColor(rand() % 255, rand() % 255, rand() % 255, 255));
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
}