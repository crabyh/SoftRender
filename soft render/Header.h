#include <iostream>
#include <vector>
#include <stdlib.h>
#include "glut/glut.h"
#include <math.h>
#include <string.h>
#include <algorithm>
#include <windows.h>
#include <math.h>

#define REAL_UNIT 1.0
#define SQR3 sqrt(3)
#define PI 3.14159265359
#define PIC 512
#define PICE 128
using namespace std;

struct Matrix44{
	double v[4][4];

	Matrix44(){}
	Matrix44(double a[4][4]){
		for (int i = 0; i < 4; ++i){
			for (int j = 0; j < 4; j++){
				v[i][j] = a[i][j];
			}
		}
	}
};

struct Matrix41{
	double v[4] = {};

	Matrix41(){}
	Matrix41(double a[4]){
		for (int i = 0; i < 4; ++i){
			v[i] = a[i];
		}
	}
	double sum(){
		return v[0] + v[1] + v[2];
	}
};

/* Position only */
struct Vertex{
	Matrix41 pos = {};

	Vertex(){}
	Vertex(double x, double y, double z){
		double a[4] = { x, y, z, 1.0f };
		pos = Matrix41(a);
	}
};

struct Facet{
	int a, b, c;
	Matrix41 vector = {};
	Facet(int a = 0, int b = 0, int c = 0, double x = 0, double y = 0, double z = 0) :a(a), b(b), c(c){
		double v[4] = { x, y, z, 1.0f };
		vector = Matrix41(v);
	}
};

Matrix44 operator *(Matrix44 m1, Matrix44 m2){
	Matrix44 res;
	for (int i = 0; i < 4; ++i){
		for (int j = 0; j < 4; ++j){
			res.v[i][j] = m1.v[i][0] * m2.v[0][j] + m1.v[i][1] * m2.v[1][j] + m1.v[i][2] * m2.v[2][j] + m1.v[i][3] * m2.v[3][j];
		}
	}
	return res;
}

Matrix41 operator *(Matrix44 m1, Matrix41 m2){
	Matrix41 res;
	for (int i = 0; i < 4; ++i){
		res.v[i] = m1.v[i][0] * m2.v[0] + m1.v[i][1] * m2.v[1] + m1.v[i][2] * m2.v[2] + m1.v[i][3] * m2.v[3];
	}
	return res;
}