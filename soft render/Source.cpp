#include "Header.h"

Vertex node[8] = { 
	Vertex(1, 1, 1), Vertex(1, 1, -1), Vertex(-1, 1, -1), Vertex(-1, 1, 1),
	Vertex(1, -1, 1), Vertex(1, -1, -1), Vertex(-1, -1, -1), Vertex(-1, -1, 1)
	 };
Facet facets[12] = {
	Facet(2, 1, 3, 0, 1, 0), Facet(0, 3, 1, 0, 1, 0), Facet(1, 5, 0, 1, 0, 0), Facet(4, 0, 5, 1, 0, 0),
	Facet(0, 4, 3, 0, 0, 1), Facet(7, 3, 4, 0, 0, 1), Facet(6, 7, 5, 0, -1, 0), Facet(4, 5, 7, 0, -1, 0),
	Facet(3, 7, 2, -1, 0, 0), Facet(6, 2, 7, -1, 0, 0), Facet(2, 6, 1, 0, 0, -1), Facet(5, 1, 6, 0, 0, -1), };
Vertex node2D[8];
Vertex node3D[8];
Facet facet3D[12];
double zprp = -4;	//ͶӰ�ο���z����
Vertex light = Vertex(0, 6, 0);	//��Դ����
Vertex position = Vertex(0, 0, 6);	//����������λ��
double fRotate = 0;	//��ת����
double zbuffer[800][800];	//���buffer
short gbuffer[800][800];	//�Դ�
short pic[PIC+1][PIC+1] = {};	//����
bool bAnim = 1;	//������ת

void Map(){
	//��ת����
	double rot[4][4] = {
		cos(fRotate), 0, -sin(fRotate), 0,
		0, 1, 0, 0,
		sin(fRotate), 0, cos(fRotate), 0,
		0, 0, 0, 1,
	};
	Matrix44 mRotate = Matrix44(rot);

	//ƽ�ƾ����ƶ���position
	double T[4][4] = {
		1, 0, 0, position.pos.v[0],
		0, 1, 0, position.pos.v[1],
		0, 0, 1, position.pos.v[2],
		0, 0, 0, 1,
	};
	Matrix44 mT = Matrix44(T);

	//��������ϵת�����۲�����ϵ
	double R[4][4] = {
		1, 0, 0, 0,
		0, SQR3 / 2, 1.0 / 2, 0,
		0, -1.0 / 2, SQR3 / 2, 0,
		0, 0, 0, 1,
	};
	Matrix44 mR = Matrix44(R);

	for (int i = 0; i < 8; ++i){
		//͸��ͶӰ�任
		node3D[i].pos = node2D[i].pos = mT*(mR*(mRotate * node[i].pos));
		double x = node2D[i].pos.v[0];
		double y = node2D[i].pos.v[1];
		double z = node2D[i].pos.v[2];
		node2D[i].pos.v[0] = x * zprp / (zprp - z) ;
		node2D[i].pos.v[1] = y * zprp / (zprp - z) ;
	}
}

void reshape(int width, int height){
	
}

double getBilinearFilteredPixelColor(double B, double C, int x, int y) {
	double u_ratio = B;
	double v_ratio = C;
	double u_opposite = 1 - u_ratio;
	double v_opposite = 1 - v_ratio;
	double result = (pic[x][y] * u_opposite + pic[x + 1][y] * u_ratio) * v_opposite +
		(pic[x][y + 1] * u_opposite + pic[x + 1][y + 1] * u_ratio) * v_ratio;
	return result;
}

//���
void fill(int facet, double angle){
	double xa = node2D[facets[facet].a].pos.v[0];
	double xb = node2D[facets[facet].b].pos.v[0];
	double xc = node2D[facets[facet].c].pos.v[0];
	double ya = node2D[facets[facet].a].pos.v[1];
	double yb = node2D[facets[facet].b].pos.v[1];
	double yc = node2D[facets[facet].c].pos.v[1];
	double za = node2D[facets[facet].a].pos.v[2];
	double zb = node2D[facets[facet].b].pos.v[2];
	double zc = node2D[facets[facet].c].pos.v[2];

	//��������������ʽ�ķ�ĸ
	double den1 = (yb - yc)*xa + (xc - xb)*ya + xb*yc - xc*yb;
	double den2 = (yc - ya)*xb + (xa - xc)*yb + xc*ya - xa*yc;
	double den3 = (ya - yb)*xc + (xb - xa)*yc + xa*yb - xb*ya;
	
	//������Ļ�ռ������ε���С����
	int xmin = min(xa, min(xb, xc)) * 400 + 399;
	int ymin = min(ya, min(yb, yc)) * 400 + 399;
	int xmax = max(xa, max(xb, xc)) * 400 + 401;
	int ymax = max(ya, max(yb, yc)) * 400 + 401;

	for (int j = xmin; j < xmax; ++j){
		double x, y;
		for (int k = ymin; k < ymax; ++k){
			double A, B, C;
			x = (j - 400) / 400.0;
			y = (k - 400) / 400.0;
			//������������
			A = ((yb - yc)*x + (xc - xb)*y + xb*yc - xc*yb) / den1;
			B = ((yc - ya)*x + (xa - xc)*y + xc*ya - xa*yc) / den2;
			C = ((ya - yb)*x + (xb - xa)*y + xa*yb - xb*ya) / den3;
			
			if (A>=0 && A <=1 && B>=0 && B <= 1 && C >= 0 && C <=1)	//�ж����ص��Ƿ�����������
			{

				//�������ص�����Ӧ��ά�ռ���������꣬ע�ⲻ�Ǽ�Ȩ�����
				double d = 1.0f / (A / (za - zprp) + B / (zb - zprp) + C / (zc - zprp));

				if (d < zbuffer[j][k]){
					zbuffer[j][k] = d;
					int tx, ty;
					double pix = 0;
					
					//����ӳ�����������������������
					tx = B * PIC *d / (zb - zprp);
					ty = C * PIC *d / (zc - zprp);
			
					//�����ε������棬һ�����������¿�ʼ����һ�������Ͽ�ʼ��
					if (facet % 2 == 0){
						pix = pic[tx][ty];
						//pix = getBilinearFilteredPixelColor( B ,  C , tx, ty);
					}
					else{
						tx = PIC -1 - tx;
						ty = PIC -1 - ty;
						pix = pic[tx][ty];
						//pix = getBilinearFilteredPixelColor((1 - B), (1 - C), tx, ty);
					}
					
					//�������ص������ά�ռ��е�����
					double xa3D = node3D[facets[facet].a].pos.v[0];
					double xb3D = node3D[facets[facet].b].pos.v[0];
					double xc3D = node3D[facets[facet].c].pos.v[0];
					double ya3D = node3D[facets[facet].a].pos.v[1];
					double yb3D = node3D[facets[facet].b].pos.v[1];
					double yc3D = node3D[facets[facet].c].pos.v[1];

					double x3D = A*xa3D + B*xb3D + C*xc3D;
					double y3D = A*ya3D + B*yb3D + C*yc3D;
					double z3D = d;

					//����������ߺ�ƽ�淨�����нǵ�COS
					double lAngle = ((x3D + light.pos.v[0])*facet3D[facet].vector.v[0] + (y3D + light.pos.v[1])*facet3D[facet].vector.v[1] + z3D*facet3D[facet].vector.v[2]) /
						(pow(x3D - light.pos.v[0], 2) + pow(y3D - light.pos.v[1], 2) + pow(z3D, 2) + 1);

					gbuffer[j][k] = 10 * (pix + 64) *lAngle + 0 * (pix + 64) * pow(lAngle, 10);

				}
			}
		}
	}
}

//��䣬��������
void fill2(int facet, double angle){}

void redraw(){
	glClear(GL_COLOR_BUFFER_BIT);
	glPointSize(1);
	Map();
	if (bAnim){
		fRotate += 0.01;
	}

	//��ʼ��buffer
	for (int j = 0; j < 800; ++j){
		for (int k = 0; k < 800; ++k){
			zbuffer[j][k] = 100;
			gbuffer[j][k] = 0;
		}
	}
	for (int i = 0; i < 12; ++i){
		//�����浽�ӽǵĽǶ�
		double inir = PI * 180 / 180;
		double rot[4][4] = {
			cos(fRotate + inir), 0, -sin(fRotate + inir), 0,
			0, 1, 0, 0,
			sin(fRotate + inir), 0, cos(fRotate + inir), 0,
			0, 0, 0, 1,
		};
		Matrix44 mRotate = Matrix44(rot);
		facet3D[i] = facets[i];
		facet3D[i].vector = mRotate * facets[i].vector;
		double angle = (mRotate * facets[i].vector).v[1] * 1 / 2.0 + (mRotate * facets[i].vector).v[2] * SQR3 / 2; 
		
		//�����޳���ֻ������������߳���ֵ����
		if (angle >= 0){
			fill(i, angle);
		}
	}
	
	//OpenGL����
	glBegin(GL_POINTS);
	for (int j = 0; j < 800; ++j){
		for (int k = 0; k < 800; ++k){
			if (gbuffer[j][k] > 0){
				glColor3f(gbuffer[j][k] / 256.0, gbuffer[j][k] / 256.0, gbuffer[j][k] / 256.0);
				glVertex2f((j - 400.0) / 400.0 , (k - 400.0) / 400.0);
				}
		}
	}
	glEnd();
	glutSwapBuffers();
}

void idle()
{
	glutPostRedisplay();
}

void key(unsigned char k, int x, int y)
{
	switch (k){
	case ' ':
	{
		bAnim = !bAnim;
		break;
	}
	case 'w':{
		position.pos.v[1] += 0.1;
		break;
	}
	case's':{
		position.pos.v[1] -= 0.1;
		break;
	}
	case'z':{
		position.pos.v[2] += 0.1;
		break;
	}
	case'c':{
		position.pos.v[2] -= 0.1;
		break;
	}
	case'i':{
		light.pos.v[1] += 0.1;
		break;
	}
	case'k':{
		light.pos.v[1] -= 0.1;
		break;
	}
	case'j':{
		light.pos.v[0] += 0.1;
		break;
	}
	case'l':{
		light.pos.v[0] -= 0.1;
		break;
	}
	default: break;
	}
}

//��ʼ������
void iniPic(){
	for (int j = 0; j < PIC; j++){
		for (int k = 0; k < PIC; k++){
			if (j / PICE % 2 ^ k / PICE % 2){
				pic[j][k] = 255;
			}
		}
	}
}

int main(int argc, char *argv[]){
	iniPic();
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowPosition(500, 200);
	glutInitWindowSize(800, 800);
	glutCreateWindow("Soft Render");

	glutDisplayFunc(redraw);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(key);
	glutIdleFunc(idle);

	glutMainLoop();
	return 0;
}
