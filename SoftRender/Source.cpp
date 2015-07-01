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
double zprp = -4;	//投影参考点z坐标
Vertex light = Vertex(0, 6, 0);	//光源坐标
Vertex position = Vertex(0, 0, 6);	//立方体中心位置
double fRotate = 0;	//旋转度数
double zbuffer[800][800];	//深度buffer
int gbuffer[800][800];	//显存
int pic[PIC+1][PIC+1];	//纹理
bool bAnim = 1;	//控制旋转
bool bLight = 0;	//控制光照
bool bTex = 0;	//控制纹理
bool bFill = 0;	//渲染方法 0:增量 1:重心坐标

void LoadBitMapFile(char *filename)
{
	FILE *filePtr;//文件指针
	BITMAPFILEHEADER bitmapfileheader;//bitmap文件头
	BITMAPINFOHEADER bitmapinfoheader;//bitmap信息头
	int index;
	int j=0,k=0;

	filePtr = fopen(filename,"rb");
	if(filePtr == NULL) return;
	//读入文件头
	fread(&bitmapfileheader,sizeof(BITMAPFILEHEADER),1,filePtr);
	//验证是否为bitmap文件
	if (bitmapfileheader.bfType != 0x4D42) return;
	//读入信息头
	fread(&bitmapinfoheader,sizeof(BITMAPINFOHEADER),1,filePtr);
	//验证是否为16位的位图文件
	if (bitmapinfoheader.biBitCount != 24) return;
	//将文件指针移至bitmap数据
	index = bitmapfileheader.bfOffBits;

	for(j = 0; j < bitmapinfoheader.biWidth; j++){
		for(k = 0; k <bitmapinfoheader.biHeight; k++){
			GLubyte B,G,R;
			fseek(filePtr,index,SEEK_SET);
			fread(&B,1,1,filePtr);
			pic[j][k] = B;
			fseek(filePtr,index+1,SEEK_SET);
			fread(&G,1,1,filePtr);
			pic[j][k] = (pic[j][k] << 8) | G;
			fseek(filePtr,index+2,SEEK_SET);
			fread(&R,1,1,filePtr);
			pic[j][k] = (pic[j][k] << 8) | R;
			index = index + 3;
		}
	}
	return;
}

void Map(){
	//旋转矩阵
	double rot[4][4] = {
		cos(fRotate), 0, -sin(fRotate), 0,
		0, 1, 0, 0,
		sin(fRotate), 0, cos(fRotate), 0,
		0, 0, 0, 1,
	};
	Matrix44 mRotate = Matrix44(rot);

	//平移矩阵，移动到position
	double T[4][4] = {
		1, 0, 0, position.pos.v[0],
		0, 1, 0, position.pos.v[1],
		0, 0, 1, position.pos.v[2],
		0, 0, 0, 1,
	};
	Matrix44 mT = Matrix44(T);

	//世界坐标系转换到观察坐标系
	double R[4][4] = {
		1, 0, 0, 0,
		0, SQR3 / 2, 1.0 / 2, 0,
		0, -1.0 / 2, SQR3 / 2, 0,
		0, 0, 0, 1,
	};
	Matrix44 mR = Matrix44(R);

	for (int i = 0; i < 8; ++i){
		//透视投影变换
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

//填充
void fill1(int facet, double angle){
	double xa = node2D[facets[facet].a].pos.v[0];
	double xb = node2D[facets[facet].b].pos.v[0];
	double xc = node2D[facets[facet].c].pos.v[0];
	double ya = node2D[facets[facet].a].pos.v[1];
	double yb = node2D[facets[facet].b].pos.v[1];
	double yc = node2D[facets[facet].c].pos.v[1];
	double za = node2D[facets[facet].a].pos.v[2];
	double zb = node2D[facets[facet].b].pos.v[2];
	double zc = node2D[facets[facet].c].pos.v[2];

	//计算重心坐标表达式的分母
	double den1 = (yb - yc)*xa + (xc - xb)*ya + xb*yc - xc*yb;
	double den2 = (yc - ya)*xb + (xa - xc)*yb + xc*ya - xa*yc;
	double den3 = (ya - yb)*xc + (xb - xa)*yc + xa*yb - xb*ya;
	
	//计算屏幕空间三角形的最小矩形
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
			//计算重心坐标
			A = ((yb - yc)*x + (xc - xb)*y + xb*yc - xc*yb) / den1;
			B = ((yc - ya)*x + (xa - xc)*y + xc*ya - xa*yc) / den2;
			C = ((ya - yb)*x + (xb - xa)*y + xa*yb - xb*ya) / den3;
			
			if (A>=0 && A <=1 && B>=0 && B <= 1 && C >= 0 && C <=1)	//判断像素点是否在三角形内
			{

				//计算像素点所对应三维空间点的深度坐标，注意不是简单权重相加
				double d = 1.0f / (A / (za - zprp) + B / (zb - zprp) + C / (zc - zprp));

				if (d < zbuffer[j][k]){
					zbuffer[j][k] = d;
					int tx, ty;
					int pix = 0;
					
					//纹理映射修正，根据深度坐标修正
					tx = B * PIC *d / (zb - zprp);
					ty = C * PIC *d / (zc - zprp);
			
					//正方形的两个面，一个从纹理左下开始画，一个从右上开始画
					if (facet % 2 != 0){
						pix = pic[tx][ty];
						//pix = getBilinearFilteredPixelColor( B ,  C , tx, ty);
					}
					else{
						tx = PIC -1 - tx;
						ty = PIC -1 - ty;
						pix = pic[tx][ty];
						//pix = getBilinearFilteredPixelColor((1 - B), (1 - C), tx, ty);
					}
					
					//计算像素点的在三维空间中的坐标
					double xa3D = node3D[facets[facet].a].pos.v[0];
					double xb3D = node3D[facets[facet].b].pos.v[0];
					double xc3D = node3D[facets[facet].c].pos.v[0];
					double ya3D = node3D[facets[facet].a].pos.v[1];
					double yb3D = node3D[facets[facet].b].pos.v[1];
					double yc3D = node3D[facets[facet].c].pos.v[1];

					double x3D = A*xa3D + B*xb3D + C*xc3D;
					double y3D = A*ya3D + B*yb3D + C*yc3D;
					double z3D = d;

					GLubyte r, g, b;
					r = pix & 0x0000ff;
					g = (pix >> 8) & 0x0000ff;
					b = (pix >> 16) & 0x0000ff;

					//计算入射光线和平面法向量夹角的COS
					if (bLight){
						double lAngle = ((x3D + light.pos.v[0])*facet3D[facet].vector.v[0] + (y3D + light.pos.v[1])*facet3D[facet].vector.v[1] + z3D*facet3D[facet].vector.v[2]) /
							(pow(x3D - light.pos.v[0], 2) + pow(y3D - light.pos.v[1], 2) + pow(z3D, 2) + 1);
						r = min(15 * (r + 0) * lAngle + 0 * (r + 64) * pow(lAngle, 10), 255);
						g = min(15 * (g + 0) * lAngle + 0 * (g + 64) * pow(lAngle, 10), 255);
						b = min(15 * (b + 0) * lAngle + 0 * (b + 64) * pow(lAngle, 10), 255);
					}
					
					gbuffer[j][k] = (r << 16) | (g << 8) | b;

				}
			}
		}
	}
}

void setPixel(int x, int y, double z, double s, double t, int facet, double angle){
	int tx, ty;
	int pix = 0;

	//z buffer
	if (z < zbuffer[x][y]){
		zbuffer[x][y] = z;
		tx = (int)(PIC*t);
		ty = (int)(PIC*s);

		//the two triangle of a face has the translation below
		//如果是对称三角形，坐标对称变换
		if (facet % 2 == 0){
			pix = pic[tx][ty];
			//pix = getBilinearFilteredPixelColor( B ,  C , tx, ty);
		}
		else{
			tx = PIC - 1 - tx;
			ty = PIC - 1 - ty;
			pix = pic[tx][ty];
			//pix = getBilinearFilteredPixelColor((1 - B), (1 - C), tx, ty);
		}

		//using the midpoint of a triangle to judge the light effect

		//提取三个顶点的3d坐标，取中点信息带入光照模型计算真实的rgb色彩
		double xa3D = node3D[facets[facet].a].pos.v[0];
		double xb3D = node3D[facets[facet].b].pos.v[0];
		double xc3D = node3D[facets[facet].c].pos.v[0];
		double ya3D = node3D[facets[facet].a].pos.v[1];
		double yb3D = node3D[facets[facet].b].pos.v[1];
		double yc3D = node3D[facets[facet].c].pos.v[1];

		double x3D = (xa3D + xb3D + xc3D) / 3;
		double y3D = (ya3D + yb3D + yc3D) / 3;
		double z3D = z;

		GLubyte r, g, b;
		r = pix & 0x0000ff;
		g = (pix >> 8) & 0x0000ff;
		b = (pix >> 16) & 0x0000ff;

		//º∆À„»Î…‰π‚œﬂ∫Õ∆Ω√Ê∑®œÚ¡øº–Ω«µƒCOS
		if (bLight){
			double lAngle = ((x3D + light.pos.v[0])*facet3D[facet].vector.v[0] + (y3D + light.pos.v[1])*facet3D[facet].vector.v[1] + z3D*facet3D[facet].vector.v[2]) /
				(pow(x3D - light.pos.v[0], 2) + pow(y3D - light.pos.v[1], 2) + pow(z3D, 2) + 1);
			r = fmin(15 * (r + 0) * lAngle + 0 * (r + 64) * pow(lAngle, 10), 255);
			g = fmin(15 * (g + 0) * lAngle + 0 * (g + 64) * pow(lAngle, 10), 255);
			b = fmin(15 * (b + 0) * lAngle + 0 * (b + 64) * pow(lAngle, 10), 255);
		}

		gbuffer[x][y] = (r << 16) | (g << 8) | b;

	}
}


//
//typedef struct tagEDGE{
//    double xi;
//    double dx;
//    int ymax;
//}EDGE;

void fill2(int facet, double angle){


	double x1, x2, x3, y1, y2, y3, z1, z2, z3;
	double point[3][5];

	//get the position of the point after propersition
	//提取透视后的坐标
	x1 = node2D[facets[facet].a].pos.v[0];
	x2 = node2D[facets[facet].b].pos.v[0];
	x3 = node2D[facets[facet].c].pos.v[0];

	y1 = node2D[facets[facet].a].pos.v[1];
	y2 = node2D[facets[facet].b].pos.v[1];
	y3 = node2D[facets[facet].c].pos.v[1];

	z1 = node2D[facets[facet].a].pos.v[2];
	z2 = node2D[facets[facet].b].pos.v[2];
	z3 = node2D[facets[facet].c].pos.v[2];


	//translate the position into screen
	//变换到屏幕坐标系
	point[0][0] = (int)(x1 * 400 + 400);
	point[0][1] = (y1 * 400 + 400);
	point[0][2] = z1 - zprp;
	point[1][0] = (int)(x2 * 400 + 400);
	point[1][1] = (y2 * 400 + 400);
	point[1][2] = z2 - zprp;
	point[2][0] = (int)(x3 * 400 + 400);
	point[2][1] = (y3 * 400 + 400);
	point[2][2] = z3 - zprp;

	//the texure cord of the three points
	//三个顶点纹理坐标赋值
	point[0][3] = 0.0;
	point[0][4] = 1.0;
	point[1][3] = 1.0;
	point[1][4] = 1.0;
	point[2][3] = 0.0;
	point[2][4] = 0.0;


	double tempx, tempy, tempz;
	double temps, tempt;
	int left;

	//SORT the point by its Y value
	//根据y的大小，三个顶点排序
	for (int i = 0; i<3; i++) {
		if (point[i][1]>point[0][1]){
			tempx = point[i][0];
			tempy = point[i][1];
			tempz = point[i][2];
			temps = point[i][3];
			tempt = point[i][4];
			point[i][0] = point[0][0];
			point[i][1] = point[0][1];
			point[i][2] = point[0][2];
			point[i][3] = point[0][3];
			point[i][4] = point[0][4];
			point[0][0] = tempx;
			point[0][1] = tempy;
			point[0][2] = tempz;
			point[0][3] = temps;
			point[0][4] = tempt;
		}
		if (point[i][1]<point[2][1]){
			tempx = point[i][0];
			tempy = point[i][1];
			tempz = point[i][2];
			temps = point[i][3];
			tempt = point[i][4];
			point[i][0] = point[2][0];
			point[i][1] = point[2][1];
			point[i][2] = point[2][2];
			point[i][3] = point[2][3];
			point[i][4] = point[2][4];
			point[2][0] = tempx;
			point[2][1] = tempy;
			point[2][2] = tempz;
			point[2][3] = temps;
			point[2][4] = tempt;
		}
	}



	//计算三条边的斜率，判断哪条边在底点左边
	//calculate the k of three edge
	double dx1 = (point[0][1] - point[2][1]) / (point[0][0] - point[2][0]);
	double dx2 = (point[1][1] - point[2][1]) / (point[1][0] - point[2][0]);
	//determine which edge is on the left side
	if (max(dx1, dx2)<0 || min(dx1, dx2)>0){
		left = (dx1>dx2) ? 1 : 0;
	}
	else{
		left = (dx1 <= dx2) ? 1 : 0;
	}

	//int the y value to avoid accuracy error
	//整数化，避免精度错误
	point[0][1] = (int)point[0][1];
	point[1][1] = (int)point[1][1];
	point[2][1] = (int)point[2][1];

	int ymin = point[2][1];
	int ymax = point[0][1];

	//INIT the edge table
	//初始化活动边表
	std::vector< std::list<EDGE> > slNet(ymax - ymin + 1);

	EDGE e;

	//底部边表项
	e.xi = point[2][0];
	e.ymax = ymax;
	e.dx = 1.0 / ((point[0][1] - point[2][1]) / (point[0][0] - point[2][0]));
	slNet[0].push_front(e);

	//中间点边表项
	if (point[1][1] == ymin) {
		e.xi = point[1][0];
		e.dx = 1.0 / ((point[0][1] - point[1][1]) / (point[0][0] - point[1][0]));
		slNet[0].push_front(e);
	}
	else{
		e.xi = point[1][0];
		e.ymax = ymax;
		e.dx = 1.0 / ((point[0][1] - point[1][1]) / (point[0][0] - point[1][0]));
		slNet[point[1][1] - ymin].push_front(e);


		e.xi = point[2][0];
		e.ymax = point[1][1];
		e.dx = 1.0 / ((point[1][1] - point[2][1]) / (point[1][0] - point[2][0]));
		slNet[0].push_front(e);

	}

	//divide the entile triangle into 2 part,both are pingjiao
	//fill
	//
	//1/z,is proved to be in the linear relations of x or y
	//s,t is linearial of 1/z

	//分成两个平角三角形进行填充
	//纹理坐标s,t和1/z线性相关

	EDGE e1, e2;
	double oneoverz_left, oneoverz_right, oneoverz_step, oneoverz_top, oneoverz_bottom, oneoverz;

	double soverz_top, soverz_bottom; // 上下顶点s/z
	double toverz_top, toverz_bottom; // 上下顶点t/z
	double soverz_left, soverz_right; // 左右线段s/z
	double toverz_left, toverz_right; // 左右线段t/z
	double soverz, soverz_step; // 插值s/z以及扫描线步长
	double toverz, toverz_step; // 插值t/z以及扫描线步长
	double s, t; // 要求的原始s和t


	//calculate the triangle bottom
	//add y & calculate the step
	//从下到上填充下面的三角形
	for (int y = ymin; y<point[1][1]; y++) {
		e1 = slNet[y - ymin].front();
		e2 = slNet[y - ymin].back();

		//swap
		//交换保证填充顺序正确
		if (e1.xi>e2.xi){
			e2 = slNet[y - ymin].front();
			e1 = slNet[y - ymin].back();
		}

		//calculate the step of 1/z, s/z,t/z
		//use its linear relation

		//分别计算1/z s/z t/z 和 线性梯度
		oneoverz_top = 1.0 / point[2][2];
		oneoverz_bottom = 1.0 / point[1][2];

		if (point[1][1] == point[2][1]) {
			oneoverz_left = (y - point[2][1]) * (oneoverz_bottom - oneoverz_top) + oneoverz_top;

		}
		else{
			oneoverz_left = (y - point[2][1]) * (oneoverz_bottom - oneoverz_top) / (point[1][1] - point[2][1]) + oneoverz_top;
		}
		oneoverz_bottom = 1.0 / point[0][2];
		if (point[0][1] == point[2][1]) {
			oneoverz_right = (y - point[2][1]) * (oneoverz_bottom - oneoverz_top) + oneoverz_top;
		}
		else{
			oneoverz_right = (y - point[2][1]) * (oneoverz_bottom - oneoverz_top) / (point[0][1] - point[2][1]) + oneoverz_top;
		}
		if ((int)e2.xi == (int)e1.xi) {
			//        if ((point[2][0] - point[1][0]) == 0) {
			oneoverz_step = 0;
		}
		else{
			oneoverz_step = (oneoverz_right - oneoverz_left) / (e2.xi - e1.xi);
			//            oneoverz_step = (oneoverz_right-oneoverz_left) / (point[2][0]-point[1][0]);
		}


		soverz_top = point[2][3] / point[2][2];
		soverz_bottom = point[1][3] / point[1][2];
		if (point[1][1] == point[2][1]) {
			soverz_left = (y - point[2][1]) * (soverz_bottom - soverz_top) + soverz_top;
		}
		else{
			soverz_left = (y - point[2][1]) * (soverz_bottom - soverz_top) / (point[1][1] - point[2][1]) + soverz_top;
		}
		soverz_bottom = point[0][3] / point[0][2];
		if (point[0][1] == point[2][1]) {
			soverz_right = (y - point[2][1]) * (soverz_bottom - soverz_top) + soverz_top;
		}
		else{
			soverz_right = (y - point[2][1]) * (soverz_bottom - soverz_top) / (point[0][1] - point[2][1]) + soverz_top;
		}
		if ((int)e2.xi == (int)e1.xi) {
			//        if ((point[2][0] - point[1][0]) == 0) {
			soverz_step = 0;
		}
		else{
			soverz_step = (soverz_right - soverz_left) / (e2.xi - e1.xi);
			//            soverz_step = (soverz_right-soverz_left) / (point[2][0]-point[1][0]);
		}


		toverz_top = point[2][4] / point[2][2];
		toverz_bottom = point[1][4] / point[1][2];
		if (point[1][1] == point[2][1]){
			toverz_left = (y - point[2][1]) * (toverz_bottom - toverz_top) + toverz_top;
		}
		else{
			toverz_left = (y - point[2][1]) * (toverz_bottom - toverz_top) / (point[1][1] - point[2][1]) + toverz_top;
		}
		toverz_bottom = point[0][4] / point[0][2];
		if (point[0][1] == point[2][1]) {
			toverz_right = (y - point[2][1]) * (toverz_bottom - toverz_top) + toverz_top;
		}
		else{
			toverz_right = (y - point[2][1]) * (toverz_bottom - toverz_top) / (point[0][1] - point[2][1]) + toverz_top;
		}
		if ((int)e2.xi == (int)e1.xi) {
			//        if ((point[2][0] - point[1][0]) == 0) {
			toverz_step = 0;
		}
		else{
			toverz_step = (toverz_right - toverz_left) / (e2.xi - e1.xi);

		}

		//调整边的左右关系，确保x递增顺序正确
		if (left == 1){
			oneoverz = oneoverz_right;
			soverz = soverz_right;
			toverz = toverz_right;
			soverz_step = -soverz_step;
			toverz_step = -toverz_step;
			oneoverz_step = -oneoverz_step;
		}
		else{
			oneoverz = oneoverz_left;
			soverz = soverz_left;
			toverz = toverz_left;
		}

		//fill the line by add X
		//递增x，计算每个点真实的s,t，进行填充
		for (int x = (int)e1.xi; x<(int)e2.xi; x++, oneoverz += (oneoverz_step), soverz += (soverz_step), toverz += (toverz_step)) {

			s = soverz / oneoverz;
			t = toverz / oneoverz;

			setPixel(x, y, 1.0 / oneoverz, s, t, facet, angle);

		}

		//更新下一行边表项
		if ((y + 1)<e1.ymax) {
			e1.xi += e1.dx;
			slNet[y - ymin + 1].push_front(e1);
		}
		if ((y + 1)<e2.ymax) {
			e2.xi += e2.dx;
			slNet[y - ymin + 1].push_front(e2);
		}

	}


	//the same way to the top triangle
	//同上，填充上面的
	for (int y = point[1][1]; y<ymax; y++) {
		e1 = slNet[y - ymin].front();
		e2 = slNet[y - ymin].back();

		if (e1.xi>e2.xi){
			e2 = slNet[y - ymin].front();
			e1 = slNet[y - ymin].back();
		}

		oneoverz_top = 1.0 / point[0][2];
		oneoverz_bottom = 1.0 / point[1][2];

		if (point[1][1] == point[0][1]) {
			oneoverz_left = (y - point[0][1]) * (oneoverz_bottom - oneoverz_top) + oneoverz_top;
		}
		else{
			oneoverz_left = (y - point[0][1]) * (oneoverz_bottom - oneoverz_top) / (point[1][1] - point[0][1]) + oneoverz_top;
		}
		oneoverz_bottom = 1.0 / point[2][2];
		if (point[2][1] == point[0][1]) {
			oneoverz_right = (y - point[0][1]) * (oneoverz_bottom - oneoverz_top) + oneoverz_top;
		}
		else{
			oneoverz_right = (y - point[0][1]) * (oneoverz_bottom - oneoverz_top) / (point[2][1] - point[0][1]) + oneoverz_top;
		}
		if (e2.xi == e1.xi) {
			//        if ((point[2][0] - point[1][0]) == 0) {
			oneoverz_step = 0;
		}
		else{
			oneoverz_step = (oneoverz_right - oneoverz_left) / (e2.xi - e1.xi);
			//            oneoverz_step = (oneoverz_right-oneoverz_left) / (point[2][0]-point[1][0]);
		}


		soverz_top = point[0][3] / point[0][2];
		soverz_bottom = point[1][3] / point[1][2];
		if (point[1][1] == point[0][1]) {
			soverz_left = (y - point[0][1]) * (soverz_bottom - soverz_top) + soverz_top;
		}
		else{
			soverz_left = (y - point[0][1]) * (soverz_bottom - soverz_top) / (point[1][1] - point[0][1]) + soverz_top;
		}
		soverz_bottom = point[2][3] / point[2][2];
		if (point[2][1] == point[0][1]) {
			soverz_right = (y - point[0][1]) * (soverz_bottom - soverz_top) + soverz_top;
		}
		else{
			soverz_right = (y - point[0][1]) * (soverz_bottom - soverz_top) / (point[2][1] - point[0][1]) + soverz_top;
		}
		if (e2.xi == e1.xi) {
			//        if ((point[2][0] - point[1][0]) == 0) {
			soverz_step = 0;
		}
		else{
			soverz_step = (soverz_right - soverz_left) / (e2.xi - e1.xi);
			//            soverz_step = (soverz_right-soverz_left) / (point[2][0]-point[1][0]);
		}


		toverz_top = point[0][4] / point[0][2];
		toverz_bottom = point[1][4] / point[1][2];
		if (point[1][1] == point[0][1]){
			toverz_left = (y - point[0][1]) * (toverz_bottom - toverz_top) + toverz_top;
		}
		else{
			toverz_left = (y - point[0][1]) * (toverz_bottom - toverz_top) / (point[1][1] - point[0][1]) + toverz_top;
		}
		toverz_bottom = point[2][4] / point[2][2];
		if (point[0][1] == point[2][1]) {
			toverz_right = (y - point[0][1]) * (toverz_bottom - toverz_top) + toverz_top;
		}
		else{
			toverz_right = (y - point[0][1]) * (toverz_bottom - toverz_top) / (point[2][1] - point[0][1]) + toverz_top;
		}
		if (e2.xi == e1.xi) {
			//        if ((point[2][0] - point[1][0]) == 0) {
			toverz_step = 0;
		}
		else{
			toverz_step = (toverz_right - toverz_left) / (e2.xi - e1.xi);
			//            toverz_step = (toverz_right-toverz_left) / (point[2][0]-point[1][0]);
		}

		//        oneoverz = oneoverz_left+(e1.xi-point[1][0])*oneoverz_step;
		//        soverz = soverz_left+(e1.xi-point[1][0])*soverz_step;
		//        toverz = toverz_left+(e1.xi-point[1][0])*toverz_step;


		if (left == 1){
			oneoverz = oneoverz_right;
			soverz = soverz_right;
			toverz = toverz_right;
			soverz_step = -soverz_step;
			toverz_step = -toverz_step;
			oneoverz_step = -oneoverz_step;
		}
		else{
			oneoverz = oneoverz_left;
			soverz = soverz_left;
			toverz = toverz_left;
		}


		for (int x = (int)e1.xi; x<(int)e2.xi; x++, oneoverz += (oneoverz_step), soverz += (soverz_step), toverz += (toverz_step)) {

			s = soverz / oneoverz;
			t = toverz / oneoverz;

			setPixel(x, y, 1.0 / oneoverz, s, t, facet, angle);

		}


		if ((y + 1)<e1.ymax) {
			e1.xi += e1.dx;
			slNet[y - ymin + 1].push_front(e1);
		}
		if ((y + 1)<e2.ymax) {
			e2.xi += e2.dx;
			slNet[y - ymin + 1].push_front(e2);
		}

	}

}


void getFPS()
{
	static int frame = 0, time, timebase = 0;
	static char buffer[256];

	frame++ ;
	time = glutGet(GLUT_ELAPSED_TIME);
	char mode[64];
	if (bLight){
		if (bFill)
			strcpy(mode, "light:\ton, fill: triangle");
		else
			strcpy(mode, "light:\ton, fill: Poly");
	}

	else{
		if (bFill)
			strcpy(mode, "light:\toff, fill: triangle");
		else
			strcpy(mode, "light:\toff, fill: Poly");
	}
	if (time - timebase > 1000) {
		sprintf(buffer,"FPS:%4.2f",frame*1000.0/(time-timebase));
		timebase = time;
		frame = 0;
	}

	char *c;
	glDisable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);  // 选择投影矩阵
	glPushMatrix();               // 保存原矩阵
	glLoadIdentity();             // 装入单位矩阵
	glOrtho(0,480,0,480,-1,1);    // 位置正投影
	glMatrixMode(GL_MODELVIEW);   // 选择Modelview矩阵
	glPushMatrix();               // 保存原矩阵
	glLoadIdentity();             // 装入单位矩阵
	glRasterPos2f(10,10);
	for (c=buffer; *c != '\0'; c++) {		
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *c);
	}
	glRasterPos2f(10,25);
	for (c=mode; *c != '\0'; c++) {		
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *c);
	}
	glMatrixMode(GL_PROJECTION);  // 选择投影矩阵
	glPopMatrix();                // 重置为原保存矩阵
	glMatrixMode(GL_MODELVIEW);   // 选择Modelview矩阵
	glPopMatrix();                // 重置为原保存矩阵
	//glEnable(GL_DEPTH_TEST);	
}

void redraw(){
	glClear(GL_COLOR_BUFFER_BIT);
	glPointSize(1);
	Map();
	if (bAnim){
		fRotate += 0.01;
	}

	//初始化buffer
	for (int j = 0; j < 800; ++j){
		for (int k = 0; k < 800; ++k){
			zbuffer[j][k] = 100;
			gbuffer[j][k] = 0;
		}
	}
	for (int i = 0; i < 12; ++i){
		//计算面到视角的角度
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
		
		//背面剔除，只填充向量和视线成正值的面
		if (angle >= 0){
			if (bFill){
				fill1(i, angle);
			}
			else{
				fill2(i, angle);
			}
		}
	}
	
	//OpenGL画点
	glBegin(GL_POINTS);
	for (int j = 0; j < 800; ++j){
		for (int k = 0; k < 800; ++k){
			if (gbuffer[j][k] > 0){
				glColor3f(((gbuffer[j][k]>>16)& 0x0000ff) / 256.0, ((gbuffer[j][k] >> 8) & 0x0000ff)/ 256.0, (gbuffer[j][k] & 0x0000ff) / 256.0);
				glVertex2f((j - 400.0) / 400.0 , (k - 400.0) / 400.0);
				}
		}
	}
	glEnd();

	glColor3f(1, 1, 1);
	getFPS();
	
	glutSwapBuffers();

}

void idle()
{
	glutPostRedisplay();
}

//初始化纹理
void iniPic(){
	if (bTex){
		LoadBitMapFile("Monet.bmp");
	}
	else{
		for (int j = 0; j < PIC; j++){
			for (int k = 0; k < PIC; k++){
				if (j / PICE % 2 ^ k / PICE % 2)
					pic[j][k] = (255 << 16) | (255 << 8) | 255;
				else
					pic[j][k] = (32 << 16) | (32 << 8) | 32;
			}
		}
	}
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
		if (zprp < 0)zprp += 0.1;
		break;
	}
	case'c':{
		if (zprp > -8)	zprp -= 0.1;
		break;
	}
	case'l':{
		bLight=!bLight;
		break;
	}
	case't':{
		bTex = !bTex;
		iniPic();
		break;
	}
	case'p':{
		bFill = !bFill;
		iniPic();
		break;
	}
	default: break;
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

