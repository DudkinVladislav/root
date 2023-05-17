#include <Windows.h>
#include <GL\glew.h>
#include <GL\freeglut.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <clocale>
#include <string>
using namespace std;
double tochnost = 0.0000000000001;
struct pixel
{
	double r;
	double g;
	double b;
};
void changeViewPort1(int w, int h)
{
	glViewport(0, 0, w, h);
}
void changeViewPort2(int w, int h)
{
	glViewport(0, 0, w, h);
}
void changeViewPort3(int w, int h)
{
	glViewport(0, 0, w, h);
}

void koefu(int* koefiout, int n)
{
	for (int i = 1; i <= n; i++)
		koefiout[i] = 0;
	koefiout[0] = 1;
	for (int j = 1; j <= n; j++)
		for (int i = j; i >= 1; i--)
			koefiout[i] = koefiout[i - 1] + koefiout[i];
}

void proizv_pol_odnoi_perem(double* koef1, double* koef2, double*otv, int m1, int m2)
{
	int m = m2 + m1 - 1;
	double*otvn = new double[m];
	for (int i = 0; i < m; i++)
		otvn[i] = 0;
	for (int i = 0; i < m2; i++)
	{
		for (int j = 0; j < m1; j++)
		{
			otvn[i + j] += (koef2[i] * koef1[j]);
		}
	}
	for (int i = 0; i < m; i++)
	{
		if (abs(otvn[i]) < tochnost)
		{
			otv[i] = 0;
		}
		else
		{
			otv[i] = otvn[i];
		}
	}
}

void proizv_pol_dvuh_perem(double** kof, double**otv, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (abs(kof[0][i] * kof[1][j]) < tochnost)
			{
				otv[i][j] = 0;
			}
			else {
				otv[i][j] = kof[0][i] * kof[1][j];
			}
		}
	}
}
int colvo;
int n1;
int m1;
int *n;
int *m;
int razm;
double ***otvr;
double ***otvg;
double ***otvb;
double** bx;
double** by;
pixel** orig;
pixel** colors;

void interpolate(int nachn, int nachm,int n, int m, double**otvr, double **otvg, double **otvb, int c)
{
	double mnog = 1;
	double mnog2 = 1;
	double mnogr = 1;
	double mnogb = 1;
	double mnogg = 1;
	double razn11 = 0;
	double razn22 = 0;
	int razm = (n - 1)*(m - 1) + 1;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			for (int l = 0; l < n; l++)
			{
				if (i != l)
				{
					razn11 = bx[c][i] - bx[c][l];
					mnog *= pow(razn11, double(m - 1));
				}
			}
			for (int r = 0; r < m; r++)
			{
				if (j != r)
				{
					razn22 = (by[c][j] - by[c][r]);
					mnog2 *= pow(razn22, double(n - 1));
				}
			}
			mnogg = orig[nachn+i][nachm+j].g / mnog / mnog2;
			mnogb = orig[nachn + i][nachm + j].b / mnog / mnog2;
			mnogr = orig[nachn + i][nachm + j].r / mnog / mnog2;
			int* koefurx = new int[m];
			int* koefury = new int[n];
			koefu(koefurx, m - 1);
			for (int er = 0; er < m; er++)
			{
				if (er % 2 == 1)
				{
					koefurx[er] *= -1;
				}
			}
			koefu(koefury, n - 1);
			for (int er = 0; er < n; er++)
			{
				if (er % 2 == 1)
				{
					koefury[er] *= -1;
				}
			}
			double** kofx = new double*[n - 1];
			for (int zz = 0; zz < n - 1; zz++)
			{
				kofx[zz] = new double[m];
			}
			for (int aa = 0; aa < n; aa++)
			{
				if (aa < i)
				{
					for (int bb = 0; bb < m; bb++)
					{
						kofx[aa][bb] = koefurx[bb] * pow(bx[c][aa], bb);
					}

				}
				if (aa > i)
				{
					for (int bb = 0; bb < m; bb++)
					{
						kofx[aa - 1][bb] = koefurx[bb] * pow(bx[c][aa], bb);
					}
				}
			}
			double** kofy = new double*[m - 1];
			for (int zz = 0; zz < m - 1; zz++)
			{
				kofy[zz] = new double[n];
			}
			for (int aa = 0; aa < m; aa++)
			{
				if (aa < j)
				{
					for (int bb = 0; bb < n; bb++)
					{
						kofy[aa][bb] = koefury[bb] * pow(by[c][aa], bb);
					}
				}
				if (aa > j)
					for (int bb = 0; bb < n; bb++)
					{
						kofy[aa - 1][bb] = koefury[bb] * pow(by[c][aa], bb);
					}
			}
			double **finalkoef = new double*[razm];
			double **newkoef = new double*[2];
			newkoef[0] = new double[razm];
			newkoef[1] = new double[razm];
			double *newkoefx = new double[razm];
			double *newkoefy = new double[razm];
			for (int qq = 0; qq < razm; qq++)
			{
				finalkoef[qq] = new double[razm];
				newkoefx[qq] = 0;
				newkoefy[qq] = 0;
				newkoef[0][qq] = 0;
				newkoef[1][qq] = 0;
			}
			for (int q1 = 0; q1 < razm; q1++)
				for (int q2 = 0; q2 < razm; q2++)
					finalkoef[q1][q2] = 0;
			int newrazx = m;
			for (int ab = 0; ab < m; ab++)
			{
				newkoefx[ab] = kofx[0][ab];
			}
			for (int ba = 0; ba < n; ba++)
			{
				newkoefy[ba] = kofy[0][ba];
			}
			for (int xx = 1; xx < n - 1; xx++)
			{
				proizv_pol_odnoi_perem(kofx[xx], newkoefx, newkoefx, m, newrazx);
				newrazx += (m - 1);
			}
			int newrazy = n;
			for (int yy = 1; yy < m - 1; yy++)
			{
				proizv_pol_odnoi_perem(kofy[yy], newkoefy, newkoefy, n, newrazy);
				newrazy += (n - 1);
			}
			newkoef[0] = newkoefx;
			newkoef[1] = newkoefy;

			proizv_pol_dvuh_perem(newkoef, finalkoef, razm);

			for (int q1 = 0; q1 < razm; q1++)
			{
				for (int q2 = 0; q2 < razm; q2++)
				{
					otvr[q1][q2] += finalkoef[q1][q2] * mnogr;
					otvg[q1][q2] += finalkoef[q1][q2] * mnogg;
					otvb[q1][q2] += finalkoef[q1][q2] * mnogb;
				}
			}
			mnogr = 1;
			mnogg = 1;
			mnogb = 1;
			mnog2 = 1;
			mnog = 1;
		}
	}
}
void render1()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	double r;
	double g;
	double b;
	double tochkax;
	double tochkay;
	double nn = double(n1) / 2;
	double mm = double(m1) / 2;
	double rx = 1 / nn;
	double ry = 1 / mm;
	int flagx = 1;
	int flagy = 1;
	for (int i = 0; i < n1; i++)
	{
		tochkax = i / nn - 1;

		for (int j = 0; j < m1; j++)
		{
			tochkay = j / mm - 1;
			r = orig[i][j].r;
			g = orig[i][j].g;
			b = orig[i][j].b;
			glColor3f(r, g, b);
			glBegin(GL_QUADS);
			glVertex2f(tochkax, tochkay);
			glVertex2f(tochkax, tochkay + ry * flagy);
			glVertex2f(tochkax + rx * flagx, tochkay + ry * flagy);
			glVertex2f(tochkax + rx * flagx, tochkay);
			glEnd();
			r = 0.0;
			g = 0.0;
			b = 0.0;

			flagy = 1;

		}
		flagx = 1;
	}
	glutSwapBuffers();
	glFinish();
}

void render2()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glRotatef(5.0, 1.0, 1.0, 0.0);
	glLineWidth(10.0);
	glPushMatrix();
	glRotatef(180.0, 1.0, 1.0, 0.0);
	glScalef(0.5, 0.5, 0.5);
	for (int c = 0; c < colvo; c++) 
	{
	for (int i = 0; i < n[c] ; i++)
	{
		for (int j =0; j < m[c]; j++)
		{

			glBegin(GL_LINES);
			glColor3f(0.0f, 0.0f, 0.0f);//black
			if (i > 0 && i < n[c] && j>0 && j < m[c])
			{
				glVertex3f(bx[c][i], by[c][j], orig[i][j].r);
				glVertex3f(bx[c][i], by[c][j + 1], orig[i][j + 1].r);
				glVertex3f(bx[c][i], by[c][j], orig[i][j].r);
				glVertex3f(bx[c][i + 1], by[c][j], orig[i + 1][j].r);
				glVertex3f(bx[c][i], by[c][j], orig[i][j].r);
				glVertex3f(bx[c][i], by[c][j - 1], orig[i][j - 1].r);
				glVertex3f(bx[c][i], by[c][j], orig[i][j].r);
				glVertex3f(bx[c][i - 1], by[c][j], orig[i - 1][j].r);
			}
			glEnd();
			glBegin(GL_LINES);
			glColor3f(1.0f, 0.0f, 0.0f);
			if (i > 0 && i < n[c] && j>0 && j < m[c])
			{
			glVertex3f(bx[c][i], by[c][j], colors[i][j].r);
			glVertex3f(bx[c][i], by[c][j + 1], colors[i][j + 1].r);
			glVertex3f(bx[c][i], by[c][j], colors[i][j].r);
			glVertex3f(bx[c][i + 1], by[c][j], colors[i + 1][j].r);

			glVertex3f(bx[c][i], by[c][j], colors[i][j].r);
			glVertex3f(bx[c][i], by[c][j - 1], colors[i][j - 1].r);
			glVertex3f(bx[c][i], by[c][j], colors[i][j].r);
			glVertex3f(bx[c][i - 1], by[c][j], colors[i - 1][j].r);
			}
			glEnd();
		}
	}
	for (int i = 1; i < n[c] - 1; i++)
	{
		glBegin(GL_LINES);
		glColor3f(0.0f, 0.0f, 0.0f);
		glVertex3f(bx[c][i], by[c][0], orig[i][0].r);
		glVertex3f(bx[c][i - 1], by[c][0], orig[i - 1][0].r);
		glVertex3f(bx[c][i], by[c][0], orig[i][0].r);
		glVertex3f(bx[c][i + 1], by[c][0], orig[i + 1][0].r);
		glVertex3f(bx[c][i], by[c][n[c] - 1], orig[i][n[c] - 1].r);
		glVertex3f(bx[c][i + 1], by[c][n[c] - 1], orig[i + 1][n[c] - 1].r);
		glVertex3f(bx[c][i], by[c][n[c] - 1], orig[i][n[c] - 1].r);
		glVertex3f(bx[c][i - 1], by[c][n[c] - 1], orig[i - 1][n[c] - 1].r);
		glEnd();

		glBegin(GL_LINES);
		glColor3f(1.0f, 0.0f, 0.0f);
		glVertex3f(bx[c][i], by[c][0], colors[i][0].r);
		glVertex3f(bx[c][i - 1], by[c][0], colors[i - 1][0].r);
		glVertex3f(bx[c][i], by[c][0], colors[i][0].r);
		glVertex3f(bx[c][i + 1], by[c][0], colors[i + 1][0].r);

		glVertex3f(bx[c][i], by[c][n[c] - 1], colors[i][n[c] - 1].r);
		glVertex3f(bx[c][i + 1], by[c][n[c] - 1], colors[i + 1][n[c] - 1].r);
		glVertex3f(bx[c][i], by[c][n[c] - 1], colors[i][n[c] - 1].r);
		glVertex3f(bx[c][i - 1], by[c][n[c] - 1], colors[i - 1][n[c] - 1].r);

		glEnd();
	}
	for (int i = 1; i < m[c] - 1; i++)
	{
		glBegin(GL_LINES);
		glColor3f(0.0f, 0.0f, 0.0f);
		glVertex3f(bx[c][0], by[c][i], orig[0][i].r);
		glVertex3f(bx[c][0], by[c][i - 1], orig[0][i - 1].r);
		glVertex3f(bx[c][0], by[c][i], orig[0][i].r);
		glVertex3f(bx[c][0], by[c][i + 1], orig[0][i + 1].r);
		glVertex3f(bx[c][n[c] - 1], by[c][i], orig[n[c] - 1][i].r);
		glVertex3f(bx[c][n[c] - 1], by[c][i - 1], orig[n[c] - 1][i - 1].r);
		glVertex3f(bx[c][n[c] - 1], by[c][i], orig[n[c] - 1][i].r);
		glVertex3f(bx[c][n[c] - 1], by[c][i + 1], orig[n[c] - 1][i + 1].r);
		glEnd();

		glBegin(GL_LINES);
		glColor3f(1.0f, 0.0f, 0.0f);

		glVertex3f(bx[c][0], by[c][i], colors[0][i].r);
		glVertex3f(bx[c][0], by[c][i - 1], colors[0][i - 1].r);
		glVertex3f(bx[c][0], by[c][i], colors[0][i].r);
		glVertex3f(bx[c][0], by[c][i + 1], colors[0][i + 1].r);

		glVertex3f(bx[c][n[c] - 1], by[c][i], colors[n[c] - 1][i].r);
		glVertex3f(bx[c][n[c] - 1], by[c][i - 1], colors[n[c] - 1][i - 1].r);
		glVertex3f(bx[c][n[c] - 1], by[c][i], colors[n[c] - 1][i].r);
		glVertex3f(bx[c][n[c] - 1], by[c][i + 1], colors[n[c] - 1][i + 1].r);
		glEnd();
	}
}
	glPopMatrix();
	glutSwapBuffers();
	glFinish();
}

void render3()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glRotatef(5.0, 1.0, 1.0, 0.0);
	glPushMatrix();
	glRotatef(180.0, 1.0, 1.0, 0.0);
	glScalef(0.5, 0.5, 0.5);
	glLineWidth(10.0);
	for (int c = 0; c < colvo; c++)
	{
		for (int i = 1; i < n[c] - 1; i++)
		{
			for (int j = 1; j < m[c] - 1; j++)
			{
				glBegin(GL_LINES);
				glColor3f(0.0f, 0.0f, 0.0f);//black
				glVertex3f(bx[c][i], by[c][j], orig[i][j].g);
				glVertex3f(bx[c][i], by[c][j + 1], orig[i][j + 1].g);
				glVertex3f(bx[c][i], by[c][j], orig[i][j].g);
				glVertex3f(bx[c][i + 1], by[c][j], orig[i + 1][j].g);
				glVertex3f(bx[c][i], by[c][j], orig[i][j].g);
				glVertex3f(bx[c][i], by[c][j - 1], orig[i][j - 1].g);
				glVertex3f(bx[c][i], by[c][j], orig[i][j].g);
				glVertex3f(bx[c][i - 1], by[c][j], orig[i - 1][j].g);
				glEnd();
				glBegin(GL_LINES);
				glColor3f(0.0f, 1.0f, 0.0f);

				glVertex3f(bx[c][i], by[c][j], colors[i][j].g);
				glVertex3f(bx[c][i], by[c][j + 1], colors[i][j + 1].g);
				glVertex3f(bx[c][i], by[c][j], colors[i][j].g);
				glVertex3f(bx[c][i + 1], by[c][j], colors[i + 1][j].g);

				glVertex3f(bx[c][i], by[c][j], colors[i][j].g);
				glVertex3f(bx[c][i], by[c][j - 1], colors[i][j - 1].g);
				glVertex3f(bx[c][i], by[c][j], colors[i][j].g);
				glVertex3f(bx[c][i - 1], by[c][j], colors[i - 1][j].g);
				glEnd();
			}
		}
		for (int i = 1; i < n[c] - 1; i++)
		{
			glBegin(GL_LINES);
			glColor3f(0.0f, 0.0f, 0.0f);
			glVertex3f(bx[c][i], by[c][0], orig[i][0].g);
			glVertex3f(bx[c][i - 1], by[c][0], orig[i - 1][0].g);
			glVertex3f(bx[c][i], by[c][0], orig[i][0].g);
			glVertex3f(bx[c][i + 1], by[c][0], orig[i + 1][0].g);
			glVertex3f(bx[c][i], by[c][n[c] - 1], orig[i][n[c] - 1].g);
			glVertex3f(bx[c][i + 1], by[c][n[c] - 1], orig[i + 1][n[c] - 1].g);
			glVertex3f(bx[c][i], by[c][n[c] - 1], orig[i][n[c] - 1].g);
			glVertex3f(bx[c][i - 1], by[c][n[c] - 1], orig[i - 1][n[c] - 1].g);
			glEnd();

			glBegin(GL_LINES);
			glColor3f(0.0f, 1.0f, 0.0f);
			glVertex3f(bx[c][i], by[c][0], colors[i][0].g);
			glVertex3f(bx[c][i - 1], by[c][0], colors[i - 1][0].g);
			glVertex3f(bx[c][i], by[c][0], colors[i][0].g);
			glVertex3f(bx[c][i + 1], by[c][0], colors[i + 1][0].g);

			glVertex3f(bx[c][i], by[c][n[c] - 1], colors[i][n[c] - 1].g);
			glVertex3f(bx[c][i + 1], by[c][n[c] - 1], colors[i + 1][n[c] - 1].g);
			glVertex3f(bx[c][i], by[c][n[c] - 1], colors[i][n[c] - 1].g);
			glVertex3f(bx[c][i - 1], by[c][n[c] - 1], colors[i - 1][n[c] - 1].g);
			glEnd();
		}
		for (int i = 1; i < m[c] - 1; i++)
		{
			glBegin(GL_LINES);
			glColor3f(0.0f, 0.0f, 0.0f);

			glVertex3f(bx[c][0], by[c][i], orig[0][i].g);
			glVertex3f(bx[c][0], by[c][i - 1], orig[0][i - 1].g);
			glVertex3f(bx[c][0], by[c][i], orig[0][i].g);
			glVertex3f(bx[c][0], by[c][i + 1], orig[0][i + 1].g);
			glVertex3f(bx[c][n[c] - 1], by[c][i], orig[n[c] - 1][i].g);
			glVertex3f(bx[c][n[c] - 1], by[c][i - 1], orig[n[c] - 1][i - 1].g);
			glVertex3f(bx[c][n[c] - 1], by[c][i], orig[n[c] - 1][i].g);
			glVertex3f(bx[c][n[c] - 1], by[c][i + 1], orig[n[c] - 1][i + 1].g);
			glEnd();

			glBegin(GL_LINES);
			glColor3f(0.0f, 1.0f, 0.0f);

			glVertex3f(bx[c][0], by[c][i], colors[0][i].g);
			glVertex3f(bx[c][0], by[c][i - 1], colors[0][i - 1].g);
			glVertex3f(bx[c][0], by[c][i], colors[0][i].g);
			glVertex3f(bx[c][0], by[c][i + 1], colors[0][i + 1].g);

			glVertex3f(bx[c][n[c] - 1], by[c][i], colors[n[c] - 1][i].g);
			glVertex3f(bx[c][n[c] - 1], by[c][i - 1], colors[n[c] - 1][i - 1].g);
			glVertex3f(bx[c][n[c] - 1], by[c][i], colors[n[c] - 1][i].g);
			glVertex3f(bx[c][n[c] - 1], by[c][i + 1], colors[n[c] - 1][i + 1].g);
			glEnd();
		}
	}
	glPopMatrix();
	glutSwapBuffers();
	glFinish();
}

void render4()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glRotatef(5.0, 1.0, 1.0, 0.0);
	glPushMatrix();
	glRotatef(180.0, 1.0, 1.0, 0.0);
	glScalef(0.5, 0.5, 0.5);
	glLineWidth(10.0);
	for (int c = 0; c < colvo; c++)
	{
		for (int i = 1; i < n[c] - 1; i++)
		{
			for (int j = 1; j < m[c] - 1; j++)
			{
				glBegin(GL_LINES);
				glColor3f(0.0f, 0.0f, 0.0f);//black
				glVertex3f(bx[c][i], by[c][j], orig[i][j].b);
				glVertex3f(bx[c][i], by[c][j + 1], orig[i][j + 1].b);
				glVertex3f(bx[c][i], by[c][j], orig[i][j].b);
				glVertex3f(bx[c][i + 1], by[c][j], orig[i + 1][j].b);
				glVertex3f(bx[c][i], by[c][j], orig[i][j].b);
				glVertex3f(bx[c][i], by[c][j - 1], orig[i][j - 1].b);
				glVertex3f(bx[c][i], by[c][j], orig[i][j].b);
				glVertex3f(bx[c][i - 1], by[c][j], orig[i - 1][j].b);
				glEnd();
				glBegin(GL_LINES);
				glColor3f(0.0f, 0.0f, 1.0f);

				glVertex3f(bx[c][i], by[c][j], colors[i][j].b);
				glVertex3f(bx[c][i], by[c][j + 1], colors[i][j + 1].b);
				glVertex3f(bx[c][i], by[c][j], colors[i][j].b);
				glVertex3f(bx[c][i + 1], by[c][j], colors[i + 1][j].b);

				glVertex3f(bx[c][i], by[c][j], colors[i][j].b);
				glVertex3f(bx[c][i], by[c][j - 1], colors[i][j - 1].b);
				glVertex3f(bx[c][i], by[c][j], colors[i][j].b);
				glVertex3f(bx[c][i - 1], by[c][j], colors[i - 1][j].b);
				glEnd();
			}
		}
		for (int i = 1; i < n[c] - 1; i++)
		{
			glBegin(GL_LINES);
			glColor3f(0.0f, 0.0f, 0.0f);
			glVertex3f(bx[c][i], by[c][0], orig[i][0].b);
			glVertex3f(bx[c][i - 1], by[c][0], orig[i - 1][0].b);
			glVertex3f(bx[c][i], by[c][0], orig[i][0].b);
			glVertex3f(bx[c][i + 1], by[c][0], orig[i + 1][0].b);

			glVertex3f(bx[c][i], by[c][n[c] - 1], orig[i][n[c] - 1].b);
			glVertex3f(bx[c][i + 1], by[c][n[c] - 1], orig[i + 1][n[c] - 1].b);
			glVertex3f(bx[c][i], by[c][n[c] - 1], orig[i][n[c] - 1].b);
			glVertex3f(bx[c][i - 1], by[c][n[c] - 1], orig[i - 1][n[c] - 1].b);
			glEnd();

			glBegin(GL_LINES);
			glColor3f(0.0f, 0.0f, 1.0f);
			glVertex3f(bx[c][i], by[c][0], colors[i][0].b);
			glVertex3f(bx[c][i - 1], by[c][0], colors[i - 1][0].b);
			glVertex3f(bx[c][i], by[c][0], colors[i][0].b);
			glVertex3f(bx[c][i + 1], by[c][0], colors[i + 1][0].b);

			glVertex3f(bx[c][i], by[c][n[c] - 1], colors[i][n[c] - 1].b);
			glVertex3f(bx[c][i + 1], by[c][n[c] - 1], colors[i + 1][n[c] - 1].b);
			glVertex3f(bx[c][i], by[c][n[c] - 1], colors[i][n[c] - 1].b);
			glVertex3f(bx[c][i - 1], by[c][n[c] - 1], colors[i - 1][n[c] - 1].b);

			glEnd();
		}
		for (int i = 1; i < m[c] - 1; i++)
		{
			glBegin(GL_LINES);
			glColor3f(0.0f, 0.0f, 0.0f);
			glVertex3f(bx[c][0], by[c][i], orig[0][i].b);
			glVertex3f(bx[c][0], by[c][i - 1], orig[0][i - 1].b);
			glVertex3f(bx[c][0], by[c][i], orig[0][i].b);
			glVertex3f(bx[c][0], by[c][i + 1], orig[0][i + 1].b);
			glVertex3f(bx[c][n[c] - 1], by[c][i], orig[n[c] - 1][i].b);
			glVertex3f(bx[c][n[c] - 1], by[c][i - 1], orig[n[c] - 1][i - 1].b);
			glVertex3f(bx[c][n[c] - 1], by[c][i], orig[n[c] - 1][i].b);
			glVertex3f(bx[c][n[c] - 1], by[c][i + 1], orig[n[c] - 1][i + 1].b);
			glEnd();

			glBegin(GL_LINES);
			glColor3f(0.0f, 0.0f, 1.0f);

			glVertex3f(bx[c][0], by[c][i], colors[0][i].b);
			glVertex3f(bx[c][0], by[c][i - 1], colors[0][i - 1].b);
			glVertex3f(bx[c][0], by[c][i], colors[0][i].b);
			glVertex3f(bx[c][0], by[c][i + 1], colors[0][i + 1].b);

			glVertex3f(bx[c][n[c] - 1], by[c][i], colors[n[c] - 1][i].b);
			glVertex3f(bx[c][n[c] - 1], by[c][i - 1], colors[n[c] - 1][i - 1].b);
			glVertex3f(bx[c][n[c] - 1], by[c][i], colors[n[c] - 1][i].b);
			glVertex3f(bx[c][n[c] - 1], by[c][i + 1], colors[n[c] - 1][i + 1].b);
			glEnd();
		}
	}
	glPopMatrix();
	glutSwapBuffers();
	glFinish();
}

void render()
{
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	double tochkax;
	double tochkay;
	double nn = double(n1) / 2;
	double mm = double(m1) / 2;
	double rx = 1 / nn;
	double ry = 1 / mm;
	for (int c = 0; c < colvo; c++)
	{
		razm = (n[c] - 1)*(m[c] - 1)+1;
		int in = 0;
		for (int i = 0; i < n[c]; i++)
		{
			int im = 0;
			tochkax = tochkax = bx[c][i];
			for (int j = 0; j < m[c]; j++)
			{
				tochkay = by[c][j];
				double clr = 0.0;
				double clg = 0.0;
				double clb = 0.0;
				double xyn = 1.0;
				for (int k = 0; k < razm; k++)
				{
					for (int l = 0; l < razm; l++)
					{
						xyn *= pow(tochkax, (n[c] - 1)*(m[c] - 1) - k);
						xyn *= pow(tochkay, (n[c] - 1)*(m[c] - 1) - l);

						if (abs(otvr[c][k][l]) > tochnost)
						{
							clr += otvr[c][k][l] * xyn;
						}
						if (abs(otvg[c][k][l]) > tochnost)
						{
							clg += otvg[c][k][l] * xyn;
						}
						if (abs(otvb[c][k][l]) > tochnost)
						{
							clb += otvb[c][k][l] * xyn;
						}
						xyn = 1.0;
					}
				}
				if (clr > 1)
					clr = 1;
				if (clr < 0)
					clr = 0;
				if (clg > 1)
					clg = 1;
				if (clg < 0)
					clg = 0;
				if (clb > 1)
					clb = 1;
				if (clb < 0)
					clb = 0;
				colors[in][im].r = clr;
				colors[in][im].g = clg;
				colors[in][im].b = clb;
				glColor3f(clr, clg, clb);
				glBegin(GL_QUADS);
				glVertex2f(tochkax, tochkay);
				glVertex2f(tochkax, tochkay + ry);
				glVertex2f(tochkax + rx, tochkay + ry);
				glVertex2f(tochkax + rx, tochkay);
				glEnd();
				clr = 0.0;
				clg = 0.0;
				clb = 0.0;
				im++;
			}
			in++;
		}
	}
	double maxnevyazr = 0;
	double maxnevyazg = 0;
	double maxnevyazb = 0;
	int xnevr = -1;
	int ynevr = -1;
	int xnevg = -1;
	int ynevg = -1;
	int xnevb = -1;
	int ynevb = -1;
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < m1; j++)
		{
			if (abs(orig[i][j].r - colors[i][j].r) > maxnevyazr)
			{
				maxnevyazr = abs(orig[i][j].r - colors[i][j].r);
				xnevr = i;
				ynevr = j;
			}
			if (abs(orig[i][j].g - colors[i][j].g) > maxnevyazg)
			{
				maxnevyazg = abs(orig[i][j].g - colors[i][j].g);
				xnevg = i;
				ynevg = j;
			}
			if (abs(orig[i][j].b - colors[i][j].b) > maxnevyazb)
			{
				maxnevyazb = abs(orig[i][j].b - colors[i][j].b);
				xnevb = i;
				ynevb = j;
			}
		}
	}
	cout << "Координата х разницы=" << xnevr << " Координата у разницы=" << ynevr << endl;
	cout << "Максимальная разница по красному цвету=" << maxnevyazr << endl;
	cout << "Координата х разницы=" << xnevg << " Координата у разницы=" << ynevg << endl;
	cout << "Максимальная разница по зелёному цвету=" << maxnevyazg << endl;
	cout << "Координата х разницы=" << xnevb << " Координата у разницы=" << ynevb << endl;
	cout << "Максимальная разница по синему цвету=" << maxnevyazb << endl;
	glutSwapBuffers();
	glFinish();
}

int main(int argc, char* argv[]) {
	setlocale(LC_ALL, "Russian");
	cout << "Введите количество пикселей по х:";
	cin >> n1;
	cout << "Введите количество пикселей по у:";
	cin >> m1;
	orig = new pixel*[n1];
	colors = new pixel*[n1];
	int cx = n1 / 3;
	int cy = m1 / 3;
	int ostn = n1 % 3;
	int ostm = m1 % 3;
	int flagx = 0;
	int flagy = 0;
	if (ostn > 0)
	{
		cx++;
	}
	if (ostm > 0)
	{
		cy++;
	}
	if (ostn < 2)
	{
		flagx++;
	}
	if (ostm < 2)
	{
		flagy++;
	}
	if (ostn == 0)
		flagx = 0;
	if (ostm == 0)
		flagy = 0;
	colvo = cx * cy;
	n = new int[colvo];
	m = new int[colvo];
	int **start = new int*[colvo];
	int **kon = new int*[colvo];
	for (int i = 0; i < colvo; i++)
	{
		start[i] = new int[2];
		kon[i] = new int[2];
	}
	start[0][0] = 0;
	start[0][1] = 0;
	kon[0][0] = 2;
	kon[0][1] = 2;
	int coln=0;
	int pribx = ostn;
	int priby = ostm;
	if (ostn == 0)
		pribx = 3;
	if (ostm == 0)
		priby = 3;
	for (int i = 1; i < colvo; i++)
	{
		if (i%cy != 0 && i < (colvo - cx))
		{
			start[i][0] = start[i-1][0];
			kon[i][0] = kon[i-1][0];
		}
		if (i%cy == 0)
		{
			start[i][0] = start[i - 1][0]+3;
			kon[i][0] = kon[i - 1][0]+3;
		}
		if (i >= (colvo - cx))
		{
			start[i][0] = start[i - 1][0];
			kon[i][0] = kon[i - 1][0];
		}
		if ((i%cy == 0) && (i >= (colvo - cx)))
		{
			start[i][0] = start[i - 1][0]+3;
			kon[i][0] = kon[i - 1][0] + pribx;
		}
	}
	for (int i = 1; i < cy; i++)
	{
		start[i][1] = start[i - 1][1]+3;
		kon[i][1] = kon[i - 1][1] + 3;
		if(i%cy==(cy-1))
		{
			kon[i][1] = kon[i - 1][1] + priby;
		}
	}
	for (int i = cy; i < colvo; i++)
	{
		start[i][1] = start[i - cy][1];
		kon[i][1] = kon[i - cy][1];
	}
	for (int i = 0; i < colvo; i++)
	{
		n[i] = kon[i][0] - start[i][0] + 1;
		m[i] = kon[i][1] - start[i][1] + 1;
	}
	bx = new double*[colvo];
	by = new double*[colvo];
	for (int i = 0; i < n1; i++)
	{
		orig[i] = new pixel[m1];
		colors[i] = new pixel[m1];
	}
	double *bxs = new double[n1];
	double *bys = new double[m1];
	double st = -1;
	double shag1 = 2 / double(n1);
	double shag2 = 2 / double(m1);
	for (int nom1 = 0; nom1 < n1; nom1++)
	{
		for (int nom2 = 0; nom2 < m1; nom2++)
		{
			//cout << "Введите значения красного, зелёного и синего цветов для пикселя " << nom1 << " " << nom2 << endl;
			/*cin >> orig[nom1][nom2].r;
			cin >> orig[nom1][nom2].g;
			cin >> orig[nom1][nom2].b;*/
			orig[nom1][nom2].r = pow(sin((-1 + nom1 * shag1)*(-1 + nom2 * shag2)), 10);
			orig[nom1][nom2].g = pow(cos((-1 + nom1 * shag1)*(-1 + nom2 * shag2)), 15);
			orig[nom1][nom2].b = pow((-1 + nom1 * shag1)*(-1 + nom2 * shag2), 5); 
			if (orig[nom1][nom2].r < 0)
				orig[nom1][nom2].r = 0;
			if (orig[nom1][nom2].g < 0)
				orig[nom1][nom2].g = 0;
			if (orig[nom1][nom2].b < 0)
				orig[nom1][nom2].b = 0;
			if (orig[nom1][nom2].r > 1)
				orig[nom1][nom2].r = 1;
			if (orig[nom1][nom2].g > 1)
				orig[nom1][nom2].g = 1;
			if (orig[nom1][nom2].b > 1)
				orig[nom1][nom2].b = 1;/*
			cout << orig[nom1][nom2].r << endl;
			cout << orig[nom1][nom2].g << endl;
			cout << orig[nom1][nom2].b << endl;*/
		}
		bxs[nom1] = -1 + nom1 * shag1;
	}
		/*orig[0][0].r = 0;  orig[0][0].g = 0.5;	orig[0][0].b = 1;
		orig[0][1].r = 1;	orig[0][1].g = 1;	orig[0][1].b = 1;
		orig[0][2].r = 1;	orig[0][2].g = 1;	orig[0][2].b = 1;
		orig[0][3].r = 0;	orig[0][3].g = 0;	orig[0][3].b = 0;
		orig[0][4].r = 0;	orig[0][4].g = 0;	orig[0][4].b = 0;
		orig[0][5].r = 0;	orig[0][5].g = 0;	orig[0][5].b = 0;
		orig[0][6].r = 1;	orig[0][6].g = 1;	orig[0][6].b = 1;
		orig[0][7].r = 0;	orig[0][7].g = 0;	orig[0][7].b = 0;
		orig[0][8].r = 0;	orig[0][8].g = 0;	orig[0][8].b = 0;
		orig[0][9].r = 0;	orig[0][9].g = 0;	orig[0][9].b = 0;
		orig[0][10].r = 0;	orig[0][10].g = 0;	orig[0][10].b = 0;
		orig[0][11].r = 1;	orig[0][11].g = 1;	orig[0][11].b = 1;
		orig[0][12].r = 1;	orig[0][12].g = 1;	orig[0][12].b = 1;
		orig[0][13].r = 1;	orig[0][13].g = 1;	orig[0][13].b = 1;
		orig[0][14].r = 1;	orig[0][14].g = 1;	orig[0][14].b = 1;
	
		orig[1][0].r = 1;	orig[1][0].g = 1;	orig[1][0].b = 1;
		orig[1][1].r = 1;	orig[1][1].g = 1;	orig[1][1].b = 1;
		orig[1][2].r = 0;	orig[1][2].g = 0;	orig[1][2].b = 0;
		orig[1][3].r = 1;	orig[1][3].g = 1;	orig[1][3].b = 1;
		orig[1][4].r = 1;	orig[1][4].g = 1;	orig[1][4].b = 1;
		orig[1][5].r = 0;	orig[1][5].g = 0;	orig[1][5].b = 0;
		orig[1][6].r = 0;	orig[1][6].g = 0;	orig[1][6].b = 0;
		orig[1][7].r = 1;	orig[1][7].g = 1;	orig[1][7].b = 1;
		orig[1][8].r = 1;	orig[1][8].g = 1;	orig[1][8].b = 1;
		orig[1][9].r = 1;	orig[1][9].g = 1;	orig[1][9].b = 1;
		orig[1][10].r = 1;	orig[1][10].g = 1;	orig[1][10].b = 1;
		orig[1][11].r = 0;	orig[1][11].g = 0;	orig[1][11].b = 0;
		orig[1][12].r = 0;	orig[1][12].g = 0;	orig[1][12].b = 0;
		orig[1][13].r = 1;	orig[1][13].g = 1;	orig[1][13].b = 1;
		orig[1][14].r = 1;	orig[1][14].g = 1;	orig[1][14].b = 1;
	
		orig[2][0].r = 1;	orig[2][0].g = 1;	orig[2][0].b = 1;
		orig[2][1].r = 0;	orig[2][1].g = 0;	orig[2][1].b = 0;
		orig[2][2].r = 0;	orig[2][2].g = 0;	orig[2][2].b = 0;
		orig[2][3].r = 1;	orig[2][3].g = 1;	orig[2][3].b = 1;
		orig[2][4].r = 1;	orig[2][4].g = 1;	orig[2][4].b = 1;
		orig[2][5].r = 1;	orig[2][5].g = 1;	orig[2][5].b = 1;
		orig[2][6].r = 1;	orig[2][6].g = 1;	orig[2][6].b = 1;
		orig[2][7].r = 1;	orig[2][7].g = 1;	orig[2][7].b = 1;
		orig[2][8].r = 1;	orig[2][8].g = 1;	orig[2][8].b = 1;
		orig[2][9].r = 1;	orig[2][9].g = 1;	orig[2][9].b = 1;
		orig[2][10].r = 1;	orig[2][10].g = 1;	orig[2][10].b = 1;
		orig[2][11].r = 1;	orig[2][11].g = 1;	orig[2][11].b = 1;
		orig[2][12].r = 1;	orig[2][12].g = 1;	orig[2][12].b = 1;
		orig[2][13].r = 0;	orig[2][13].g = 0;	orig[2][13].b = 0;
		orig[2][14].r = 1;	orig[2][14].g = 1;	orig[2][14].b = 1;
	
		orig[3][0].r= 0;	orig[3][0].g = 0;	orig[3][0].b = 0;
		orig[3][1].r = 0;	orig[3][1].g = 0;	orig[3][1].b = 0;
		orig[3][2].r = 1;	orig[3][2].g = 1;	orig[3][2].b = 1;
		orig[3][3].r = 1;	orig[3][3].g = 1;	orig[3][3].b = 1;
		orig[3][4].r = 0;	orig[3][4].g = 0;	orig[3][4].b = 0;
		orig[3][5].r = 0;	orig[3][5].g = 0;	orig[3][5].b = 0;
		orig[3][6].r = 1;	orig[3][6].g = 1;	orig[3][6].b = 1;
		orig[3][7].r = 0;	orig[3][7].g = 0;	orig[3][7].b = 0;
		orig[3][8].r = 0;	orig[3][8].g = 0;	orig[3][8].b = 0;
		orig[3][9].r = 0;	orig[3][9].g = 0;	orig[3][9].b = 0;
		orig[3][10].r = 1;	orig[3][10].g = 1;	orig[3][10].b = 1;
		orig[3][11].r = 1;	orig[3][11].g = 1;	orig[3][11].b = 1;
		orig[3][12].r = 1;	orig[3][12].g = 1;	orig[3][12].b = 1;
		orig[3][13].r = 0;	orig[3][13].g = 0;	orig[3][13].b = 0;
		orig[3][14].r = 1;	orig[3][14].g = 1;	orig[3][14].b = 1;
	
		orig[4][0].r = 0;	orig[4][0].g = 0;	orig[4][0].b = 0;
		orig[4][1].r = 1;	orig[4][1].g = 1;	orig[4][1].b = 1;
		orig[4][2].r = 1;	orig[4][2].g = 1;	orig[4][2].b = 1;
		orig[4][3].r = 0;	orig[4][3].g = 0;	orig[4][3].b = 0;
		orig[4][4].r = 0;	orig[4][4].g = 0;	orig[4][4].b = 0;
		orig[4][5].r = 1;	orig[4][5].g = 1;	orig[4][5].b = 1;
		orig[4][6].r = 1;	orig[4][6].g = 1;	orig[4][6].b = 1;
		orig[4][7].r = 0;	orig[4][7].g = 0;	orig[4][7].b = 0;
		orig[4][8].r = 1;	orig[4][8].g = 0;	orig[4][8].b = 0;
		orig[4][9].r = 0;	orig[4][9].g = 0;	orig[4][9].b = 0;
		orig[4][10].r = 1;	orig[4][10].g = 1;	orig[4][10].b = 1;
		orig[4][11].r = 1;	orig[4][11].g = 1;	orig[4][11].b= 1;
		orig[4][12].r = 1;	orig[4][12].g = 1;	orig[4][12].b = 1;
		orig[4][13].r = 1;	orig[4][13].g = 1;	orig[4][13].b = 1;
		orig[4][14].r = 0;	orig[4][14].g = 0;	orig[4][14].b = 0;
	
		orig[5][0].r = 0;	orig[5][0].g = 0;	orig[5][0].b = 0;
		orig[5][1].r = 1;	orig[5][1].g = 1;	orig[5][1].b = 1;
		orig[5][2].r = 0;	orig[5][2].g = 0;	orig[5][2].b = 0;
		orig[5][3].r = 1;	orig[5][3].g = 1;	orig[5][3].b = 1;
		orig[5][4].r = 0;	orig[5][4].g = 0;	orig[5][4].b = 0;
		orig[5][5].r = 1;	orig[5][5].g = 1;	orig[5][5].b = 1;
		orig[5][6].r = 1;	orig[5][6].g = 1;	orig[5][6].b = 1;
		orig[5][7].r = 0;	orig[5][7].g = 0;	orig[5][7].b = 0;
		orig[5][8].r = 0;	orig[5][8].g = 0;	orig[5][8].b = 0;
		orig[5][9].r = 0;	orig[5][9].g = 0;	orig[5][9].b = 0;
		orig[5][10].r = 1;	orig[5][10].g = 1;	orig[5][10].b = 1;
		orig[5][11].r = 1;	orig[5][11].g = 1;	orig[5][11].b = 1;
		orig[5][12].r = 1;	orig[5][12].g = 1;	orig[5][12].b = 1;
		orig[5][13].r = 1;	orig[5][13].g = 1;	orig[5][13].b = 1;
		orig[5][14].r = 0;	orig[5][14].g = 0;	orig[5][14].b = 0;
	
		orig[6][0].r = 0;	orig[6][0].g = 0;	orig[6][0].b = 0;
		orig[6][1].r = 1;	orig[6][1].g = 1;	orig[6][1].b = 1;
		orig[6][2].r = 0;	orig[6][2].g = 0;	orig[6][2].b = 0;
		orig[6][3].r = 0;	orig[6][3].g = 0;	orig[6][3].b = 0;
		orig[6][4].r = 0;	orig[6][4].g = 0;	orig[6][4].b = 0;
		orig[6][5].r = 1;	orig[6][5].g = 1;	orig[6][5].b = 1;
		orig[6][6].r = 1;	orig[6][6].g = 1;	orig[6][6].b = 1;
		orig[6][7].r = 1;	orig[6][7].g = 1;	orig[6][7].b = 1;
		orig[6][8].r = 1;	orig[6][8].g = 1;	orig[6][8].b = 1;
		orig[6][9].r = 1;	orig[6][9].g = 1;	orig[6][9].b = 1;
		orig[6][10].r = 1;	orig[6][10].g = 1;	orig[6][10].b = 1;
		orig[6][11].r = 1;	orig[6][11].g = 1;	orig[6][11].b = 1;
		orig[6][12].r = 1;	orig[6][12].g = 1;	orig[6][12].b = 1;
		orig[6][13].r = 1;	orig[6][13].g = 1;	orig[6][13].b = 1;
		orig[6][14].r = 0;	orig[6][14].g = 0;	orig[6][14].b = 0;
	
		orig[7][0].r = 0;	orig[7][0].g = 0;	orig[7][0].b = 0;
		orig[7][1].r = 1;	orig[7][1].g = 1;	orig[7][1].b = 1;
		orig[7][2].r = 0;	orig[7][2].g = 0;	orig[7][2].b = 0;
		orig[7][3].r = 1;	orig[7][3].g = 1;	orig[7][3].b = 1;
		orig[7][4].r = 0;	orig[7][4].g = 0;	orig[7][4].b = 0;
		orig[7][5].r = 1;	orig[7][5].g = 1;	orig[7][5].b = 1;
		orig[7][6].r = 0;	orig[7][6].g = 0;	orig[7][6].b = 0;
		orig[7][7].r = 1;	orig[7][7].g = 1;	orig[7][7].b = 1;
		orig[7][8].r = 1;	orig[7][8].g = 1;	orig[7][8].b = 1;
		orig[7][9].r = 1;	orig[7][9].g = 1;	orig[7][9].b = 1;
		orig[7][10].r = 1;	orig[7][10].g = 1;	orig[7][10].b = 1;
		orig[7][11].r = 1;	orig[7][11].g = 1;	orig[7][11].b = 1;
		orig[7][12].r = 1;	orig[7][12].g = 1;	orig[7][12].b = 1;
		orig[7][13].r = 1;	orig[7][13].g = 1;	orig[7][13].b = 1;
		orig[7][14].r = 0;	orig[7][14].g = 0;	orig[7][14].b = 0;
	
		orig[8][0].r = 0;	orig[8][0].g = 0;	orig[8][0].b = 0;
		orig[8][1].r = 1;	orig[8][1].g = 1;	orig[8][1].b = 1;
		orig[8][2].r = 0;	orig[8][2].g = 0;	orig[8][2].b = 0;
		orig[8][3].r = 0;	orig[8][3].g = 0;	orig[8][3].b = 0;
		orig[8][4].r = 0;	orig[8][4].g = 0;	orig[8][4].b = 0;
		orig[8][5].r = 1;	orig[8][5].g = 1;	orig[8][5].b = 1;
		orig[8][6].r = 0;	orig[8][6].g = 0;	orig[8][6].b = 0;
		orig[8][7].r = 0;	orig[8][7].g = 0;	orig[8][7].b = 0;
		orig[8][8].r = 1;	orig[8][8].g = 1;	orig[8][8].b = 1;
		orig[8][9].r = 1;	orig[8][9].g = 1;	orig[8][9].b = 1;
		orig[8][10].r = 1;	orig[8][10].g = 1;	orig[8][10].b = 1;
		orig[8][11].r = 1;	orig[8][11].g = 1;	orig[8][11].b = 1;
		orig[8][12].r = 1;	orig[8][12].g = 1;	orig[8][12].b = 1;
		orig[8][13].r = 1;	orig[8][13].g = 1;	orig[8][13].b = 1;
		orig[8][14].r = 0;	orig[8][14].g = 0;	orig[8][14].b = 0;
	
		orig[9][0].r = 0;	orig[9][0].g = 0;	orig[9][0].b = 0;
		orig[9][1].r = 1;	orig[9][1].g = 1;	orig[9][1].b = 1;
		orig[9][2].r = 0;	orig[9][2].g = 0;	orig[9][2].b = 0;
		orig[9][3].r = 1;	orig[9][3].g = 1;	orig[9][3].b = 1;
		orig[9][4].r = 0;	orig[9][4].g = 0;	orig[9][4].b = 0;
		orig[9][5].r = 1;	orig[9][5].g = 1;	orig[9][5].b = 1;
		orig[9][6].r = 0;	orig[9][6].g = 0;	orig[9][6].b = 0;
		orig[9][7].r = 1;	orig[9][7].g = 1;	orig[9][7].b = 1;
		orig[9][8].r = 1;	orig[9][8].g = 1;	orig[9][8].b = 1;
		orig[9][9].r = 1;	orig[9][9].g = 1;	orig[9][9].b = 1;
		orig[9][10].r = 1;	orig[9][10].g = 1;	orig[9][10].b = 1;
		orig[9][11].r = 1;	orig[9][11].g = 1;	orig[9][11].b = 1;
		orig[9][12].r = 1;	orig[9][12].g = 1;	orig[9][12].b = 1;
		orig[9][13].r = 1;	orig[9][13].g = 1;	orig[9][13].b = 1;
		orig[9][14].r = 0;	orig[9][14].g = 0;	orig[9][14].b = 0;
	
		orig[10][0].r = 0;		orig[10][0].g = 0;		orig[10][0].b = 0;
		orig[10][1].r = 1;		orig[10][1].g = 1;		orig[10][1].b = 1;
		orig[10][2].r = 0;		orig[10][2].g = 0;		orig[10][2].b = 0;
		orig[10][3].r = 0;		orig[10][3].g = 0;		orig[10][3].b = 0;
		orig[10][4].r = 0;		orig[10][4].g = 0;		orig[10][4].b = 0;
		orig[10][5].r = 1;		orig[10][5].g = 1;		orig[10][5].b = 1;
		orig[10][6].r = 1;		orig[10][6].g = 1;		orig[10][6].b = 1;
		orig[10][7].r = 1;		orig[10][7].g = 1;		orig[10][7].b = 1;
		orig[10][8].r = 1;		orig[10][8].g = 1;		orig[10][8].b = 1;
		orig[10][9].r = 1;		orig[10][9].g = 1;		orig[10][9].b = 1;
		orig[10][10].r = 1;	orig[10][10].g = 1;	orig[10][10].b = 1;
		orig[10][11].r = 1;	orig[10][11].g = 1;	orig[10][11].b = 1;
		orig[10][12].r = 1;	orig[10][12].g = 1;	orig[10][12].b = 1;
		orig[10][13].r = 1;	orig[10][13].g = 1;	orig[10][13].b = 1;
		orig[10][14].r = 0;	orig[10][14].g = 0;	orig[10][14].b = 0;
	
		orig[11][0].r = 0;		orig[11][0].g = 0;		orig[11][0].b = 0;
		orig[11][1].r = 1;		orig[11][1].g = 1;		orig[11][1].b = 1;
		orig[11][2].r = 0;		orig[11][2].g = 0;		orig[11][2].b = 0;
		orig[11][3].r = 1;		orig[11][3].g = 1;		orig[11][3].b = 1;
		orig[11][4].r = 0;		orig[11][4].g = 0;		orig[11][4].b = 0;
		orig[11][5].r = 1;		orig[11][5].g = 1;		orig[11][5].b = 1;
		orig[11][6].r = 1;		orig[11][6].g = 1;		orig[11][6].b = 1;
		orig[11][7].r = 0;		orig[11][7].g = 0;		orig[11][7].b = 0;
		orig[11][8].r = 0;		orig[11][8].g = 0;		orig[11][8].b = 0;
		orig[11][9].r = 0;		orig[11][9].g = 0;		orig[11][9].b = 0;
		orig[11][10].r = 1;	orig[11][10].g = 1;	orig[11][10].b = 1;
		orig[11][11].r = 1;	orig[11][11].g = 1;	orig[11][11].b = 1;
		orig[11][12].r = 1;	orig[11][12].g = 1;	orig[11][12].b = 1;
		orig[11][13].r = 1;	orig[11][13].g = 1;	orig[11][13].b = 1;
		orig[11][14].r = 0;	orig[11][14].g = 0;	orig[11][14].b = 0;
	
	
		orig[12][0].r = 0;		orig[12][0].g = 0;		orig[12][0].b = 0;
		orig[12][1].r = 1;		orig[12][1].g = 1;		orig[12][1].b = 1;
		orig[12][2].r = 1;		orig[12][2].g = 1;		orig[12][2].b = 1;
		orig[12][3].r = 0;		orig[12][3].g = 0;		orig[12][3].b = 0;
		orig[12][4].r = 0;		orig[12][4].g = 0;		orig[12][4].b = 0;
		orig[12][5].r = 1;		orig[12][5].g = 1;		orig[12][5].b = 1;
		orig[12][6].r = 1;		orig[12][6].g = 1;		orig[12][6].b = 1;
		orig[12][7].r = 0;		orig[12][7].g = 0;		orig[12][7].b = 0;
		orig[12][8].r = 1;		orig[12][8].g = 0;		orig[12][8].b = 0;
		orig[12][9].r = 0;		orig[12][9].g = 0;		orig[12][9].b = 0;
		orig[12][10].r = 1;	orig[12][10].g = 1;	orig[12][10].b = 1;
		orig[12][11].r = 1;	orig[12][11].g = 1;	orig[12][11].b = 1;
		orig[12][12].r = 1;	orig[12][12].g = 1;	orig[12][12].b = 1;
		orig[12][13].r = 1;	orig[12][13].g = 1;	orig[12][13].b = 1;
		orig[12][14].r = 0;	orig[12][14].g = 0;	orig[12][14].b = 0;
	
		orig[13][0].r = 0;		orig[13][0].g = 0;		orig[13][0].b = 0;
		orig[13][1].r = 0;		orig[13][1].g = 0;		orig[13][1].b = 0;
		orig[13][2].r = 1;		orig[13][2].g = 1;		orig[13][2].b = 1;
		orig[13][3].r = 1;		orig[13][3].g = 1;		orig[13][3].b = 1;
		orig[13][4].r = 0;		orig[13][4].g = 0;		orig[13][4].b = 0;
		orig[13][5].r = 0;		orig[13][5].g = 0;		orig[13][5].b = 0;
		orig[13][6].r = 1;		orig[13][6].g = 1;		orig[13][6].b = 1;
		orig[13][7].r = 0;		orig[13][7].g = 0;		orig[13][7].b = 0;
		orig[13][8].r = 0;		orig[13][8].g = 0;		orig[13][8].b = 0;
		orig[13][9].r = 0;		orig[13][9].g = 0;		orig[13][9].b = 0;
		orig[13][10].r = 1;	orig[13][10].g = 1;	orig[13][10].b = 1;
		orig[13][11].r = 1;	orig[13][11].g = 1;	orig[13][11].b = 1;
		orig[13][12].r = 1;	orig[13][12].g = 1;	orig[13][12].b = 1;
		orig[13][13].r = 0;	orig[13][13].g = 0;	orig[13][13].b = 0;
		orig[13][14].r = 1;	orig[13][14].g = 1;	orig[13][14].b = 1;
	
		orig[14][0].r = 1;		orig[14][0].g = 1;		orig[14][0].b = 1;
		orig[14][1].r = 0;		orig[14][1].g = 0;		orig[14][1].b = 0;
		orig[14][2].r = 0;		orig[14][2].g = 0;		orig[14][2].b = 0;
		orig[14][3].r = 1;		orig[14][3].g = 1;		orig[14][3].b = 1;
		orig[14][4].r = 1;		orig[14][4].g = 1;		orig[14][4].b = 1;
		orig[14][5].r = 1;		orig[14][5].g = 1;		orig[14][5].b = 1;
		orig[14][6].r = 1;		orig[14][6].g = 1;		orig[14][6].b = 1;
		orig[14][7].r = 1;		orig[14][7].g = 1;		orig[14][7].b = 1;
		orig[14][8].r = 1;		orig[14][8].g = 1;		orig[14][8].b = 1;
		orig[14][9].r = 1;		orig[14][9].g = 1;		orig[14][9].b = 1;
		orig[14][10].r = 1;	orig[14][10].g = 1;	orig[14][10].b = 1;
		orig[14][11].r = 1;	orig[14][11].g = 1;	orig[14][11].b = 1;
		orig[14][12].r = 1;	orig[14][12].g = 1;	orig[14][12].b = 1;
		orig[14][13].r = 0;	orig[14][13].g = 0;	orig[14][13].b = 0;
		orig[14][14].r = 1;	orig[14][14].g = 1;	orig[14][14].b = 1;
	
		orig[15][0].r = 1;		orig[15][0].g = 1;		orig[15][0].b = 1;
		orig[15][1].r = 1;		orig[15][1].g = 1;		orig[15][1].b = 1;
		orig[15][2].r = 0;		orig[15][2].g = 0;		orig[15][2].b = 0;
		orig[15][3].r = 1;		orig[15][3].g = 1;		orig[15][3].b = 1;
		orig[15][4].r = 1;		orig[15][4].g = 1;		orig[15][4].b = 1;
		orig[15][5].r = 0;		orig[15][5].g = 0;		orig[15][5].b = 0;
		orig[15][6].r = 0;		orig[15][6].g = 0;		orig[15][6].b = 0;
		orig[15][7].r = 1;		orig[15][7].g = 1;		orig[15][7].b = 1;
		orig[15][8].r = 1;		orig[15][8].g = 1;		orig[15][8].b = 1;
		orig[15][9].r = 1;		orig[15][9].g = 1;		orig[15][9].b = 1;
		orig[15][10].r = 1;	orig[15][10].g = 1;	orig[15][10].b = 1;
		orig[15][11].r = 0;	orig[15][11].g = 0;	orig[15][11].b = 0;
		orig[15][12].r = 0;	orig[15][12].g = 0;	orig[15][12].b = 0;
		orig[15][13].r = 1;	orig[15][13].g = 1;	orig[15][13].b = 1;
		orig[15][14].r = 1;	orig[15][14].g = 1;	orig[15][14].b = 1;
	
		orig[16][0].r = 1;		orig[16][0].g = 1;		orig[16][0].b = 1;
		orig[16][1].r = 1;		orig[16][1].g = 1;		orig[16][1].b = 1;
		orig[16][2].r = 1;		orig[16][2].g = 1;		orig[16][2].b = 1;
		orig[16][3].r = 0;		orig[16][3].g = 0;		orig[16][3].b = 0;
		orig[16][4].r = 0;		orig[16][4].g = 0;		orig[16][4].b = 0;
		orig[16][5].r = 0;		orig[16][5].g = 0;		orig[16][5].b = 0;
		orig[16][6].r = 1;		orig[16][6].g = 1;		orig[16][6].b = 1;
		orig[16][7].r = 0;		orig[16][7].g = 0;		orig[16][7].b = 0;
		orig[16][8].r = 0;		orig[16][8].g = 0;		orig[16][8].b = 0;
		orig[16][9].r = 0;		orig[16][9].g = 0;		orig[16][9].b = 0;
		orig[16][10].r = 0;	orig[16][10].g = 0;	orig[16][10].b = 0;
		orig[16][11].r = 1;	orig[16][11].g = 1;	orig[16][11].b = 1;
		orig[16][12].r = 1;	orig[16][12].g = 1;	orig[16][12].b = 1;
		orig[16][13].r = 1;	orig[16][13].g = 1;	orig[16][13].b = 1;
		orig[16][14].r = 1;	orig[16][14].g = 1;	orig[16][14].b = 1;  */
orig[0][0].r = 0.5;  	orig[0][0].g = 0.5;	orig[0][0].b = 0.5;
orig[0][1].r = 0.5;	orig[0][1].g = 0.5;	orig[0][1].b = 0.5;
orig[0][2].r = 0.5;	orig[0][2].g = 0.5;	orig[0][2].b = 0.5;
orig[0][3].r = 0.5;	orig[0][3].g = 0.5;	orig[0][3].b = 0.5;
orig[0][4].r = 0.5;	orig[0][4].g = 0.5;	orig[0][4].b = 0.5;
orig[0][5].r = 0.5;	orig[0][5].g = 0.5;	orig[0][5].b = 0.5;
orig[0][6].r = 0.5;	orig[0][6].g = 0.5;	orig[0][6].b = 0.5;
orig[0][7].r = 0.5;	orig[0][7].g = 0.5;	orig[0][7].b = 0.5;
orig[0][8].r = 0.5;	orig[0][8].g = 0.5;	orig[0][8].b = 0.5;
orig[0][9].r = 0.5;	orig[0][9].g = 0.5;	orig[0][9].b = 0.5;
orig[0][10].r = 0.5;	orig[0][10].g = 0.5;	orig[0][10].b = 0.5;
orig[0][11].r = 0.5;	orig[0][11].g = 0.5;	orig[0][11].b = 0.5;
orig[0][12].r = 0;	orig[0][12].g = 0;	orig[0][12].b = 0;
orig[0][13].r = 0;	orig[0][13].g = 0.5;	orig[0][13].b = 1;
orig[0][14].r = 0;	orig[0][14].g = 0.5;	orig[0][14].b = 1;
orig[0][15].r = 0;	orig[0][15].g = 0.5;	orig[0][15].b = 1;
orig[0][16].r = 0;	orig[0][16].g = 0.5;	orig[0][16].b = 1;
orig[0][17].r = 0;	orig[0][17].g = 0.5;	orig[0][17].b = 1;
orig[0][18].r = 0.75;	orig[0][18].g = 0.75;	orig[0][18].b = 0.75;
orig[0][19].r = 0.75;	orig[0][19].g = 0.75;	orig[0][19].b = 0.75;
orig[0][20].r = 0.75;	orig[0][20].g = 0.75;	orig[0][20].b = 0.75;
orig[0][21].r = 0;	orig[0][21].g = 0;	orig[0][21].b = 0;
orig[0][22].r = 0;	orig[0][22].g = 0;	orig[0][22].b = 0;
orig[0][23].r = 0;	orig[0][23].g = 0;	orig[0][23].b = 0;
orig[0][24].r = 0;	orig[0][24].g = 0;	orig[0][24].b = 0;
orig[0][25].r = 0;	orig[0][25].g = 0;	orig[0][25].b = 0;
orig[0][26].r = 0;	orig[0][26].g = 0;	orig[0][26].b = 0;
orig[0][27].r = 0;	orig[0][27].g = 0;	orig[0][27].b = 0;
orig[0][28].r = 0;	orig[0][28].g = 0;	orig[0][28].b = 0;
orig[0][29].r = 0;	orig[0][29].g = 0;	orig[0][29].b = 0;

orig[1][0].r = 0.5;  	orig[1][0].g = 0.5;	orig[1][0].b = 0.5;
orig[1][1].r = 0.5;	orig[1][1].g = 0.5;	orig[1][1].b = 0.5;
orig[1][2].r = 0.5;	orig[1][2].g = 0.5;	orig[1][2].b = 0.5;
orig[1][3].r = 0.5;	orig[1][3].g = 0.5;	orig[1][3].b = 0.5;
orig[1][4].r = 0.5;	orig[1][4].g = 0.5;	orig[1][4].b = 0.5;
orig[1][5].r = 0.5;	orig[1][5].g = 0.5;	orig[1][5].b = 0.5;
orig[1][6].r = 0.5;	orig[1][6].g = 0.5;	orig[1][6].b = 0.5;
orig[1][7].r = 0.5;	orig[1][7].g = 0.5;	orig[1][7].b = 0.5;
orig[1][8].r = 0.5;	orig[1][8].g = 0.5;	orig[1][8].b = 0.5;
orig[1][9].r = 0.5;	orig[1][9].g = 0.5;	orig[1][9].b = 0.5;
orig[1][10].r = 0.5;	orig[1][10].g = 0.5;	orig[1][10].b = 0.5;
orig[1][11].r = 0.5;	orig[1][11].g = 0.5;	orig[1][11].b = 0.5;
orig[1][12].r = 0;	orig[1][12].g = 0;	orig[1][12].b = 0;
orig[1][13].r = 0.5;	orig[1][13].g = 0.5;	orig[1][13].b = 0.5;
orig[1][14].r = 0.5;	orig[1][14].g = 0.5;	orig[1][14].b = 0.5;
orig[1][15].r = 0.5;	orig[1][15].g = 0.5;	orig[1][15].b = 0.5;
orig[1][16].r = 0.5;	orig[1][16].g = 0.5;	orig[1][16].b = 0.5;
orig[1][17].r = 0.5;	orig[1][17].g = 0.5;	orig[1][17].b = 0.5;
orig[1][18].r = 0.5;	orig[1][18].g = 0.5;	orig[1][18].b = 0.5;
orig[1][19].r = 0.5;	orig[1][19].g = 0.5;	orig[1][19].b = 0.5;
orig[1][20].r = 0.5;	orig[1][20].g = 0.5;	orig[1][20].b = 0.5;
orig[1][21].r = 0.5;	orig[1][21].g = 0.5;	orig[1][21].b = 0.5;
orig[1][22].r = 0.5;	orig[1][22].g = 0.5;	orig[1][22].b = 0.5;
orig[1][23].r = 0.5;	orig[1][23].g = 0.5;	orig[1][23].b = 0.5;
orig[1][24].r = 0.5;	orig[1][24].g = 0.5;	orig[1][24].b = 0.5;
orig[1][25].r = 0.5;	orig[1][25].g = 0.5;	orig[1][25].b = 0.5;
orig[1][26].r = 0.5;	orig[1][26].g = 0.5;	orig[1][26].b = 0.5;
orig[1][27].r = 0.5;	orig[1][27].g = 0.5;	orig[1][27].b = 0.5;
orig[1][28].r = 0;	orig[1][28].g = 0;	orig[1][28].b = 0;
orig[1][29].r = 0;	orig[1][29].g = 0;	orig[1][29].b = 0;

orig[2][0].r = 0.5;  	orig[2][0].g = 0.5;	orig[2][0].b = 0.5;
orig[2][1].r = 0.5;	orig[2][1].g = 0.5;	orig[2][1].b = 0.5;
orig[2][2].r = 0;	orig[2][2].g = 0;	orig[2][2].b = 1;
orig[2][3].r = 0;	orig[2][3].g = 0;	orig[2][3].b = 1;
orig[2][4].r = 0.5;	orig[2][4].g = 0.5;	orig[2][4].b = 0.5;
orig[2][5].r = 0.5;	orig[2][5].g = 0.5;	orig[2][5].b = 0.5;
orig[2][6].r = 0.5;	orig[2][6].g = 0.5;	orig[2][6].b = 0.5;
orig[2][7].r = 0.5;	orig[2][7].g = 0.5;	orig[2][7].b = 0.5;
orig[2][8].r = 0.5;	orig[2][8].g = 0.5;	orig[2][8].b = 0.5;
orig[2][9].r = 0.5;	orig[2][9].g = 0.5;	orig[2][9].b = 0.5;
orig[2][10].r = 0.5;	orig[2][10].g = 0.5;	orig[2][10].b = 0.5;
orig[2][11].r = 0.5;	orig[2][11].g = 0.5;	orig[2][11].b = 0.5;
orig[2][12].r = 0;	orig[2][12].g = 0;	orig[2][12].b = 0;
orig[2][13].r = 0;	orig[2][13].g = 0.5;	orig[2][13].b = 1;
orig[2][14].r = 0;	orig[2][14].g = 0.5;	orig[2][14].b = 1;
orig[2][15].r = 0;	orig[2][15].g = 0.5;	orig[2][15].b = 1;
orig[2][16].r = 0;	orig[2][16].g = 0.5;	orig[2][16].b = 1;
orig[2][17].r = 0;	orig[2][17].g = 0.5;	orig[2][17].b = 1;
orig[2][18].r = 0;	orig[2][18].g = 0.5;	orig[2][18].b = 1;
orig[2][19].r = 0;	orig[2][19].g = 0.5;	orig[2][19].b = 1;
orig[2][20].r = 0;	orig[2][20].g = 0.5;	orig[2][20].b = 1;
orig[2][21].r = 0;	orig[2][21].g = 0.5;	orig[2][21].b = 1;
orig[2][22].r = 0.75;	orig[2][22].g = 0.75;	orig[2][22].b = 0.75;
orig[2][23].r = 0.75;	orig[2][23].g = 0.75;	orig[2][23].b = 0.75;
orig[2][24].r = 0.75;	orig[2][24].g = 0.75;	orig[2][24].b = 0.75;
orig[2][25].r = 0;	orig[2][25].g = 0;	orig[2][25].b = 0;
orig[2][26].r = 0;	orig[2][26].g = 0;	orig[2][26].b = 0;
orig[2][27].r = 0.5;	orig[2][27].g = 0.5;	orig[2][27].b = 0.5;
orig[2][28].r = 0;	orig[2][28].g = 0;	orig[2][28].b = 0;
orig[2][29].r = 1;	orig[2][29].g = 1;	orig[2][29].b = 1;

orig[3][0].r = 0.5;  	orig[3][0].g = 0.5;	orig[3][0].b = 0.5;
orig[3][1].r = 0;	orig[3][1].g = 0;	orig[3][1].b = 1;
orig[3][2].r = 0;	orig[3][2].g = 0;	orig[3][2].b = 1;
orig[3][3].r = 0;	orig[3][3].g = 0;	orig[3][3].b = 1;
orig[3][4].r = 0;	orig[3][4].g = 0;	orig[3][4].b = 1;
orig[3][5].r = 0.5;	orig[3][5].g = 0.5;	orig[3][5].b = 0.5;
orig[3][6].r = 0.5;	orig[3][6].g = 0.5;	orig[3][6].b = 0.5;
orig[3][7].r = 0.5;	orig[3][7].g = 0.5;	orig[3][7].b = 0.5;
orig[3][8].r = 0.5;	orig[3][8].g = 0.5;	orig[3][8].b = 0.5;
orig[3][9].r = 0.5;	orig[3][9].g = 0.5;	orig[3][9].b = 0.5;
orig[3][10].r = 0.5;	orig[3][10].g = 0.5;	orig[3][10].b = 0.5;
orig[3][11].r = 0.5;	orig[3][11].g = 0.5;	orig[3][11].b = 0.5;
orig[3][12].r = 0;	orig[3][12].g = 0;	orig[3][12].b = 0;
orig[3][13].r = 0;	orig[3][13].g = 0.5;	orig[3][13].b = 1;
orig[3][14].r = 0;	orig[3][14].g = 0.5;	orig[3][14].b = 1;
orig[3][15].r = 0;	orig[3][15].g = 0.5;	orig[3][15].b = 1;
orig[3][16].r = 0;	orig[3][16].g = 0.5;	orig[3][16].b = 1;
orig[3][17].r = 0;	orig[3][17].g = 0.5;	orig[3][17].b = 1;
orig[3][18].r = 0;	orig[3][18].g = 0.5;	orig[3][18].b = 1;
orig[3][19].r = 0;	orig[3][19].g = 0.5;	orig[3][19].b = 1;
orig[3][20].r = 0;	orig[3][20].g = 0.5;	orig[3][20].b = 1;
orig[3][21].r = 0;	orig[3][21].g = 0.5;	orig[3][21].b = 1;
orig[3][22].r = 0;	orig[3][22].g = 0.5;	orig[3][22].b = 1;
orig[3][23].r = 0;	orig[3][23].g = 0.5;	orig[3][23].b = 1;
orig[3][24].r = 0;	orig[3][24].g = 0;	orig[3][24].b = 1;
orig[3][25].r = 0.5;	orig[3][25].g = 0.5;	orig[3][25].b = 0.5;
orig[3][26].r = 0;	orig[3][26].g = 0;	orig[3][26].b = 0;
orig[3][27].r = 0.5;	orig[3][27].g = 0.5;	orig[3][27].b = 0.5;
orig[3][28].r = 0;	orig[3][28].g = 0;	orig[3][28].b = 0;
orig[3][29].r = 0;	orig[3][29].g = 0;	orig[3][29].b = 0;

orig[4][0].r = 0.5;  	orig[4][0].g = 0.5;	orig[4][0].b = 0.5;
orig[4][1].r = 0;	orig[4][1].g = 0;	orig[4][1].b = 1;
orig[4][2].r = 0;	orig[4][2].g = 0;	orig[4][2].b = 1;
orig[4][3].r = 0;	orig[4][3].g = 0;	orig[4][3].b = 1;
orig[4][4].r = 0;	orig[4][4].g = 0;	orig[4][4].b = 1;
orig[4][5].r = 0.5;	orig[4][5].g = 0.5;	orig[4][5].b = 0.5;
orig[4][6].r = 0.5;	orig[4][6].g = 0.5;	orig[4][6].b = 0.5;
orig[4][7].r = 0.5;	orig[4][7].g = 0.5;	orig[4][7].b = 0.5;
orig[4][8].r = 0.5;	orig[4][8].g = 0.5;	orig[4][8].b = 0.5;
orig[4][9].r = 0.5;	orig[4][9].g = 0.5;	orig[4][9].b = 0.5;
orig[4][10].r = 0.5;	orig[4][10].g = 0.5;	orig[4][10].b = 0.5;
orig[4][11].r = 0.5;	orig[4][11].g = 0.5;	orig[4][11].b = 0.5;
orig[4][12].r = 0;	orig[4][12].g = 0;	orig[4][12].b = 0;
orig[4][13].r = 0;	orig[4][13].g = 0.5;	orig[4][13].b = 1;
orig[4][14].r = 0;	orig[4][14].g = 0.5;	orig[4][14].b = 1;
orig[4][15].r = 0;	orig[4][15].g = 0.5;	orig[4][15].b = 1;
orig[4][16].r = 0;	orig[4][16].g = 0.5;	orig[4][16].b = 1;
orig[4][17].r = 0;	orig[4][17].g = 0.5;	orig[4][17].b = 1;
orig[4][18].r = 0;	orig[4][18].g = 0.5;	orig[4][18].b = 1;
orig[4][19].r = 0;	orig[4][19].g = 0.5;	orig[4][19].b = 1;
orig[4][20].r = 0;	orig[4][20].g = 0.5;	orig[4][20].b = 1;
orig[4][21].r = 0;	orig[4][21].g = 0.5;	orig[4][21].b = 1;
orig[4][22].r = 0;	orig[4][22].g = 0.5;	orig[4][22].b = 1;
orig[4][23].r = 0;	orig[4][23].g = 0;	orig[4][23].b = 1;
orig[4][24].r = 0;	orig[4][24].g = 0;	orig[4][24].b = 1;
orig[4][25].r = 0.5;	orig[4][25].g = 0.5;	orig[4][25].b = 0.5;
orig[4][26].r = 0.5;	orig[4][26].g = 0.5;	orig[4][26].b = 0;
orig[4][27].r = 0.5;	orig[4][27].g = 0.5;	orig[4][27].b = 0.5;
orig[4][28].r = 0;	orig[4][28].g = 0;	orig[4][28].b = 0;
orig[4][29].r = 0;	orig[4][29].g = 0;	orig[4][29].b = 0;

orig[5][0].r = 0.5;  	orig[5][0].g = 0.5;	orig[5][0].b = 0.5;
orig[5][1].r = 0;	orig[5][1].g = 0;	orig[5][1].b = 1;
orig[5][2].r = 0;	orig[5][2].g = 0;	orig[5][2].b = 1;
orig[5][3].r = 0;	orig[5][3].g = 0;	orig[5][3].b = 1;
orig[5][4].r = 0;	orig[5][4].g = 0;	orig[5][4].b = 1;
orig[5][5].r = 0.5;	orig[5][5].g = 0.5;	orig[5][5].b = 0.5;
orig[5][6].r = 0.5;	orig[5][6].g = 0.5;	orig[5][6].b = 0.5;
orig[5][7].r = 0.5;	orig[5][7].g = 0.5;	orig[5][7].b = 0.5;
orig[5][8].r = 0.5;	orig[5][8].g = 0.5;	orig[5][8].b = 0.5;
orig[5][9].r = 0.5;	orig[5][9].g = 0.5;	orig[5][9].b = 0.5;
orig[5][10].r = 0.5;	orig[5][10].g = 0.5;	orig[5][10].b = 0.5;
orig[5][11].r = 0.5;	orig[5][11].g = 0.5;	orig[5][11].b = 0.5;
orig[5][12].r = 0;	orig[5][12].g = 0;	orig[5][12].b = 0;
orig[5][13].r = 0;	orig[5][13].g = 0.5;	orig[5][13].b = 1;
orig[5][14].r = 0;	orig[5][14].g = 0.5;	orig[5][14].b = 1;
orig[5][15].r = 0;	orig[5][15].g = 0.5;	orig[5][15].b = 1;
orig[5][16].r = 0;	orig[5][16].g = 0.5;	orig[5][16].b = 1;
orig[5][17].r = 0;	orig[5][17].g = 0.5;	orig[5][17].b = 1;
orig[5][18].r = 0;	orig[5][18].g = 0.5;	orig[5][18].b = 1;
orig[5][19].r = 0;	orig[5][19].g = 0.5;	orig[5][19].b = 1;
orig[5][20].r = 0;	orig[5][20].g = 0.5;	orig[5][20].b = 1;
orig[5][21].r = 0;	orig[5][21].g = 0.5;	orig[5][21].b = 1;
orig[5][22].r = 0;	orig[5][22].g = 0.5;	orig[5][22].b = 1;
orig[5][23].r = 0;	orig[5][23].g = 0.5;	orig[5][23].b = 1;
orig[5][24].r = 0;	orig[5][24].g = 0;	orig[5][24].b = 1;
orig[5][25].r = 0.5;	orig[5][25].g = 0.5;	orig[5][25].b = 0.5;
orig[5][26].r = 0;	orig[5][26].g = 0;	orig[5][26].b = 0;
orig[5][27].r = 0;	orig[5][27].g = 0;	orig[5][27].b = 0;
orig[5][28].r = 0;	orig[5][28].g = 0;	orig[5][28].b = 0;
orig[5][29].r = 0;	orig[5][29].g = 0;	orig[5][29].b = 0;

orig[6][0].r = 0.5;  	orig[6][0].g = 0.5;	orig[6][0].b = 0.5;
orig[6][1].r = 0;	orig[6][1].g = 0;	orig[6][1].b = 1;
orig[6][2].r = 0;	orig[6][2].g = 0;	orig[6][2].b = 1;
orig[6][3].r = 0;	orig[6][3].g = 0;	orig[6][3].b = 1;
orig[6][4].r = 0;	orig[6][4].g = 0;	orig[6][4].b = 1;
orig[6][5].r = 0.5;	orig[6][5].g = 0.5;	orig[6][5].b = 0.5;
orig[6][6].r = 0.5;	orig[6][6].g = 0.5;	orig[6][6].b = 0.5;
orig[6][7].r = 0.5;	orig[6][7].g = 0.5;	orig[6][7].b = 0.5;
orig[6][8].r = 0.5;	orig[6][8].g = 0.5;	orig[6][8].b = 0.5;
orig[6][9].r = 0.5;	orig[6][9].g = 0.5;	orig[6][9].b = 0.5;
orig[6][10].r = 0.5;	orig[6][10].g = 0.5;	orig[6][10].b = 0.5;
orig[6][11].r = 0.5;	orig[6][11].g = 0.5;	orig[6][11].b = 0.5;
orig[6][12].r = 0;	orig[6][12].g = 0;	orig[6][12].b = 0;
orig[6][13].r = 0;	orig[6][13].g = 0.5;	orig[6][13].b = 1;
orig[6][14].r = 0;	orig[6][14].g = 0.5;	orig[6][14].b = 1;
orig[6][15].r = 0;	orig[6][15].g = 0.5;	orig[6][15].b = 1;
orig[6][16].r = 0;	orig[6][16].g = 0.5;	orig[6][16].b = 1;
orig[6][17].r = 0;	orig[6][17].g = 0.5;	orig[6][17].b = 1;
orig[6][18].r = 0;	orig[6][18].g = 0.5;	orig[6][18].b = 1;
orig[6][19].r = 0;	orig[6][19].g = 0.5;	orig[6][19].b = 1;
orig[6][20].r = 0;	orig[6][20].g = 0.5;	orig[6][20].b = 1;
orig[6][21].r = 0;	orig[6][21].g = 0.5;	orig[6][21].b = 1;
orig[6][22].r = 0.75;	orig[6][22].g = 0.75;	orig[6][22].b = 0.75;
orig[6][23].r = 0.75;	orig[6][23].g = 0.75;	orig[6][23].b = 0.75;
orig[6][24].r = 0.75;	orig[6][24].g = 0.75;	orig[6][24].b = 0.75;
orig[6][25].r = 0;	orig[6][25].g = 0;	orig[6][25].b = 0;
orig[6][26].r = 0;	orig[6][26].g = 0;	orig[6][26].b = 0;
orig[6][27].r = 0;	orig[6][27].g = 0;	orig[6][27].b = 0;
orig[6][28].r = 0;	orig[6][28].g = 0;	orig[6][28].b = 0;
orig[6][29].r = 1;	orig[6][29].g = 1;	orig[6][29].b = 1;

orig[7][0].r = 0.5;  	orig[7][0].g = 0.5;	orig[7][0].b = 0.5;
orig[7][1].r = 0.5;	orig[7][1].g = 0.5;	orig[7][1].b = 0.5;
orig[7][2].r = 0;	orig[7][2].g = 0;	orig[7][2].b = 1;
orig[7][3].r = 0;	orig[7][3].g = 0;	orig[7][3].b = 1;
orig[7][4].r = 0;	orig[7][4].g = 0;	orig[7][4].b = 1;
orig[7][5].r = 0.5;	orig[7][5].g = 0.5;	orig[7][5].b = 0.5;
orig[7][6].r = 0.5;	orig[7][6].g = 0.5;	orig[7][6].b = 0.5;
orig[7][7].r = 0;	orig[7][7].g = 0;	orig[7][7].b = 1;
orig[7][8].r = 0;	orig[7][8].g = 0;	orig[7][8].b = 1;
orig[7][9].r = 0.5;	orig[7][9].g = 0.5;	orig[7][9].b = 0.5;
orig[7][10].r = 0.5;	orig[7][10].g = 0.5;	orig[7][10].b = 0.5;
orig[7][11].r = 0.5;	orig[7][11].g = 0.5;	orig[7][11].b = 0.5;
orig[7][12].r = 0;	orig[7][12].g = 0;	orig[7][12].b = 0;
orig[7][13].r = 0;	orig[7][13].g = 0.5;	orig[7][13].b = 1;
orig[7][14].r = 0;	orig[7][14].g = 0.5;	orig[7][14].b = 1;
orig[7][15].r = 0;	orig[7][15].g = 0.5;	orig[7][15].b = 1;
orig[7][16].r = 0;	orig[7][16].g = 0.5;	orig[7][16].b = 1;
orig[7][17].r = 0;	orig[7][17].g = 0.5;	orig[7][17].b = 1;
orig[7][18].r = 0;	orig[7][18].g = 0.5;	orig[7][18].b = 1;
orig[7][19].r = 0;	orig[7][19].g = 0.5;	orig[7][19].b = 1;
orig[7][20].r = 0.75;	orig[7][20].g = 0.75;	orig[7][20].b = 0.75;
orig[7][21].r = 0.75;	orig[7][21].g = 0.75;	orig[7][21].b = 0.75;
orig[7][22].r = 0.75;	orig[7][22].g = 0.75;	orig[7][22].b = 0.75;
orig[7][23].r = 0;	orig[7][23].g = 0;	orig[7][23].b = 0;
orig[7][24].r = 0;	orig[7][24].g = 0;	orig[7][24].b = 0;
orig[7][25].r = 0;	orig[7][25].g = 0;	orig[7][25].b = 0;
orig[7][26].r = 0;	orig[7][26].g = 0;	orig[7][26].b = 0;
orig[7][27].r = 0;	orig[7][27].g = 0;	orig[7][27].b = 0;
orig[7][28].r = 0;	orig[7][28].g = 0;	orig[7][28].b = 0;
orig[7][29].r = 0;	orig[7][29].g = 0;	orig[7][29].b = 0;

orig[8][0].r = 0.5;  	orig[8][0].g = 0.5;	orig[8][0].b = 0.5;
orig[8][1].r = 0.5;	orig[8][1].g = 0.5;	orig[8][1].b = 0.5;
orig[8][2].r = 0.5;	orig[8][2].g = 0.5;	orig[8][2].b = 0.5;
orig[8][3].r = 0;	orig[8][3].g = 0;	orig[8][3].b = 1;
orig[8][4].r = 0.5;	orig[8][4].g = 0.5;	orig[8][4].b = 0.5;
orig[8][5].r = 0.5;	orig[8][5].g = 0.5;	orig[8][5].b = 0.5;
orig[8][6].r = 0.5;	orig[8][6].g = 0.5;	orig[8][6].b = 0.5;
orig[8][7].r = 0;	orig[8][7].g = 0;	orig[8][7].b = 1;
orig[8][8].r = 0;	orig[8][8].g = 0;	orig[8][8].b = 1;
orig[8][9].r = 0.5;	orig[8][9].g = 0.5;	orig[8][9].b = 0.5;
orig[8][10].r = 0.5;	orig[8][10].g = 0.5;	orig[8][10].b = 0.5;
orig[8][11].r = 0.5;	orig[8][11].g = 0.5;	orig[8][11].b = 0.5;
orig[8][12].r = 0;	orig[8][12].g = 0;	orig[8][12].b = 0;
orig[8][13].r = 0;	orig[8][13].g = 0.5;	orig[8][13].b = 1;
orig[8][14].r = 0;	orig[8][14].g = 0.5;	orig[8][14].b = 1;
orig[8][15].r = 0;	orig[8][15].g = 0.5;	orig[8][15].b = 1;
orig[8][16].r = 0;	orig[8][16].g = 0.5;	orig[8][16].b = 1;
orig[8][17].r = 0;	orig[8][17].g = 0.5;	orig[8][17].b = 1;
orig[8][18].r = 0.75;	orig[8][18].g = 0.75;	orig[8][18].b = 0.75;
orig[8][19].r = 0.75;	orig[8][19].g = 0.75;	orig[8][19].b = 0.75;
orig[8][20].r = 0.75;	orig[8][20].g = 0.75;	orig[8][20].b = 0.75;
orig[8][21].r = 0;	orig[8][21].g = 0;	orig[8][21].b = 0;
orig[8][22].r = 0;	orig[8][22].g = 0;	orig[8][22].b = 0;
orig[8][23].r = 0;	orig[8][23].g = 0;	orig[8][23].b = 0;
orig[8][24].r = 0;	orig[8][24].g = 0;	orig[8][24].b = 0;
orig[8][25].r = 0;	orig[8][25].g = 0;	orig[8][25].b = 0;
orig[8][26].r = 0;	orig[8][26].g = 0;	orig[8][26].b = 0;
orig[8][27].r = 0;	orig[8][27].g = 0;	orig[8][27].b = 0;
orig[8][28].r = 0;	orig[8][28].g = 0;	orig[8][28].b = 0;
orig[8][29].r = 0;	orig[8][29].g = 0;	orig[8][29].b = 0;

orig[9][0].r = 0.5;  	orig[9][0].g = 0.5;	orig[9][0].b = 0.5;
orig[9][1].r = 0.5;	orig[9][1].g = 0.5;	orig[9][1].b = 0.5;
orig[9][2].r = 0.5;	orig[9][2].g = 0.5;	orig[9][2].b = 0.5;
orig[9][3].r = 0.5;	orig[9][3].g = 0.5;	orig[9][3].b = 0.5;
orig[9][4].r = 0.5;	orig[9][4].g = 0.5;	orig[9][4].b = 0.5;
orig[9][5].r = 0.5;	orig[9][5].g = 0.5;	orig[9][5].b = 0.5;
orig[9][6].r = 0;	orig[9][6].g = 0;	orig[9][6].b = 1;
orig[9][7].r = 0;	orig[9][7].g = 0;	orig[9][7].b = 1;
orig[9][8].r = 0;	orig[9][8].g = 0;	orig[9][8].b = 1;
orig[9][9].r = 0;	orig[9][9].g = 0;	orig[9][9].b = 1;
orig[9][10].r = 0;	orig[9][10].g = 0;	orig[9][10].b = 1;
orig[9][11].r = 0.5;	orig[9][11].g = 0.5;	orig[9][11].b = 0.5;
orig[9][12].r = 0;	orig[9][12].g = 0;	orig[9][12].b = 0;
orig[9][13].r = 0;	orig[9][13].g = 0.5;	orig[9][13].b = 1;
orig[9][14].r = 0;	orig[9][14].g = 0.5;	orig[9][14].b = 1;
orig[9][15].r = 0;	orig[9][15].g = 0.5;	orig[9][15].b = 1;
orig[9][16].r = 0.75;	orig[9][16].g = 0.75;	orig[9][16].b = 0.75;
orig[9][17].r = 0.75;	orig[9][17].g = 0.75;	orig[9][17].b = 0.75;
orig[9][18].r = 0.75;	orig[9][18].g = 0.75;	orig[9][18].b = 0.75;
orig[9][19].r = 0;	orig[9][19].g = 0;	orig[9][19].b = 0;
orig[9][20].r = 0;	orig[9][20].g = 0;	orig[9][20].b = 0;
orig[9][21].r = 0;	orig[9][21].g = 0;	orig[9][21].b = 0;
orig[9][22].r = 0;	orig[9][22].g = 0;	orig[9][22].b = 0;
orig[9][23].r = 0;	orig[9][23].g = 0;	orig[9][23].b = 0;
orig[9][24].r = 0;	orig[9][24].g = 0;	orig[9][24].b = 0;
orig[9][25].r = 0;	orig[9][25].g = 0;	orig[9][25].b = 0;
orig[9][26].r = 0;	orig[9][26].g = 0;	orig[9][26].b = 0;
orig[9][27].r = 0;	orig[9][27].g = 0;	orig[9][27].b = 0;
orig[9][28].r = 1;	orig[9][28].g = 1;	orig[9][28].b = 1;
orig[9][29].r = 0;	orig[9][29].g = 0;	orig[9][29].b = 0;

orig[10][0].r = 0.5;  	orig[10][0].g = 0.5;	orig[10][0].b = 0.5;
orig[10][1].r = 0.5;	orig[10][1].g = 0.5;	orig[10][1].b = 0.5;
orig[10][2].r = 0.5;	orig[10][2].g = 0.5;	orig[10][2].b = 0.5;
orig[10][3].r = 0.5;	orig[10][3].g = 0.5;	orig[10][3].b = 0.5;
orig[10][4].r = 0.5;	orig[10][4].g = 0.5;	orig[10][4].b = 0.5;
orig[10][5].r = 0.5;	orig[10][5].g = 0.5;	orig[10][5].b = 0.5;
orig[10][6].r = 0;	orig[10][6].g = 0;	orig[10][6].b = 1;
orig[10][7].r = 0;	orig[10][7].g = 0;	orig[10][7].b = 1;
orig[10][8].r = 0;	orig[10][8].g = 0;	orig[10][8].b = 1;
orig[10][9].r = 0;	orig[10][9].g = 0;	orig[10][9].b = 1;
orig[10][10].r = 0;	orig[10][10].g = 0;	orig[10][10].b = 1;
orig[10][11].r = 0.5;	orig[10][11].g = 0.5;	orig[10][11].b = 0.5;
orig[10][12].r = 0;	orig[10][12].g = 0;	orig[10][12].b = 0;
orig[10][13].r = 0;	orig[10][13].g = 0.5;	orig[10][13].b = 1;
orig[10][14].r = 0.75;	orig[10][14].g = 0.75;	orig[10][14].b = 0.75;
orig[10][15].r = 0.75;	orig[10][15].g = 0.75;	orig[10][15].b = 0.75;
orig[10][16].r = 0.75;	orig[10][16].g = 0.75;	orig[10][16].b = 0.75;
orig[10][17].r = 0;	orig[10][17].g = 0;	orig[10][17].b = 0;
orig[10][18].r = 0;	orig[10][18].g = 0;	orig[10][18].b = 0;
orig[10][19].r = 0;	orig[10][19].g = 0;	orig[10][19].b = 0;
orig[10][20].r = 0;	orig[10][20].g = 0;	orig[10][20].b = 0;
orig[10][21].r = 0;	orig[10][21].g = 0;	orig[10][21].b = 0;
orig[10][22].r = 0;	orig[10][22].g = 0;	orig[10][22].b = 0;
orig[10][23].r = 0;	orig[10][23].g = 0;	orig[10][23].b = 0;
orig[10][24].r = 1;	orig[10][24].g = 1;	orig[10][24].b = 1;
orig[10][25].r = 0;	orig[10][25].g = 0;	orig[10][25].b = 0;
orig[10][26].r = 0;	orig[10][26].g = 0;	orig[10][26].b = 0;
orig[10][27].r = 0;	orig[10][27].g = 0;	orig[10][27].b = 0;
orig[10][28].r = 0;	orig[10][28].g = 0;	orig[10][28].b = 0;
orig[10][29].r = 0;	orig[10][29].g = 0;	orig[10][29].b = 0;

orig[11][0].r = 0.5;  	orig[11][0].g = 0.5;	orig[11][0].b = 0.5;
orig[11][1].r = 0.5;	orig[11][1].g = 0.5;	orig[11][1].b = 0.5;
orig[11][2].r = 0.5;	orig[11][2].g = 0.5;	orig[11][2].b = 0.5;
orig[11][3].r = 0.5;	orig[11][3].g = 0.5;	orig[11][3].b = 0.5;
orig[11][4].r = 0.5;	orig[11][4].g = 0.5;	orig[11][4].b = 0.5;
orig[11][5].r = 0.5;	orig[11][5].g = 0.5;	orig[11][5].b = 0.5;
orig[11][6].r = 0;	orig[11][6].g = 0;	orig[11][6].b = 1;
orig[11][7].r = 0;	orig[11][7].g = 0;	orig[11][7].b = 1;
orig[11][8].r = 0;	orig[11][8].g = 0;	orig[11][8].b = 1;
orig[11][9].r = 0;	orig[11][9].g = 0;	orig[11][9].b = 1;
orig[11][10].r = 0;	orig[11][10].g = 0;	orig[11][10].b = 1;
orig[11][11].r = 0.5;	orig[11][11].g = 0.5;	orig[11][11].b = 0.5;
orig[11][12].r = 0;	orig[11][12].g = 0;	orig[11][12].b = 0;
orig[11][13].r = 0.75;	orig[11][13].g = 0.75;	orig[11][13].b = 0.75;
orig[11][14].r = 0;	orig[11][14].g = 0;	orig[11][14].b = 0;
orig[11][15].r = 0;	orig[11][15].g = 0;	orig[11][15].b = 0;
orig[11][16].r = 0;	orig[11][16].g = 0;	orig[11][16].b = 0;
orig[11][17].r = 0;	orig[11][17].g = 0;	orig[11][17].b = 0;
orig[11][18].r = 0;	orig[11][18].g = 0;	orig[11][18].b = 0;
orig[11][19].r = 0;	orig[11][19].g = 0;	orig[11][19].b = 0;
orig[11][20].r = 0;	orig[11][20].g = 0;	orig[11][20].b = 0;
orig[11][21].r = 0;	orig[11][21].g = 0;	orig[11][21].b = 0;
orig[11][22].r = 0;	orig[11][22].g = 0;	orig[11][22].b = 0;
orig[11][23].r = 0;	orig[11][23].g = 0;	orig[11][23].b = 0;
orig[11][24].r = 0;	orig[11][24].g = 0;	orig[11][24].b = 0;
orig[11][25].r = 0;	orig[11][25].g = 0;	orig[11][25].b = 0;
orig[11][26].r = 0;	orig[11][26].g = 0;	orig[11][26].b = 0;
orig[11][27].r = 0;	orig[11][27].g = 0;	orig[11][27].b = 0;
orig[11][28].r = 0;	orig[11][28].g = 0;	orig[11][28].b = 0;
orig[11][29].r = 0;	orig[11][29].g = 0;	orig[11][29].b = 0;

orig[12][0].r = 0.5;  	orig[12][0].g = 0.5;	orig[12][0].b = 0.5;
orig[12][1].r = 0.5;	orig[12][1].g = 0.5;	orig[12][1].b = 0.5;
orig[12][2].r = 0.5;	orig[12][2].g = 0.5;	orig[12][2].b = 0.5;
orig[12][3].r = 0.5;	orig[12][3].g = 0.5;	orig[12][3].b = 0.5;
orig[12][4].r = 0.5;	orig[12][4].g = 0.5;	orig[12][4].b = 0.5;
orig[12][5].r = 0.5;	orig[12][5].g = 0.5;	orig[12][5].b = 0.5;
orig[12][6].r = 0;	orig[12][6].g = 0;	orig[12][6].b = 1;
orig[12][7].r = 0;	orig[12][7].g = 0;	orig[12][7].b = 1;
orig[12][8].r = 0;	orig[12][8].g = 0;	orig[12][8].b = 1;
orig[12][9].r = 0;	orig[12][9].g = 0;	orig[12][9].b = 1;
orig[12][10].r = 0;	orig[12][10].g = 0;	orig[12][10].b = 1;
orig[12][11].r = 0.5;	orig[12][11].g = 0.5;	orig[12][11].b = 0.5;
orig[12][12].r = 0;	orig[12][12].g = 0;	orig[12][12].b = 0;
orig[12][13].r = 0.5;	orig[12][13].g = 0.5;	orig[12][13].b = 0.5;
orig[12][14].r = 0;	orig[12][14].g = 0;	orig[12][14].b = 0;
orig[12][15].r = 0;	orig[12][15].g = 0;	orig[12][15].b = 0;
orig[12][16].r = 0;	orig[12][16].g = 0;	orig[12][16].b = 0;
orig[12][17].r = 0;	orig[12][17].g = 0;	orig[12][17].b = 0;
orig[12][18].r = 0;	orig[12][18].g = 0;	orig[12][18].b = 0;
orig[12][19].r = 0;	orig[12][19].g = 0;	orig[12][19].b = 0;
orig[12][20].r = 0;	orig[12][20].g = 0;	orig[12][20].b = 0;
orig[12][21].r = 0;	orig[12][21].g = 0;	orig[12][21].b = 0;
orig[12][22].r = 0;	orig[12][22].g = 0;	orig[12][22].b = 0;
orig[12][23].r = 0;	orig[12][23].g = 0;	orig[12][23].b = 0;
orig[12][24].r = 0;	orig[12][24].g = 0;	orig[12][24].b = 0;
orig[12][25].r = 0;	orig[12][25].g = 0;	orig[12][25].b = 0;
orig[12][26].r = 0;	orig[12][26].g = 0;	orig[12][26].b = 0;
orig[12][27].r = 1;	orig[12][27].g = 1;	orig[12][27].b = 1;
orig[12][28].r = 1;	orig[12][28].g = 1;	orig[12][28].b = 1;
orig[12][29].r = 0;	orig[12][29].g = 0;	orig[12][29].b = 0;

orig[13][0].r = 0.5;  	orig[13][0].g = 0.5;	orig[13][0].b = 0.5;
orig[13][1].r = 0.5;	orig[13][1].g = 0.5;	orig[13][1].b = 0.5;
orig[13][2].r = 0.5;	orig[13][2].g = 0.5;	orig[13][2].b = 0.5;
orig[13][3].r = 0.5;	orig[13][3].g = 0.5;	orig[13][3].b = 0.5;
orig[13][4].r = 0.5;	orig[13][4].g = 0.5;	orig[13][4].b = 0.5;
orig[13][5].r = 0.5;	orig[13][5].g = 0.5;	orig[13][5].b = 0.5;
orig[13][6].r = 0;	orig[13][6].g = 0;	orig[13][6].b = 1;
orig[13][7].r = 0;	orig[13][7].g = 0;	orig[13][7].b = 1;
orig[13][8].r = 0;	orig[13][8].g = 0;	orig[13][8].b = 1;
orig[13][9].r = 0;	orig[13][9].g = 0;	orig[13][9].b = 1;
orig[13][10].r = 0;	orig[13][10].g = 0;	orig[13][10].b = 1;
orig[13][11].r = 0.5;	orig[13][11].g = 0.5;	orig[13][11].b = 0.5;
orig[13][12].r = 0;	orig[13][12].g = 0;	orig[13][12].b = 0;
orig[13][13].r = 0.25;	orig[13][13].g = 0.25;	orig[13][13].b = 0.25;
orig[13][14].r = 0;	orig[13][14].g = 0;	orig[13][14].b = 0;
orig[13][15].r = 0;	orig[13][15].g = 0;	orig[13][15].b = 0;
orig[13][16].r = 0;	orig[13][16].g = 0;	orig[13][16].b = 0;
orig[13][17].r = 0;	orig[13][17].g = 0;	orig[13][17].b = 0;
orig[13][18].r = 0;	orig[13][18].g = 0;	orig[13][18].b = 0;
orig[13][19].r = 0;	orig[13][19].g = 0;	orig[13][19].b = 0;
orig[13][20].r = 0;	orig[13][20].g = 0;	orig[13][20].b = 0;
orig[13][21].r = 1;	orig[13][21].g = 1;	orig[13][21].b = 1;
orig[13][22].r = 0;	orig[13][22].g = 0;	orig[13][22].b = 0;
orig[13][23].r = 0;	orig[13][23].g = 0;	orig[13][23].b = 0;
orig[13][24].r = 0;	orig[13][24].g = 0;	orig[13][24].b = 0;
orig[13][25].r = 0;	orig[13][25].g = 0;	orig[13][25].b = 0;
orig[13][26].r = 1;	orig[13][26].g = 1;	orig[13][26].b = 1;
orig[13][27].r = 1;	orig[13][27].g = 1;	orig[13][27].b = 1;
orig[13][28].r = 1;	orig[13][28].g = 1;	orig[13][28].b = 1;
orig[13][29].r = 1;	orig[13][29].g = 1;	orig[13][29].b = 1;

orig[14][0].r = 0.5;  	orig[14][0].g = 0.5;	orig[14][0].b = 0.5;
orig[14][1].r = 0.5;	orig[14][1].g = 0.5;	orig[14][1].b = 0.5;
orig[14][2].r = 0.5;	orig[14][2].g = 0.5;	orig[14][2].b = 0.5;
orig[14][3].r = 0.5;	orig[14][3].g = 0.5;	orig[14][3].b = 0.5;
orig[14][4].r = 0.5;	orig[14][4].g = 0.5;	orig[14][4].b = 0.5;
orig[14][5].r = 0.5;	orig[14][5].g = 0.5;	orig[14][5].b = 0.5;
orig[14][6].r = 0;	orig[14][6].g = 0;	orig[14][6].b = 1;
orig[14][7].r = 0;	orig[14][7].g = 0;	orig[14][7].b = 1;
orig[14][8].r = 0;	orig[14][8].g = 0;	orig[14][8].b = 1;
orig[14][9].r = 0;	orig[14][9].g = 0;	orig[14][9].b = 1;
orig[14][10].r = 0;	orig[14][10].g = 0;	orig[14][10].b = 1;
orig[14][11].r = 0.5;	orig[14][11].g = 0.5;	orig[14][11].b = 0.5;
orig[14][12].r = 0;	orig[14][12].g = 0;	orig[14][12].b = 0;
orig[14][13].r = 0.5;	orig[14][13].g = 0.5;	orig[14][13].b = 0.5;
orig[14][14].r = 0;	orig[14][14].g = 0;	orig[14][14].b = 0;
orig[14][15].r = 0;	orig[14][15].g = 0;	orig[14][15].b = 0;
orig[14][16].r = 0;	orig[14][16].g = 0;	orig[14][16].b = 0;
orig[14][17].r = 0;	orig[14][17].g = 0;	orig[14][17].b = 0;
orig[14][18].r = 0;	orig[14][18].g = 0;	orig[14][18].b = 0;
orig[14][19].r = 0;	orig[14][19].g = 0;	orig[14][19].b = 0;
orig[14][20].r = 0;	orig[14][20].g = 0;	orig[14][20].b = 0;
orig[14][21].r = 0;	orig[14][21].g = 0;	orig[14][21].b = 0;
orig[14][22].r = 0;	orig[14][22].g = 0;	orig[14][22].b = 0;
orig[14][23].r = 0;	orig[14][23].g = 0;	orig[14][23].b = 0;
orig[14][24].r = 0;	orig[14][24].g = 0;	orig[14][24].b = 0;
orig[14][25].r = 0;	orig[14][25].g = 0;	orig[14][25].b = 0;
orig[14][26].r = 1;	orig[14][26].g = 1;	orig[14][26].b = 1;
orig[14][27].r = 1;	orig[14][27].g = 1;	orig[14][27].b = 1;
orig[14][28].r = 1;	orig[14][28].g = 1;	orig[14][28].b = 1;
orig[14][29].r = 1;	orig[14][29].g = 1;	orig[14][29].b = 1;

orig[15][0].r = 0.5;  	orig[15][0].g = 0.5;	orig[15][0].b = 0.5;
orig[15][1].r = 0.5;	orig[15][1].g = 0.5;	orig[15][1].b = 0.5;
orig[15][2].r = 0.5;	orig[15][2].g = 0.5;	orig[15][2].b = 0.5;
orig[15][3].r = 0.5;	orig[15][3].g = 0.5;	orig[15][3].b = 0.5;
orig[15][4].r = 0.5;	orig[15][4].g = 0.5;	orig[15][4].b = 0.5;
orig[15][5].r = 0.5;	orig[15][5].g = 0.5;	orig[15][5].b = 0.5;
orig[15][6].r = 0;	orig[15][6].g = 0;	orig[15][6].b = 1;
orig[15][7].r = 0;	orig[15][7].g = 0;	orig[15][7].b = 1;
orig[15][8].r = 0;	orig[15][8].g = 0;	orig[15][8].b = 1;
orig[15][9].r = 0;	orig[15][9].g = 0;	orig[15][9].b = 1;
orig[15][10].r = 0;	orig[15][10].g = 0;	orig[15][10].b = 1;
orig[15][11].r = 0.5;	orig[15][11].g = 0.5;	orig[15][11].b = 0.5;
orig[15][12].r = 0;	orig[15][12].g = 0;	orig[15][12].b = 0;
orig[15][13].r = 0.75;	orig[15][13].g = 0.75;	orig[15][13].b = 0.75;
orig[15][14].r = 0;	orig[15][14].g = 0;	orig[15][14].b = 0;
orig[15][15].r = 0;	orig[15][15].g = 0;	orig[15][15].b = 0;
orig[15][16].r = 0;	orig[15][16].g = 0;	orig[15][16].b = 0;
orig[15][17].r = 0;	orig[15][17].g = 0;	orig[15][17].b = 0;
orig[15][18].r = 0;	orig[15][18].g = 0;	orig[15][18].b = 0;
orig[15][19].r = 0;	orig[15][19].g = 0;	orig[15][19].b = 0;
orig[15][20].r = 0;	orig[15][20].g = 0;	orig[15][20].b = 0;
orig[15][21].r = 0;	orig[15][21].g = 0;	orig[15][21].b = 0;
orig[15][22].r = 0;	orig[15][22].g = 0;	orig[15][22].b = 0;
orig[15][23].r = 0;	orig[15][23].g = 0;	orig[15][23].b = 0;
orig[15][24].r = 0;	orig[15][24].g = 0;	orig[15][24].b = 0;
orig[15][25].r = 0;	orig[15][25].g = 0;	orig[15][25].b = 0;
orig[15][26].r = 0;	orig[15][26].g = 0;	orig[15][26].b = 0;
orig[15][27].r = 1;	orig[15][27].g = 1;	orig[15][27].b = 1;
orig[15][28].r = 1;	orig[15][28].g = 1;	orig[15][28].b = 1;
orig[15][29].r = 0;	orig[15][29].g = 0;	orig[15][29].b = 0;

orig[16][0].r = 0.5;  	orig[16][0].g = 0.5;	orig[16][0].b = 0.5;
orig[16][1].r = 0.5;	orig[16][1].g = 0.5;	orig[16][1].b = 0.5;
orig[16][2].r = 0.5;	orig[16][2].g = 0.5;	orig[16][2].b = 0.5;
orig[16][3].r = 0.5;	orig[16][3].g = 0.5;	orig[16][3].b = 0.5;
orig[16][4].r = 0.5;	orig[16][4].g = 0.5;	orig[16][4].b = 0.5;
orig[16][5].r = 0.5;	orig[16][5].g = 0.5;	orig[16][5].b = 0.5;
orig[16][6].r = 0;	orig[16][6].g = 0;	orig[16][6].b = 1;
orig[16][7].r = 0;	orig[16][7].g = 0;	orig[16][7].b = 1;
orig[16][8].r = 0;	orig[16][8].g = 0;	orig[16][8].b = 1;
orig[16][9].r = 0;	orig[16][9].g = 0;	orig[16][9].b = 1;
orig[16][10].r = 0;	orig[16][10].g = 0;	orig[16][10].b = 1;
orig[16][11].r = 0.5;	orig[16][11].g = 0.5;	orig[16][11].b = 0.5;
orig[16][12].r = 0;	orig[16][12].g = 0;	orig[16][12].b = 0;
orig[16][13].r = 0;	orig[16][13].g = 0.5;	orig[16][13].b = 1;
orig[16][14].r = 0.75;	orig[16][14].g = 0.75;	orig[16][14].b = 0.75;
orig[16][15].r = 0.75;	orig[16][15].g = 0.75;	orig[16][15].b = 0.75;
orig[16][16].r = 0.75;	orig[16][16].g = 0.75;	orig[16][16].b = 0.75;
orig[16][17].r = 0;	orig[16][17].g = 0;	orig[16][17].b = 0;
orig[16][18].r = 0;	orig[16][18].g = 0;	orig[16][18].b = 0;
orig[16][19].r = 0;	orig[16][19].g = 0;	orig[16][19].b = 0;
orig[16][20].r = 0;	orig[16][20].g = 0;	orig[16][20].b = 0;
orig[16][21].r = 0;	orig[16][21].g = 0;	orig[16][21].b = 0;
orig[16][22].r = 0;	orig[16][22].g = 0;	orig[16][22].b = 0;
orig[16][23].r = 1;	orig[16][23].g = 1;	orig[16][23].b = 1;
orig[16][24].r = 0;	orig[16][24].g = 0;	orig[16][24].b = 0;
orig[16][25].r = 0;	orig[16][25].g = 0;	orig[16][25].b = 0;
orig[16][26].r = 0;	orig[16][26].g = 0;	orig[16][26].b = 0;
orig[16][27].r = 0;	orig[16][27].g = 0;	orig[16][27].b = 0;
orig[16][28].r = 0;	orig[16][28].g = 0;	orig[16][28].b = 0;
orig[16][29].r = 0;	orig[16][29].g = 0;	orig[16][29].b = 0;

orig[17][0].r = 0.5;  	orig[17][0].g = 0.5;	orig[17][0].b = 0.5;
orig[17][1].r = 0.5;	orig[17][1].g = 0.5;	orig[17][1].b = 0.5;
orig[17][2].r = 0.5;	orig[17][2].g = 0.5;	orig[17][2].b = 0.5;
orig[17][3].r = 0.5;	orig[17][3].g = 0.5;	orig[17][3].b = 0.5;
orig[17][4].r = 0.5;	orig[17][4].g = 0.5;	orig[17][4].b = 0.5;
orig[17][5].r = 0.5;	orig[17][5].g = 0.5;	orig[17][5].b = 0.5;
orig[17][6].r = 0;	orig[17][6].g = 0;	orig[17][6].b = 1;
orig[17][7].r = 0;	orig[17][7].g = 0;	orig[17][7].b = 1;
orig[17][8].r = 0;	orig[17][8].g = 0;	orig[17][8].b = 1;
orig[17][9].r = 0;	orig[17][9].g = 0;	orig[17][9].b = 1;
orig[17][10].r = 0;	orig[17][10].g = 0;	orig[17][10].b = 1;
orig[17][11].r = 0.5;	orig[17][11].g = 0.5;	orig[17][11].b = 0.5;
orig[17][12].r = 0;	orig[17][12].g = 0;	orig[17][12].b = 0;
orig[17][13].r = 0;	orig[17][13].g = 0.5;	orig[17][13].b = 1;
orig[17][14].r = 0;	orig[17][14].g = 0.5;	orig[17][14].b = 1;
orig[17][15].r = 0;	orig[17][15].g = 0.5;	orig[17][15].b = 1;
orig[17][16].r = 0.75;	orig[17][16].g = 0.75;	orig[17][16].b = 0.75;
orig[17][17].r = 0.75;	orig[17][17].g = 0.75;	orig[17][17].b = 0.75;
orig[17][18].r = 0.75;	orig[17][18].g = 0.75;	orig[17][18].b = 0.75;
orig[17][19].r = 0;	orig[17][19].g = 0;	orig[17][19].b = 0;
orig[17][20].r = 0;	orig[17][20].g = 0;	orig[17][20].b = 0;
orig[17][21].r = 0;	orig[17][21].g = 0;	orig[17][21].b = 0;
orig[17][22].r = 0;	orig[17][22].g = 0;	orig[17][22].b = 0;
orig[17][23].r = 0;	orig[17][23].g = 0;	orig[17][23].b = 0;
orig[17][24].r = 0;	orig[17][24].g = 0;	orig[17][24].b = 0;
orig[17][25].r = 0;	orig[17][25].g = 0;	orig[17][25].b = 0;
orig[17][26].r = 0;	orig[17][26].g = 0;	orig[17][26].b = 0;
orig[17][27].r = 0;	orig[17][27].g = 0;	orig[17][27].b = 0;
orig[17][28].r = 0;	orig[17][28].g = 0;	orig[17][28].b = 0;
orig[17][29].r = 0;	orig[17][29].g = 0;	orig[17][29].b = 0;

orig[18][0].r = 0.5;  	orig[18][0].g = 0.5;	orig[18][0].b = 0.5;
orig[18][1].r = 0.5;	orig[18][1].g = 0.5;	orig[18][1].b = 0.5;
orig[18][2].r = 0.5;	orig[18][2].g = 0.5;	orig[18][2].b = 0.5;
orig[18][3].r = 0.5;	orig[18][3].g = 0.5;	orig[18][3].b = 0.5;
orig[18][4].r = 0.5;	orig[18][4].g = 0.5;	orig[18][4].b = 0.5;
orig[18][5].r = 0.5;	orig[18][5].g = 0.5;	orig[18][5].b = 0.5;
orig[18][6].r = 0;	orig[18][6].g = 0;	orig[18][6].b = 1;
orig[18][7].r = 0;	orig[18][7].g = 0;	orig[18][7].b = 1;
orig[18][8].r = 0;	orig[18][8].g = 0;	orig[18][8].b = 1;
orig[18][9].r = 0;	orig[18][9].g = 0;	orig[18][9].b = 1;
orig[18][10].r = 0;	orig[18][10].g = 0;	orig[18][10].b = 1;
orig[18][11].r = 0.5;	orig[18][11].g = 0.5;	orig[18][11].b = 0.5;
orig[18][12].r = 0;	orig[18][12].g = 0;	orig[18][12].b = 0;
orig[18][13].r = 0;	orig[18][13].g = 0.5;	orig[18][13].b = 1;
orig[18][14].r = 0;	orig[18][14].g = 0.5;	orig[18][14].b = 1;
orig[18][15].r = 0;	orig[18][15].g = 0.5;	orig[18][15].b = 1;
orig[18][16].r = 0;	orig[18][16].g = 0.5;	orig[18][16].b = 1;
orig[18][17].r = 0;	orig[18][17].g = 0.5;	orig[18][17].b = 1;
orig[18][18].r = 0.75;	orig[18][18].g = 0.75;	orig[18][18].b = 0.75;
orig[18][19].r = 0.75;	orig[18][19].g = 0.75;	orig[18][19].b = 0.75;
orig[18][20].r = 0.75;	orig[18][20].g = 0.75;	orig[18][20].b = 0.75;
orig[18][21].r = 0;	orig[18][21].g = 0;	orig[18][21].b = 0;
orig[18][22].r = 0;	orig[18][22].g = 0;	orig[18][22].b = 0;
orig[18][23].r = 0;	orig[18][23].g = 0;	orig[18][23].b = 0;
orig[18][24].r = 0;	orig[18][24].g = 0;	orig[18][24].b = 0;
orig[18][25].r = 0;	orig[18][25].g = 0;	orig[18][25].b = 0;
orig[18][26].r = 1;	orig[18][26].g = 1;	orig[18][26].b = 1;
orig[18][27].r = 0;	orig[18][27].g = 0;	orig[18][27].b = 0;
orig[18][28].r = 0;	orig[18][28].g = 0;	orig[18][28].b = 0;
orig[18][29].r = 0;	orig[18][29].g = 0;	orig[18][29].b = 0;

orig[19][0].r = 0.5;  	orig[19][0].g = 0.5;	orig[19][0].b = 0.5;
orig[19][1].r = 0.5;	orig[19][1].g = 0.5;	orig[19][1].b = 0.5;
orig[19][2].r = 0.5;	orig[19][2].g = 0.5;	orig[19][2].b = 0.5;
orig[19][3].r = 0;	orig[19][3].g = 0;	orig[19][3].b = 1;
orig[19][4].r = 0.5;	orig[19][4].g = 0.5;	orig[19][4].b = 0.5;
orig[19][5].r = 0.5;	orig[19][5].g = 0.5;	orig[19][5].b = 0.5;
orig[19][6].r = 0.5;	orig[19][6].g = 0.5;	orig[19][6].b = 0.5;
orig[19][7].r = 0;	orig[19][7].g = 0;	orig[19][7].b = 1;
orig[19][8].r = 0;	orig[19][8].g = 0;	orig[19][8].b = 1;
orig[19][9].r = 0;	orig[19][9].g = 0;	orig[19][9].b = 1;
orig[19][10].r = 0;	orig[19][10].g = 0;	orig[19][10].b = 1;
orig[19][11].r = 0.5;	orig[19][11].g = 0.5;	orig[19][11].b = 0.5;
orig[19][12].r = 0;	orig[19][12].g = 0;	orig[19][12].b = 0;
orig[19][13].r = 0;	orig[19][13].g = 0.5;	orig[19][13].b = 1;
orig[19][14].r = 0;	orig[19][14].g = 0.5;	orig[19][14].b = 1;
orig[19][15].r = 0;	orig[19][15].g = 0.5;	orig[19][15].b = 1;
orig[19][16].r = 0;	orig[19][16].g = 0.5;	orig[19][16].b = 1;
orig[19][17].r = 0;	orig[19][17].g = 0.5;	orig[19][17].b = 1;
orig[19][18].r = 0;	orig[19][18].g = 0.5;	orig[19][18].b = 1;
orig[19][19].r = 0;	orig[19][19].g = 0.5;	orig[19][19].b = 1;
orig[19][20].r = 0.75;	orig[19][20].g = 0.75;	orig[19][20].b = 0.75;
orig[19][21].r = 0.75;	orig[19][21].g = 0.75;	orig[19][21].b = 0.75;
orig[19][22].r = 0.75;	orig[19][22].g = 0.75;	orig[19][22].b = 0.75;
orig[19][23].r = 0;	orig[19][23].g = 0;	orig[19][23].b = 0;
orig[19][24].r = 0;	orig[19][24].g = 0;	orig[19][24].b = 0;
orig[19][25].r = 0;	orig[19][25].g = 0;	orig[19][25].b = 0;
orig[19][26].r = 0;	orig[19][26].g = 0;	orig[19][26].b = 0;
orig[19][27].r = 0;	orig[19][27].g = 0;	orig[19][27].b = 0;
orig[19][28].r = 0;	orig[19][28].g = 0;	orig[19][28].b = 0;
orig[19][29].r = 0;	orig[19][29].g = 0;	orig[19][29].b = 0;

orig[20][0].r = 0.5;  	orig[20][0].g = 0.5;	orig[20][0].b = 0.5;
orig[20][1].r = 0.5;	orig[20][1].g = 0.5;	orig[20][1].b = 0.5;
orig[20][2].r = 0;	orig[20][2].g = 0;	orig[20][2].b = 1;
orig[20][3].r = 0;	orig[20][3].g = 0;	orig[20][3].b = 1;
orig[20][4].r = 0;	orig[20][4].g = 0;	orig[20][4].b = 1;
orig[20][5].r = 0.5;	orig[20][5].g = 0.5;	orig[20][5].b = 0.5;
orig[20][6].r = 0.5;	orig[20][6].g = 0.5;	orig[20][6].b = 0.5;
orig[20][7].r = 0.5;	orig[20][7].g = 0.5;	orig[20][7].b = 0.5;
orig[20][8].r = 0;	orig[20][8].g = 0;	orig[20][8].b = 1;
orig[20][9].r = 0;	orig[20][9].g = 0;	orig[20][9].b = 1;
orig[20][10].r = 0.5;	orig[20][10].g = 0.5;	orig[20][10].b = 0.5;
orig[20][11].r = 0.5;	orig[20][11].g = 0.5;	orig[20][11].b = 0.5;
orig[20][12].r = 0;	orig[20][12].g = 0;	orig[20][12].b = 0;
orig[20][13].r = 0;	orig[20][13].g = 0.5;	orig[20][13].b = 1;
orig[20][14].r = 0;	orig[20][14].g = 0.5;	orig[20][14].b = 1;
orig[20][15].r = 0;	orig[20][15].g = 0.5;	orig[20][15].b = 1;
orig[20][16].r = 0;	orig[20][16].g = 0.5;	orig[20][16].b = 1;
orig[20][17].r = 0;	orig[20][17].g = 0.5;	orig[20][17].b = 1;
orig[20][18].r = 0;	orig[20][18].g = 0.5;	orig[20][18].b = 1;
orig[20][19].r = 0;	orig[20][19].g = 0.5;	orig[20][19].b = 1;
orig[20][20].r = 0;	orig[20][20].g = 0.5;	orig[20][20].b = 1;
orig[20][21].r = 0;	orig[20][21].g = 0.5;	orig[20][21].b = 1;
orig[20][22].r = 0.75;	orig[20][22].g = 0.75;	orig[20][22].b = 0.75;
orig[20][23].r = 0.75;	orig[20][23].g = 0.75;	orig[20][23].b = 0.75;
orig[20][24].r = 0.75;	orig[20][24].g = 0.75;	orig[20][24].b = 0.75;
orig[20][25].r = 0;	orig[20][25].g = 0;	orig[20][25].b = 0;
orig[20][26].r = 0;	orig[20][26].g = 0;	orig[20][26].b = 0;
orig[20][27].r = 0;	orig[20][27].g = 0;	orig[20][27].b = 0;
orig[20][28].r = 0;	orig[20][28].g = 0;	orig[20][28].b = 0;
orig[20][29].r = 0;	orig[20][29].g = 0;	orig[20][29].b = 0;

orig[21][0].r = 0.5;  	orig[21][0].g = 0.5;	orig[21][0].b = 0.5;
orig[21][1].r = 0.5;	orig[21][1].g = 0.5;	orig[21][1].b = 0.5;
orig[21][2].r = 0;	orig[21][2].g = 0;	orig[21][2].b = 1;
orig[21][3].r = 0;	orig[21][3].g = 0;	orig[21][3].b = 1;
orig[21][4].r = 0;	orig[21][4].g = 0;	orig[21][4].b = 1;
orig[21][5].r = 0.5;	orig[21][5].g = 0.5;	orig[21][5].b = 0.5;
orig[21][6].r = 0.5;	orig[21][6].g = 0.5;	orig[21][6].b = 0.5;
orig[21][7].r = 0.5;	orig[21][7].g = 0.5;	orig[21][7].b = 0.5;
orig[21][8].r = 0.5;	orig[21][8].g = 0.5;	orig[21][8].b = 0.5;
orig[21][9].r = 0.5;	orig[21][9].g = 0.5;	orig[21][9].b = 0.5;
orig[21][10].r = 0.5;	orig[21][10].g = 0.5;	orig[21][10].b = 0.5;
orig[21][11].r = 0.5;	orig[21][11].g = 0.5;	orig[21][11].b = 0.5;
orig[21][12].r = 0;	orig[21][12].g = 0;	orig[21][12].b = 0;
orig[21][13].r = 0;	orig[21][13].g = 0.5;	orig[21][13].b = 1;
orig[21][14].r = 0;	orig[21][14].g = 0.5;	orig[21][14].b = 1;
orig[21][15].r = 0;	orig[21][15].g = 0.5;	orig[21][15].b = 1;
orig[21][16].r = 0;	orig[21][16].g = 0.5;	orig[21][16].b = 1;
orig[21][17].r = 0;	orig[21][17].g = 0.5;	orig[21][17].b = 1;
orig[21][18].r = 0;	orig[21][18].g = 0.5;	orig[21][18].b = 1;
orig[21][19].r = 0;	orig[21][19].g = 0.5;	orig[21][19].b = 1;
orig[21][20].r = 0;	orig[21][20].g = 0.5;	orig[21][20].b = 1;
orig[21][21].r = 0;	orig[21][21].g = 0.5;	orig[21][21].b = 1;
orig[21][22].r = 0;	orig[21][22].g = 0.5;	orig[21][22].b = 1;
orig[21][23].r = 0;	orig[21][23].g = 0.5;	orig[21][23].b = 1;
orig[21][24].r = 0;	orig[21][24].g = 0;	orig[21][24].b = 1;
orig[21][25].r = 0.5;	orig[21][25].g = 0.5;	orig[21][25].b = 0.5;
orig[21][26].r = 0;	orig[21][26].g = 0;	orig[21][26].b = 0;
orig[21][27].r = 0;	orig[21][27].g = 0;	orig[21][27].b = 0;
orig[21][28].r = 1;	orig[21][28].g = 1;	orig[21][28].b = 1;
orig[21][29].r = 0;	orig[21][29].g = 0;	orig[21][29].b = 0;

orig[22][0].r = 0.5;  	orig[22][0].g = 0.5;	orig[22][0].b = 0.5;
orig[22][1].r = 0.5;	orig[22][1].g = 0.5;	orig[22][1].b = 0.5;
orig[22][2].r = 0;	orig[22][2].g = 0;	orig[22][2].b = 1;
orig[22][3].r = 0;	orig[22][3].g = 0;	orig[22][3].b = 1;
orig[22][4].r = 0;	orig[22][4].g = 0;	orig[22][4].b = 1;
orig[22][5].r = 0;	orig[22][5].g = 0;	orig[22][5].b = 1;
orig[22][6].r = 0.5;	orig[22][6].g = 0.5;	orig[22][6].b = 0.5;
orig[22][7].r = 0.5;	orig[22][7].g = 0.5;	orig[22][7].b = 0.5;
orig[22][8].r = 0.5;	orig[22][8].g = 0.5;	orig[22][8].b = 0.5;
orig[22][9].r = 0.5;	orig[22][9].g = 0.5;	orig[22][9].b = 0.5;
orig[22][10].r = 0.5;	orig[22][10].g = 0.5;	orig[22][10].b = 0.5;
orig[22][11].r = 0.5;	orig[22][11].g = 0.5;	orig[22][11].b = 0.5;
orig[22][12].r = 0;	orig[22][12].g = 0;	orig[22][12].b = 0;
orig[22][13].r = 0;	orig[22][13].g = 0.5;	orig[22][13].b = 1;
orig[22][14].r = 0;	orig[22][14].g = 0.5;	orig[22][14].b = 1;
orig[22][15].r = 0;	orig[22][15].g = 0.5;	orig[22][15].b = 1;
orig[22][16].r = 0;	orig[22][16].g = 0.5;	orig[22][16].b = 1;
orig[22][17].r = 0;	orig[22][17].g = 0.5;	orig[22][17].b = 1;
orig[22][18].r = 0;	orig[22][18].g = 0.5;	orig[22][18].b = 1;
orig[22][19].r = 0;	orig[22][19].g = 0.5;	orig[22][19].b = 1;
orig[22][20].r = 0;	orig[22][20].g = 0.5;	orig[22][20].b = 1;
orig[22][21].r = 0;	orig[22][21].g = 0.5;	orig[22][21].b = 1;
orig[22][22].r = 0;	orig[22][22].g = 0.5;	orig[22][22].b = 1;
orig[22][23].r = 0;	orig[22][23].g = 0;	orig[22][23].b = 1;
orig[22][24].r = 0;	orig[22][24].g = 0;	orig[22][24].b = 1;
orig[22][25].r = 0.5;	orig[22][25].g = 0.5;	orig[22][25].b = 0.5;
orig[22][26].r = 0.5;	orig[22][26].g = 0.5;	orig[22][26].b = 0;
orig[22][27].r = 0.5;	orig[22][27].g = 0.5;	orig[22][27].b = 0.5;
orig[22][28].r = 0;	orig[22][28].g = 0;	orig[22][28].b = 0;
orig[22][29].r = 0;	orig[22][29].g = 0;	orig[22][29].b = 0;

orig[23][0].r = 0.5;  	orig[23][0].g = 0.5;	orig[23][0].b = 0.5;
orig[23][1].r = 0.5;	orig[23][1].g = 0.5;	orig[23][1].b = 0.5;
orig[23][2].r = 0;	orig[23][2].g = 0;	orig[23][2].b = 1;
orig[23][3].r = 0;	orig[23][3].g = 0;	orig[23][3].b = 1;
orig[23][4].r = 0;	orig[23][4].g = 0;	orig[23][4].b = 1;
orig[23][5].r = 0;	orig[23][5].g = 0;	orig[23][5].b = 1;
orig[23][6].r = 0.5;	orig[23][6].g = 0.5;	orig[23][6].b = 0.5;
orig[23][7].r = 0.5;	orig[23][7].g = 0.5;	orig[23][7].b = 0.5;
orig[23][8].r = 0.5;	orig[23][8].g = 0.5;	orig[23][8].b = 0.5;
orig[23][9].r = 0.5;	orig[23][9].g = 0.5;	orig[23][9].b = 0.5;
orig[23][10].r = 0.5;	orig[23][10].g = 0.5;	orig[23][10].b = 0.5;
orig[23][11].r = 0.5;	orig[23][11].g = 0.5;	orig[23][11].b = 0.5;
orig[23][12].r = 0;	orig[23][12].g = 0;	orig[23][12].b = 0;
orig[23][13].r = 0;	orig[23][13].g = 0.5;	orig[23][13].b = 1;
orig[23][14].r = 0;	orig[23][14].g = 0.5;	orig[23][14].b = 1;
orig[23][15].r = 0;	orig[23][15].g = 0.5;	orig[23][15].b = 1;
orig[23][16].r = 0;	orig[23][16].g = 0.5;	orig[23][16].b = 1;
orig[23][17].r = 0;	orig[23][17].g = 0.5;	orig[23][17].b = 1;
orig[23][18].r = 0;	orig[23][18].g = 0.5;	orig[23][18].b = 1;
orig[23][19].r = 0;	orig[23][19].g = 0.5;	orig[23][19].b = 1;
orig[23][20].r = 0;	orig[23][20].g = 0.5;	orig[23][20].b = 1;
orig[23][21].r = 0;	orig[23][21].g = 0.5;	orig[23][21].b = 1;
orig[23][22].r = 0;	orig[23][22].g = 0.5;	orig[23][22].b = 1;
orig[23][23].r = 0;	orig[23][23].g = 0.5;	orig[23][23].b = 1;
orig[23][24].r = 0;	orig[23][24].g = 0;	orig[23][24].b = 1;
orig[23][25].r = 0.5;	orig[23][25].g = 0.5;	orig[23][25].b = 0.5;
orig[23][26].r = 0;	orig[23][26].g = 0;	orig[23][26].b = 0;
orig[23][27].r = 0.5;	orig[23][27].g = 0.5;	orig[23][27].b = 0.5;
orig[23][28].r = 0;	orig[23][28].g = 0;	orig[23][28].b = 0;
orig[23][29].r = 0;	orig[23][29].g = 0;	orig[23][29].b = 0;

orig[24][0].r = 0.5;  	orig[24][0].g = 0.5;	orig[24][0].b = 0.5;
orig[24][1].r = 0.5;	orig[24][1].g = 0.5;	orig[24][1].b = 0.5;
orig[24][2].r = 0.5;	orig[24][2].g = 0.5;	orig[24][2].b = 0.5;
orig[24][3].r = 0;	orig[24][3].g = 0;	orig[24][3].b = 1;
orig[24][4].r = 0;	orig[24][4].g = 0;	orig[24][4].b = 1;
orig[24][5].r = 0.5;	orig[24][5].g = 0.5;	orig[24][5].b = 0.5;
orig[24][6].r = 0.5;	orig[24][6].g = 0.5;	orig[24][6].b = 0.5;
orig[24][7].r = 0.5;	orig[24][7].g = 0.5;	orig[24][7].b = 0.5;
orig[24][8].r = 0.5;	orig[24][8].g = 0.5;	orig[24][8].b = 0.5;
orig[24][9].r = 0.5;	orig[24][9].g = 0.5;	orig[24][9].b = 0.5;
orig[24][10].r = 0.5;	orig[24][10].g = 0.5;	orig[24][10].b = 0.5;
orig[24][11].r = 0.5;	orig[24][11].g = 0.5;	orig[24][11].b = 0.5;
orig[24][12].r = 0;	orig[24][12].g = 0;	orig[24][12].b = 0;
orig[24][13].r = 0;	orig[24][13].g = 0.5;	orig[24][13].b = 1;
orig[24][14].r = 0;	orig[24][14].g = 0.5;	orig[24][14].b = 1;
orig[24][15].r = 0;	orig[24][15].g = 0.5;	orig[24][15].b = 1;
orig[24][16].r = 0;	orig[24][16].g = 0.5;	orig[24][16].b = 1;
orig[24][17].r = 0;	orig[24][17].g = 0.5;	orig[24][17].b = 1;
orig[24][18].r = 0;	orig[24][18].g = 0.5;	orig[24][18].b = 1;
orig[24][19].r = 0;	orig[24][19].g = 0.5;	orig[24][19].b = 1;
orig[24][20].r = 0;	orig[24][20].g = 0.5;	orig[24][20].b = 1;
orig[24][21].r = 0;	orig[24][21].g = 0.5;	orig[24][21].b = 1;
orig[24][22].r = 0.75;	orig[24][22].g = 0.75;	orig[24][22].b = 0.75;
orig[24][23].r = 0.75;	orig[24][23].g = 0.75;	orig[24][23].b = 0.75;
orig[24][24].r = 0.75;	orig[24][24].g = 0.75;	orig[24][24].b = 0.75;
orig[24][25].r = 0;	orig[24][25].g = 0;	orig[24][25].b = 0;
orig[24][26].r = 0;	orig[24][26].g = 0;	orig[24][26].b = 0;
orig[24][27].r = 0.5;	orig[24][27].g = 0.5;	orig[24][27].b = 0.5;
orig[24][28].r = 0;	orig[24][28].g = 0;	orig[24][28].b = 0;
orig[24][29].r = 0;	orig[24][29].g = 0;	orig[24][29].b = 0;

orig[25][0].r = 0.5;  	orig[25][0].g = 0.5;	orig[25][0].b = 0.5;
orig[25][1].r = 0.5;	orig[25][1].g = 0.5;	orig[25][1].b = 0.5;
orig[25][2].r = 0.5;	orig[25][2].g = 0.5;	orig[25][2].b = 0.5;
orig[25][3].r = 0.5;	orig[25][3].g = 0.5;	orig[25][3].b = 0.5;
orig[25][4].r = 0;	orig[25][4].g = 0;	orig[25][4].b = 1;
orig[25][5].r = 0.5;	orig[25][5].g = 0.5;	orig[25][5].b = 0.5;
orig[25][6].r = 0.5;	orig[25][6].g = 0.5;	orig[25][6].b = 0.5;
orig[25][7].r = 0.5;	orig[25][7].g = 0.5;	orig[25][7].b = 0.5;
orig[25][8].r = 0.5;	orig[25][8].g = 0.5;	orig[25][8].b = 0.5;
orig[25][9].r = 0.5;	orig[25][9].g = 0.5;	orig[25][9].b = 0.5;
orig[25][10].r = 0.5;	orig[25][10].g = 0.5;	orig[25][10].b = 0.5;
orig[25][11].r = 0.5;	orig[25][11].g = 0.5;	orig[25][11].b = 0.5;
orig[25][12].r = 0;	orig[25][12].g = 0;	orig[25][12].b = 0;
orig[25][13].r = 0.5;	orig[25][13].g = 0.5;	orig[25][13].b = 0.5;
orig[25][14].r = 0.5;	orig[25][14].g = 0.5;	orig[25][14].b = 0.5;
orig[25][15].r = 0.5;	orig[25][15].g = 0.5;	orig[25][15].b = 0.5;
orig[25][16].r = 0.5;	orig[25][16].g = 0.5;	orig[25][16].b = 0.5;
orig[25][17].r = 0.5;	orig[25][17].g = 0.5;	orig[25][17].b = 0.5;
orig[25][18].r = 0.5;	orig[25][18].g = 0.5;	orig[25][18].b = 0.5;
orig[25][19].r = 0.5;	orig[25][19].g = 0.5;	orig[25][19].b = 0.5;
orig[25][20].r = 0.5;	orig[25][20].g = 0.5;	orig[25][20].b = 0.5;
orig[25][21].r = 0.5;	orig[25][21].g = 0.5;	orig[25][21].b = 0.5;
orig[25][22].r = 0.5;	orig[25][22].g = 0.5;	orig[25][22].b = 0.5;
orig[25][23].r = 0.5;	orig[25][23].g = 0.5;	orig[25][23].b = 0.5;
orig[25][24].r = 0.5;	orig[25][24].g = 0.5;	orig[25][24].b = 0.5;
orig[25][25].r = 0.5;	orig[25][25].g = 0.5;	orig[25][25].b = 0.5;
orig[25][26].r = 0.5;	orig[25][26].g = 0.5;	orig[25][26].b = 0.5;
orig[25][27].r = 0.5;	orig[25][27].g = 0.5;	orig[25][27].b = 0.5;
orig[25][28].r = 0;	orig[25][28].g = 0;	orig[25][28].b = 0;
orig[25][29].r = 1;	orig[25][29].g = 1;	orig[25][29].b = 1;

orig[26][0].r = 0.5;  	orig[26][0].g = 0.5;	orig[26][0].b = 0.5;
orig[26][1].r = 0.5;	orig[26][1].g = 0.5;	orig[26][1].b = 0.5;
orig[26][2].r = 0.5;	orig[26][2].g = 0.5;	orig[26][2].b = 0.5;
orig[26][3].r = 0.5;	orig[26][3].g = 0.5;	orig[26][3].b = 0.5;
orig[26][4].r = 0.5;	orig[26][4].g = 0.5;	orig[26][4].b = 0.5;
orig[26][5].r = 0.5;	orig[26][5].g = 0.5;	orig[26][5].b = 0.5;
orig[26][6].r = 0.5;	orig[26][6].g = 0.5;	orig[26][6].b = 0.5;
orig[26][7].r = 0.5;	orig[26][7].g = 0.5;	orig[26][7].b = 0.5;
orig[26][8].r = 0.5;	orig[26][8].g = 0.5;	orig[26][8].b = 0.5;
orig[26][9].r = 0.5;	orig[26][9].g = 0.5;	orig[26][9].b = 0.5;
orig[26][10].r = 0.5;	orig[26][10].g = 0.5;	orig[26][10].b = 0.5;
orig[26][11].r = 0.5;	orig[26][11].g = 0.5;	orig[26][11].b = 0.5;
orig[26][12].r = 0;	orig[26][12].g = 0;	orig[26][12].b = 0;
orig[26][13].r = 0;	orig[26][13].g = 0.5;	orig[26][13].b = 1;
orig[26][14].r = 0;	orig[26][14].g = 0.5;	orig[26][14].b = 1;
orig[26][15].r = 0;	orig[26][15].g = 0.5;	orig[26][15].b = 1;
orig[26][16].r = 0;	orig[26][16].g = 0.5;	orig[26][16].b = 1;
orig[26][17].r = 0;	orig[26][17].g = 0.5;	orig[26][17].b = 1;
orig[26][18].r = 0.75;	orig[26][18].g = 0.75;	orig[26][18].b = 0.75;
orig[26][19].r = 0.75;	orig[26][19].g = 0.75;	orig[26][19].b = 0.75;
orig[26][20].r = 0.75;	orig[26][20].g = 0.75;	orig[26][20].b = 0.75;
orig[26][21].r = 0;	orig[26][21].g = 0;	orig[26][21].b = 0;
orig[26][22].r = 0;	orig[26][22].g = 0;	orig[26][22].b = 0;
orig[26][23].r = 0;	orig[26][23].g = 0;	orig[26][23].b = 0;
orig[26][24].r = 0;	orig[26][24].g = 0;	orig[26][24].b = 0;
orig[26][25].r = 0;	orig[26][25].g = 0;	orig[26][25].b = 0;
orig[26][26].r = 0;	orig[26][26].g = 0;	orig[26][26].b = 0;
orig[26][27].r = 0;	orig[26][27].g = 0;	orig[26][27].b = 0;
orig[26][28].r = 0;	orig[26][28].g = 0;	orig[26][28].b = 0;
orig[26][29].r = 0;	orig[26][29].g = 0;	orig[26][29].b = 0;
	for (int nom = 0; nom < m1; nom++)
	{
		bys[nom] = -1 + nom * shag2;
	}
	
	int op1 = 0;
	int op2 = n[0];
	for (int i = 0; i < colvo; i++)
	{
		bx[i] = new double[n[i]];
	}
	int no = 0;
	for (int i = 0; i < colvo; i++)
	{
		for (int j = start[i][0]; j <= kon[i][0]; j++)
		{
			bx[i][no] = bxs[j];
			no++;
		}
		no = 0;
	}
	for (int i = 0; i < colvo; i++)
	{

		by[i] = new double[m[i]];
		for (int j = 0; j < m[i]; j++)
		{
			by[i][j] = -1 + shag2 * op1;
			op1++;
			if (op1 == m1)
				op1 = 0;
		}
		otvr = new double**[colvo];
		otvg = new double**[colvo];
		otvb = new double**[colvo];
		for (int i = 0; i < colvo; i++)
		{
			int razm = (n[i] - 1)*(m[i] - 1) + 1;
			otvr[i] = new double*[razm];
			otvg[i] = new double*[razm];
			otvb[i] = new double*[razm];

			for (int j = 0; j < razm; j++)
			{
				otvr[i][j] = new double[razm];
				otvg[i][j] = new double[razm];
				otvb[i][j] = new double[razm];
				for (int k = 0; k < razm; k++)
				{

					otvr[i][j][k] = 0;
					otvg[i][j][k] = 0;
					otvb[i][j][k] = 0;
				}
			}
		}
	}
	for (int i = 0; i < colvo; i++) {
		interpolate(start[i][0], start[i][1], n[i], m[i], otvr[i], otvg[i], otvb[i], i);
	}
	/*for (int i = 0; i < colvo; i++)
	{
		int razm = (n[i] - 1)*(m[i] - 1) + 1;
		for (int j = 0; j < razm; j++)
		{
			for (int k = 0; k < razm; k++)
				cout << "(" << otvr[i][j][k] << ")*x^" << razm - j - 1 << "*y^" << razm - k - 1 << "+";
			cout << endl;
		}
		cout << endl;
		cout << endl;
	}
	cout << endl;
	cout << endl;
	cout << endl;
	for (int i = 0; i < colvo; i++)
	{
		int razm = (n[i] - 1)*(m[i] - 1) + 1;
		for (int j = 0; j < razm; j++)
		{
			for (int k = 0; k < razm; k++)
				cout << "(" << otvg[i][j][k] << ")*x^" << razm - j - 1 << "*y^" << razm - k - 1 << "+";
			cout << endl;
		}
		cout << endl;
		cout << endl;
	}
	cout << endl;
	cout << endl;
	cout << endl;
	for (int i = 0; i < colvo; i++)
	{
		int razm = (n[i] - 1)*(m[i] - 1) + 1;
		for (int j = 0; j < razm; j++)
		{
			for (int k = 0; k < razm; k++)
				cout << "(" << otvb[i][j][k] << ")*x^" << razm - j - 1 << "*y^" << razm - k - 1 << "+";
			cout << endl;
		}
		cout << endl;
		cout << endl;
	}
	cout << endl;
	cout << endl;
	cout << endl;*/
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 0);
	GLint win1 = glutCreateWindow("image1");
	glutReshapeFunc(changeViewPort1);
	glEnable(GL_DEPTH_TEST);
	glutDisplayFunc(render);
	GLenum err = glewInit();
	if (GLEW_OK != err) {
		fprintf(stderr, "GLEW error");
		return 1;
	}
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(550, 0);
	GLint win2 = glutCreateWindow("image2");
	glutReshapeFunc(changeViewPort2);
	glEnable(GL_DEPTH_TEST);
	glutDisplayFunc(render1);
	GLenum err2 = glewInit();
	if (GLEW_OK != err2) {
		fprintf(stderr, "GLEW error");
		return 1;
	}
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(1100, 0);
	GLint win3 = glutCreateWindow("graphic1");
	glutReshapeFunc(changeViewPort3);
	glEnable(GL_DEPTH_TEST);
	glutDisplayFunc(render2);
	GLenum err3 = glewInit();
	if (GLEW_OK != err3) {
		fprintf(stderr, "GLEW error");
		return 1;
	}
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(0, 550);
	GLint win4 = glutCreateWindow("graphic2");
	glutReshapeFunc(changeViewPort2);
	glEnable(GL_DEPTH_TEST);
	glutDisplayFunc(render3);
	GLenum err4 = glewInit();
	if (GLEW_OK != err4) {
		fprintf(stderr, "GLEW error");
		return 1;
	}
	glutInitWindowSize(500, 500);
	glutInitWindowPosition(550, 550);
	GLint win5 = glutCreateWindow("graphic3");
	glutReshapeFunc(changeViewPort2);
	glEnable(GL_DEPTH_TEST);
	glutDisplayFunc(render4);
	GLenum err5 = glewInit();
	if (GLEW_OK != err5) {
		fprintf(stderr, "GLEW error");
		return 1;
	}
	glutMainLoop();
	return 0;
}