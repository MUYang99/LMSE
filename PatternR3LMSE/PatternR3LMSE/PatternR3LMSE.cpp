// PatternR3LMSE.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include "math.h"
#include<iostream>

using namespace std;

void mt(float *a, float *b, int row, int column)//矩阵转置函数
{
	int i, j;//a是转置前矩阵 b是转置后的矩阵 row为矩阵的行数 column是矩阵的列数
	float temp;
	for (i = 0; i<row; i++)
		for (j = 0; j<column; j++)
		{
			temp = *(a + i*column + j);
			*(b + j*row + i) = temp;
		}
}

void bmul(float *a, float *b, int m, int n, int k, float *result)//矩阵相乘矩阵
{
	int i, j, l, u;              //矩阵result=矩阵a*矩阵b
	for (i = 0; i <= m - 1; i++)    //m是矩阵a的行数，n是矩阵a的列数 k是矩阵b的列数
		for (j = 0; j <= k - 1; j++)
		{
			u = i*k + j;
			result[u] = 0.0;
			for (l = 0; l <= n - 1; l++)
				result[u] = result[u] + a[i*n + l] * b[l*k + j];
		}
	return;
}

int inv(float *a, int n) //矩阵求逆函数 n是矩阵a的阶数
{
	int *is, *js, i, j, k, l, u, v;
	float d, p;

	is = new int[n];
	js = new int[n];

	for (k = 0; k <= n - 1; k++)
	{
		d = 0.0;
		for (i = k; i <= n - 1; i++)
			for (j = k; j <= n - 1; j++)
			{
				l = i*n + j;
				p = fabs(a[l]);
				if (p>d) { d = p; is[k] = i; js[k] = j; }
			}
		if (d + 1.0 == 1.0)
		{
			delete is;
			delete js;
			cout << "error matrix can not inv" << endl;
			return(0);
		}
		if (is[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k*n + j; v = is[k] * n + j;
				p = a[u]; a[u] = a[v]; a[v] = p;
			}
		if (js[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i*n + k; v = i*n + js[k];
				p = a[u]; a[u] = a[v]; a[v] = p;
			}
		l = k*n + k;
		a[l] = 1.0 / a[l];
		for (j = 0; j <= n - 1; j++)
			if (j != k)
			{
				u = k*n + j; a[u] = a[u] * a[l];
			}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
				for (j = 0; j <= n - 1; j++)
					if (j != k)
					{
						u = i*n + j;
						a[u] = a[u] - a[i*n + k] * a[k*n + j];
					}
		for (i = 0; i <= n - 1; i++)
			if (i != k)
			{
				u = i*n + k; a[u] = -a[u] * a[l];
			}
	}
	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
			for (j = 0; j <= n - 1; j++)
			{
				u = k*n + j; v = js[k] * n + j;
				p = a[u]; a[u] = a[v]; a[v] = p;
			}
		if (is[k] != k)
			for (i = 0; i <= n - 1; i++)
			{
				u = i*n + k; v = i*n + is[k];
				p = a[u]; a[u] = a[v]; a[v] = p;
			}
	}
	delete is; delete js;
	return(1);
}

void subb(float *a, float *b, int n, float *c)//矩阵相减函数
{                                          //c中存放差矩阵
	int i;
	for (i = 0; i<n; i++)
	{
		*(c + i) = *(a + i) - *(b + i);
	}
}

void kmul(float *a, int n, int c, float *b)//矩阵数乘函数
{
	int i;
	for (i = 0; i<n; i++)
	{
		*(b + i) = (*(a + i))*c;
	}
}

void madd(float *a, float *b, int n, float *c)//矩阵相加函数
{                                          //c中存放和矩阵
	int i;
	for (i = 0; i<n; i++)
	{
		*(c + i) = *(a + i) + *(b + i);//fabs是浮点数绝对值函数
	}
}

void main()
{
	float x[8][4] = { { 0, 0, 0, 1 }, { 1, 0, 0, 1 }, { 1, 0, 1, 1 }, { 1, 1, 0, 1 }, { 0, 0, -1, -1 }, { 0, -1, -1, -1 }, { 0, -1, 0, -1 }, { -1, -1, -1, -1 } };//正规化
	float xt[4][8];
	mt(&x[0][0], &xt[0][0], 8, 4);
	float xm[4][4];
	bmul(&xt[0][0], &x[0][0], 4, 8, 4, &xm[0][0]);
	inv(&xm[0][0], 4);
	float xn[4][8];
	bmul(&xm[0][0], &xt[0][0], 4, 4, 8, &xn[0][0]);//求得规范矩阵x+

	float b[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
	float w[4];
	bmul(&xn[0][0], &b[0], 4, 8, 1, &w[0]);//求得权向量w
	
	cout << "第一题 \n" << endl;
	cout << "解向量为:(" << w[0] << "," << w[1] << "," << w[2] << "," << w[3] << ")\n" << endl;//输出解向量
	
	
	float x1[4][3] = { { 0, 1, 1 }, { 0, -1, 1 }, { -1, 0, -1 }, { 1, 0, -1 } };//正规化
	float xt1[3][4];
	mt(&x1[0][0], &xt1[0][0], 4, 3);
	float xm1[3][3];
	bmul(&xt1[0][0], &x1[0][0], 3, 4, 3, &xm1[0][0]);
	inv(&xm1[0][0], 3);
	float xn1[3][4];
	bmul(&xm1[0][0], &xt1[0][0], 3, 3, 4, &xn1[0][0]);//求得规范矩阵x+

	float b1[4] = { 1, 1, 1, 1 };
	float w1[3];
	bmul(&xn1[0][0], &b1[0], 3, 4, 1, &w1[0]);//求得权向量w
	cout << "第二题 \n" << endl;
	cout << "解向量为:(" << w1[0] << "," << w1[1] << "," <<  w1[2] << ")" << endl;//输出解向量
	
	float xw[4];
	float e[4];
	float ef[4];
	float wa[3];
	float ba[4];
	int i,f0, f1, f2;
	int c = 2;
	
	while (1)
	{
		bmul(&x1[0][0], &w1[0], 4, 3, 1, &xw[0]);
		subb(&xw[0], &b1[0], 4, &e[0]); //e的值
		//根据e的值进行判别
		f0 = 0; f1 = 0; f2 = 0;
		for (i = 0; i<4; i++)
		{
			if (e[i] == 0) f0++;
			else if (e[i]>0) f1++;
			else f2++;
			ef[i] = fabs(e[i]);
		}
		if (f0 == 4) //e=0得解
			break;
		if (f1>0) //e>0继续迭代
		{
			bmul(&xn1[0][0], &ef[0], 3, 4, 1, &wa[0]);
			kmul(&wa[0], 3, c, &wa[0]);
			madd(&w1[0], &wa[0], 3, &w1[0]); //修正w

			madd(&e[0], &ef[0], 4, &ba[0]);
			kmul(&ba[0], 4, c, &ba[0]);
			madd(&b1[0], &ba[0], 4, &b1[0]); //修正b
		}
		if (f2 == 4) //e的全部分量为负值，模式线性不可分
		{
			cout << "模式线性不可分!" << endl;
			break;
		}
	}
	getchar();
}