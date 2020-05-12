// PatternR3LMSE.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "stdafx.h"
#include "math.h"
#include<iostream>

using namespace std;

void mt(float *a, float *b, int row, int column)//����ת�ú���
{
	int i, j;//a��ת��ǰ���� b��ת�ú�ľ��� rowΪ��������� column�Ǿ��������
	float temp;
	for (i = 0; i<row; i++)
		for (j = 0; j<column; j++)
		{
			temp = *(a + i*column + j);
			*(b + j*row + i) = temp;
		}
}

void bmul(float *a, float *b, int m, int n, int k, float *result)//������˾���
{
	int i, j, l, u;              //����result=����a*����b
	for (i = 0; i <= m - 1; i++)    //m�Ǿ���a��������n�Ǿ���a������ k�Ǿ���b������
		for (j = 0; j <= k - 1; j++)
		{
			u = i*k + j;
			result[u] = 0.0;
			for (l = 0; l <= n - 1; l++)
				result[u] = result[u] + a[i*n + l] * b[l*k + j];
		}
	return;
}

int inv(float *a, int n) //�������溯�� n�Ǿ���a�Ľ���
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

void subb(float *a, float *b, int n, float *c)//�����������
{                                          //c�д�Ų����
	int i;
	for (i = 0; i<n; i++)
	{
		*(c + i) = *(a + i) - *(b + i);
	}
}

void kmul(float *a, int n, int c, float *b)//�������˺���
{
	int i;
	for (i = 0; i<n; i++)
	{
		*(b + i) = (*(a + i))*c;
	}
}

void madd(float *a, float *b, int n, float *c)//������Ӻ���
{                                          //c�д�ź;���
	int i;
	for (i = 0; i<n; i++)
	{
		*(c + i) = *(a + i) + *(b + i);//fabs�Ǹ���������ֵ����
	}
}

void main()
{
	float x[8][4] = { { 0, 0, 0, 1 }, { 1, 0, 0, 1 }, { 1, 0, 1, 1 }, { 1, 1, 0, 1 }, { 0, 0, -1, -1 }, { 0, -1, -1, -1 }, { 0, -1, 0, -1 }, { -1, -1, -1, -1 } };//���滯
	float xt[4][8];
	mt(&x[0][0], &xt[0][0], 8, 4);
	float xm[4][4];
	bmul(&xt[0][0], &x[0][0], 4, 8, 4, &xm[0][0]);
	inv(&xm[0][0], 4);
	float xn[4][8];
	bmul(&xm[0][0], &xt[0][0], 4, 4, 8, &xn[0][0]);//��ù淶����x+

	float b[8] = { 1, 1, 1, 1, 1, 1, 1, 1 };
	float w[4];
	bmul(&xn[0][0], &b[0], 4, 8, 1, &w[0]);//���Ȩ����w
	
	cout << "��һ�� \n" << endl;
	cout << "������Ϊ:(" << w[0] << "," << w[1] << "," << w[2] << "," << w[3] << ")\n" << endl;//���������
	
	
	float x1[4][3] = { { 0, 1, 1 }, { 0, -1, 1 }, { -1, 0, -1 }, { 1, 0, -1 } };//���滯
	float xt1[3][4];
	mt(&x1[0][0], &xt1[0][0], 4, 3);
	float xm1[3][3];
	bmul(&xt1[0][0], &x1[0][0], 3, 4, 3, &xm1[0][0]);
	inv(&xm1[0][0], 3);
	float xn1[3][4];
	bmul(&xm1[0][0], &xt1[0][0], 3, 3, 4, &xn1[0][0]);//��ù淶����x+

	float b1[4] = { 1, 1, 1, 1 };
	float w1[3];
	bmul(&xn1[0][0], &b1[0], 3, 4, 1, &w1[0]);//���Ȩ����w
	cout << "�ڶ��� \n" << endl;
	cout << "������Ϊ:(" << w1[0] << "," << w1[1] << "," <<  w1[2] << ")" << endl;//���������
	
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
		subb(&xw[0], &b1[0], 4, &e[0]); //e��ֵ
		//����e��ֵ�����б�
		f0 = 0; f1 = 0; f2 = 0;
		for (i = 0; i<4; i++)
		{
			if (e[i] == 0) f0++;
			else if (e[i]>0) f1++;
			else f2++;
			ef[i] = fabs(e[i]);
		}
		if (f0 == 4) //e=0�ý�
			break;
		if (f1>0) //e>0��������
		{
			bmul(&xn1[0][0], &ef[0], 3, 4, 1, &wa[0]);
			kmul(&wa[0], 3, c, &wa[0]);
			madd(&w1[0], &wa[0], 3, &w1[0]); //����w

			madd(&e[0], &ef[0], 4, &ba[0]);
			kmul(&ba[0], 4, c, &ba[0]);
			madd(&b1[0], &ba[0], 4, &b1[0]); //����b
		}
		if (f2 == 4) //e��ȫ������Ϊ��ֵ��ģʽ���Բ��ɷ�
		{
			cout << "ģʽ���Բ��ɷ�!" << endl;
			break;
		}
	}
	getchar();
}