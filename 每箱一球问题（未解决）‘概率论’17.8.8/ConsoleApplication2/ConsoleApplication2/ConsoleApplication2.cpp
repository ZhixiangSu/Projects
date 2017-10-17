// ConsoleApplication2.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
using namespace std;
int a[10000], b[10000], c[10000];
int main()
{//	cout << n << ": ";
	int n;
	cin >> n;
	cout << endl;
	memset(c, 0, sizeof(c));
	for (int i = 0; i < 10000; i++)
	{
		srand(i);
		memset(b, 0, sizeof(b));
		int j = 0, k = 0;
		for (; !k; j++)
		{
			b[rand() % n] = 1; k = 1;
			for (int l = 0; l < n; l++)
			{
				if (!b[l])k = 0;
			}
		}
		a[i] = j;
		c[j]++;
	}
	double sum = 0;
	for (int i = 0; i < 10000; i++)
	{
		sum += a[i];
	}
	sum /= 10000;
	cout << sum << endl;
	for (int i = 1; i < 10000; i++)
		if(c[i])cout << c[i] << " ";
	cout << endl;

	return 0;
}

