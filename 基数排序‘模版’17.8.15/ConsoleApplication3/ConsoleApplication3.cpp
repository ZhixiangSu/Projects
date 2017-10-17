// ConsoleApplication3.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
#include<algorithm>
using namespace std;
int a[10000],b[10000],c[11],d[10000];
int jsort(int n)
{
	memset(c, 0, sizeof(c));
	memset(d, 0, sizeof(d));
	for (int i = 0; i < n; i++)
	{
		c[b[i]]++;
	}
	for (int i = 0; i < 10; i++)
	{
		c[i] += c[i - 1];
	}
	for (int i = n - 1; i >= 0; i--)
	{
		d[c[b[i]]] = a[i];
		c[b[i]]--;
	}
	for (int i = 0; i < n; i++)
	{
		a[i] = d[i+1];
	}

	return 0;
}
int main()
{
	int n;
	cin >> n;
	srand(n);
	int max0 = 0;
	for (int i = 0; i < n; i++)
	{
		a[i] = rand();
		max0 = max(max0, a[i]);
	}
	int m = 0;
	while (max0)
	{
		max0 /= 10;
		m++;
		for (int i = 0; i < n; i++)
		{
			b[i] =a[i]% (int)pow(10, m);
			b[i] /=(int) pow(10, m - 1);
		}
		jsort(n);
	}
    return 0;
}