// ConsoleApplication1.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
using namespace std;
int a[10000];
int dsort(int i, int n)
{
	int max = i;
	int l = i * 2 + 1, r = i * 2 + 2;
	if (l<n&&a[l] > a[max])
	{
		max = l;
	}
	if (r<n&&a[r] > a[max])
	{
		max = r;
	}
	if (max != i)
	{
		swap(a[i], a[max]);
		dsort(max, n);
	}
	return 0;
}
int main()
{
	int n;
	cin >> n;
	for (int i = 0; i < n; i++)
	{
		a[i]=rand();
	}
	for (int i = n / 2; i >= 0; i--)
	{
		dsort(i, n);
	}
	for (int i = n - 1; i > 0; i--)
	{
		swap(a[0], a[i]);
		dsort(0, i);
	}
    return 0;
}

