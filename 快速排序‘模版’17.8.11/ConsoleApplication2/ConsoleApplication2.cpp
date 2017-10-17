// ConsoleApplication2.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
using namespace std;
int a[10000];
int check(int s, int e);
int qsort(int s,int e)
{
	if (s < e)
	{
		int m = check(s, e);
		qsort(s, m - 1);
		qsort(m + 1, e);
	}
	
	return 0;
}
int check(int s, int e)
{
	int x = rand() % (e - s+1);
	swap(a[e], a[e - x]);
	int i = s-1;
	for (int j = s; j < e; j++)
	{
		if (a[j] <a[e])
		{
			i++;
			swap(a[j], a[i]);
		}
	}
	swap(a[i + 1], a[e]);
	return i + 1;
}
int main()
{
	int n;
	cin >> n;
	srand(n);
	for (int i = 0; i < n; i++)
	{
		a[i] = rand();
	}
	qsort(0,n-1);
    return 0;
}

