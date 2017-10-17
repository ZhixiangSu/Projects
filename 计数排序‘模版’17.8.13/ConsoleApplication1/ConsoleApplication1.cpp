// ConsoleApplication1.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
#include<algorithm>
using namespace std;
int a[1000], b[100000], c[1000];
int main()
{
	int n;
	cin >> n;
	srand(n);
	int max1 = 0;
	for (int i = 0; i < n; i++)
	{
		a[i]=rand();
		b[a[i]]++;
		max1 = max(max1, a[i]);
	}
	for (int i = 0; i <= max1; i++)
	{
		b[i] += b[i - 1];
	}
	memset(c, 0, sizeof(c));
	for (int i = n-1; i >=0; i--)
	{
		c[b[a[i]]-1] = a[i];
		b[a[i]]--;
	}
    return 0;
}

