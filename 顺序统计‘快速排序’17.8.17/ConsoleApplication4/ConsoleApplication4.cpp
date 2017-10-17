// ConsoleApplication4.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
using namespace std;
int a[1000];
int search(int s, int e)
{
	int j = s - 1;
	for (int i = s; i < e; i++)
	{
		if (a[i] < a[e])
		{
			j++;
			swap(a[i], a[j]);
		}
	}
	swap(a[e], a[j + 1]);
	return j + 1;
}
int qsort(int s,int e,int m)
{

	int mid = search(s,e);
	if (mid == m)
	{
		cout << a[m] << endl;
		return 0;
	}
	if(mid>m)qsort(s, mid - 1,m);
	else if(mid<m)qsort(mid + 1, e,m);
	return 0;
}
int main()
{
	int n, m;
	cin >> n >> m;
	srand(n);
	for (int i = 0; i < n; i++)
	{
		a[i] = rand();
	}
	qsort(0,n-1,m-1);
    return 0;
}

