// ConsoleApplication8.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
#include<algorithm>
using namespace std;
int x[1000], y[1000];int n;
int less_x(int num)
{
	int s = 0, e = n - 1;
	while (e > s+1)
	{
		int m = (s + e) / 2;
		if (y[m] > num)
		{
			e = m;
		}
		else if (y[m] < num)
		{
			s = m;
		}
	}
	if (y[0] > num)return -1;
	return s;
}
int dsort_x(int s, int e, int m)
{
	int mid = (s + e) / 2;
	int mid2 = less_x(x[mid]);
	if (mid + mid2 == m-2)return x[mid];
	else if (s == e)return 0;
	else if (mid + mid2 > m-2)return dsort_x(s, mid, m);
	else return dsort_x(mid + 1, e, m);
}
int main()
{
	
	cin >> n;
	for (int i = 0; i < n; i++)
	{
		x[i] = rand();
		y[i] = rand();
	}
	sort(x, x+n);
	sort(y, y+n);
	int w=dsort_x(0, n - 1, n);
	if (w)cout << w;
	else
	{
		for (int i = 0; i < n; i++)
		{
			swap(x[i], y[i]);
		}
		cout << dsort_x(0, n - 1, n);
	}
    return 0;
}

