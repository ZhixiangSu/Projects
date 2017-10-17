// ConsoleApplication1.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
#include<algorithm>
using namespace std;
int d[1000];
int dfs(int s, int e)
{
	int dmax = 0;
	if (s < e)
	{
		int mid = (s + e) / 2;
		dmax=max(dfs(s, mid),dfs(mid + 1, e));
		int sum = 0,da=0,db=0;
		for (int i = mid; i >= s; i--)
		{
			sum += d[i];
			da = max(da, sum);
		}
		sum = 0;
		for (int j = mid + 1; j <= e; j++)
		{
			sum += d[j];
			db = max(db, sum);
		}
		dmax = max(dmax, da + db);
		return dmax;
	}
	else return d[s];
}
int main()
{
	int n;
	cin >> n;
	int a=0, b=0;
	for (int i = 0; i < n; i++)
	{
		a = b;
		cin >> b;
		d[i] = b - a;
	}
	cout<<dfs(1, n-1);
    return 0;
}

