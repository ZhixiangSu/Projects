// ConsoleApplication1.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
using namespace std;
class young
{
public:
	int size, q[100][100];
	void dsort(int i, int j)
	{
		int imin = i, jmin = j;
		if (i + 1 < size&&q[i + 1][j] < q[imin][jmin])
		{
			imin = i + 1; jmin = j;
		}
		if (j + 1 < size&&q[i][j + 1] < q[imin][jmin])
		{
			imin = i; jmin = j + 1;
		}
		if (imin != i || jmin != j)
		{
			swap(q[i][j], q[imin][jmin]);
			dsort(imin, jmin);
		}
	}
	void gouzao()
	{
		for (int i = size - 1; i >= 0; i--)
		{
			for (int j = size - 1; j >= 0; j--)
			{
				dsort(i, j);
			}
		}
	}
	int check()
	{

	}
};
int main()
{
	young a;
	cin >> a.size;
	srand(a.size);
	for (int i = 0; i < a.size; i++)
	{
		for(int j = 0; j < a.size; j++)
		{
			a.q[i][j] = rand();
		}
	}
	a.gouzao();

	/*for(int i = 0; i < a.size; i++)
	{
		for(int j = 0; j < a.size; j++)
		{
			cout << a.q[i][j] << " ";
		}
		cout << endl;
	}*/

    return 0;
}

