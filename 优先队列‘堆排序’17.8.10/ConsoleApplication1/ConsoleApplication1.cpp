// ConsoleApplication1.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
using namespace std;

class queue_int
{
public:
	int length=0,q[1000];
	void dsort(int i)
	{
		int l = i * 2 + 1, r = i * 2 + 2;
		int max = i;
		if (l<length&&q[l] > q[max])
		{
			max = l;
		}
		if (r<length&&q[r] > q[max])
		{
			max = r;
		}
		if (max != i)
		{
			swap(q[i], q[max]);
			dsort(max);
		}
	}
	void gouzao()
	{
		for (int i = length / 2; i >= 0; i--)
		{
			dsort(i);
		}
	}
	void increase(int i, int num)
	{
		q[i] = num;
		while (i&&q[i / 2] < q[i])
		{
			swap(q[i], q[i / 2]);
			i /= 2;
		}
	}
	void insert(int i)
	{
		length++;
		increase(length - 1, i);
		
	}
	int maximum()
	{
		return q[0];
	}
	int extract_max()
	{
		int m = q[0];
		q[0] = q[length - 1];
		length--;
		dsort(0);
		return m;
	}
	void delete1(int i)
	{
		q[i] = q[length - 1];
		length--;
		dsort(i);
	}
	void check()
	{
		for (int i = 0; i < length / 2; i++)
		{
			int l = i * 2 + 1, r = i * 2 + 2;
			if (l<length&&q[l]>q[i])
			{
				cout << "error!"<<endl;
				cout << q[i] << " " << q[l] << endl;
			}
			if (r<length&&q[r]>q[i])
			{
				cout << "error!" << endl;
				cout << q[i] << " " << q[r] << endl;
			}
		}
	}
};
int main()
{
	queue_int a;
	cin >> a.length;
	srand(a.length);
	for (int i = 0; i <a.length; i++)
	{
		a.q[i] = rand();
	}
	a.gouzao();
	a.insert(10000);
	a.extract_max();
	a.delete1(5);
	a.check();
    return 0;
}

