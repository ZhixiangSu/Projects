// ConsoleApplication12.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
using namespace std;
int prev0[1000], num0[1000], next0[1000];
int key[1005], next1[1005];
int nw = 0; int nx = 1;
int h(int num)
{
	int k[1000], a[1000];
	memset(k, 0, sizeof(k));
	memset(a, 0, sizeof(a));
	int num1 = num,num2=rand();
	for (int i = 0; num1; i++)
	{
		k[i] = num1 % 103;
		a[i] = num2 % 103;
		num1 /= 103;
		num2 /= 103;
	}
	int sum = 0;
	for (int i = 0; k[i]; i++)
	{
		sum += a[i] * k[i]%1005;
	}
	return sum%1005;
}
int push(int num1)
{
	if (next0[nx] == nx)
	{
		cout << "full!" << endl;
		return 0;
	}
	next0[nw] = nx;
	prev0[nx] = nw;
	nw = nx;
	nx = next0[nx];
	num0[nw] = num1;
	int x = h(num1);
	if (key[x])
	{
		while (key[x])
		{
			next1[x] = h(x + num1 * 103 + 1);
			x = h(x + num1 * 103+1);
			
		}
	}
	key[x] = nw;
	return 0;
}
int search(int x)
{
	while (next1[x] != x)
	{
		cout << key[x]<<endl;
		x = next1[x];
	}
	cout << key[x]<<endl;
	return 0;
}
int del(int x)
{
	int num = x;
	while (next1[x] != x)
	{
		if (num0[key[x]] == num)
		{
			next0[prev0[key[x]]] = next0[key[x]];
			prev0[next0[key[x]]] = prev0[key[x]];
			next0[key[x]] = 0; prev0[key[x]] = 0; num0[key[x]] = 0;
			key[x] = 0;
			int x2 = x;
			x = next1[x];
			next1[x2] = x2;
			next0[x] = nx;
			next0[nw] = x;
			nx = x;
		}
	}
	if (num0[key[x]] == num)
	{
		next0[prev0[key[x]]] = next0[key[x]];
		prev0[next0[key[x]]] = prev0[key[x]];
		next0[key[x]] = 0; prev0[key[x]] = 0; num0[key[x]] = 0;
		key[x] = 0;
		key[x] = 0; next1[x] = x;
		next0[x] = nx;
		next0[nw] = x;
		nx = x;
	}
	
	return 0;
}
int main()
{
	memset(prev0, 0, sizeof(prev0));
	memset(num0, 0, sizeof(num0));
	for (int i = 0; i < 999; i++)
	{
		next0[i] = i + 1;
	}
	next0[999] = 999;
	memset(key, 0, sizeof(key));
	for (int i = 0; i <= 999; i++)
	{
		next1[i] = i;
	}
	int i;
	while (cin >> i)
	{
		push(i);
		search(i);
	}
	del(10);
	del(11);
	return 0;
}

