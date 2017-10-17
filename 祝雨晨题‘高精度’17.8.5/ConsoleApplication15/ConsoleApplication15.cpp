// ConsoleApplication12.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
#include<string>
#include<cstring>
#include<cmath>
#include<stack>
#include<algorithm>
#include<cstdio>
using namespace std;
class lll
{
public:
	int wei[5000];
	
	int length = 0; int point = 0;
	void set()
	{
		memset(wei, 0, sizeof(wei));
		length = 0;point =0;
	}
	lll()
	{
		memset(wei, 0, sizeof(wei));
	}
	lll(int a)
	{
		memset(wei, 0, sizeof(wei));
		for (int i = 0; a > 0; i++)
		{
			wei[i] = a % 10;
			a /= 10;
			length++;
		}
	}
	lll(string a)
	{
		memset(wei, 0, sizeof(wei));
		for (int i = 0; i < a.length(); i++)
		{
			if (a[a.length() - i - 1] != '.')
			{
				length++;
				wei[length - 1] = a[a.length() - i - 1] - '0';
			}
			else point = i;
		}
	}
	int out()
	{
		for (int i = 0; i < length; i++)
		{
			if (length - i == point)cout << ".";
			cout << wei[length - i - 1];
		}
		return 0;
	}
};

int jian(lll a, lll b, lll& c);
int jia(lll a, lll b, lll& c)
{
	memset(c.wei, 0, sizeof(c.wei));
	if (b.wei[b.length - 1] < 0)
	{
		b.wei[b.length - 1] = -b.wei[b.length - 1];
		jian(a, b, c);
		return 0;
	}
	int dpoint = a.point - b.point;
	c.length = max(a.point, b.point) + max((a.length - a.point), (b.length - b.point)) + 1;
	c.point = max(a.point, b.point);
	for (int i = min(0, dpoint); i < c.length; i++)
	{
		if (i >= 0)c.wei[i] = a.wei[i];
		if (i - dpoint >= 0)c.wei[i] += b.wei[i - dpoint];
	}
	for (int i = 0; i < c.length; i++)
	{
		if (c.wei[i] >= 10)
		{
			c.wei[i + 1] += c.wei[i] / 10;
			c.wei[i] %= 10;
		}
	}
	while (c.wei[c.length - 1] == 0 && c.length - c.point > 1)c.length--;
	return 0;
}
int jian(lll a, lll b, lll& c)
{
	memset(c.wei, 0, sizeof(c.wei));
	if (b.wei[b.length - 1] < 0)
	{
		b.wei[b.length - 1] = -b.wei[b.length - 1];
		jia(a, b, c);
		return 0;
	}
	int dpoint = a.point - b.point;
	c.length = max(a.point, b.point) + max((a.length - a.point), (b.length - b.point)) + 1;
	c.point = max(a.point, b.point);
	for (int i = min(0, dpoint); i < c.length; i++)
	{
		if (i >= 0)c.wei[i] = a.wei[i];
		if (i - dpoint >= 0)c.wei[i] -= b.wei[i - dpoint];
	}
	for (int i = 0; i < c.length; i++)
	{
		if (c.wei[i] < 0&&i!=c.length-1)
		{
			c.wei[i] += 10;
			c.wei[i + 1]--;
		}
	}
	if (c.wei[c.length - 1] < 0)
	{
		jian(b, a, c);
		c.wei[c.length - 1] = -c.wei[c.length - 1];
	}
	while (c.wei[c.length - 1] == 0 && c.length - c.point > 1)c.length--;
	return 0;
}
int cheng(lll a, lll b, lll &c)
{
	memset(c.wei, 0, sizeof(c.wei));
	c.length = a.length + b.length;
	for (int i = 0; i < a.length; i++)
	{
		for (int j = 0; j < b.length; j++)
		{
			int sum = a.wei[i] * b.wei[j];
			c.wei[i + j] += sum % 10;
			c.wei[i + j + 1] += sum / 10;
		}
	}
	for (int i = 0; i < c.length; i++)
	{
		if (c.wei[i] >= 10)
		{
			c.wei[i + 1] += c.wei[i] / 10;
			c.wei[i] %= 10;
		}
	}
	c.point = a.point + b.point;
	while (c.wei[c.length - 1] == 0 && c.length - c.point > 1)c.length--;
	return 0;
}
int chu(int jingdu, lll a, lll b, lll &c)
{
	c.point =jingdu-(a.length - a.point - b.length + b.point);
	int dpointa = jingdu - a.length, dpointb = jingdu - b.length;
	for (int i = 0; i < a.length; i++)
	{
		if (i + dpointa >= 0)
		{
			a.wei[i + dpointa] = a.wei[i];
			a.wei[i] = 0;
		}
	}
	for (int i = 0; i < b.length; i++)
	{
		if (i + dpointa >= 0)
		{
			b.wei[i + dpointa] = b.wei[i];
			b.wei[i] = 0;
		}
	}
	stack<int>s;
	lll d(1),n(1),e;
	while (s.size() < jingdu)
	{
		lll i(0);
		for (; d.wei[max(d.length-1,0)] > 0; jia(i, n, i))
		{
			cheng(i, b, e);
			jian(a, e, d);
		}
		s.push(i.wei[0]);
		lll shi(10);
		jian(a, e, d);
		cheng(a, shi, a);
	}
	c.length = jingdu;
	for (int i = 0; s.size(); i++)
	{
		c.wei[i] = s.top();
		s.pop();
	}
	return 0;
}
int kaifang(int jingdu, lll a, lll &b)
{
	stack<int> s;
	if ((a.length - a.point) % 2)a.length++;
	b.point = a.point / 2;
	lll num, k;
	lll bai(100),shi(10);
	for (int i = a.length - 1; i >= 0; i -= 2)
	{
		lll shiwei(a.wei[i] * 10);
		if (i > 0)
		{

			cheng(num, bai, num);
			lll gewei(a.wei[i - 1]);
			jia(num, shiwei, num);
			jia(num, gewei, num);
			//num = num * 100 + a.wei[i] * 10 + a.wei[i - 1];
		}
		else
		{
			cheng(num, bai, num);
			jia(num, shiwei, num);
			b.point++;
			//num = num * 100 + a.wei[i] * 10;
		}
		lll j(0); lll n(1); lll o;lll m;
		//jia(k, j, m); cheng(m, j, m); jian(num, m, o);
		while (o.wei[max(0,o.length-1)]>=0)
		{ 
			jia(j, n,j);
			jia(k, j, m);
			cheng(m, j, m);
			jian(num, m, o);
			
		}
		jian(j, n, j);
		jia(k, j, m); cheng(m, j, m); jian(num, m, o);
		jian(num, m, num);
		//num -= (k + j)*j;
		s.push(j.wei[0]);
		cheng(k, shi, k);
		lll ershi(20);
		cheng(j, ershi, j);
		jia(k, j, k);
	}
	while (s.size() < jingdu)
	{
		cheng(num, 100, num);
		lll j(0); lll m; lll n(1); lll o;
		jia(k, j, m); cheng(m, j, m); jian(num, m, o);
		while (o.wei[max(0, o.length - 1)] >= 0)
		{
			jia(j, n, j);
			jia(k, j, m);
			cheng(m, j, m);
			jian(num, m, o);
		}
		jian(j, n, j);
		jia(k, j, m); cheng(m, j, m); jian(num, m, o);
		jian(num, m, num);
		//num -= (k + j)*j;
		s.push(j.wei[0]);
		cheng(k, shi, k);
		lll ershi(20);
		cheng(j, ershi, j);
		jia(k, j, k);
		b.point++;
	}
	b.length = s.size();
	for (int i = 0; s.size(); i++)
	{
		b.wei[i] = s.top();
		s.pop();
	}
	return 0;
}
int main()
{
	string a, b;
	int n;
	cin >> n;	
	//cout << "a:";
	cin >> a;
	//cout << "b:";
	cin >> b;
	lll c(a);
	lll d(b);
	lll e;
	//c.out();
	//d.out();
	//cout << "ans:";
	//cheng(c, d,e);
	//jian(c, d, e);
	/*for (int i = 1; i <= 10000; i++)
	{
		if (i % 2)
		{
			cheng(c, d, c);
			kaifang(n, c, c);
			c.out();
			cout << endl;
		}
		else
		{
			cheng(c, d, d);
			kaifang(n, d, d);         //zyc的题
			d.out();
			cout << endl;
		}
	}*/
	chu(n,c, d, e);
	e.out();
	return 0;
}

