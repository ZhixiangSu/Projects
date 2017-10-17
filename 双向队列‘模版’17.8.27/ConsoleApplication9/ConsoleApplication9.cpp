// ConsoleApplication9.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
using namespace std;
class sq
{
public:
	 sq()
	{
		 s = 0;
		 e = 0;
		 length = 0;
		 memset(l, 0, sizeof(l));
	}
	int spush(int x)
	{
		if (length > 1000)cout << "full!" << endl;
		s = (s - 1) % 1000;
		l[s] = x;
		length++;
		return 0;
	}
	int spop()
	{
		if (!length)
		{
			cout << "empty!" << endl;
			return 0;
		}
		l[s] = 0;
		s = (s + 1) % 1000;
		length--;
		return 0;
	}
	int qpush(int x)
	{
		if (length > 1000)cout << "full!" << endl;
		l[e] = x;
		e = (e + 1) % 1000;
		length++;
		return 0;
	}
	int qpop()
	{
		spop();
		return 0;
	}
private:
	int l[1000];
	int s;
	int e;
	int length;
};
int main()
{
	sq a;
	a.qpop();
	a.qpush(1);
	a.spop();
	a.spush(2);
	a.qpop();
    return 0;
}

