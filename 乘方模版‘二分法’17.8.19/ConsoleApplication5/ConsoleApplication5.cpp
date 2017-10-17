// ConsoleApplication5.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include<iostream>
using namespace std;
int power(int a, int n)
{
	int sum = 1; if (n)
	{
		sum = pow(a, n / 2);
		sum *= sum;
		if (n % 2)sum *= a;
	}

	return sum;
}
int main()
{
	int a, n;
	cin >> a >> n;
	cout << power(a, n);
	return 0;
}

