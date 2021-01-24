/*
最小二乘法求线性回归方程实现原理：           相关指数计算原理:                                          均方误差(MSE)计算原理:               ^
公式:                                               ESS           RSS                                            ^               SSE    Σ(y-y) ^2 
        nΣxy —ΣxΣy                       R^2=---------= 1 —-------                                 SSE=Σ(y-y) ^2      MSE=-----=--------------
    k=-------------------     b=ky-x                TSS           TSS                                                             n         n 
      nΣ(x^2)—(Σx)^2
      
	  A=Σxy   B=Σx   C=Σy   D=Σ(x^2)     TSS:总偏差平方和   ESS:回归平方和   RSS:残差平方和         SSE:和方差   MSE:均方误差 

	    nA—BC                                        -              ^ -              - ^
	k=-----------     b=k(C/n)-(B/n)         TSS=Σ(y-y) ^2   ESS=Σ(y-y) ^2   RSS=Σ(y-y) ^2 
	    nD—BB
	n为拟合点数量
	
	经计算，结果(保留20位小数)为:   k=2.83354284078336520000   b=6.69231425482736950000   R^2=0.96139995028534708000   MSE=2.66125669622954540000 
*/

#include <stdio.h>

struct point  //定义一个结构体point(点)，包含x坐标和y坐标 
{
	double x;
	double y;
};

void least_square_apporach(struct point p[],double* target_k,double* target_b);  //最小二乘法求k,b函数 
double R2_apporach(struct point p[],int n,double k,double b);  //R^2函数
double MSE_apporach(struct point p[],int n,double k,double b);  //MSE函数 

int main()
{
	struct point p[20];  //定义一个结构数组，含十个元素，代表输入的十个点，而每个点又都含有一个x坐标，一个y坐标 
	int n=10;  //拟合点数量
	int i=0;
	double k=0;  //线性回归直线斜率
	double b=0;  //线性回归直线纵截距
	double* target_k=&k;  //利用指针target_k储存k的地址
	double* target_b=&b;  //利用指针target_b储存b的地址 
	double R2=0;  //相关指数R^2
	double MSE=0;  //均方误差(MSE)
	freopen("data.txt","r",stdin);  //输入(data)文本文档的数据 
	for(i=0;i<n;i++)
	{
		scanf("%lf",&p[i].x);
	}
	for(i=0;i<n;i++)
	{
		scanf("%lf",&p[i].y);
	}
    /*测试点数据(检验程序是否正常计算)
	p[0].x=111;		p[0].y=124;
	p[1].x=123;     p[1].y=126;
	p[2].x=134;     p[2].y=138;
	p[3].x=147;     p[3].y=154;
	p[4].x=155;     p[4].y=165;
	p[5].x=167;     p[5].y=178;
	p[6].x=178;     p[6].y=185;
	p[7].x=186;     p[7].y=198;
	p[8].x=194;     p[8].y=206;
	p[9].x=207;     p[9].y=215;
	*/
	least_square_apporach(p,&k,&b);
	R2=R2_apporach(p,n,k,b);
	MSE=MSE_apporach(p,n,k,b);
	if(b>0)  //输出保留20位小数 
	{
		printf("y=%.20fx+%.20f   R2=%.20f   MSE=%.20f\n",k,b,R2,MSE);
	}
	else
	{
		if(b<0)
		{
			printf("y=%.20fx%.20f   R2=%.20f   MSE=%.20f\n",k,b,R2,MSE);
		}
		if(b=0)
		{
			printf("y=%.20fx   R2=%.20f   MSE=%.20f\n",k,R2,MSE);
		}
	}
	return 0;
}

void least_square_apporach(struct point p[],double* target_k,double* target_b)
{
	int n=10;  //拟合点数量 
	int i=0;
	double k=0;
	double b=0;
	double A=0;
	double B=0;
	double C=0; 
	double D=0;
	for(i=0;i<n;i++)
	{
		A+=p[i].x*p[i].y;  //Σxy 
		B+=p[i].x;  //Σx
		C+=p[i].y;  //Σy  
		D+=p[i].x*p[i].x;  //Σ(x^2) 
	}
	k=(n*A-B*C)/(n*D-B*B);
	b=(C/n)-k*(B/n);
	*target_k=k;  //通过指针传回k 
	*target_b=b;  //通过指针传回b
}

double R2_apporach(struct point p[],int n,double k,double b)
{
	int i=0;
	double R2=0;
	double TSS=0;  //总偏差平方和 
	double ESS=0;  //回归平方和 
	double C=0;  //Σy 
	double average_C=0;  //Σy/n; 
	double yi=0;  //根据线性回归方程得出的估计值 
	for(i=0;i<n;i++)
	{
		C+=p[i].y;
	}
	average_C=C/n;
	for(i=0;i<n;i++)
	{
		TSS+=(p[i].y-average_C)*(p[i].y-average_C);
	}
	for(i=0;i<n;i++)
	{
		yi=0;  //重新初始化yi; 
		yi=k*p[i].x+b;
		ESS+=(yi-average_C)*(yi-average_C);
	}
	R2=ESS/TSS; 
	return R2;
}

double MSE_apporach(struct point p[],int n,double k,double b)
{
	double SSE=0;  //和方差 
	double MSE=0;  //均方误差 
	int i=0;
	double yi=0;  //根据线性回归方程得出的估计值 
	for(i=0;i<n;i++)
	{
		yi=0;
		yi=k*p[i].x+b;
		SSE+=(p[i].y-yi)*(p[i].y-yi);
	}
	MSE=SSE/n;
	return MSE;
}











