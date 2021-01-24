/*
��С���˷������Իع鷽��ʵ��ԭ��           ���ָ������ԭ��:                                          �������(MSE)����ԭ��:               ^
��ʽ:                                               ESS           RSS                                            ^               SSE    ��(y-y) ^2 
        n��xy ����x��y                       R^2=---------= 1 ��-------                                 SSE=��(y-y) ^2      MSE=-----=--------------
    k=-------------------     b=ky-x                TSS           TSS                                                             n         n 
      n��(x^2)��(��x)^2
      
	  A=��xy   B=��x   C=��y   D=��(x^2)     TSS:��ƫ��ƽ����   ESS:�ع�ƽ����   RSS:�в�ƽ����         SSE:�ͷ���   MSE:������� 

	    nA��BC                                        -              ^ -              - ^
	k=-----------     b=k(C/n)-(B/n)         TSS=��(y-y) ^2   ESS=��(y-y) ^2   RSS=��(y-y) ^2 
	    nD��BB
	nΪ��ϵ�����
	
	�����㣬���(����20λС��)Ϊ:   k=2.83354284078336520000   b=6.69231425482736950000   R^2=0.96139995028534708000   MSE=2.66125669622954540000 
*/

#include <stdio.h>

struct point  //����һ���ṹ��point(��)������x�����y���� 
{
	double x;
	double y;
};

void least_square_apporach(struct point p[],double* target_k,double* target_b);  //��С���˷���k,b���� 
double R2_apporach(struct point p[],int n,double k,double b);  //R^2����
double MSE_apporach(struct point p[],int n,double k,double b);  //MSE���� 

int main()
{
	struct point p[20];  //����һ���ṹ���飬��ʮ��Ԫ�أ����������ʮ���㣬��ÿ�����ֶ�����һ��x���꣬һ��y���� 
	int n=10;  //��ϵ�����
	int i=0;
	double k=0;  //���Իع�ֱ��б��
	double b=0;  //���Իع�ֱ���ݽؾ�
	double* target_k=&k;  //����ָ��target_k����k�ĵ�ַ
	double* target_b=&b;  //����ָ��target_b����b�ĵ�ַ 
	double R2=0;  //���ָ��R^2
	double MSE=0;  //�������(MSE)
	freopen("data.txt","r",stdin);  //����(data)�ı��ĵ������� 
	for(i=0;i<n;i++)
	{
		scanf("%lf",&p[i].x);
	}
	for(i=0;i<n;i++)
	{
		scanf("%lf",&p[i].y);
	}
    /*���Ե�����(��������Ƿ���������)
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
	if(b>0)  //�������20λС�� 
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
	int n=10;  //��ϵ����� 
	int i=0;
	double k=0;
	double b=0;
	double A=0;
	double B=0;
	double C=0; 
	double D=0;
	for(i=0;i<n;i++)
	{
		A+=p[i].x*p[i].y;  //��xy 
		B+=p[i].x;  //��x
		C+=p[i].y;  //��y  
		D+=p[i].x*p[i].x;  //��(x^2) 
	}
	k=(n*A-B*C)/(n*D-B*B);
	b=(C/n)-k*(B/n);
	*target_k=k;  //ͨ��ָ�봫��k 
	*target_b=b;  //ͨ��ָ�봫��b
}

double R2_apporach(struct point p[],int n,double k,double b)
{
	int i=0;
	double R2=0;
	double TSS=0;  //��ƫ��ƽ���� 
	double ESS=0;  //�ع�ƽ���� 
	double C=0;  //��y 
	double average_C=0;  //��y/n; 
	double yi=0;  //�������Իع鷽�̵ó��Ĺ���ֵ 
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
		yi=0;  //���³�ʼ��yi; 
		yi=k*p[i].x+b;
		ESS+=(yi-average_C)*(yi-average_C);
	}
	R2=ESS/TSS; 
	return R2;
}

double MSE_apporach(struct point p[],int n,double k,double b)
{
	double SSE=0;  //�ͷ��� 
	double MSE=0;  //������� 
	int i=0;
	double yi=0;  //�������Իع鷽�̵ó��Ĺ���ֵ 
	for(i=0;i<n;i++)
	{
		yi=0;
		yi=k*p[i].x+b;
		SSE+=(p[i].y-yi)*(p[i].y-yi);
	}
	MSE=SSE/n;
	return MSE;
}











