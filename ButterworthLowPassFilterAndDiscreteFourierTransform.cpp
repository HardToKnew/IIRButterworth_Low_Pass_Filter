// Discrete Fourier Transform.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <string.h>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <bitset>
#include <complex>
#include <math.h>
using namespace std;
#define low 1
#define high 2


struct BUTTORD{
	int N;
	double  Wn;
};
struct BUTTER {
	double  A;
	double  B;
};

struct BUTTAP{
	complex < double > z;
	complex < double > p;
	complex < double > k;
};

#define pi 3.1415926
void  fft(complex < double >* Data, int  Log2N, int  flag);
void  openfile(string filename, double* Arraty,complex < double >* Data);
BUTTORD buttord(double, double, double, double);
//BUTTER* butter(int n, double Wn, int filter, double *az,double *bz);//filter=0 lowpass 1 highpass
BUTTER* butter(int N, double Wn, int filter,double* az, double* bz);
int Complex_Multiple(complex<double> a, complex<double> b
	, complex<double>* Res_Real);
int Bilinear(int N,
	double *as, double *bs,
	double *az, double *bz);
int Filter(double *pdAz,
	double *pdBz,	//滤波器参数表2
	double *dDataIn,//输入数据
	int	nDataLen,
	int nABLen,
	double *DataOut//数据输出
);
double FiltButter(double *pdAz,	//滤波器参数表1
	double *pdBz,	//滤波器参数表2
	int	nABLen,	//参数序列的长度
	double dDataIn,//输入数据
	double *pdBuf);




int main()
{
	string Filename = "C:\\Users\\lenovo\\Desktop\\zizaodata.txt";
	double Arraty[20480];
	complex < double > Data[20480];
	
	openfile(Filename, Arraty, Data);
	fft( Data,1024, -1);
	//输出为文件
	FILE *FourierBackgroundfile;
	FourierBackgroundfile = fopen("C:\\Users\\lenovo\\Desktop\\Data\\cpptemp.txt", "w+");
	for (int i = 0; i < 1024; i++) {
		fprintf(FourierBackgroundfile, "%d\t%f %f\n", i, Data[i].real(), Data[i].imag());
	}
	fclose(FourierBackgroundfile);
	double  Stopband_attenuation = 30;
	double	Passband_attenuation = 1;
	double Stopband = 0.4 , Passband = 0.3 ;
	Passband = 0.2;
	Stopband = 0.4;
	Passband_attenuation = 1;
	Stopband_attenuation = 30;
	int N;
	
	char num;
	cout << "c 继续"<<endl;
	cin >> num;
	while (num != 'q') {

		cout << "输入0-1之间的数字频率，如0.1代表0.1pi" << endl;
		cout << "Passband:";
		cin >> Passband;
		cout << "Stopband:";
		cin >> Stopband;
		cout << "Passband_attenuation :";
		cin >>Passband_attenuation;
		cout << "Stopband_attenuation: ";
		cin >> Stopband_attenuation;
		BUTTORD Buttord = buttord(Passband, Stopband, Passband_attenuation, Stopband_attenuation);
		BUTTER* ButterZ = new BUTTER[Buttord.N + 1];
		double *az = new double[5];
		double *bz = new double[5];
		ButterZ=butter(Buttord.N, Buttord.Wn, low,az,bz);
		double* pdBuf = new double[1025];
		//double
		/*for (int i = 0; i < 1024; i++) {*/
		double* az1 = new double[5];
			*az1=1; *(az1+1)= 0.982564079238018; *(az1 + 2)= 0.792712081405628; 
			*(az1 + 3) = 0.243240475352550; *(az1 + 4) = 0.037920407142041;
		double *bz1 = new double[5]; 
		*bz1 = 0.191027315196140; *(bz1 + 1) = 0.764109260784559; *(bz1 + 2) = 1.146163891176839;
		*(bz1 + 3) = 0.764109260784559; *(bz1 + 4) = 0.191027315196140;
			//cout << (Filter(az1, bz1, Buttord.N, Arraty[0], pdBuf)) << " ";
			Filter(az1, bz1,  Arraty,1025, Buttord.N, pdBuf);
		//}*/
		cout << endl<<"c 继续" << endl;
		cin >> num;
		if(num == 'q')
			break;
	}
	
	
}


void   openfile(string filename, double* Arraty, complex < double >* Data)//数据读出
{
	int datastart, dataend;
	string strbuf[4096];
	ifstream fin;
	fin.open(filename, ios::in);
	int j = 0, count = 0;
	while (!fin.eof() && j < 40960)//文件是否读完，行数是否小于数组
	{
		string inbuf;
		getline(fin, inbuf, '\n');
		strbuf[j] = inbuf;
		j = j + 1;
		if (inbuf.substr(0, 8) == "<<DATA>>") {
			datastart = j;
		}
		if (inbuf.substr(0, 7) == "<<END>>")
		{
			dataend = j;
		}
	}
	if (datastart&&dataend) {
		for (int i = datastart; i < dataend - 1; i++) {
			Arraty[count] = stof(strbuf[i]);
			Data[count] = complex < double >(Arraty[count]);
			count++;
		}
	}
}

void  fft(complex < double >* Data,//输入数据返回（逆）傅里叶变换结果
	int  Log2N,//数据长度 
	int  flag)//1 反变换
{
	int i, j, length;
	complex<double> wn;
	length = Log2N;
	complex<double>* temp= new complex<double>[length];
	for (i = 0; i < length; i++)
	{
		temp[i] = 0;

		for (j = 0; j < length; j++)
		{
			wn = complex<double>(cos(2.0*pi / length * i*j), flag*sin(2.0*pi / length * i*j));
			temp[i] += Data[j] * wn;
		}
	}
	/*if(flag==1)
	{
		for(i=0;i<length;i++)
			Data[i]=temp[i]/complex<double> (length);
	}*/

	for (i = 0; i < length; i++) {
		if (flag == 1)
			Data[i] = temp[i] / complex<double>(length);
		else
			Data[i] = temp[i];
	}
}
BUTTORD buttord(double wp, double ws, double rp, double rs)
{
	BUTTORD Temp;
	double WA;
	double WN,wn;
	double W0;
	double order;
	int ftype = 0;
	if (wp < ws)
		ftype = ftype + 1;	// low (1) or reject (3)
	else
		ftype = ftype + 2;	// high (2) or pass (4)
	double WP = wp;
	double WS = ws;
	WP = tan(pi*wp / 2);
	WS = tan(pi*ws / 2);
	if (ftype == 1)	// low
		WA = WS / WP;
	else	// high
		WA = WP / WS;
	order = ceil(log10((pow(10 , (0.1*abs(rs))) - 1) /
		(pow(10 ,(0.1000*abs(rp))) - 1) ) / (2*log10(WA)));
	W0 = WA / pow(pow(10, 0.1000*abs(rs)) - 1, 1 / (2 * abs(order)));
	if (ftype == 1) // low
		WN = W0 * WP;
	else	// high
		WN = WP / W0;



	wn = (2 / pi)*atan(WN);//数字频率
	//wn = WN;
	cout << "N:" << order << " wn:" <<wn << endl;
	printf("WN:%f\n",WN);
	Temp.N = order;
	Temp.Wn = WN;//模拟频率
	return Temp;
}

BUTTER* butter(int N, double Wn, int filter,double *az,double *bz )
{
	BUTTER* Butter = new BUTTER[N+1];
	BUTTER* Butterz = new BUTTER[N + 1];
	complex<double>* poles = new complex<double>[N];
	complex<double>* Res = new complex<double>[N + 1];
	complex<double>* Res_Save = new complex<double>[N + 1];
	double delta = 0.5;
	double *a = new double[N + 1];
	double *b = new double[N + 1];
	int count = 0, count_1 = 0;
	for (int k = 0; k < N; k++)
	{
		poles[k] = complex<double>(Wn*cos(pi*(delta + (2.0 * k + 1) / (2.0 * N))), Wn*sin(pi*(delta + (2.0 * k + 1) / (2.0 * N)))); //
	}
	Res[0] = poles[0];
	Res[1] = complex<double>(1, 0);
	for (count_1 = 0; count_1 < N - 1; count_1++)
	{
		for (count = 0; count <= count_1 + 2; count++)
		{
			if (0 == count)
			{
				Complex_Multiple(Res[count], poles[count_1 + 1],
					&Res_Save[count]);
			}

			else if ((count_1 + 2) == count)
			{
				Res_Save[count] += Res[count - 1];
			}
			else
			{
				Complex_Multiple(Res[count], poles[count_1 + 1],
					&Res_Save[count]);
				Res_Save[count] += Res[count - 1];
			}

		}


		for (count = 0; count <= N; count++)
		{
			Res[count] = Res_Save[count];
			*(a + N - count) = abs(Res[count].real());
			*(b + N - count) = 0;
			(*(Butter + N - count)).A = abs(Res[count].real());
			(*(Butter + N - count)).B = 0;

		}
		(*(Butter + N )).B = (*(Butter + N)).A;
		*(b + N) = *(a + N);

	}
	printf("a:[");
	for (count = 0; count <= N; count++)
	{
		printf("%lf ", *(a + count));
	}
	printf(" ] \nb:[");
	for (count = 0; count <= N; count++)
	{
		printf("%lf ", *(b + count));
	}
	printf(" ] \n");
	Bilinear(N,a,b,az,bz);
	for (count = 0; count <= N; count++)
	{
		 
		(*(Butterz + N - count)).A = *(az + N - count);
		(*(Butterz + N - count)).B = *(bz + N - count);

	}
	return Butterz;

}


int Complex_Multiple(complex<double> a, complex<double> b
	,complex<double> *Res)
{
	*(Res) = a * b;
	return (int)1;
}


int Bilinear(int N,//滤波器阶数
	double *as, double *bs,
	double *az, double *bz)
{

	int Count = 0, Count_1 = 0, Count_2 = 0, Count_Z = 0;
	double *Res, *Res_Save;
	Res = new double[N + 1]();
	Res_Save = new double[N + 1]();
	memset(Res, 0, sizeof(double)*(N + 1));
	memset(Res_Save, 0, sizeof(double)*(N + 1)); 
	double A = 0;
	for (Count_Z = N; Count_Z >= 0; Count_Z--)
	{
		A += *(as + Count_Z);
	}

	for (Count_Z = 0; Count_Z <= N; Count_Z++)
	{
		*(az + Count_Z) = 0;
		*(bz + Count_Z) = 0;
	}


	for (Count = 0; Count <= N; Count++)
	{
		for (Count_Z = 0; Count_Z <= N; Count_Z++)
		{
			Res[Count_Z] = 0;
			Res_Save[Count_Z] = 0;
		}
		Res_Save[0] = 1;

		for (Count_1 = 0; Count_1 < N - Count; Count_1++)//计算（1-Z^-1）^N-Count的系数,
		{												//Res_Save[]=Z^-1多项式的系数，从常数项开始
			for (Count_2 = 0; Count_2 <= Count_1 + 1; Count_2++)
			{
				if (Count_2 == 0)
				{
					Res[Count_2] += Res_Save[Count_2];
				}

				else if ((Count_2 == (Count_1 + 1)) && (Count_1 != 0))
				{
					Res[Count_2] += -Res_Save[Count_2 - 1];
				}

				else
				{
					Res[Count_2] += Res_Save[Count_2] - Res_Save[Count_2 - 1];
				}
				//printf( "Res[%d] %lf  \n" , Count_2 ,Res[Count_2]);

			}

			//printf( "Res : ");
			for (Count_Z = 0; Count_Z <= N; Count_Z++)
			{
				Res_Save[Count_Z] = Res[Count_Z];
				Res[Count_Z] = 0;
				//printf("[%d]  %lf  ", Count_Z, Res_Save[Count_Z]);
			}
			//printf(" \n");
		}

		for (Count_1 = (N - Count); Count_1 < N; Count_1++)//计算(1-Z^-1)^N-Count*（1+Z^-1）^Count的系数,
		{												//Res_Save[]=Z^-1多项式的系数，从常数项开始
			for (Count_2 = 0; Count_2 <= Count_1 + 1; Count_2++)
			{
				if (Count_2 == 0)
				{
					Res[Count_2] += Res_Save[Count_2];
				}

				else if ((Count_2 == (Count_1 + 1)) && (Count_1 != 0))
				{
					Res[Count_2] += Res_Save[Count_2 - 1];
				}

				else
				{
					Res[Count_2] += Res_Save[Count_2] + Res_Save[Count_2 - 1];
				}
				//printf("Res[%d] %lf  \n", Count_2, Res[Count_2]);

			}

			//printf("Res : ");
			for (Count_Z = 0; Count_Z <= N; Count_Z++)
			{
				Res_Save[Count_Z] = Res[Count_Z];
				Res[Count_Z] = 0;
				//printf("[%d]  %lf  ", Count_Z, Res_Save[Count_Z]);
			}
			//printf(" \n" );
		}
		//printf( "Res : ");
		for (Count_Z = 0; Count_Z <= N; Count_Z++)
		{
			//*(az + Count_Z) += pow(2, N - Count)  *  (*(as + Count)) * Res_Save[Count_Z];
			*(az + Count_Z) += (*(as + Count)) * Res_Save[Count_Z];
			*(bz + Count_Z) += (*(bs + Count)) * Res_Save[Count_Z];
			//printf( " az: %lf=%lf*%lf*%lf  " ,*(az+Count_Z), pow(2, N - Count), (*(as + Count)), Res_Save[Count_Z]);
		}
		//printf(" \n" );

	}//最外层for循环
	
	for (Count_Z = N; Count_Z >= 0; Count_Z--)
	{
		*(bz + Count_Z) = (*(bz + Count_Z)) / A;//(*(az + 0));
		*(az + Count_Z) = (*(az + Count_Z)) / A;// (*(az + 0));

	}


	//------------------------display---------------------------------//
	printf("bz =  [");
	for (Count_Z = 0; Count_Z <= N; Count_Z++)
	{
		printf("%lf ", *(bz + Count_Z));
	}
	printf(" ] \n");
	printf("az =  [");
	for (Count_Z = 0; Count_Z <= N; Count_Z++)
	{
		printf("%lf ", *(az + Count_Z));
	}
	printf(" ] \n");
	printf("--------------------------------------------------------\n");



	return (int)1;
}
int Filter(double *pdAz,//滤波器分母系数
	double *pdBz,//滤波器分子系数
	double *dDataIn,//输入数据
	int	nDataLen,//数据长度
	int nABLen,//滤波器阶数
	double *dDataOut//数据输出
) 
{
	int nALen = nABLen;
	int nBLen = nABLen;
	double *YBuf = new double[nALen+1];//
	double *XBuf = new double[nBLen+1];

	cout << "nALen:" << nALen << "nBLen:" << nBLen <<endl;
	for (int i = 0; i < nDataLen; i++)
	{
		dDataOut[i] = 0;
		for (int j = 0; j <= nABLen; j++) {
			if (i - j < 0) {
				YBuf[j] = 0;
				XBuf[j] = 0;
			}
			else {
				XBuf[j] = dDataIn[i - j];
				YBuf[j]	= dDataOut[i - j];
			}
			
		}
		dDataOut[i] = pdBz[0] * XBuf[0];
		for (int j = 1; j <= nABLen; j++) {
			dDataOut[i] += pdBz[j] * XBuf[j] - pdAz[j] * YBuf[j];// b(j)*x(n-j)-a(j)*y(n-j)
		}
		//cout << i+1<<"\t"<<dDataOut[i]  << endl;;

	}
	//cout << endl;
	delete[] YBuf;
	delete[] XBuf;

	return 0;


}

/*double FiltButter(double *pdAz,	//滤波器参数表1
	double *pdBz,	//滤波器参数表2
	int	nABLen,	//参数序列的长度
	double dDataIn,//输入数据
	double *pdBuf)	//数据缓冲区
{
	int	i;
	int nALen;
	int nBLen;
	int	nBufLen;
	double	dOut;

	if (nABLen < 1)return 0.0;
	//根据参数,自动求取序列有效长度
	nALen = nABLen;
	for (i = nABLen - 1; i; --i)
	{
		if (*(pdAz + i) != 0.0)//从最后一个系数判断是否为0
		{
			nALen = i + 1;
			break;
		}
	}
	//printf("%lf ", nALen);
	if (i == 0) nALen = 0;

	nBLen = nABLen;
	for (i = nABLen - 1; i; --i)
	{
		if (*(pdBz + i) != 0.0)
		{
			nBLen = i + 1;

			break;
		}
	}
	//printf("%lf ", nBLen);
	if (i == 0) nBLen = 0;
	//计算缓冲区有效长度
	nBufLen = nALen;
	if (nALen < nBLen)
		nBufLen = nBLen;

	//滤波: 与系数a卷乘
	dOut = (*pdAz) * dDataIn;  // a(0) * x(i)   

	for (i = 1; i < nALen; i++)	// a(i) * w(n-i),i=1toN
	{
		dOut -= *(pdAz + i) * *(pdBuf + (nBufLen - 1) - i);
	}

	//卷乘结果保存为缓冲序列的最后一个
	*(pdBuf + nBufLen - 1) = dOut;

	//滤波: 与系数b卷乘
	dOut = 0.0;
	for (i = 0; i < nBLen; i++)	// b(i) * w(n-i)
	{
		dOut += *(pdBz + i) * *(pdBuf + (nBufLen - 1) - i);
	}

	//丢弃缓冲序列中最早的一个数, 最后一个数清零
	for (i = 0; i < nBufLen - 1; i++)
	{
		*(pdBuf + i) = *(pdBuf + i + 1);
	}
	*(pdBuf + nBufLen - 1) = 0;

	//返回输出值
	return dOut;
}*/








// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
