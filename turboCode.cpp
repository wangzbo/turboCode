//��������Խ��еı����볤��Ϊ8+2*n;18+3*n,32+4*n,50+5*n,72+6*n,98+7*n,128+8*n,162+9*n,200+10*n,242+11*n��
#define Ea 0.8862
#include <iostream>
#include <fstream>                        //new added
#include <bitset>                         //new added
#include <iomanip>
#include<conio.h>
#include <ctime>
#include <cmath>
using namespace std;
//��������
void A(void);//����A����
void B(void);//����B����
void R(void);//����R����
char* juanjima(char* str);   //���뺯��
char* interleave(char* str); //��֯����
char* decoder(int n);        //���뺯��
int* suijichongpai (int a);  //������к���
int* SuiJiChongPai (int a);
double*AWGN(int C,double N0);//AWGN��������
double* rayleigh(double fm,int M,int N,double sigma);
double max0(int s00,int s01,int s,int k);//��һ��max*()
double max1(int k);                      //�ڶ���max*()
double max2(int s00,int s01,int s,int k);//������max*()
inline double modifyfunction(double a,double b);//��������
//ȫ�ֱ�������
int*array1;
int *array2;
int arr[3000],Arr[3000];
int len;      //֡��
double dB;    //�����(dB��)
double Eb_N0;//��dB����Ϊ������
double noise[3000];//�洢�����������
double ray[3000];
double decodersaver0[1000],decodersaver1[1000],decodersaver2[1000];//����ʱ��ŷ������������
double M[8],Le[1000]={0};                                                                        
double AA[1000][8],BB[1000][8],RR[1000][8][8];//�洢A,B,R������ֵ
char ee1[1000],ee2[1000],ee3[1000];
clock_t start,start_,finish,finish_;
double duration,duration_;
int malvflag,decoderflag,channelflag,iteration;

void main()
{  
	srand((unsigned)time(NULL));
	char information[1000];
	int jishu=0;
	double zhongjian[3000];
	char info[100];                          //new added
	char receive[100];                          //new added
	bitset<8> bitvect;                       //new added
	
	
	ifstream infile;//("test.txt");             //new added
	infile.open("image009.png",ios::binary);
	
	infile.seekg(0,ios::end);//�����ļ�ָ�뵽�ļ�����β��  //new added
	streampos ps = infile.tellg();  //��ȡ�ļ�ָ���λ��   //new added
    
	cout<<"File size is "<<ps<<endl;                       //new added
	infile.seekg(0,ios::beg);//�����ļ�ָ�뵽�ļ������ײ�  //new added
	
	ofstream outfile;//("1.txt");                             //new added
	outfile.open("1.png",ios::binary);
	if (outfile.fail()){//������                         //new added
		cout<<"�ļ��򿪴���"<<endl;                        //new added
		return;                                            //new added
	}                                                      //new added
	
	
	
	//for(int i=0;i<1000;i++)
	   //information[i]='!';
	//cout<<"������Ҫ���б������Ϣ���У�"<<'\n'<<endl;
	//cin>>information;
	//while (information[len+1]!='!')
	//len++;                       //ȷ������
	cout<<"Please enter Eb/N0 in dB :"<<"    ";
	cin>>dB;  cout<<'\n';
	Eb_N0=pow(10,dB/10.0);
	cout<<"Please enter the frame size :"<<"    ";
	cin>>len;  cout<<'\n';
	cout<<"Please choose punctured or unpunctured :  0: 1/2 rate  1:1/3 rate"<<"    ";
	cin>>malvflag;  cout<<'\n';
	cout<<"Please enter the decoding algorithm :  0: LOG-MAP  1:MAX-LOG-MAP"<<"    ";
	cin>>decoderflag;  cout<<'\n';
	cout<<"Please enter the maximum number of iterations for each frame :"<<"    ";
	cin>>iteration;  cout<<'\n';
	cout<<"Please choose the type of channel :0: AWGN  1:Rayleigh"<<"    ";
	cin>>channelflag;  cout<<'\n';
	   
	start_=clock();
	for(int zhenshu=1;zhenshu<=ps/100;zhenshu++) //֡��
	{
		start=clock();
		/*for(int i=0;i<len;i++)
		information[i]=rand()%2+48;*/
		infile.read(info,100);                              //new added
		for(int j=0;j<100;j++)   //��info��ÿһ�ַ���ÿһ���ش��һ���ַ�
		{
			bitvect=info[j];                                //new added
			for(int i=7;i>=0;i--)                           //new added
				information[7-i+8*j]=bitvect[i]+48;         //new added
		}
		
		char*str=information;
		char*str1=juanjima(str);
		//cout<<'\n'<<"������1���������Ϊ��"<<'\n'<<'\n'<<str1<<endl;
		char*str2=new char[1000];
		for(int i=0;i<len;i++)
			*(str2+i)=*(str1+i);            //����str1
		char*str3=juanjima(interleave(str));
		//cout<<'\n'<<"������2���������Ϊ��"<<'\n'<<'\n'<<str3<<endl;
		
		
		if(malvflag==0)
		{
			char*str4=new char[2000];
			for(i=0;i<2*len;i++)
			{if(i%2==0)
			*(str4+i)=*(str+i/2);
			else if(i%4==1)
				*(str4+i)=*(str2+(i-1)/2);
			else 
				*(str4+i)=*(str3+(i-1)/2);
			}
			*(str4+2*len)='\0';             //���и���
			//cout<<'\n'<<"1/2���ʵ�Turbo�����Ϊ��"<<'\n'<<'\n'<<str4<<endl;
			
			delete []str2;
			double guodu[2000];
			for(i=0;i<2*len;i++)
			{
				guodu[i]=double(*(str4+i)-48);
			}
			delete []str4;
			for(i=0;i<2*len;i++)
			{
				if(*(guodu+i)==0)
					*(guodu+i)=-1;           //ӳ��
			}
			
			
			if(channelflag==0)
			{
				AWGN(2*len,1/Eb_N0);                                      //biaoji
				for(i=0;i<2*len;i++)
					*(guodu+i)+=noise[i];         //ģ��AWGN�ŵ�
			}
			else
			{
				array2=SuiJiChongPai(2*len);
				
				for(i=0;i<2*len;i++)
				{
					*(zhongjian+*(array2+i)-1)=*(guodu+i);
				}
				for(i=0;i<2*len;i++)
				{
					*(guodu+i)=*(zhongjian+i);
				}
				rayleigh(41.665,10,2*len,1); 
				for(i=0;i<2*len;i++)
					*(guodu+i)*=ray[i]; 
				AWGN(2*len,1/Eb_N0); 
				
				for(i=0;i<2*len;i++)
					*(guodu+i)+=noise[i];         //ģ��AWGN�ŵ�
				
				for(i=0;i<2*len;i++)
				{
					*(zhongjian+i)=*(guodu+array2[i]-1);
				}
				for(i=0;i<2*len;i++)
				{
					*(guodu+i)=*(zhongjian+i);
				}
				
			}
			
			
			
			for(i=0;i<2*len;i++)
			{if(i%2==0)
			*(decodersaver0+i/2+1)=*(guodu+i);
			else if(i%4==1)
				*(decodersaver1+(i-1)/2+1)=*(guodu+i);
			else 
				*(decodersaver2+(i-1)/2+1)=*(guodu+i);
			}                                //�����������
			for(i=0;i<len;i++)
			{
				if(i%2==0)
					*(decodersaver2+i+1)=0;
				else
					*(decodersaver1+i+1)=0;
			}                                //��ɾ��λ����������
			
		}
		
		else
		{
			char*str4=new char[3000];
			for(i=0;i<3*len;i++)
			{if(i%3==0)
			*(str4+i)=*(str+i/3);
			else if(i%3==1)
				*(str4+i)=*(str2+(i-1)/3);
			else 
				*(str4+i)=*(str3+(i-2)/3);
			}
			*(str4+3*len)='\0';             //���и���
			//cout<<'\n'<<"1/3���ʵ�Turbo�����Ϊ��"<<'\n'<<'\n'<<str4<<endl;
			
			delete []str2;
			double guodu[3000];
			for(i=0;i<3*len;i++)
			{
				guodu[i]=double(*(str4+i)-48);
			}
			delete []str4;
			for(i=0;i<3*len;i++)
			{
				if(*(guodu+i)==0)
					*(guodu+i)=-1;           //ӳ��
			}
			
			if(channelflag==0)
			{
				AWGN(3*len,1.5/Eb_N0);                                     //biaoji
				for(i=0;i<3*len;i++)
					*(guodu+i)+=noise[i];         //ģ��AWGN�ŵ�
			}
			else
			{
				
				array2=SuiJiChongPai(3*len);
				
				for(i=0;i<3*len;i++)
				{
					*(zhongjian+*(array2+i)-1)=*(guodu+i);
				}
				for(i=0;i<3*len;i++)
				{
					*(guodu+i)=*(zhongjian+i);
				}
				rayleigh(41.665,10,3*len,1); 
				for(i=0;i<3*len;i++)
					*(guodu+i)*=ray[i]; 
				AWGN(3*len,1.5/Eb_N0); 
				
				
				for(i=0;i<3*len;i++)
					*(guodu+i)+=noise[i];         //ģ��AWGN�ŵ�
				
				for(i=0;i<3*len;i++)
				{
					*(zhongjian+i)=*(guodu+array2[i]-1);
				}
				for(i=0;i<3*len;i++)
				{
					*(guodu+i)=*(zhongjian+i);
				}
			}
			
			for(i=0;i<3*len;i++)
			{
			if(i%3==0)
			*(decodersaver0+i/3+1)=*(guodu+i);
			else if(i%3==1)
				*(decodersaver1+(i-1)/3+1)=*(guodu+i);
			else 
				*(decodersaver2+(i-2)/3+1)=*(guodu+i);
			}                                //�����������
		}
		//cout<<'\n'<<"LOG-MAP�㷨�µ�������Ϊ��"<<'\n'<<'\n'<<decoder(len)<<'\n'<<endl;
		//start=clock();
		decoder(len);
		finish=clock();
		duration=(double)(finish-start)/CLOCKS_PER_SEC;
		
		cout<</*"��"<<zhenshu<<"֡�������,��ʱ"<<*/duration<<"��."<<endl;
		//cout<<setfill('*')<<setw(161)<<endl;
		for(i=0;i<len;i++)
		{  
			if(information[i]!=ee3[i])
				jishu++;
		}
		
		
		for(j=0;j<100;j++)                       //new added
		{
			for(int i=7;i>=0;i--)                //new added
				bitvect[i]=ee3[8*j+7-i]-48;      //new added
			receive[j]=(char)bitvect.to_ulong(); //new added
			
		}                                        //new added
		outfile.write(receive,100);              //new added
		
}
outfile.close();                         //new added
finish_=clock();
duration_=(double)(finish_-start_)/CLOCKS_PER_SEC;

cout<</*'\n'<<"����֡�Ѿ��������,�ܺ�ʱ"<<*/duration_<<"��"<<"����������Ϊ"<<jishu<<'\n'<<endl;
getch();

}
/////////////////////////////////////////////////////////////////////////////

//ѭ��ϵͳ�������루RSC�������ɶ���ʽΪg1=��1101����g2=��1011��

char* juanjima(char* str)
{   
	int array[1000];
	int u[1000];
	int arra1[1000];
    for(int i=0;i<len;i++)
	{
		array[i]=*(str+i)-48;
	}
    //���Ͻ��ַ���ת��Ϊ����
	u[0]=array[0];
	u[1]=(array[1]+u[0])%2;
	u[2]=(array[2]+u[1])%2;
	for(i=3;i<len;i++)
	{
        u[i]=(array[i]+u[i-1]+u[i-3])%2;
	}
	
    arra1[0]=u[0];
	arra1[1]=u[1];
	arra1[2]=(u[2]+u[0])%2;
	for(i=3;i<len;i++)
	{
        arra1[i]=(u[i]+u[i-2]+u[i-3])%2;
	}
    //���Ͻ���ñ�������
	
	for(i=0;i<len;i++)
	{
		*(ee1+i)=char(arra1[i]+48);
	}
	*(ee1+len)='\0';
	return ee1;
	//������ת��Ϊ�ַ���
}

//Sα�����֯������
char* interleave(char* str)
{
	//�ַ���ת��Ϊ����
	int array[1000];//�洢ԭʼ����
	
	for( int i=0;i<len;i++)
	{
		array[i]=*(str+i)-48;
	}
	//���ź��������ò���S��
	/*	int s;
	if (int(sqrt(len/2))==0)
	{
	s=1;
	}
    else s=int(sqrt(len/2))-1;
	
	  array1=suijichongpai(len);//���Ȼ��һ�������
	  
		int*array11=new int[len+s];
		
		  for(i=0;i<s;i++)
		  { 
		  *(array11+i)=1000;	
		  }
		  for(i=s;i<len+s;i++)
		  { 
		  *(array11+i)=*(array1+i-s);	
		  }
		  for(i=s;i<len+s;i++)	{
		  for (int j=i-1;j>=i-s;j--)
		  {
		  while(fabs(array11[i]-array11[j])<s)
		  {
		  array1=suijichongpai(len);
		  for(i=s;i<len+s;i++)
		  { 
		  *(array11+i)=*(array1+i-s);	
		  }
		  i=s;
		  j=i-1;
		  }
		  }
		  }
	delete[]array11;*/
	//���ϼ�����������Ƿ�����������ֱ���ó�����������Sα�漴���У�����array1��
	
    array1=suijichongpai(len);
	
    //���µõ��������������
    for(i=0;i<len;i++)
	{
		*(ee2+*(array1+i)-1)=char(*(array+i)+48);
	}
	*(ee2+len)='\0';
	return ee2;
}

int* suijichongpai (int a)
{
	
/*int c,p;
for(unsigned int i=0;i<a;i++)
{
arr[i]=i+1;
}
for(i=0;i<a-1;i++)
{

  p=rand()%(a-(i+1))+(i+1);
  c=arr[i];
  arr[i]=arr[p];
  arr[p]=c;
  }
  
	//ע�͵����������*/
	
	int zu=int(sqrt(a/2));
	int arr1[100];
	int arr2[1000];//
	int c,p;
	for(int i=0;i<zu;i++)
	{
		arr1[i]=i+1;
	}
	//srand((unsigned)time(NULL));
    for(i=0;i<zu-1;i++)
	{
		
		p=rand()%(zu-(i+1))+(i+1);
		c=arr1[i];
		arr1[i]=arr1[p];
		arr1[p]=c;
	}
	//���϶Է����������������
	for(i=0;i<zu;i++)
		
	{
		
		for(int j=0;j<a/zu;j++)
		{
			arr[j]=j+1+i*(a/zu);
		}
		//srand((unsigned)time(NULL));
		for(j=0;j<a/zu-1;j++)
		{
			
			p=rand()%(a/zu-(j+1))+(j+1);
			c=arr[j];
			arr[j]=arr[p];
			arr[p]=c;
			
		}
		
		for(j=0;j<a/zu;j++)
		{
			arr2[j+i*(a/zu)]=arr[j];
		}
	}
	//���϶�ÿ�������������������У�����arr2��
	
    for(i=0;i<zu;i++)
	{
		for(int j=0;j<a/zu;j++)
		{
			arr[arr1[i]-1+j*zu]=arr2[j+i*(a/zu)];
		}
	}
	
	//���ϰ���ӳ���ϵ����arr2�ж�Ӧ��arr�еõ��������*/
	
	
	
	return arr;//�洢���ɵ������
}


//���²���C������(0,N0/2)��̬�ֲ��ĸ�˹���������ֵ�������noise��
double*   AWGN(int C,double sigma2)   
{ 
	for (int j=0;j<C;j++)
	{
		double   n=0;    
		for(int   i=0;i<120;i++)   
		{   
			n+=(double)rand()/RAND_MAX;   
		}   
		
		n=(n-60)/sqrt(10); //��׼��   
		
		noise[j]=sqrt(sigma2)*n; 
	} 
	return   noise; 
}   

//����ģ��:

char* decoder(int n)
{  
	double L[1000],Le1[1000],x[1000],y[1000],W[16];
	char ee4[1000];
	int error_num;
	for(int q=0;q<iteration;q++)//����10��
	{  //��ʼ�Ա�����1��������Ϣ
		error_num=0;
		R();
		A();//��������AA
		B();//��������BB
		for(int k=1;k<=n;k++)
		{
			W[0]=AA[k-1][0]+RR[k][0][4]+BB[k][4];
			W[1]=AA[k-1][1]+RR[k][1][0]+BB[k][0];
			W[2]=AA[k-1][2]+RR[k][2][5]+BB[k][5];
			W[3]=AA[k-1][3]+RR[k][3][1]+BB[k][1];
			W[4]=AA[k-1][4]+RR[k][4][2]+BB[k][2];
			W[5]=AA[k-1][5]+RR[k][5][6]+BB[k][6];
			W[6]=AA[k-1][6]+RR[k][6][3]+BB[k][3];
			W[7]=AA[k-1][7]+RR[k][7][7]+BB[k][7];
			for(int j=1;j<8;j++)
			{ if(W[0]<W[j])
			W[0]=W[j]+modifyfunction(W[0],W[j]);
			else 
				W[0]=W[0]+modifyfunction(W[0],W[j]);
			}
			W[8]=AA[k-1][0]+RR[k][0][0]+BB[k][0];
			W[9]=AA[k-1][1]+RR[k][1][4]+BB[k][4];
			W[10]=AA[k-1][2]+RR[k][2][1]+BB[k][1];
			W[11]=AA[k-1][3]+RR[k][3][5]+BB[k][5];
			W[12]=AA[k-1][4]+RR[k][4][6]+BB[k][6];
			W[13]=AA[k-1][5]+RR[k][5][2]+BB[k][2];
			W[14]=AA[k-1][6]+RR[k][6][7]+BB[k][7];
			W[15]=AA[k-1][7]+RR[k][7][3]+BB[k][3];
			for(j=9;j<16;j++)
			{ if(W[8]<W[j])
			W[8]=W[j]+modifyfunction(W[8],W[j]);
			else 
				W[8]=W[8]+modifyfunction(W[8],W[j]);
			}
			L[k]=W[0]-W[8];
			Le[k]=L[k]-Le[k]-decodersaver0[k]*Eb_N0*((6-2*malvflag)/3.0);
		}
		
		//����Le
		
		for(int i=1;i<=len;i++)
		{
			Le1[i]=Le[i];
		}
		
		//���¶�decodersaver0��Le���н�֯
		for(i=0;i<len;i++)
		{
			*(x+*(array1+i)-1)=*(decodersaver0+i+1);
			*(y+*(array1+i)-1)=*(Le+i+1);
		}
		for(i=0;i<len;i++)
		{
			*(decodersaver0+i+1)=*(x+i);
			*(Le+i+1)=*(y+i);
		}
		
		//����Ϊ��ʹ������2�ܹ�����R����������decodersaver2��decodersaver1
		for(i=1;i<=len;i++)
		{
			*(x+i-1)=*(decodersaver1+i);
			*(decodersaver1+i)=*(decodersaver2+i);
			*(decodersaver2+i)=*(x+i-1);
		}
		
		//��ʼ�Ա�����2��������Ϣ
		R();
		A();
		B();
		for(k=1;k<=n;k++)
		{
			W[0]=AA[k-1][0]+RR[k][0][4]+BB[k][4];
			W[1]=AA[k-1][1]+RR[k][1][0]+BB[k][0];
			W[2]=AA[k-1][2]+RR[k][2][5]+BB[k][5];
			W[3]=AA[k-1][3]+RR[k][3][1]+BB[k][1];
			W[4]=AA[k-1][4]+RR[k][4][2]+BB[k][2];
			W[5]=AA[k-1][5]+RR[k][5][6]+BB[k][6];
			W[6]=AA[k-1][6]+RR[k][6][3]+BB[k][3];
			W[7]=AA[k-1][7]+RR[k][7][7]+BB[k][7];
			for(int j=1;j<8;j++)
			{ 
				if(W[0]<W[j])
					W[0]=W[j]+modifyfunction(W[0],W[j]);
				else 
					W[0]=W[0]+modifyfunction(W[0],W[j]);
			}
			W[8]=AA[k-1][0]+RR[k][0][0]+BB[k][0];
			W[9]=AA[k-1][1]+RR[k][1][4]+BB[k][4];
			W[10]=AA[k-1][2]+RR[k][2][1]+BB[k][1];
			W[11]=AA[k-1][3]+RR[k][3][5]+BB[k][5];
			W[12]=AA[k-1][4]+RR[k][4][6]+BB[k][6];
			W[13]=AA[k-1][5]+RR[k][5][2]+BB[k][2];
			W[14]=AA[k-1][6]+RR[k][6][7]+BB[k][7];
			W[15]=AA[k-1][7]+RR[k][7][3]+BB[k][3];
			for(j=9;j<16;j++)
			{ if(W[8]<W[j])
			W[8]=W[j]+modifyfunction(W[8],W[j]);
			else 
				W[8]=W[8]+modifyfunction(W[8],W[j]);
			}
			L[k]=W[0]-W[8];
			Le[k]=L[k]-Le[k]-decodersaver0[k]*Eb_N0*((6-2*malvflag)/3.0);
		}
		//����decodersaver2��decodersaver1
		for(i=1;i<=len;i++)
		{
			*(x+i-1)=*(decodersaver1+i);
			*(decodersaver1+i)=*(decodersaver2+i);
			*(decodersaver2+i)=*(x+i-1);
		}
		
		//��decodersaver0��Le�⽻֯
		for(i=0;i<len;i++)
		{
			*(x+i)=*(decodersaver0+array1[i]);
			*(y+i)=*(Le+array1[i]);
		}
		for(i=0;i<len;i++)
		{
			*(decodersaver0+i+1)=*(x+i);
			*(Le+i+1)=*(y+i);
		}
		
		
		for(k=1;k<=n;k++)
		{
			L[k]=decodersaver0[k]*Eb_N0*((6-2*malvflag)/3.0)+Le[k]+Le1[k];
			if(L[k]>0)
				*(ee3+k-1)='1';
			else 
				*(ee3+k-1)='0';
		}
        
		
		
		for(i=0;i<len;i++)
		{  
			if(ee4[i]!=ee3[i])
				error_num++;
			ee4[i]=ee3[i];
		}
		
		if (error_num==0)
			break;
}

//for(int k=1;k<=n;k++)
//{
//L[k]=decodersaver0[k]*Eb_N0*((6-2*malvflag)/3.0)+Le[k]+Le1[k];
//if(L[k]>0)
//*(ee3+k-1)='1';
//else 
//*(ee3+k-1)='0';
//}
//*(ee3+n)='\0';
return ee3;

}

/*double R(unsigned int k,int s0,int s)  //����R��k��s0��s����ֵ
{   
if((s0==0)&(s==4))
return 0.5*Le[k]+(decodersaver0[k]+decodersaver1[k])/(1/Eb_N0);
else if((s0==1)&(s==0))
return 0.5*Le[k]+(decodersaver0[k]+decodersaver1[k])/(1/Eb_N0);
else if((s0==6)&(s==3))
return 0.5*Le[k]+(decodersaver0[k]+decodersaver1[k])/(1/Eb_N0);
else if((s0==7)&(s==7))
return 0.5*Le[k]+(decodersaver0[k]+decodersaver1[k])/(1/Eb_N0);
else if((s0==0)&(s==0))
return -0.5*Le[k]-(decodersaver0[k]+decodersaver1[k])/(1/Eb_N0);
else if((s0==1)&(s==4))
return -0.5*Le[k]-(decodersaver0[k]+decodersaver1[k])/(1/Eb_N0);
else if((s0==6)&(s==7))
return -0.5*Le[k]-(decodersaver0[k]+decodersaver1[k])/(1/Eb_N0);
else if((s0==7)&(s==3))
return -0.5*Le[k]-(decodersaver0[k]+decodersaver1[k])/(1/Eb_N0);
else if((s0==2)&(s==5))
return 0.5*Le[k]+(decodersaver0[k]-decodersaver1[k])/(1/Eb_N0);
else if((s0==3)&(s==1))
return 0.5*Le[k]+(decodersaver0[k]-decodersaver1[k])/(1/Eb_N0);
else if((s0==4)&(s==2))
return 0.5*Le[k]+(decodersaver0[k]-decodersaver1[k])/(1/Eb_N0);
else if((s0==5)&(s==6))
return 0.5*Le[k]+(decodersaver0[k]-decodersaver1[k])/(1/Eb_N0);
else 
return -0.5*Le[k]+(decodersaver1[k]-decodersaver0[k])/(1/Eb_N0);
}*/

void R(void)  //����R��k��s0��s����ֵ                             //biaoji 
{   
	for(int i=1;i<=len;i++)
	{
		RR[i][0][4]=0.5*Le[i]+((Ea-1)*channelflag+1)*(decodersaver0[i]+decodersaver1[i])/((malvflag+2)/(2*Eb_N0));
		RR[i][1][0]=0.5*Le[i]+((Ea-1)*channelflag+1)*(decodersaver0[i]+decodersaver1[i])/((malvflag+2)/(2*Eb_N0));
		RR[i][6][3]=0.5*Le[i]+((Ea-1)*channelflag+1)*(decodersaver0[i]+decodersaver1[i])/((malvflag+2)/(2*Eb_N0));
		RR[i][7][7]=0.5*Le[i]+((Ea-1)*channelflag+1)*(decodersaver0[i]+decodersaver1[i])/((malvflag+2)/(2*Eb_N0));
		
		RR[i][0][0]=-0.5*Le[i]-((Ea-1)*channelflag+1)*(decodersaver0[i]+decodersaver1[i])/((malvflag+2)/(2*Eb_N0));
		RR[i][1][4]=-0.5*Le[i]-((Ea-1)*channelflag+1)*(decodersaver0[i]+decodersaver1[i])/((malvflag+2)/(2*Eb_N0));
		RR[i][6][7]=-0.5*Le[i]-((Ea-1)*channelflag+1)*(decodersaver0[i]+decodersaver1[i])/((malvflag+2)/(2*Eb_N0));
		RR[i][7][3]=-0.5*Le[i]-((Ea-1)*channelflag+1)*(decodersaver0[i]+decodersaver1[i])/((malvflag+2)/(2*Eb_N0));
		
		RR[i][2][5]=0.5*Le[i]+((Ea-1)*channelflag+1)*(decodersaver0[i]-decodersaver1[i])/((malvflag+2)/(2*Eb_N0));
		RR[i][3][1]=0.5*Le[i]+((Ea-1)*channelflag+1)*(decodersaver0[i]-decodersaver1[i])/((malvflag+2)/(2*Eb_N0));
		RR[i][4][2]=0.5*Le[i]+((Ea-1)*channelflag+1)*(decodersaver0[i]-decodersaver1[i])/((malvflag+2)/(2*Eb_N0));
		RR[i][5][6]=0.5*Le[i]+((Ea-1)*channelflag+1)*(decodersaver0[i]-decodersaver1[i])/((malvflag+2)/(2*Eb_N0));
		
		RR[i][2][1]=-0.5*Le[i]+((Ea-1)*channelflag+1)*(decodersaver1[i]-decodersaver0[i])/((malvflag+2)/(2*Eb_N0));
		RR[i][3][5]=-0.5*Le[i]+((Ea-1)*channelflag+1)*(decodersaver1[i]-decodersaver0[i])/((malvflag+2)/(2*Eb_N0));
		RR[i][4][6]=-0.5*Le[i]+((Ea-1)*channelflag+1)*(decodersaver1[i]-decodersaver0[i])/((malvflag+2)/(2*Eb_N0));
		RR[i][5][2]=-0.5*Le[i]+((Ea-1)*channelflag+1)*(decodersaver1[i]-decodersaver0[i])/((malvflag+2)/(2*Eb_N0));
	}
	
}


void A(void)  //����A��k��s����ֵ
{
    AA[0][0]=0;
	AA[0][1]=-10000;
	AA[0][2]=-10000;
	AA[0][3]=-10000;
	AA[0][4]=-10000;
	AA[0][5]=-10000;
	AA[0][6]=-10000;
	AA[0][7]=-10000;
	for(int i=1;i<=len;i++)
	{
		AA[i][0]=max0(0,1,0,i)-max1(i);
		AA[i][1]=max0(2,3,1,i)-max1(i);
		AA[i][2]=max0(4,5,2,i)-max1(i);
		AA[i][3]=max0(6,7,3,i)-max1(i);
		AA[i][4]=max0(0,1,4,i)-max1(i);
		AA[i][5]=max0(2,3,5,i)-max1(i);
		AA[i][6]=max0(4,5,6,i)-max1(i);
		AA[i][7]=max0(6,7,7,i)-max1(i);
	}	
}

void B(void)  //����B��k��s0����ֵ
{
    for(int i=0;i<8;i++)
		BB[len][i]=-2.08;
	for(i=len-1;i>=0;i--)
	{
		BB[i][0]=max2(0,0,4,i+1)-max1(i+1);
		BB[i][1]=max2(1,4,0,i+1)-max1(i+1);
		BB[i][2]=max2(2,5,1,i+1)-max1(i+1);
		BB[i][3]=max2(3,5,1,i+1)-max1(i+1);
		BB[i][4]=max2(4,2,6,i+1)-max1(i+1);
		BB[i][5]=max2(5,6,2,i+1)-max1(i+1);
		BB[i][6]=max2(6,3,7,i+1)-max1(i+1);
		BB[i][7]=max2(7,7,3,i+1)-max1(i+1);
	}
}

double max0(int s00,int s01,int s,int k)//��һȷ����s��������һ״̬�����ֵ
{
	if(AA[k-1][s00]+RR[k][s00][s]>AA[k-1][s01]+RR[k][s01][s])
		return AA[k-1][s00]+RR[k][s00][s]+modifyfunction(AA[k-1][s00]+RR[k][s00][s],AA[k-1][s01]+RR[k][s01][s]);
	else
		return AA[k-1][s01]+RR[k][s01][s]+modifyfunction(AA[k-1][s00]+RR[k][s00][s],AA[k-1][s01]+RR[k][s01][s]);
}

double max1(int k)//������s��������������һ״̬�����ֵ
{  
	M[0]=max0(0,1,0,k);
    M[1]=max0(2,3,1,k);
	M[2]=max0(4,5,2,k);
	M[3]=max0(6,7,3,k);
	M[4]=max0(0,1,4,k);
	M[5]=max0(2,3,5,k);
	M[6]=max0(4,5,6,k);
	M[7]=max0(6,7,7,k);
	for(int i=1;i<8;i++)
	{
		if(M[0]<M[i])
			M[0]=M[i]+modifyfunction(M[0],M[i]);
        else
            M[0]=M[0]+modifyfunction(M[0],M[i]);
	}
	return M[0];
}

double max2(int s0,int s,int s1,int k)//��һȷ����s0��������һ״̬�����ֵ
{  
    if(BB[k][s]+RR[k][s0][s]>BB[k][s1]+RR[k][s0][s1])
		return BB[k][s]+RR[k][s0][s]+modifyfunction(BB[k][s]+RR[k][s0][s],BB[k][s1]+RR[k][s0][s1]);
	else
		return BB[k][s1]+RR[k][s0][s1]+modifyfunction(BB[k][s]+RR[k][s0][s],BB[k][s1]+RR[k][s0][s1]);
}


inline double modifyfunction(double a,double b)
{
	switch(decoderflag)
	{
	case 0 :
		if ((a=fabs(a-b))<0.182)
			return 0.65;
		else if (a<0.384)
			return 0.563;
		else if (a<0.613)
			return 0.476;
		else if (a<0.882)
			return 0.39;
		else if (a<1.215)
			return 0.303;
		else if (a<1.666)
			return 0.22;
		else if (a<2.4)
			return 0.13;
		else if (a<5)
			return 0.0433;
		else 
			return 0;
		
		
			/*if ((a=fabs(a-b))<0.625)
			return 0.55;
			else if (a<1.25)
			return 0.33;
			else if (a<1.875)
			return 0.19;
			else if (a<2.5)
			return 0.106;
			else if (a<3.125)
			return 0.058;
			else if (a<3.75)
			return 0.032;
			else if (a<4.375)
			return 0.017;
			else if (a<5)
			return 0.0092;
			else 
		return 0;*/
		break;
	case 1:
		return 0;
		break;
	}	
	//return log(1+exp(-fabs(a-b)));
}




double* rayleigh(double fm,int M,int N,double sigma)
{
	double *t=new double[N];
	double x[3000]={0} ;
	double y[3000]={0} ;
	double alpha[100],ph1[100],ph2[100];
	for(int i=0;i<N;i++)
		t[i]=0.00024*i;
	double C=sqrt(2/(double)M)*sigma;
	double W=6.28*fm;
	for(i=0;i<M;i++)
	{
		alpha[i]=(6.28*i-3.14+(rand()/(double)RAND_MAX*6.28-3.14))/(double)(4*M);
		ph1[i]=rand()/(double)RAND_MAX*6.28-3.14;
		ph2[i]=rand()/(double)RAND_MAX*6.28-3.14;
	}
	
	for(i=0;i<M;i++)
	{
		for(int j=0;j<N;j++) 
		{	
			x[j]=x[j]+C*cos(W*t[j]*cos(alpha[i])+ph1[i]);
			y[j]=y[j]+C*cos(W*t[j]*sin(alpha[i])+ph2[i]);
		}
	}
	for(int j=0;j<N;j++)
		ray[j]=sqrt(pow(x[j],2)+pow(y[j],2))/sqrt(2.0);
	delete []t;
	return ray;
}



int* SuiJiChongPai (int a)
{
	
    int c,p;
    for(unsigned int i=0;i<a;i++)
	{
		Arr[i]=i+1;
	}
    for(i=0;i<a-1;i++)
	{
		
		p=rand()%(a-(i+1))+(i+1);
		c=Arr[i];
		Arr[i]=Arr[p];
		Arr[p]=c;
	}
	return Arr;
}
