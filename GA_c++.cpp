#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>
#include <list>
#include <algorithm>
#include <fstream>
using namespace std;

ofstream out;
short int APWakeTime[4][4] = {{0, 0, 0, 0}, {0, 80, 80, 80}, {0, 107, 80, 80}, {0, 160, 133, 107}}; // ��1��ʼ

class AP
{
public:
	string Number;	  // �ɻ����
	string TS;		  // ��ͨ��
	short int APKind; // �ͺ�
	int delay;		  // �ִ�ʱ��
	AP() {}
	AP(string n, short int k, string s, int d) : Number(n), APKind(k), TS(s), delay(d) {}
	friend ostream &operator<<(ostream &os, const AP &ap)
	{
		os << ap.Number << ":" << ap.APKind << ":" << ap.TS << ":" << ap.delay;
		return os;
	}
	friend bool operator==(const AP &A, const AP &B)
	{
		return A.Number == B.Number;
	}
};

class individual
{
public:
	vector<AP *> chromosome; // Ⱦɫ��
	double fitness;			 // ��Ӧ��
	double val;				 // ����ֵ

	void print()
	{
		typename vector<AP *>::iterator it = this->chromosome.begin();
		cout << "������" << chromosome.size() << "������";
		for (; it != this->chromosome.end(); it++)
			cout << *(*it) << "|";
		cout << endl;
	}
	void write()
	{
		typename vector<AP *>::iterator it = this->chromosome.begin();
		cout << "������" << chromosome.size() << "������\n";
		for (; it != this->chromosome.end(); it++)
			out << **it << " ";
		out << endl;
		out << val << endl;
	}

	bool isAfrontofB() {}
};

static void mySwap(AP *&A, AP *&B)
{
	AP *temp = A;
	A = B;
	B = temp;
}

// GA�Ŵ��㷨��ģ��
class GAalg
{

public:
	// �㷨��ʼ��
	GAalg(int sg, int sc, int mg, double pc, double pm) : sizeof_generation(sg), sizeof_chrom(sc), max_generation(mg), p_crossover(pc), p_mutation(pm)
	{

		// �ļ���ʼ��
		out.open("data.txt");
		cout << "Start GA algorithm...\n";
		generationbox.clear();
		// ���´���
		cur_generation = 0;
		// ���������
		srand(time(NULL));
	}
	// �㷨����
	virtual ~GAalg()
	{
		cout << "End GA algorithm...\n";
		cout << endl
			 << endl
			 << "��ʷ���Ÿ���Ϊ:\n";
		bestone.print();
		generationbox.clear();
		system("pause");
	}
	// ��ȡ��������
	void ReadAPData()
	{
		ifstream in("input.txt");
		while (!in.eof())
		{
			AP ap;
			in >> ap.Number >> ap.APKind >> ap.TS >> ap.delay;
			apbox.push_back(ap);
		}
		in.close();
	}
	// ��ӡ���и���
	virtual void print()
	{
		cout << "����" << sizeof_generation << "�����壬ÿ������������£�\n";
		typename vector<individual>::iterator it = generationbox.begin();
		for (; it != generationbox.end(); it++)
			(*it).print();
	}
	// ��ʼ����Ⱥ�еĸ���
	virtual void InitIndividual(individual &my)
	{
		my.chromosome.clear();
		// ʹ��ϴ���㷨
		for (auto &i : apbox)
		{
			my.chromosome.push_back(&i);
		}
		for (int i = my.chromosome.size() - 1; i > 0; i--)
		{
			mySwap(my.chromosome[i], my.chromosome[rand() % i]);
		}
		return;
	}
	// ��Ⱥ��ʼ������ʼ������
	virtual void GenerateInit()
	{
		ReadAPData();
		for (int i = 0; i < sizeof_generation; i++)
		{

			cout << "���ڳ�ʼ����" << i + 1 << "������" << endl;
			individual obj;
			InitIndividual(obj);
			generationbox.push_back(obj);
		}

		cur_generation = 0;
		return;
	}
	// ������һ��
	virtual void NextPopulation()
	{
		cur_generation++;
		cout << endl
			 << endl
			 << "���ڴ����" << cur_generation << "��...\n";
		SelectOperator(); // ѡ������
		// print();
		// system("pause");
		cout << "���ڽ���" << endl;
		CrossoverOperator(); // ��������
		cout << "���ڱ���" << endl;
		MutationOperator(); // ��������
	}
	// ���۸���
	virtual void ValuatePopulation()
	{
		CalculateFitnessValue();   // ������Ӧ��
		FindBestWorstIndividual(); // ������Ѹ����������
	}
	// ������Ѹ����������
	virtual void FindBestWorstIndividual()
	{
		// ��ʼ����ǰ��Ⱥ��ֵ
		cur_bestone = *generationbox.begin();
		cur_worstone = *generationbox.begin();
		// ��������iterator
		typename vector<individual>::iterator it = generationbox.begin();
		int i = 0;
		// Ѱ�ҵ�ǰ������
		for (it = generationbox.begin(); it != generationbox.end(); it++)
		{
			if (it->fitness >= cur_bestone.fitness)
			{
				cur_bestone = *it;
				best_index = i;
			}
			else if (it->fitness <= cur_worstone.fitness)
			{
				cur_worstone = *it;
				worst_index = i;
			}
			i++;
		}
		// ������ʷ������
		if (cur_generation == 0)
		{
			bestone = cur_bestone;
			worstone = cur_worstone;
		}
		else
		{
			if (cur_bestone.fitness >= bestone.fitness)
				bestone = cur_bestone;
			if (cur_worstone.fitness <= worstone.fitness)
				worstone = cur_worstone;
		}
		return;
	}
	// ���䷽���������Դ�λ��pos֮���x��y��Ⱦɫ��
	void cross(int x, int y, int pos, int end)
	{
		cout << "���䣺" << x << " " << y << " " << pos << " " << end << endl;
		// ��ʱ����
		typename vector<individual>::iterator itx = generationbox.begin();
		typename vector<individual>::iterator ity = generationbox.begin();
		individual valx;
		individual valy;
		// Ѱ���ݴ�
		for (int i = 0; i != x; i++)
			itx++;
		for (int i = 0; i != y; i++)
			ity++;
		valx = *itx;
		valy = *ity;
		map<AP *, AP *> mapa, mapb;
		// cout<<"����"<<x+1<<"��"<<y+1<<"�ڻ���λ��"<<pos+1<<"��ʼ����ֱ��λ��"<<end+1<<endl;
		// ��ʼ����x
		typename vector<AP *>::iterator from, targe;
		from = valy.chromosome.begin();
		targe = (*itx).chromosome.begin();
		for (int i = 0; i != pos; i++)
		{
			from++;
			targe++;
		}
		for (int p = pos; p <= end; p++)
		{
			mapa[*targe] = *from;
			*targe = *from;
			from++;
			targe++;
		}
		// ��ʼ����y
		from = valx.chromosome.begin();
		targe = (*ity).chromosome.begin();
		for (int i = 0; i != pos; i++)
		{
			from++;
			targe++;
		}
		for (int p = pos; p <= end; p++)
		{
			mapb[*targe] = *from;
			*targe = *from;
			from++;
			targe++;
		}
		// �����ظ�Ԫ��
		for (auto pair : mapa)
		{
			if (mapa.find(pair.second) != mapa.end())
			{
				mapa[pair.first] = mapa[pair.second];
			}
		}
		for (auto pair : mapb)
		{
			if (mapb.find(pair.second) != mapb.end())
			{
				mapb[pair.first] = mapb[pair.second];
			}
		}
		// �м�ӳ��
		from = (*itx).chromosome.begin();
		targe = (*ity).chromosome.begin();
		for (int i = 0; i != pos; i++)
		{
			if (mapb.find(*from) != mapb.end())
				*from = mapb[*from];
			if (mapa.find(*targe) != mapa.end())
				*targe = mapa[*targe];
			from++;
			targe++;
		}
		for (int p = pos; p <= end; p++)
		{
			from++;
			targe++;
		}
		for (int e = end + 1; e < valx.chromosome.size(); e++)
		{
			if (mapb.find(*from) != mapb.end())
				*from = mapb[*from];
			if (mapa.find(*targe) != mapa.end())
				*targe = mapa[*targe];
			from++;
			targe++;
		}
		// �ж��Ƿ�Ϊ��Ч��
		if (!IsOk(*itx))
		{
			*itx = valx;
			*ity = valy;
			cout << "fail\n";
			return;
		}
		if (!IsOk(*ity))
		{
			*itx = valx;
			*ity = valy;
			cout << "fail\n";
			return;
		}
		return;
	}
	// ѡ������*���̶�ѡ�㷨ROUTE WHEEL SELECTION*
	virtual void SelectOperator()
	{
		// cout<<"ѡ������*���̶�ѡ�㷨�����ڸ�����Ⱥ...\n";
		// cout<<"����ǰ��\n";
		// this->print();
		int index, i = 0;
		// ��������Լ�����Ӧ��
		double p, sum = 0.0;
		// ����Ⱥtemp
		vector<individual> newgenerationbox;
		// ǰ׺����Ӧ��
		double *addfit = new double[sizeof_generation];
		// ��������iterator
		typename vector<individual>::iterator it = generationbox.begin();
		// ������Ӧ��
		for (it = generationbox.begin(); it != generationbox.end(); it++)
			sum += it->fitness;
		// �������������
		for (it = generationbox.begin(); it != generationbox.end(); it++)
			addfit[i++] = it->fitness / sum;

		// �����������������ǰ׺��
		i = 1;
		for (it = ++generationbox.begin(); it != generationbox.end(); it++, i++)
			addfit[i] = addfit[i - 1] + addfit[i];
		/*//���ǰ׺������
		for(int j=0;j<sizeof_generation;j++)
			cout<<addfit[j]<<" ";
		cout<<"\n"<<p<<"\n";*/
		// ģ�����̶�ѡ�㷨
		for (i = 0; i < sizeof_generation; i++)
		{
			p = rand() % 1000 / 1000.0;
			index = 0;
			while (p >= addfit[index])
				index++;
			it = generationbox.begin();
			for (int j = 0; j < index; j++)
				it++;
			newgenerationbox.push_back(*it);
		}
		// ���µ�����Ⱥ
		generationbox.clear();
		it = newgenerationbox.begin();
		for (i = 0; i < sizeof_generation; i++)
		{
			generationbox.push_back(*it);
			it++;
		}
		// cout<<"���º�\n";
		// this->print();
		return;
	}
	// ��������
	virtual void CrossoverOperator()
	{
		// cout<<"�������ӣ����ڸ�����Ⱥ...\n";
		// ����������ʣ��ж��Ƿ���
		bool *pi = new bool[sizeof_generation];
		for (int i = 0; i < sizeof_generation; i++)
		{
			double p = rand() % 1000 / 1000.0;
			pi[i] = (p < p_crossover);
		}
		// ����һ���������
		list<int> cl;
		for (int i = 0; i < sizeof_generation; i++)
			if (pi[i])
				cl.push_back(i);
		// �ж��Ƿ�Ϊż��
		if (cl.size() % 2 != 0)
			cl.pop_back();
		// ��ʼ����
		while (!cl.empty())
		{
			int x = cl.front();
			cl.pop_front();
			int y = cl.front();
			cl.pop_front();
			double p1 = rand() % sizeof_chrom, p2 = rand() % sizeof_chrom;
			if (p2 < p1)
				swap(p1, p2);
			cross(x, y, p1, p2);
		}
		return;
	}
	// ��������
	virtual void MutationOperator()
	{
		// cout<<"�������ӣ����ڸ�����Ⱥ...\n";
		int i, j;
		double p;
		// ��ȡÿһ������
		typename vector<individual>::iterator iti;
		for (iti = generationbox.begin(); iti != generationbox.end(); iti++)
		{
			// ��ȡ�����ÿһ������
			individual temp = *iti;
			i = rand() % temp.chromosome.size(), j = rand() % temp.chromosome.size();
			if (i > j)
				swap(i, j);
			auto itx = temp.chromosome.begin();
			while (i > 0)
			{
				++itx;
				--i;
				--j;
			}
			auto ity = itx;
			while (j > 0)
			{
				++ity;
				--j;
			}
			mySwap(*itx, *ity);
			// �Ƿ񱣴����
			if (IsOk(temp))
				*iti = temp;
		}
		return;
	}
	// �������չʾ����
	virtual void OutputPrint()
	{
		out << bestone.fitness << endl;
		/*
				cout<<"����չʾ��"<<cur_generation<<"��...\n";
				cout<<"�������Ÿ���Ϊ:\n";
				cur_bestone.print();
				out<<bestone.val<<endl;
				cout<<"����ֵΪ��"<<cur_bestone.val<<endl;
				cout<<"��Ӧ��Ϊ��"<<cur_bestone.fitness<<endl;
				cout<<"����ƽ����Ӧ�ȣ�";
				double sum=0;
				for(vector<individual<T> >::iterator it=generationbox.begin();it!=generationbox.end();it++)
					sum+=it->fitness;
				cout<<sum/(double)sizeof_generation<<endl;
				cout<<"��ʷ���Ÿ���Ϊ:\n";
				bestone.print();
				cout<<"����ֵΪ��"<<bestone.val<<endl;
				cout<<"��Ӧ��Ϊ��"<<bestone.fitness<<endl;*/
	}
	// ������մ�
	virtual void ans()
	{
		out << endl
			<< endl
			<< "��ʷ���Ÿ���Ϊ:\n";
		bestone.write();
		out << "���ô�����" << cur_generation << endl;
		out << "��Ӧֵ��" << bestone.fitness << endl;
		out << "����ֵ��" << bestone.val << endl;
		out.close();
	}

	// ���㺯��ֵ
	virtual void CalculateVal(individual &my)
	{
		// ����ֻ�����˸��ݻ��Ͳ�ͬβ��ʱ��õ������Ŷ�ʱ��
		auto i = my.chromosome.begin();
		auto j = i;
		++j;
		double time = 0;
		// my.print();
		// cout<<"time=";
		while (j != my.chromosome.end())
		{
			// cout<<"+"<<APWakeTime[(*i)->APKind][(*j)->APKind];
			time += APWakeTime[(*i)->APKind][(*j)->APKind];
			++i;
			++j;
		}
		// cout<<"="<<time<<endl;
		// system("pause");
		my.val = time;
	}
	// ������Ӧ��
	virtual void CalculateFitnessValue()
	{
		double totaltime = apbox.size() * 160;
		auto it = generationbox.begin();
		// ����ֵ
		for (it = generationbox.begin(); it != generationbox.end(); it++)
		{
			CalculateVal(*it);
			it->fitness = pow(1 / (it->val / totaltime), 2);
		}
	}
	// ��ֹ����
	virtual bool IsEnd()
	{
		return cur_generation >= max_generation;
	}
	// �Ƿ�Ϊ��Ч��
	virtual bool IsOk(individual cur) { return 1; }

public:
	// �ɻ�����
	vector<AP> apbox;
	// Ⱦɫ������
	vector<individual> generationbox;
	// ��Ѹ���
	individual cur_bestone;
	// ������
	individual cur_worstone;
	// ��ʷ��Ѹ���
	individual bestone;
	// ��ʷ������
	individual worstone;

	// ���ݳ�Ա
	int sizeof_generation; // ��Ⱥ��С
	int sizeof_chrom;	   // �������
	int max_generation;	   // ����Ŵ�������
	int cur_generation;	   // ��ǰ������
	double p_crossover;	   // �������
	double p_mutation;	   // �������
	int best_index;		   // ��Ѹ���λ��
	int worst_index;	   // ������λ��
};

int main()
{
	int sg, sc, mg; // ��Ⱥ���������������������
	double pc, pm;	// ������ʣ��������
	/*printf("��ʼ��ȫ�ֱ���:\n");
	printf("\t��Ⱥ��С:");
	cin>>sg;
	if((sg%2) != 0){
		printf("��Ⱥ��С������Ϊż��\n");
		sg++;
	};
	printf("\t�����������");
	cin>>mg;
	printf("\t�����ʣ�");
	cin>>pc;
	printf("\t�����ʣ�");
	cin>>pm;
	printf("\t�������������");
	cin>>_max_sum;*/
	GAalg a(30, 4, 8, 0.60, 0.005); // ��Ⱥ�����������������������������ʣ��������
	a.GenerateInit();
	// a.print();
	cout << "��ʼ�Ŵ��㷨��" << endl;
	a.ValuatePopulation();
	while (!a.IsEnd())
	{
		a.NextPopulation();
		a.ValuatePopulation();
		a.OutputPrint();
	}
	a.ans();
}