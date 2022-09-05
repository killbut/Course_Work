#include<cmath>
#include<fstream>
#include<iostream>
#include<vector>
#include<algorithm>
#include <set>

using namespace std;

int numVert; // ���������� ������ �� ��������� �������
int numTriag; // ���������� �������������
int numFinitElem; // ���������� �������� ��������� �� ���� ��������� �������
int timeIntervals; // ���������� ��������� ����������
vector<double> timeT; // ������ �������
int countTime = 0; // ������� ��������� ����

vector<int> ig; // ������ ���������� ������� � ����������� �������
vector<int> jg; // ������� ���������� ������� � ����������� �������


struct Point
{
	int glNumber; // ���������� �����
	double x, y; // ���������� ����
};

struct Range
{
	vector<double> diffusionCoff; // �������� ������������ ������
	double gammaCoff; // �������� ������������ �����
	vector<int> functionValue; // �������� �������
	vector<int> functionForAppr; // �������� ������� ���������� �����������
	double FormulaForFunction(double x, double y, int numTriag) // �������� ������� f
	{
		switch (functionValue[numTriag])
		{
		case 0:
			return 0; // ���� 1, 2, 3
		case 1:
			return -6; // ���� 4
		case 2:
			return -8; // ���� 5
		case 3:
			return -18 * x; // ���� 6
		case 4:
			return -6 * y; // ���� 7
		case 5:
			return -12 * x * x; // ���� 8
		case 6:
			return -48 * y * y; // ���� 9
		case 7:
			return (3 * y * y) - (2 * timeT[countTime]); // ���� 11 (u = y * y * t)
		case 8:
			return (4 * y * y * y * timeT[countTime]) - (6 * timeT[countTime] * timeT[countTime] * y); // ���� 12 (u = y * y * y * t * t)
		case 9:
			return 7; // ���� 10 (u = 7t)

		case 10: // ���� 13 (u = y * y * y * y * t * t * t)
			return (-12 * y * y * timeT[countTime] * timeT[countTime] * timeT[countTime]) + (3 * y * y * y * y * timeT[countTime] * timeT[countTime]);
		case 11:
			return -2;
		case 12:
			return y * y * y - (6 * y * timeT[countTime]);
		case 13:
			return 49 * y ;

		default:
			break;
		}

	}

	double FormulaForApprox(double x, double y, int numTriag) // �������� ������� u
	{
		switch (functionForAppr[numTriag])
		{
		case 0:
			return 0;
		case 1:
			return 7 * timeT[countTime]; // ���� 10 (u = 7 * t)
		case 2:
			return y * y * timeT[countTime]; // ���� 11 (u = y * y * t)
		case 3:
			return y * y * y * timeT[countTime] * timeT[countTime]; // ���� 12 (u = y * y * y * t * t)
		case 4:
			return 7; // ���� 1 (u = 7)
		case 5:
			return x; // ���� 2 (u = x)
		case 6:
			return y; // ���� 3 (u = y) 
		case 7:
			return x * x; // ���� 4 (u = x * x)
		case 8:
			return y * y; // ���� 5 (u = y * y)
		case 9:
			return x * x * x; // ���� 6 (u = x * x * x)
		case 10:
			return y * y * y; // ���� 7 (u = y * y * y)
		case 11:
			return x * x * x * x; // ���� 8 (u = x * x * x * x)
		case 12:
			return y * y * y * y; // ���� 9 (u = y * y * y * y)
		case 13:
			return y * y * y * y * timeT[countTime] * timeT[countTime]; // ���� 12 (u = y * y * y * y * t * t * t)
		case 14:
			return 7 * y; // 
		case 15:
			return y * y * y * timeT[countTime]; // 
		case 16:
			return 49 * y * timeT[countTime]; // 

		default:
			break;
		}

	}
};

// �������� �� L
struct gradL
{
	double comp1; // ������ �������� ��� �������� ���������
	double comp2; // ������ �������� ��� �������� ���������
};

// ������ ������� �������
struct BoundaryCondition1
{
	Point node; // ���� �� ������� �������� ������� �������
	int formulaNum; // ����� �������� ������� �������� �������

	double FormulaForCondition1() // ����� �������� ������� �������� �������
	{
		switch (formulaNum)
		{
		case 0:
			return 7; // ��� ����� 1
		case 1:
			return 1; // ��� ������ 2, 3, 4, 5, 6, 7, 8, 9
		case 2:
			return 7; // ��� ������ 2, 3
		case 3:
			return 49; // ��� ������ 4, 5
		case 4:
			return 343; // ��� ������ 6, 7
		case 5:
			return 2401; // ��� ������ 8, 9
		case 6:
			return timeT[countTime]; // ��� ����� 11
		case 7:
			return 49 * timeT[countTime]; // ��� ����� 11
		case 8:
			return timeT[countTime] * timeT[countTime]; // ��� ����� 12
		case 9:
			return 343 * timeT[countTime] * timeT[countTime]; // ��� ����� 12
		case 10:
			return 7 * timeT[countTime]; // ��� ����� 10
		case 11:
			return timeT[countTime] * timeT[countTime] * timeT[countTime]; // ��� ����� 13
		case 12:
			return 2401 * timeT[countTime] * timeT[countTime] * timeT[countTime]; // ��� ����� 13
		case 13:
			return 13;
		case 14:
			return 13 * 7;
		case 15:
			return 13 * 13 * timeT[countTime];
		case 16:
			return 49 * 13 * timeT[countTime];

		default:
			break;
		}
	}
};

// ������ ������� �������
struct BoundaryCondition2
{
	Point node[4]; // ����� �� ������� �������� ������� �������
	int formulaNum; // ����� �������� ������� �������� �������

	double FormulaForCondition2() // ����� �������� ������� �������� �������
	{
		switch (formulaNum)
		{
		case(0):
			return 0;
		default:
			break;
		}
	}
};


// ������� ����� ��������� ������ �������������, ���������� ������� ���� ����� ������������, ���������� ��������� ������� � �������� � �������
// ���������: 
//            point - ���������� ����� ������� ������������
//            triangle - ���������� ������ ����� ������������
//            range - ��������� ��������� �������
void Input(vector<Point>& point, vector<vector<int>>& triangle, Range& range)
{
	// �������� �� �����
	ifstream finXY("coordxy.txt");
	finXY >> numVert >> numFinitElem;
	point.resize(numVert);
	for (int i = 0; i < numVert; i++)
	{
		finXY >> point[i].x >> point[i].y;
		point[i].glNumber = i;
	}
	finXY.close();

	// �������� � �������������
	ifstream finParam("m.txt");
	finParam >> numTriag;
	triangle.resize(numTriag, vector<int>(10, 0));
	for (int i = 0; i < numTriag; i++)
		for (int j = 0; j < 10; j++)
			finParam >> triangle[i][j];

	finParam.close();

	// �������� � ���������� ����
	ifstream fieldIn("field.txt");
	fieldIn >> range.gammaCoff;
	range.diffusionCoff.resize(numTriag);
	range.functionValue.resize(numTriag);
	range.functionForAppr.resize(numTriag);
	for (int i = 0; i < numTriag; i++)
		fieldIn >> range.diffusionCoff[i];
	for (int i = 0; i < numTriag; i++)
		fieldIn >> range.functionValue[i];
	for (int i = 0; i < numTriag; i++)
		fieldIn >> range.functionForAppr[i];
	fieldIn.close();

	// �������� � �������
	ifstream finTime("time.txt");
	double timeEnd, oneInter;
	finTime >> timeIntervals;
	timeT.resize(timeIntervals + 1);
	finTime >> timeT[0] >> timeEnd;
	oneInter = (timeEnd - timeT[0]) / timeIntervals;
	for (int i = 1; i < timeIntervals + 1; i++)
		timeT[i] = timeT[i - 1] + oneInter;
	finTime.close();
}

// ������� �������������� ��������� �������� �������, �.�. ig � jg
// ���������: 
//            triangle - ���������� ������ ����� ������������
//            n - ���������� ���� ����� �� ��������� �������
void Portrait(vector<vector<int>> triangle, int n)
{
	vector<set<int>> connection(n);

	for (auto& t : triangle)
	{
		for (int i = 1; i < 10; i++)
			for (int j = 0; j < i; j++)
			{
				int a = t[i];
				int b = t[j];
				if (a < b)
					swap(a, b);

				connection[a].insert(b);
			}
	}


	ig.resize(n + 1);
	int* IA = &ig[0];
	IA[0] = IA[1] = 1;

	for (int i = 2; i <= n; i++)
	{
		int col = IA[i - 1];
		IA[i] = col + connection[i - 1].size();
	}

	jg.resize(IA[n]);
	int* JA = &jg[0];
	for (int i = 1, k = 0; i < n; i++)
		for (int j : connection[i])
		{
			JA[k] = j;
			k++;
		}

}

// ������� �������� ���������� (��� ������ ��������� �������)
// ���������: 
//            num - �����, ���������� � ���������
//            res - ��������� ���������� � ���������
void Factorial(int num, double& res)
{
	res = 1;
	for (int i = 2; i <= num; i++)
		res *= i;
}

// ������� �������� ���������� (��� ������ ��������� �������)
// ���������: 
//            sign - ���� ���������
//            val - ������, � ������� �������� L1, L2, L3
//            res - ��������� �������� ���������
//            detD - ������������ �������
void CalculationIntegralL(double sign, vector<int> val, double& res, double detD)
{

	vector<double> rs(4, 0);
	Factorial(val[0], rs[0]);
	Factorial(val[1], rs[1]);
	Factorial(val[2], rs[2]);
	Factorial(val[0] + val[1] + val[2] + 2, rs[3]);
	res = (rs[0] * rs[1] * rs[2]) / (rs[3]) * sign * detD;
}

// ������� �������� ���������� ������������ gladL(n) �� gladL(k)
// ���������: 
//            g1 - gladL(n)
//            g2 - gladL(k)
double Scalar(gradL g1, gradL g2)
{
	return g1.comp1 * g2.comp1 + g1.comp2 * g2.comp2;
}

// ������� ���������� ��������� �������
// ���������: 
//            point - ���������� ����� ������� ������������
//            triangle - ���������� ������ ����� ������������
//            localMatrixHard - ��������� ������� ���������
//            localMatrixMass - ��������� ������� ����
//            numLocalTriag - ����� ���������� ������������
void BuildLocalMatrix(vector<Point> point, vector<vector<int>> triangle, vector<vector<double>>& localMatrixHard, vector<vector<double>>& localMatrixMass, int numLocalTriag)
{
	double detD = ((point[1].x - point[0].x) * (point[2].y - point[0].y)) - ((point[2].x - point[0].x) * (point[1].y - point[0].y));


	double alpha[3][3] =
	{
		{(point[1].x * point[2].y - point[2].x * point[1].y) / detD, (point[1].y - point[2].y) / detD, (point[2].x - point[1].x) / detD},
		{(point[2].x * point[0].y - point[0].x * point[2].y) / detD, (point[2].y - point[0].y) / detD, (point[0].x - point[2].x) / detD},
		{(point[0].x * point[1].y - point[1].x * point[0].y) / detD, (point[0].y - point[1].y) / detD, (point[1].x - point[0].x) / detD},
	}; // ������� �������� ����

	detD = abs(detD); // ������� ������������ � ������� ��� ���������� ��������

	// ������� �������� gradL1, gradL2, gradL3
	vector<gradL> gradLN(3);
	for (int i = 0; i < 3; i++)
	{
		gradLN[i].comp1 = alpha[i][1];
		gradLN[i].comp2 = alpha[i][2];
	}

	// ������ ������ �������� L-��������� � gradPsi
	// ������������ ��� ���������� ������� ���������
	vector<vector<vector<int>>> powL =
	{
		{ { 2, 0, 0 }, { 1, 0, 0 }, { 0, 0, 0 } }, // gradPsi1
		{ { 0, 2, 0 }, { 0, 1, 0 }, { 0, 0, 0 } }, // gradPsi2
		{ { 0, 0, 2 }, { 0, 0, 1 }, { 0, 0, 0 } }, // gradPsi3

		{ { 1, 1, 0 }, { 2, 0, 0 }, { 0, 1, 0 }, { 1, 0, 0 } }, // gradPsi4
		{ { 1, 1, 0 }, { 0, 2, 0 }, { 1, 0, 0 }, { 0, 1, 0 } }, // gradPsi5

		{ { 0, 1, 1 }, { 0, 2, 0 }, { 0, 0, 1 }, { 0, 1, 0 } }, // gradPsi6
		{ { 0, 1, 1 }, { 0, 0, 2 }, { 0, 1, 0 }, { 0, 0, 1 } }, // gradPsi7

		{ { 1, 0, 1 }, { 0, 0, 2 }, { 1, 0, 0 }, { 0, 0, 1 } }, // gradPsi8
		{ { 1, 0, 1 }, { 2, 0, 0 }, { 0, 0, 1 }, { 1, 0, 0 } }, // gradPsi9

		{ { 0, 1, 1 }, { 1, 0, 1 }, { 1, 1, 0 } } // gradPsi10
	};

	// ������ ������ �������� L-��������� � Psi
	// ������������ ��� ���������� ������� ����
	vector<vector<vector<int>>> basis =
	{
		{ { 3, 0, 0}, { 2, 0, 0}, { 1, 0, 0} }, // Psi1
		{ { 0, 3, 0}, { 0, 2, 0}, { 0, 1, 0} }, // Psi2
		{ { 0, 0, 3}, { 0, 0, 2}, { 0, 0, 1} }, // Psi3

		{ { 2, 1, 0}, { 1, 1, 0} }, // Psi4
		{ { 1, 2, 0}, { 1, 1, 0} }, // Psi5

		{ { 0, 2, 1}, { 0, 1, 1} }, // Psi6
		{ { 0, 1, 2}, { 0, 1, 1} }, // Psi7

		{ { 1, 0, 2}, { 1, 0, 1} }, // Psi8
		{ { 2, 0, 1}, { 1, 0, 1} }, // Psi9

		{ { 1, 1, 1} } // Psi10
	};

	// ������ ������ ������������� � gradPsi
	// ������������ ��� ������� ���������
	vector<vector<gradL>> gradPsi =
	{
		{{gradLN[0]}, {gradLN[0]}, {gradLN[0]}}, // gradPsi1
		{{gradLN[1]}, {gradLN[1]}, {gradLN[1]}}, // gradPsi2
		{{gradLN[2]}, {gradLN[2]}, {gradLN[2]}}, // gradPsi3

		{{gradLN[0]}, {gradLN[1]}, {gradLN[0]}, {gradLN[1]}}, // gradPsi4
		{{gradLN[1]}, {gradLN[0]}, {gradLN[1]}, {gradLN[0]}}, // gradPsi5

		{{gradLN[1]}, {gradLN[2]}, {gradLN[1]}, {gradLN[2]}}, // gradPsi6
		{{gradLN[2]}, {gradLN[1]}, {gradLN[2]}, {gradLN[1]}}, // gradPsi7

		{{gradLN[2]}, {gradLN[0]}, {gradLN[2]}, {gradLN[0]}}, // gradPsi8
		{{gradLN[0]}, {gradLN[2]}, {gradLN[0]}, {gradLN[2]}}, // gradPsi9

		{{gradLN[0]}, {gradLN[1]}, {gradLN[2]}} // gradPsi10
	};

	vector<vector<double>> num =
	{
		{{13.5}, {-9}, {1}}, // gradPsi1
		{{13.5}, {-9}, {1}}, // gradPsi2
		{{13.5}, {-9}, {1}}, // gradPsi3

		{{27}, {13.5}, {-4.5}, {-4.5}}, // gradPsi4
		{{27}, {13.5}, {-4.5}, {-4.5}}, // gradPsi5

		{{27}, {13.5}, {-4.5}, {-4.5}}, // gradPsi6
		{{27}, {13.5}, {-4.5}, {-4.5}}, // gradPsi7

		{{27}, {13.5}, {-4.5}, {-4.5}}, // gradPsi8
		{{27}, {13.5}, {-4.5}, {-4.5}}, // gradPsi9

		{{27}, {27}, {27}} // gradPsi10
	};

	// ������ ������ ������������� � Psi
	// ������������ � ������� ����
	vector<vector<double>> basisCoff =
	{
		{ { 4.5}, { -4.5}, { 1} }, // Psi1
		{ { 4.5}, { -4.5}, { 1} }, // Psi2
		{ { 4.5}, { -4.5}, { 1} }, // Psi3
		{ { 13.5}, { -4.5} }, // Psi4
		{ { 13.5}, { -4.5} }, // Psi5
		{ { 13.5}, { -4.5} }, // Psi6
		{ { 13.5}, { -4.5} }, // Psi7
		{ { 13.5}, { -4.5} }, // Psi8
		{ { 13.5}, { -4.5} }, // Psi9
		{ { 27} } // Psi10
	};

	vector<int> sumV(3, 0);
	double valIntegralDop, valIntegral = 0; // ���������� ��� �������� ���������� 
	// ������� ������� ���������
	for (int i = 0; i < 10; i++)
	{
		for (int j = i; j < 10; j++)
		{
			for (int v = 0; v < powL[i].size(); v++)
			{
				for (int k = 0; k < powL[j].size(); k++)
				{
					sumV[0] = powL[i][v][0] + powL[j][k][0];
					sumV[1] = powL[i][v][1] + powL[j][k][1];
					sumV[2] = powL[i][v][2] + powL[j][k][2];
					CalculationIntegralL(num[i][v] * num[j][k], sumV, valIntegralDop, detD);
					valIntegral += valIntegralDop * Scalar(gradPsi[i][v], gradPsi[j][k]);
				}
			}
			localMatrixHard[i][j] = valIntegral;
			localMatrixHard[j][i] = localMatrixHard[i][j];

			valIntegral = 0;
		}
	}

	// ������� ������� ����
	for (int i = 0; i < 10; i++)
	{
		for (int j = i; j < 10; j++)
		{
			for (int v = 0; v < basis[i].size(); v++)
			{
				for (int k = 0; k < basis[j].size(); k++)
				{
					sumV[0] = basis[i][v][0] + basis[j][k][0];
					sumV[1] = basis[i][v][1] + basis[j][k][1];
					sumV[2] = basis[i][v][2] + basis[j][k][2];
					CalculationIntegralL(basisCoff[i][v] * basisCoff[j][k], sumV, valIntegralDop, detD);
					valIntegral += valIntegralDop;
				}
			}
			localMatrixMass[i][j] = valIntegral;
			localMatrixMass[j][i] = localMatrixMass[i][j];

			valIntegral = 0;
		}
	}

}

// ��������� ���������� �������
// ���������: 
//            node - ���������� ����� ������� ������������
//            localVector - ��������� ������ ������ �����
//            range - �������� ��������� �������
void BuildLocalVector(vector<Point> node, vector<double>& localVector, Range range, int numTriag)
{
	// ����������
	// L1 = p1, L2 = p2, L3 = 1 - L2 - L1
	double p1[21] =
	{ 0.0451890097844, 0.0451890097844, 0.9096219804312, 0.7475124727339, 0.2220631655373, 0.7475124727339,
		0.2220631655373, 0.0304243617288, 0.0304243617288, 0.1369912012649, 0.6447187277637, 0.1369912012649, 0.2182900709714,
		0.2182900709714, 0.6447187277637, 0.0369603304334, 0.4815198347833, 0.4815198347833, 0.4036039798179, 0.4036039798179,
		0.1927920403641
	};
	double p2[21] =
	{ 0.0451890097844, 0.9096219804312, 0.0451890097844, 0.0304243617288, 0.0304243617288, 0.2220631655373,
		0.7475124727339, 0.7475124727339, 0.2220631655373, 0.2182900709714, 0.2182900709714, 0.6447187277637, 0.6447187277637,
		0.1369912012649, 0.1369912012649, 0.4815198347833, 0.0369603304334, 0.4815198347833, 0.1927920403641, 0.4036039798179,
		0.4036039798179
	};
	// ����
	double w[21] =
	{ 0.0519871420646, 0.0519871420646, 0.0519871420646, 0.0707034101784, 0.0707034101784, 0.0707034101784,
		0.0707034101784, 0.0707034101784, 0.0707034101784, 0.0909390760952, 0.0909390760952, 0.0909390760952, 0.0909390760952,
		0.0909390760952, 0.0909390760952, 0.1032344051380, 0.1032344051380, 0.1032344051380, 0.1881601469167, 0.1881601469167,
		0.1881601469167
	};

	// ��������� ���������� ��� ������ ����������
	double LQ[3];
	// �������� x � y � ������ ����������
	double xQ = 0, yQ = 0;
	// ��������������� ����������
	double s = 0;
	double detD = abs((node[1].x - node[0].x) * (node[2].y - node[0].y) - (node[2].x - node[0].x) * (node[1].y - node[0].y));

	for (int i = 0; i < 21; i++)
	{
		LQ[0] = p1[i];
		LQ[1] = p2[i];
		LQ[2] = 1 - LQ[1] - LQ[0];
		xQ = node[0].x * LQ[0] + node[1].x * LQ[1] + node[2].x * LQ[2];
		yQ = node[0].y * LQ[0] + node[1].y * LQ[1] + node[2].y * LQ[2];

		if (range.functionForAppr[numTriag] == 0)
			s = 0.25 * w[i] * range.FormulaForFunction(xQ, yQ, numTriag);
		else
			s = 0.25 * w[i] * range.FormulaForApprox(xQ, yQ, numTriag);


		// �������� ��� � ������ ����������
		vector<double> valuePsi =
		{
			{ 4.5 * LQ[0] * LQ[0] * LQ[0] - 4.5 * LQ[0] * LQ[0] + LQ[0] }, // Psi1
			{ 4.5 * LQ[1] * LQ[1] * LQ[1] - 4.5 * LQ[1] * LQ[1] + LQ[1] }, // Psi2
			{ 4.5 * LQ[2] * LQ[2] * LQ[2] - 4.5 * LQ[2] * LQ[2] + LQ[2] }, // Psi3

			{ 13.5 * LQ[0] * LQ[0] * LQ[1] - 4.5 * LQ[0] * LQ[1] }, // Psi4
			{ 13.5 * LQ[0] * LQ[1] * LQ[1] - 4.5 * LQ[0] * LQ[1] }, // Psi5
			{ 13.5 * LQ[1] * LQ[1] * LQ[2] - 4.5 * LQ[1] * LQ[2] }, // Psi6

			{ 13.5 * LQ[1] * LQ[2] * LQ[2] - 4.5 * LQ[1] * LQ[2] }, // Psi7
			{ 13.5 * LQ[2] * LQ[2] * LQ[0] - 4.5 * LQ[2] * LQ[0] }, // Psi8
			{ 13.5 * LQ[2] * LQ[0] * LQ[0] - 4.5 * LQ[2] * LQ[0] }, // Psi9

			{ 27 * LQ[0] * LQ[1] * LQ[2] } // Psi10
		};
		for (int j = 0; j < 10; j++)
			localVector[j] += valuePsi[j] * s * detD;
	}

}

// ������� ��� ��������� ��������� ������� � ����������
// ���������: 
//            num - ���������� ����� � ����������� �������
//            vectorNode - ������ ���� ���������� ������� ����� � ����� ������������
//            n - ���������� ����� ����
bool isInGlobal(int num, vector<int> vectorNode, int& n)
{
	for (int i = 0; i < vectorNode.size(); i++)
	{
		if (vectorNode[i] == num)
		{
			n = i;
			return true;
		}
	}
	return false;
}

// ������� ��� ��������� ��������� ������� � ���������� ������� ������ ����� � ���������� ������� � ���������� ������ ��������������
// ���������: 
//            gglGlobalMatrixHard - ������ ����������� ���������� ������� ���������
//            gglGlobalMatrixMass - ������ ����������� ���������� ������� ����
//            diGlobalMatrixHard - ��������� ���������� ������� ���������
//            diGlobalMatrixMass - ��������� ���������� ������� ����
//            globalVector - ���������� ������ ������ �����
//            point - ���������� ����� ������� ������������
//            triangle - ���������� ������ ����� ������������
//            range - ��������� ��������� �������
void BuildGlobalMatrixAndVector(vector<double>& gglGlobalMatrixHard, vector<double>& gglGlobalMatrixMass, vector<double>& diGlobalMatrixHard, vector<double>& diGlobalMatrixMass, vector<double>& globalVector, vector<Point> point, vector<vector<int>> triangle, Range range)
{
	vector<double> localVector(10, 0);
	vector<vector<double>> localMatrixHard(10, vector<double>(10, 0));
	vector<vector<double>> localMatrixMass(10, vector<double>(10, 0));

	/*vector<vector<double>> globalMatrixMassPL(16, vector<double>(16, 0));
	vector<vector<double>> globalMatrixMassRAZ(16, vector<double>(16, 0));*/


	// numLocalTriag - ���-�� �������������
	for (int numLocalTriag = 0; numLocalTriag < numTriag; numLocalTriag++)
	{
		// node - ���� ��� ��������� ������� � �������(3 �������)
		vector<Point> node(3);
		for (int i = 0; i < 3; i++)
		{
			node[i].glNumber = triangle[numLocalTriag][i];
			node[i].x = point[node[i].glNumber].x;
			node[i].y = point[node[i].glNumber].y;
		}

		// �������� ��������� �������
		BuildLocalMatrix(node, triangle, localMatrixHard, localMatrixMass, numLocalTriag);
		for (int i = 0; i < localMatrixHard.size(); i++)
			for (int j = 0; j < localMatrixHard[i].size(); j++)
				localMatrixHard[i][j] *= range.diffusionCoff[numLocalTriag];

		// �������� ��������� ������
		BuildLocalVector(node, localVector, range, numLocalTriag);

		vector<int> globalNode(10, 0);
		// ������ ��� �������� � ���������� ������
		for (int i = 0; i < 10; i++)
			globalNode[i] = triangle[numLocalTriag][i];


		for (int i = 0; i < ig.size() - 1; i++)
		{
			int x, y;
			int num = 0;
			if (isInGlobal(i, globalNode, y))
			{
				for (int j = ig[i] - 1; num < ig[i + 1] - ig[i]; j++, num++)
				{
					if (isInGlobal(jg[j], globalNode, x))
					{
						gglGlobalMatrixHard[j] += localMatrixHard[x][y];
						gglGlobalMatrixMass[j] += localMatrixMass[x][y];
					}
				}
			}
		}

		for (int i = 0; i < 10; i++)
		{
			diGlobalMatrixHard[globalNode[i]] += localMatrixHard[i][i];
			diGlobalMatrixMass[globalNode[i]] += localMatrixMass[i][i];
		}


		// ������� ��������� ������ � ����������
		for (int i = 0; i < 10; i++)
		{
			globalVector[globalNode[i]] += localVector[i];
		}


		// ��������� ��������� ������
		for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 10; j++)
			{
				localMatrixHard[i][j] = 0;
				localMatrixMass[i][j] = 0;
			}
		}


		// ��������� ���������� �������
		for (int i = 0; i < 10; i++)
			localVector[i] = 0;

	}

}

// ������� ������ ������� �������
// ���������: 
//            nodes - ���������� ����� ������� ������������, �� ������� ����������� ������� �������
//            condition1 - ������� ������� ������� ����
//            condition2 - ������� ������� ������� ����
void InputBoundaryCondition(vector<Point>& nodes, vector<BoundaryCondition1>& condition1, vector<BoundaryCondition2>& condition2)
{
	ifstream inputBoundary1("boundaryCondition1.txt");
	int numBound1;
	inputBoundary1 >> numBound1;
	condition1.resize(numBound1);
	for (int i = 0; i < numBound1; i++)
		inputBoundary1 >> condition1[i].node.glNumber >> condition1[i].formulaNum;
	inputBoundary1.close();
	ifstream inputBoundary2("boundaryCondition2.txt");
	int numBound2;
	inputBoundary2 >> numBound2;
	condition2.resize(numBound2);
	for (int i = 0; i < numBound2; i++)
	{
		for (int j = 0; j < 4; j++)
			inputBoundary2 >> condition2[i].node[j].glNumber;
		inputBoundary2 >> condition2[i].formulaNum;
		for (int j = 0; j < 2; j++)
		{
			condition2[i].node[j].x = nodes[condition2[i].node[j].glNumber].x;
			condition2[i].node[j].y = nodes[condition2[i].node[j].glNumber].y;
		}
	}
	inputBoundary2.close();
}

// ������� ����� ������� ������� 1 ����
// ���������: 
//            condition - ������� ������� ������� ����
//            di - ��������� �������
//            globalVector - ������ ������ �����
void EntryBoundaryCondition1(vector<BoundaryCondition1> condition, vector<double>& di, vector<double>& globalVector)
{

	for (int i = 0; i < condition.size(); i++)
	{
		di[condition[i].node.glNumber] = 1.0E50;
		globalVector[condition[i].node.glNumber] = 1.0E50 * condition[i].FormulaForCondition1();
	}
}

// ���� ������� ������� 2 ����
// ������� ����� ������� ������� 2 ����
// ���������: 
//            condition - ������� ������� ������� ����
//            globalVector - ������ ������ �����
void EntryBoundaryCondition2(vector<BoundaryCondition2> condition, vector<double>& globalVector)
{
	double h = 0;
	for (int i = 0; i < condition.size(); i++)
	{
		h = condition[i].FormulaForCondition2() * sqrt((condition[i].node[1].x - condition[i].node[0].x) * (condition[i].node[1].x - condition[i].node[0].x)
			+ (condition[i].node[1].y - condition[i].node[0].y) * (condition[i].node[1].y - condition[i].node[0].y));
		globalVector[condition[i].node[0].glNumber] += h * 0.125;
		globalVector[condition[i].node[1].glNumber] += h * 0.125;
		globalVector[condition[i].node[2].glNumber] += h * 0.375;
		globalVector[condition[i].node[3].glNumber] += h * 0.375;
	}
}

// ������� ���������� ������������ ��������
// ���������: 
//            x - ������ ������
//            y - ������ ������
double ScalarProduct(vector<double>& x, vector<double>& y)
{
	double scalarProduct = 0;
	for (int i = 0; i < x.size(); i++)
		scalarProduct += x[i] * y[i];
	return scalarProduct;
}

// ������� ������������ ������� �� ������
// ���������: 
//            di - ��������� �������
//            ggl - ������ ����������� �������
//            x - ������
vector<double> MatrixVectorProduct(vector<double> di, vector<double> ggl, vector<double>& x)
{
	vector<double> resultVector(x.size());
	for (int i = 0; i < di.size(); i++)
	{
		resultVector[i] = di[i] * x[i];
		for (int j = ig[i] - 1; j < ig[i + 1] - 1; j++)
		{
			resultVector[i] += ggl[j] * x[jg[j]];
			resultVector[jg[j]] += ggl[j] * x[i];
		}
	}
	return resultVector;
}

// ������� ��� ������� ���� ������� ����������� ����������
// ���������: 
//            di - ��������� �������
//            ggl - ������ ����������� �������
//            B - ������ ������ ����� (��������� ����������������)
void MCG(vector<double> di, vector<double> ggl, vector<double>& B)
{
	int n = di.size();
	vector<double> composition_Az(n);
	vector<double> z = B;
	vector<double> r = B;
	vector<double> approximation(n);
	double scalarProduct_rr, normVectorF = 0, normVectorR = 0;
	double relativeDiscrepancy;
	double a, b, eps = 1.0E-35;
	int iteration = 0, maxIteration = 100;
	scalarProduct_rr = ScalarProduct(r, r);
	for (int i = 0; i < B.size(); i++)
		normVectorF += B[i] * B[i];
	normVectorF = sqrt(normVectorF);
	do
	{
		composition_Az = MatrixVectorProduct(di, ggl, z);
		a = scalarProduct_rr / ScalarProduct(composition_Az, z);
		for (int i = 0; i < n; i++)
			approximation[i] += a * z[i];
		for (int i = 0; i < n; i++)
			r[i] -= a * composition_Az[i];
		b = scalarProduct_rr;
		scalarProduct_rr = ScalarProduct(r, r);
		if (abs(b) < 1.0e-60)
			break;
		b = scalarProduct_rr / b;
		for (int i = 0; i < n; i++)
			z[i] = z[i] * b + r[i];
		for (int i = 0; i < r.size(); i++)
			normVectorR += r[i] * r[i];
		normVectorR = sqrt(normVectorR);
		relativeDiscrepancy = normVectorR / normVectorF;
		iteration++;
	} while (iteration < maxIteration);
	B = approximation;
}

// ������� ��� ���������� Q1 c ������� ������� ����������� �����
// ���������: 
//            resultQ0 - ��������� ����������� Q0
//            point - ���������� ����� ������� ������������
//            triangle - ���������� ������ ����� ������������
//            range - ��������� ��������� �������
//            condition1 - ������� ������� ������� ����
//            condition2 - ������� ������� ������� ����
// ���������: 
//            ������ Q1
vector<double> ImplicitTwolayerScheme(vector<double> resultQ0, vector<Point> point, vector<vector<int>> triangle, Range range, vector<BoundaryCondition1> cond1, vector<BoundaryCondition2> cond2)
{

	// ������ �����������
	vector<double> gglGlobalMatrix(ig[numFinitElem] - 1, 0);
	vector<double> gglGlobalMatrixHard(ig[numFinitElem] - 1, 0);
	vector<double> gglGlobalMatrixMass(ig[numFinitElem] - 1, 0);

	// ������������ �������
	vector<double> diGlobalMatrix(numFinitElem, 0);
	vector<double> diGlobalMatrixHard(numFinitElem, 0);
	vector<double> diGlobalMatrixMass(numFinitElem, 0);

	vector<double> globalVector(numFinitElem, 0);
	vector<double> resultQ1(numFinitElem, 0);
	vector<double> temp(numFinitElem, 0);

	double deltaTime = timeT[countTime] - timeT[countTime - 1]; // �������� deltaT

	// ��������� ���������� ������ ��������� � ���� � ����������� �������
	BuildGlobalMatrixAndVector(gglGlobalMatrixHard, gglGlobalMatrixMass, diGlobalMatrixHard, diGlobalMatrixMass, globalVector, point, triangle, range);

	// A = G + (1/deltaT * M)
	for (int i = 0; i < diGlobalMatrix.size(); i++)
		diGlobalMatrix[i] = (diGlobalMatrixHard[i]) + (diGlobalMatrixMass[i] * (range.gammaCoff / deltaTime));
	for (int i = 0; i < gglGlobalMatrix.size(); i++)
		gglGlobalMatrix[i] = (gglGlobalMatrixHard[i]) + (gglGlobalMatrixMass[i] * (range.gammaCoff / deltaTime));

	// M * q0 
	for (int i = 0; i < numFinitElem; i++)
	{
		temp[i] = diGlobalMatrixMass[i] * resultQ0[i];
		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			temp[i] += gglGlobalMatrixMass[j - 1] * resultQ0[jg[j - 1]];
			temp[jg[j - 1]] += gglGlobalMatrixMass[j - 1] * resultQ0[i];
		}
		
	}

	// d = b + 1/deltaT * M * q0
	for (int i = 0; i < resultQ1.size(); i++)
	{
		temp[i] *= (range.gammaCoff / deltaTime);
		resultQ1[i] = globalVector[i] + temp[i]; 
	}

	// ���� ������� �������
	EntryBoundaryCondition2(cond2, resultQ1);
	EntryBoundaryCondition1(cond1, diGlobalMatrix, resultQ1);

	MCG(diGlobalMatrix, gglGlobalMatrix, resultQ1);

	// �������� q1
	return resultQ1;

}

// ������� ��� ���������� Q2 c ������� ������� ����������� �����
// ���������: 
//            resultQ0 - ��������� ����������� Q0
//            resultQ1 - �������� Q1, ����������� � ������� ������� ����������� ������
//            point - ���������� ����� ������� ������������
//            triangle - ���������� ������ ����� ������������
//            range - ��������� ��������� �������
//            condition1 - ������� ������� ������� ����
//            condition2 - ������� ������� ������� ����
// ���������: 
//            ������ Q2
vector<double> ImplicitThreelayerScheme(vector<double> resultQ0, vector<double> resultQ1, vector<Point> point, vector<vector<int>> triangle, Range range, vector<BoundaryCondition1> cond1, vector<BoundaryCondition2> cond2)
{

	// ������ �����������
	vector<double> gglGlobalMatrix(ig[numFinitElem] - 1, 0);
	vector<double> gglGlobalMatrixHard(ig[numFinitElem] - 1, 0);
	vector<double> gglGlobalMatrixMass(ig[numFinitElem] - 1, 0);

	// ������������ �������
	vector<double> diGlobalMatrix(numFinitElem, 0);
	vector<double> diGlobalMatrixHard(numFinitElem, 0);
	vector<double> diGlobalMatrixMass(numFinitElem, 0);

	vector<double> globalVector(numFinitElem, 0);
	vector<double> resultQ2(numFinitElem, 0);
	vector<double> tempQ0(numFinitElem, 0);
	vector<double> tempQ1(numFinitElem, 0);

	BuildGlobalMatrixAndVector(gglGlobalMatrixHard, gglGlobalMatrixMass, diGlobalMatrixHard, diGlobalMatrixMass, globalVector, point, triangle, range);

	double deltaTime = timeT[countTime] - timeT[countTime - 2];
	double deltaTime1 = timeT[countTime - 1] - timeT[countTime - 2];
	double deltaTime0 = timeT[countTime] - timeT[countTime - 1];

	// A = ((deltaT + deltaT0) / (detlaT * deltaT0)) * M + G
	for (int i = 0; i < diGlobalMatrix.size(); i++)
		diGlobalMatrix[i] = (diGlobalMatrixHard[i]) + (diGlobalMatrixMass[i] * ((deltaTime + deltaTime0) / (deltaTime * deltaTime0)) * range.gammaCoff);
	for (int i = 0; i < gglGlobalMatrix.size(); i++)
		gglGlobalMatrix[i] = (gglGlobalMatrixHard[i]) + (gglGlobalMatrixMass[i] * ((deltaTime + deltaTime0) / (deltaTime * deltaTime0)) * range.gammaCoff);

	// ������� �������� M * q0 � M * q1
	for (int i = 0; i < numFinitElem; i++)
	{
		tempQ0[i] = diGlobalMatrixMass[i] * resultQ0[i];
		tempQ1[i] = diGlobalMatrixMass[i] * resultQ1[i];
		for (int j = ig[i]; j < ig[i + 1]; j++)
		{
			tempQ0[i] += gglGlobalMatrixMass[j - 1] * resultQ0[jg[j - 1]];
			tempQ0[jg[j - 1]] += gglGlobalMatrixMass[j - 1] * resultQ0[i];

			tempQ1[i] += gglGlobalMatrixMass[j - 1] * resultQ1[jg[j - 1]];
			tempQ1[jg[j - 1]] += gglGlobalMatrixMass[j - 1] * resultQ1[i];
		}
		
	}

	// d = b - (deltaT0 / (deltaT * deltaT1)) * M * q0 + (deltaT / (deltaT1 * deltaT0)) * M * q1
	for (int i = 0; i < resultQ2.size(); i++)
	{
		tempQ0[i] *= (range.gammaCoff * (deltaTime0 / (deltaTime * deltaTime1))); 
		tempQ1[i] *= (range.gammaCoff * (deltaTime / (deltaTime1 * deltaTime0))); 
		resultQ2[i] = globalVector[i] - tempQ0[i] + tempQ1[i];
	}

	// ���� ������� �������
	EntryBoundaryCondition2(cond2, resultQ2);
	EntryBoundaryCondition1(cond1, diGlobalMatrix, resultQ2);

	MCG(diGlobalMatrix, gglGlobalMatrix, resultQ2);

	// �������� ������� q2
	return resultQ2;
}

int main()
{
	countTime = 0;
	vector<Point> node;
	vector<vector<int>> triangle;
	vector<double> time;
	vector<BoundaryCondition1> boundCond1; // ������ ������� �������
	vector<BoundaryCondition2> boundCond2; // ������ ������� �������
	Range coeff;
	Input(node, triangle, coeff);
	Portrait(triangle, numFinitElem);

	// ������ �����������
	vector<double> gglGlobalMatrix(ig[numFinitElem] - 1, 0);
	vector<double> gglGlobalMatrixHard(ig[numFinitElem] - 1, 0);
	vector<double> gglGlobalMatrixMass(ig[numFinitElem] - 1, 0);

	// ������������ �������
	vector<double> diGlobalMatrix(numFinitElem, 0);
	vector<double> diGlobalMatrixHard(numFinitElem, 0);
	vector<double> diGlobalMatrixMass(numFinitElem, 0);

	vector<double> globalVector(numFinitElem, 0); // ���������� ������ ������ �����

	BuildGlobalMatrixAndVector(gglGlobalMatrixHard, gglGlobalMatrixMass, diGlobalMatrixHard, diGlobalMatrixMass, globalVector, node, triangle, coeff);

	vector <double> resultQ0(numFinitElem, 0);
	InputBoundaryCondition(node, boundCond1, boundCond2);
	EntryBoundaryCondition1(boundCond1, diGlobalMatrixMass, globalVector);

	MCG(diGlobalMatrixMass, gglGlobalMatrixMass, globalVector);
	cout << "Time layer: " << 1 << " (" << timeT[countTime] << ")" << endl;
	for (int i = 0; i < numFinitElem; i++)
	{
		resultQ0[i] = globalVector[i];
		cout << resultQ0[i] << " (" << i << ") " << endl;
	}
	
	countTime++;
	vector<double> resultQ1(numFinitElem, 0);

	// ���������� ������� �������� ��� t * t
	/*globalVector.clear(); globalVector.resize(numFinitElem);
	gglGlobalMatrixMass.clear(); gglGlobalMatrixMass.resize(ig[numFinitElem] - 1);
	diGlobalMatrixMass.clear(); diGlobalMatrixMass.resize(numFinitElem);
	BuildGlobalMatrixAndVector(gglGlobalMatrixHard, gglGlobalMatrixMass, diGlobalMatrixHard, diGlobalMatrixMass, globalVector, node, triangle, coeff);
	EntryBoundaryCondition1(boundCond1, diGlobalMatrixMass, globalVector);
	MCG(diGlobalMatrixMass, gglGlobalMatrixMass, globalVector);
	resultQ1 = globalVector;

	cout << endl << "Time layer: " << 2 << endl;
	for (int i = 0; i < numFinitElem; i++)
		cout << resultQ1[i] << " (" << i << ") " << endl;*/

	for (int i = 0; i < numTriag; i++)
		coeff.functionForAppr[i] = 0;
	resultQ1 = ImplicitTwolayerScheme(resultQ0, node, triangle, coeff, boundCond1, boundCond2);
	cout << endl << "Time layer: " << 2 << " (" << timeT[countTime] << ")"<<endl;
	for (int i = 0; i < numFinitElem; i++)
		cout << resultQ1[i] << " (" << i << ") " << endl;

	vector<double> resultQ2(numFinitElem, 0);
	for (countTime = 2; countTime <= timeIntervals; countTime++)
	{
		resultQ2 = ImplicitThreelayerScheme(resultQ0, resultQ1, node, triangle, coeff, boundCond1, boundCond2);
		cout << endl << "Time layer: " << countTime + 1 << " (" << timeT[countTime] << ")" << endl;
		for (int j = 0; j < numFinitElem; j++)
			cout << resultQ2[j] << " (" << j << ") " << endl;
		resultQ0 = resultQ1;
		resultQ1 = resultQ2;
	}

	system("pause");
	return 0;
}