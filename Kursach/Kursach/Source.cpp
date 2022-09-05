#include<cmath>
#include<fstream>
#include<iostream>
#include<vector>
#include<algorithm>
#include <set>

using namespace std;

int numVert; // количество вершин на расчётной области
int numTriag; // количество треугольников
int numFinitElem; // количество конечных элементов на всей расчетной области
int timeIntervals; // количество временных интервалов
vector<double> timeT; // вектор времени
int countTime = 0; // счетчик временных слоёв

vector<int> ig; // строки глобальной матрицы в разреженном формате
vector<int> jg; // столбцы глобальной матрицы в разреженном формате


struct Point
{
	int glNumber; // глобальный номер
	double x, y; // координаты узла
};

struct Range
{
	vector<double> diffusionCoff; // значение коэффициента лямбда
	double gammaCoff; // значение коэффициента гамма
	vector<int> functionValue; // значение функции
	vector<int> functionForAppr; // значение функции начального приближения
	double FormulaForFunction(double x, double y, int numTriag) // значения функции f
	{
		switch (functionValue[numTriag])
		{
		case 0:
			return 0; // тест 1, 2, 3
		case 1:
			return -6; // тест 4
		case 2:
			return -8; // тест 5
		case 3:
			return -18 * x; // тест 6
		case 4:
			return -6 * y; // тест 7
		case 5:
			return -12 * x * x; // тест 8
		case 6:
			return -48 * y * y; // тест 9
		case 7:
			return (3 * y * y) - (2 * timeT[countTime]); // тест 11 (u = y * y * t)
		case 8:
			return (4 * y * y * y * timeT[countTime]) - (6 * timeT[countTime] * timeT[countTime] * y); // тест 12 (u = y * y * y * t * t)
		case 9:
			return 7; // тест 10 (u = 7t)

		case 10: // тест 13 (u = y * y * y * y * t * t * t)
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

	double FormulaForApprox(double x, double y, int numTriag) // значения функции u
	{
		switch (functionForAppr[numTriag])
		{
		case 0:
			return 0;
		case 1:
			return 7 * timeT[countTime]; // тест 10 (u = 7 * t)
		case 2:
			return y * y * timeT[countTime]; // тест 11 (u = y * y * t)
		case 3:
			return y * y * y * timeT[countTime] * timeT[countTime]; // тест 12 (u = y * y * y * t * t)
		case 4:
			return 7; // тест 1 (u = 7)
		case 5:
			return x; // тест 2 (u = x)
		case 6:
			return y; // тест 3 (u = y) 
		case 7:
			return x * x; // тест 4 (u = x * x)
		case 8:
			return y * y; // тест 5 (u = y * y)
		case 9:
			return x * x * x; // тест 6 (u = x * x * x)
		case 10:
			return y * y * y; // тест 7 (u = y * y * y)
		case 11:
			return x * x * x * x; // тест 8 (u = x * x * x * x)
		case 12:
			return y * y * y * y; // тест 9 (u = y * y * y * y)
		case 13:
			return y * y * y * y * timeT[countTime] * timeT[countTime]; // тест 12 (u = y * y * y * y * t * t * t)
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

// градиент от L
struct gradL
{
	double comp1; // первый параметр для подсчёта градиента
	double comp2; // второй параметр для подсчёта градиента
};

// Первое краевое условие
struct BoundaryCondition1
{
	Point node; // узел на котором задаются краевые условия
	int formulaNum; // выбор значения первого краевого условия

	double FormulaForCondition1() // задаёт значение первого краевого условия
	{
		switch (formulaNum)
		{
		case 0:
			return 7; // для теста 1
		case 1:
			return 1; // для тестов 2, 3, 4, 5, 6, 7, 8, 9
		case 2:
			return 7; // для тестов 2, 3
		case 3:
			return 49; // для тестов 4, 5
		case 4:
			return 343; // для тестов 6, 7
		case 5:
			return 2401; // для тестов 8, 9
		case 6:
			return timeT[countTime]; // для теста 11
		case 7:
			return 49 * timeT[countTime]; // для теста 11
		case 8:
			return timeT[countTime] * timeT[countTime]; // для теста 12
		case 9:
			return 343 * timeT[countTime] * timeT[countTime]; // для теста 12
		case 10:
			return 7 * timeT[countTime]; // для теста 10
		case 11:
			return timeT[countTime] * timeT[countTime] * timeT[countTime]; // для теста 13
		case 12:
			return 2401 * timeT[countTime] * timeT[countTime] * timeT[countTime]; // для теста 13
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

// Второе краевое условие
struct BoundaryCondition2
{
	Point node[4]; // ребро на котором задаются краевые условия
	int formulaNum; // выбор значения второго краевого условия

	double FormulaForCondition2() // задаёт значение второго краевого условия
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


// Функция ввода координат вершин треугольников, глобальных номеров всех узлов треугольника, параметров расчетной области и сведений о времени
// параметры: 
//            point - глобальный номер вершины треугольника
//            triangle - глобальные номера узлов треугольника
//            range - параметры расчётной области
void Input(vector<Point>& point, vector<vector<int>>& triangle, Range& range)
{
	// сведения об узлах
	ifstream finXY("coordxy.txt");
	finXY >> numVert >> numFinitElem;
	point.resize(numVert);
	for (int i = 0; i < numVert; i++)
	{
		finXY >> point[i].x >> point[i].y;
		point[i].glNumber = i;
	}
	finXY.close();

	// сведения о треугольниках
	ifstream finParam("m.txt");
	finParam >> numTriag;
	triangle.resize(numTriag, vector<int>(10, 0));
	for (int i = 0; i < numTriag; i++)
		for (int j = 0; j < 10; j++)
			finParam >> triangle[i][j];

	finParam.close();

	// сведения о параметрах поля
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

	// сведения о времени
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

// Функция автоматической генерации портрета матрицы, т.е. ig и jg
// параметры: 
//            triangle - глобальные номера узлов треугольника
//            n - количество всех узлов на расчётной области
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

// Функция подсчёта факториала (для сборки локальной матрицы)
// параметры: 
//            num - число, возводимое в факториал
//            res - результат возведение в факториал
void Factorial(int num, double& res)
{
	res = 1;
	for (int i = 2; i <= num; i++)
		res *= i;
}

// Функция подсчёта интегралов (для сборки локальной матрицы)
// параметры: 
//            sign - знак интеграла
//            val - вектор, в котором хранятся L1, L2, L3
//            res - результат подсчёта интеграла
//            detD - определитель матрицы
void CalculationIntegralL(double sign, vector<int> val, double& res, double detD)
{

	vector<double> rs(4, 0);
	Factorial(val[0], rs[0]);
	Factorial(val[1], rs[1]);
	Factorial(val[2], rs[2]);
	Factorial(val[0] + val[1] + val[2] + 2, rs[3]);
	res = (rs[0] * rs[1] * rs[2]) / (rs[3]) * sign * detD;
}

// Функция подсчёта скалярного произведения gladL(n) на gladL(k)
// параметры: 
//            g1 - gladL(n)
//            g2 - gladL(k)
double Scalar(gradL g1, gradL g2)
{
	return g1.comp1 * g2.comp1 + g1.comp2 * g2.comp2;
}

// Функция построение локальных матрицы
// параметры: 
//            point - глобальный номер вершины треугольника
//            triangle - глобальные номера узлов треугольника
//            localMatrixHard - локальная матрица жесткости
//            localMatrixMass - локальная матрица масс
//            numLocalTriag - номер локального треугольника
void BuildLocalMatrix(vector<Point> point, vector<vector<int>> triangle, vector<vector<double>>& localMatrixHard, vector<vector<double>>& localMatrixMass, int numLocalTriag)
{
	double detD = ((point[1].x - point[0].x) * (point[2].y - point[0].y)) - ((point[2].x - point[0].x) * (point[1].y - point[0].y));


	double alpha[3][3] =
	{
		{(point[1].x * point[2].y - point[2].x * point[1].y) / detD, (point[1].y - point[2].y) / detD, (point[2].x - point[1].x) / detD},
		{(point[2].x * point[0].y - point[0].x * point[2].y) / detD, (point[2].y - point[0].y) / detD, (point[0].x - point[2].x) / detD},
		{(point[0].x * point[1].y - point[1].x * point[0].y) / detD, (point[0].y - point[1].y) / detD, (point[1].x - point[0].x) / detD},
	}; // считаем значения альф

	detD = abs(detD); // считаем определитель и находим его абсолютное значение

	// считаем значения gradL1, gradL2, gradL3
	vector<gradL> gradLN(3);
	for (int i = 0; i < 3; i++)
	{
		gradLN[i].comp1 = alpha[i][1];
		gradLN[i].comp2 = alpha[i][2];
	}

	// строим массив степеней L-координат в gradPsi
	// используется для построения матрицы жесткости
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

	// строим массив степеней L-координат в Psi
	// используется для построения матрицы масс
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

	// строим массив коэффициентов в gradPsi
	// используется для матрицы жесткости
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

	// строим массив коэффициентов в Psi
	// используется в матрице масс
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
	double valIntegralDop, valIntegral = 0; // переменные для подсчёта интегралов 
	// считаем матрицу жесткости
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

	// считаем матрицу масс
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

// Генерация локального вектора
// параметры: 
//            node - глобальный номер вершины треугольника
//            localVector - локальный вектор правой части
//            range - параметр расчётной области
void BuildLocalVector(vector<Point> node, vector<double>& localVector, Range range, int numTriag)
{
	// квадратуры
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
	// веса
	double w[21] =
	{ 0.0519871420646, 0.0519871420646, 0.0519871420646, 0.0707034101784, 0.0707034101784, 0.0707034101784,
		0.0707034101784, 0.0707034101784, 0.0707034101784, 0.0909390760952, 0.0909390760952, 0.0909390760952, 0.0909390760952,
		0.0909390760952, 0.0909390760952, 0.1032344051380, 0.1032344051380, 0.1032344051380, 0.1881601469167, 0.1881601469167,
		0.1881601469167
	};

	// локальные координаты для каждой квадратуры
	double LQ[3];
	// значения x и y в каждой квадратуре
	double xQ = 0, yQ = 0;
	// вспомогательная переменная
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


		// значения пси в каждой квадратуре
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

// Функция для занесения локальной матрицы в глобальную
// параметры: 
//            num - глобальное место в разреженном формате
//            vectorNode - вектор всех глобальных номеров узлов в одном треугольнике
//            n - глобальный номер узла
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

// Функция для занесения локальной матрицы и локального вектора правой части в глобальную матрицу и глобальный вектор соответственно
// параметры: 
//            gglGlobalMatrixHard - нижний треугольник глобальной матрицы жесткости
//            gglGlobalMatrixMass - нижний треугольник глобальной матрицы масс
//            diGlobalMatrixHard - диагональ глобальной матрицы жесткости
//            diGlobalMatrixMass - диагональ глобальной матрицы масс
//            globalVector - глобальный вектор правой части
//            point - глобальный номер вершины треугольника
//            triangle - глобальные номера узлов треугольника
//            range - параметры расчётной области
void BuildGlobalMatrixAndVector(vector<double>& gglGlobalMatrixHard, vector<double>& gglGlobalMatrixMass, vector<double>& diGlobalMatrixHard, vector<double>& diGlobalMatrixMass, vector<double>& globalVector, vector<Point> point, vector<vector<int>> triangle, Range range)
{
	vector<double> localVector(10, 0);
	vector<vector<double>> localMatrixHard(10, vector<double>(10, 0));
	vector<vector<double>> localMatrixMass(10, vector<double>(10, 0));

	/*vector<vector<double>> globalMatrixMassPL(16, vector<double>(16, 0));
	vector<vector<double>> globalMatrixMassRAZ(16, vector<double>(16, 0));*/


	// numLocalTriag - кол-во треугольников
	for (int numLocalTriag = 0; numLocalTriag < numTriag; numLocalTriag++)
	{
		// node - узлы для локальной матрицы и вектора(3 вершины)
		vector<Point> node(3);
		for (int i = 0; i < 3; i++)
		{
			node[i].glNumber = triangle[numLocalTriag][i];
			node[i].x = point[node[i].glNumber].x;
			node[i].y = point[node[i].glNumber].y;
		}

		// получаем локальную матрицу
		BuildLocalMatrix(node, triangle, localMatrixHard, localMatrixMass, numLocalTriag);
		for (int i = 0; i < localMatrixHard.size(); i++)
			for (int j = 0; j < localMatrixHard[i].size(); j++)
				localMatrixHard[i][j] *= range.diffusionCoff[numLocalTriag];

		// получаем локальный вектор
		BuildLocalVector(node, localVector, range, numLocalTriag);

		vector<int> globalNode(10, 0);
		// вектор для перевода в глобальную матриц
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


		// заносим локальный вектор в глобальный
		for (int i = 0; i < 10; i++)
		{
			globalVector[globalNode[i]] += localVector[i];
		}


		// обнуление локальных матриц
		for (int i = 0; i < 10; i++)
		{
			for (int j = 0; j < 10; j++)
			{
				localMatrixHard[i][j] = 0;
				localMatrixMass[i][j] = 0;
			}
		}


		// обнуление локального вектора
		for (int i = 0; i < 10; i++)
			localVector[i] = 0;

	}

}

// Функция чтения краевых условий
// параметры: 
//            nodes - глобальный номер вершины треугольника, на которых учитываются краевые условия
//            condition1 - краевые условия первого рода
//            condition2 - краевые условия второго рода
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

// Функция учета краевых условий 1 рода
// параметры: 
//            condition - краевые условия первого рода
//            di - диагональ матрицы
//            globalVector - вектор правой части
void EntryBoundaryCondition1(vector<BoundaryCondition1> condition, vector<double>& di, vector<double>& globalVector)
{

	for (int i = 0; i < condition.size(); i++)
	{
		di[condition[i].node.glNumber] = 1.0E50;
		globalVector[condition[i].node.glNumber] = 1.0E50 * condition[i].FormulaForCondition1();
	}
}

// Учет краевых условий 2 рода
// Функция учета краевых условий 2 рода
// параметры: 
//            condition - краевые условия второго рода
//            globalVector - вектор правой части
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

// Функция скалярного происведения векторов
// параметры: 
//            x - первый вектор
//            y - второй вектор
double ScalarProduct(vector<double>& x, vector<double>& y)
{
	double scalarProduct = 0;
	for (int i = 0; i < x.size(); i++)
		scalarProduct += x[i] * y[i];
	return scalarProduct;
}

// Функция произведения матрицы на вектор
// параметры: 
//            di - диагональ матрицы
//            ggl - нижний треугольник матрицы
//            x - вектор
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

// Функция для решения СЛАУ методом сопряженных градиентов
// параметры: 
//            di - диагональ матрицы
//            ggl - нижний треугольник матрицы
//            B - вектор правой части (результат перезаписывается)
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

// Функция для нахождения Q1 c помощью неявной двухслойной схемы
// параметры: 
//            resultQ0 - начальное приближение Q0
//            point - глобальный номер вершины треугольника
//            triangle - глобальные номера узлов треугольника
//            range - параметры расчётной области
//            condition1 - краевые условия первого рода
//            condition2 - краевые условия второго рода
// результат: 
//            вектор Q1
vector<double> ImplicitTwolayerScheme(vector<double> resultQ0, vector<Point> point, vector<vector<int>> triangle, Range range, vector<BoundaryCondition1> cond1, vector<BoundaryCondition2> cond2)
{

	// нижний треугольник
	vector<double> gglGlobalMatrix(ig[numFinitElem] - 1, 0);
	vector<double> gglGlobalMatrixHard(ig[numFinitElem] - 1, 0);
	vector<double> gglGlobalMatrixMass(ig[numFinitElem] - 1, 0);

	// диагональные элемнты
	vector<double> diGlobalMatrix(numFinitElem, 0);
	vector<double> diGlobalMatrixHard(numFinitElem, 0);
	vector<double> diGlobalMatrixMass(numFinitElem, 0);

	vector<double> globalVector(numFinitElem, 0);
	vector<double> resultQ1(numFinitElem, 0);
	vector<double> temp(numFinitElem, 0);

	double deltaTime = timeT[countTime] - timeT[countTime - 1]; // значение deltaT

	// Получение глобальных матриц жесткости и масс в разреженном формате
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

	// Учёт краевых условий
	EntryBoundaryCondition2(cond2, resultQ1);
	EntryBoundaryCondition1(cond1, diGlobalMatrix, resultQ1);

	MCG(diGlobalMatrix, gglGlobalMatrix, resultQ1);

	// Получаем q1
	return resultQ1;

}

// Функция для нахождения Q2 c помощью неявной трехслойной схемы
// параметры: 
//            resultQ0 - начальное приближение Q0
//            resultQ1 - значение Q1, посчитанное в помощью неявной двухслойной схемой
//            point - глобальный номер вершины треугольника
//            triangle - глобальные номера узлов треугольника
//            range - параметры расчётной области
//            condition1 - краевые условия первого рода
//            condition2 - краевые условия второго рода
// результат: 
//            вектор Q2
vector<double> ImplicitThreelayerScheme(vector<double> resultQ0, vector<double> resultQ1, vector<Point> point, vector<vector<int>> triangle, Range range, vector<BoundaryCondition1> cond1, vector<BoundaryCondition2> cond2)
{

	// нижний треугольник
	vector<double> gglGlobalMatrix(ig[numFinitElem] - 1, 0);
	vector<double> gglGlobalMatrixHard(ig[numFinitElem] - 1, 0);
	vector<double> gglGlobalMatrixMass(ig[numFinitElem] - 1, 0);

	// диагональные элемнты
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

	// Находим значения M * q0 и M * q1
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

	// Учёт краевых условий
	EntryBoundaryCondition2(cond2, resultQ2);
	EntryBoundaryCondition1(cond1, diGlobalMatrix, resultQ2);

	MCG(diGlobalMatrix, gglGlobalMatrix, resultQ2);

	// Получаем решение q2
	return resultQ2;
}

int main()
{
	countTime = 0;
	vector<Point> node;
	vector<vector<int>> triangle;
	vector<double> time;
	vector<BoundaryCondition1> boundCond1; // первые краевые условия
	vector<BoundaryCondition2> boundCond2; // вторые краевые условия
	Range coeff;
	Input(node, triangle, coeff);
	Portrait(triangle, numFinitElem);

	// нижний треугольник
	vector<double> gglGlobalMatrix(ig[numFinitElem] - 1, 0);
	vector<double> gglGlobalMatrixHard(ig[numFinitElem] - 1, 0);
	vector<double> gglGlobalMatrixMass(ig[numFinitElem] - 1, 0);

	// диагональные элемнты
	vector<double> diGlobalMatrix(numFinitElem, 0);
	vector<double> diGlobalMatrixHard(numFinitElem, 0);
	vector<double> diGlobalMatrixMass(numFinitElem, 0);

	vector<double> globalVector(numFinitElem, 0); // глобальный вектор правой части

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

	// нахождение точного значения для t * t
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