#include "Matrix.h"

int* Matrix::operator[](int i) const
{

	return this->mat[i];

}

int* Matrix::operator[](int i)
{

	return this->mat[i];

}

Matrix Matrix::operator=(const Matrix& other)
{

	if (this != &other) 
	{

		n = other.n;

		mat = new int* [n];

		for (int i = 0; i < n; i++)
		{

			mat[i] = new int[n];

			for (int j = 0; j < n; j++) mat[i][j] = other.mat[i][j];

		}

	}

	return this;
}

Matrix Matrix::operator=(Matrix&& other) noexcept
{

	if (mat != nullptr)
	{

		for (int i = 0; i < n; i++) delete[] mat[i];
		delete[] mat;

	}

	n = other.n;
	mat = other.mat;
	other.mat = nullptr;

	return this;

}

Matrix Matrix::operator+(const Matrix& other)
{

	if (this->n != other.n) throw "Unable because of different sizes: " + to_string(this->n) + " and " + to_string(other.n);

	Matrix temp = Matrix(n, 0);
	Matrix ths = Matrix(this);
	int i, j;

	for (i = 0; i < n; i++)
	{

		for (j = 0; j < n; j++) temp.mat[i][j] = this->mat[i][j] + other.mat[i][j];

	}

	return temp;

}

Matrix Matrix::operator-(const Matrix& other)
{

	if (this->n != other.n) throw "Unable because of different sizes: " + to_string(this->n) + " and " + to_string(other.n);

	Matrix temp = Matrix(n, 0);

	for (int i = 0; i < n; i++)
	{

		for (int j = 0; j < n; j++) temp.mat[i][j] = this->mat[i][j] - other.mat[i][j];

	}

	return temp;

}

bool Matrix::operator==(const Matrix& other)
{

	if (this->n != other.n) return false;

	for (int i = 0; i < n; i++)
	{

		for (int j = 0; j < n; j++)
		{

			if (abs(this->mat[i][j] != other.mat[i][j])) return false;

		}

	}

	return true;

}

ostream& operator<<(ostream& out, const Matrix& a)
{

	for (int i = 0; i < a.n; i++)
	{

		for (int j = 0; j < a.n; j++) out << setw(10) << a[i][j];
		out << endl;

	}
	
	return out;

}

void Matrix::fill_random()
{

	srand(time(NULL));

	for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) mat[i][j] = rand()* 0.001;

}

// ���������

Matrix Matrix::multiply_usual(const Matrix& other)
{

	// �������� �� ������������ ��������
	if (this->n != other.n) throw "Unable because of different sizes: " + to_string(this->n) + " and " + to_string(other.n);

	// ������� ������������
	Matrix temp = Matrix(n, 0);

	// ������ ������� ������������
	for (int i = 0; i < n; i++)
	{

		// ������� ������� ������������
		for (int j = 0; j < n; j++)
		{

			//���� ��� ���������� ������������ ������ �� �������
			for (int k = 0; k < n; k++)
			{

				temp[i][j] += this->mat[i][k] * other.mat[k][j];

			}

		}

	}

	// �����
	return temp;

}

Matrix Matrix::multiply_usual_omp(const Matrix& other)
{
	// �������� �� ������������ ��������
	if (this->n != other.n) throw "Unable because of different sizes: " + to_string(this->n) + " and " + to_string(other.n);

	// ������� ������������
	Matrix temp = Matrix(n, 0);
	// ����� this ��� �������
	Matrix ths = Matrix(this);

	// ���������
	int i, j, k;

	// ����������������� �� ������ ������� ������������
#pragma omp parallel for shared(ths, other, temp) private(i, j, k)
	for (i = 0; i < n; i++)
	{

		// ������� ������� ������������
		for (j = 0; j < n; j++)
		{

			//���� ��� ���������� ������������ ������ �� �������
			for (k = 0; k < n; k++)
			{

				temp[i][j] += ths.mat[i][k] * other.mat[k][j];

			}

		}

	}

	// �����
	return temp;

}

Matrix Matrix::multiply_vinograd(const Matrix& other) const
{

	// �������� �� ������������ ��������
	if (this->n != other.n) throw "Unable because of different sizes: " + to_string(this->n) + " and " + to_string(other.n);

	Matrix temp = Matrix(n, 0); // ������� ������������
	double* mulH = new double[n]; // ��������� ������
	double* mulV = new double[n]; // ��������� �������

	// ���������� ���������� ������ � �������
	for (int i = 0; i < n; ++i) {

		mulH[i] = 0;
		mulV[i] = 0;

		for (int j = 0; j < n / 2; ++j) {

			mulH[i] += this->mat[i][2 * j] * this->mat[i][2 * j + 1];
			mulV[i] += other.mat[2 * j][i] * other.mat[2 * j + 1][i];

		}

	}

	

	// ������������ ��������� �� ������� ���������
	for (int i = 0; i < n; ++i) {

		for (int j = 0; j < n; ++j) {

			temp.mat[i][j] = -mulH[i] - mulV[j];

			for (int k = 0; k < n / 2; ++k) {

				temp.mat[i][j] += (this->mat[i][2 * k] + other.mat[2 * k + 1][j]) * (this->mat[i][2 * k + 1] + other.mat[2 * k][j]);
			}
		}
	}

	// ���������� ������ n-�������� (���� ���������)
	if (n % 2) {
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				temp.mat[i][j] += this->mat[i][n - 1] * other.mat[n - 1][j];
			}
		}
	}

	// ������������ ������
	delete[] mulH;
	delete[] mulV;

	// �����
	return temp;

}

Matrix Matrix::multiply_vinograd_omp(const Matrix& other) const
{

	// �������� �� ������������ ��������
	if (this->n != other.n) throw "Unable because of different sizes: " + to_string(this->n) + " and " + to_string(other.n);

	Matrix temp = Matrix(n, 0);// ������� ������������
	Matrix ths = Matrix(this); // ����� this ��� �������
	double* mulH = new double[n]; // ��������� ������
	double* mulV = new double[n]; // ��������� �������

	// ���������
	int i, j, k;

	// ���������������� ���������� ���������� ������ � �������
#pragma omp parallel for shared(ths, other, temp) private(i, j, k)
	for (i = 0; i < n; ++i) {

		mulH[i] = 0;
		mulV[i] = 0;

		for (j = 0; j < n / 2; ++j) {

			mulH[i] += this->mat[i][2 * j] * this->mat[i][2 * j + 1];
			mulV[i] += other.mat[2 * j][i] * other.mat[2 * j + 1][i];

		}

	}

	// ���������������� ������������ ���������
#pragma omp parallel for shared(ths, other, temp) private(i, j, k)
	for (i = 0; i < n; ++i) {

		for (j = 0; j < n; ++j) {

			temp.mat[i][j] = -mulH[i] - mulV[j];

			for (k = 0; k < n / 2; ++k) {

				temp.mat[i][j] += (this->mat[i][2 * k] + other.mat[2 * k + 1][j]) * (this->mat[i][2 * k + 1] + other.mat[2 * k][j]);
			}
		}
	}

	// ���������������� ���������� ������ n-�������� (���� ���������)
	if (n % 2) {
#pragma omp parallel for shared(ths, other, temp) private(i, j, k)
		for (i = 0; i < n; ++i) {
			for (j = 0; j < n; ++j) {
				temp.mat[i][j] += this->mat[i][n - 1] * other.mat[n - 1][j];
			}
		}
	}

	// ������������ ������
	delete[] mulH;
	delete[] mulV;

	// �����
	return temp;

}

Matrix Matrix::multiply_shtrassen(const Matrix& other, int n1)
{
	// �������� �� ������������ ��������
	if (this->n != other.n) throw "Unable because of different sizes: " + to_string(this->n) + " and " + to_string(other.n);

	int i, j, k; // ���������
	Matrix a11 = Matrix(n, 0), a12 = Matrix(n, 0), a21 = Matrix(n, 0), a22 = Matrix(n, 0); // ���������� ��� this
	Matrix b11 = Matrix(n, 0), b12 = Matrix(n, 0), b21 = Matrix(n, 0), b22 = Matrix(n, 0); // ���������� ��� other
	Matrix c11 = Matrix(n, 0), c12 = Matrix(n, 0), c21 = Matrix(n, 0), c22 = Matrix(n, 0); // ���������� ��� ������
	Matrix p1 = Matrix(n, 0), p2 = Matrix(n, 0), p3 = Matrix(n, 0), p4 = Matrix(n, 0), p5 = Matrix(n, 0), p6 = Matrix(n, 0), p7 = Matrix(n, 0);



	// ��������� ���������� ������ �������
	for (i = 0; i < n1 / 2; i++)
	{
		for (j = 0; j < n1 / 2; j++)
		{
			a11[i][j] = this->mat[i][j];
			b11[i][j] = other.mat[i][j];

			a12[i][j] = this->mat[i][j + n1 / 2];
			b12[i][j] = other.mat[i][j + n1 / 2];

			a21[i][j] = this->mat[i + n1 / 2][j];
			b21[i][j] = other.mat[i + n1 / 2][j];

			a22[i][j] = this->mat[i + n1 / 2][j + n1 / 2];
			b22[i][j] = other.mat[i + n1 / 2][j + n1 / 2];
		}
	}

		p1 = (a11 + a22).multiply_vinograd(b11 + b22);
		p2 = (a21 + a22).multiply_vinograd(b11);
		p3 = a11.multiply_vinograd(b12 - b22);
		p4 = a22.multiply_vinograd(b21 - b11);
		p5 = (a11 + a12).multiply_vinograd(b22);
		p6 = (a21 - a11).multiply_vinograd(b11 + b12);
		p7 = (a12 - a22).multiply_vinograd(b21 + b22);


	c11 = p1 + p4 - p5 + p7;
	c12 = p3 + p5;
	c21 = p2 + p4;
	c22 = p1 + p3 - p2 + p6;

	//cout << c11 << endl;

	Matrix temp = Matrix(n, 0);

	for (int i = 0; i < n1 / 2; i++)
	{
		for (int j = 0; j < n1 / 2; j++)
		{
			temp.mat[i][j] = c11[i][j];
		}
	}

	for (int i = 0; i < n1 / 2; i++)
	{
		for (int j = n1 / 2; j < n1; j++)
		{
			temp.mat[i][j] = c12[i][j - n1 / 2];
		}
	}

	for (int i = n1 / 2; i < n1; i++)
	{
		for (int j = 0; j < n1 / 2; j++)
		{
			temp.mat[i][j] = c21[i - n1 / 2][j];
		}
	}

	for (int i = n1 / 2; i < n1; i++)
	{
		for (int j = n1 / 2; j < n1; j++)
		{
			temp.mat[i][j] = c22[i - n1 / 2][j - n1 / 2];
		}

	}

	return temp;
}

Matrix Matrix::multiply_shtrassen_omp(const Matrix& other, int n1)
{
	if (this->n != other.n) throw "Unable because of different sizes: " + to_string(this->n) + " and " + to_string(other.n);

	int i, j, k;
	Matrix a11 = Matrix(n, 0), a12 = Matrix(n, 0), a21 = Matrix(n, 0), a22 = Matrix(n, 0);
	Matrix b11 = Matrix(n, 0), b12 = Matrix(n, 0), b21 = Matrix(n, 0), b22 = Matrix(n, 0);
	Matrix c11 = Matrix(n, 0), c12 = Matrix(n, 0), c21 = Matrix(n, 0), c22 = Matrix(n, 0);
	Matrix ths = Matrix(this);
	Matrix p1 = Matrix(n, 0), p2 = Matrix(n, 0), p3 = Matrix(n, 0), p4 = Matrix(n, 0), p5 = Matrix(n, 0), p6 = Matrix(n, 0), p7 = Matrix(n, 0);


	// ��������� ���������� ������ �������
//#pragma omp parallel for shared(ths, other, a11,a12,a21,a22,b11,b12,b21,b22,n1) private(i, j, k)
	for (i = 0; i < n1 / 2; i++)
	{
		for (j = 0; j < n1 / 2; j++)
		{
			a11[i][j] = ths.mat[i][j];
			b11[i][j] = other.mat[i][j];

			a12[i][j] = ths.mat[i][j + n1 / 2];
			b12[i][j] = other.mat[i][j + n1 / 2];

			a21[i][j] = ths.mat[i + n1 / 2][j];
			b21[i][j] = other.mat[i + n1 / 2][j];

			a22[i][j] = ths.mat[i + n1 / 2][j + n1 / 2];
			b22[i][j] = other.mat[i + n1 / 2][j + n1 / 2];
		}
	}

	
#pragma omp parallel shared(a11,a12,a21,a22,b11,b12,b21,b22,n1)
#pragma omp sections
		{

#pragma omp section
			{
				p1 = (a11 + a22).multiply_vinograd(b11 + b22);
			}
#pragma omp section
			{
				p2 = (a21 + a22).multiply_vinograd(b11);
			}
#pragma omp section
			{
				p3 = a11.multiply_vinograd(b12 - b22);
			}
#pragma omp section
			{
				p4 = a22.multiply_vinograd(b21 - b11);
			}
#pragma omp section
			{
				p5 = (a11 + a12).multiply_vinograd(b22);
			}
#pragma omp section
			{
				p6 = (a21 - a11).multiply_vinograd(b11 + b12);
			}
#pragma omp section
			{
				p7 = (a12 - a22).multiply_vinograd(b21 + b22);
			}

		}

	c11 = p1 + p4 - p5 + p7;
	c12 = p3 + p5;
	c21 = p2 + p4;
	c22 = p1 + p3 - p2 + p6;

	Matrix temp = Matrix(n, 0);

	for (int i = 0; i < n1 / 2; i++)
	{
		for (int j = 0; j < n1 / 2; j++)
		{
			temp.mat[i][j] = c11[i][j];
		}
	}

	for (int i = 0; i < n1 / 2; i++)
	{
		for (int j = n1 / 2; j < n1; j++)
		{
			temp.mat[i][j] = c12[i][j - n1 / 2];
		}
	}

	for (int i = n1 / 2; i < n1; i++)
	{
		for (int j = 0; j < n1 / 2; j++)
		{
			temp.mat[i][j] = c21[i - n1 / 2][j];
		}
	}

	for (int i = n1 / 2; i < n1; i++)
	{
		for (int j = n1 / 2; j < n1; j++)
		{
			temp.mat[i][j] = c22[i - n1 / 2][j - n1 / 2];
		}

	}

	return temp;
}