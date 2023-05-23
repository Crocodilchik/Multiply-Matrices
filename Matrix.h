#pragma once

// ����������� ���� ����������� ���������
#include <iostream>
#include <iomanip>
#include <random>
#include <time.h>
#include <string>
#include <omp.h>

using namespace std;

//����� ������
class Matrix
{

// ���� ������
private:

	int n; // ������� �������
	int** mat; // ������

public:

	// ������������ � ����������
	Matrix(int n = 0, double x = 0) : n(n)
	{

		if (n != 0)
		{

			mat = new int* [n];

			for (int i = 0; i < n; i++)
			{

				mat[i] = new int[n];

				for (int j = 0; j < n; j++) mat[i][j] = x;

			}

		}
		else mat = nullptr;

	}

	Matrix(const Matrix* other): n(other->n)
	{

		if (n != 0)
		{

			mat = new int* [n];

			for (int i = 0; i < n; i++)
			{

				mat[i] = new int[n];

				for (int j = 0; j < n; j++) mat[i][j] = other->mat[i][j];

			}

		}
		else mat = nullptr;

	}

	Matrix(Matrix&& other): n(other.n)
	{

		mat = other.mat;
		other.mat = nullptr;

	}

	~Matrix()
	{

		if (mat)
		{
			for (int i = 0; i < n; i++) delete[] mat[i];
			delete[] mat;
		}

	}

	// ���������
	int* operator[](int i) const; // Rvalue �����
	int* operator[](int i); // Lvalue �����
	Matrix operator=(const Matrix& other); // ������������
	Matrix operator=(Matrix&& other) noexcept; // ������������ ������������
	Matrix operator+(const Matrix& other); // ��������
	Matrix operator-(const Matrix& other); // ���������
	bool operator==(const Matrix& other); // ��������� �� ���������

	// ���������� ���������� ������ �������
	void fill_random();

	//���������
	Matrix multiply_usual(const Matrix& other); // �������
	Matrix multiply_usual_omp(const Matrix& other); // ������� omp
	Matrix multiply_vinograd(const Matrix& other) const; // ���������
	Matrix multiply_vinograd_omp(const Matrix& other) const; // ��������� omp
	Matrix multiply_shtrassen(const Matrix& other, int n1); // �������
	Matrix multiply_shtrassen_omp(const Matrix& other, int n1); // ������� omp

	friend ostream& operator<<(ostream& out,const Matrix& a); // ����� �������

};
