#pragma once

// Подключение всех необходимых библиотек
#include <iostream>
#include <iomanip>
#include <random>
#include <time.h>
#include <string>
#include <omp.h>

using namespace std;

//класс матриц
class Matrix
{

// Поля класса
private:

	int n; // порядок матрицы
	int** mat; // мтрица

public:

	// Конструкторы и деструктор
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

	// Операторы
	int* operator[](int i) const; // Rvalue вызов
	int* operator[](int i); // Lvalue вызов
	Matrix operator=(const Matrix& other); // присваивание
	Matrix operator=(Matrix&& other) noexcept; // перемещающее присваивание
	Matrix operator+(const Matrix& other); // сложение
	Matrix operator-(const Matrix& other); // вычитание
	bool operator==(const Matrix& other); // сравнение на равенство

	// Заполнение рандомными целыми числами
	void fill_random();

	//умножение
	Matrix multiply_usual(const Matrix& other); // обычное
	Matrix multiply_usual_omp(const Matrix& other); // обычное omp
	Matrix multiply_vinograd(const Matrix& other) const; // Винограда
	Matrix multiply_vinograd_omp(const Matrix& other) const; // Винограда omp
	Matrix multiply_shtrassen(const Matrix& other, int n1); // Штрасса
	Matrix multiply_shtrassen_omp(const Matrix& other, int n1); // Штрасса omp

	friend ostream& operator<<(ostream& out,const Matrix& a); // вывод матрицы

};
