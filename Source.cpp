#include "Matrix.h"


constexpr int SIZE = 1500;

int main()
{

	cout << "Count of elements: " << SIZE << endl << endl;
	cout  << setiosflags(ios::left) << setw(15) << "Kind of mult" << setw(15) << "Time" << setw(15) << "Correct" << endl;
	cout << endl;

	Matrix a = Matrix(SIZE, 0);
	Matrix b = Matrix(SIZE, 0);

	a.fill_random();
	b.fill_random();

	// Обычное умножение

	double tstart, tend, ftime, stime;

	tstart = omp_get_wtime();
	Matrix standart = a.multiply_usual(b);
	tend = omp_get_wtime();

	ftime = tend - tstart;

	cout << setiosflags(ios::left) << setw(15) << "Usual" << setw(15) << ftime << setw(15) << "Standart" << endl;

	tstart = omp_get_wtime();
	Matrix c1 = a.multiply_usual_omp(b);
	tend = omp_get_wtime();

	stime = tend - tstart;

	cout << setiosflags(ios::left) << setw(15) << "Usual omp" << setw(15) << stime << setw(15) << (c1 == standart)  << endl;
	cout << endl << "Efficiency: " << ftime / stime << endl << endl;

	// Умножение Винограда

	tstart = omp_get_wtime();
	Matrix c2 = a.multiply_vinograd(b);
	tend = omp_get_wtime();

	ftime = tend - tstart;

	cout << setiosflags(ios::left) << setw(15) << "Vinograd" << setw(15) << ftime << setw(15) << (c2 == standart) << endl;

	tstart = omp_get_wtime();
	c1 = a.multiply_vinograd_omp(b);
	tend = omp_get_wtime();

	stime = tend - tstart;

	cout << setiosflags(ios::left) << setw(15) << "Vinograd omp" << setw(15) << stime << setw(15) << (c1 == standart) << endl;
	cout << endl << "Efficiency: " << ftime / stime << endl << endl;

	// Умножение Штрассена
	
	tstart = omp_get_wtime();
	c2 = a.multiply_shtrassen(b, SIZE);
	tend = omp_get_wtime();

	ftime = tend - tstart;

	cout << setiosflags(ios::left) << setw(15) << "Strassen" << setw(15) << ftime << setw(15) << (c2 == standart) << endl;

	tstart = omp_get_wtime();
	c1 = a.multiply_shtrassen_omp(b, SIZE);
	tend = omp_get_wtime();

	stime = tend - tstart;

	cout << setiosflags(ios::left) << setw(15) << "Strassen omp" << setw(15) << stime << setw(15) << (c1 == standart) << endl;
	cout << endl << "Efficiency: " << ftime / stime << endl << endl;

	return 0;

}