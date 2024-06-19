#include <stdio.h>
#include <unistd.h>
#include <math.h>
#if __unix__
#include "gnuplot_i.h"
#endif

// Исходная производная функция dy/dx = Nx
double G(double const x, double const u)
{
    return 8*u;
}

// Точное решение ДУ первого порядка
double exactSolution(double x)
{
    return exp(8*x);
}

// Метод Эйлера для численного приближения ДУ перого порядка ломаными 
void EulerMethod(double *x, double *y, int n, double h, double x0, double y0)
{
    for (int i = 0; i <= n; i++)
    {
	x[i] = x0 + i * h;
    }
    
    y[0] = 1;
    
    // Применяем разностную схему для расчета аппроксимации в узлах
    for (int i = 0; i < n; i++)
    {
	y[i+1] = y[i] + h * G(x[i], y[i]);
    }
}

//  Метод Рунге-Кутты 2-го порядка для численного решения ДУ перого порядка
void rungeKutta2nd(double *x, double *y, int n, double h, double x0, double y0)
{
    double k1, k2;

    for (int i = 0; i <= n; i++)
    {
	x[i] = x0 + i * h;
    }
    
    y[0] = 1;

    // Вычисляем коэффициенты и применяем разностную схему для расчета аппроксимации в узлах
    for (int i = 0; i < n; i++)
    {
        k1 = G(x[i], y[i]);
        k2 = G(x[i] + h, y[i] + h * k1);
 
        y[i+1] = y[i] + 0.5 * h * (k1 + k2);
 
        x[i+1] = x[i] + h;
    }
}

//  Метод Рунге-Кутты 3-го порядка для численного решения ДУ перого порядка
void rungeKutta3nd(double *x, double *y, int n, double h, double x0, double y0)
{
    double k1, k2, k3;

    for (int i = 0; i <= n; i++)
    {
	x[i] = x0 + i * h;
    }
    
    y[0] = 1;

    // Вычисляем коэффициенты и применяем разностную схему для расчета аппроксимации в узлах
    for (int i = 0; i < n; i++)
    {
        k1 = G(x[i], y[i]);
        k2 = G(x[i] + h/3.0, y[i] + h * k1/3.0);
	k3 = G(x[i] + h/3.0, y[i] + 2 * h * k2/3.0);
	
        y[i+1] = y[i] + (h / 4.0) * (k1 + 3*k3);
 
        x[i+1] = x[i] + h;
    }
}

//  Метод Рунге-Кутты 4-го порядка для численного решения ДУ перого порядка
void rungeKutta4nd(double *x, double *y, int n, double h, double x0, double y0)
{
    double k1, k2, k3, k4;

    for (int i = 0; i <= n; i++)
    {
	x[i] = x0 + i * h;
    }

    y[0] = 1;

    // Вычисляем коэффициенты и применяем разностную схему для расчета аппроксимации в узлах
    for (int i = 0; i < n; i++)
    {
        k1 = G(x[i], y[i]);
        k2 = G(x[i] + h/2.0, y[i] + h * k1/2.0);
	k3 = G(x[i] + h/2.0, y[i] + h * k2/2.0);
	k4 = G(x[i] + h, y[i] + h*k3);
 
        y[i+1] = y[i] + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
 
        x[i+1] = x[i] + h;
    }
}

// Значения точного решения функции в узлах аппроксимации
void myfunc(double *x, double *u, int n)
{
    for (int i = 0; i < n; i++)
    {
	u[i] = exactSolution(x[i]);
    }
}

// Параметры для вывода графиков
#define NPOINTS 20
#define SECONDS 60

// Вывод графиков в gnuplot
void plot_graphs(double *x, double *u, int n)
{
    #if __unix__
    gnuplot_ctrl *h1, *h2;
    int i;

    h1 = gnuplot_init();
    h2 = gnuplot_init();
    
    gnuplot_set_axislabel(h1, "x", "X");
    gnuplot_set_axislabel(h1, "y", "Y");

    gnuplot_cmd(h1, "set key top left");
    gnuplot_cmd(h1, "set title 'Сравнение методов решения дифференциальных уравнений'");

    gnuplot_cmd(h1, "set xrange[0:0.2]");
    gnuplot_cmd(h1, "set yrange[1:5]");

    gnuplot_cmd(h1, "y(x) = exp(8*x)");

    gnuplot_cmd(h1, "set grid xtics ytics mxtics mytics");
    
    
    gnuplot_cmd(h1, "plot 'euler.csv' using 1:2 w l dt 4 title 'Метод Эйлера', \
    'runge2.csv' using 1:2 w l dt 2 title 'Рунге-Кутта 2-го порядка', \
    'runge3.csv' using 1:2 w l dt 3 title 'Рунге-Кутта 3-го порядка', \
    'runge4.csv' using 1:2 w l dt 4 title 'Рунге-Кутта 4-го порядка', \
    y(x) with lines title 'Точное решение'");


    gnuplot_set_axislabel(h2, "x", "X");
    gnuplot_set_axislabel(h2, "y", "Y");

    gnuplot_cmd(h2, "set grid xtics ytics mxtics mytics");

    gnuplot_cmd(h2, "set key bottom right");
    gnuplot_cmd(h2, "set title 'График погрешностей методов аппроксимации'");

    gnuplot_cmd(h2, "set xrange[-0.5:0.5]");
    gnuplot_cmd(h2, "set yrange[0:2.5]");

    gnuplot_cmd(h2, "y(x) = exp(8*x)");

    gnuplot_cmd(h2, "plot 'euler_rate.csv' using 1:2:3 with xerrorlines dt 5 title 'Метод Эйлера', \
    'runge2_rate.csv' using 1:2:3 with xerrorlines dt 5 title 'Рунге-Кутта 2-го порядка', \
    'runge3_rate.csv' using 1:2:3 with xerrorlines dt 5 title 'Рунге-Кутта 3-го порядка', \
    'runge4_rate.csv' using 1:2:3 with xerrorlines dt 5 title 'Рунге-Кутта 4-го порядка', \
    y(x) with lines title 'Точное решение'");

    sleep(SECONDS);
    gnuplot_close(h1);
    gnuplot_close(h2);
    #endif
}

// Записываем значения аппроксимации в файлы
void write_plots(double *x, double *u, int n, double h, double x0, double y0)
{
    // Открываем файл для записи координат графиков
    FILE *f = fopen("euler.csv", "wb");
    
    if (f == NULL)
	perror("Не открывается!");
    else
    {
	// Применяем метод Эйлера для решения задачи Коши
	EulerMethod(x, u, n, h, x0, y0);

	// Записываем значения сетки
	for (int i = 0; i < n; i++)
	{
	    fprintf(f, "%f, %f\n", x[i], u[i]);
	}

	fclose(f);
    }

    f = fopen("runge2.csv", "wb");

    if (f == NULL)
	perror("Не открывается!");
    else
    {
	// Применяем метод Рунге-Кутты 2-го порядка для решения задачи Коши
	rungeKutta2nd(x, u, n, h, x0, y0);

	// Записываем значения сетки
	for (int i = 0; i < n; i++)
	{
	    fprintf(f, "%f, %f\n", x[i], u[i]);
	}

	fclose(f);
    };

    f = fopen("runge3.csv", "wb");

    if (f == NULL)
	perror("Не открывается!");
    else
    {
	// Применяем метод Рунге-Кутты 3-го порядка для решения задачи Коши
	rungeKutta3nd(x, u, n, h, x0, y0);
	
	// Записываем значения сетки
	for (int i = 0; i < n; i++)
	{
	    fprintf(f, "%f, %f\n", x[i], u[i]);
	}

	fclose(f);
    };

    f = fopen("runge4.csv", "wb");

    if (f == NULL)
	perror("Не открывается!");
    else
    {
	// Применяем метод Рунге-Кутты 4-го порядка для решения задачи Коши
	rungeKutta4nd(x, u, n, h, x0, y0);
	
	// Записываем значения сетки
	for (int i = 0; i < n; i++)
	{
	    fprintf(f, "%f, %f\n", x[i], u[i]);
	}

	fclose(f);
    };

    f = fopen("myfunc.csv", "wb");

    if (f == NULL)
	perror("not opened!");
    else
    {
	// Считаем значения исходной функции в узлах шага 
	myfunc(x, u, n);
	
	// Записываем значения сетки
	for (int i = 0; i < n; i++)
	{
	    fprintf(f, "%f, %f\n", x[i], u[i]);
	}

	fclose(f);
    }; 
}

// Реализация правила Рунге для расчёта погрешностей аппроксимации для каждого из четырех методов
void error_rate(double *x, double *u, int n, double h, double x0, double y0, double epsilon)
{
    FILE *f = fopen("euler_rate.csv", "wb");
    double runge_rule = 0.0;
    double new_h = h;
    
    double h_2 = h / 2, u_h_2[n+1];

    // Расчёт погрешностей метода Рунге-Кутты 1-го порядка (метода Эйлера)
    EulerMethod(x, u, n, h, x0, y0);
    EulerMethod(x, u_h_2, n, h_2, x0, y0);

    int cnt = 0;
    for (int i = 0; i < n; i++)
    {
	double syg = ((u[i] - u_h_2[i]) / (pow(2, 1) - 1)) + pow(new_h, 2);
	cnt++;
	
	if (syg > epsilon)
	{
	    new_h /= 2;
	    
	    i = 0;
	}

	if (syg <= 10 * epsilon)
	    new_h *= 2;

	if (cnt == n)
	    break;

	printf("%f\n", syg);
	fprintf(f, "%f, %f, %f\n", x[cnt], exactSolution(x[cnt]), syg);
    }

    fclose(f);

    f = fopen("runge2_rate.csv", "wb");

    printf("\n");

    cnt = 0;
    new_h = h;

    // Расчёт погрешностей метода Рунге-Кутты 2-го порядка
    rungeKutta2nd(x, u, n, h, x0, y0);
    rungeKutta2nd(x, u_h_2, n, h_2, x0, y0);

    for (int i = 0; i < n; i++)
    {
	double syg = ((u[i] - u_h_2[i]) / (pow(2, 2) - 1)) + pow(new_h, 3);
	
	if (syg > epsilon)
	{
	    new_h /= 2;
	    i = 0;
	}

	if (syg * 10 <= epsilon)
	    new_h *= 2;

	if (cnt == n)
	    break;
	
	printf("%f\n", syg);
	fprintf(f, "%f, %f, %f\n", x[cnt], exactSolution(x[cnt]), syg);
	cnt++;
    }

    fclose(f);

    cnt = 0;
    
    f = fopen("runge3_rate.csv", "wb");

    printf("\n");

    new_h = h;

    // Расчёт погрешностей метода Рунге-Кутты 3-го порядка
    rungeKutta3nd(x, u, n, h, x0, y0);
    rungeKutta3nd(x, u_h_2, n, h_2, x0, y0);

    for (int i = 0; i < n; i++)
    {
	double syg = ((u[i] - u_h_2[i]) / (pow(2, 3) - 1)) + pow(new_h, 4);
	cnt++;
	if (syg > epsilon)
	{
	    new_h /= 2;
	    i = 0;
	}

	if (syg * 10 <= epsilon)
	    new_h *= 2;

	if (cnt == n)
	    break;

	printf("%f\n", syg);
	fprintf(f, "%f, %f, %f\n", x[cnt], exactSolution(x[cnt]), syg);
    }

    fclose(f);

    f = fopen("runge4_rate.csv", "wb");

    printf("\n");

    new_h = h;

    cnt = 0;

    // Расчёт погрешностей метода Рунге-Кутты 4-го порядка
    rungeKutta4nd(x, u, n, h, x0, y0);
    rungeKutta4nd(x, u_h_2, n, h_2, x0, y0);

    for (int i = 0; i < n; i++)
    {
	double syg = ((u[i] - u_h_2[i]) / (pow(2, 4) - 1)) + pow(new_h, 5);
	cnt++;
	if (syg > epsilon)
	{
	    new_h /= 2;
	    i = 0;
	}

	if (syg * 10 <= epsilon)
	    new_h *= 2;

	if (cnt == n)
	    break;
	
	printf("%f\n", syg);
	fprintf(f, "%f, %f, %f\n", x[cnt], exactSolution(x[cnt]), syg);
    }
    fclose(f);

    printf("\n");
}

int main()
{
    // Исходное дифференциальное уравнение u'x - 8u(x) = 0 при u(0) = 1
    // Его решением является функция u = e^(8x)

    // Задаем сетку для численного интегрирования
    // Здесь 1 - это y0 = b
    double x0 = 0.0, y0 = 0.2, end = 0.2;
    double h = 0.05;
    int n = y0/h + 1;
    
    double x[n+1], u[n+1];
    u[0] = 1.0;

    write_plots(x, u, n, h, x0, y0);
    plot_graphs(x, u, n);
    error_rate(x, u, n, h, x0, y0, 0.01);
}
