# Аппроксимация с помощью методов Рунге-Кутты 1-4 порядков
Задача аппроксимации состоит в построении приближенной (аппроксимирующей) функции, в целом наиболее близко проходящей около данных точек или около данной непрерывной функции.
Используя расчётные формулы для построения аппроксимирующих кривых получим численное решение задачи Коши для ДУ первого порядка:

![image](https://github.com/denadorite/runge_kutta_approximation/assets/52493471/dc64907d-ff8a-4dae-85a7-292ff2914840)

График погрешностей аппроксимации:

![image](https://github.com/denadorite/runge_kutta_approximation/assets/52493471/066d7923-0ee7-4416-a8ae-271eaa64e184)

# Инструкция по использованию:
Линковать исходный файл runge_kutta_approximation.c необходимо с файлами gnuplot UNIX API для библиотеки отрисовки графиков gnuplot: gnuplot_i.c; gnuplot_i.h, которые можно найти по следующей ссылке:

gnuplot C interface library: https://github.com/longradix/gnuplot_i

Отрисовка графиков для системы Windows не предусмотрена.
