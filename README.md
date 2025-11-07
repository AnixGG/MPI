# Отчёт по лабораторной работе (MPI)

### Установка нужных библиотек
```bash
sudo apt update  
sudo apt install openmpi-bin openmpi-doc libopenmpi-dev python3-pandas python3-matplotlib python3-numpy -y
```
### Запуск экспериментов (пример)
в корне проекта:
```bash
./utils/run_experiments_1.sh
./utils/run_experiments_2.sh
./utils/run_experiments_3.sh
```

## Задание 1 — вычисление π методом Монте-Карло

**Условия**
- процессы: `1 4 8`  
- выборки: `10^6 15*10^5 10^7 15*10^6 2*10^7`  
- повторов на конфигурацию: `3`

**Графики**  
- Ускорение
  
![Иллюстрация к проекту](https://github.com/TypicalCode0/distributed_computing/blob/main/utils/graphics_1/speedup.png)    
- Эффективность (ускорение/количество процессоров)
  
![Иллюстрация к проекту](https://github.com/TypicalCode0/distributed_computing/blob/main/utils/graphics_1/efficiency.png)  

**Наблюдения**
- Ускорение при увеличении числа процессов растёт почти линейно
- Эффективность высока: при 4 процессах она почти равна 1, при 8 процессах в среднем 0.9
- Для больших объёмов выборкипараллельная реализация демонстрирует ближе к идеальному масштабированию
- Эффективность в одном измерении превышает 1, но это можно списать на нестабильность данных в силу небольшого количества данных


**Вывод**
- Метод Монте-Карло для оценки π хорошо распараллеливается, емееет высокие эффективность и ускорение

## Задание 2 — умножение матрицы на вектор

**Режимы тестирования**
- `row` (разбиение по строкам)  
- `col` (разбиение по столбцам)  
- `block` (блочное разбиение / сетка процессов)  
**Условия**
- процессы: `1 4 8`  
- выборки: `10^3 5*10^3 10^4`  
- повторов на конфигурацию: `3`

**Графики**  
- Ускорение

![Иллюстрация к проекту](https://github.com/TypicalCode0/distributed_computing/blob/main/utils/graphics_2/speedup_row.png)    
![Иллюстрация к проекту](https://github.com/TypicalCode0/distributed_computing/blob/main/utils/graphics_2/speedup_col.png)  
![Иллюстрация к проекту](https://github.com/TypicalCode0/distributed_computing/blob/main/utils/graphics_2/speedup_block.png)  
- Эффективность (ускорение/количество процессоров)
  
![Иллюстрация к проекту](https://github.com/TypicalCode0/distributed_computing/blob/main/utils/graphics_2/efficiency_row.png)  
![Иллюстрация к проекту](https://github.com/TypicalCode0/distributed_computing/blob/main/utils/graphics_2/efficiency_col.png)  
![Иллюстрация к проекту](https://github.com/TypicalCode0/distributed_computing/blob/main/utils/graphics_2/efficiency_block.png)  


**Наблюдения**
- ?

**Вывод**
- ?


## Задание 3 — умножение матриц (алгоритм Кэннона)

**Условия**
- число процессов: полный квадрат (`1 4 9`)  
- `n` должно делиться на `sqrt(p)`
- выборки: `6*10^2 12*10^2 15*10^2`
- повторов на конфигурацию: `3`

**Графики**  
- Ускорение
    
![Иллюстрация к проекту](https://github.com/TypicalCode0/distributed_computing/blob/main/utils/graphics_3/speedup.png)    
- Эффективность (ускорение/количество процессоров)
  
![Иллюстрация к проекту](https://github.com/TypicalCode0/distributed_computing/blob/main/utils/graphics_3/efficiency.png)  

**Наблюдения**
- ?

**Вывод**
- ?


## Итоговые выводы

- 

