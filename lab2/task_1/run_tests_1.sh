#!/bin/bash


THREADS_LIST="1 2 4 8"
POINTS_LIST="500 1000 2000"
EXECUTABLE="./task_1"
RESULTS_FILE="results.csv"


if [ ! -f "$EXECUTABLE" ]; then
    echo "Ошибка: Исполняемый файл '$EXECUTABLE' не найден."
    echo "Пожалуйста, скомпилируйте программу: gcc task_1.c -o task_1 -fopenmp -lm"
    exit 1
fi

echo "Начинаю замеры производительности..."

echo "threads,points,time" > $RESULTS_FILE

for points in $POINTS_LIST; do
    echo "----------------------------------------"
    echo "Тестирование для сетки: ${points}x${points}"
    
    for threads in $THREADS_LIST; do
        echo "  -> Запуск с $threads потоками..."
        
        output_time=$( $EXECUTABLE $threads $points | grep "Время вычислений" | grep -oP '[\d\.]+' )
        
        if [ -n "$output_time" ]; then
            echo "     Время: $output_time секунд"
            echo "$threads,$points,$output_time" >> $RESULTS_FILE
        else
            echo "     Не удалось измерить время."
        fi
    done
done

echo "----------------------------------------"
echo "Все тесты завершены. Результаты сохранены в $RESULTS_FILE"