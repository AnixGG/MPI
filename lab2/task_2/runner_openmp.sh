#!/bin/bash

PROGRAM="./task2_openmp"
ARGS="1000.0 input.txt"
RUNS=10

if [ ! -f "$PROGRAM" ]; then
    echo "Ошибка: Файл $PROGRAM не найден."
    exit 1
fi

total_time=0

echo "Запускаем тест ($RUNS раз)..."
echo "------------------------------"

for ((i=1; i<=RUNS; i++)); do
    current_time=$($PROGRAM $ARGS)

    echo "Run #$i: $current_time сек."

    total_time=$(echo "$total_time + $current_time" | bc -l)
done

avg_time=$(echo "scale=6; $total_time / $RUNS" | bc -l)

echo "------------------------------"
echo "Среднее время: $avg_time сек."