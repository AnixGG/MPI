import subprocess
import matplotlib.pyplot as plt
import sys
import os
import re

SOURCE_FILE = "main.c"
OUTPUT_STD = "bench_std"
OUTPUT_MY = "bench_my"

PARAMS = {
    "initial_inserts": 1000,
    "total_ops": 100000,
    "search_percent": 0.80,
    "insert_percent": 0.10
}

THREADS_LIST = [1, 2, 4, 8, 12, 16]

def compile_programs():
    print("=== Компиляция ===")
    
    cmd_std = ["gcc", "-O2", "-pthread", SOURCE_FILE, "-o", OUTPUT_STD]
    print(f"Компиляция Standard: {' '.join(cmd_std)}")
    if subprocess.call(cmd_std) != 0:
        print("Ошибка компиляции стандартной версии!")
        sys.exit(1)

    cmd_my = ["gcc", "-O2", "-pthread", SOURCE_FILE, "-o", OUTPUT_MY, "-DUSE_MY_RWLOCK"]
    print(f"Компиляция Custom:   {' '.join(cmd_my)}")
    if subprocess.call(cmd_my) != 0:
        print("Ошибка компиляции моей версии!")
        sys.exit(1)
        
    print("Компиляция успешна.\n")

def run_benchmark(executable, threads):
    input_str = f"{PARAMS['initial_inserts']}\n{PARAMS['total_ops']}\n{PARAMS['search_percent']}\n{PARAMS['insert_percent']}\n"
    
    try:
        # Запуск процесса
        process = subprocess.Popen(
            [f"./{executable}", str(threads)],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        stdout, stderr = process.communicate(input=input_str)
        
        if process.returncode != 0:
            print(f"Ошибка выполнения {executable}: {stderr}")
            return 0.0

        match = re.search(r"Elapsed time\s*=\s*([0-9\.e\-\+]+)", stdout)
        if match:
            return float(match.group(1))
        else:
            print(f"Не удалось найти время в выводе {executable}")
            print("Вывод программы:\n", stdout)
            return 0.0

    except Exception as e:
        print(f"Exception при запуске {executable}: {e}")
        return 0.0

def main():
    if not os.path.exists(SOURCE_FILE):
        print(f"Файл {SOURCE_FILE} не найден.")
        return

    compile_programs()

    times_std = []
    times_my = []

    print(f"Запуск тестов (Total Ops: {PARAMS['total_ops']}, Read: {PARAMS['search_percent']*100}%)...")
    print(f"{'Threads':<10} | {'Std Time (s)':<15} | {'My Time (s)':<15}")
    print("-" * 45)

    for t in THREADS_LIST:
        t_std = run_benchmark(OUTPUT_STD, t)
        times_std.append(t_std)

        t_my = run_benchmark(OUTPUT_MY, t)
        times_my.append(t_my)

        print(f"{t:<10} | {t_std:<15.6f} | {t_my:<15.6f}")

    plt.figure(figsize=(10, 6))
    
    plt.plot(THREADS_LIST, times_std, marker='o', label='Pthread Library RWLock', linewidth=2)
    plt.plot(THREADS_LIST, times_my, marker='s', label='Custom RWLock (Ваш код)', linewidth=2, linestyle='--')

    plt.title(f'Сравнение RWLock\nOps: {PARAMS["total_ops"]}, Read: {PARAMS["search_percent"]*100}%, Write: {(1.0-PARAMS["search_percent"])*100:.0f}%')
    plt.xlabel('Количество потоков')
    plt.ylabel('Время выполнения (сек)')
    plt.grid(True, linestyle=':', alpha=0.6)
    plt.legend()
    
    output_png = "rwlock_comparison.png"
    plt.savefig(output_png)
    print(f"\nГрафик сохранен в файл: {output_png}")
    plt.show()

    if os.path.exists(OUTPUT_STD): os.remove(OUTPUT_STD)
    if os.path.exists(OUTPUT_MY): os.remove(OUTPUT_MY)

if __name__ == "__main__":
    main()