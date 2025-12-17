import pandas as pd
import matplotlib.pyplot as plt

RESULTS_FILE = 'results.csv'

try:
    df = pd.read_csv(RESULTS_FILE)
except FileNotFoundError:
    print(f"Ошибка: Файл '{RESULTS_FILE}' не найден. Сначала запустите скрипт run_tests.sh")
    exit()

points_sizes = df['points'].unique()

for size in points_sizes:
    print(f"Создание графиков для сетки {size}x{size}...")
    
    subset = df[df['points'] == size].copy()
    
    try:
        t1 = subset[subset['threads'] == 1]['time'].iloc[0]
    except IndexError:
        print(f"Предупреждение: Не найдены данные для 1 потока для сетки {size}x{size}. Пропускаю этот график.")
        continue

    subset['speedup'] = t1 / subset['time']
    subset['efficiency'] = subset['speedup'] / subset['threads']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle(f'Анализ производительности для сетки {size}x{size}', fontsize=16)

    ax1.plot(subset['threads'], subset['speedup'], marker='o', linestyle='-', label='Практическое ускорение')
    ax1.plot(subset['threads'], subset['threads'], linestyle='--', color='gray', label='Идеальное ускорение')
    ax1.set_xlabel('Количество потоков')
    ax1.set_ylabel('Ускорение (раз)')
    ax1.set_title('График ускорения')
    ax1.legend()
    ax1.grid(True)
    ax1.set_xticks(subset['threads'])


    ax2.plot(subset['threads'], subset['efficiency'], marker='o', linestyle='-', color='r', label='Практическая эффективность')
    ax2.axhline(y=1.0, linestyle='--', color='gray', label='Идеальная эффективность (100%)')
    ax2.set_xlabel('Количество потоков')
    ax2.set_ylabel('Эффективность')
    ax2.set_title('График эффективности')
    ax2.legend()
    ax2.grid(True)
    ax2.set_ylim(0, 1.2)
    ax2.set_xticks(subset['threads'])
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    output_filename = f'performance_{size}x{size}.png'
    plt.savefig(output_filename)
    print(f"Графики сохранены в файл: {output_filename}")

