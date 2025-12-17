import random

counts = [16, 64, 128, 256, 512, 1024, 2048, 4096, 8192]

for n in counts:
    filename = f"input_{n}.txt"
    with open(filename, "w") as f:
        f.write(f"{n}\n")

        for _ in range(n):
            m = random.uniform(1000.0, 100000.0)

            x = random.uniform(-100.0, 100.0)
            y = random.uniform(-100.0, 100.0)
            z = random.uniform(-100.0, 100.0)

            vx = random.uniform(-2.0, 2.0)
            vy = random.uniform(-2.0, 2.0)
            vz = random.uniform(-2.0, 2.0)

            f.write(f"{m:.2f} {x:.2f} {y:.2f} {z:.2f} {vx:.2f} {vy:.2f} {vz:.2f}\n")

    print(f"Создан файл: {filename}")