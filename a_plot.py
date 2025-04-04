import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
# === 1. Compilation du programme C++ ===
cpp_file = "Ex3_2025.cpp"
executable = "./Ex3_2025"

if not os.path.exists(executable):
    print("Compilation du programme C++...")
    compile_cmd = ["g++", "-O2", cpp_file, "-o", executable, "-std=c++17"]
    result = subprocess.run(compile_cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(" Erreur de compilation :")
        print(result.stderr)
        exit(1)
    else:
        print(" Compilation réussie.")

# === 2. Exécution du programme ===
print("Exécution du programme C++...")
run_cmd = [executable, "configuration.in.example"]
result = subprocess.run(run_cmd)
if result.returncode != 0:
    print("Erreur à l'exécution.")
    exit(1)
else:
    print("Simulation terminée.")

# === 3. Chargement des données ===
def load_data(filename):
    data = np.loadtxt(filename)
    t = data[:, 0]
    x = data[:, 1]
    y = data[:, 2]
    vx = data[:, 3]
    vy = data[:, 4]
    E = data[:, 5]
    return t, x, y, vx, vy, E

filename = "Ex3.txt"
if not os.path.exists(filename):
    print(f"Fichier de sortie {filename} introuvable.")
    exit(1)

t, x, y, vx, vy, E = load_data(filename)

# === 4. Tracer la trajectoire ===
def plot_trajectory(x, y):
    plt.figure()
    ax = plt.gca()
    sun = patches.Circle((0, 0), 6.957e8, color='black', alpha=0.5, label="Soleil")
    ax.add_patch(sun)
    plt.plot(x, y, label="Astéroïde", color="blue")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.axis('equal')
    plt.legend()
    plt.title("Trajectoire de l'astéroïde")
    plt.grid()
    plt.show()

# === 5. Tracer l'énergie ===
def plot_energy(t, E):
    plt.figure(figsize=(8, 6))
    plt.plot(t, E, label="Énergie mécanique", color="red")
    plt.xlabel("Temps (s)")
    plt.ylabel("Énergie (J)")
    plt.title("Évolution de l'énergie mécanique")
    plt.legend()
    plt.grid()
    plt.show()

# === 6. Affichage final ===
plot_trajectory(x, y)
plot_energy(t, E)
print("vmin : ", np.min(np.sqrt(vx*vx+vy*vy)))
print("vmax : ", np.max(np.sqrt(vx*vx+vy*vy)))
print("xmin : ", np.min(np.sqrt(x*x+y*y)))
print("xmax : ", np.max(np.sqrt(x*x+y*y)))
