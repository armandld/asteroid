import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import os

# Parameters
# TODO adapt to what you need (folder path executable input filename)
executable = 'Ex3_2025'  # Name of the executable (NB: .exe extension is required on Windows)
repertoire = r"/Users/Sayu/Desktop/asteroid"

os.chdir(repertoire)


CONFIG_FILE = os.path.join(os.path.dirname(__file__), "configuration.in.example")

input_filename = 'configuration.in.example'  # Name of the input file

outputs = []
errors = []

def lire_configuration():
    config_path = os.path.join(os.path.dirname(__file__), "configuration.in.example")
    configuration = {}
    
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Le fichier {config_path} n'existe pas.")
    
    with open(config_path, "r", encoding="utf-8") as fichier:
        for ligne in fichier:
            ligne = ligne.strip()
            if ligne and "=" in ligne and not ligne.startswith("#"):
                cle, valeur = ligne.split("=", 1)
                configuration[cle.strip()] = valeur.strip()
    
    return configuration

def ecrire_configuration(nouvelles_valeurs):
    """Écrit les nouvelles valeurs dans le fichier de configuration."""
    if not os.path.exists(CONFIG_FILE):
        raise FileNotFoundError(f"Le fichier {CONFIG_FILE} n'existe pas.")

    lignes_modifiees = []
    
    with open(CONFIG_FILE, "r", encoding="utf-8") as fichier:
        for ligne in fichier:
            ligne_strippée = ligne.strip()
            if ligne_strippée and "=" in ligne_strippée and not ligne_strippée.startswith("#"):
                cle, _ = ligne_strippée.split("=", 1)
                cle = cle.strip()
                if cle in nouvelles_valeurs:
                    ligne = f"{cle} = {nouvelles_valeurs[cle]}\n"
            lignes_modifiees.append(ligne)

    with open(CONFIG_FILE, "w", encoding="utf-8") as fichier:
        fichier.writelines(lignes_modifiees)


tFin = 3153600000
m1 = 1.98892e30
m2 = 0
msat = 1
a = 778.479e9
adapt = 1
tol = 1
nsteps = 300000
sampling = 1
nsel_physics = 1

valeurs = lire_configuration()

def actualise_valeur():
    global tFin, m1, m2, msat, a, adapt, tol, nsteps, sampling, nsel_physics
    tFin = float(valeurs.get("tFin"))
    m1 = float(valeurs.get("m1"))
    m2 = float(valeurs.get("m2"))
    msat = float(valeurs.get("msat"))
    a = float(valeurs.get("a"))
    adapt = bool(valeurs.get("adapt"))
    tol = float(valeurs.get("tol"))
    nsteps = float(valeurs.get("nsteps"))
    sampling = float(valeurs.get("sampling"))
    nsel_physics = int(valeurs.get("nsel_physics"))


def ecrire_valeur(nom,valeur):
    global valeurs
    valeurs[nom] = valeur
    ecrire_configuration(valeurs)
    actualise_valeur()

def lancer_simulation(theta0, output_file):
    ecrire_configuration({"theta0": theta0})
    cmd = f"./{executable} {input_filename} output={output_file}"
    subprocess.run(cmd, shell=True)


adapt = [0,1]  # Nombre de pas par période

paramstr = 'adapt'  # Paramètre à scanner
param = adapt


# Question 1
xjup = a*m1/(m1-m2)
xsol = a - xjup

l = []
plt.figure()
ax = plt.gca()
for i, adapt in enumerate(param):
    output_file = f"{paramstr}={adapt}.out"
    outputs.append(output_file)
    cmd = f"./{executable} {input_filename} {paramstr}={adapt} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Simulation terminée.')

    # Chargement des données
    data = np.loadtxt(output_file)
    t = data[:, 0]
    x = data[:, 1]
    y = data[:, 2]
    vx = data[:, 3]
    vy = data[:, 4]
    E = data[:, 5]

    coul = "blue"
    sun = plt.Circle((xsol, 0), 696340000, color='red')
    jupiter = plt.Circle((xjup, 0), 1737100, color='orange')
    ax.add_patch(sun)
    ax.add_patch(jupiter)
    q = '--'
    if (adapt == True):
        coul = "red"
        #q = '-'
        dt = data[:, 6]
        #x = np.zeros(len(y))
    p = plt.plot(x, y,color = coul,linestyle=q,label=f"Adapt={adapt}")
    l.append(p)
    
    # Solution analytique


# Tracé de l'étude de convergence
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.grid(True, linestyle="--", alpha=0.3)
plt.legend()
plt.title("Trajectoire")


plt.figure()
plt.xlabel("t [s]")
plt.ylabel("E [J]")
plt.plot(t, E,color = coul,linestyle='-',label=f"Adapt={adapt}")
plt.title("energy")

plt.figure()
plt.xlabel("t [s]")
plt.ylabel("d [m]")
plt.plot(t, np.sqrt(x*x + y*y),color = coul,linestyle='-',label=f"Adapt={adapt}")
plt.axhline(np.max(np.sqrt(x*x + y*y)), color='black', linestyle='--',  label=f"Max Distance={np.max(np.sqrt(x*x + y*y))}")
plt.axhline(np.min(np.sqrt(x*x + y*y)), color='purple', linestyle='--', label=f"Min Distance={np.min(np.sqrt(x*x + y*y))}")
plt.title("Distance par rapport au Soleil en fonction du temps")
plt.legend(loc = 'center')

plt.figure()
plt.xlabel("t [s]")
plt.ylabel("v [m/s]")
plt.plot(t, np.sqrt(vx*vx + vy*vy),color = coul,linestyle='-',label=f"Adapt={adapt}")
plt.axhline(np.max(np.sqrt(vx*vx + vy*vy)), color='black', linestyle='--', label=f"Max Velocity={np.max(np.sqrt(vx*vx + vy*vy))}")
plt.axhline(np.min(np.sqrt(vx*vx + vy*vy)), color='purple', linestyle='--',label=f"Min Velocity={np.min(np.sqrt(vx*vx + vy*vy))}")
plt.title("Vitesse en fonction du temps")
plt.legend()

plt.figure()
plt.xlabel("t [s]")
plt.ylabel("dt [s]")
plt.plot(t, dt,color = coul,linestyle='-',label=f"Adapt={adapt}")
plt.title("Pas de temps en fonction du temps")
plt.legend()



adapt = [0,1]
paramstr = 'adapt'  # Paramètre à scanner
param = adapt

l = []
plt.figure()
ax = plt.gca()
for i, adapt in enumerate(param):
    output_file = f"{paramstr}={adapt}.out"
    outputs.append(output_file)
    cmd = f"./{executable} {input_filename} {paramstr}={adapt} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Simulation terminée.')

    # Chargement des données
    data = np.loadtxt(output_file)
    t = data[:, 0]
    #countsteps = data[:,1]
    x = data[:, 1]
    y = data[:, 2]
    vx = data[:, 3]
    vy = data[:, 4]
    E = data[:, 5]
    coul = "blue"
    q = '-'
    if (adapt == True):
        coul = "red"
        #q = '-'
        #x = np.zeros(len(y))
    plt.plot(t, E,color = coul,linestyle=q,label=f"Adapt={adapt}")
    
    
    # Solution analytique
# Tracé de l'étude de convergence
plt.xlabel("t [s]")
plt.ylabel("E [J]")
plt.grid(True, linestyle="--", alpha=0.3)
plt.legend()
plt.title("Énergie")


ecrire_valeur("adapt",1)
tol=[1e7, 1e6, 1e5, 1e4,1e3, 1e2, 1e1, 1e0, 1e-1,1e-2]  # Nombre de pas par période

paramstr = 'tol'  # Paramètre à scanner
param = tol

l = []
K = []
j = []
plt.figure()
ax = plt.gca()
for i, tol0 in enumerate(param):
    output_file = f"{paramstr}={tol0}.out"
    outputs.append(output_file)
    cmd = f"./{executable} {input_filename} {paramstr}={tol0} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Simulation terminée.')

    # Chargement des données
    data = np.loadtxt(output_file)
    t = data[:, 0]
    E = data[:, 5]
    Emin = np.min(E)
    Emax = np.max(E)
    coul = "blue"
    q = '-'
    if (adapt == True):
        coul = "red"
    K.append(Emax-Emin)
    j.append(len(t))
    # Solution analytique
# Tracé de l'étude de convergence
plt.loglog(j, K,color = coul,linestyle=q)
plt.xlabel("nsteps")
plt.ylabel("$\\Delta E$ [J]")
plt.grid(True, linestyle="--", alpha=0.3)
plt.legend()
plt.title("Convergence énergie")


tol=[1e2, 1e1, 1e0, 1e-1,1e-2]  # Nombre de pas par période

l = []
K = []
j = []
plt.figure()
ax = plt.gca()
for i, tol0 in enumerate(param):
    output_file = f"{paramstr}={tol0}.out"
    outputs.append(output_file)
    cmd = f"./{executable} {input_filename} {paramstr}={tol0} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Simulation terminée.')

    # Chargement des données
    data = np.loadtxt(output_file)
    t = data[:, 0]
    x = data[:,1]
    y = data[:,2]
    r = x[-1]
    E = data[:, 5]
    coul = "blue"
    q = '-'
    if (adapt == True):
        coul = "red"
    K.append(r)

    j.append(len(t))
    # Solution analytique
# Tracé de l'étude de convergence
j = np.array(j)
j = 1/j
plt.loglog(j, K,color = coul,linestyle=q)
plt.xlabel("nsteps")
plt.ylabel("$\\Delta x$ [J]")
plt.grid(True, linestyle="--", alpha=0.3)
plt.legend()
plt.title("Convergence position finale")

plt.show()