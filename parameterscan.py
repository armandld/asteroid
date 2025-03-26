import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import os

# Parameters
# TODO adapt to what you need (folder path executable input filename)
executable = 'Ex3_2024'  # Name of the executable (NB: .exe extension is required on Windows)
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


tFin = 63072000 
msol = 1.98892e30
mjup = 5.9736e24
msat = 1
a = 149598023e3
adapt = True
tol = 30
nsteps = 3000000
sampling = 1
nsel_physics=2
m1 = 1.98892e30
m2 = 5.9736e24

valeurs = lire_configuration()

def actualise_valeur():
    global tFin, msol, mjup, msat, a, adapt, tol, output, nsteps, sampling, nsel_physics, m1, m2
    tFin = float(valeurs.get("tFin"))
    msol = float(valeurs.get("msol"))
    mjup = float(valeurs.get("mjup"))
    msat = float(valeurs.get("msat"))
    a = float(valeurs.get("a"))
    adapt = bool(valeurs.get("adapt"))
    tol = float(valeurs.get("tol"))
    nsteps = float(valeurs.get("nsteps"))
    sampling = float(valeurs.get("sampling"))
    nsel_physics = float(valeurs.get("nsel_physics"))
    m1 = float(valeurs.get("m1"))
    m2 = float(valeurs.get("m2"))

def ecrire_valeur(nom,valeur):
    global valeurs
    valeurs[nom] = valeur
    ecrire_configuration(valeurs)
    actualise_valeur()

def lancer_simulation(theta0, output_file):
    ecrire_configuration({"theta0": theta0})
    cmd = f"./{executable} {input_filename} output={output_file}"
    subprocess.run(cmd, shell=True)

adapt = [True, False]  # Nombre de pas par période

paramstr = 'adapt'  # Paramètre à scanner
param = adapt


# Question 1
ecrire_valeur("adapt",True)
ecrire_valeur("nsteps",300000)
ecrire_valeur("mjup",0)
ecrire_valeur("tol",1e-11)

xjup = a*msol/(msol-mjup)
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
    #t = data[:, 0]
    #countsteps = data[:,1]
    x = data[:, 2]
    y = data[:, 3]
    #vx = data[:, 4]
    #vy = data[:, 5]
    #E = data[:, 6]
    coul = "blue"
    sun = plt.Circle((xsol, 0), 696340000, color='red')
    jupiter = plt.Circle((xjup, 0), 1737100, color='orange')
    ax.add_patch(sun)
    ax.add_patch(jupiter)
    if (adapt == False):
        coul = "red"
        #x = np.zeros(len(y))
    p = plt.plot(x, y,color = coul,linestyle='-',label=f"Adapt={adapt}")
    l.append(p)
    
    # Solution analytique


# Tracé de l'étude de convergence
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.grid(True, linestyle="--", alpha=0.3)
plt.legend()
plt.title("Trajectoire")

plt.show()