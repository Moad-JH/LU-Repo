import numpy as np
import tkinter as tk
from tkinter import ttk

def principale(n, root):
    root.title("Décomposition LU")

    a_entries = [[None for _ in range(n)] for _ in range(n)]
    b_entries = [None for _ in range(n)]

    # Création des étiquettes pour la matrice A
    a_label = ttk.Label(root, text="Matrice A")
    a_label.grid(row=0, column=0, columnspan=n, padx=5, pady=5)

    # Création de la grille pour les entrées de la matrice A
    for i in range(n):
        for j in range(n):
            a_entries[i][j] = tk.Entry(root, width=10)
            a_entries[i][j].grid(row=i+1, column=j, padx=5, pady=5)

    # Création des étiquettes pour le vecteur B
    b_label = ttk.Label(root, text="Vecteur B")
    b_label.grid(row=0, column=n, padx=5, pady=5)

    # Création des entrées pour le vecteur B
    for i in range(n):
        b_entries[i] = tk.Entry(root, width=10)
        b_entries[i].grid(row=i+1, column=n, padx=5, pady=5)

    def calculer():
        A = np.array([[float(entry.get()) for entry in row] for row in a_entries])
        B = np.array([float(entry.get()) for entry in b_entries])

        if not test_A_et_Ak(A, n):
            return  # Arrêter si la matrice n'est pas inversible

        U, new_B, L = elimination_Gauss_et_remplissage_de_L(A, B, np.eye(n))
        Y = Calcule_Y(L, np.zeros(n), new_B)
        X = Calcule_X(U, np.zeros(n), Y)

        print("Matrice U :\n", U)
        print("Matrice L :\n", L)
        print("Vecteur Y :", Y)
        print("Vecteur X :", X)

    def test_A_et_Ak(A,n):
        for k in range(n):
            Ak=A[:k+1,:k+1]
            if np.linalg.det(Ak)==0:
                print("A"+str(k)+" n'est pas inversible, donc la decomposition LU n'est pas applicable")
                return False
        return True

    def elimination_Gauss_et_remplissage_de_L(A,B,In):
        for  i in range(n-1):
            for j in range(i+1, n):
                C=A[j,i]/A[i,i]
                A[j] = A[j] - C*A[i]
                B[j] = B[j] - C*B[i]
                In[i,j]=C
        return A, B, In

    def Calcule_Y(L,Y,new_B):
        Y[0]=new_B[0]/L[0,0]
        for i in range(1,n):
            sum_term=0
            for j in range(i):
                sum_term+=L[i,j]*Y[j]
                Y[i]=(new_B[i] - sum_term)/L[i,i]
        return Y

    def Calcule_X(U,X,Y):
        X[n-1]=Y[n-1]/U[n-1,n-1]
        for i in range(n-2,-1,-1):
            sum_term=0
            for j in range(i+1,n):
                sum_term+=U[i,j]*X[j]
            X[i]=(Y[i] - sum_term)/U[i,i]
        return X

    # Bouton pour déclencher le calcul
    calculate_button = ttk.Button(root, text="Calculer", command=calculer)
    calculate_button.grid(row=n+1, column=n, pady=10)

# Demander à l'utilisateur de saisir la taille n
n = int(input('Entrer la taille n :'))

# Créer la fenêtre principale Tkinter
root = tk.Tk()

# Appeler la fonction principale avec la taille n et la fenêtre Tkinter
principale(n, root)

# Lancer la boucle principale Tkinter
root.mainloop()
