#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

// Fonction pour calculer le déterminant d'une matrice carrée
double determinant(double **mat, int n) {
    if (n == 1) {
        return mat[0][0];
    } else {
        double det = 0;
        for (int i = 0; i < n; i++) {
            double **submat = (double **)malloc((n - 1) * sizeof(double *));
            for (int j = 0; j < n - 1; j++) {
                submat[j] = (double *)malloc((n - 1) * sizeof(double));
                for (int k = 0; k < n - 1; k++) {
                    submat[j][k] = mat[j + 1][(k >= i) ? k + 1 : k];
                }
            }
            det += (i % 2 == 0 ? 1 : -1) * mat[0][i] * determinant(submat, n - 1);

            // Free the dynamically allocated memory for submat
            for (int j = 0; j < n - 1; j++) {
                free(submat[j]);
            }
            free(submat);
        }
        return det;
    }
}

// Fonction pour tester si A_k est inversible pour chaque k de 0 à n-1
int test_A_et_Ak(double **A, int n) {
    for (int k = 0; k < n; k++) {
        double **Ak = (double **)malloc((k + 1) * sizeof(double *));
        for (int i = 0; i <= k; i++) {
            Ak[i] = (double *)malloc((k + 1) * sizeof(double));
            for (int j = 0; j <= k; j++) {
                Ak[i][j] = A[i][j];
            }
        }

        if (determinant(Ak, k + 1) == 0) {
            printf("A%d n'est pas inversible, donc la decomposition LU n'est pas applicable\n", k + 1);
            
            // Free the allocated memory for Ak
            for (int i = 0; i <= k; i++) {
                free(Ak[i]);
            }
            free(Ak);

            return 0;
        }

        // Free the allocated memory for Ak
        for (int i = 0; i <= k; i++) {
            free(Ak[i]);
        }
        free(Ak);
    }
    return 1;
}


void lu_decomposition(double **A, double **L, double **U, int n)
{
    // Initialisation de L comme une matrice identité 
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                L[i][j] = 1.0;
            }
            else
            {
                L[i][j] = 0.0;
            }
        }
    }

    // Initialisation de U comme la matrice A d'origine 
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            U[i][j] = A[i][j];
        }
    }

    // Étapes de la décomposition LU 
    for (int k = 0; k < n - 1; k++)
    {
        for (int i = k + 1; i < n; i++)
        {
            L[i][k] = U[i][k] / U[k][k];
            for (int j = k; j < n; j++)
            {
                U[i][j] = U[i][j] - L[i][k] * U[k][j];
            }
        }
    }
}

void solve_lu(double **L, double **U, double *b, double *x, int n)
{
    double *y = (double *)malloc(n * sizeof(double));

    // Résolution de Ly = b
    for (int i = 0; i < n; i++)
    {
        double sum = 0.0;
        for (int j = 0; j < i; j++)
        {
            sum = sum + L[i][j] * y[j];
        }
            y[i] = (b[i] - sum) / L[i][i];
    }

    // Résolution de Ux = y
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++)
        {
            sum = sum + U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }

    free(y);
}


// Function to print a fractional representation of a double
void print_fraction(double value) {
    double intPart;
    double fracPart = modf(value, &intPart);
    int numerator = (int)round(fracPart * 1000); // Multiplying by 1000 to get 3 decimal places
    int denominator = 1000; // 3 decimal places

    if (numerator == 0) {
        printf("%.0lf", intPart);
    } else {
        printf("%d/%d", numerator, denominator);
    }
}

void print_matrix(double **matrix, int n, const char *name)
{
    printf("\nMatrice %s :\n", name);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            print_fraction(matrix[i][j]);
            printf("\t");
        }
        printf("\n");
    }
}

int main()
{
    int n;
    printf("Entrer la taille de la matrice : ");
    scanf("%d", &n);

    double **A, **L, **U;
    double *b, *x;

    // Allocation de mémoire pour A, L, U, b, x et y
    A = (double **)malloc(n * sizeof(double *));
    L = (double **)malloc(n * sizeof(double *));
    U = (double **)malloc(n * sizeof(double *));
    b = (double *)malloc(n * sizeof(double));
    x = (double *)malloc(n * sizeof(double));

    for (int i = 0; i < n; i++)
    {
        A[i] = (double *)malloc(n * sizeof(double));
        L[i] = (double *)malloc(n * sizeof(double));
        U[i] = (double *)malloc(n * sizeof(double));
    }

    // Entrer les valeurs de la matrice A
    printf("Entrer les valeurs de la matrice A :\n");
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("A[%d][%d] = ", i, j);
            scanf("%lf", &A[i][j]);
        }
    }
    printf("\nMatrice A :\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            print_fraction(A[i][j]);
            printf("\t"); // Afficher chaque élément suivi d'une tabulation
        }
        printf("\n"); // Nouvelle ligne après chaque ligne de la matrice
    }
    if (test_A_et_Ak(A, n)) {
        printf("Toutes les sous-matrices A_k sont inversibles.\n");
        // Entrer les valeurs du vecteur b
        printf("Entrer les valeurs du vecteur b :\n");
        for (int i = 0; i < n; i++)
        {
            printf("b[%d] = ", i);
            scanf("%lf", &b[i]);
        }

        // Résolution du système Ax = b
        lu_decomposition(A, L, U, n);

        // Afficher les matrices L, U, et les solutions
        print_matrix(L, n, "L");
        print_matrix(U, n, "U");

        solve_lu(L, U, b, x, n);

        // Afficher les solutions
        printf("\nSolutions du système :\n");
        for (int i = 0; i < n; i++)
        {
            printf("x[%d] = ",i);
            print_fraction(x[i]);
            printf("\n");
        }

    }
    else {
        printf("\n");
    }
    
    // Libérer la mémoire
    for (int i = 0; i < n; i++)
    {
        free(A[i]);
        free(L[i]);
        free(U[i]);
    }
    free(A);
    free(L);
    free(U);
    free(b);
    free(x);

    return 0;
}
