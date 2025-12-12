#ifndef PENDULO_H
#define PENDULO_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// Estrutura para armazenar o estado do pêndulo
typedef struct {
    double theta;  // ângulo
    double omega;  // velocidade angular
} Estado;

// Estrutura para armazenar resultados da simulação
typedef struct {
    double *tempo;
    double *theta;
    double *omega;
    int n_pontos;
    int capacidade;
} Resultado;

// Parâmetros do pêndulo
typedef struct {
    double g;  // aceleração da gravidade
    double l;  // comprimento do pêndulo
} Parametros;

// Funções principais
void derivadas_pendulo(double t, Estado *estado, Estado *derivada, Parametros *params);
void rk4_passo_sistema(double t, double h, Estado *y, Estado *y_novo, Parametros *params);
Resultado* simular_pendulo_fixo(double theta0, double h, int n_periodos, Parametros *params);
Resultado* simular_pendulo_adaptativo(double theta0, double tol, int n_periodos, Parametros *params);
double calcular_periodo_numerico(Resultado *res, int n_inversoes);
double periodo_analitico(Parametros *params);
double theta_analitico(double t, double theta0, Parametros *params);

// Funções auxiliares
Resultado* criar_resultado(int capacidade_inicial);
void adicionar_ponto(Resultado *res, double t, double theta, double omega);
void liberar_resultado(Resultado *res);
void salvar_dados(Resultado *res, const char *filename);

#endif
