#include "pendulo.h"

// Calcula as derivadas do sistema: theta' = omega, omega' = -(g/l)*sin(theta)
void derivadas_pendulo(double t, Estado *estado, Estado *derivada, Parametros *params) {
    derivada->theta = estado->omega;
    derivada->omega = -(params->g / params->l) * sin(estado->theta);
}

// Executa um passo do método Runge-Kutta de ordem 4 para sistemas
void rk4_passo_sistema(double t, double h, Estado *y, Estado *y_novo, Parametros *params) {
    Estado k1, k2, k3, k4;
    Estado temp;
    
    // k1 = f(t, y)
    derivadas_pendulo(t, y, &k1, params);
    
    // k2 = f(t + h/2, y + h*k1/2)
    temp.theta = y->theta + 0.5 * h * k1.theta;
    temp.omega = y->omega + 0.5 * h * k1.omega;
    derivadas_pendulo(t + 0.5*h, &temp, &k2, params);
    
    // k3 = f(t + h/2, y + h*k2/2)
    temp.theta = y->theta + 0.5 * h * k2.theta;
    temp.omega = y->omega + 0.5 * h * k2.omega;
    derivadas_pendulo(t + 0.5*h, &temp, &k3, params);
    
    // k4 = f(t + h, y + h*k3)
    temp.theta = y->theta + h * k3.theta;
    temp.omega = y->omega + h * k3.omega;
    derivadas_pendulo(t + h, &temp, &k4, params);
    
    // y_novo = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    y_novo->theta = y->theta + (h / 6.0) * (k1.theta + 2*k2.theta + 2*k3.theta + k4.theta);
    y_novo->omega = y->omega + (h / 6.0) * (k1.omega + 2*k2.omega + 2*k3.omega + k4.omega);
}

// Simula o pêndulo com passo fixo
Resultado* simular_pendulo_fixo(double theta0, double h, int n_periodos, Parametros *params) {
    Resultado *res = criar_resultado(100000);
    
    Estado estado_atual, estado_novo;
    estado_atual.theta = theta0;
    estado_atual.omega = 0.0;
    
    double t = 0.0;
    int inversoes = 0;
    double omega_anterior = 0.0;
    
    adicionar_ponto(res, t, estado_atual.theta, estado_atual.omega);
    
    // Simula até completar n_periodos (2*n_periodos inversões)
    while (inversoes < 2 * n_periodos) {
        rk4_passo_sistema(t, h, &estado_atual, &estado_novo, params);
        
        t += h;
        
        // Detecta inversão de velocidade
        if (omega_anterior * estado_novo.omega < 0) {
            inversoes++;
        }
        
        omega_anterior = estado_novo.omega;
        estado_atual = estado_novo;
        
        adicionar_ponto(res, t, estado_atual.theta, estado_atual.omega);
    }
    
    return res;
}

// Simula o pêndulo com passo adaptativo
Resultado* simular_pendulo_adaptativo(double theta0, double tol, int n_periodos, Parametros *params) {
    Resultado *res = criar_resultado(100000);
    
    Estado estado_atual, estado_novo, estado_meio, estado_duplo;
    estado_atual.theta = theta0;
    estado_atual.omega = 0.0;
    
    double t = 0.0;
    double h = 0.01;  // passo inicial
    int inversoes = 0;
    double omega_anterior = 0.0;
    
    adicionar_ponto(res, t, estado_atual.theta, estado_atual.omega);
    
    while (inversoes < 2 * n_periodos) {
        // Um passo com h
        rk4_passo_sistema(t, h, &estado_atual, &estado_novo, params);
        
        // Dois passos com h/2
        rk4_passo_sistema(t, h/2.0, &estado_atual, &estado_meio, params);
        rk4_passo_sistema(t + h/2.0, h/2.0, &estado_meio, &estado_duplo, params);
        
        // Estima o erro local
        double erro_theta = fabs(estado_duplo.theta - estado_novo.theta) / 15.0;
        double erro_omega = fabs(estado_duplo.omega - estado_novo.omega) / 15.0;
        double erro = fmax(erro_theta, erro_omega);
        
        if (erro < tol || erro < 1e-12) {
            // Aceita o passo
            estado_novo.theta = estado_duplo.theta + (estado_duplo.theta - estado_novo.theta) / 15.0;
            estado_novo.omega = estado_duplo.omega + (estado_duplo.omega - estado_novo.omega) / 15.0;
            
            t += h;
            
            // Detecta inversão
            if (omega_anterior * estado_novo.omega < 0) {
                inversoes++;
            }
            
            omega_anterior = estado_novo.omega;
            estado_atual = estado_novo;
            
            adicionar_ponto(res, t, estado_atual.theta, estado_atual.omega);
            
            // Aumenta o passo se possível
            if (erro > 0) {
                double fator = pow(tol / erro, 0.2);
                h = h * fmin(fator, 2.0);
            }
        } else {
            // Rejeita o passo e diminui h
            double fator = pow(tol / erro, 0.2);
            h = h * fmax(fator, 0.5);
        }
        
        // Limita o passo
        if (h > 0.1) h = 0.1;
        if (h < 1e-6) h = 1e-6;
    }
    
    return res;
}

// Calcula o período numericamente usando interpolação linear
double calcular_periodo_numerico(Resultado *res, int n_inversoes) {
    if (n_inversoes < 2) n_inversoes = 2;
    
    int inversoes_encontradas = 0;
    double tempos_inversao[100];
    
    for (int i = 1; i < res->n_pontos && inversoes_encontradas < 100; i++) {
        // Detecta mudança de sinal (inversão de velocidade)
        if (res->omega[i-1] * res->omega[i] < 0) {
            double t1 = res->tempo[i-1];
            double t2 = res->tempo[i];
            double v1 = res->omega[i-1];
            double v2 = res->omega[i];
            
            // Interpolação linear para encontrar t exato da inversão
            double t_inversao = t1 + fabs(v1) / (fabs(v1) + fabs(v2)) * (t2 - t1);
            
            tempos_inversao[inversoes_encontradas] = t_inversao;
            inversoes_encontradas++;
            
            if (inversoes_encontradas >= n_inversoes) {
                break;
            }
        }
    }
    
    if (inversoes_encontradas >= n_inversoes) {
        // Calcula o período médio usando n_inversoes
        double tempo_total = tempos_inversao[n_inversoes-1] - tempos_inversao[0];
        double n_meios_periodos = n_inversoes - 1;
        return 2.0 * tempo_total / n_meios_periodos;
    } else if (inversoes_encontradas >= 2) {
        // Se não chegou a n_inversoes, usa o que tem
        double tempo_total = tempos_inversao[inversoes_encontradas-1] - tempos_inversao[0];
        double n_meios_periodos = inversoes_encontradas - 1;
        return 2.0 * tempo_total / n_meios_periodos;
    }
    
    return 0.0;
}

// Calcula o período analítico (aproximação linear)
double periodo_analitico(Parametros *params) {
    return 2.0 * M_PI * sqrt(params->l / params->g);
}

// Calcula theta analítico (aproximação linear)
double theta_analitico(double t, double theta0, Parametros *params) {
    return theta0 * cos(sqrt(params->g / params->l) * t);
}

// Funções auxiliares para gerenciamento de resultados
Resultado* criar_resultado(int capacidade_inicial) {
    Resultado *res = (Resultado*)malloc(sizeof(Resultado));
    res->capacidade = capacidade_inicial;
    res->n_pontos = 0;
    res->tempo = (double*)malloc(capacidade_inicial * sizeof(double));
    res->theta = (double*)malloc(capacidade_inicial * sizeof(double));
    res->omega = (double*)malloc(capacidade_inicial * sizeof(double));
    return res;
}

void adicionar_ponto(Resultado *res, double t, double theta, double omega) {
    if (res->n_pontos >= res->capacidade) {
        res->capacidade *= 2;
        res->tempo = (double*)realloc(res->tempo, res->capacidade * sizeof(double));
        res->theta = (double*)realloc(res->theta, res->capacidade * sizeof(double));
        res->omega = (double*)realloc(res->omega, res->capacidade * sizeof(double));
    }
    
    res->tempo[res->n_pontos] = t;
    res->theta[res->n_pontos] = theta;
    res->omega[res->n_pontos] = omega;
    res->n_pontos++;
}

void liberar_resultado(Resultado *res) {
    free(res->tempo);
    free(res->theta);
    free(res->omega);
    free(res);
}

void salvar_dados(Resultado *res, const char *filename) {
    FILE *f = fopen(filename, "w");
    if (f == NULL) {
        printf("Erro ao abrir arquivo %s\n", filename);
        return;
    }
    
    fprintf(f, "tempo,theta,omega\n");
    for (int i = 0; i < res->n_pontos; i++) {
        fprintf(f, "%.10f,%.10f,%.10f\n", res->tempo[i], res->theta[i], res->omega[i]);
    }
    
    fclose(f);
}
