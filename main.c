#include "pendulo.h"

void imprimir_cabecalho() {
    printf("========================================\n");
    printf("   SIMULADOR DE PÊNDULO - INF1608\n");
    printf("========================================\n\n");
}

void analise_comparativa(double theta0, Parametros *params) {
    printf("\n--- ANÁLISE PARA θ₀ = %.4f rad (%.2f°) ---\n", theta0, theta0 * 180.0 / M_PI);
    
    // Período analítico
    double T_analitico = periodo_analitico(params);
    printf("\nPeríodo Analítico (linearizado): %.6f s\n", T_analitico);
    
    // Testa diferentes passos fixos
    double passos[] = {0.01, 0.001, 0.0001};
    int n_passos = 3;
    
    printf("\n%-15s %-15s %-15s %-15s\n", "Método", "Período (s)", "Nº Passos", "Erro (s)");
    printf("----------------------------------------------------------------\n");
    
    for (int i = 0; i < n_passos; i++) {
        double h = passos[i];
        clock_t inicio = clock();
        Resultado *res = simular_pendulo_fixo(theta0, h, 1, params);
        clock_t fim = clock();
        
        double T_numerico = calcular_periodo_numerico(res, 10);
        double erro = fabs(T_numerico - T_analitico);
        double tempo_exec = ((double)(fim - inicio)) / CLOCKS_PER_SEC;
        
        printf("h = %-9.4f   %-15.6f %-15d %-15.6f\n", 
               h, T_numerico, res->n_pontos, erro);
        
        liberar_resultado(res);
    }
    
    // Passo adaptativo
    clock_t inicio = clock();
    Resultado *res_adapt = simular_pendulo_adaptativo(theta0, 1e-5, 1, params);
    clock_t fim = clock();
    
    double T_adaptativo = calcular_periodo_numerico(res_adapt, 10);
    double erro_adapt = fabs(T_adaptativo - T_analitico);
    double tempo_adapt = ((double)(fim - inicio)) / CLOCKS_PER_SEC;
    
    printf("Adaptativo       %-15.6f %-15d %-15.6f\n", 
           T_adaptativo, res_adapt->n_pontos, erro_adapt);
    
    // Salva dados para plotagem
    char filename[100];
    sprintf(filename, "dados_theta_%.2f.csv", theta0);
    salvar_dados(res_adapt, filename);
    printf("\nDados salvos em: %s\n", filename);
    
    liberar_resultado(res_adapt);
}

void encontrar_angulo_maximo(Parametros *params) {
    printf("\n\n--- DETERMINAÇÃO DO ÂNGULO MÁXIMO ---\n");
    printf("(para erro < 0.001 s na fórmula simplificada)\n\n");
    
    double T_analitico_base = periodo_analitico(params);
    
    printf("%-15s %-15s %-15s\n", "θ₀ (rad)", "θ₀ (°)", "Erro (s)");
    printf("-----------------------------------------------\n");
    
    double theta_max = 0;
    for (double theta = 0.1; theta <= M_PI; theta += 0.05) {
        Resultado *res = simular_pendulo_adaptativo(theta, 1e-5, 1, params);
        double T_real = calcular_periodo_numerico(res, 10);
        double erro = fabs(T_real - T_analitico_base);
        
        printf("%-15.4f %-15.2f %-15.6f", theta, theta * 180.0 / M_PI, erro);
        
        if (erro < 0.001) {
            theta_max = theta;
            printf(" ✓\n");
        } else {
            printf("\n");
            liberar_resultado(res);
            break;
        }
        
        liberar_resultado(res);
    }
    
    printf("\nÂngulo máximo encontrado: %.4f rad (%.2f°)\n", 
           theta_max, theta_max * 180.0 / M_PI);
}

void teste_tempo_real(double theta0, Parametros *params) {
    printf("\n\n--- TESTE DE EXECUÇÃO EM TEMPO REAL ---\n");
    printf("Simulando 10 períodos...\n\n");
    
    // Passo fixo h = 0.01
    clock_t inicio = clock();
    Resultado *res1 = simular_pendulo_fixo(theta0, 0.01, 10, params);
    clock_t fim = clock();
    double tempo_exec1 = ((double)(fim - inicio)) / CLOCKS_PER_SEC;
    double tempo_fisico1 = res1->tempo[res1->n_pontos - 1];
    
    printf("Passo fixo h=0.01:\n");
    printf("  Tempo de execução: %.6f s\n", tempo_exec1);
    printf("  Tempo físico simulado: %.6f s\n", tempo_fisico1);
    printf("  Razão (exec/físico): %.4f\n", tempo_exec1 / tempo_fisico1);
    printf("  Tempo real? %s\n\n", tempo_exec1 < tempo_fisico1 ? "SIM ✓" : "NÃO");
    
    liberar_resultado(res1);
    
    // Passo adaptativo
    inicio = clock();
    Resultado *res2 = simular_pendulo_adaptativo(theta0, 1e-5, 10, params);
    fim = clock();
    double tempo_exec2 = ((double)(fim - inicio)) / CLOCKS_PER_SEC;
    double tempo_fisico2 = res2->tempo[res2->n_pontos - 1];
    
    printf("Passo adaptativo:\n");
    printf("  Tempo de execução: %.6f s\n", tempo_exec2);
    printf("  Tempo físico simulado: %.6f s\n", tempo_fisico2);
    printf("  Razão (exec/físico): %.4f\n", tempo_exec2 / tempo_fisico2);
    printf("  Tempo real? %s\n", tempo_exec2 < tempo_fisico2 ? "SIM ✓" : "NÃO");
    
    liberar_resultado(res2);
}

void gerar_dados_plotagem(Parametros *params) {
    printf("\n\n--- GERAÇÃO DE DADOS PARA PLOTAGEM ---\n");
    
    double angulos[] = {0.1, 0.3, 0.5, 1.0, 1.5, 2.0};
    int n_angulos = 6;
    
    for (int i = 0; i < n_angulos; i++) {
        double theta0 = angulos[i];
        
        // Simula com passo adaptativo
        Resultado *res = simular_pendulo_adaptativo(theta0, 1e-5, 1, params);
        
        // Salva dados numéricos
        char filename_num[100];
        sprintf(filename_num, "plot_numerico_%.2f.csv", theta0);
        
        FILE *f = fopen(filename_num, "w");
        fprintf(f, "tempo,theta_numerico,theta_analitico\n");
        
        for (int j = 0; j < res->n_pontos; j++) {
            double t = res->tempo[j];
            double theta_num = res->theta[j];
            double theta_ana = theta_analitico(t, theta0, params);
            fprintf(f, "%.10f,%.10f,%.10f\n", t, theta_num, theta_ana);
        }
        
        fclose(f);
        printf("Gerado: %s\n", filename_num);
        
        liberar_resultado(res);
    }
}

int main() {
    imprimir_cabecalho();
    
    // Parâmetros do pêndulo
    Parametros params;
    params.g = 9.81;  // m/s²
    params.l = 1.0;   // m
    
    printf("Parâmetros da simulação:\n");
    printf("  g = %.2f m/s²\n", params.g);
    printf("  l = %.2f m\n", params.l);
    printf("  Tolerância (passo adaptativo): ε = 10⁻⁵\n");
    
    // Análises para diferentes ângulos
    double angulos_teste[] = {0.1, 0.5, 1.0, 1.5, 2.0, 2.5};
    int n_angulos = 6;
    
    for (int i = 0; i < n_angulos; i++) {
        analise_comparativa(angulos_teste[i], &params);
    }
    
    // Encontra ângulo máximo
    encontrar_angulo_maximo(&params);
    
    // Teste de tempo real
    teste_tempo_real(0.5, &params);
    
    // Gera dados para plotagem
    gerar_dados_plotagem(&params);
    
    printf("\n========================================\n");
    printf("   SIMULAÇÃO CONCLUÍDA\n");
    printf("========================================\n");
    
    return 0;
}
