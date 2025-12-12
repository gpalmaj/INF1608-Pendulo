#!/usr/bin/env python3
"""
Script para plotar os resultados da simulação do pêndulo
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

def plotar_comparacao(theta0_graus):
    """Plota a comparação entre solução numérica e analítica"""
    
    filename = f"plot_numerico_{theta0_graus:.2f}.csv"
    
    if not os.path.exists(filename):
        print(f"Arquivo {filename} não encontrado")
        return
    
    # Lê os dados
    dados = pd.read_csv(filename)
    
    # Cria o gráfico
    plt.figure(figsize=(12, 6))
    
    plt.plot(dados['tempo'], dados['theta_numerico'], 'b-', 
             label='Solução Numérica (RK4)', linewidth=2)
    plt.plot(dados['tempo'], dados['theta_analitico'], 'r--', 
             label='Solução Analítica (linearizada)', linewidth=2, alpha=0.7)
    
    plt.xlabel('Tempo (s)', fontsize=12)
    plt.ylabel('θ (rad)', fontsize=12)
    plt.title(f'Comparação: Pêndulo com θ₀ = {theta0_graus:.2f} rad ({theta0_graus*180/np.pi:.1f}°)', 
              fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Salva o gráfico
    output_file = f"grafico_theta_{theta0_graus:.2f}.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Gráfico salvo: {output_file}")
    plt.close()

def plotar_todos():
    """Plota comparações para todos os ângulos"""
    
    # Lista todos os arquivos CSV gerados
    arquivos = glob.glob("plot_numerico_*.csv")
    
    if not arquivos:
        print("Nenhum arquivo de dados encontrado")
        return
    
    print(f"\nEncontrados {len(arquivos)} arquivos de dados")
    print("Gerando gráficos...\n")
    
    for arquivo in sorted(arquivos):
        # Extrai o valor de theta0 do nome do arquivo
        theta0_str = arquivo.replace("plot_numerico_", "").replace(".csv", "")
        theta0 = float(theta0_str)
        
        plotar_comparacao(theta0)
    
    print("\nTodos os gráficos foram gerados!")

def plotar_comparacao_multipla():
    """Plota múltiplos ângulos em um único gráfico"""
    
    arquivos = glob.glob("plot_numerico_*.csv")
    
    if not arquivos:
        print("Nenhum arquivo de dados encontrado")
        return
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    axes = axes.flatten()
    
    for idx, arquivo in enumerate(sorted(arquivos)[:6]):
        theta0_str = arquivo.replace("plot_numerico_", "").replace(".csv", "")
        theta0 = float(theta0_str)
        
        dados = pd.read_csv(arquivo)
        
        axes[idx].plot(dados['tempo'], dados['theta_numerico'], 'b-', 
                      label='Numérica', linewidth=2)
        axes[idx].plot(dados['tempo'], dados['theta_analitico'], 'r--', 
                      label='Analítica', linewidth=2, alpha=0.7)
        
        axes[idx].set_xlabel('Tempo (s)', fontsize=10)
        axes[idx].set_ylabel('θ (rad)', fontsize=10)
        axes[idx].set_title(f'θ₀ = {theta0:.2f} rad ({theta0*180/np.pi:.1f}°)', 
                           fontsize=11, fontweight='bold')
        axes[idx].legend(fontsize=9)
        axes[idx].grid(True, alpha=0.3)
    
    plt.suptitle('Comparação Numérico vs Analítico para Diferentes Ângulos Iniciais', 
                 fontsize=16, fontweight='bold', y=1.00)
    plt.tight_layout()
    
    output_file = "comparacao_multipla.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nGráfico comparativo salvo: {output_file}")
    plt.close()

def plotar_erro_periodo():
    """Plota o erro no período em função do ângulo inicial"""
    
    # Dados extraídos da saída do programa
    theta0_rad = np.array([0.1, 0.5, 1.0, 1.5, 2.0, 2.5])
    theta0_graus = theta0_rad * 180 / np.pi
    
    # Períodos numéricos (adaptativo)
    T_numerico = np.array([2.007264, 2.037868, 2.138830, 2.330940, 2.665632, 3.291926])
    
    # Período analítico
    T_analitico = 2.006067
    
    # Calcula erros
    erro_absoluto = np.abs(T_numerico - T_analitico)
    erro_relativo = (erro_absoluto / T_analitico) * 100
    
    # Cria o gráfico
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Erro absoluto
    ax1.plot(theta0_graus, erro_absoluto, 'bo-', linewidth=2, markersize=8)
    ax1.axhline(y=0.001, color='r', linestyle='--', label='Limite de erro (0.001 s)')
    ax1.set_xlabel('Ângulo Inicial (°)', fontsize=12)
    ax1.set_ylabel('Erro Absoluto (s)', fontsize=12)
    ax1.set_title('Erro no Período vs Ângulo Inicial', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    
    # Erro relativo
    ax2.plot(theta0_graus, erro_relativo, 'go-', linewidth=2, markersize=8)
    ax2.set_xlabel('Ângulo Inicial (°)', fontsize=12)
    ax2.set_ylabel('Erro Relativo (%)', fontsize=12)
    ax2.set_title('Erro Relativo no Período vs Ângulo Inicial', fontsize=13, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    output_file = "erro_periodo.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Gráfico de erro salvo: {output_file}")
    plt.close()

if __name__ == "__main__":
    print("=" * 60)
    print("  GERADOR DE GRÁFICOS - SIMULAÇÃO DO PÊNDULO")
    print("=" * 60)
    
    # Plota gráficos individuais
    plotar_todos()
    
    # Plota comparação múltipla
    plotar_comparacao_multipla()
    
    # Plota análise de erro
    plotar_erro_periodo()
    
    print("\n" + "=" * 60)
    print("  PROCESSAMENTO CONCLUÍDO")
    print("=" * 60)
