# ILP_2burning

Programa Linear Inteiro (PLI) para o **número de 2-queima** (*2-burning number*, $b_2(G)$) de grafos, implementado em Julia com [JuMP](https://jump.dev/) e [CPLEX](https://www.ibm.com/products/ilog-cplex-optimization-studio). O modelo recebe como solução inicial (*warm start*) o resultado de um algoritmo guloso, que também fornece um limitante superior para reduzir o horizonte de tempo.

## O processo de 2-queima

O processo de 2-queima é uma variante do processo de queima de grafos (*graph burning*) em que um vértice precisa de **dois** vizinhos queimados para queimar por propagação. O processo ocorre em rodadas $t = 1, 2, \dots$:

1. **Propagação:** todo vértice não queimado com pelo menos 2 vizinhos queimados (no estado da rodada $t-1$) passa a ficar queimado na rodada $t$.
2. **Nova fonte:** se ainda restarem vértices não queimados, escolhe-se um deles como nova fonte, que é queimado.

O processo termina na primeira rodada em que todos os vértices estão queimados. O número de 2-queima $b_2(G)$ é o menor número de rodadas necessário, minimizado sobre todas as escolhas de fontes.

## Modelo PLI

Variáveis (com $T$ o horizonte de tempo):

| Variável | Significado |
|---|---|
| `x[v,t]` | 1 se o vértice `v` está queimado até o instante `t` ($t = 0,\dots,T$) |
| `s[v,t]` | 1 se `v` é a fonte escolhida no instante `t` ($t = 1,\dots,T$) |
| `f[t]` | 1 se o processo termina exatamente no instante `t` |
| `z` | instante final (valor objetivo) |

O objetivo é minimizar `z`, sujeito às restrições R1–R11 (monotonicidade de `x`, no máximo uma fonte por instante, cada vértice é fonte no máximo uma vez, regras de propagação, instante final único e todos queimados no instante final). A descrição detalhada de cada restrição está comentada no código-fonte, em `build_2burning_model`.

### Algoritmo guloso e warm start

A cada rodada, o guloso escolhe como fonte o vértice não queimado que maximiza o tamanho do *fecho de propagação* (quantos vértices ficariam queimados se a propagação rodasse até estabilizar), com desempate pelo maior número de vizinhos ainda não queimados. O valor obtido, $z_g$, é um limitante superior válido para $b_2(G)$ e é usado de duas formas:

- **`WARM_START`**: a solução gulosa completa (`x`, `s`, `f`, `z`) é fornecida ao CPLEX como solução inicial (MIP start).
- **`SHRINK_HORIZON`**: o horizonte passa a ser $T = z_g$ em vez de $T = n$, o que reduz bastante o número de variáveis e restrições.

## Requisitos

- [Julia](https://julialang.org/downloads/) 1.9 ou superior
- Pacotes: `JuMP`, `Graphs`, `CPLEX`
- Uma instalação do **IBM ILOG CPLEX Optimization Studio** com licença válida (há licença acadêmica gratuita). O pacote `CPLEX.jl` precisa encontrar a biblioteca; veja as [instruções do CPLEX.jl](https://github.com/jump-dev/CPLEX.jl#installation) sobre a variável de ambiente `CPLEX_STUDIO_BINARIES`.

Instalação dos pacotes Julia:

```julia
using Pkg
Pkg.add(["JuMP", "Graphs", "CPLEX"])
```

## Uso

```bash
julia burning2_pli_warmstart.jl <dir_entrada> <dir_saida>
```

Exemplo (é o que o script `run.sh` executa):

```bash
julia burning2_pli_warmstart.jl input output
```

O programa percorre `<dir_entrada>` recursivamente, resolve cada arquivo `.txt` e grava **um CSV por grafo** em `<dir_saida>`, espelhando a estrutura de subpastas da entrada:

```
input/                      output/
├── hb/                     ├── hb/
│   └── can__96.txt         │   └── can__96.csv
└── paths/                  └── paths/
    ├── p5.txt                  ├── p5.csv
    ├── p6.txt                  ├── p6.csv
    └── p15.txt                 └── p15.csv
```

Todas as instâncias são resolvidas em uma única sessão Julia, o que evita pagar o tempo de carga e compilação do JuMP/CPLEX a cada grafo. Grafos que já possuem CSV completo (cabeçalho e pelo menos uma linha) são **pulados**, então é possível interromper e retomar a execução. Se a resolução de um grafo falhar, o erro é registrado e o lote continua com os demais.

## Formato de entrada

Arquivos `.txt` com uma aresta por linha, dada por dois inteiros separados por espaço:

```
0 1
1 2
2 3
3 4
```

- Os rótulos dos vértices podem ser quaisquer inteiros; internamente são normalizados para índices consecutivos.
- Laços e arestas repetidas são descartados (o grafo é simples e não direcionado).
- Linhas vazias e linhas iniciadas por `#` ou `%` são ignoradas, assim como colunas além da segunda.

## Formato de saída

Cada CSV tem o cabeçalho e uma linha:

| Coluna | Descrição |
|---|---|
| `graph` | nome do grafo (nome do arquivo sem extensão) |
| `greedy_burning_number` | valor obtido pelo algoritmo guloso |
| `pli_burning_number` | melhor valor encontrado pelo PLI (vazio se nenhuma solução viável foi encontrada) |
| `time` | tempo do solver em segundos (não inclui o tempo do guloso) |
| `#sources` | número de fontes na solução do PLI |
| `status` | status de terminação do solver (`OPTIMAL`, `TIME_LIMIT`, ...) |

Exemplo (`output/paths/p15.csv`):

```csv
graph,greedy_burning_number,pli_burning_number,time,#sources,status
p15,9,9,5.522,9,OPTIMAL
```

> Quando o `status` é `TIME_LIMIT`, o valor em `pli_burning_number` é a melhor solução encontrada dentro do limite de tempo e **não** é necessariamente ótimo. No exemplo incluído, `can__96` termina com `TIME_LIMIT`, com guloso 10 e PLI 8.

## Configuração

Os parâmetros ficam no início de `burning2_pli_warmstart.jl`:

| Constante | Padrão | Descrição |
|---|---|---|
| `TIME_LIMIT` | `600.0` | limite de tempo do PLI por grafo, em segundos |
| `WARM_START` | `true` | fornece a solução gulosa como solução inicial ao CPLEX |
| `SHRINK_HORIZON` | `true` | usa $T = z_g$ em vez de $T = n$ |
| `SILENT` | `true` | `false` exibe o log do CPLEX (útil para conferir o MIP start) |

Para rodar vários processos em paralelo, convém limitar as threads do CPLEX para que as instâncias não disputem CPU e os tempos fiquem comparáveis.

## Estrutura do repositório

```
ILP_2burning/
├── burning2_pli_warmstart.jl   # modelo PLI, guloso, leitura de grafos e execução em lote
├── run.sh                      # exemplo de execução
├── input/                      # grafos de exemplo (listas de arestas)
├── output/                     # CSVs gerados
├── LICENSE
└── README.md
```

## Licença

Distribuído sob a licença [GNU GPL v3](LICENSE).

## Como citar

Se você usar este código em um trabalho acadêmico, considere citar o repositório (adicione aqui a referência ou o DOI quando houver).
