"""
Modelo de Programação Linear Inteira (PLI) para o Número de 2-Queima
(2-Burning Number, b2(G)) via variável z que codifica o instante final
usando a variável indicadora f_t.

Variáveis:
  x[v,t]  -> 1 se o vértice v está queimado até o instante t (t = 0,...,T)
  s[v,t]  -> 1 se v é a fonte de queima escolhida no instante t (t = 1,...,T)
  f[t]    -> 1 se o processo termina exatamente no instante t
  z       -> instante final (valor objetivo)

Requer: JuMP, CPLEX, Graphs (todos instaláveis via Pkg).
"""

using JuMP
using CPLEX
using Graphs

# ---------------------------------------------------------------------
# Construção do modelo
# ---------------------------------------------------------------------

"""
    build_2burning_model(g::AbstractGraph, T::Int; silent::Bool=false)

Constrói o modelo PLI (R1-R15) para o grafo `g` com horizonte de tempo `T`.
`T` deve ser um limitante superior válido para o número de queima
(por exemplo, T = nv(g), ou um limitante mais justo baseado no diâmetro).

Retorna o `Model` do JuMP pronto para ser otimizado.
"""
function build_2burning_model(g::AbstractGraph, T::Int; silent::Bool=false)
    V = collect(vertices(g))

    model = Model(CPLEX.Optimizer)
    silent && set_silent(model)

    d(v) = degree(g, v)
    N(v) = neighbors(g, v)

    # ---------------- Variáveis ----------------
    @variable(model, x[v in V, t in 0:T], Bin)
    @variable(model, s[v in V, t in 1:T], Bin)
    @variable(model, f[t in 1:T], Bin)
    @variable(model, z >= 1, Int)

    # ---------------- R1: x_{v,0} = 0 ----------------
    @constraint(model, R1[v in V], x[v, 0] == 0)

    # ---------------- R2: x_{v,t} >= x_{v,t-1} (monotonicidade) ----------------
    @constraint(model, R2[v in V, t in 1:T], x[v, t] >= x[v, t-1])

    # ---------------- R3: x_{v,t} >= s_{v,t} ----------------
    @constraint(model, R3[v in V, t in 1:T], x[v, t] >= s[v, t])

    # ---------------- R4: cada vértice é fonte no máximo uma vez ----------------
    @constraint(model, R4[v in V], sum(s[v, t] for t in 1:T) <= 1)

    # ---------------- R5: no máximo uma fonte por instante ----------------
    @constraint(model, R5[t in 1:T], sum(s[v, t] for v in V) <= 1)

    # ---------------- R6: só pode ser fonte se ainda não queimado ----------------
    @constraint(model, R6[v in V, t in 1:T], s[v, t] <= 1 - x[v, t-1])

    # ---------------- R7: propagação (upper bound de queima em 2 saltos) --------
    @constraint(model, R7[v in V, t in 1:T],
        x[v, t] <= x[v, t-1] + s[v, t] + 0.5 * sum(x[u, t-1] for u in N(v)))

    # ---------------- R8: só p/ vértices com grau >= 2 ---------------------------
    @constraint(model, R8[v in V, t in 1:T; d(v) >= 2],
        0.5 + d(v) * x[v, t] >=
            0.5 * sum(x[u, t-1] for u in N(v)) - d(v) * (x[v, t-1] + s[v, t]))

    # ---------------- R9: exatamente um instante final ----------------
    @constraint(model, R9, sum(f[t] for t in 1:T) == 1)

    # ---------------- R10: z = instante final ----------------
    @constraint(model, R10, z == sum(t * f[t] for t in 1:T))

    # ---------------- R11: todo vértice deve estar queimado no instante final ---
    @constraint(model, R11[v in V, t in 1:T], x[v, t] >= f[t])

    # ---------------- Objetivo ----------------
    @objective(model, Min, z)

    return model
end

# ---------------------------------------------------------------------
# Extração da solução
# ---------------------------------------------------------------------

"""
    extract_solution(model, g)

Após `optimize!(model)`, extrai:
  - b2 : valor ótimo de z (número de 2-queima)
  - sources : vetor de (t, v) com as fontes escolhidas em ordem de t
"""
function extract_solution(model, g::AbstractGraph)
    b2 = round(Int, value(model[:z]))

    s = model[:s]
    sources = Tuple{Int,Int}[]
    for t in 1:last(axes(s, 2))
        for v in vertices(g)
            if value(s[v, t]) > 0.5
                push!(sources, (t, v))
            end
        end
    end
    sort!(sources, by = first)

    return b2, sources
end

# ---------------------------------------------------------------------
# Exemplo de uso
# ---------------------------------------------------------------------
#=
function example()
    # Grafo de exemplo: caminho P5 (troque por sua leitura via graph6)
    g = path_graph(5)

    n = nv(g)
    T = n  # limitante superior simples; pode ser refinado (ex.: baseado no diâmetro)

    model = build_2burning_model(g, T; silent = true)
    optimize!(model)

    status = termination_status(model)
    println("Status: ", status)

    if status == MOI.OPTIMAL || status == MOI.TIME_LIMIT
        b2, sources = extract_solution(model, g)
        println("Número de 2-queima (b2): ", b2)
        println("Sequência de fontes (t, v): ", sources)
    end

    return model
end

# Descomente para rodar o exemplo diretamente:
example()
=#