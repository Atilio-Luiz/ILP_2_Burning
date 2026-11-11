"""
Modelo de Programação Linear Inteira (PLI) para o Número de 2-Queima
(2-Burning Number, b2(G)) via variável z que codifica o instante final
usando a variável indicadora f_t, com WARM START vindo de um algoritmo guloso.

Variáveis:
  x[v,t]  -> 1 se o vértice v está queimado até o instante t (t = 0,...,T)
  s[v,t]  -> 1 se v é a fonte de queima escolhida no instante t (t = 1,...,T)
  f[t]    -> 1 se o processo termina exatamente no instante t
  z       -> instante final (valor objetivo)

Requer: JuMP, CPLEX, Graphs (todos instaláveis via Pkg).
"""

using JuMP
using Graphs
using CPLEX

# ---------------------------------------------------------------------
# Parâmetros de execução (ajuste aqui)
# ---------------------------------------------------------------------
const TIME_LIMIT     = 600.0  # limite de tempo do PLI (segundos)
const WARM_START     = true   # fornece a solução gulosa como solução inicial ao CPLEX
const SHRINK_HORIZON = true   # usa T = valor guloso (limitante superior válido) em vez de T = n
const SILENT         = true   # false mostra o log do CPLEX (útil p/ conferir o MIP start)


# ---------------------------------------------------------------------
# Função que lê um grafo simples de um arquivo texto no formato de lista de arestas 
# ---------------------------------------------------------------------
function read_simple_graph(filepath::String)::SimpleGraph
    unique_edges = Set{Tuple{Int,Int}}()
    for line in eachline(filepath)
        line = strip(line)
        (isempty(line) || startswith(line, '#') || startswith(line, '%')) && continue
        parts = split(line)
        length(parts) < 2 && continue
        u, v = parse(Int, parts[1]), parse(Int, parts[2])
        u != v && push!(unique_edges, u < v ? (u, v) : (v, u))
    end

    all_vertices = sort(unique(vcat([[u, v] for (u, v) in unique_edges]...)))
    vertex_map = Dict(v => i for (i, v) in enumerate(all_vertices))

    g = SimpleGraph(length(all_vertices))
    for (u, v) in unique_edges
        add_edge!(g, vertex_map[u], vertex_map[v])
    end
    return g
end


# ---------------------------------------------------------------------
# Algoritmo guloso para o 2-queima
# ---------------------------------------------------------------------
#
# Regra do processo (a mesma codificada nas restrições R1-R8):
#   * Rodada t: primeiro ocorre a propagação a partir do estado em t-1
#     (um vértice não queimado queima se tiver >= 2 vizinhos queimados);
#     depois, se ainda restar vértice não queimado, escolhe-se uma nova
#     fonte entre os não queimados.
#   * O processo termina na primeira rodada t em que todos estão queimados.
#
# Escolha gulosa da fonte: o vértice não queimado que maximiza o tamanho do
# "fecho" de propagação (quantos vértices ficam queimados se, a partir dele,
# a propagação rodasse até estabilizar, sem novas fontes). Desempate: maior
# número de vizinhos ainda não queimados.

"Um passo de propagação: retorna o novo estado de queima a partir de `burned`."
function propagate_step(g::AbstractGraph, burned::AbstractVector{Bool})
    new = copy(burned)
    for v in vertices(g)
        burned[v] && continue
        c = 0
        for u in neighbors(g, v)
            if burned[u]
                c += 1
                c >= 2 && break
            end
        end
        c >= 2 && (new[v] = true)
    end
    return new
end

"""
Quantidade de vértices queimados após adicionar `u` como fonte e deixar a
propagação estabilizar. `base_cnt[v]` = nº de vizinhos queimados de v em `burned`.
"""
function closure_size(g::AbstractGraph, burned::AbstractVector{Bool},
                      base_cnt::Vector{Int}, u::Int)
    b = copy(burned)
    cnt = copy(base_cnt)
    b[u] = true
    total = count(b)
    queue = Int[u]
    while !isempty(queue)
        w = pop!(queue)
        for y in neighbors(g, w)
            cnt[y] += 1
            if !b[y] && cnt[y] >= 2
                b[y] = true
                total += 1
                push!(queue, y)
            end
        end
    end
    return total
end

"""
    greedy_2burning(g)

Retorna `(z, sources, states)`:
  - z       : instante final obtido (limitante superior para b2(g))
  - sources : vetor de (t, v) com a fonte escolhida em cada rodada t
  - states  : states[t+1][v] == true se v está queimado no instante t (t = 0,...,z)
"""
function greedy_2burning(g::AbstractGraph)
    n = nv(g)
    burned = falses(n)
    states = [copy(burned)]               # instante t = 0
    sources = Tuple{Int,Int}[]
    t = 0

    while true
        t += 1
        new = propagate_step(g, burned)

        if all(new)                        # terminou só com a propagação
            push!(states, copy(new))
            break
        end

        base_cnt = [count(u -> new[u], neighbors(g, v)) for v in 1:n]
        best, best_score = 0, (-1, -1)
        for u in 1:n
            new[u] && continue
            score = (closure_size(g, new, base_cnt, u),
                     count(w -> !new[w], neighbors(g, u)))
            if score > best_score
                best, best_score = u, score
            end
        end

        new[best] = true
        push!(sources, (t, best))
        push!(states, copy(new))
        burned = new
        all(new) && break
    end

    return t, sources, states
end


# ---------------------------------------------------------------------
# Construção do modelo
# ---------------------------------------------------------------------

"""
    build_2burning_model(g::AbstractGraph, T::Int;
                         silent::Bool=false,
                         time_limit::Union{Nothing,Real}=nothing)

Constrói o modelo PLI (R1-R11) para o grafo `g` com horizonte de tempo `T`.
`T` deve ser um limitante superior válido para o número de queima
(por exemplo, T = nv(g), ou o valor obtido pelo algoritmo guloso).
Se `time_limit` for informado, o solver para após esse número de segundos.

Retorna o `Model` do JuMP pronto para ser otimizado.
"""
function build_2burning_model(g::AbstractGraph, T::Int;
                              silent::Bool=false,
                              time_limit::Union{Nothing,Real}=nothing)
    V = collect(vertices(g))

    model = Model(CPLEX.Optimizer)
    silent && set_silent(model)
    time_limit !== nothing && set_time_limit_sec(model, Float64(time_limit))

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
# Warm start
# ---------------------------------------------------------------------

"""
    set_warm_start!(model, g, T, z_g, sources, states)

Fornece ao modelo a solução gulosa como solução inicial (MIP start),
preenchendo TODAS as variáveis (x, s, f, z). Requer T >= z_g.
Para t > z_g, o estado fica todo queimado (x = 1) e s = 0.
"""
function set_warm_start!(model, g::AbstractGraph, T::Int, z_g::Int,
                         sources::Vector{Tuple{Int,Int}},
                         states::Vector{<:AbstractVector{Bool}})
    T >= z_g || error("Horizonte T=$T menor que o valor guloso z=$z_g")

    x, s, f, z = model[:x], model[:s], model[:f], model[:z]

    for v in vertices(g), t in 0:T
        tt = min(t, z_g)
        set_start_value(x[v, t], states[tt+1][v] ? 1.0 : 0.0)
    end

    for v in vertices(g), t in 1:T
        set_start_value(s[v, t], 0.0)
    end
    for (t, v) in sources
        set_start_value(s[v, t], 1.0)
    end

    for t in 1:T
        set_start_value(f[t], t == z_g ? 1.0 : 0.0)
    end

    set_start_value(z, float(z_g))
    return nothing
end

# ---------------------------------------------------------------------
# Extração da solução
# ---------------------------------------------------------------------

"""
    extract_solution(model, g)

Após `optimize!(model)`, extrai:
  - b2 : valor de z (ótimo, ou melhor solução encontrada se houve limite de tempo)
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
# Saída em CSV
# ---------------------------------------------------------------------

"""
    init_csv(path)

Cria o arquivo CSV com o cabeçalho, caso ele ainda não exista.
Se já existir, as novas linhas são anexadas ao final.
"""
function init_csv(path::AbstractString)
    if !isfile(path)
        open(path, "w") do io
            println(io, "graph,greedy_burning_number,pli_burning_number,time,#sources,status")
        end
    end
end

"""
    append_csv_row(path, graph_name, burning_number, time, status, greedy_bn, sources)

Anexa uma linha ao CSV. `burning_number` pode ser `missing` quando o solver
não encontrou nenhuma solução viável (campo fica vazio). `greedy_bn` é o valor
obtido pelo algoritmo guloso.
"""
function append_csv_row(path::AbstractString, graph_name::AbstractString,
                        burning_number, time::Real, status, greedy_bn, sources)
    bn = ismissing(burning_number) ? "" : string(burning_number)
    open(path, "a") do io
        println(io, string(graph_name, ",", greedy_bn, ",", bn, ",", round(time; digits = 3), ",", length(sources), ",",
                           status))
    end
end

# ---------------------------------------------------------------------
# Resolução de uma instância
# ---------------------------------------------------------------------

"""
    resolve_instancia(name, g; time_limit, csv_path, warm_start, shrink_horizon)

Resolve o PLI para `g`, opcionalmente com warm start guloso e horizonte
reduzido, e grava uma linha no CSV.
Retorna (b2, sources) ou (missing, []) se não houve solução viável.

Obs.: a coluna `time` do CSV é o tempo do solver (`solve_time`); o tempo do
guloso não está incluído (é desprezível perto do PLI, mas pode ser medido à parte).
"""
function resolve_instancia(name::AbstractString, g::AbstractGraph;
                           time_limit::Real = TIME_LIMIT,
                           csv_path::AbstractString = CSV_PATH,
                           warm_start::Bool = WARM_START,
                           shrink_horizon::Bool = SHRINK_HORIZON)
    n = nv(g)
    T = n  # limitante superior simples

    # O guloso sempre roda: seu valor é gravado no CSV mesmo com as flags desligadas.
    z_g, src_g, states_g = greedy_2burning(g)
    println("[", name, "] guloso: z = ", z_g)
    shrink_horizon && (T = z_g)

    model = build_2burning_model(g, T; silent = SILENT, time_limit = time_limit)
    warm_start && set_warm_start!(model, g, T, z_g, src_g, states_g)

    optimize!(model)

    status = termination_status(model)
    tempo  = solve_time(model)

    b2, sources = missing, Tuple{Int,Int}[]
    if has_values(model)
        b2, sources = extract_solution(model, g)
    end

    #append_csv_row(csv_path, name, z_g, b2, tempo, status)
    append_csv_row(csv_path, name, b2, tempo, status, z_g, sources)
    println("[", name, "] status = ", status, " | b2 = ", b2,
            " | tempo = ", round(tempo; digits = 2), " s")

    return b2, sources
end


"""
Percorre `input_dir` recursivamente, resolve cada .txt e grava um CSV por grafo
em `output_dir`, espelhando a estrutura de subpastas. Pula grafos já processados.
"""
function processa_diretorio(input_dir::AbstractString, output_dir::AbstractString)
    arquivos = String[]
    for (root, _, files) in walkdir(input_dir)
        for f in files
            endswith(f, ".txt") && push!(arquivos, joinpath(root, f))
        end
    end
    sort!(arquivos)  # ordem determinística

    total = length(arquivos)
    for (i, filepath) in enumerate(arquivos)
        rel_dir = relpath(dirname(filepath), input_dir)
        out_sub = rel_dir == "." ? output_dir : joinpath(output_dir, rel_dir)
        name    = splitext(basename(filepath))[1]
        out_csv = joinpath(out_sub, name * ".csv")

        println("[", i, "/", total, "] ", filepath)

        # Só considera "processado" se tiver cabeçalho + pelo menos uma linha
        if isfile(out_csv) && countlines(out_csv) >= 2
            println("  Já processado, pulando -> ", out_csv)
            continue
        end

        mkpath(out_sub)
        try
            g = read_simple_graph(filepath)
            isfile(out_csv) && rm(out_csv)
            init_csv(out_csv)
            resolve_instancia(name, g; csv_path = out_csv)
        catch e
            @error "Falha em $filepath" exception = (e, catch_backtrace())
            # remove CSV incompleto (só cabeçalho) para ser refeito depois
            isfile(out_csv) && countlines(out_csv) < 2 && rm(out_csv)
        end

        GC.gc()  # libera memória entre instâncias grandes
    end
end

# ---------------------------------------------------------------------
# Função principal
# ---------------------------------------------------------------------
function main()
    if length(ARGS) < 2
        println(stderr, "Uso: julia burning2_pli_warmstart.jl <dir_entrada> <dir_saida>")
        exit(1)
    end
    input_dir, output_dir = ARGS[1], ARGS[2]
    isdir(input_dir) || (println(stderr, "Diretório não encontrado: ", input_dir); exit(1))
    processa_diretorio(input_dir, output_dir)
end

# ---------------------------------------------------------------------
# Chamada
# ---------------------------------------------------------------------
main()